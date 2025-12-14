## =============================================================================
## m6A prediction + survival + matched mRNA survival + m6A–mRNA correlation
## (clean, English-only, copy/paste ready)
## =============================================================================

.libPaths("/rsrch5/home/bcb/yzhao19/R/ubuntu/4.4.1")
setwd("/home/yzhao19/m6a/Model/")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(survival)
  library(broom)
  library(ggplot2)
})

set.seed(42)

## =============================================================================
## 0) Inputs
## =============================================================================
m6a_matrix_norm <- readRDS("Data/m6A_matrix_226.rds")
colnames(m6a_matrix_norm)[1:2] <- c("gene_id", "gene_name")

TCGA_mRNA <- readRDS("Data/TCGA_mRNA_QN_zscore.rds")  # must contain column "gene_id"

survival_all <- read.csv(
  "survival/survival.csv",
  header = TRUE, check.names = FALSE, stringsAsFactors = FALSE
)[, 1:4]

external_prediction_mlp <- fread(
  "Prediction/New_1/external_3485_predictions_v16_10152025.csv",
  select = c("peak_id", "sample_id", "pred_Meta_final", "CancerType")
) %>%
  as.data.frame() %>%
  mutate(predicted_m6A = pred_Meta_final)

## =============================================================================
## 1) peak_id -> gene_name mapping
## =============================================================================
gene_map <- m6a_matrix_norm %>%
  tibble::rownames_to_column("peak_id") %>%
  dplyr::select(peak_id, gene_name)

## =============================================================================
## 2) Merge m6A prediction with survival (this will be used for m6A KM + plotting)
## =============================================================================
merged_data <- external_prediction_mlp %>%
  mutate(bcr_patient_barcode = substr(sample_id, 1, 12)) %>%
  left_join(
    survival_all %>% dplyr::select(bcr_patient_barcode, type, OS, OS.time),
    by = "bcr_patient_barcode"
  ) %>%
  mutate(
    OS.time = as.numeric(as.character(OS.time)),
    OS      = as.numeric(as.character(OS))
  ) %>%
  filter(!is.na(OS.time), !is.na(OS))

## =============================================================================
## 3) mRNA long table matched to the SAME samples as m6A survival cohort
## =============================================================================
stopifnot("gene_id" %in% colnames(TCGA_mRNA))

mrna_all_samples <- setdiff(colnames(TCGA_mRNA), "gene_id")
common_samples   <- intersect(mrna_all_samples, unique(merged_data$sample_id))
if (length(common_samples) == 0) stop("No overlapping samples between TCGA_mRNA and m6A sample_id.")

TCGA_mRNA_common <- TCGA_mRNA[, c("gene_id", common_samples), drop = FALSE] %>%
  mutate(gene_name = gene_id)

target_genes <- unique(gene_map$gene_name)
TCGA_mRNA_common <- TCGA_mRNA_common %>% filter(gene_name %in% target_genes)

# Collapse duplicated gene_name rows by mean (recommended)
TCGA_mRNA_common <- TCGA_mRNA_common %>%
  group_by(gene_name) %>%
  summarise(across(all_of(common_samples), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

mrna_long <- TCGA_mRNA_common %>%
  pivot_longer(
    cols = all_of(common_samples),
    names_to = "sample_id",
    values_to = "mrna_expr"
  ) %>%
  mutate(bcr_patient_barcode = substr(sample_id, 1, 12)) %>%
  left_join(
    survival_all %>% dplyr::select(bcr_patient_barcode, type, OS, OS.time),
    by = "bcr_patient_barcode"
  ) %>%
  mutate(
    OS.time = as.numeric(as.character(OS.time)),
    OS      = as.numeric(as.character(OS))
  ) %>%
  filter(!is.na(type), !is.na(OS.time), !is.na(OS))

## =============================================================================
## 4) mRNA KM per gene per cancer type (median split within cancer)
## =============================================================================
min_n_per_arm <- 5

mrna_split <- mrna_long %>%
  group_by(gene_name, type) %>%
  mutate(group = ifelse(mrna_expr > median(mrna_expr, na.rm = TRUE), "High", "Low")) %>%
  ungroup()

safe_km_row <- function(df, min_n_per_arm = 5) {
  grp_tab <- table(df$group)
  valid_arms <- length(grp_tab) >= 2
  enough_n   <- all(grp_tab >= min_n_per_arm)
  any_event  <- any(df$OS == 1, na.rm = TRUE)
  
  if (!valid_arms || !enough_n || !any_event) {
    return(tibble(
      chisq = NA_real_, p.value = NA_real_,
      med_H = NA_real_, med_L = NA_real_,
      mrna_better = NA_character_
    ))
  }
  
  fit <- tryCatch(survfit(Surv(OS.time, OS) ~ group, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble(
      chisq = NA_real_, p.value = NA_real_,
      med_H = NA_real_, med_L = NA_real_,
      mrna_better = NA_character_
    ))
  }
  
  chisq <- tryCatch(survdiff(Surv(OS.time, OS) ~ group, data = df)$chisq, error = function(e) NA_real_)
  pval  <- if (is.na(chisq)) NA_real_ else pchisq(chisq, df = 1, lower.tail = FALSE)
  
  tb <- tryCatch(summary(fit)$table, error = function(e) NULL)
  
  get_median <- function(tb, label) {
    if (is.null(tb) || is.null(dim(tb))) return(NA_real_)
    idx <- grep(label, rownames(tb), ignore.case = TRUE)
    if (length(idx) == 0) return(NA_real_)
    as.numeric(tb[idx, "median"])
  }
  
  med_H <- get_median(tb, "High")
  med_L <- get_median(tb, "Low")
  
  mrna_better <- case_when(
    !is.na(med_H) & !is.na(med_L) & med_H > med_L ~ "High",
    !is.na(med_H) & !is.na(med_L) & med_L > med_H ~ "Low",
    TRUE ~ NA_character_
  )
  
  tibble(chisq = chisq, p.value = pval, med_H = med_H, med_L = med_L, mrna_better = mrna_better)
}

km_mrna <- mrna_split %>%
  group_by(gene_name, type) %>%
  nest() %>%
  mutate(stats = map(data, ~ safe_km_row(.x, min_n_per_arm = min_n_per_arm))) %>%
  unnest(stats) %>%
  ungroup() %>%
  group_by(gene_name) %>%
  mutate(fdr_gene = p.adjust(p.value, method = "fdr")) %>%
  ungroup() %>%
  transmute(
    gene_name, type,
    p.value, fdr_gene,
    High = med_H, Low = med_L,
    mrna_better
  )

saveRDS(km_mrna, "survival/km_mrna.rds")

## =============================================================================
## 5) m6A–mRNA correlation per (peak, cancer type), using matched samples
## =============================================================================
m6a_with_gene <- merged_data %>%
  left_join(gene_map, by = "peak_id") %>%
  filter(sample_id %in% common_samples) %>%
  filter(!is.na(gene_name))

merged_for_corr <- m6a_with_gene %>%
  left_join(
    mrna_long %>% dplyr::select(sample_id, gene_name, type, mrna_expr),
    by = c("sample_id", "gene_name", "type")
  ) %>%
  filter(!is.na(mrna_expr), !is.na(predicted_m6A))

corr_df <- merged_for_corr %>%
  group_by(peak_id, gene_name, type) %>%
  summarise(
    n    = n(),
    corr = suppressWarnings(cor(predicted_m6A, mrna_expr, use = "pairwise.complete.obs", method = "pearson")),
    .groups = "drop"
  ) %>%
  filter(is.finite(corr), n >= 10)

saveRDS(corr_df, "survival/m6a_mrna_corr_matched_samples.rds")

## =============================================================================
## 6) Direction-consistent significant pairs (mRNA + m6A KM + correlation)
## =============================================================================
if (!exists("km_peak")) {
  km_peak <- readRDS("survival/km_peak_results_10162025.rds")
}
stopifnot(all(c("peak_id", "gene_name", "type", "fdr_peak", "better") %in% colnames(km_peak)))

sig_mrna <- km_mrna %>%
  filter(!is.na(mrna_better), fdr_gene <= 0.1) %>%
  rename(mrna_p = p.value)

sig_m6a <- km_peak %>%
  filter(!is.na(better), fdr_peak <= 0.1) %>%
  transmute(peak_id, gene_name, type, m6a_fdr = fdr_peak, m6a_better = better)

opposite_dir <- function(x) ifelse(x == "High", "Low", ifelse(x == "Low", "High", NA_character_))

final_hits <- corr_df %>%
  inner_join(sig_mrna, by = c("gene_name", "type")) %>%
  inner_join(sig_m6a,  by = c("peak_id", "gene_name", "type")) %>%
  mutate(
    inferred_m6a_better = ifelse(corr >= 0, mrna_better, opposite_dir(mrna_better)),
    direction_match     = (m6a_better == inferred_m6a_better)
  ) %>%
  filter(direction_match) %>%
  arrange(gene_name, type, desc(abs(corr)))

final_summary <- final_hits %>%
  transmute(
    peak_id, gene_name, type,
    n, corr,
    mrna_p,
    mrna_fdr = fdr_gene, mrna_better,
    m6a_fdr,  m6a_better,
    inferred_m6a_better
  )

write.csv(final_summary, "survival/final_consistent_hits_summary.csv", row.names = FALSE)
print(head(final_summary, 20))

## =============================================================================
## 7) KM plot helper for a given peak + cancer type (prints group sizes)
## =============================================================================
plot_peak_km_simple2 <- function(peak_id, gene, cancer_type) {
  df <- merged_data %>%
    filter(peak_id == !!peak_id, type == !!cancer_type) %>%
    mutate(group = ifelse(predicted_m6A > median(predicted_m6A, na.rm = TRUE), "High", "Low"))
  
  cat("\n==============================\n")
  cat("Peak:", peak_id, "\nGene:", gene, "\nCancer:", cancer_type, "\n")
  cat("Sample size:\n")
  print(df %>% count(group))
  cat("==============================\n\n")
  
  fit     <- survfit(Surv(OS.time, OS) ~ group, data = df)
  surv_df <- tidy(fit) %>% mutate(group = sub("^group=", "", strata))
  
  cross_df <- surv_df %>%
    group_by(group) %>%
    summarise(
      cross_time = suppressWarnings(min(time[estimate <= 0.5])),
      estimate   = 0.5,
      .groups    = "drop"
    ) %>%
    filter(is.finite(cross_time))
  
  ggplot(surv_df, aes(x = time, y = estimate, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, colour = NA) +
    geom_step(linewidth = 1.2) +
    geom_segment(
      data = cross_df,
      aes(x = cross_time, xend = cross_time, y = 0, yend = estimate),
      colour = "black", linetype = "dashed"
    ) +
    geom_segment(
      data = cross_df,
      aes(x = 0, xend = cross_time, y = estimate, yend = estimate),
      colour = "black", linetype = "dashed"
    ) +
    scale_colour_manual(values = c(Low = "#4E669A", High = "#B8474D")) +
    scale_fill_manual(values  = c(Low = "#4E669A", High = "#B8474D")) +
    labs(
      title  = paste0(gene, " - ", cancer_type),
      x      = "Time (days)",
      y      = "Survival Probability",
      colour = NULL,
      fill   = NULL
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.95, 0.95),
      legend.justification = c("right", "top"),
      legend.background = element_blank()
    )
}

## Example:
# plot_peak_km_simple2("consensus_peak-chrX_101162911_101163112_+", gene = "CENPI", cancer_type = "LUAD")
