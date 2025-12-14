# ==============================================================================
# RPPA regulator: add CancerType, QC plots, and multi-omic variation analyses
# English-only + paste-ready
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(data.table)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(gridExtra)
  library(ggrepel)
  library(patchwork)
})

# ------------------------------------------------------------------------------
# Paths (edit if needed)
# ------------------------------------------------------------------------------
PATH_RPPA_REG_IN   <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/Regulator/rppa_regulator.rds"
PATH_SURVIVAL_CSV  <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/TCGA_data/RPPA/survival.csv"
PATH_RPPA_REG_OUT  <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/Regulator/rppa_regulator_10022025.rds"

PATH_TCGA_COLOR    <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/TCGA_Cancer_Type_color.csv"

PATH_GISTIC_CN     <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/Regulator/20250304_gistic_reg_gene_level_thresh_rppa_samples.rds"

PATH_MAF           <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/TCGA_data/mc3.v0.2.8.PUBLIC.GRCh38_converted_09022025.maf"

PATH_METH_TSV      <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/Regulator/selected_regulators_silencing_status.tsv"

BASE_GISTIC_AMP    <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/1107/Gistic2/Amplifications/"
BASE_GISTIC_DEL    <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/1107/Gistic2/Deletions/"

BASE_MUTSIG2CV     <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/1107/MutSig2CV/"

# ------------------------------------------------------------------------------
# Regulator gene list
# ------------------------------------------------------------------------------
REGULATORS <- c(
  "ALKBH5", "FTO", "IGF2BP1", "IGF2BP2", "IGF2BP3", "KIAA1429",
  "METTL14", "METTL3", "RBM15", "WTAP", "YTHDC1", "YTHDC2",
  "YTHDF1", "YTHDF2", "YTHDF3"
)

# Allow alias mapping for mutation symbols (e.g., VIRMA for KIAA1429)
GENE_ALIASES <- list(
  KIAA1429 = c("KIAA1429", "VIRMA")
)
REGULATORS_WITH_ALIASES <- unlist(lapply(REGULATORS, function(g) {
  if (g %in% names(GENE_ALIASES)) GENE_ALIASES[[g]] else g
}))

# ==============================================================================
# 1) Load RPPA regulator matrix and map CancerType
# ==============================================================================
rppa_regulator <- readRDS(PATH_RPPA_REG_IN)

# Drop any existing CancerType columns if present (safe)
if ("CancerType" %in% colnames(rppa_regulator)) rppa_regulator$CancerType <- NULL
if ("type" %in% colnames(rppa_regulator))       rppa_regulator$type <- NULL

overlap_samples <- rownames(rppa_regulator)
cat("RPPA regulator samples:", length(overlap_samples), "\n")
print(head(rppa_regulator))

# ---- clinical survival table: 12-char barcode -> type (CancerType) ----
clinical_data <- read.csv(PATH_SURVIVAL_CSV)
clinical_data <- clinical_data[, 1:4]
clinical_data$bcr_patient_barcode <- substr(clinical_data$bcr_patient_barcode, 1, 12)

barcode_map <- clinical_data %>%
  filter(!is.na(bcr_patient_barcode), bcr_patient_barcode != "") %>%
  select(bcr_patient_barcode, type) %>%
  distinct(bcr_patient_barcode, .keep_all = TRUE)

# Add CancerType to rppa_regulator by 12-char patient barcode
rppa_regulator_with_ct <- rppa_regulator %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  mutate(bcr_patient_barcode = substr(sample_id, 1, 12)) %>%
  left_join(barcode_map, by = "bcr_patient_barcode") %>%
  mutate(CancerType = type) %>%
  select(-bcr_patient_barcode, -type) %>%
  column_to_rownames("sample_id")

cat("Unmatched CancerType after clinical join:", sum(is.na(rppa_regulator_with_ct$CancerType)), "\n")

# ---- Optional fallback mapping using Sample_Source (requires `filtered_data`) ----
# Expected: filtered_data has columns: sample_id, Sample_Source
# This block only runs if `filtered_data` exists in your environment.
if (exists("filtered_data")) {
  test_na <- rppa_regulator_with_ct %>% as.data.frame() %>% filter(is.na(CancerType))
  test_sample <- rownames(test_na)
  
  test_filtered <- filtered_data %>% filter(sample_id %in% test_sample)
  
  src2ct <- c(
    "TCGA (BRCA)"     = "BRCA",
    "TCGA (LGG)"      = "LGG",
    "TCGA (Melanoma)" = "SKCM",
    "TCGA (OVCA)"     = "OV",
    "TCGA (STAD)"     = "STAD",
    "TCGA (TGCT)"     = "TGCT"
  )
  
  sample_src_map <- test_filtered %>%
    as.data.frame() %>%
    select(sample_id, Sample_Source) %>%
    mutate(CancerType_from_src = unname(src2ct[Sample_Source])) %>%
    distinct(sample_id, .keep_all = TRUE)
  
  rppa_regulator_filled <- rppa_regulator_with_ct %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    left_join(sample_src_map, by = "sample_id") %>%
    mutate(CancerType = coalesce(CancerType, CancerType_from_src)) %>%
    select(-CancerType_from_src) %>%
    column_to_rownames("sample_id")
  
  filled_now <- sum(is.na(rppa_regulator_with_ct$CancerType) & !is.na(rppa_regulator_filled$CancerType))
  cat("Filled CancerType from Sample_Source for samples:", filled_now, "\n")
  cat("Remaining NA CancerType:", sum(is.na(rppa_regulator_filled$CancerType)), "\n")
  
  rppa_regulator <- rppa_regulator_filled
} else {
  rppa_regulator <- rppa_regulator_with_ct
}

# Save updated
saveRDS(rppa_regulator, PATH_RPPA_REG_OUT)
cat("Saved:", PATH_RPPA_REG_OUT, "\n")

# Reload (optional)
rppa_regulator <- readRDS(PATH_RPPA_REG_OUT)
overlap_samples <- rownames(rppa_regulator)
cat("Final RPPA regulator samples:", length(overlap_samples), "\n")
print(head(rppa_regulator))

# ==============================================================================
# 2) Figure 5A: RPPA samples by CancerType (bar plots)
# ==============================================================================
TCGA_Cancer_Type_color <- read.csv(PATH_TCGA_COLOR)

cancer_counts <- rppa_regulator %>% count(CancerType, name = "n")

cancer_counts_color <- cancer_counts %>%
  left_join(TCGA_Cancer_Type_color, by = c("CancerType" = "cancer_type"))

# Horizontal bar plot
ggplot(cancer_counts_color, aes(x = reorder(CancerType, n), y = n, fill = cancer_type_color)) +
  geom_col(alpha = 0.85) +
  scale_fill_identity() +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

# Vertical bar plot
ggplot(cancer_counts_color, aes(x = reorder(CancerType, -n), y = n, fill = cancer_type_color)) +
  geom_col(alpha = 0.85) +
  scale_fill_identity() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

# ==============================================================================
# 3) Data preparation: GISTIC CN (gene-level threshold), Mutation MAF, Methylation
# ==============================================================================
# ---- GISTIC CN ----
gistic_cn <- readRDS(PATH_GISTIC_CN)
gistic_cn$sample_abbrv <- gsub("_", "-", gistic_cn$sample_abbrv)

gistic_cn <- gistic_cn %>%
  filter(sample_abbrv %in% overlap_samples)

cat("GISTIC CN samples:", length(unique(gistic_cn$sample_abbrv)), "\n")

# ---- Mutation MAF ----
mutation <- fread(PATH_MAF, sep = "\t", quote = "", showProgress = FALSE) %>% as.data.frame()

mutation_1 <- mutation[, c(
  "Hugo_Symbol", "Chromosome", "Variant_Classification", "Tumor_Sample_Barcode",
  "HGVSp_Short", "t_ref_count", "t_alt_count"
)]
mutation_1$sample_id <- substr(mutation_1$Tumor_Sample_Barcode, 1, 16)

# read depth filter
mutation_1 <- mutation_1 %>% filter(t_ref_count + t_alt_count > 60)

# keep only regulators (with aliases), remove silent
filtered_mutation <- mutation_1 %>%
  filter(Hugo_Symbol %in% REGULATORS_WITH_ALIASES) %>%
  filter(Variant_Classification != "Silent") %>%
  distinct()

cat("Filtered mutation rows:", nrow(filtered_mutation), "\n")
print(table(filtered_mutation$Hugo_Symbol))

# Add CancerType to mutation: 12-char from clinical first, then fill by rppa_regulator
clinical_data_renamed <- clinical_data
clinical_data_renamed$sample_id_12 <- clinical_data_renamed$bcr_patient_barcode

filtered_mutation$sample_id_12 <- substr(filtered_mutation$sample_id, 1, 12)
filtered_mutation <- filtered_mutation %>%
  left_join(clinical_data_renamed[, c("sample_id_12", "type")], by = "sample_id_12") %>%
  rename(CancerType = type) %>%
  select(-sample_id_12)

na_idx <- which(is.na(filtered_mutation$CancerType))
if (length(na_idx) > 0) {
  sample2type_rppa <- rppa_regulator$CancerType
  names(sample2type_rppa) <- rownames(rppa_regulator)
  filtered_mutation$CancerType[na_idx] <- sample2type_rppa[filtered_mutation$sample_id[na_idx]]
}
cat("Unmatched mutation CancerType:", sum(is.na(filtered_mutation$CancerType)), "\n")

# ---- Methylation silencing status ----
methylation <- read.table(PATH_METH_TSV, sep = "\t", header = TRUE, check.names = FALSE) %>%
  filter(sample_type != "Solid Tissue Normal") %>%
  select(-meth_sample_id, -exp_sample_id) %>%
  distinct() %>%
  group_by(sample) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  as.data.frame()

rownames(methylation) <- methylation$sample

overlap_samples_met <- intersect(overlap_samples, rownames(methylation))
methylation <- methylation[overlap_samples_met, , drop = FALSE]

cat("Methylation overlap samples:", length(overlap_samples_met), "\n")
print(head(methylation))

# ==============================================================================
# 4) Methylation silencing proportion per cancer (stacked vertical panels)
# ==============================================================================
plot_silenced_prop <- function(gene_status_col, gene_name, color_df = TCGA_Cancer_Type_color) {
  color_map <- color_df %>% select(cancer_type, cancer_type_color)
  
  df_prop <- methylation %>%
    mutate(cancer_type = sub("^TCGA-", "", as.character(project_id))) %>%
    filter(!is.na(.data[[gene_status_col]])) %>%
    group_by(cancer_type) %>%
    summarise(
      total    = n(),
      silenced = sum(.data[[gene_status_col]] == "Silenced", na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(prop = silenced / total) %>%
    left_join(color_map, by = "cancer_type") %>%
    mutate(fill_color = ifelse(is.na(cancer_type_color), "#CCCCCC", cancer_type_color))
  
  ggplot(df_prop, aes(x = reorder(cancer_type, -prop), y = prop, fill = fill_color)) +
    geom_col(width = 0.75) +
    scale_fill_identity() +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
    labs(title = NULL, x = NULL, y = gene_name) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid  = element_blank(),
      axis.line   = element_line(color = "#555555")
    )
}

p_m1 <- plot_silenced_prop("IGF2BP1_status", "IGF2BP1")
p_m2 <- plot_silenced_prop("IGF2BP2_status", "IGF2BP2")
p_m3 <- plot_silenced_prop("IGF2BP3_status", "IGF2BP3")
grid.arrange(p_m1, p_m2, p_m3, ncol = 1)

# ==============================================================================
# 5) Build variation dot plot inputs: Amp/Del from GISTIC peaks, Mut from MutSig2CV
# ==============================================================================
read_gistic_peak_table <- function(file_path, cancer_label, genes_keep = REGULATORS) {
  if (!file.exists(file_path) || file.info(file_path)$size <= 0) {
    return(data.frame())
  }
  raw <- read.table(file_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(raw) == 0) return(data.frame())
  
  # Transpose so first row becomes header-like row
  dt <- as.data.frame(t(raw), stringsAsFactors = FALSE)
  colnames(dt) <- dt[1, ]
  dt <- dt[-1, , drop = FALSE]
  
  # Standardize first 5 columns if present
  if (ncol(dt) >= 5) {
    colnames(dt)[1:5] <- c("cytoband", "q_value", "residual_q_value", "wide_peak_boundaries", "genes_in_wide_peak")
  } else {
    return(data.frame())
  }
  
  dt$q_value <- as.numeric(dt$q_value)
  
  out <- dt %>%
    select(cytoband, q_value, residual_q_value, wide_peak_boundaries, genes_in_wide_peak) %>%
    mutate(gene = str_split(genes_in_wide_peak, pattern = "\n")) %>%
    unnest(gene) %>%
    mutate(gene = str_trim(gene)) %>%
    filter(!is.na(gene), gene != "") %>%
    filter(gene %in% genes_keep) %>%
    distinct(cytoband, q_value, residual_q_value, wide_peak_boundaries, gene) %>%
    mutate(CancerType = cancer_label) %>%
    rename(genes_in_wide_peak = gene)
  
  as.data.frame(out)
}

cancer_types_tcga <- unique(methylation$project_id)

# Amplifications
amplifications_cancers <- lapply(cancer_types_tcga, function(cancer) {
  fp <- file.path(BASE_GISTIC_AMP, paste0(cancer, ".txt"))
  read_gistic_peak_table(fp, cancer_label = cancer, genes_keep = REGULATORS)
})
names(amplifications_cancers) <- cancer_types_tcga

amplifications_cancers_df <- bind_rows(amplifications_cancers)
amplifications_cancers_df <- amplifications_cancers_df %>%
  filter(!is.na(q_value)) %>%
  filter(q_value < 0.2)

cat("Amplifications rows:", nrow(amplifications_cancers_df), "\n")

# Deletions
deletions_cancers <- lapply(cancer_types_tcga, function(cancer) {
  fp <- file.path(BASE_GISTIC_DEL, paste0(cancer, ".txt"))
  read_gistic_peak_table(fp, cancer_label = cancer, genes_keep = REGULATORS)
})
names(deletions_cancers) <- cancer_types_tcga

deletions_cancers_df <- bind_rows(deletions_cancers)
deletions_cancers_df <- deletions_cancers_df %>%
  filter(!is.na(q_value)) %>%
  filter(q_value < 0.2)

cat("Deletions rows:", nrow(deletions_cancers_df), "\n")

# MutSig2CV
mutation_MutSig2CV_list <- lapply(cancer_types_tcga, function(cancer) {
  fp <- file.path(BASE_MUTSIG2CV, paste0(cancer, ".csv"))
  if (!file.exists(fp)) return(data.frame())
  read.csv(fp)
})
names(mutation_MutSig2CV_list) <- cancer_types_tcga

mutation_MutSig2CV_df <- bind_rows(lapply(names(mutation_MutSig2CV_list), function(cancer) {
  df <- mutation_MutSig2CV_list[[cancer]]
  if (nrow(df) == 0) return(NULL)
  df <- df %>% filter(gene %in% REGULATORS)
  if (nrow(df) == 0) return(NULL)
  df$CancerType <- cancer
  df
}))

if (!is.null(mutation_MutSig2CV_df) && nrow(mutation_MutSig2CV_df) > 0) {
  mutation_MutSig2CV_df <- mutation_MutSig2CV_df %>%
    group_by(CancerType) %>%
    mutate(FDR = p.adjust(p, method = "fdr")) %>%
    ungroup() %>%
    filter(FDR < 0.2)
}

cat("MutSig2CV rows after FDR<0.2:", ifelse(is.null(mutation_MutSig2CV_df), 0, nrow(mutation_MutSig2CV_df)), "\n")

# ==============================================================================
# 6) Variation dot plot (Mut / Amp / Del): size=q bin, alpha=percent (capped)
# ==============================================================================
amp_df <- amplifications_cancers_df %>%
  mutate(Type = "Amp") %>%
  select(CancerType, genes_in_wide_peak, q_value, Type)

del_df <- deletions_cancers_df %>%
  mutate(Type = "Del") %>%
  select(CancerType, genes_in_wide_peak, q_value, Type)

mut_df <- mutation_MutSig2CV_df %>%
  mutate(Type = "Mut") %>%
  select(CancerType, gene, FDR, Type) %>%
  rename(genes_in_wide_peak = gene, q_value = FDR)

combined_df0 <- bind_rows(amp_df, del_df, mut_df) %>%
  filter(genes_in_wide_peak %in% REGULATORS) %>%
  mutate(genes_in_wide_peak = factor(genes_in_wide_peak, levels = REGULATORS))

# Count samples per cancer for denominators
gistic_sample_sizes <- gistic_cn %>%
  distinct(tt, sample_abbrv) %>%
  count(tt, name = "n_gistic")

mut_sample_sizes <- rppa_regulator %>%
  count(CancerType, name = "n_mut")

# Counts per event
df_Amp_count <- combined_df0 %>% filter(Type == "Amp") %>%
  rowwise() %>%
  mutate(count = sum(gistic_cn$gene_symbol == genes_in_wide_peak &
                       gistic_cn$tt == CancerType &
                       gistic_cn$value == 2, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(gistic_sample_sizes, by = c("CancerType" = "tt")) %>%
  mutate(percentage = count / n_gistic)

df_Del_count <- combined_df0 %>% filter(Type == "Del") %>%
  rowwise() %>%
  mutate(count = sum(gistic_cn$gene_symbol == genes_in_wide_peak &
                       gistic_cn$tt == CancerType &
                       gistic_cn$value == -2, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(gistic_sample_sizes, by = c("CancerType" = "tt")) %>%
  mutate(percentage = count / n_gistic)

# Use unique mutated samples as count for mutation percentage (recommended)
df_Mut_count <- combined_df0 %>% filter(Type == "Mut") %>%
  rowwise() %>%
  mutate(count = dplyr::n_distinct(filtered_mutation$sample_id[
    filtered_mutation$Hugo_Symbol == as.character(genes_in_wide_peak) &
      filtered_mutation$CancerType == CancerType
  ])) %>%
  ungroup() %>%
  left_join(mut_sample_sizes, by = "CancerType") %>%
  mutate(percentage = count / n_mut)

combined_df <- bind_rows(df_Amp_count, df_Del_count, df_Mut_count) %>%
  mutate(
    q_bin = case_when(
      q_value <= 0.01 ~ "0.01",
      q_value <= 0.05 ~ "0.05",
      q_value <= 0.10 ~ "0.10",
      TRUE            ~ "0.20"
    ),
    perc_cap = pmin(percentage, 0.15)
  ) %>%
  filter(count >= 5) %>%
  distinct()

color_map_var <- c("Mut" = "#187016", "Amp" = "#911E20", "Del" = "#383783")

p_dot <- ggplot(
  combined_df,
  aes(x = CancerType, y = genes_in_wide_peak, color = Type, alpha = perc_cap, size = q_bin)
) +
  geom_point(shape = 16) +
  scale_color_manual(name = "Variation Type", values = color_map_var) +
  scale_alpha_continuous(
    name   = "Variation %",
    limits = c(0, 0.15),
    range  = c(0.5, 1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_size_manual(
    name   = "q-value",
    breaks = c("0.20", "0.10", "0.05", "0.01"),
    values = c(2, 4, 6, 8)
  ) +
  scale_y_discrete(limits = rev(REGULATORS)) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 10),
    panel.grid   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) +
  labs(x = "Cancer Type", y = "Regulator Gene")

print(p_dot)

# ==============================================================================
# 7) Hotspot vs LoF mutation fractions (restricted to significant Mut events)
# ==============================================================================
sig_mut <- combined_df %>% filter(Type == "Mut") %>% distinct(CancerType, genes_in_wide_peak)

matched_mutations <- filtered_mutation %>%
  semi_join(sig_mut, by = c("CancerType" = "CancerType", "Hugo_Symbol" = "genes_in_wide_peak"))

hotspot_types <- c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins")
lof_types <- c("Frame_Shift_Ins", "Frame_Shift_Del", "Nonsense_Mutation",
               "Nonstop_Mutation", "Splice_Site", "Translation_Start_Site")

hotspot_table <- matched_mutations %>%
  filter(Variant_Classification %in% hotspot_types) %>%
  group_by(CancerType, Hugo_Symbol, HGVSp_Short) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 2) %>%
  select(CancerType, Hugo_Symbol, HGVSp_Short)

filtered_labeled <- matched_mutations %>%
  mutate(
    is_hotspot = as.integer(paste(CancerType, Hugo_Symbol, HGVSp_Short) %in%
                              paste(hotspot_table$CancerType, hotspot_table$Hugo_Symbol, hotspot_table$HGVSp_Short)),
    is_lof = as.integer(Variant_Classification %in% lof_types),
    is_nonsilent = as.integer(Variant_Classification %in% c(hotspot_types, lof_types))
  )

mutation_summary <- filtered_labeled %>%
  filter(is_nonsilent == 1) %>%
  group_by(CancerType, Hugo_Symbol) %>%
  summarise(
    n_hotspot = sum(is_hotspot),
    n_lof     = sum(is_lof),
    n_total   = n(),
    fraction_hotspot = n_hotspot / n_total,
    fraction_lof     = n_lof / n_total,
    .groups = "drop"
  ) %>%
  left_join(
    filtered_labeled %>%
      filter(is_hotspot == 1) %>%
      group_by(CancerType, Hugo_Symbol) %>%
      summarise(n_hotspot_sites = n_distinct(HGVSp_Short), .groups = "drop"),
    by = c("CancerType", "Hugo_Symbol")
  ) %>%
  mutate(
    n_hotspot_sites = ifelse(is.na(n_hotspot_sites), 0, n_hotspot_sites),
    cancer_gene = paste(CancerType, Hugo_Symbol, sep = "_"),
    category = case_when(
      fraction_hotspot > 0.20 & fraction_lof < 0.30 & n_hotspot_sites >= 2 ~ "Hotspot enriched",
      fraction_lof > 0.20 & fraction_hotspot < 0.40 & n_lof >= 3 ~ "LoF enriched",
      fraction_hotspot > 0.20 & fraction_lof > 0.20 ~ "Dual enriched",
      TRUE ~ "Non-enriched"
    ),
    size_group = case_when(
      n_total < 10 ~ "n < 10",
      n_total < 25 ~ "10 ≤ n < 25",
      n_total < 50 ~ "25 ≤ n < 50",
      TRUE         ~ "n ≥ 50"
    ),
    size_group = factor(size_group, levels = c("n < 10", "10 ≤ n < 25", "25 ≤ n < 50", "n ≥ 50"))
  )

size_values <- c("n < 10" = 3, "10 ≤ n < 25" = 5, "25 ≤ n < 50" = 7, "n ≥ 50" = 9)

p_hotspot <- ggplot(mutation_summary, aes(x = fraction_lof, y = fraction_hotspot)) +
  geom_point(aes(color = category, alpha = category, size = size_group)) +
  geom_text_repel(
    data = mutation_summary %>% filter(category != "Non-enriched"),
    aes(label = cancer_gene, color = category),
    size = 2.5,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.5,
    segment.alpha = 0.7,
    show.legend = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "Hotspot enriched" = alpha("#911E20", 0.7),
    "LoF enriched"     = alpha("#383783", 0.7),
    "Dual enriched"    = alpha("orange", 0.8),
    "Non-enriched"     = alpha("grey70", 0.8)
  )) +
  scale_alpha_manual(values = c(
    "Hotspot enriched" = 0.9,
    "LoF enriched"     = 0.9,
    "Dual enriched"    = 1,
    "Non-enriched"     = 0.5
  )) +
  scale_size_manual(name = "Non-silent mutation count", values = size_values) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_line(color = "black")) +
  labs(
    x = "Fraction of LoF mutations",
    y = "Fraction of hotspot mutations",
    color = "Enrichment",
    alpha = "Enrichment"
  )

print(p_hotspot)

# ==============================================================================
# 8) RPPA expression comparisons (Wilcoxon) for Mut / Amp / Del and boxplots
# ==============================================================================
compare_rppa <- function(var_samples, wt_samples, gene, rppa_data) {
  common_samples <- intersect(rownames(rppa_data), c(var_samples, wt_samples))
  rppa_sub <- rppa_data[common_samples, , drop = FALSE]
  
  var_values <- rppa_sub[rownames(rppa_sub) %in% var_samples, gene, drop = TRUE]
  wt_values  <- rppa_sub[rownames(rppa_sub) %in% wt_samples,  gene, drop = TRUE]
  
  var_n <- length(var_values)
  wt_n  <- length(wt_values)
  
  if (var_n > 1 && wt_n > 1) {
    p_value <- tryCatch(wilcox.test(var_values, wt_values)$p.value, error = function(e) NA_real_)
    median_diff <- median(var_values, na.rm = TRUE) - median(wt_values, na.rm = TRUE)
    var_med <- median(var_values, na.rm = TRUE)
    wt_med  <- median(wt_values, na.rm = TRUE)
  } else {
    p_value <- NA_real_
    median_diff <- NA_real_
    var_med <- NA_real_
    wt_med  <- NA_real_
  }
  
  list(
    p_value = p_value,
    median_difference = median_diff,
    var_sample_size = var_n,
    wt_sample_size = wt_n,
    var_mean = var_med,
    wt_mean  = wt_med
  )
}

analyze_variation <- function(df, var_condition_func, wt_condition_func, gene_col, rppa_data) {
  df %>%
    rowwise() %>%
    mutate(
      gene = as.character(.data[[gene_col]]),
      cancer = as.character(CancerType),
      var_samples = list(var_condition_func(gene, cancer)),
      wt_samples  = list(wt_condition_func(gene, cancer)),
      result = list(compare_rppa(
        var_samples = unlist(var_samples),
        wt_samples  = unlist(wt_samples),
        gene = gene,
        rppa_data = rppa_data
      )),
      p_value = result$p_value,
      median_difference = result$median_difference,
      var_sample_size = result$var_sample_size,
      wt_sample_size  = result$wt_sample_size,
      var_mean = result$var_mean,
      wt_mean  = result$wt_mean
    ) %>%
    ungroup() %>%
    select(-result, -var_samples, -wt_samples, -gene, -cancer)
}

adjust_fdr <- function(results) {
  results %>% mutate(FDR = p.adjust(p_value, method = "fdr"))
}

# Split dot-plot combined_df into Amp/Del/Mut event tables
df_Amp <- combined_df %>% filter(Type == "Amp")
df_Del <- combined_df %>% filter(Type == "Del")
df_Mut <- combined_df %>% filter(Type == "Mut")

df_Amp$genes_in_wide_peak <- as.character(df_Amp$genes_in_wide_peak)
df_Del$genes_in_wide_peak <- as.character(df_Del$genes_in_wide_peak)
df_Mut$genes_in_wide_peak <- as.character(df_Mut$genes_in_wide_peak)

amp_results <- analyze_variation(
  df_Amp,
  function(gene, cancer) gistic_cn %>% filter(gene_symbol == gene, value == 2,  tt == cancer) %>% pull(sample_abbrv),
  function(gene, cancer) gistic_cn %>% filter(gene_symbol == gene, value == 0,  tt == cancer) %>% pull(sample_abbrv),
  "genes_in_wide_peak",
  rppa_regulator
)

del_results <- analyze_variation(
  df_Del,
  function(gene, cancer) gistic_cn %>% filter(gene_symbol == gene, value == -2, tt == cancer) %>% pull(sample_abbrv),
  function(gene, cancer) gistic_cn %>% filter(gene_symbol == gene, value == 0,  tt == cancer) %>% pull(sample_abbrv),
  "genes_in_wide_peak",
  rppa_regulator
)

mut_results <- analyze_variation(
  df_Mut,
  function(gene, cancer) filtered_mutation %>% filter(Hugo_Symbol == gene, CancerType == cancer) %>% pull(sample_id) %>% unique(),
  function(gene, cancer) {
    setdiff(
      rownames(rppa_regulator[rppa_regulator$CancerType == cancer, , drop = FALSE]),
      filtered_mutation %>% filter(Hugo_Symbol == gene, CancerType == cancer) %>% pull(sample_id) %>% unique()
    )
  },
  "genes_in_wide_peak",
  rppa_regulator
)

amp_results <- adjust_fdr(as.data.frame(amp_results))
del_results <- adjust_fdr(as.data.frame(del_results))
mut_results <- adjust_fdr(as.data.frame(mut_results))

amp_results <- amp_results %>% select(CancerType, genes_in_wide_peak, p_value, FDR, median_difference, var_mean, wt_mean, var_sample_size, wt_sample_size)
del_results <- del_results %>% select(CancerType, genes_in_wide_peak, p_value, FDR, median_difference, var_mean, wt_mean, var_sample_size, wt_sample_size)
mut_results <- mut_results %>% select(CancerType, genes_in_wide_peak, p_value, FDR, median_difference, var_mean, wt_mean, var_sample_size, wt_sample_size)

# ---- Boxplots with custom axis-line style (minimal panels) ----
color_map_box <- c("Mut" = "#187016", "Amp" = "#911E20", "Del" = "#383783", "WT" = "#B0B0B0")

.darkener <- function(cols, factor = 0.75) {
  vapply(cols, function(cl) {
    rgb_val <- grDevices::col2rgb(cl) / 255
    rgb_new <- pmax(0, pmin(1, rgb_val * factor))
    grDevices::rgb(rgb_new[1], rgb_new[2], rgb_new[3])
  }, character(1))
}

generate_mut_plots <- function(results, rppa_data) {
  sig <- results %>% filter(p_value < 0.05)
  if (nrow(sig) == 0) { message("No significant mutation results under p < 0.05"); return(invisible(NULL)) }
  
  plots <- vector("list", nrow(sig))
  for (i in seq_len(nrow(sig))) {
    gene <- sig$genes_in_wide_peak[i]
    cancer <- sig$CancerType[i]
    
    var_samples <- filtered_mutation %>%
      filter(Hugo_Symbol == gene, CancerType == cancer) %>%
      pull(sample_id) %>% unique()
    
    wt_samples <- setdiff(
      rownames(rppa_data[rppa_data$CancerType == cancer, , drop = FALSE]),
      var_samples
    )
    
    var_values <- rppa_data[var_samples, gene, drop = TRUE]
    wt_values  <- rppa_data[wt_samples,  gene, drop = TRUE]
    
    df <- data.frame(
      SampleType = factor(rep(c("Mut", "WT"), c(length(var_values), length(wt_values))),
                          levels = c("Mut", "WT")),
      Expression = c(var_values, wt_values)
    ) %>% filter(!is.na(Expression))
    
    used_levels <- levels(droplevels(df$SampleType))
    cm <- color_map_box[used_levels]
    cm_dark <- .darkener(cm, factor = 0.75)
    
    rng_y <- range(df$Expression, na.rm = TRUE)
    y_breaks <- scales::pretty_breaks(n = 4)(rng_y)
    y_span <- diff(rng_y)
    
    yticks_df <- data.frame(
      y = y_breaks,
      x0 = 0.5 - 0.08,
      x1 = 0.5 - 0.14
    )
    xticks_df <- data.frame(
      x = seq_along(used_levels),
      y0 = rng_y[1] - 0.015 * y_span,
      y1 = rng_y[1] - 0.045 * y_span
    )
    
    plots[[i]] <- ggplot(df, aes(x = SampleType, y = Expression)) +
      geom_boxplot(aes(fill = SampleType), width = 0.3, alpha = 0.8,
                   outlier.shape = NA, linewidth = 0.6, color = "black") +
      geom_point(aes(color = SampleType), size = 1.8, alpha = 0.9,
                 position = position_jitter(width = 0.06, height = 0)) +
      geom_segment(data = yticks_df, aes(x = x0, xend = x1, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = xticks_df, aes(x = x, xend = x, y = y0, yend = y1),
                   inherit.aes = FALSE, linewidth = 0.4) +
      scale_fill_manual(values = cm, drop = FALSE) +
      scale_color_manual(values = cm_dark, drop = FALSE) +
      scale_y_continuous(breaks = y_breaks) +
      coord_cartesian(clip = "off") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 20, 30)
      ) +
      labs(y = "Mean RPPA") +
      ggtitle(paste(gene, "-", cancer))
  }
  
  gridExtra::grid.arrange(grobs = plots, ncol = min(5, length(plots)))
}

generate_cna_plots <- function(results, rppa_data, gistic_cn) {
  sig <- results %>% filter(p_value < 0.05)
  if (nrow(sig) == 0) { message("No significant CNA results under p < 0.05"); return(invisible(NULL)) }
  
  plots <- vector("list", nrow(sig))
  for (i in seq_len(nrow(sig))) {
    gene <- sig$genes_in_wide_peak[i]
    cancer <- sig$CancerType[i]
    
    amp_samples <- gistic_cn %>% filter(gene_symbol == gene, value == 2,  tt == cancer) %>% pull(sample_abbrv)
    del_samples <- gistic_cn %>% filter(gene_symbol == gene, value == -2, tt == cancer) %>% pull(sample_abbrv)
    wt_samples  <- gistic_cn %>% filter(gene_symbol == gene, value == 0,  tt == cancer) %>% pull(sample_abbrv)
    
    amp_values <- rppa_data[amp_samples, gene, drop = TRUE]
    del_values <- rppa_data[del_samples, gene, drop = TRUE]
    wt_values  <- rppa_data[wt_samples,  gene, drop = TRUE]
    
    df <- data.frame(
      SampleType = factor(rep(c("Amp", "WT", "Del"),
                              c(length(amp_values), length(wt_values), length(del_values))),
                          levels = c("Amp", "WT", "Del")),
      Expression = c(amp_values, wt_values, del_values)
    ) %>% filter(!is.na(Expression))
    
    used_levels <- levels(droplevels(df$SampleType))
    cm <- color_map_box[used_levels]
    cm_dark <- .darkener(cm, factor = 0.75)
    
    rng_y <- range(df$Expression, na.rm = TRUE)
    y_breaks <- scales::pretty_breaks(n = 4)(rng_y)
    y_span <- diff(rng_y)
    
    yticks_df <- data.frame(
      y = y_breaks,
      x0 = 0.5 - 0.08,
      x1 = 0.5 - 0.14
    )
    xticks_df <- data.frame(
      x = seq_along(used_levels),
      y0 = rng_y[1] - 0.015 * y_span,
      y1 = rng_y[1] - 0.045 * y_span
    )
    
    plots[[i]] <- ggplot(df, aes(x = SampleType, y = Expression)) +
      geom_boxplot(aes(fill = SampleType), width = 0.3, alpha = 0.8,
                   outlier.shape = NA, linewidth = 0.6, color = "black") +
      geom_point(aes(color = SampleType), size = 1.8, alpha = 0.9,
                 position = position_jitter(width = 0.06, height = 0)) +
      geom_segment(data = yticks_df, aes(x = x0, xend = x1, y = y, yend = y),
                   inherit.aes = FALSE, linewidth = 0.4) +
      geom_segment(data = xticks_df, aes(x = x, xend = x, y = y0, yend = y1),
                   inherit.aes = FALSE, linewidth = 0.4) +
      scale_fill_manual(values = cm, drop = FALSE) +
      scale_color_manual(values = cm_dark, drop = FALSE) +
      scale_y_continuous(breaks = y_breaks) +
      coord_cartesian(clip = "off") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.line.x = element_line(color = "black", linewidth = 0.6),
        axis.line.y = element_line(color = "black", linewidth = 0.6),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.margin = margin(10, 20, 20, 30)
      ) +
      labs(y = "Mean RPPA") +
      ggtitle(paste(gene, "-", cancer))
    
    plots[[i]] <- plots[[i]]
  }
  
  gridExtra::grid.arrange(grobs = plots, ncol = min(4, length(plots)))
}

cat("Generating mutation plots...\n")
generate_mut_plots(mut_results, rppa_regulator)

cat("Generating deletion plots...\n")
generate_cna_plots(del_results, rppa_regulator, gistic_cn)

cat("Generating amplification plots...\n")
generate_cna_plots(amp_results, rppa_regulator, gistic_cn)

# ==============================================================================
# 9) Paired t-test across cancer types: silenced vs not silenced (by cancer mean)
# ==============================================================================
paired_t_test_and_plot <- function(gene_status_col, gene_name, show_y = TRUE) {
  cancer_types <- unique(methylation$project_id)
  
  silenced_means <- numeric(0)
  not_silenced_means <- numeric(0)
  kept_cancers <- character(0)
  
  for (ct in cancer_types) {
    ct_df <- methylation[methylation$project_id == ct, , drop = FALSE]
    sil_samples <- ct_df$sample[ct_df[[gene_status_col]] == "Silenced"]
    ns_samples  <- ct_df$sample[ct_df[[gene_status_col]] == "Not Silenced"]
    
    if (length(sil_samples) >= 3 && length(ns_samples) >= 3) {
      sil_expr <- rppa_regulator[sil_samples, gene_name, drop = TRUE]
      ns_expr  <- rppa_regulator[ns_samples,  gene_name, drop = TRUE]
      
      silenced_means <- c(silenced_means, mean(as.numeric(sil_expr), na.rm = TRUE))
      not_silenced_means <- c(not_silenced_means, mean(as.numeric(ns_expr), na.rm = TRUE))
      kept_cancers <- c(kept_cancers, ct)
    }
  }
  
  if (length(kept_cancers) == 0) return(NULL)
  
  p_val <- t.test(silenced_means, not_silenced_means, paired = TRUE)$p.value
  
  plot_data <- data.frame(
    Expression = c(silenced_means, not_silenced_means),
    Status = factor(rep(c("Silenced", "Not Silenced"), each = length(kept_cancers)),
                    levels = c("Not Silenced", "Silenced")),
    Pair = rep(kept_cancers, 2),
    gene = gene_name
  )
  
  color_map_sil <- c("Silenced" = "#7830A0", "Not Silenced" = "gray70")
  
  ggplot(plot_data, aes(x = Status, y = Expression)) +
    geom_boxplot(aes(fill = Status), width = 0.3, outlier.shape = NA, alpha = 0.7) +
    geom_point(aes(color = Status), size = 2, alpha = 0.8, position = position_jitter(width = 0.05)) +
    geom_line(aes(group = Pair), color = "gray60", linewidth = 0.4, linetype = "dotted") +
    scale_fill_manual(values = color_map_sil) +
    scale_color_manual(values = color_map_sil) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 10)
    ) +
    labs(y = ifelse(show_y, "Mean RPPA", NULL)) +
    facet_grid(. ~ gene) +
    annotate(
      "text", x = 1.5, y = max(plot_data$Expression, na.rm = TRUE) * 1.05,
      label = paste0("p = ", signif(p_val, 3)),
      size = 4, fontface = "bold"
    )
}

plots <- list(
  paired_t_test_and_plot("IGF2BP1_status", "IGF2BP1", show_y = TRUE),
  paired_t_test_and_plot("IGF2BP2_status", "IGF2BP2", show_y = TRUE),
  paired_t_test_and_plot("IGF2BP3_status", "IGF2BP3", show_y = TRUE)
)
plots <- Filter(Negate(is.null), plots)
wrap_plots(plots, ncol = 3)
