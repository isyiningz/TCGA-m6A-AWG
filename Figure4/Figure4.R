# ==============================================================================
# Pan-cancer m6A–mRNA correlation + enrichment + example downstream analyses
# (English-only, formatted, ready to paste)
# ==============================================================================

suppressPackageStartupMessages({
  library(preprocessCore)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(scales)
  library(viridis)
  library(grid)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(msigdbr)
  library(pheatmap)
  library(reshape2)
  library(Matrix)
  library(tidyr)
  library(purrr)
  library(ggpubr)
  library(ggExtra)
  library(gridExtra)
})

# ------------------------------------------------------------------------------
# 0) Paths (edit if needed)
# ------------------------------------------------------------------------------
base_path <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0811/"

rna_tsv_path <- "Figure4/Data/EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv"
rna_qn_z_rds  <- file.path(base_path, "TCGA_mRNA_QN_zscore.rds")

m6a_matrix_rds <- file.path(base_path, "m6A_matrix_norm.rds")
barcode_txt    <- file.path(base_path, "TCGA.m6A.sample_barcodes.Synapse.11302022.txt")

m6a_list_rds    <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/m6A_matrix_norm_list_mrna.rds"
filtered_list_rds <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/filtered_results_list.rds"
filtered_list_pc_rds <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/filtered_results_list_Pancancer.rds"
mean_corr_rds   <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/mean_corr_summary.rds"

gsea_out_rds    <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/pancancer_all_gsea_results_enrich_10012025.rds"
shuffle_null_rds <- "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/shuffle_pearson_1000.rds"

# ------------------------------------------------------------------------------
# 1) RNA preprocessing: quantile normalization + row-wise z-score
#     (If you already have rna_qn_z_rds, you can skip to 1.3)
# ------------------------------------------------------------------------------

# 1.1 Load raw RNA table
TCGA_RNA <- read.table(rna_tsv_path, header = TRUE, sep = "\t", check.names = FALSE)
TCGA_RNA$gene_id <- sub("\\|.*", "", TCGA_RNA$gene_id)
gene_id <- as.character(TCGA_RNA$gene_id)

# 1.2 Quantile normalize + z-score (by gene)
TCGA_mRNA_raw <- TCGA_RNA[, -1, drop = FALSE]
data_matrix <- as.matrix(TCGA_mRNA_raw)
qn <- normalize.quantiles(data_matrix)
rownames(qn) <- rownames(data_matrix)
colnames(qn) <- colnames(data_matrix)
z <- t(scale(t(qn)))
TCGA_mRNA2_normalized <- as.data.frame(z)

# 1.3 Standardize sample IDs to 16 chars with "-"
sample_id_rna <- colnames(TCGA_RNA)[-1]
sample_id_rna <- gsub("\\.", "-", sample_id_rna)
sample_id_rna <- substr(sample_id_rna, 1, 16)

TCGA_mRNA2 <- cbind(gene_id = gene_id, TCGA_mRNA2_normalized)
colnames(TCGA_mRNA2) <- c("gene_id", sample_id_rna)

# Save if desired
# saveRDS(TCGA_mRNA2, rna_qn_z_rds)

# Or load the saved matrix
TCGA_mRNA2 <- readRDS(rna_qn_z_rds)

# ------------------------------------------------------------------------------
# 2) m6A matrix + common genes/samples + m6A list by gene
# ------------------------------------------------------------------------------
m6A_matrix_norm <- readRDS(m6a_matrix_rds)

m6A_genename <- m6A_matrix_norm$`as.character(PeakToGene[rownames(m6A_matrix_norm)])`
m6A_sample16 <- substr(colnames(m6A_matrix_norm), 1, 16)
colnames(m6A_matrix_norm) <- m6A_sample16

common_genes <- intersect(TCGA_mRNA2$gene_id, m6A_genename)
common_id16  <- intersect(substr(colnames(m6A_matrix_norm), 1, 16), colnames(TCGA_mRNA2))

# Subset RNA to common genes and samples
TCGA_mRNA_common_1 <- TCGA_mRNA2[, c("gene_id", common_id16), drop = FALSE]
TCGA_mRNA_common_1 <- TCGA_mRNA_common_1[TCGA_mRNA_common_1$gene_id %in% common_genes, , drop = FALSE]
rownames(TCGA_mRNA_common_1) <- TCGA_mRNA_common_1$gene_id
TCGA_mRNA_common_1 <- TCGA_mRNA_common_1[, -1, drop = FALSE]
TCGA_mRNA_common_1 <- TCGA_mRNA_common_1[common_genes, , drop = FALSE]

# Build (or load) per-gene m6A list (rows=peaks for that gene, cols=samples)
if (file.exists(m6a_list_rds)) {
  m6A_matrix_norm_list <- readRDS(m6a_list_rds)
} else {
  m6A_matrix_norm_list <- vector("list", length(common_genes))
  names(m6A_matrix_norm_list) <- common_genes
  for (g in common_genes) {
    m6A_matrix_norm_list[[g]] <- m6A_matrix_norm[m6A_genename == g, common_id16, drop = FALSE]
  }
  m6A_matrix_norm_list <- m6A_matrix_norm_list[common_genes]
  # saveRDS(m6A_matrix_norm_list, m6a_list_rds)
}

# ------------------------------------------------------------------------------
# 3) Cancer type one-hot encoding (by 16-char sample IDs)
# ------------------------------------------------------------------------------
m6A_sample_barcodes <- read.table(barcode_txt, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
m6A_sample_barcodes <- m6A_sample_barcodes %>%
  filter(!grepl("A933-53$", IP)) %>%
  mutate(sample_16 = substr(IP, 1, 16)) %>%
  select(Cancer_type_abbv, sample_16) %>%
  distinct()

rownames(m6A_sample_barcodes) <- m6A_sample_barcodes$sample_16
m6A_sample_barcodes <- m6A_sample_barcodes[common_id16, , drop = FALSE]

cancer_types <- m6A_sample_barcodes$Cancer_type_abbv
sample_ids   <- m6A_sample_barcodes$sample_16
unique_cancer_types <- unique(cancer_types)

cancer_types_encoded <- matrix(0, nrow = length(sample_ids), ncol = length(unique_cancer_types))
colnames(cancer_types_encoded) <- unique_cancer_types
rownames(cancer_types_encoded) <- sample_ids
for (i in seq_along(sample_ids)) {
  cancer_types_encoded[i, cancer_types[i]] <- 1
}
cancer_types_encoded <- as.data.frame(cancer_types_encoded)

# ------------------------------------------------------------------------------
# 4) Barplot: sample counts by cancer type
# ------------------------------------------------------------------------------
cancer_counts <- colSums(cancer_types_encoded)
cancer_data <- data.frame(
  CancerType = names(cancer_counts),
  CohortSize = as.numeric(cancer_counts),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Color = ifelse(CohortSize >= 3, ">=3", "<3"),
    CancerType = factor(CancerType, levels = CancerType[order(CohortSize, decreasing = TRUE)])
  )

ggplot(cancer_data, aes(x = CancerType, y = CohortSize, fill = Color)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = CohortSize), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    ">=3" = adjustcolor("#0099FF", alpha.f = 0.6),
    "<3"  = adjustcolor("#FF6666", alpha.f = 0.6)
  )) +
  labs(x = NULL, y = "Cohort Size") +
  theme_light() +
  theme(
    axis.text.x = element_text(size = 15, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

# ------------------------------------------------------------------------------
# 5) Per-cancer: compute peak–mRNA correlations; keep common pairs across cancers
# ------------------------------------------------------------------------------
subset_list_by_samples <- function(data_list, keep_ids) {
  lapply(data_list, function(df) df[, colnames(df) %in% keep_ids, drop = FALSE])
}

cancers_keep <- cancer_data %>%
  filter(CohortSize >= 3) %>%
  pull(CancerType) %>%
  as.character()

results_list_correlations <- list()

for (cancer in cancers_keep) {
  cat("Processing cancer:", cancer, "\n")
  
  cancer_subset <- cancer_types_encoded[cancer_types_encoded[[cancer]] == 1, , drop = FALSE]
  common_id_cancer <- rownames(cancer_subset)
  
  m6A_list_sub <- subset_list_by_samples(m6A_matrix_norm_list, common_id_cancer)
  RNA_sub <- TCGA_mRNA_common_1[, common_id_cancer, drop = FALSE]
  
  correlations_df <- data.frame(Gene_m6A_Pair = character(), Correlation = numeric(), stringsAsFactors = FALSE)
  
  for (i in seq_along(common_genes)) {
    gene_name <- common_genes[i]
    if (is.null(m6A_list_sub[[gene_name]]) || nrow(m6A_list_sub[[gene_name]]) == 0) next
    
    mRNA_expression <- as.numeric(RNA_sub[i, , drop = TRUE])
    
    for (m in seq_len(nrow(m6A_list_sub[[gene_name]]))) {
      m6A_modification <- as.numeric(m6A_list_sub[[gene_name]][m, , drop = TRUE])
      
      dat <- data.frame(mRNA_expression = mRNA_expression, m6A_modification = m6A_modification)
      dat <- na.omit(dat)
      
      if (nrow(dat) < 3) next
      
      corr <- suppressWarnings(cor(dat$mRNA_expression, dat$m6A_modification, method = "pearson"))
      m6A_name <- rownames(m6A_list_sub[[gene_name]])[m]
      id <- paste(gene_name, m6A_name, sep = "_")
      
      correlations_df <- rbind(
        correlations_df,
        data.frame(Gene_m6A_Pair = id, Correlation = corr, stringsAsFactors = FALSE)
      )
    }
  }
  
  results_list_correlations[[cancer]] <- correlations_df
}

# Common pairs across all cancers
common_pairs <- lapply(results_list_correlations, function(df) df$Gene_m6A_Pair)
common_gene_m6A_pairs <- Reduce(intersect, common_pairs)

filtered_results_list <- lapply(results_list_correlations, function(df) {
  df[df$Gene_m6A_Pair %in% common_gene_m6A_pairs, , drop = FALSE]
})

# Save/load as needed
# saveRDS(filtered_results_list, "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/filtered_results_list_pearson_QN_Z.rds")
if (file.exists(filtered_list_rds)) filtered_results_list <- readRDS(filtered_list_rds)

# ------------------------------------------------------------------------------
# 6) Pan-cancer mean correlation across cancers (for common pairs)
# ------------------------------------------------------------------------------
all_pairs_correlations <- lapply(filtered_results_list, function(df) df[, c("Gene_m6A_Pair", "Correlation")])
combined_correlations <- Reduce(function(x, y) full_join(x, y, by = "Gene_m6A_Pair"), all_pairs_correlations)
combined_correlations$Mean_Correlation <- rowMeans(combined_correlations[, -1, drop = FALSE], na.rm = TRUE)

mean_corr_summary <- combined_correlations[, c("Gene_m6A_Pair", "Mean_Correlation")]
mean_corr_summary$Gene <- sub("_.*", "", mean_corr_summary$Gene_m6A_Pair)
colnames(mean_corr_summary) <- c("Gene_m6A_Pair", "Correlation", "Gene")

filtered_results_list[["PanCancer"]] <- mean_corr_summary

# Save/load as needed
# saveRDS(filtered_results_list, "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/filtered_results_list_Pancancer.rds")
# saveRDS(mean_corr_summary, "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/0919/mean_corr_summary.rds")
if (file.exists(filtered_list_pc_rds)) filtered_results_list <- readRDS(filtered_list_pc_rds)
if (file.exists(mean_corr_rds)) mean_corr_summary <- readRDS(mean_corr_rds)

# ------------------------------------------------------------------------------
# 7) Figure 4A: ranked genes (PanCancer)
# ------------------------------------------------------------------------------
top_n <- 10
df_rank <- filtered_results_list[["PanCancer"]] %>%
  arrange(desc(Correlation)) %>%
  mutate(rank = row_number())

N <- nrow(df_rank)
label_df <- df_rank %>% filter(rank <= top_n | rank > (N - top_n))

p_rank <- ggplot(df_rank, aes(rank, Correlation)) +
  geom_segment(aes(xend = rank, y = 0, yend = Correlation),
               linewidth = 0.3, colour = "#E0E0E0", alpha = 0.6) +
  geom_point(aes(colour = Correlation),
             size = 3, shape = 21, stroke = 0.4, fill = "white") +
  scale_colour_viridis_c(
    option = "C", end = 0.9,
    guide = guide_colorbar(
      title = "ρ", title.position = "top",
      barwidth = unit(0.4, "cm"),
      barheight = unit(5, "cm"),
      frame.colour = "grey80",
      frame.linewidth = 0.5
    )
  ) +
  geom_text_repel(
    data = label_df, aes(label = Gene),
    size = 4.5,
    box.padding = 0.3, point.padding = 0.2,
    segment.size = 0.25, max.overlaps = Inf
  ) +
  scale_x_reverse(expand = expansion(mult = c(0.02, 0.02)),
                  breaks = scales::pretty_breaks()) +
  labs(
    x = "Gene Rank",
    y = "Mean correlation (ρ)",
    title = "Pan-Cancer **m<sup>6</sup>A** site – gene correlation ranking"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = ggtext::element_markdown(face = "bold", size = 14),
    axis.line = element_line(color = "grey70"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90"),
    axis.text = element_text(color = "grey30"),
    axis.title = element_text(color = "grey20"),
    legend.position = "right",
    legend.title.align = 0.5,
    legend.key.width = unit(0.3, "cm"),
    legend.key.height = unit(3, "cm")
  )

print(p_rank)

# ------------------------------------------------------------------------------
# 8) Figure 4B / S4: Hallmark GSEA per cancer + dotplot summary
# ------------------------------------------------------------------------------
hallmark_genes <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_term <- hallmark_genes[, c("gs_name", "gene_symbol")]

process_for_gsea <- function(df) {
  df$Gene <- sub("_.*", "", df$Gene_m6A_Pair)
  df$Correlation <- -df$Correlation
  df <- df[order(df$Correlation, decreasing = TRUE), , drop = FALSE]
  df
}

generate_gsea_input <- function(df) {
  df_mean <- aggregate(Correlation ~ Gene, data = df, FUN = mean)
  df_mean <- df_mean[order(df_mean$Correlation, decreasing = TRUE), , drop = FALSE]
  setNames(df_mean$Correlation, df_mean$Gene)
}

filtered_for_gsea <- lapply(filtered_results_list, process_for_gsea)
gsea_gene_lists <- lapply(filtered_for_gsea, generate_gsea_input)

cancer_names <- names(gsea_gene_lists)
all_gsea_results <- data.frame()

for (i in seq_along(gsea_gene_lists)) {
  cat("GSEA:", cancer_names[i], "\n")
  tryCatch({
    set.seed(123)
    gsea_res <- GSEA(
      geneList = gsea_gene_lists[[i]],
      TERM2GENE = hallmark_term,
      minGSSize = 10,
      maxGSSize = 5000,
      pvalueCutoff = 0.2,
      verbose = FALSE
    )
    res_df <- as.data.frame(gsea_res)
    res_df$cancer_type <- cancer_names[i]
    all_gsea_results <- rbind(all_gsea_results, res_df)
  }, error = function(e) {
    message("GSEA error for ", cancer_names[i], ": ", e$message)
  })
}

# Save/load as needed
# saveRDS(all_gsea_results, gsea_out_rds)
if (file.exists(gsea_out_rds)) all_gsea_results <- readRDS(gsea_out_rds)

# Dotplot-style heatmap (NES color, p.adjust as size)
cancer_levels <- unique(all_gsea_results$cancer_type)
cancer_levels <- c("PanCancer", setdiff(cancer_levels, "PanCancer"))
all_gsea_results$cancer_type <- factor(all_gsea_results$cancer_type, levels = cancer_levels)

all_gsea_results <- all_gsea_results %>%
  mutate(p_adjust_range = case_when(
    p.adjust < 0.01 ~ "0.01",
    p.adjust < 0.05 ~ "0.05",
    p.adjust < 0.1  ~ "0.1",
    TRUE            ~ "0.2"
  ))

pal_c2 <- c("#053061", "#134b87", "#327db7", "#6fafd2", "#c7e0ed",
            "#ffffff", "#fbd2bc", "#feab88", "#b71c2c", "#8b0824", "#6a0624")

all_gsea_results$p_adjust_range <- factor(all_gsea_results$p_adjust_range, levels = c("0.01", "0.05", "0.1", "0.2"))
x_levels <- cancer_levels
y_levels <- unique(all_gsea_results$ID)

all_gsea_results <- all_gsea_results %>%
  mutate(
    cancer_type = factor(cancer_type, levels = x_levels),
    ID = factor(ID, levels = y_levels)
  )

grid_df <- expand.grid(
  cancer_type = x_levels,
  ID = y_levels,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    cancer_type = factor(cancer_type, levels = x_levels),
    ID = factor(ID, levels = y_levels)
  )

lim <- max(abs(all_gsea_results$NES), na.rm = TRUE)
vals <- scales::rescale(seq(-lim, lim, length.out = length(pal_c2)), to = c(0, 1))

p_gsea_dot <- ggplot() +
  geom_tile(
    data = grid_df,
    aes(x = cancer_type, y = ID),
    width = 1, height = 1,
    fill = NA, color = "#D9D9D9", linewidth = 0.6
  ) +
  geom_point(
    data = all_gsea_results,
    aes(x = cancer_type, y = ID, size = p_adjust_range, color = NES)
  ) +
  scale_color_gradientn(
    colours = pal_c2,
    values = vals,
    limits = c(-lim, lim),
    name = "NES"
  ) +
  scale_size_manual(
    values = c("0.01" = 4, "0.05" = 3, "0.1" = 2, "0.2" = 1),
    name = expression(p.adjust)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(
    expand = c(0, 0),
    labels = function(x) gsub("_", " ", sub("^HALLMARK[_ ]?", "", x))
  ) +
  coord_cartesian(clip = "on") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Hallmark GSEA summary across cancer types",
    x = "Cancer Type",
    y = "Hallmark Pathway"
  )

print(p_gsea_dot)

# ------------------------------------------------------------------------------
# 9) Null shuffling (pan-cancer mean correlation) + empirical p + BH q + cutoffs
# ------------------------------------------------------------------------------
filtered_list_for_null <- readRDS(filtered_list_rds)

set.seed(123)
num_shuffles <- 1000L
num_pairs    <- nrow(filtered_list_for_null[[1L]])
num_cancers  <- length(filtered_list_for_null)

corr_mat <- do.call(cbind, lapply(filtered_list_for_null, function(df) as.numeric(df$Correlation)))
storage.mode(corr_mat) <- "double"

mean_corr_across_cancers <- matrix(NA_real_, nrow = num_pairs, ncol = num_shuffles)

for (s in seq_len(num_shuffles)) {
  idx_mat <- matrix(
    sample.int(num_pairs, size = num_pairs * num_cancers, replace = TRUE),
    nrow = num_pairs, ncol = num_cancers
  )
  col_idx <- rep(seq_len(num_cancers), each = num_pairs)
  picked <- corr_mat[cbind(as.vector(idx_mat), col_idx)]
  dim(picked) <- c(num_pairs, num_cancers)
  mean_corr_across_cancers[, s] <- rowMeans(picked)
  cat("Shuffle:", s, "\n")
}

null_vec <- as.vector(mean_corr_across_cancers)
# saveRDS(null_vec, shuffle_null_rds)

null_vec <- readRDS(shuffle_null_rds)
obs_df <- readRDS(mean_corr_rds)
rownames(obs_df) <- obs_df$Gene_m6A_Pair
obs <- as.numeric(obs_df$Correlation)
names(obs) <- rownames(obs_df)

null_sorted <- sort(null_vec)
N0 <- length(null_sorted)
k_i <- findInterval(obs, null_sorted, all.inside = TRUE)

p_left  <- (k_i + 1) / (N0 + 1)
p_right <- (N0 - k_i + 1) / (N0 + 1)

p_one <- ifelse(obs > 0, p_right,
                ifelse(obs < 0, p_left, pmin(p_left, p_right) * 2))

qval <- p.adjust(p_one, method = "BH")

alpha <- 0.20
sig_pos <- obs >= 0 & qval <= alpha
sig_neg <- obs <= 0 & qval <= alpha

cutoff_pos <- if (any(sig_pos)) min(obs[sig_pos]) else NA_real_
cutoff_neg <- if (any(sig_neg)) max(obs[sig_neg]) else NA_real_

cat("FDR <= 0.2 positive cutoff:", cutoff_pos, "\n")
cat("FDR <= 0.2 negative cutoff:", cutoff_neg, "\n")

out_tbl <- data.frame(
  pair = names(obs),
  mean_corr = obs,
  p_one = p_one,
  q = qval,
  direction = ifelse(obs > 0, "pos", ifelse(obs < 0, "neg", "zero")),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 10) Figure 4C: density plot with cutoffs
# ------------------------------------------------------------------------------
df_pc <- readRDS(filtered_list_pc_rds)[["PanCancer"]] %>%
  dplyr::select(Correlation)

neg_cut <- cutoff_neg
pos_cut <- cutoff_pos

dens <- density(df_pc$Correlation, adjust = 1)
densDF <- data.frame(x = dens$x, y = dens$y)

low_col  <- "#003070"
mid_col  <- "#E5E5E5"
high_col <- "#C22121"

p_density <- ggplot(densDF, aes(x, y)) +
  geom_area(aes(fill = x), alpha = 0.8, colour = NA) +
  scale_fill_gradient2(
    midpoint = 0,
    low = low_col,
    mid = mid_col,
    high = high_col,
    space = "Lab",
    guide = "none"
  ) +
  geom_area(data = subset(densDF, x <= neg_cut), aes(x, y), fill = low_col, alpha = 0.5) +
  geom_area(data = subset(densDF, x >= pos_cut), aes(x, y), fill = high_col, alpha = 0.5) +
  geom_line(colour = "grey20", linewidth = 0.8) +
  geom_vline(xintercept = neg_cut, linetype = "dashed", colour = low_col, linewidth = 0.6) +
  geom_vline(xintercept = pos_cut, linetype = "dashed", colour = high_col, linewidth = 0.6) +
  annotate("text", x = neg_cut, y = max(densDF$y) * 1.05,
           label = paste0("r = ", signif(neg_cut, 4)), hjust = 0,
           colour = low_col, size = 5, fontface = "bold") +
  annotate("text", x = pos_cut, y = max(densDF$y) * 1.05,
           label = paste0("r = ", signif(pos_cut, 4)), hjust = 1,
           colour = high_col, size = 5, fontface = "bold") +
  geom_rug(data = df_pc, mapping = aes(x = Correlation), inherit.aes = FALSE,
           sides = "b", alpha = 0.25, colour = "grey40") +
  labs(x = "Mean correlation (ρ)", y = "Density") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    axis.title = element_text(colour = "grey20"),
    axis.text  = element_text(colour = "grey30"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "grey70", fill = NA, linewidth = 0.5),
    plot.margin = margin(t = 8, r = 12, b = 8, l = 8)
  )

print(p_density)

# ------------------------------------------------------------------------------
# 11) Figure 4D: heatmap for negative/positive pair sets across cancer types
# ------------------------------------------------------------------------------
mean_correlations <- readRDS(mean_corr_rds)
colnames(mean_correlations) <- c("Pair", "Mean_Correlation", "Gene_name")

pos_pairs <- unique(mean_correlations$Pair[mean_correlations$Mean_Correlation >= cutoff_pos])
neg_pairs <- unique(mean_correlations$Pair[mean_correlations$Mean_Correlation <= cutoff_neg])

filtered_list_core <- readRDS(filtered_list_rds)
combined_data <- do.call(rbind, lapply(names(filtered_list_core), function(ct) {
  df <- filtered_list_core[[ct]]
  df$cancer_type <- ct
  df
}))

heatmap_data <- combined_data %>%
  filter(Gene_m6A_Pair %in% c(neg_pairs, pos_pairs))

heatmap_matrix <- reshape2::acast(heatmap_data, Gene_m6A_Pair ~ cancer_type, value.var = "Correlation")

neg_rows <- rownames(heatmap_matrix)[rownames(heatmap_matrix) %in% neg_pairs]
pos_rows <- rownames(heatmap_matrix)[rownames(heatmap_matrix) %in% pos_pairs]
heatmap_matrix <- heatmap_matrix[c(neg_rows, pos_rows), , drop = FALSE]

pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#3A53A4", "#5C59A7", "#8E7BB9", "#EFE4F1",
                             "#FFFFFF",
                             "#FCCCBC", "#F69175", "#EE3324"))(201),
  na_col = "grey",
  show_rownames = FALSE,
  show_colnames = TRUE,
  fontsize_col = 12
)

# ------------------------------------------------------------------------------
# 12) TF overlap enrichment test (Negative/Positive vs Overall)
# ------------------------------------------------------------------------------
target_go_ids <- c(
  "GO:0001227",
  "GO:0001228",
  "GO:0140297",
  "GO:0003714",
  "GO:0003713",
  "GO:0001217",
  "GO:0001216",
  "GO:0003700"
)

annotated_tf <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = target_go_ids,
  columns = c("SYMBOL", "GENENAME", "ENTREZID"),
  keytype = "GO"
)

TF_list <- unique(annotated_tf$SYMBOL)

positive_genes <- unique(mean_correlations$Gene_name[mean_correlations$Mean_Correlation >= cutoff_pos])
negative_genes <- unique(mean_correlations$Gene_name[mean_correlations$Mean_Correlation <= cutoff_neg])
overall_genes  <- unique(common_genes)

neg_TF <- length(intersect(negative_genes, TF_list))
pos_TF <- length(intersect(positive_genes, TF_list))
all_TF <- length(intersect(overall_genes,  TF_list))

n_neg <- length(negative_genes)
n_pos <- length(positive_genes)
n_all <- length(overall_genes)

counts <- c(Negative = neg_TF, Positive = pos_TF, Overall = all_TF)
totals <- c(Negative = n_neg,  Positive = n_pos,  Overall = n_all)

overall_test <- prop.test(x = counts, n = totals, correct = TRUE)
print(overall_test)

# ------------------------------------------------------------------------------
# 13) Helper utilities for 3-panel example plotting (RPPA/mRNA/m6A/residual)
# ------------------------------------------------------------------------------
split_gene_peak <- function(gene_peak) {
  parts <- strsplit(gene_peak, "_")[[1]]
  gene  <- parts[1]
  peak  <- paste(parts[-1], collapse = "_")
  list(gene = gene, peak = peak)
}

format_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 1e-4) return("P < 1e-4")
  paste0("P = ", signif(p, 3))
}

# NOTE: This function assumes you already built:
#   m6A_matrix_norm_list_selected, rppa_norm_list_selected, TCGA_mRNA_common_1_selected, observed_expected_list
# It will stop if any are missing.
plot_correlation_triplet <- function(gene, peak,
                                     m6A_matrix_norm_list_selected,
                                     rppa_norm_list_selected,
                                     TCGA_mRNA_common_1_selected,
                                     observed_expected_list) {
  stopifnot(gene %in% names(m6A_matrix_norm_list_selected))
  stopifnot(gene %in% names(rppa_norm_list_selected))
  stopifnot(gene %in% rownames(TCGA_mRNA_common_1_selected))
  stopifnot(gene %in% names(observed_expected_list))
  
  m6A_data  <- m6A_matrix_norm_list_selected[[gene]]
  rppa_data <- rppa_norm_list_selected[[gene]]
  mrna_data <- TCGA_mRNA_common_1_selected[gene, , drop = FALSE]
  obs_df    <- observed_expected_list[[gene]]
  
  stopifnot(peak %in% rownames(m6A_data))
  
  common_samples <- Reduce(intersect, list(colnames(m6A_data), colnames(rppa_data), colnames(mrna_data)))
  stopifnot(length(common_samples) >= 3)
  
  m6A_vec  <- as.numeric(m6A_data[peak, common_samples])
  rppa_vec <- as.numeric(rppa_data[1, common_samples])
  mrna_vec <- as.numeric(mrna_data[1, common_samples])
  
  obs_samples <- intersect(obs_df$SampleID, common_samples)
  stopifnot(length(obs_samples) >= 3)
  
  obs_vec <- obs_df$ObservedMinusExpected[match(obs_samples, obs_df$SampleID)]
  m6A_obs <- as.numeric(m6A_data[peak, obs_samples])
  
  ct1 <- suppressWarnings(cor.test(rppa_vec, mrna_vec, use = "complete.obs"))
  ct2 <- suppressWarnings(cor.test(mrna_vec, m6A_vec,  use = "complete.obs"))
  ct3 <- suppressWarnings(cor.test(m6A_obs, obs_vec,   use = "complete.obs"))
  
  mk_plot <- function(x, y, xlab, ylab, title, ct, col_point, col_smooth) {
    ggplot(data.frame(x, y), aes(x = x, y = y)) +
      geom_point(color = col_point, alpha = 0.6, size = 2) +
      geom_smooth(method = "lm",
                  color = col_smooth, fill = col_smooth,
                  alpha = 0.2, linewidth = 1) +
      annotate("text",
               x = min(x, na.rm = TRUE),
               y = max(y, na.rm = TRUE),
               label = paste0("cor = ", round(unname(ct$estimate), 3), "\n", format_p(ct$p.value)),
               hjust = 0, vjust = 1.2, size = 4) +
      labs(title = title, x = xlab, y = ylab) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "grey70", fill = NA),
        axis.title = element_text(face = "bold")
      )
  }
  
  p1 <- mk_plot(mrna_vec, rppa_vec, "mRNA expression", "Protein (RPPA)",
                paste0(gene, ": Protein vs mRNA"), ct1, "#4E669A", "#B8474D")
  p2 <- mk_plot(m6A_vec, mrna_vec, "m6A level", "mRNA expression",
                paste0(gene, " / ", peak, ": m6A vs mRNA"), ct2, "#2A9D8F", "#E76F51")
  p3 <- mk_plot(m6A_obs, obs_vec, "m6A level", "Observed - Expected",
                paste0(gene, " / ", peak, ": m6A vs residual"), ct3, "#264653", "#F4A261")
  
  p1m <- ggMarginal(p1, type = "histogram", margins = "both", size = 5, fill = "#4E669A", color = "white", alpha = 0.6)
  p2m <- ggMarginal(p2, type = "histogram", margins = "both", size = 5, fill = "#2A9D8F", color = "white", alpha = 0.6)
  p3m <- ggMarginal(p3, type = "histogram", margins = "both", size = 5, fill = "#264653", color = "white", alpha = 0.6)
  
  grid.arrange(p1m, p2m, p3m, ncol = 3)
}

# ==============================================================================
# End of script
# ==============================================================================
