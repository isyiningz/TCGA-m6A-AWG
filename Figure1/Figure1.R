################################################################################
# Setup
################################################################################
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# ---- Input: cancer type colors (CSV must contain: cancer_type, cancer_type_color)
TCGA_Cancer_Type_color <- read.csv(
  "/Users/yiningzhao/OneDrive - Baylor College of Medicine/M6A/TCGA_Cancer_Type_color.csv",
  stringsAsFactors = FALSE
)

################################################################################
# Figure 1A: Sample counts by cancer type (colored bars)
################################################################################

# Order shown in the plot (top-to-bottom after coord_flip)
cancer_order <- c(
  "GBM","HNSC","DLBC","THYM","THCA","BRCA","LUAD","LUSC","STAD",
  "LIHC","ACC","PCPG","CHOL","KICH","KIRC","KIRP","COAD","UCEC",
  "CESC","BLCA","TGCT","SARC","SKCM"
)

# Counts aligned to cancer_order
counts <- data.frame(
  cancer_type = cancer_order,
  count = c(
    2,2,15,8,16,38,13,4,4,
    10,4,3,1,5,2,11,3,27,
    6,13,5,7,27
  )
)

plot_df <- counts %>%
  left_join(
    TCGA_Cancer_Type_color %>% dplyr::select(cancer_type, cancer_type_color),
    by = "cancer_type"
  ) %>%
  mutate(cancer_type = factor(cancer_type, levels = rev(cancer_order)))

ref_lines <- c(10, 20, 30)

p_fig1a <- ggplot(plot_df, aes(x = cancer_type, y = count, fill = cancer_type)) +
  geom_hline(
    yintercept = ref_lines,
    linetype = "dashed",
    color = "grey70",
    linewidth = 0.4
  ) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(
    values = setNames(plot_df$cancer_type_color, as.character(plot_df$cancer_type))
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Number of samples", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  )

print(p_fig1a)

################################################################################
# Technical replicate / null distribution test (fixed 8 same-cancer-type pairs)
################################################################################
# Assumptions (must already exist in your environment):
#   - m: matrix [peaks x samples], colnames(m) are sample names
#   - m6A_sample_barcodes: data.frame with columns sample_16, Cancer_type_abbv
#   - score_pair(vec1, vec2): function returning a numeric score (e.g., correlation)
#   - obs_median: observed median score across the true replicate pairs
#   - USE_CORR: TRUE => right-tail test (larger is more extreme), FALSE => left-tail
#   - plot_null_vs_obs_gg(): your plotting helper (optional)

# ---- Map samples to cancer types
sample_names <- colnames(m)
orig_ids     <- substr(sample_names, 1, 16)

ct_df <- m6A_sample_barcodes
stopifnot(all(c("sample_16", "Cancer_type_abbv") %in% names(ct_df)))

mt <- match(orig_ids, ct_df$sample_16)
cancer_type_vec <- ct_df$Cancer_type_abbv[mt]

idx_by_ct <- split(seq_along(sample_names), cancer_type_vec)
idx_by_ct <- idx_by_ct[!is.na(names(idx_by_ct))]

# ---- Controls
N_PAIRS_PER_ITER  <- 8L
B                 <- 1000L
EXCLUDE_SAME_ORIG <- TRUE
EXCLUDE_REPL      <- FALSE              # set TRUE to drop replicate columns from pools
REPL_PATTERN      <- "A933-53$"         # used only when EXCLUDE_REPL = TRUE
NO_REUSE_SAMPLES  <- TRUE               # prefer disjoint pairs per iteration when possible

# ---- Single-iteration sampler: draw 8 pairs within the same cancer type and return median(score)
make_sameCT_pairs_median_fixed <- function(n_pairs = N_PAIRS_PER_ITER) {
  
  # Build per-cancer pools; optionally exclude specific replicate columns
  pools <- lapply(idx_by_ct, function(v) {
    if (EXCLUDE_REPL) v <- v[!grepl(REPL_PATTERN, sample_names[v])]
    v
  })
  
  # Keep cancer types with at least 2 samples
  pools <- pools[sapply(pools, length) >= 2]
  if (!length(pools)) return(NA_real_)
  
  # Check if strict no-reuse is feasible
  noreuse <- NO_REUSE_SAMPLES
  if (noreuse) {
    cap <- sum(floor(sapply(pools, length) / 2))
    if (cap < n_pairs) noreuse <- FALSE
  }
  
  pairs <- matrix(NA_integer_, nrow = 0, ncol = 2)
  avail <- pools
  
  trials     <- 0L
  max_trials <- 5000L
  
  while (nrow(pairs) < n_pairs && trials < max_trials) {
    
    elig_cts <- names(avail)[sapply(avail, length) >= 2]
    if (!length(elig_cts)) break
    
    # Weight by number of possible pairs in each cancer type
    ws <- sapply(avail[elig_cts], function(v) choose(length(v), 2))
    ws <- ws / sum(ws)
    
    ct <- sample(elig_cts, 1, prob = ws)
    v  <- avail[[ct]]
    
    pr <- sample(v, 2, replace = FALSE)
    
    # Optional: forbid pairing two columns from the same 16-char original ID
    if (EXCLUDE_SAME_ORIG && orig_ids[pr[1]] == orig_ids[pr[2]]) {
      trials <- trials + 1L
      next
    }
    
    pairs <- rbind(pairs, pr)
    
    if (noreuse) {
      avail[[ct]] <- setdiff(v, pr)
    }
    
    trials <- trials + 1L
  }
  
  # If not enough pairs, relax no-reuse to finish the round
  if (nrow(pairs) < n_pairs) {
    need   <- n_pairs - nrow(pairs)
    tries2 <- 0L
    
    while (need > 0L && tries2 < 5000L) {
      elig_cts <- names(pools)[sapply(pools, length) >= 2]
      if (!length(elig_cts)) break
      
      ws <- sapply(pools[elig_cts], function(v) choose(length(v), 2))
      ws <- ws / sum(ws)
      
      ct <- sample(elig_cts, 1, prob = ws)
      v  <- pools[[ct]]
      
      pr <- sample(v, 2, replace = FALSE)
      
      if (EXCLUDE_SAME_ORIG && orig_ids[pr[1]] == orig_ids[pr[2]]) {
        tries2 <- tries2 + 1L
        next
      }
      
      pairs <- rbind(pairs, pr)
      need  <- n_pairs - nrow(pairs)
      tries2 <- tries2 + 1L
    }
  }
  
  if (nrow(pairs) < n_pairs) return(NA_real_)
  
  sc <- apply(pairs, 1, function(v) score_pair(m[, v[1]], m[, v[2]]))
  median(sc, na.rm = TRUE)
}

# ---- Generate null distribution
set.seed(2025)
null_medians <- replicate(B, make_sameCT_pairs_median_fixed(N_PAIRS_PER_ITER))
null_medians <- null_medians[is.finite(null_medians)]

if (!length(null_medians)) {
  stop("Failed to build a null distribution. Check sample counts per cancer type.")
}

# ---- Empirical p-value vs observed median
p_empirical <- if (USE_CORR) {
  (sum(null_medians >= obs_median) + 1) / (length(null_medians) + 1)
} else {
  (sum(null_medians <= obs_median) + 1) / (length(null_medians) + 1)
}

cat(sprintf(
  "Same-CT fixed-%d pairs: mean=%.4f, sd=%.4f, empirical p=%.4g, n_valid=%d\n",
  N_PAIRS_PER_ITER,
  mean(null_medians),
  sd(null_medians),
  p_empirical,
  length(null_medians)
))

# ---- Optional plot (requires your helper)
if (exists("plot_null_vs_obs_gg")) {
  xlab_txt <- if (USE_CORR) {
    sprintf("Median correlation across %d pairs", N_PAIRS_PER_ITER)
  } else {
    sprintf("Median |Δ| across %d pairs", N_PAIRS_PER_ITER)
  }
  
  p_null <- plot_null_vs_obs_gg(
    null_vals  = null_medians,
    obs_val    = obs_median,
    p          = p_empirical,
    title      = if (USE_CORR) {
      sprintf("Null (same cancer type, %d pairs): median(correlation)", N_PAIRS_PER_ITER)
    } else {
      sprintf("Null (same cancer type, %d pairs): median(|Δ|)", N_PAIRS_PER_ITER)
    },
    right_tail = USE_CORR,
    xlab       = xlab_txt
  )
  print(p_null)
}

# ---- Optional: append to summary table if present
if (exists("summary_df")) {
  summary_df2 <- rbind(
    summary_df,
    data.frame(
      metric          = if (USE_CORR) "correlation" else "median_abs_diff",
      n_true_pairs    = if (exists("obs_scores")) length(obs_scores) else NA_integer_,
      observed_median = obs_median,
      nullA_mean      = mean(null_medians),
      nullA_sd        = sd(null_medians),
      nullA_p         = p_empirical,
      nullB_mean      = NA_real_,
      nullB_sd        = NA_real_,
      nullB_p         = NA_real_,
      row.names = "Free_sameCancerType_fixed8"
    )
  )
  print(summary_df2)
}

################################################################################
# Z-score + normal-approx p-value from the null distribution
################################################################################

z_p_from_null <- function(obs_val, null_vals, right_tail = TRUE, two_sided = FALSE) {
  null_vals <- null_vals[is.finite(null_vals)]
  mu  <- mean(null_vals)
  sdv <- sd(null_vals)
  
  if (!is.finite(sdv) || sdv == 0) {
    return(list(
      z = NA_real_, p = NA_real_, mean = mu, sd = sdv, n = length(null_vals),
      note = "Null distribution SD is 0 or not finite; cannot compute z-score."
    ))
  }
  
  z <- (obs_val - mu) / sdv
  
  if (two_sided) {
    p <- 2 * min(pnorm(z), 1 - pnorm(z))
  } else {
    p <- if (right_tail) 1 - pnorm(z) else pnorm(z)
  }
  
  list(z = z, p = p, mean = mu, sd = sdv, n = length(null_vals))
}

# One-sided (direction matches USE_CORR)
res_z_one_sided <- z_p_from_null(
  obs_val    = obs_median,
  null_vals  = null_medians,
  right_tail = USE_CORR,
  two_sided  = FALSE
)

cat(sprintf(
  "Normal approx (one-sided): mean=%.6f, sd=%.6f, z=%.3f, p=%.4g, n=%d\n",
  res_z_one_sided$mean, res_z_one_sided$sd,
  res_z_one_sided$z, res_z_one_sided$p,
  res_z_one_sided$n
))

# Two-sided
res_z_two_sided <- z_p_from_null(
  obs_val    = obs_median,
  null_vals  = null_medians,
  right_tail = USE_CORR,
  two_sided  = TRUE
)

cat(sprintf("Normal approx (two-sided): p=%.3e\n", res_z_two_sided$p))

# Echo empirical p-value if computed
if (exists("p_empirical")) {
  cat(sprintf("Empirical p (one-sided) = %.4g\n", p_empirical))
}
