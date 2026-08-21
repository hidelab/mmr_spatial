# 14_01_knn_admixture.R
# Quantify tumor spatial isolation vs. admixture using k-NN purity scores.
# For each tumor cell, computes the fraction of its k spatial nearest neighbors
# that are also tumor cells. High purity = isolated tumor (WT pattern);
# low purity = admixed with other cell types (DKO pattern).

library(Seurat)
library(ggplot2)
library(patchwork)
library(FNN)
library(dplyr)
library(tidyr)
library(glmmTMB)

# Source centralized parameters
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "14_knn_admixture")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters
# =============================================================================
K_NEIGHBORS <- 30
TARGET_CELLTYPE <- "Tumor"

# =============================================================================
# Load data
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# =============================================================================
# Compute k-NN tumor purity per sample
# =============================================================================
# We build spatial k-NN graphs per sample to avoid cross-sample neighbors.

all_coords <- Embeddings(seu, "spatial")
celltypes <- seu$Tier1_celltype
samples <- seu$sample
is_tumor <- as.integer(celltypes == TARGET_CELLTYPE)

cat(sprintf("Computing %d-NN purity for %s cells...\n", K_NEIGHBORS, TARGET_CELLTYPE))

purity_scores <- rep(NA_real_, ncol(seu))

for (s in levels(samples)) {
  idx <- which(samples == s)
  if (length(idx) < K_NEIGHBORS + 1) next

  coords_s <- all_coords[idx, , drop = FALSE]
  knn_result <- get.knn(coords_s, k = K_NEIGHBORS)

  is_tumor_s <- is_tumor[idx]

  # For each cell, fraction of k neighbors that are tumor
  neighbor_is_tumor <- matrix(is_tumor_s[knn_result$nn.index],
                              nrow = nrow(knn_result$nn.index),
                              ncol = K_NEIGHBORS)
  purity_scores[idx] <- rowMeans(neighbor_is_tumor)

  n_tumor <- sum(is_tumor_s)
  cat(sprintf("  %s: %d cells (%d tumor)\n", s, length(idx), n_tumor))
}

seu$knn_tumor_purity <- purity_scores

# Subset to tumor cells only for downstream analysis
tumor_idx <- which(celltypes == TARGET_CELLTYPE)
df_tumor <- data.frame(
  cell = colnames(seu)[tumor_idx],
  sample = samples[tumor_idx],
  condition = seu$condition[tumor_idx],
  purity = purity_scores[tumor_idx],
  stringsAsFactors = FALSE
)
df_tumor <- df_tumor[!is.na(df_tumor$purity), ]

cat(sprintf("\nTumor cells with purity scores: %d\n", nrow(df_tumor)))

# =============================================================================
# Summary statistics per sample
# =============================================================================
sample_summary <- df_tumor %>%
  group_by(sample, condition) %>%
  summarise(
    n_cells = n(),
    mean_purity = mean(purity),
    median_purity = median(purity),
    sd_purity = sd(purity),
    q25_purity = quantile(purity, 0.25),
    q75_purity = quantile(purity, 0.75),
    .groups = "drop"
  )

cat("\nPer-sample summary:\n")
print(as.data.frame(sample_summary))

write.csv(sample_summary, file.path(OUTPUT_DIR, "tumor_purity_sample_summary.csv"),
          row.names = FALSE)

# =============================================================================
# Statistical test: Beta GLMM (logit link) — appropriate for bounded [0,1] response
# =============================================================================
# Beta regression requires open interval (0, 1). Squeeze exact 0/1 values inward.
eps <- 1e-4
df_tumor$purity_beta <- pmin(pmax(df_tumor$purity, eps), 1 - eps)

cat("\nFitting beta GLMM: purity ~ condition + (1 | sample) ...\n")
fit_beta <- glmmTMB(purity_beta ~ condition + (1 | sample),
                    data = df_tumor,
                    family = beta_family(link = "logit"))

beta_summary <- summary(fit_beta)
cat("\n--- Beta GLMM Summary ---\n")
print(beta_summary)

# Extract condition effect on logit scale
beta_coefs <- fixef(fit_beta)$cond
beta_se <- sqrt(diag(vcov(fit_beta)$cond))
beta_z <- beta_coefs / beta_se
beta_pval <- 2 * pnorm(abs(beta_z), lower.tail = FALSE)

cat(sprintf("\nCondition effect (%s vs %s):\n", TEST_CONDITION, REFERENCE_CONDITION))
cat(sprintf("  Logit-scale estimate: %.4f (SE = %.4f)\n", beta_coefs[2], beta_se[2]))
cat(sprintf("  z = %.3f, p = %.2e\n", beta_z[2], beta_pval[2]))

# Marginal means on response scale
mean_wt_logit <- beta_coefs[1]
mean_dko_logit <- beta_coefs[1] + beta_coefs[2]
mean_wt_prob <- plogis(mean_wt_logit)
mean_dko_prob <- plogis(mean_dko_logit)
cat(sprintf("  Predicted mean purity: WT = %.3f, DKO = %.3f\n", mean_wt_prob, mean_dko_prob))

# Save model summary
sink(file.path(OUTPUT_DIR, "beta_glmm_summary.txt"))
cat("Beta GLMM: purity ~ condition + (1 | sample)\n")
cat("Family: beta (logit link)\n")
cat(sprintf("N cells = %d, N samples = %d\n\n", nrow(df_tumor), length(unique(df_tumor$sample))))
print(beta_summary)
cat("\n\nMarginal predicted means (response scale):\n")
cat(sprintf("  WT:  %.4f\n", mean_wt_prob))
cat(sprintf("  DKO: %.4f\n", mean_dko_prob))
sink()

# Also run Wilcoxon on sample medians as a simple non-parametric check
wt_medians <- sample_summary$median_purity[sample_summary$condition == REFERENCE_CONDITION]
dko_medians <- sample_summary$median_purity[sample_summary$condition == TEST_CONDITION]
wilcox_result <- wilcox.test(wt_medians, dko_medians, alternative = "two.sided", exact = FALSE)
cat(sprintf("\nWilcoxon on sample medians (supplementary): W = %.1f, p = %.4f\n",
            wilcox_result$statistic, wilcox_result$p.value))

# =============================================================================
# Plot 1: Violin/box plot of purity by condition
# =============================================================================
p_violin <- ggplot(df_tumor, aes(x = condition, y = purity, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "k-NN Tumor Purity",
       title = sprintf("Tumor Cell Neighborhood Purity (k = %d)", K_NEIGHBORS),
       subtitle = sprintf("Beta GLMM p = %.2e", beta_pval[2])) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "violin_tumor_purity_by_condition.pdf"),
       p_violin, width = 4, height = 5)

# =============================================================================
# Plot 2: Per-sample box plot
# =============================================================================
p_sample <- ggplot(df_tumor, aes(x = sample, y = purity, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "k-NN Tumor Purity",
       title = sprintf("Tumor Purity per Sample (k = %d)", K_NEIGHBORS)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "violin_tumor_purity_by_sample.pdf"),
       p_sample, width = 7, height = 5)

# =============================================================================
# Plot 3: Sample-level median dot plot (the "n=7" summary view)
# =============================================================================
p_dots <- ggplot(sample_summary, aes(x = condition, y = median_purity,
                                      color = condition)) +
  geom_jitter(width = 0.1, size = 3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = "black", linewidth = 0.5) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "Median Tumor Purity (per sample)",
       title = "Sample-Level Tumor Purity",
       subtitle = sprintf("Wilcoxon p = %.3f", wilcox_result$p.value)) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "dotplot_sample_median_purity.pdf"),
       p_dots, width = 3.5, height = 5)

# =============================================================================
# Plot 4: Spatial scatter colored by purity score (tumor cells only)
# =============================================================================
sample_ids <- levels(seu$sample)

# Max span for uniform scaling
tumor_spatial <- all_coords[tumor_idx, ]
max_x_span <- max(tapply(tumor_spatial[, 1], samples[tumor_idx], function(x) diff(range(x))))
max_y_span <- max(tapply(tumor_spatial[, 2], samples[tumor_idx], function(x) diff(range(x))))
half_span <- max(max_x_span, max_y_span) / 2 * 1.05

df_spatial <- data.frame(
  x = all_coords[tumor_idx, 1],
  y = all_coords[tumor_idx, 2],
  sample = samples[tumor_idx],
  purity = purity_scores[tumor_idx]
)
df_spatial <- df_spatial[!is.na(df_spatial$purity), ]

plots_spatial <- lapply(sample_ids, function(s) {
  sub <- df_spatial[df_spatial$sample == s, ]
  if (nrow(sub) == 0) return(NULL)
  ggplot(sub, aes(x = x, y = y, color = purity)) +
    geom_point(size = 0.3, alpha = 0.5) +
    scale_color_viridis_c(option = "viridis", limits = c(0, 1), name = "Purity") +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    ggtitle(s) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
})
plots_spatial <- plots_spatial[!sapply(plots_spatial, is.null)]

p_spatial <- wrap_plots(plots_spatial, nrow = 1) +
  plot_annotation(
    title = sprintf("Tumor k-NN Purity (k = %d)", K_NEIGHBORS),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_tumor_purity.pdf"),
       p_spatial, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_purity.png"),
       p_spatial, width = 4 * length(sample_ids), height = 5, dpi = 300)

# =============================================================================
# Plot 5: Density plot overlay (WT vs DKO distributions)
# =============================================================================
p_density <- ggplot(df_tumor, aes(x = purity, fill = condition, color = condition)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = "k-NN Tumor Purity", y = "Density",
       title = "Distribution of Tumor Neighborhood Purity") +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "density_tumor_purity_by_condition.pdf"),
       p_density, width = 6, height = 4)

# =============================================================================
# Save full per-cell results
# =============================================================================
write.csv(df_tumor, file.path(OUTPUT_DIR, "tumor_cell_purity_scores.csv"),
          row.names = FALSE)

cat("\nDone. Outputs saved to:", OUTPUT_DIR, "\n")
