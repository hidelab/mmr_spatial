# 14_02_moranI.R
# Compute global and local Moran's I on the binary "is_tumor" variable.
# Global Moran's I: one value per sample — high positive = tumor cells cluster
# together (WT pattern); near zero = random mixing (DKO pattern).
# Local Moran's I (LISA): per-cell statistic identifying cells in significant
# tumor clusters vs. mixed zones.

library(Seurat)
library(ggplot2)
library(patchwork)
library(FNN)
library(spdep)
library(dplyr)

# Source centralized parameters
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "14_moranI")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters
# =============================================================================
K_NEIGHBORS <- 20
N_PERMUTATIONS <- 999
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

all_coords <- Embeddings(seu, "spatial")
celltypes <- seu$Tier1_celltype
samples <- seu$sample
is_tumor <- as.numeric(celltypes == TARGET_CELLTYPE)

sample_ids <- levels(samples)

# =============================================================================
# Global Moran's I per sample
# =============================================================================
cat(sprintf("Computing global Moran's I (k = %d, %d permutations)...\n",
            K_NEIGHBORS, N_PERMUTATIONS))

global_results <- list()

for (s in sample_ids) {
  idx <- which(samples == s)
  coords_s <- all_coords[idx, , drop = FALSE]
  is_tumor_s <- is_tumor[idx]

  # Build k-NN spatial weights
  knn_result <- get.knn(coords_s, k = K_NEIGHBORS)
  nb <- lapply(seq_len(nrow(knn_result$nn.index)), function(i) {
    as.integer(knn_result$nn.index[i, ])
  })
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_along(nb))
  lw <- nb2listw(nb, style = "W")

  # Global Moran's I with permutation test
  moran_test <- moran.mc(is_tumor_s, lw, nsim = N_PERMUTATIONS, alternative = "greater")

  global_results[[s]] <- data.frame(
    sample = s,
    condition = as.character(seu$condition[idx[1]]),
    n_cells = length(idx),
    n_tumor = sum(is_tumor_s),
    tumor_frac = mean(is_tumor_s),
    moran_I = moran_test$statistic,
    moran_expected = -1 / (length(idx) - 1),
    moran_p = moran_test$p.value,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %s: I = %.4f, p = %.4f (n = %d, tumor_frac = %.3f)\n",
              s, moran_test$statistic, moran_test$p.value,
              length(idx), mean(is_tumor_s)))
}

df_global <- bind_rows(global_results)
df_global$condition <- factor(df_global$condition,
                              levels = c(REFERENCE_CONDITION, TEST_CONDITION))

write.csv(df_global, file.path(OUTPUT_DIR, "global_moranI_per_sample.csv"),
          row.names = FALSE)

cat("\nGlobal Moran's I summary:\n")
print(df_global[, c("sample", "condition", "moran_I", "moran_p", "tumor_frac")])

# =============================================================================
# Plot: Global Moran's I by condition
# =============================================================================
p_global <- ggplot(df_global, aes(x = condition, y = moran_I, color = condition)) +
  geom_jitter(width = 0.1, size = 3.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = "black", linewidth = 0.5) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "Global Moran's I",
       title = "Spatial Autocorrelation of Tumor Identity",
       subtitle = sprintf("k = %d neighbors, %d permutations", K_NEIGHBORS, N_PERMUTATIONS)) +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "dotplot_global_moranI.pdf"),
       p_global, width = 3.5, height = 5)

# =============================================================================
# Local Moran's I (LISA) per sample
# =============================================================================
cat("\nComputing local Moran's I (LISA)...\n")

seu$lisa_I <- NA_real_
seu$lisa_cluster <- NA_character_

for (s in sample_ids) {
  idx <- which(samples == s)
  coords_s <- all_coords[idx, , drop = FALSE]
  is_tumor_s <- is_tumor[idx]

  knn_result <- get.knn(coords_s, k = K_NEIGHBORS)
  nb <- lapply(seq_len(nrow(knn_result$nn.index)), function(i) {
    as.integer(knn_result$nn.index[i, ])
  })
  class(nb) <- "nb"
  attr(nb, "region.id") <- as.character(seq_along(nb))
  lw <- nb2listw(nb, style = "W")

  # Local Moran's I
  lisa <- localmoran(is_tumor_s, lw)
  seu$lisa_I[idx] <- lisa[, "Ii"]

  # Classify LISA quadrants
  z_x <- scale(is_tumor_s)[, 1]
  lag_x <- lag.listw(lw, is_tumor_s)
  z_lag <- scale(lag_x)[, 1]
  p_lisa <- lisa[, "Pr(z != E(Ii))"]

  quad <- rep("Not significant", length(idx))
  sig <- p_lisa < 0.05
  quad[sig & z_x > 0 & z_lag > 0] <- "High-High (Tumor cluster)"
  quad[sig & z_x < 0 & z_lag < 0] <- "Low-Low (Non-tumor cluster)"
  quad[sig & z_x > 0 & z_lag < 0] <- "High-Low (Tumor outlier)"
  quad[sig & z_x < 0 & z_lag > 0] <- "Low-High (Non-tumor outlier)"

  seu$lisa_cluster[idx] <- quad
  cat(sprintf("  %s: HH = %d, LL = %d, HL = %d, LH = %d, NS = %d\n",
              s,
              sum(quad == "High-High (Tumor cluster)"),
              sum(quad == "Low-Low (Non-tumor cluster)"),
              sum(quad == "High-Low (Tumor outlier)"),
              sum(quad == "Low-High (Non-tumor outlier)"),
              sum(quad == "Not significant")))
}

# =============================================================================
# LISA summary: fraction of tumor cells in HH clusters by condition
# =============================================================================
tumor_cells <- which(celltypes == TARGET_CELLTYPE)
df_lisa_tumor <- data.frame(
  sample = samples[tumor_cells],
  condition = seu$condition[tumor_cells],
  lisa_cluster = seu$lisa_cluster[tumor_cells],
  stringsAsFactors = FALSE
)

lisa_summary <- df_lisa_tumor %>%
  group_by(sample, condition) %>%
  summarise(
    n_tumor = n(),
    n_HH = sum(lisa_cluster == "High-High (Tumor cluster)"),
    frac_HH = n_HH / n_tumor,
    .groups = "drop"
  )

cat("\nFraction of tumor cells in HH (tumor cluster) per sample:\n")
print(as.data.frame(lisa_summary))

write.csv(lisa_summary, file.path(OUTPUT_DIR, "lisa_tumor_HH_fraction.csv"),
          row.names = FALSE)

# =============================================================================
# Plot: LISA HH fraction by condition
# =============================================================================
lisa_summary$condition <- factor(lisa_summary$condition,
                                 levels = c(REFERENCE_CONDITION, TEST_CONDITION))

p_lisa_frac <- ggplot(lisa_summary, aes(x = condition, y = frac_HH, color = condition)) +
  geom_jitter(width = 0.1, size = 3.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = "black", linewidth = 0.5) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "Fraction of Tumor Cells in HH Cluster",
       title = "Local Moran's I: Tumor Cells in Spatial Clusters",
       subtitle = "High-High = tumor cell surrounded by tumor") +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "dotplot_lisa_HH_fraction.pdf"),
       p_lisa_frac, width = 3.5, height = 5)

# =============================================================================
# Plot: Spatial maps of LISA clusters
# =============================================================================
LISA_COLORS <- c(
  "High-High (Tumor cluster)"     = "#B71C1C",
  "Low-Low (Non-tumor cluster)"   = "#1565C0",
  "High-Low (Tumor outlier)"      = "#F57C00",
  "Low-High (Non-tumor outlier)"  = "#7B1FA2",
  "Not significant"               = "#E0E0E0"
)

max_x_span <- max(tapply(all_coords[, 1], samples, function(x) diff(range(x))))
max_y_span <- max(tapply(all_coords[, 2], samples, function(x) diff(range(x))))
half_span <- max(max_x_span, max_y_span) / 2 * 1.05

plots_lisa <- lapply(sample_ids, function(s) {
  idx <- which(samples == s)
  df_plot <- data.frame(
    x = all_coords[idx, 1],
    y = all_coords[idx, 2],
    cluster = factor(seu$lisa_cluster[idx], levels = names(LISA_COLORS))
  )
  # Plot non-significant first, significant on top
  df_plot <- df_plot[order(df_plot$cluster == "Not significant", decreasing = TRUE), ]

  ggplot(df_plot, aes(x = x, y = y, color = cluster)) +
    geom_point(size = 0.2, alpha = 0.6) +
    scale_color_manual(values = LISA_COLORS, drop = FALSE) +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    ggtitle(s) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10)) +
    NoLegend()
})

p_lisa_spatial <- wrap_plots(plots_lisa, nrow = 1) +
  plot_annotation(
    title = "Local Moran's I — LISA Clusters",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_lisa_clusters.pdf"),
       p_lisa_spatial, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_lisa_clusters.png"),
       p_lisa_spatial, width = 4 * length(sample_ids), height = 5, dpi = 300)

# Shared legend
p_legend <- ggplot(data.frame(cluster = names(LISA_COLORS)),
                   aes(x = 1, y = seq_along(cluster), color = cluster)) +
  geom_point(size = 4) +
  scale_color_manual(values = LISA_COLORS, name = "LISA Cluster") +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "right")
ggsave(file.path(OUTPUT_DIR, "spatial_lisa_legend.pdf"),
       p_legend, width = 4, height = 3)

cat("\nDone. Outputs saved to:", OUTPUT_DIR, "\n")
