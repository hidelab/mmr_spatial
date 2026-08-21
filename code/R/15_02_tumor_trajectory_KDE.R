# 15_02_tumor_trajectory_KDE.R
# Tumor density trajectory using 2D kernel density estimation (KDE).
# Estimates a smooth tumor density field per sample from tumor cell positions,
# then evaluates it at every cell's location. This gives a continuous, smooth
# gradient compared to the discrete k-NN approach in 15_01.
# Downstream analysis: same as 15_01 (binned heatmaps, GAM curves, WT vs DKO).

library(Seurat)
library(ggplot2)
library(patchwork)
library(MASS)
library(dplyr)
library(tidyr)
library(mgcv)
library(viridis)

# Source centralized parameters
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "15_tumor_trajectory_kde")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters (modifiable)
# =============================================================================
KDE_BANDWIDTH <- NULL
TARGET_CELLTYPE <- "Tumor"
N_BINS <- 10
TOP_N_GENES <- 50
GAM_K <- 4
CELLTYPES_TRAJECTORY <- IMMUNE_CELLTYPES_DE
KDE_GRID_N <- 200

# =============================================================================
# Load data
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

DefaultAssay(seu) <- "RNA"

# =============================================================================
# Compute KDE-based tumor density for ALL cells (per sample)
# =============================================================================
all_coords <- Embeddings(seu, "spatial")
celltypes <- seu$Tier1_celltype
samples <- seu$sample
is_tumor <- celltypes == TARGET_CELLTYPE

cat("Computing KDE tumor density per sample...\n")

tumor_density <- rep(NA_real_, ncol(seu))

for (s in levels(samples)) {
  idx <- which(samples == s)
  coords_s <- all_coords[idx, , drop = FALSE]
  tumor_idx_s <- which(is_tumor[idx])

  if (length(tumor_idx_s) < 10) {
    cat(sprintf("  %s: too few tumor cells (%d), skipping\n", s, length(tumor_idx_s)))
    next
  }

  tumor_coords <- coords_s[tumor_idx_s, , drop = FALSE]

  # Estimate bandwidth from tumor cells (or use provided)
  bw <- KDE_BANDWIDTH
  if (is.null(bw)) {
    bw_est <- bandwidth.nrd(tumor_coords[, 1])
    bw_est2 <- bandwidth.nrd(tumor_coords[, 2])
    bw <- c(bw_est, bw_est2)
  }

  # Build KDE on a grid from tumor cell positions
  x_range <- range(coords_s[, 1])
  y_range <- range(coords_s[, 2])

  kde_fit <- kde2d(tumor_coords[, 1], tumor_coords[, 2],
                   h = bw,
                   n = KDE_GRID_N,
                   lims = c(x_range, y_range))

  # Interpolate KDE density at each cell's position
  # Use bilinear interpolation from the grid
  x_grid <- kde_fit$x
  y_grid <- kde_fit$y
  z_grid <- kde_fit$z

  cell_density <- numeric(length(idx))
  for (i in seq_along(idx)) {
    cx <- coords_s[i, 1]
    cy <- coords_s[i, 2]

    # Find nearest grid indices
    xi <- findInterval(cx, x_grid, all.inside = TRUE)
    yi <- findInterval(cy, y_grid, all.inside = TRUE)

    # Bilinear interpolation
    x1 <- x_grid[xi]; x2 <- x_grid[xi + 1]
    y1 <- y_grid[yi]; y2 <- y_grid[yi + 1]

    fQ11 <- z_grid[xi, yi]
    fQ21 <- z_grid[xi + 1, yi]
    fQ12 <- z_grid[xi, yi + 1]
    fQ22 <- z_grid[xi + 1, yi + 1]

    denom <- (x2 - x1) * (y2 - y1)
    cell_density[i] <- (fQ11 * (x2 - cx) * (y2 - cy) +
                         fQ21 * (cx - x1) * (y2 - cy) +
                         fQ12 * (x2 - cx) * (cy - y1) +
                         fQ22 * (cx - x1) * (cy - y1)) / denom
  }

  # Normalize to [0, 1] within sample
  d_min <- min(cell_density)
  d_max <- max(cell_density)
  if (d_max > d_min) {
    cell_density <- (cell_density - d_min) / (d_max - d_min)
  } else {
    cell_density <- rep(0.5, length(cell_density))
  }

  tumor_density[idx] <- cell_density
  cat(sprintf("  %s: %d cells, %d tumor, bw = (%.1f, %.1f)\n",
              s, length(idx), length(tumor_idx_s), bw[1], bw[2]))
}

seu$tumor_density_kde <- tumor_density
cat("KDE tumor density computed for all cells.\n")

# =============================================================================
# Helper: bin cells into density deciles and compute mean expression
# =============================================================================
compute_binned_expression <- function(expr_mat, density_vec, n_bins = N_BINS) {
  bin_breaks <- quantile(density_vec, probs = seq(0, 1, length.out = n_bins + 1),
                         na.rm = TRUE)
  bin_breaks[1] <- bin_breaks[1] - 1e-6
  if (length(unique(bin_breaks)) < length(bin_breaks)) {
    bin_breaks <- seq(min(density_vec, na.rm = TRUE) - 1e-6,
                      max(density_vec, na.rm = TRUE),
                      length.out = n_bins + 1)
  }
  bins <- cut(density_vec, breaks = bin_breaks, labels = seq_len(n_bins))
  bin_midpoints <- (bin_breaks[-1] + bin_breaks[-(n_bins + 1)]) / 2

  binned_expr <- matrix(NA, nrow = nrow(expr_mat), ncol = n_bins)
  rownames(binned_expr) <- rownames(expr_mat)
  colnames(binned_expr) <- paste0("bin_", seq_len(n_bins))

  for (b in seq_len(n_bins)) {
    cells_in_bin <- which(bins == b)
    if (length(cells_in_bin) > 1) {
      binned_expr[, b] <- rowMeans(expr_mat[, cells_in_bin, drop = FALSE])
    } else if (length(cells_in_bin) == 1) {
      binned_expr[, b] <- expr_mat[, cells_in_bin]
    }
  }

  list(binned_expr = binned_expr, midpoints = bin_midpoints, bins = bins)
}

# =============================================================================
# Helper: select top variable genes along the trajectory
# =============================================================================
select_trajectory_genes <- function(binned_expr, top_n = TOP_N_GENES) {
  row_vars <- apply(binned_expr, 1, var, na.rm = TRUE)
  row_means <- rowMeans(binned_expr, na.rm = TRUE)
  keep <- row_means > 0.05 & !is.na(row_vars)
  row_vars_filt <- row_vars[keep]
  top_genes <- names(sort(row_vars_filt, decreasing = TRUE))[seq_len(min(top_n, length(row_vars_filt)))]
  top_genes
}

# =============================================================================
# Analysis per cell type: WT, DKO, and combined
# =============================================================================

for (ct in CELLTYPES_TRAJECTORY) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
  ct_dir <- file.path(OUTPUT_DIR, ct_safe)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n=== %s ===\n", ct))

  ct_idx <- which(celltypes == ct & !is.na(tumor_density))
  if (length(ct_idx) < 50) {
    cat(sprintf("  Skipping %s: too few cells (%d)\n", ct, length(ct_idx)))
    next
  }

  ct_density <- tumor_density[ct_idx]
  ct_condition <- seu$condition[ct_idx]
  ct_sample <- seu$sample[ct_idx]

  expr_mat <- GetAssayData(seu, assay = "RNA", layer = "data")[, ct_idx]

  cat(sprintf("  Total cells: %d (WT: %d, DKO: %d)\n",
              length(ct_idx),
              sum(ct_condition == REFERENCE_CONDITION),
              sum(ct_condition == TEST_CONDITION)))

  # ---------------------------------------------------------------------------
  # Combined: binned heatmap
  # ---------------------------------------------------------------------------
  binned_all <- compute_binned_expression(expr_mat, ct_density)
  top_genes <- select_trajectory_genes(binned_all$binned_expr)

  if (length(top_genes) < 5) {
    cat(sprintf("  Skipping %s: too few variable genes\n", ct))
    next
  }

  heatmap_mat <- binned_all$binned_expr[top_genes, ]
  heatmap_scaled <- t(scale(t(heatmap_mat)))

  hm_df <- as.data.frame(heatmap_scaled) %>%
    mutate(gene = rownames(heatmap_scaled)) %>%
    pivot_longer(-gene, names_to = "bin", values_to = "zscore") %>%
    mutate(bin_num = as.integer(gsub("bin_", "", bin)),
           gene = factor(gene, levels = rev(top_genes)))

  p_heatmap <- ggplot(hm_df, aes(x = bin_num, y = gene, fill = zscore)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno", name = "Z-score") +
    scale_x_continuous(breaks = c(1, N_BINS), labels = c("Low", "High")) +
    labs(x = "Tumor Density (KDE)", y = NULL,
         title = sprintf("%s — Gene Expression along Tumor Density (KDE)", ct)) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 7))

  ggsave(file.path(ct_dir, "heatmap_trajectory_combined.pdf"),
         p_heatmap, width = 6, height = max(4, length(top_genes) * 0.25))

  # ---------------------------------------------------------------------------
  # Per-condition: separate heatmaps
  # ---------------------------------------------------------------------------
  for (cond in levels(ct_condition)) {
    cond_idx <- which(ct_condition == cond)
    if (length(cond_idx) < 30) next

    binned_cond <- compute_binned_expression(expr_mat[, cond_idx, drop = FALSE],
                                             ct_density[cond_idx])

    hm_cond <- binned_cond$binned_expr[top_genes, , drop = FALSE]
    hm_cond_scaled <- t(scale(t(hm_cond)))

    hm_cond_df <- as.data.frame(hm_cond_scaled) %>%
      mutate(gene = rownames(hm_cond_scaled)) %>%
      pivot_longer(-gene, names_to = "bin", values_to = "zscore") %>%
      mutate(bin_num = as.integer(gsub("bin_", "", bin)),
             gene = factor(gene, levels = rev(top_genes)))

    p_cond <- ggplot(hm_cond_df, aes(x = bin_num, y = gene, fill = zscore)) +
      geom_tile() +
      scale_fill_viridis_c(option = "inferno", name = "Z-score") +
      scale_x_continuous(breaks = c(1, N_BINS), labels = c("Low", "High")) +
      labs(x = "Tumor Density (KDE)", y = NULL,
           title = sprintf("%s — %s (KDE)", ct, cond)) +
      theme_publication() +
      theme(axis.text.y = element_text(size = 7))

    ggsave(file.path(ct_dir, sprintf("heatmap_trajectory_%s.pdf", cond)),
           p_cond, width = 6, height = max(4, length(top_genes) * 0.25))
  }

  # ---------------------------------------------------------------------------
  # GAM fits: top genes, WT vs DKO smooth curves
  # ---------------------------------------------------------------------------
  gam_dir <- file.path(ct_dir, "gam_curves")
  dir.create(gam_dir, recursive = TRUE, showWarnings = FALSE)

  gam_genes <- top_genes
  gam_results <- list()

  for (gene in gam_genes) {
    gene_expr <- as.numeric(expr_mat[gene, ])
    df_gam <- data.frame(
      expr = gene_expr,
      density = ct_density,
      condition = ct_condition,
      sample = ct_sample
    )

    fit_gam <- tryCatch(
      gam(expr ~ s(density, by = condition, k = GAM_K) + condition,
          data = df_gam, method = "REML"),
      error = function(e) NULL
    )

    if (is.null(fit_gam)) next

    newdata_pred <- expand.grid(
      density = seq(min(ct_density), max(ct_density), length.out = 100),
      condition = levels(ct_condition)
    )
    pred <- predict(fit_gam, newdata = newdata_pred, se.fit = TRUE)
    newdata_pred$fit <- pred$fit
    newdata_pred$se <- pred$se.fit
    newdata_pred$ci_lo <- pred$fit - 1.96 * pred$se.fit
    newdata_pred$ci_hi <- pred$fit + 1.96 * pred$se.fit

    fit_null <- tryCatch(
      gam(expr ~ s(density, k = GAM_K), data = df_gam, method = "REML"),
      error = function(e) NULL
    )
    anova_p <- NA
    if (!is.null(fit_null)) {
      aov_test <- anova(fit_null, fit_gam, test = "Chisq")
      anova_p <- aov_test[["Pr(>Chi)"]][2]
    }

    gam_results[[gene]] <- data.frame(gene = gene, p_condition = anova_p)

    p_gam <- ggplot(newdata_pred, aes(x = density, y = fit, color = condition, fill = condition)) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = CONDITION_COLORS) +
      scale_fill_manual(values = CONDITION_COLORS) +
      labs(x = "Tumor Density (KDE)", y = "Expression",
           title = sprintf("%s — %s", gene, ct),
           subtitle = ifelse(!is.na(anova_p),
                             sprintf("Condition effect p = %.3g", anova_p), "")) +
      theme_publication()

    ggsave(file.path(gam_dir, sprintf("gam_%s.pdf", gene)),
           p_gam, width = 5, height = 4)
  }

  if (length(gam_results) > 0) {
    gam_summary <- bind_rows(gam_results)
    gam_summary <- gam_summary[order(gam_summary$p_condition), ]
    write.csv(gam_summary, file.path(ct_dir, "gam_condition_test.csv"), row.names = FALSE)
    cat(sprintf("  GAM results saved for %d genes\n", nrow(gam_summary)))
  }

  # ---------------------------------------------------------------------------
  # Density distribution: where along the gradient do cells of this type sit?
  # ---------------------------------------------------------------------------
  df_dist <- data.frame(density = ct_density, condition = ct_condition)

  p_dist <- ggplot(df_dist, aes(x = density, fill = condition, color = condition)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(x = "Tumor Density (KDE)", y = "Density of cells",
         title = sprintf("%s — Position along Tumor Gradient (KDE)", ct)) +
    theme_publication()

  ggsave(file.path(ct_dir, "density_distribution_along_gradient.pdf"),
         p_dist, width = 5, height = 4)

  cat(sprintf("  Done: %s\n", ct))
}

# =============================================================================
# Summary: spatial map of KDE tumor density (all cells)
# =============================================================================
sample_ids <- levels(seu$sample)
max_x_span <- max(tapply(all_coords[, 1], samples, function(x) diff(range(x))))
max_y_span <- max(tapply(all_coords[, 2], samples, function(x) diff(range(x))))
half_span <- max(max_x_span, max_y_span) / 2 * 1.05

plots_density <- lapply(sample_ids, function(s) {
  idx <- which(samples == s)
  df_plot <- data.frame(x = all_coords[idx, 1], y = all_coords[idx, 2],
                        density = tumor_density[idx])
  ggplot(df_plot, aes(x = x, y = y, color = density)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_viridis_c(option = "magma", limits = c(0, 1), name = "Tumor\nDensity") +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    ggtitle(s) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
})

p_spatial_density <- wrap_plots(plots_density, nrow = 1) +
  plot_annotation(title = "Tumor Density Field (KDE)",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_tumor_density_field_kde.pdf"),
       p_spatial_density, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_density_field_kde.png"),
       p_spatial_density, width = 4 * length(sample_ids), height = 5, dpi = 300)

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")