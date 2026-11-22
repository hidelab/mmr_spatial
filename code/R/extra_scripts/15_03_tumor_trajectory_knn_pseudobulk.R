# 15_01_tumor_trajectory_knn_pseudobulk.R
# Tumor density trajectory with pseudobulk GAM fits.
# Same k-NN tumor density as 15_01, but gene expression is aggregated into
# pseudobulk (per sample x density bin) before GAM fitting. This accounts for
# sample-level correlation and avoids the overpowered single-cell approach.

library(Seurat)
library(ggplot2)
library(patchwork)
library(FNN)
library(dplyr)
library(tidyr)
library(mgcv)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# Source centralized parameters
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "15_tumor_trajectory_knn_pseudobulk")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters (modifiable)
# =============================================================================
K_NEIGHBORS <- 30
TARGET_CELLTYPE <- "Tumor"
N_BINS <- 10
TOP_N_GENES <- 100
MIN_CELLS_PER_BIN <- 5
GAM_K <- 4
CELLTYPES_TRAJECTORY <- IMMUNE_CELLTYPES_DE

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
# Compute tumor density score for ALL cells (per sample k-NN)
# =============================================================================
all_coords <- Embeddings(seu, "spatial")
celltypes <- seu$Tier1_celltype
samples <- seu$sample
is_tumor <- as.integer(celltypes == TARGET_CELLTYPE)

cat(sprintf("Computing %d-NN tumor density for all cells...\n", K_NEIGHBORS))

tumor_density <- rep(NA_real_, ncol(seu))

for (s in levels(samples)) {
  idx <- which(samples == s)
  if (length(idx) < K_NEIGHBORS + 1) next

  coords_s <- all_coords[idx, , drop = FALSE]
  knn_result <- get.knn(coords_s, k = K_NEIGHBORS)
  is_tumor_s <- is_tumor[idx]

  neighbor_is_tumor <- matrix(is_tumor_s[knn_result$nn.index],
                              nrow = nrow(knn_result$nn.index),
                              ncol = K_NEIGHBORS)
  tumor_density[idx] <- rowMeans(neighbor_is_tumor)
  cat(sprintf("  %s: %d cells\n", s, length(idx)))
}

seu$tumor_density <- tumor_density
cat("Tumor density computed for all cells.\n")

# =============================================================================
# Helper: create pseudobulk per sample x density bin
# =============================================================================
create_pseudobulk <- function(expr_mat, density_vec, sample_vec, condition_vec,
                              n_bins = N_BINS, min_cells = MIN_CELLS_PER_BIN) {
  # Bin density into quantiles (across all cells for consistent bins)
  bin_breaks <- quantile(density_vec, probs = seq(0, 1, length.out = n_bins + 1),
                         na.rm = TRUE)
  bin_breaks[1] <- bin_breaks[1] - 1e-6
  if (length(unique(bin_breaks)) < length(bin_breaks)) {
    bin_breaks <- seq(min(density_vec, na.rm = TRUE) - 1e-6,
                      max(density_vec, na.rm = TRUE),
                      length.out = n_bins + 1)
  }
  bins <- as.integer(cut(density_vec, breaks = bin_breaks, labels = seq_len(n_bins)))
  bin_midpoints <- (bin_breaks[-1] + bin_breaks[-(n_bins + 1)]) / 2

  # Aggregate per sample x bin
  groups <- paste(sample_vec, bins, sep = "__")
  unique_groups <- unique(groups)

  pb_list <- list()
  pb_meta <- list()

  for (g in unique_groups) {
    cell_idx <- which(groups == g)
    if (length(cell_idx) < min_cells) next

    parts <- strsplit(g, "__")[[1]]
    samp <- parts[1]
    bin_id <- as.integer(parts[2])
    if (is.na(bin_id)) next

    pb_list[[g]] <- rowMeans(expr_mat[, cell_idx, drop = FALSE])
    pb_meta[[g]] <- data.frame(
      group = g,
      sample = samp,
      condition = as.character(condition_vec[cell_idx[1]]),
      bin = bin_id,
      density_midpoint = bin_midpoints[bin_id],
      n_cells = length(cell_idx),
      stringsAsFactors = FALSE
    )
  }

  if (length(pb_list) == 0) return(NULL)

  pb_mat <- do.call(cbind, pb_list)
  pb_meta_df <- bind_rows(pb_meta)
  pb_meta_df$condition <- factor(pb_meta_df$condition,
                                 levels = c(REFERENCE_CONDITION, TEST_CONDITION))

  list(expr = pb_mat, meta = pb_meta_df, bin_midpoints = bin_midpoints)
}

# =============================================================================
# Helper: select top variable genes from pseudobulk
# =============================================================================
select_trajectory_genes_pb <- function(pb_expr, pb_meta, top_n = TOP_N_GENES) {
  # Compute variance of bin-level means (collapse samples within bins)
  bins <- sort(unique(pb_meta$bin))
  bin_means <- matrix(NA, nrow = nrow(pb_expr), ncol = length(bins))
  rownames(bin_means) <- rownames(pb_expr)

  for (i in seq_along(bins)) {
    idx <- which(pb_meta$bin == bins[i])
    if (length(idx) > 1) {
      bin_means[, i] <- rowMeans(pb_expr[, idx, drop = FALSE])
    } else if (length(idx) == 1) {
      bin_means[, i] <- pb_expr[, idx]
    }
  }

  row_vars <- apply(bin_means, 1, var, na.rm = TRUE)
  row_means <- rowMeans(bin_means, na.rm = TRUE)
  keep <- row_means > 0.05 & !is.na(row_vars)
  row_vars_filt <- row_vars[keep]
  top_genes <- names(sort(row_vars_filt, decreasing = TRUE))[seq_len(min(top_n, length(row_vars_filt)))]
  top_genes
}

# =============================================================================
# Analysis per cell type
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
  # Create pseudobulk (sample x density bin)
  # ---------------------------------------------------------------------------
  pb <- create_pseudobulk(expr_mat, ct_density, ct_sample, ct_condition)

  if (is.null(pb) || nrow(pb$meta) < 10) {
    cat(sprintf("  Skipping %s: too few pseudobulk groups (%s)\n",
                ct, ifelse(is.null(pb), "0", nrow(pb$meta))))
    next
  }

  cat(sprintf("  Pseudobulk groups: %d (samples x bins with >= %d cells)\n",
              nrow(pb$meta), MIN_CELLS_PER_BIN))

  # ---------------------------------------------------------------------------
  # Select top variable genes
  # ---------------------------------------------------------------------------
  top_genes <- select_trajectory_genes_pb(pb$expr, pb$meta)

  if (length(top_genes) < 5) {
    cat(sprintf("  Skipping %s: too few variable genes\n", ct))
    next
  }

  # ---------------------------------------------------------------------------
  # Heatmap: bin-averaged expression (combined)
  # ---------------------------------------------------------------------------
  bins <- sort(unique(pb$meta$bin))
  bin_means_all <- matrix(NA, nrow = length(top_genes), ncol = length(bins))
  rownames(bin_means_all) <- top_genes
  colnames(bin_means_all) <- paste0("bin_", bins)

  for (i in seq_along(bins)) {
    idx <- which(pb$meta$bin == bins[i])
    if (length(idx) > 1) {
      bin_means_all[, i] <- rowMeans(pb$expr[top_genes, idx, drop = FALSE])
    } else if (length(idx) == 1) {
      bin_means_all[, i] <- pb$expr[top_genes, idx]
    }
  }

  heatmap_scaled <- t(scale(t(bin_means_all)))

  hm_df <- as.data.frame(heatmap_scaled) %>%
    mutate(gene = rownames(heatmap_scaled)) %>%
    pivot_longer(-gene, names_to = "bin", values_to = "zscore") %>%
    mutate(bin_num = as.integer(gsub("bin_", "", bin)),
           gene = factor(gene, levels = rev(top_genes)))

  p_heatmap <- ggplot(hm_df, aes(x = bin_num, y = gene, fill = zscore)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno", name = "Z-score") +
    scale_x_continuous(breaks = c(1, N_BINS), labels = c("Low", "High")) +
    labs(x = "Tumor Density", y = NULL,
         title = sprintf("%s — Pseudobulk Trajectory (k-NN, k=%d)", ct, K_NEIGHBORS)) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 7))

  ggsave(file.path(ct_dir, "heatmap_trajectory_combined.pdf"),
         p_heatmap, width = 6, height = max(4, length(top_genes) * 0.25))

  # ---------------------------------------------------------------------------
  # Per-condition heatmaps
  # ---------------------------------------------------------------------------
  for (cond in levels(pb$meta$condition)) {
    cond_idx <- which(pb$meta$condition == cond)
    if (length(cond_idx) < 5) next

    bins_cond <- sort(unique(pb$meta$bin[cond_idx]))
    bin_means_cond <- matrix(NA, nrow = length(top_genes), ncol = length(bins_cond))
    rownames(bin_means_cond) <- top_genes
    colnames(bin_means_cond) <- paste0("bin_", bins_cond)

    for (i in seq_along(bins_cond)) {
      idx_b <- which(pb$meta$bin == bins_cond[i] & pb$meta$condition == cond)
      if (length(idx_b) > 1) {
        bin_means_cond[, i] <- rowMeans(pb$expr[top_genes, idx_b, drop = FALSE])
      } else if (length(idx_b) == 1) {
        bin_means_cond[, i] <- pb$expr[top_genes, idx_b]
      }
    }

    hm_cond_scaled <- t(scale(t(bin_means_cond)))

    hm_cond_df <- as.data.frame(hm_cond_scaled) %>%
      mutate(gene = rownames(hm_cond_scaled)) %>%
      pivot_longer(-gene, names_to = "bin", values_to = "zscore") %>%
      mutate(bin_num = as.integer(gsub("bin_", "", bin)),
             gene = factor(gene, levels = rev(top_genes)))

    p_cond <- ggplot(hm_cond_df, aes(x = bin_num, y = gene, fill = zscore)) +
      geom_tile() +
      scale_fill_viridis_c(option = "inferno", name = "Z-score") +
      scale_x_continuous(breaks = range(bins_cond), labels = c("Low", "High")) +
      labs(x = "Tumor Density", y = NULL,
           title = sprintf("%s — %s (Pseudobulk)", ct, cond)) +
      theme_publication() +
      theme(axis.text.y = element_text(size = 7))

    ggsave(file.path(ct_dir, sprintf("heatmap_trajectory_%s.pdf", cond)),
           p_cond, width = 6, height = max(4, length(top_genes) * 0.25))
  }

  # ---------------------------------------------------------------------------
  # GAM fits on pseudobulk: expression ~ density, weighted by n_cells
  # ---------------------------------------------------------------------------
  gam_dir <- file.path(ct_dir, "gam_curves")
  dir.create(gam_dir, recursive = TRUE, showWarnings = FALSE)

  gam_genes <- top_genes
  gam_results <- list()

  for (gene in gam_genes) {
    df_gam <- data.frame(
      expr = as.numeric(pb$expr[gene, ]),
      density = pb$meta$density_midpoint,
      condition = pb$meta$condition,
      sample = pb$meta$sample,
      weight = pb$meta$n_cells
    )

    # Weighted GAM with condition-specific smooths
    fit_gam <- tryCatch(
      gam(expr ~ s(density, by = condition, k = GAM_K) + condition,
          data = df_gam, weights = weight, method = "REML"),
      error = function(e) NULL
    )

    if (is.null(fit_gam)) next

    # Predict smooth curves
    newdata_pred <- expand.grid(
      density = seq(min(df_gam$density), max(df_gam$density), length.out = 100),
      condition = levels(df_gam$condition)
    )
    pred <- predict(fit_gam, newdata = newdata_pred, se.fit = TRUE)
    newdata_pred$fit <- pred$fit
    newdata_pred$se <- pred$se.fit
    newdata_pred$ci_lo <- pred$fit - 1.96 * pred$se.fit
    newdata_pred$ci_hi <- pred$fit + 1.96 * pred$se.fit

    # Test condition effect
    fit_null <- tryCatch(
      gam(expr ~ s(density, k = GAM_K), data = df_gam, weights = weight, method = "REML"),
      error = function(e) NULL
    )
    anova_p <- NA
    if (!is.null(fit_null)) {
      aov_test <- anova(fit_null, fit_gam, test = "Chisq")
      anova_p <- aov_test[["Pr(>Chi)"]][2]
    }

    gam_results[[gene]] <- data.frame(gene = gene, p_condition = anova_p)

    # Plot with pseudobulk points overlaid
    p_gam <- ggplot(newdata_pred, aes(x = density, y = fit, color = condition, fill = condition)) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      geom_point(data = df_gam, aes(x = density, y = expr, color = condition, size = weight),
                 alpha = 0.5, inherit.aes = FALSE) +
      scale_color_manual(values = CONDITION_COLORS) +
      scale_fill_manual(values = CONDITION_COLORS) +
      scale_size_continuous(range = c(1, 4), guide = "none") +
      labs(x = "Tumor Density (k-NN)", y = "Mean Expression (pseudobulk)",
           title = sprintf("%s — %s", gene, ct),
           subtitle = ifelse(!is.na(anova_p),
                             sprintf("Condition effect p = %.3g", anova_p), "")) +
      theme_publication()

    ggsave(file.path(gam_dir, sprintf("gam_%s.pdf", gene)),
           p_gam, width = 5, height = 4)
  }

  # Save GAM summary
  if (length(gam_results) > 0) {
    gam_summary <- bind_rows(gam_results)
    gam_summary <- gam_summary[order(gam_summary$p_condition), ]
    write.csv(gam_summary, file.path(ct_dir, "gam_condition_test.csv"), row.names = FALSE)
    cat(sprintf("  GAM results saved for %d genes\n", nrow(gam_summary)))
  }

  # ---------------------------------------------------------------------------
  # ComplexHeatmap: WT and DKO side by side, rows sorted by condition p-value
  # ---------------------------------------------------------------------------
  if (length(gam_results) > 0) {
    gam_summary_valid <- gam_summary[!is.na(gam_summary$p_condition), ]
    hm_genes <- gam_summary_valid$gene
    hm_genes <- hm_genes[hm_genes %in% rownames(pb$expr)]

    if (length(hm_genes) >= 5) {
      bins_sorted <- sort(unique(pb$meta$bin))

      # Build per-condition bin-averaged matrices
      build_bin_matrix <- function(cond) {
        mat <- matrix(NA, nrow = length(hm_genes), ncol = length(bins_sorted))
        rownames(mat) <- hm_genes
        colnames(mat) <- bins_sorted
        for (i in seq_along(bins_sorted)) {
          idx_b <- which(pb$meta$bin == bins_sorted[i] & pb$meta$condition == cond)
          if (length(idx_b) > 1) {
            mat[, i] <- rowMeans(pb$expr[hm_genes, idx_b, drop = FALSE])
          } else if (length(idx_b) == 1) {
            mat[, i] <- pb$expr[hm_genes, idx_b]
          }
        }
        mat
      }

      mat_wt <- build_bin_matrix(REFERENCE_CONDITION)
      mat_dko <- build_bin_matrix(TEST_CONDITION)

      # Row z-score across BOTH conditions together (shared scale)
      mat_combined <- cbind(mat_wt, mat_dko)
      row_means <- rowMeans(mat_combined, na.rm = TRUE)
      row_sds <- apply(mat_combined, 1, sd, na.rm = TRUE)
      row_sds[row_sds == 0] <- 1

      mat_wt_z <- (mat_wt - row_means) / row_sds
      mat_dko_z <- (mat_dko - row_means) / row_sds

      # P-value annotation as color
      pvals <- gam_summary_valid$p_condition[match(hm_genes, gam_summary_valid$gene)]
      neglog10p <- -log10(pmax(pvals, 1e-300))

      pval_col_fun <- colorRamp2(c(0, -log10(0.05), max(neglog10p)),
                                  c("#E0E0E0", "#FFCDD2", "#B71C1C"))

      row_ha <- rowAnnotation(
        `p-value` = neglog10p,
        col = list(`p-value` = pval_col_fun),
        annotation_name_rot = 90,
        width = unit(0.4, "cm"),
        show_legend = TRUE,
        annotation_legend_param = list(`p-value` = list(title = "-log10(p)"))
      )

      # Same inferno palette as the ggplot heatmaps
      col_fun <- colorRamp2(c(-2, -1, 0, 1, 2),
                            viridis::inferno(5))

      ht_wt <- Heatmap(mat_wt_z,
                        cluster_rows = FALSE,
                        cluster_columns = FALSE,
                        col = col_fun,
                        na_col = "#F5F5F5",
                        name = "Z-score",
                        column_title = REFERENCE_CONDITION,
                        column_names_rot = 45,
                        show_row_names = TRUE,
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                        width = unit(3, "cm"),
                        border = TRUE,
                        show_heatmap_legend = TRUE)

      ht_dko <- Heatmap(mat_dko_z,
                         cluster_rows = FALSE,
                         cluster_columns = FALSE,
                         col = col_fun,
                         na_col = "#F5F5F5",
                         name = " ",
                         column_title = TEST_CONDITION,
                         column_names_rot = 45,
                         show_row_names = FALSE,
                         width = unit(3, "cm"),
                         border = TRUE,
                         show_heatmap_legend = FALSE)

      ht_combined <- ht_wt + ht_dko + row_ha

      pdf(file.path(ct_dir, "heatmap_trajectory_sidebyside.pdf"),
          width = 7, height = max(5, length(hm_genes) * 0.2 + 2))
      draw(ht_combined,
           column_title = sprintf("%s — Trajectory (sorted by condition p-value)", ct),
           column_title_gp = gpar(fontsize = 11, fontface = "bold"))
      dev.off()

      cat(sprintf("  Side-by-side heatmap saved (%d genes)\n", length(hm_genes)))
    }
  }

  # ---------------------------------------------------------------------------
  # Density distribution: where along the gradient do cells sit?
  # ---------------------------------------------------------------------------
  df_dist <- data.frame(density = ct_density, condition = ct_condition)

  p_dist <- ggplot(df_dist, aes(x = density, fill = condition, color = condition)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(x = "Tumor Density (k-NN)", y = "Density of cells",
         title = sprintf("%s — Position along Tumor Gradient", ct)) +
    theme_publication()

  ggsave(file.path(ct_dir, "density_distribution_along_gradient.pdf"),
         p_dist, width = 5, height = 4)

  cat(sprintf("  Done: %s\n", ct))
}

# =============================================================================
# Summary: spatial map of tumor density (all cells)
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
  plot_annotation(title = sprintf("Tumor Density Field (k = %d)", K_NEIGHBORS),
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_tumor_density_field.pdf"),
       p_spatial_density, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_density_field.png"),
       p_spatial_density, width = 4 * length(sample_ids), height = 5, dpi = 300)

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")