# 15_05_spatial_gradient_masks.R
# Spatial visualization of the tumor-proximity trajectory genes (from 15_04).
#
# Motivation: 15_04 found genes whose CD8+ T expression trajectory *reshapes* in
# DKO — flat/low in WT but surging toward the tumor core in DKO. This script shows
# that in real tissue. Each gene's figure has one column per sample (WT samples
# first, then DKO) and three rows:
#   1. tumor cells only, highlighted in a single color (tissue context)
#   2. CD8+ T cells only, colored by expression
#   3. overlay: tumor cells underneath, CD8+ T cells on top
# The CD8+ T color scale is shared across all samples, so WT-vs-DKO is directly
# comparable. Violin plots of the same genes in CD8+ T cells (one violin per
# sample, filled by condition) go to the violin_plots/ subfolder.
#
# Genes are read from the 15_04 primary-spec (pseudobulk_nb_fs) condition-test CSV
# and filtered to those significant on EITHER the shape OR the level test
# (BH-FDR < GENE_FDR), so the plotted set tracks the reported hits automatically.
#
# Adapted from 08_01_spatial_masks_08lambda.R.

# --- Thread limits (shared server) ------------------------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS)

library(Seurat)
library(ggplot2)
library(patchwork)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Parameters (modifiable)
# =============================================================================
# Cell type(s) to color by expression. The 15_04 trajectory was on CD8+ T; add
# "Treg" / "NK" here to include more T/lymphoid populations.
MASK_CELLTYPE  <- c("CD8+ T")
TUMOR_CELLTYPE <- "Tumor"
TUMOR_COLOR    <- TIER1_COLORS[["Tumor"]]   # tumor cells (context only, not expression)
# CD8+ T expression gradient: light blue (low) -> dark blue (high)
CD8_GRADIENT   <- c("#9ECAE1", "#08306B")

# Genes to plot: read from the 15_04 condition test and keep those significant on
# EITHER the shape OR the level test (BH-FDR < GENE_FDR). This tracks the primary
# spec's hit list directly instead of a hardcoded set — re-run 15_04 and this
# updates automatically.
CONDITION_TEST_CSV <- file.path(PROJECT_DIR, "output", "R", "15_tumor_trajectory_working",
                                "CD8pos_T", "gam_condition_test_pseudobulk_nb.csv")
GENE_FDR <- 0.1        # keep genes with fdr_shape < GENE_FDR OR fdr_level < GENE_FDR

ct_res <- read.csv(CONDITION_TEST_CSV, stringsAsFactors = FALSE)
sig_gene <- (!is.na(ct_res$fdr_shape) & ct_res$fdr_shape < GENE_FDR) |
            (!is.na(ct_res$fdr_level) & ct_res$fdr_level < GENE_FDR)
GENES <- ct_res$gene[sig_gene]
cat(sprintf("Read %d genes from %s (FDR < %.2f on shape OR level)\n",
            length(GENES), basename(CONDITION_TEST_CSV), GENE_FDR))
if (length(GENES) == 0) stop("No genes passed the FDR cutoff; nothing to plot.")

EXPR_CLIP_Q  <- 0.99   # cap color scale at this quantile (per gene) to tame outliers
PT_SIZE_MASK <- 0.8    # masked (T) cell point size
PT_SIZE_BG   <- 0.2    # tumor cell point size

# Selected-samples figure (PNG only, in selected/): one row per sample,
# columns = tumor / CD8+ T expression / overlay
SELECTED_SAMPLES <- c("WT1", "DKO4")

# Outputs go in a subfolder named after the mask cell type(s), e.g. CD8pos_T/
mask_dir_name <- paste(gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", MASK_CELLTYPE)), collapse = "_")
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "15_spatial_gradient_masks", mask_dir_name)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
SELECTED_DIR <- file.path(OUTPUT_DIR, "selected")
dir.create(SELECTED_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Load merged-domain Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$sample    <- factor(seu$sample,
                        levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))
DefaultAssay(seu) <- "RNA"

# QC filter intentionally OFF (reviewed): plot all cells, matching 15_04.

# =============================================================================
# Pull coordinates, cell-type masks, and expression
# =============================================================================
coords     <- Embeddings(seu, "spatial")
samples    <- seu$sample
celltypes  <- seu$Tier1_celltype
is_mask    <- celltypes %in% MASK_CELLTYPE
is_tumor   <- celltypes == TUMOR_CELLTYPE

cat(sprintf("Masked cell type(s): %s -> %d cells\n",
            paste(MASK_CELLTYPE, collapse = ", "), sum(is_mask)))

sample_ids <- levels(samples)[levels(samples) %in% unique(as.character(samples))]

# Uniform, centered coordinate frame across samples (per-sample centroid removed)
max_x_span <- max(tapply(coords[, 1], samples, function(x) diff(range(x))))
max_y_span <- max(tapply(coords[, 2], samples, function(x) diff(range(x))))
half_span  <- max(max_x_span, max_y_span) / 2 * 1.05

available_genes <- intersect(GENES, rownames(seu))
missing_genes   <- setdiff(GENES, rownames(seu))
if (length(missing_genes) > 0)
  cat("Genes not found (skipping):", paste(missing_genes, collapse = ", "), "\n")
cat(sprintf("Plotting %d genes across %d samples\n",
            length(available_genes), length(sample_ids)))

expr_layer <- GetAssayData(seu, assay = "RNA", layer = "data")

# Center a sample's coords about its own midrange so coord_fixed(-half, half) fits.
center_coords <- function(idx) {
  x <- coords[idx, 1]; y <- coords[idx, 2]
  cbind(x = x - mean(range(x)), y = y - mean(range(y)))
}

mask_label  <- paste(MASK_CELLTYPE, collapse = " + ")
panel_frame <- list(
  coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)),
  theme_void(),
  theme(plot.title = element_text(hjust = 0.5, size = 10), legend.position = "right")
)
row_label <- function(label) {
  ggplot() + annotate("text", x = 0, y = 0, label = label, angle = 90,
                      fontface = "bold", size = 4.5) + theme_void()
}

# =============================================================================
# Per gene: rows = tumor only / CD8+ T only / overlay; one column per sample
# =============================================================================
for (gene in available_genes) {
  cat(sprintf("  %s ...\n", gene))
  gvec <- as.numeric(expr_layer[gene, ])

  # Shared color limits across all samples (masked cells only) for comparability
  mask_vals <- gvec[is_mask]
  hi <- as.numeric(quantile(mask_vals[mask_vals > 0], EXPR_CLIP_Q, na.rm = TRUE))
  if (!is.finite(hi) || hi <= 0) hi <- max(mask_vals, na.rm = TRUE)
  if (!is.finite(hi) || hi <= 0) hi <- 1
  mask_scale <- scale_color_gradient(low = CD8_GRADIENT[1], high = CD8_GRADIENT[2],
                                     name = paste0(mask_label, "\n", gene),
                                     limits = c(0, hi), oob = scales::squish)

  panels <- lapply(sample_ids, function(s) {
    s_idx  <- which(samples == s)
    cc     <- center_coords(s_idx)
    m_tum  <- is_tumor[s_idx]
    m_msk  <- is_mask[s_idx]
    df_tum <- data.frame(x = cc[m_tum, 1], y = cc[m_tum, 2])
    df_msk <- data.frame(x = cc[m_msk, 1], y = cc[m_msk, 2], expr = pmin(gvec[s_idx][m_msk], hi))
    df_msk <- df_msk[order(df_msk$expr), ]   # draw high-expressers on top
    title  <- sprintf("%s (%s)", s, unique(as.character(seu$condition[s_idx])))

    tumor_layer <- geom_point(data = df_tum, aes(x, y), color = TUMOR_COLOR,
                              size = PT_SIZE_BG, alpha = 0.5)
    mask_layer  <- geom_point(data = df_msk, aes(x, y, color = expr), size = PT_SIZE_MASK, alpha = 0.9)

    list(
      tumor   = ggplot() + tumor_layer + ggtitle(title) + panel_frame,
      mask    = ggplot() + mask_layer + mask_scale + panel_frame,
      overlay = ggplot() + tumor_layer + mask_layer + mask_scale + panel_frame
    )
  })
  names(panels) <- sample_ids

  rows <- c(tumor = paste(TUMOR_CELLTYPE, "cells"), mask = paste(mask_label, "cells"),
            overlay = "Overlay")
  grid <- unlist(lapply(names(rows), function(r) {
    c(list(row_label(rows[[r]])), lapply(panels, `[[`, r))
  }), recursive = FALSE)

  p_all <- wrap_plots(grid, ncol = length(sample_ids) + 1,
                      widths = c(0.15, rep(1, length(sample_ids))), guides = "collect") +
    plot_annotation(
      title = sprintf("%s in %s and %s cells", gene, TUMOR_CELLTYPE, mask_label),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  fig_w <- 3.6 * length(sample_ids) + 1.5
  fig_h <- 3 * 4.2
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", gene, ".pdf")), p_all, width = fig_w, height = fig_h)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", gene, ".png")), p_all, width = fig_w, height = fig_h, dpi = 300)

  # Selected samples: rows = samples, columns = tumor / CD8+ T / overlay.
  # Same color limits as the full figure, so the two are directly comparable.
  sel <- intersect(SELECTED_SAMPLES, sample_ids)
  if (length(sel) > 0) {
    sel_grid <- unlist(lapply(sel, function(s) {
      s_label <- sprintf("%s (%s)", s, unique(as.character(seu$condition[samples == s])))
      c(list(row_label(s_label)),
        lapply(names(rows), function(r) panels[[s]][[r]] + labs(title = if (s == sel[1]) rows[[r]] else NULL)))
    }), recursive = FALSE)

    p_sel <- wrap_plots(sel_grid, ncol = length(rows) + 1,
                        widths = c(0.15, rep(1, length(rows))), guides = "collect") +
      plot_annotation(
        title = sprintf("%s in %s and %s cells", gene, TUMOR_CELLTYPE, mask_label),
        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

    ggsave(file.path(SELECTED_DIR, paste0("spatial_", gene, ".png")), p_sel,
           width = 3.6 * length(rows) + 1.5, height = 4.2 * length(sel), dpi = 300)
  }
}

# =============================================================================
# Violin plots: expression in CD8+ T cells, one violin per sample
# =============================================================================
VIOLIN_DIR <- file.path(OUTPUT_DIR, "violin_plots")
dir.create(VIOLIN_DIR, recursive = TRUE, showWarnings = FALSE)

mask_idx <- which(is_mask)
vln_df <- data.frame(
  sample    = factor(as.character(samples[mask_idx]), levels = sample_ids),
  condition = factor(as.character(seu$condition[mask_idx]),
                     levels = c(REFERENCE_CONDITION, TEST_CONDITION))
)

violin_layers <- list(
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, alpha = 0.8),
  stat_summary(fun = median, geom = "point", size = 1.2, color = "black", show.legend = FALSE),
  scale_fill_manual(values = CONDITION_COLORS, name = NULL),
  labs(x = NULL, y = "Log-normalized expression"),
  theme_publication(),
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

cat("Violin plots...\n")
for (gene in available_genes) {
  vln_df$expr <- as.numeric(expr_layer[gene, mask_idx])
  p_vln <- ggplot(vln_df, aes(x = sample, y = expr, fill = condition)) +
    violin_layers +
    ggtitle(sprintf("%s in %s cells", gene, mask_label))
  ggsave(file.path(VIOLIN_DIR, paste0("violin_", gene, ".pdf")), p_vln, width = 5, height = 4)
  ggsave(file.path(VIOLIN_DIR, paste0("violin_", gene, ".png")), p_vln, width = 5, height = 4, dpi = 300)
}

# All genes in one panel (free y-scales)
vln_all <- do.call(rbind, lapply(available_genes, function(gene) {
  cbind(vln_df[, c("sample", "condition")], gene = gene,
        expr = as.numeric(expr_layer[gene, mask_idx]))
}))
vln_all$gene <- factor(vln_all$gene, levels = available_genes)
p_vln_all <- ggplot(vln_all, aes(x = sample, y = expr, fill = condition)) +
  violin_layers +
  facet_wrap(~ gene, scales = "free_y", ncol = PANEL_NCOL) +
  ggtitle(sprintf("Trajectory genes in %s cells", mask_label))
vln_h <- ceiling(length(available_genes) / PANEL_NCOL) * 3
ggsave(file.path(VIOLIN_DIR, "violin_all_genes.pdf"), p_vln_all, width = 16, height = vln_h,
       limitsize = FALSE)
ggsave(file.path(VIOLIN_DIR, "violin_all_genes.png"), p_vln_all, width = 16, height = vln_h,
       dpi = 300, limitsize = FALSE)

# =============================================================================
# Dot plot: all genes x condition (WT | DKO) in CD8+ T cells
# =============================================================================
# % expressing (expr > 0) and mean log-normalized expression are computed per
# tumour and then averaged within condition, so each tumour counts equally
# rather than tumours with more CD8+ T cells dominating.
vln_all$expressed <- vln_all$expr > 0
per_tumour <- aggregate(cbind(pct_expressing = expressed * 100, mean_expression = expr) ~
                          gene + condition + sample, data = vln_all, FUN = mean)
dot_df <- aggregate(cbind(pct_expressing, mean_expression) ~ gene + condition,
                    data = per_tumour, FUN = mean)
dot_df$n_tumours <- aggregate(sample ~ gene + condition, data = per_tumour, FUN = length)$sample
write.csv(dot_df[order(dot_df$gene, dot_df$condition), ],
          file.path(VIOLIN_DIR, "dotplot_values.csv"), row.names = FALSE)

# Order genes by the DKO - WT difference in % expressing (higher in DKO on top)
pct_wide  <- tapply(dot_df$pct_expressing, list(dot_df$gene, dot_df$condition), mean)
pct_diff  <- pct_wide[, TEST_CONDITION] - pct_wide[, REFERENCE_CONDITION]
dot_df$gene <- factor(dot_df$gene, levels = names(sort(pct_diff)))

p_dot <- ggplot(dot_df, aes(x = condition, y = gene, size = pct_expressing, color = mean_expression)) +
  geom_point() +
  scale_size_area(max_size = 7, name = "% expressing") +
  scale_color_gradient(low = CD8_GRADIENT[1], high = CD8_GRADIENT[2],
                       name = "Mean log-norm.\nexpression") +
  labs(x = NULL, y = NULL, title = sprintf("Trajectory genes in %s cells", mask_label),
       subtitle = "Per-tumour values averaged within condition") +
  theme_publication() +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3))
dot_h <- 0.22 * length(available_genes) + 1.8
ggsave(file.path(VIOLIN_DIR, "dotplot_all_genes.pdf"), p_dot, width = 5, height = dot_h)
ggsave(file.path(VIOLIN_DIR, "dotplot_all_genes.png"), p_dot, width = 5, height = dot_h, dpi = 300)

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")
