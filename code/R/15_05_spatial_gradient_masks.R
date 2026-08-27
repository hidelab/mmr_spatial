# 15_05_spatial_gradient_masks.R
# Spatial visualization of the tumor-proximity trajectory genes (from 15_04).
#
# Motivation: 15_04 found genes whose CD8+ T expression trajectory *reshapes* in
# DKO — flat/low in WT but surging toward the tumor core in DKO (p_shape hits:
# Nt5e, Bcl2, Icos, Irf1, Gzmb, Prf1, Cd274, Fasl, Il7r, ...). This script shows
# that in real tissue: T cells colored by expression, overlaid on the tumor mass
# (grey) as spatial context, one panel per sample (WT samples first, then DKO).
# A shared per-gene color scale across all panels makes WT-vs-DKO comparable, so
# the DKO-tumor-core-specific pickup is visible directly.
#
# Adapted from 08_01_spatial_masks_08lambda.R.

# --- Thread limits (shared server) ------------------------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS)

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "15_spatial_gradient_masks")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters (modifiable)
# =============================================================================
# Cell type(s) to color by expression. The 15_04 trajectory was on CD8+ T; add
# "Treg" / "NK" here to include more T/lymphoid populations.
MASK_CELLTYPE  <- c("CD8+ T")
TUMOR_CELLTYPE <- "Tumor"          # drawn as grey background context
SHOW_TUMOR_BACKGROUND <- TRUE      # FALSE = only draw the masked cells

# DKO-tumor-core-surge genes (bottom cluster of the 15_04 significant heatmap).
GENES <- c("Nt5e", "Bcl2", "Icos", "Irf1", "Gzmb", "Prf1", "Cd274", "Fasl",
           "Il7r", "Cx3cl1", "Vegfa", "Nfatc1", "Lat2", "Grem1", "Rnf213")

EXPR_CLIP_Q  <- 0.99   # cap color scale at this quantile (per gene) to tame outliers
PT_SIZE_MASK <- 1.1    # masked (T) cell point size
PT_SIZE_BG   <- 0.25   # tumor background point size

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

# =============================================================================
# One panel per sample; grey tumor background + masked cells colored by gene
# =============================================================================
for (gene in available_genes) {
  cat(sprintf("  %s ...\n", gene))
  gvec <- as.numeric(expr_layer[gene, ])

  # Shared color limits across all samples (masked cells only) for comparability
  mask_vals <- gvec[is_mask]
  hi <- as.numeric(quantile(mask_vals[mask_vals > 0], EXPR_CLIP_Q, na.rm = TRUE))
  if (!is.finite(hi) || hi <= 0) hi <- max(mask_vals, na.rm = TRUE)
  if (!is.finite(hi) || hi <= 0) hi <- 1
  lims <- c(0, hi)

  plots <- lapply(sample_ids, function(s) {
    s_idx <- which(samples == s)
    cc    <- center_coords(s_idx)
    m_bg  <- is_tumor[s_idx]
    m_fg  <- is_mask[s_idx]

    df_bg <- data.frame(x = cc[m_bg, 1], y = cc[m_bg, 2])
    df_fg <- data.frame(x = cc[m_fg, 1], y = cc[m_fg, 2],
                        expr = pmin(gvec[s_idx][m_fg], hi))
    df_fg <- df_fg[order(df_fg$expr), ]   # draw high-expressers on top

    p <- ggplot()
    if (SHOW_TUMOR_BACKGROUND && nrow(df_bg) > 0)
      p <- p + geom_point(data = df_bg, aes(x, y), color = "grey80",
                          size = PT_SIZE_BG, alpha = 0.5)
    p <- p +
      geom_point(data = df_fg, aes(x, y, color = expr), size = PT_SIZE_MASK, alpha = 0.9) +
      scale_color_viridis(option = "inferno", name = gene, limits = lims,
                          oob = scales::squish) +
      coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
      ggtitle(sprintf("%s (%s)", s, unique(as.character(seu$condition[s_idx])))) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 10),
            legend.position = "right")
    p
  })

  p_all <- wrap_plots(plots, nrow = 1, guides = "collect") +
    plot_annotation(
      title = sprintf("%s — %s (grey = %s)", gene,
                      paste(MASK_CELLTYPE, collapse = ", "), TUMOR_CELLTYPE),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", gene, ".pdf")),
         p_all, width = 3.6 * length(sample_ids), height = 4.5)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", gene, ".png")),
         p_all, width = 3.6 * length(sample_ids), height = 4.5, dpi = 300)
}

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")
