# 08_01_spatial_masks_08lambda.R
# Spatial feature plots masked to cells in specific domain(s) AND cell type.
# For each gene, produces a multi-panel figure (one panel per sample),
# showing only the masked cells colored by expression level.

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

PROJECT_DIR <- "/data/work/Pourya/mmr_spatial"
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "08_spatial_features_masked")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Mask parameters
# Leave empty (c()) or NULL to use all domains / all cell types
# =============================================================================
DOMAIN_COL <- "banksy_domain_merged"
MASK_DOMAIN <- c(2)
MASK_CELLTYPE <- c()

# Cell-level QC filters
MIN_TOTAL_COUNTS <- 50
MIN_GENES_DETECTED <- 10

# =============================================================================
# Selected genes of interest
# =============================================================================
GENES <- c(
  # Cytotoxicity
  "Klrk1", "Gzmb", "Prf1",
  # Chemokine axis
  "Cxcl16", "Cxcr6",
  # Cytokine signaling
  "Il15",
  # Immune markers
  "Cd8a", "Cd4", "Foxp3", "Nkg7",
  # Macrophage polarization
  "Mrc1", "Cd163", "Nos2",
  # Tumor / proliferation
  "Mki67", "Top2a", "Tigit", "Arg1",
  # Some DEG genes
  "Osm", "Osmr", "Il24",
  "Pdgfra", "Fcgr4", "Lig1"
)

# =============================================================================
# Load merged-domain Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", "DKO", seu$condition)
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# Use RNA assay for expression
DefaultAssay(seu) <- "RNA"

# =============================================================================
# Apply mask: specified domain(s) AND/OR cell type(s)
# If either is NULL or empty, no filtering on that axis (use all).
# =============================================================================
mask <- rep(TRUE, ncol(seu))

if (!is.null(MASK_DOMAIN) && length(MASK_DOMAIN) > 0) {
  mask <- mask & as.character(seu@meta.data[[DOMAIN_COL]]) %in% MASK_DOMAIN
}
if (!is.null(MASK_CELLTYPE) && length(MASK_CELLTYPE) > 0) {
  mask <- mask & seu$Tier1_celltype %in% MASK_CELLTYPE
}

seu_masked <- subset(seu, cells = colnames(seu)[mask])

# Apply QC filter
n_before_qc <- ncol(seu_masked)
keep_qc <- seu_masked$total_counts >= MIN_TOTAL_COUNTS &
            seu_masked$n_genes_by_counts >= MIN_GENES_DETECTED
seu_masked <- subset(seu_masked, cells = colnames(seu_masked)[keep_qc])
cat(sprintf("After QC (total_counts >= %d, n_genes >= %d): %d cells (removed %d)\n",
            MIN_TOTAL_COUNTS, MIN_GENES_DETECTED, ncol(seu_masked),
            n_before_qc - ncol(seu_masked)))

domain_label <- if (!is.null(MASK_DOMAIN) && length(MASK_DOMAIN) > 0) {
  paste("Domain", paste(MASK_DOMAIN, collapse = "/"))
} else {
  "All domains"
}
celltype_label <- if (!is.null(MASK_CELLTYPE) && length(MASK_CELLTYPE) > 0) {
  paste(MASK_CELLTYPE, collapse = ", ")
} else {
  "All cell types"
}
cat(sprintf("Masked to %s + %s: %d cells (from %d)\n",
            domain_label, celltype_label, ncol(seu_masked), ncol(seu)))

sample_ids <- levels(seu_masked$sample)
sample_ids <- sample_ids[sample_ids %in% unique(seu_masked$sample)]

# Filter to genes present in the dataset
available_genes <- intersect(GENES, rownames(seu_masked))
missing_genes <- setdiff(GENES, rownames(seu_masked))
if (length(missing_genes) > 0) {
  cat("Genes not found in dataset (skipping):", paste(missing_genes, collapse = ", "), "\n")
}
cat(sprintf("Plotting %d genes across %d samples\n", length(available_genes), length(sample_ids)))

# =============================================================================
# Plot each gene across all samples (masked cells only)
# =============================================================================
for (gene in available_genes) {
  cat(sprintf("  %s ...\n", gene))

  plots <- lapply(sample_ids, function(s) {
    sub <- subset(seu_masked, subset = sample == s)
    FeaturePlot(sub, features = gene, reduction = "spatial", pt.size = 3,
                alpha = 0.3, raster = TRUE, raster.dpi = c(300, 300)) +
      scale_color_viridis(option = "inferno", name = gene) +
      ggtitle(s) +
      NoAxes() +
      coord_fixed()
  })

  p <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = sprintf("%s — %s, %s", gene, domain_label, celltype_label),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  ggsave(file.path(OUTPUT_DIR, paste0("spatial_masked_", gene, ".pdf")),
         p, width = 4 * length(sample_ids), height = 5)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_masked_", gene, ".png")),
         p, width = 4 * length(sample_ids), height = 5, dpi = 300)
}

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")