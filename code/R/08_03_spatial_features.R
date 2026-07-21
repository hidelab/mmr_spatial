# 08_03_spatial_features.R
# Spatial feature plots showing expression of selected genes across all samples.
# For each gene, produces a multi-panel figure (one panel per sample),
# colored by expression level on spatial coordinates.

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

PROJECT_DIR <- "/data/work/Pourya/mmr_spatial"
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "08_spatial_features_masked", "all_cells")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

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
  "Osm", "Osmr", "Il24"
)

# =============================================================================
# Load Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", "DKO", seu$condition)
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# Use RNA assay for expression
DefaultAssay(seu) <- "RNA"

sample_ids <- levels(seu$sample)

# Filter to genes present in the dataset
available_genes <- intersect(GENES, rownames(seu))
missing_genes <- setdiff(GENES, rownames(seu))
if (length(missing_genes) > 0) {
  cat("Genes not found in dataset (skipping):", paste(missing_genes, collapse = ", "), "\n")
}
cat(sprintf("Plotting %d genes across %d samples\n", length(available_genes), length(sample_ids)))

# =============================================================================
# Plot each gene across all samples
# =============================================================================
for (gene in available_genes) {
  cat(sprintf("  %s ...\n", gene))

  plots <- lapply(sample_ids, function(s) {
    sub <- subset(seu, subset = sample == s)
    FeaturePlot(sub, features = gene, reduction = "spatial", pt.size = 1,
                order = TRUE, raster = TRUE, raster.dpi = c(300, 300)) +
      scale_color_viridis(option = "magma", name = gene) +
      ggtitle(s) +
      NoAxes() +
      coord_fixed()
  })

  p <- wrap_plots(plots, nrow = 1) +
    plot_annotation(title = gene,
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  ggsave(file.path(OUTPUT_DIR, paste0("spatial_feature_", gene, ".pdf")),
         p, width = 4 * length(sample_ids), height = 5)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_feature_", gene, ".png")),
         p, width = 4 * length(sample_ids), height = 5, dpi = 300)
}

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")