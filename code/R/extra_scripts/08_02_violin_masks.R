# 08_02_violin_masks.R
# Violin plots of selected genes in masked cells (specific domain + cell type),
# split by condition (WT vs DKO).

library(Seurat)
library(ggplot2)
library(patchwork)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "08_spatial_features_masked", "violins")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Mask parameters (from 00_parameters.R)
# =============================================================================
DOMAIN_COL <- COL_DOMAIN
MASK_DOMAIN <- MASK_DOMAIN_VIOLIN
MASK_CELLTYPE <- MASK_CELLTYPE_VIOLIN

# Genes of interest (combined list from 00_parameters.R)
GENES <- GENES_SPATIAL_VIS

# =============================================================================
# Load merged-domain Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
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

# Set condition as identity for VlnPlot
Idents(seu_masked) <- "condition"

# Filter to genes present in the dataset
available_genes <- intersect(GENES, rownames(seu_masked))
missing_genes <- setdiff(GENES, rownames(seu_masked))
if (length(missing_genes) > 0) {
  cat("Genes not found in dataset (skipping):", paste(missing_genes, collapse = ", "), "\n")
}
cat(sprintf("Plotting %d genes\n", length(available_genes)))

# =============================================================================
# Violin plots split by condition
# =============================================================================

# --- Individual gene plots ---
for (gene in available_genes) {
  cat(sprintf("  %s ...\n", gene))

  p <- VlnPlot(seu_masked, features = gene, group.by = "condition",
               pt.size = 0, cols = CONDITION_COLORS) +
    ggtitle(sprintf("%s — %s, %s", gene, domain_label, celltype_label)) +
    theme_publication() +
    theme(legend.position = "none")

  ggsave(file.path(OUTPUT_DIR, paste0("violin_masked_", gene, ".pdf")),
         p, width = 4, height = 4)
}

# --- Combined multi-panel figure ---
p_combined <- VlnPlot(seu_masked, features = available_genes, group.by = "condition",
                      pt.size = 0, cols = CONDITION_COLORS, ncol = PANEL_NCOL)

ggsave(file.path(OUTPUT_DIR, "violin_masked_all_genes.pdf"),
       p_combined, width = 16, height = ceiling(length(available_genes) / PANEL_NCOL) * 4)
ggsave(file.path(OUTPUT_DIR, "violin_masked_all_genes.png"),
       p_combined, width = 16, height = ceiling(length(available_genes) / PANEL_NCOL) * 4, dpi = 300)

# =============================================================================
# DotPlot split by condition
# =============================================================================
cat("Generating DotPlot...\n")

p_dot <- DotPlot(seu_masked, features = available_genes, group.by = "condition",
                 cols = c(CONDITION_COLORS[[REFERENCE_CONDITION]], CONDITION_COLORS[[TEST_CONDITION]])) +
  ggtitle(sprintf("%s, %s — %s vs %s", domain_label, celltype_label, REFERENCE_CONDITION, TEST_CONDITION)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "dotplot_masked_all_genes.pdf"),
       p_dot, width = 12, height = 4)
ggsave(file.path(OUTPUT_DIR, "dotplot_masked_all_genes.png"),
       p_dot, width = 12, height = 4, dpi = 300)

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")