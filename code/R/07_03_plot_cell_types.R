# 06_plot_cell_types.R
# Spatial plots highlighting each Tier1 cell type individually.
# For each cell type, cells belonging to that type are colored red,
# all other cells are gray. One multi-panel figure per cell type (panels = samples).

library(Seurat)
library(ggplot2)
library(patchwork)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "celltype_plots")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load Seurat object ---
seu <- readRDS(SEURAT_RDS)
cat("Loaded:", ncol(seu), "cells\n")

celltypes <- sort(unique(seu$Tier1_celltype))
sample_ids <- unique(seu$sample)
coords <- Embeddings(seu, "spatial")

cat(sprintf("Plotting %d cell types across %d samples\n", length(celltypes), length(sample_ids)))

for (ct in celltypes) {
  cat(sprintf("  %s ...\n", ct))

  plots <- lapply(sample_ids, function(s) {
    idx <- which(seu$sample == s)
    df <- data.frame(
      x = coords[idx, 1],
      y = coords[idx, 2],
      highlight = ifelse(seu$Tier1_celltype[idx] == ct, "CellType", "Other")
    )
    # Plot "Other" first (behind), then highlighted on top
    df$highlight <- factor(df$highlight, levels = c("Other", "CellType"))
    df <- df[order(df$highlight), ]

    ggplot(df, aes(x, y, color = highlight)) +
      geom_point(size = 0.1) +
      scale_color_manual(values = c("Other" = "grey80", "CellType" = "red")) +
      coord_fixed() +
      theme_minimal() +
      ggtitle(s) +
      theme(legend.position = "none",
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank())
  })

  p <- wrap_plots(plots, ncol = 4) +
    plot_annotation(title = ct,
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))

  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_celltype_", ct_safe, ".pdf")),
         p, width = 16, height = 12)
}

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")