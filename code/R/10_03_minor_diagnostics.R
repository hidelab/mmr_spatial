# 10_03_minor_diagnostics.R
# Check if spatial coordinates are at the same scale across samples
# Focus on one marker and two samples

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)

# Source centralized parameters
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "10_minor_diagnostics")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Load Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

DefaultAssay(seu) <- "RNA"

# =============================================================================
# Select one marker and two samples for detailed coordinate inspection
# =============================================================================
MARKER_NAME <- "Cd4"
MARKER_EXPR <- "cd4_expr > 0"
SAMPLES_TO_CHECK <- c("WT1", "DKO2")

cat("\n=== Checking coordinates for marker:", MARKER_NAME, "in samples:", paste(SAMPLES_TO_CHECK, collapse=", "), "===\n")

# Get expression matrix
counts_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")
cd4_expr <- counts_mat["Cd4", ]

# Subset to selected samples
seu_subset <- subset(seu, sample %in% SAMPLES_TO_CHECK)

# Add marker expression to metadata
seu_subset$marker_positive <- cd4_expr[colnames(seu_subset)] > 0

# =============================================================================
# Extract and compare coordinate ranges
# =============================================================================
cat("\nCoordinate ranges by sample:\n")

coord_summary <- NULL

for (samp in SAMPLES_TO_CHECK) {
  samp_subset <- subset(seu_subset, sample == samp)
  coords <- Embeddings(samp_subset, reduction = "spatial")

  x_range <- range(coords[, 1])
  y_range <- range(coords[, 2])

  cat(sprintf("\n%s:\n", samp))
  cat(sprintf("  X range: %.2f to %.2f (span: %.2f)\n", x_range[1], x_range[2], x_range[2] - x_range[1]))
  cat(sprintf("  Y range: %.2f to %.2f (span: %.2f)\n", y_range[1], y_range[2], y_range[2] - y_range[1]))

  coord_summary <- rbind(coord_summary, data.frame(
    sample = samp,
    x_min = x_range[1],
    x_max = x_range[2],
    x_span = x_range[2] - x_range[1],
    y_min = y_range[1],
    y_max = y_range[2],
    y_span = y_range[2] - y_range[1],
    n_cells = ncol(samp_subset),
    n_marker_positive = sum(samp_subset$marker_positive)
  ))
}

print(coord_summary)

# Save coordinate summary
write.csv(coord_summary, file.path(OUTPUT_DIR, "coordinate_ranges_sample_comparison.csv"), row.names = FALSE)

# =============================================================================
# Create spatial plots using Seurat with matched unit scale
# =============================================================================
cat("\nGenerating spatial plots with matched unit scale...\n")

# Determine the max span across both samples so both plots use the same scale
max_x_span <- max(coord_summary$x_span)
max_y_span <- max(coord_summary$y_span)

cat(sprintf("\nMax X span: %.2f, Max Y span: %.2f\n", max_x_span, max_y_span))

plot_list <- list()

for (samp in SAMPLES_TO_CHECK) {
  samp_obj <- subset(seu_subset, sample == samp)

  # Get this sample's coordinate center
  coords <- Embeddings(samp_obj, reduction = "spatial")
  x_center <- mean(range(coords[, 1]))
  y_center <- mean(range(coords[, 2]))

  # Set axis limits centered on the sample but with the max span (+ 5% padding)
  pad <- 0.05
  x_lim <- c(x_center - max_x_span * (0.5 + pad), x_center + max_x_span * (0.5 + pad))
  y_lim <- c(y_center - max_y_span * (0.5 + pad), y_center + max_y_span * (0.5 + pad))

  # Use Seurat DimPlot on spatial reduction, coloring by marker_positive
  samp_obj$marker_status <- ifelse(samp_obj$marker_positive, paste0(MARKER_NAME, "+"), "Negative")
  samp_obj$marker_status <- factor(samp_obj$marker_status, levels = c("Negative", paste0(MARKER_NAME, "+")))

  p <- DimPlot(samp_obj, reduction = "spatial", group.by = "marker_status",
               cols = c("Negative" = "lightgray", "Cd4+" = "#e74c3c"),
               pt.size = 1.5) +
    coord_fixed(xlim = x_lim, ylim = y_lim) +
    ggtitle(paste0(samp, " - ", MARKER_NAME, " Expression")) +
    xlab("Spatial X") + ylab("Spatial Y") +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

  plot_list[[samp]] <- p
}

# Combine plots
combined_plot <- wrap_plots(plot_list, ncol = 2)

# Save as JPG
ggsave(file.path(OUTPUT_DIR, "spatial_marker_comparison_cd4_wt1_dko2.jpg"),
       combined_plot, width = 14, height = 7, dpi = 300, units = "in")

cat("\nDone. Output saved to:", OUTPUT_DIR, "\n")