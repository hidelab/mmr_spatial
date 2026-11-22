# 10_2_minor_diagnostics.R
# Generate marker expression summary table with user-defined markers
# Markers are specified at the beginning of the script

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "10_minor_diagnostics")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# USER-DEFINED MARKERS
# =============================================================================
# Define markers as a list of expressions to evaluate
# Each element should be a named list with:
#   - name: character string for the marker name (display name)
#   - expr: character string of R expression to evaluate (e.g., "cd4_expr > 0")
#
# Example:
#   MARKERS <- list(
#     list(name = "Cd4", expr = "cd4_expr > 0"),
#     list(name = "Cd19", expr = "cd19_expr > 0"),
#     list(name = "Cd4+Cd3d+", expr = "cd4_expr > 0 & cd3d_expr > 0")
#   )

MARKERS <- list(
  list(name = "Cd4", expr = "cd4_expr > 0"),
  list(name = "Cd19", expr = "cd19_expr > 0"),
  list(name = "Cd8a", expr = "cd8a_expr > 0"),
  list(name = "Cd4+Cd3d+", expr = "cd4_expr > 0 & cd3d_expr > 0"),
  list(name = "Cd4+Ptprc+", expr = "cd4_expr > 0 & ptprc_expr > 0"),
  list(name = "Cd4+Cd3d+Ptprc+", expr = "cd4_expr > 0 & cd3d_expr > 0 & ptprc_expr > 0"),
  list(name = "Cd8+Cd3d+", expr = "cd8a_expr > 0 & cd3d_expr > 0"),
  list(name = "Cd8+Ptprc+", expr = "cd8a_expr > 0 & ptprc_expr > 0"),
  list(name = "Cd8+Cd3d+Ptprc+", expr = "cd8a_expr > 0 & cd3d_expr > 0 & ptprc_expr > 0"),
  list(name = "Cd11b+ (Macrophage)", expr = "cd11b_expr > 0"),
  list(name = "Itgax+ (Dendritic Cell)", expr = "itgax_expr > 0"),
  list(name = "Csf3r+ (Neurtophil)", expr = "csf3r_expr > 0"),
  list(name = "Fbn1+ (Fibroblast)", expr = "fbn1_expr > 0"),
  list(name = "Klrb1c+ (NK cell)", expr = "klrb1c_expr > 0"),
  list(name = "tnfrsf12a+ (Tumor)", expr = "tnfrsf12a_expr > 0")
)

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
# Get expression matrix and define gene expressions
# =============================================================================
counts_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")

# Define gene expressions for all markers used in MARKERS list
cd4_expr <- counts_mat["Cd4", ]
cd19_expr <- counts_mat["Cd19", ]
cd8a_expr <- counts_mat["Cd8a", ]
cd3d_expr <- counts_mat["Cd3d", ]
ptprc_expr <- counts_mat["Ptprc", ]
cd11b_expr <- counts_mat["Itgam", ]
itgax_expr <- counts_mat["Itgax", ]
csf3r_expr <- counts_mat["Csf3r",]
fbn1_expr  <- counts_mat["Fbn1",]
klrb1c_expr <- counts_mat["Klrb1c",]
tnfrsf12a_expr <- counts_mat["Tnfrsf12a",]
# =============================================================================
# 5. Generate Summary Table
# =============================================================================
cat("\n=== Generating Summary Table ===\n")

summary_df <- data.frame(
  marker = character(length(MARKERS)),
  n_cells = integer(length(MARKERS)),
  pct_total = numeric(length(MARKERS)),
  stringsAsFactors = FALSE
)

for (i in seq_along(MARKERS)) {
  marker_name <- MARKERS[[i]]$name
  marker_expr <- MARKERS[[i]]$expr

  # Evaluate the expression in the current environment
  pos <- eval(parse(text = marker_expr))
  n_pos <- sum(pos)
  pct_pos <- round(100 * n_pos / ncol(seu), 3)

  summary_df$marker[i] <- marker_name
  summary_df$n_cells[i] <- n_pos
  summary_df$pct_total[i] <- pct_pos

  cat(sprintf("  %s: %d cells (%.2f%%)\n", marker_name, n_pos, pct_pos))
}

# =============================================================================
# Cell Type Composition from Tier1_celltype annotation
# =============================================================================
cat("\n\n=== Cell Type Composition (Tier1_celltype) ===\n")

celltype_table <- table(seu$Tier1_celltype)
celltype_pct <- round(100 * celltype_table / ncol(seu), 3)

celltype_df <- data.frame(
  marker = names(celltype_table),
  n_cells = as.integer(celltype_table),
  pct_total = as.numeric(celltype_pct),
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Sort by n_cells descending
celltype_df <- celltype_df[order(celltype_df$n_cells, decreasing = TRUE), ]

print(celltype_df)

# =============================================================================
# Combine tables for comparison
# =============================================================================
cat("\n\n=== Combined Summary (Markers + Cell Types) ===\n")

combined_df <- rbind(summary_df, celltype_df)

print(combined_df)

# =============================================================================
# Confusion Matrix: Markers vs Tier1_celltype
# =============================================================================
cat("\n\n=== Confusion Matrix (Markers vs Tier1_celltype) ===\n")

# Get unique cell types
celltypes <- sort(unique(seu$Tier1_celltype))

# Initialize confusion matrix
confusion_matrix <- matrix(0, nrow = length(MARKERS), ncol = length(celltypes))
rownames(confusion_matrix) <- sapply(MARKERS, function(x) x$name)
colnames(confusion_matrix) <- celltypes

# Fill confusion matrix
for (i in seq_along(MARKERS)) {
  marker_expr <- MARKERS[[i]]$expr
  pos <- eval(parse(text = marker_expr))

  for (j in seq_along(celltypes)) {
    celltype_match <- seu$Tier1_celltype == celltypes[j]
    confusion_matrix[i, j] <- sum(pos & celltype_match)
  }
}

print(confusion_matrix)

# Convert to data frame for saving
confusion_df <- as.data.frame(confusion_matrix)
confusion_df$marker <- rownames(confusion_matrix)
confusion_df <- confusion_df[, c(ncol(confusion_df), 1:(ncol(confusion_df)-1))]

write.csv(confusion_df, file.path(OUTPUT_DIR, "confusion_matrix_markers_vs_celltypes.csv"), row.names = FALSE)

# Save both individual and combined tables
write.csv(summary_df, file.path(OUTPUT_DIR, "marker_positive_summary_custom.csv"), row.names = FALSE)
write.csv(celltype_df, file.path(OUTPUT_DIR, "celltype_composition_tier1.csv"), row.names = FALSE)
write.csv(combined_df, file.path(OUTPUT_DIR, "combined_marker_celltype_summary.csv"), row.names = FALSE)

# =============================================================================
# Save list of all expressed genes
# =============================================================================
cat("\n\n=== Extracting Expressed Genes ===\n")

# Get all genes from the counts matrix
all_genes <- rownames(counts_mat)

# Calculate expression statistics per gene
gene_stats <- data.frame(
  gene = all_genes,
  n_cells_expressed = rowSums(counts_mat > 0),
  mean_expression = Matrix::rowMeans(counts_mat),
  max_expression = apply(counts_mat, 1, max),
  stringsAsFactors = FALSE
)

# Sort by number of cells expressing the gene
gene_stats <- gene_stats[order(gene_stats$n_cells_expressed, decreasing = TRUE), ]

cat(sprintf("Total genes: %d\n", nrow(gene_stats)))
cat(sprintf("Genes expressed in at least 1 cell: %d\n", sum(gene_stats$n_cells_expressed > 0)))
cat(sprintf("Genes expressed in at least 10 cells: %d\n", sum(gene_stats$n_cells_expressed >= 10)))

write.csv(gene_stats, file.path(OUTPUT_DIR, "all_expressed_genes.csv"), row.names = FALSE)

cat("\nDone. Output saved to:", OUTPUT_DIR, "\n")