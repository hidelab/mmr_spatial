# 08_04_violin_x_cells.R
# Violin plots of selected genes with cell types on x-axis, grouped/colored by domain.
# Each gene produces a separate plot showing expression across cell types, faceted by domain.

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

PROJECT_DIR <- "/data/work/Pourya/mmr_spatial"
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "08_violin_x_cells")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Parameters
# Leave MASK_DOMAIN empty (c()) or NULL to use all domains
# =============================================================================
DOMAIN_COL <- "banksy_domain_merged"
DOMAIN_NAMES <- c(
  "2" = "Tumor_Core",
  "3" = "Immune_Engulfing",
  "4" = "Stroma"
)

# Domains to include (use domain IDs). NULL or c() = all domains.
MASK_DOMAIN <- c("2", "4")

# Cell types to include. NULL or c() = all cell types.
MASK_CELLTYPE <- c(
  "CD8+ T",
  "Treg",
  "NK",
  "cDC1",
  "pDC",
  "DC (Ccr7+)",
  "Macrophage (Cxcl16+)",
  "Neutrophil"
)

# Cell-level QC filters (for pseudobulk)
MIN_TOTAL_COUNTS <- 50
MIN_GENES_DETECTED <- 10
MIN_CELLS_PER_SAMPLE <- 10

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
  "Pdgfra", "Fcgr4", "Lig1",
  # DEG Neutrophil
  "Cxcl2","Ccl3","Nlrp3","Cd274", "Il1b",
  # DEG CD8 T 
  "Ifng","Lag3"
)

# =============================================================================
# Load merged-domain Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", "DKO", seu$condition)
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# Use RNA assay for expression
DefaultAssay(seu) <- "RNA"

# Create named domain column
seu$domain_id <- as.character(seu@meta.data[[DOMAIN_COL]])
seu$domain_name <- ifelse(seu$domain_id %in% names(DOMAIN_NAMES),
                          DOMAIN_NAMES[seu$domain_id],
                          paste0("Domain_", seu$domain_id))

# =============================================================================
# Apply mask
# =============================================================================
mask <- rep(TRUE, ncol(seu))

if (!is.null(MASK_DOMAIN) && length(MASK_DOMAIN) > 0) {
  mask <- mask & seu$domain_id %in% MASK_DOMAIN
}
if (!is.null(MASK_CELLTYPE) && length(MASK_CELLTYPE) > 0) {
  mask <- mask & seu$Tier1_celltype %in% MASK_CELLTYPE
}

seu_masked <- subset(seu, cells = colnames(seu)[mask])

domain_label <- if (!is.null(MASK_DOMAIN) && length(MASK_DOMAIN) > 0) {
  paste("Domains:", paste(DOMAIN_NAMES[MASK_DOMAIN], collapse = ", "))
} else {
  "All domains"
}
celltype_label <- if (!is.null(MASK_CELLTYPE) && length(MASK_CELLTYPE) > 0) {
  paste(length(MASK_CELLTYPE), "cell types")
} else {
  "All cell types"
}
cat(sprintf("Masked to %s + %s: %d cells (from %d)\n",
            domain_label, celltype_label, ncol(seu_masked), ncol(seu)))

# Set factor levels for domain_name
domain_levels <- if (!is.null(MASK_DOMAIN) && length(MASK_DOMAIN) > 0) {
  DOMAIN_NAMES[MASK_DOMAIN]
} else {
  sort(unique(seu_masked$domain_name))
}
seu_masked$domain_name <- factor(seu_masked$domain_name, levels = domain_levels)

# Set factor levels for cell type (follow TIER1_ORDER)
ct_levels <- if (!is.null(MASK_CELLTYPE) && length(MASK_CELLTYPE) > 0) {
  intersect(TIER1_ORDER, MASK_CELLTYPE)
} else {
  intersect(TIER1_ORDER, unique(seu_masked$Tier1_celltype))
}
seu_masked$Tier1_celltype <- factor(seu_masked$Tier1_celltype, levels = ct_levels)

# Domain colors
n_domains <- length(domain_levels)
domain_colors <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")[seq_len(n_domains)],
  domain_levels
)

# Filter to genes present in the dataset
available_genes <- intersect(GENES, rownames(seu_masked))
missing_genes <- setdiff(GENES, rownames(seu_masked))
if (length(missing_genes) > 0) {
  cat("Genes not found in dataset (skipping):", paste(missing_genes, collapse = ", "), "\n")
}
cat(sprintf("Plotting %d genes across %d cell types and %d domains\n",
            length(available_genes), length(ct_levels), n_domains))

# Report cell counts
cat("\nCells per cell type x domain:\n")
print(table(seu_masked$Tier1_celltype, seu_masked$domain_name))
cat("\n")

# =============================================================================
# Extract expression data and build combined data frame
# =============================================================================
expr_data <- FetchData(seu_masked, vars = available_genes)
expr_data$celltype <- seu_masked$Tier1_celltype
expr_data$domain <- seu_masked$domain_name
expr_data$condition <- seu_masked$condition

# =============================================================================
# Violin plots: x = cell type, fill = domain
# =============================================================================
cat("Generating violin plots...\n")

for (gene in available_genes) {
  cat(sprintf("  %s ...\n", gene))

  p <- ggplot(expr_data, aes(x = celltype, y = .data[[gene]], fill = domain)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, alpha = 0.8) +
    scale_fill_manual(values = domain_colors, name = "Domain") +
    labs(title = gene,
         x = NULL,
         y = "Expression") +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "top"
    )

  ggsave(file.path(OUTPUT_DIR, paste0("violin_xcell_", gene, ".pdf")),
         p, width = max(6, length(ct_levels) * 1.2), height = 5)
}

# =============================================================================
# Combined multi-panel figure
# =============================================================================
cat("Generating combined panel...\n")

plot_list <- lapply(available_genes, function(gene) {
  ggplot(expr_data, aes(x = celltype, y = .data[[gene]], fill = domain)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, alpha = 0.8) +
    scale_fill_manual(values = domain_colors, name = "Domain") +
    labs(title = gene, x = NULL, y = NULL) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
})

# Create a shared legend
legend_plot <- ggplot(expr_data, aes(x = celltype, y = .data[[available_genes[1]]], fill = domain)) +
  geom_violin() +
  scale_fill_manual(values = domain_colors, name = "Domain") +
  theme_publication() +
  theme(legend.position = "bottom")
shared_legend <- cowplot::get_legend(legend_plot)

n_cols <- 4
n_rows <- ceiling(length(available_genes) / n_cols)

p_combined <- wrap_plots(plot_list, ncol = n_cols) +
  plot_annotation(
    title = sprintf("Gene expression by cell type and domain"),
    subtitle = sprintf("%s | %s", domain_label, celltype_label),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  )

plot_width <- max(16, length(ct_levels) * n_cols * 0.8)
plot_height <- n_rows * 4 + 1

ggsave(file.path(OUTPUT_DIR, "violin_xcell_all_genes.pdf"),
       p_combined, width = plot_width, height = plot_height)
ggsave(file.path(OUTPUT_DIR, "violin_xcell_all_genes.png"),
       p_combined, width = plot_width, height = plot_height, dpi = 300)

# =============================================================================
# Pseudobulk boxplots: x = cell type, fill = domain
# Aggregate per sample x domain x celltype, then plot per-sample means
# =============================================================================
cat("\n--- Pseudobulk boxplots ---\n")

pb_dir <- file.path(OUTPUT_DIR, "boxplots_pseudobulk")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)

# Apply QC filter on cells
keep_qc <- seu_masked$total_counts >= MIN_TOTAL_COUNTS &
            seu_masked$n_genes_by_counts >= MIN_GENES_DETECTED
seu_qc <- subset(seu_masked, cells = colnames(seu_masked)[keep_qc])
cat(sprintf("After QC (total_counts >= %d, n_genes >= %d): %d cells (removed %d)\n",
            MIN_TOTAL_COUNTS, MIN_GENES_DETECTED, ncol(seu_qc),
            ncol(seu_masked) - ncol(seu_qc)))

# Get counts matrix
counts_mat <- GetAssayData(seu_qc, assay = "RNA", layer = "counts")

# Build pseudobulk groups: sample x domain x celltype
seu_qc$pb_group <- paste(seu_qc$sample, seu_qc$domain_name, seu_qc$Tier1_celltype, sep = "__")
group_table <- table(seu_qc$pb_group)

# Aggregate only groups with enough cells
valid_groups <- names(group_table[group_table >= MIN_CELLS_PER_SAMPLE])
cat(sprintf("Pseudobulk groups with >= %d cells: %d (of %d total)\n",
            MIN_CELLS_PER_SAMPLE, length(valid_groups), length(group_table)))

pb_list <- list()
for (g in valid_groups) {
  cells_g <- colnames(seu_qc)[seu_qc$pb_group == g]
  pb_list[[g]] <- Matrix::rowSums(counts_mat[, cells_g, drop = FALSE])
}

pb_mat <- do.call(cbind, pb_list)

# Normalize to CPM
lib_sizes <- colSums(pb_mat)
pb_cpm <- t(t(pb_mat) / lib_sizes * 1e6)
pb_log <- log1p(pb_cpm)

# Build pseudobulk metadata
pb_meta <- data.frame(
  group = valid_groups,
  sample = sapply(strsplit(valid_groups, "__"), `[`, 1),
  domain = sapply(strsplit(valid_groups, "__"), `[`, 2),
  celltype = sapply(strsplit(valid_groups, "__"), `[`, 3),
  n_cells = as.numeric(group_table[valid_groups]),
  stringsAsFactors = FALSE
)
pb_meta$condition <- ifelse(grepl("^WT", pb_meta$sample), "WT", "DKO")
pb_meta$domain <- factor(pb_meta$domain, levels = domain_levels)
pb_meta$celltype <- factor(pb_meta$celltype, levels = ct_levels)

cat("\nPseudobulk samples per celltype x domain:\n")
print(table(pb_meta$celltype, pb_meta$domain))
cat("\n")

# Build long-format expression data for available genes
pb_genes <- intersect(available_genes, rownames(pb_log))
pb_expr <- as.data.frame(t(pb_log[pb_genes, , drop = FALSE]))
pb_expr$group <- rownames(pb_expr)
pb_long <- merge(pb_expr, pb_meta, by = "group")

# --- Individual gene boxplots ---
cat("Generating pseudobulk boxplots...\n")

for (gene in pb_genes) {
  cat(sprintf("  %s ...\n", gene))

  p <- ggplot(pb_long, aes(x = celltype, y = .data[[gene]], fill = domain)) +
    geom_boxplot(outlier.size = 0.8, alpha = 0.8, linewidth = 0.4) +
    geom_point(aes(group = domain), position = position_dodge(width = 0.75),
               size = 1.2, alpha = 0.6) +
    scale_fill_manual(values = domain_colors, name = "Domain") +
    labs(title = sprintf("%s (pseudobulk)", gene),
         x = NULL,
         y = "log(CPM + 1)") +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "top"
    )

  ggsave(file.path(pb_dir, paste0("boxplot_pb_", gene, ".pdf")),
         p, width = max(6, length(ct_levels) * 1.2), height = 5)
}

# --- Combined multi-panel boxplot figure ---
cat("Generating combined pseudobulk panel...\n")

pb_plot_list <- lapply(pb_genes, function(gene) {
  ggplot(pb_long, aes(x = celltype, y = .data[[gene]], fill = domain)) +
    geom_boxplot(outlier.size = 0.5, alpha = 0.8, linewidth = 0.3) +
    geom_point(aes(group = domain), position = position_dodge(width = 0.75),
               size = 0.8, alpha = 0.5) +
    scale_fill_manual(values = domain_colors, name = "Domain") +
    labs(title = gene, x = NULL, y = NULL) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      plot.title = element_text(size = 10, face = "bold"),
      legend.position = "none"
    )
})

n_cols <- 4
n_rows <- ceiling(length(pb_genes) / n_cols)

p_pb_combined <- wrap_plots(pb_plot_list, ncol = n_cols) +
  plot_annotation(
    title = "Pseudobulk gene expression by cell type and domain",
    subtitle = sprintf("%s | %s", domain_label, celltype_label),
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5)
    )
  )

plot_width <- max(16, length(ct_levels) * n_cols * 0.8)
plot_height <- n_rows * 4 + 1

ggsave(file.path(pb_dir, "boxplot_pb_all_genes.pdf"),
       p_pb_combined, width = plot_width, height = plot_height)
ggsave(file.path(pb_dir, "boxplot_pb_all_genes.png"),
       p_pb_combined, width = plot_width, height = plot_height, dpi = 300)

cat("\nDone. All plots saved to:", OUTPUT_DIR, "\n")