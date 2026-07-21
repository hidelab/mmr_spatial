# 01_data_import.R
# Import DKO_Tier1.h5ad into Seurat, inspect, and produce overview plots.

library(Seurat)
library(anndataR)
library(ggplot2)
library(patchwork)
library(cowplot)

PROJECT_DIR <- "/data/work/Pourya/mmr_spatial"
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "01_summary_exploratory")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- 1. Load h5ad via anndataR -------------------------------------------------
h5ad_path <- file.path(PROJECT_DIR, "input", "DKO_results", "DKO_Tier1.h5ad")
adata <- read_h5ad(h5ad_path)

# --- 2. Convert to Seurat ------------------------------------------------------
sobj <- adata$as_Seurat()

# Add spatial coordinates
coords <- adata$obsm[["spatial"]]
colnames(coords) <- c("spatial_1", "spatial_2")
sobj[["spatial"]] <- CreateDimReducObject(
  embeddings = coords,
  key = "spatial_",
  assay = DefaultAssay(sobj)
)

# Transfer UMAP if present
if ("X_umap" %in% names(adata$obsm)) {
  umap_coords <- adata$obsm[["X_umap"]]
  colnames(umap_coords) <- c("umap_1", "umap_2")
  sobj[["umap"]] <- CreateDimReducObject(
    embeddings = umap_coords,
    key = "umap_",
    assay = DefaultAssay(sobj)
  )
}

# Recode condition: "KO" → "DKO" for display consistency
sobj$condition <- gsub("^KO$", "DKO", sobj$condition)

# Set factor levels for consistent legend ordering across all plots
# Only include levels that actually exist in the data to avoid NAs
sobj$Tier1_celltype <- factor(sobj$Tier1_celltype,
                              levels = intersect(TIER1_ORDER, unique(sobj$Tier1_celltype)))
sobj$condition <- factor(sobj$condition,
                         levels = intersect(names(CONDITION_COLORS), unique(sobj$condition)))
sobj$sample <- factor(sobj$sample,
                      levels = intersect(names(SAMPLE_COLORS), unique(sobj$sample)))

# --- 3. Basic inspection -------------------------------------------------------
cat("== Object overview ==\n")
print(sobj)

cat("\n== Metadata columns ==\n")
print(colnames(sobj@meta.data))

cat("\n== Cells per sample ==\n")
print(table(sobj$sample))

cat("\n== Cells per condition ==\n")
print(table(sobj$condition))

cat("\n== Cell types (Tier1) ==\n")
print(sort(table(sobj$Tier1_celltype), decreasing = TRUE))

cat("\n== QC summary ==\n")
print(summary(sobj$total_counts))
print(summary(sobj$n_genes_by_counts))

# --- 4. Visualizations ----------------------------------------------------------

if ("umap" %in% Reductions(sobj)) {

  # UMAP by condition
  p1 <- DimPlot(sobj, reduction = "umap", group.by = "condition", pt.size = 1,
                cols = CONDITION_COLORS, shuffle = TRUE,
                raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle("Condition") +
    NoAxes() +
    guides(color = guide_legend(override.aes = list(size = 3)))

  # UMAP by cell type (labeled, no legend)
  p2 <- DimPlot(sobj, reduction = "umap", group.by = "Tier1_celltype", pt.size = 1,
                cols = TIER1_COLORS, label = TRUE, label.size = 3, repel = TRUE,
                shuffle = TRUE, raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle("Cell Type") +
    NoAxes() +
    NoLegend()

  # UMAP by cell type with legend (for reference)
  p2_legend <- DimPlot(sobj, reduction = "umap", group.by = "Tier1_celltype", pt.size = 1,
                       cols = TIER1_COLORS, shuffle = TRUE,
                       raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle("Cell Type") +
    NoAxes() +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

  # UMAP by sample
  p3 <- DimPlot(sobj, reduction = "umap", group.by = "sample", pt.size = 1,
                cols = SAMPLE_COLORS, shuffle = TRUE,
                raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle("Sample") +
    NoAxes() +
    guides(color = guide_legend(override.aes = list(size = 3)))

  ggsave(file.path(OUTPUT_DIR, "umap_condition.pdf"), p1, width = 6, height = 5)
  ggsave(file.path(OUTPUT_DIR, "umap_celltype_labeled.pdf"), p2, width = 8, height = 6)
  ggsave(file.path(OUTPUT_DIR, "umap_celltype_legend.pdf"), p2_legend, width = 9, height = 6)
  ggsave(file.path(OUTPUT_DIR, "umap_sample.pdf"), p3, width = 7, height = 5)

  # Combined panel (2x2 grid, equal sizes)
  p_combined <- (p1 | p3) / (p2_legend | plot_spacer()) +
    plot_annotation(title = "4T1 Spatial Profiling",
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
  ggsave(file.path(OUTPUT_DIR, "umap_overview_panel.pdf"), p_combined, width = 14, height = 12)
}

# Spatial plot of cell types per sample (one row)
sample_ids <- levels(sobj$sample)
spatial_plots <- lapply(sample_ids, function(s) {
  sub <- subset(sobj, subset = sample == s)
  DimPlot(sub, reduction = "spatial", group.by = "Tier1_celltype", pt.size = 1,
          cols = TIER1_COLORS, shuffle = TRUE,
          raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle(s) +
    NoAxes() +
    NoLegend() +
    coord_fixed()
})

# Shared legend from full object
p_legend <- DimPlot(sobj, reduction = "spatial", group.by = "Tier1_celltype",
                    cols = TIER1_COLORS, pt.size = 0) +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
  theme_void()
legend_grob <- cowplot::get_legend(p_legend)

p_spatial_ct <- wrap_plots(spatial_plots, nrow = 1) +
  plot_annotation(title = "Spatial — Cell Types",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_celltype_all_samples.pdf"),
       p_spatial_ct, width = 4 * length(sample_ids), height = 5)

# Version with legend
p_spatial_ct_legend <- wrap_plots(c(spatial_plots, list(legend_grob)), nrow = 1,
                                  widths = c(rep(1, length(sample_ids)), 0.4))
ggsave(file.path(OUTPUT_DIR, "spatial_celltype_all_samples_legend.pdf"),
       p_spatial_ct_legend, width = 4 * length(sample_ids) + 2, height = 5)

# Cell type proportions by condition (stacked bar)
prop_df <- as.data.frame(prop.table(table(sobj$condition, sobj$Tier1_celltype), margin = 1))
colnames(prop_df) <- c("condition", "celltype", "proportion")
prop_df$celltype <- factor(prop_df$celltype, levels = TIER1_ORDER)

p4 <- ggplot(prop_df, aes(x = condition, y = proportion, fill = celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = TIER1_COLORS, name = "Cell Type", breaks = TIER1_ORDER) +
  labs(x = NULL, y = "Proportion") +
  theme_publication() +
  theme(legend.key.size = unit(0.4, "cm"))

ggsave(file.path(OUTPUT_DIR, "celltype_proportions.pdf"), p4, width = 5, height = 6)

# --- 5. Summary statistics tables -----------------------------------------------

# Cells per cell type (overall and by condition)
celltype_counts <- as.data.frame.matrix(table(sobj$Tier1_celltype, sobj$condition))
celltype_counts$Total <- rowSums(celltype_counts)
celltype_counts$Cell_Type <- rownames(celltype_counts)
celltype_counts <- celltype_counts[, c("Cell_Type", names(CONDITION_COLORS), "Total")]
celltype_counts <- celltype_counts[match(levels(sobj$Tier1_celltype), celltype_counts$Cell_Type), ]
write.csv(celltype_counts, file.path(OUTPUT_DIR, "table_cells_per_celltype.csv"), row.names = FALSE)

# Cells per sample (overall and by cell type)
sample_counts <- as.data.frame.matrix(table(sobj$sample, sobj$Tier1_celltype))
sample_counts$Total <- rowSums(sample_counts)
sample_counts$Sample <- rownames(sample_counts)
sample_counts$Condition <- ifelse(grepl("^WT", sample_counts$Sample), "WT", "DKO")
sample_counts <- sample_counts[, c("Sample", "Condition", levels(sobj$Tier1_celltype), "Total")]
write.csv(sample_counts, file.path(OUTPUT_DIR, "table_cells_per_sample.csv"), row.names = FALSE)

# QC summary per sample (median counts, median genes)
qc_summary <- aggregate(
  cbind(total_counts, n_genes_by_counts) ~ sample + condition,
  data = sobj@meta.data,
  FUN = function(x) round(median(x, na.rm = TRUE), 1)
)
colnames(qc_summary) <- c("Sample", "Condition", "Median_Counts", "Median_Genes")
qc_summary$N_Cells <- as.numeric(table(sobj$sample)[qc_summary$Sample])
qc_summary <- qc_summary[, c("Sample", "Condition", "N_Cells", "Median_Counts", "Median_Genes")]
write.csv(qc_summary, file.path(OUTPUT_DIR, "table_qc_per_sample.csv"), row.names = FALSE)

cat("\nSummary tables saved:\n")
cat("  - table_cells_per_celltype.csv\n")
cat("  - table_cells_per_sample.csv\n")
cat("  - table_qc_per_sample.csv\n")

# --- 6. Save Seurat object ------------------------------------------------------
saveRDS(sobj, file.path(PROJECT_DIR, "output", "R", "DKO_Tier1_seurat.rds"))
cat("\nSeurat object saved to:", file.path(PROJECT_DIR, "output", "R", "DKO_Tier1_seurat.rds"), "\n")
cat("Plots/tables saved to:", OUTPUT_DIR, "\n")