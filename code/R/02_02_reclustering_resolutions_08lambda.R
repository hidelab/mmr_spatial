# 02_02_reclustering_resolutions_08lambda.R
# Re-run FindClusters at multiple resolutions on the Banksy+Harmony object
# from 02_01 to obtain coarser/finer clusterings for domain identification.

library(Seurat)
library(ggplot2)
library(patchwork)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- 1. Load the Banksy+Harmony Seurat object -----------------------------------
BANKSY_RDS <- file.path(PROJECT_DIR, "output", "R", "DKO_banksy_seurat.rds")
seu <- readRDS(BANKSY_RDS)
cat("Loaded Seurat object:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# --- 2. Re-run FindClusters at multiple resolutions --------------------------------
# The SNN graph from FindNeighbors (pc30) is already stored in the object.
resolutions <- BANKSY_RESOLUTIONS_SWEEP

for (res in resolutions) {
  col_name <- paste0("banksy_res_", res)
  seu <- FindClusters(seu, resolution = res, cluster.name = col_name)
  n_clust <- length(unique(seu@meta.data[[col_name]]))
  cat(sprintf("Resolution %.2f -> %d clusters\n", res, n_clust))
}

# Keep a reference to the original primary resolution column
orig_col <- paste0("banksy_clust_pc", BANKSY_NPC_USE, "_res.", BANKSY_RESOLUTION)
if (!orig_col %in% colnames(seu@meta.data)) {
  orig_col <- grep(paste0("banksy_clust_pc", BANKSY_NPC_USE), colnames(seu@meta.data), value = TRUE)[1]
}
seu[[paste0("banksy_res_", BANKSY_RESOLUTION)]] <- seu@meta.data[[orig_col]]

# --- 3. Summary table of cluster counts per resolution ---------------------------
all_res <- c(resolutions, BANKSY_RESOLUTION)
res_summary <- data.frame(
  resolution = all_res,
  n_clusters = sapply(paste0("banksy_res_", all_res), function(col) {
    length(unique(seu@meta.data[[col]]))
  })
)
print(res_summary)
write.csv(res_summary, file.path(OUTPUT_DIR, "resolution_sweep_summary.csv"),
          row.names = FALSE)

# --- 4. UMAP plots for each resolution ------------------------------------------
plots_umap <- lapply(all_res, function(res) {
  col_name <- paste0("banksy_res_", res)
  n_cl <- length(unique(seu@meta.data[[col_name]]))
  DimPlot(seu, reduction = "umap_banksy", group.by = col_name,
          pt.size = 1, shuffle = TRUE,
          raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle(sprintf("res=%.2f (%d clusters)", res, n_cl)) +
    NoAxes()
})

p_res_umap <- wrap_plots(plots_umap, ncol = 2) +
  plot_annotation(title = "Banksy Resolution Sweep",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUTPUT_DIR, "umap_resolution_sweep.pdf"), p_res_umap, width = 14, height = 12)
cat("UMAP resolution sweep saved.\n")

# --- 5. Spatial plots for each resolution ----------------------------------------
sample_ids <- levels(seu$sample)

for (res in all_res) {
  col_name <- paste0("banksy_res_", res)

  plots_spatial <- lapply(sample_ids, function(s) {
    sub <- subset(seu, subset = sample == s)
    coords_s <- Embeddings(sub, "spatial")
    df <- data.frame(x = coords_s[, 1], y = coords_s[, 2],
                     cluster = sub@meta.data[[col_name]])
    ggplot(df, aes(x, y, color = cluster)) +
      geom_point(size = 0.3) +
      coord_fixed() +
      ggtitle(s) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")
  })

  p_sp <- wrap_plots(plots_spatial, nrow = 1) +
    plot_annotation(title = sprintf("Banksy res=%.2f", res),
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
  ggsave(file.path(OUTPUT_DIR, sprintf("spatial_banksy_res%.2f.pdf", res)),
         p_sp, width = 4 * length(sample_ids), height = 5)
}
cat("Spatial cluster plots saved.\n")

# --- 6. Cluster composition: clusters vs cell type heatmaps ---------------------
for (res in all_res) {
  col_name <- paste0("banksy_res_", res)
  ct <- table(seu@meta.data[[col_name]], seu$Tier1_celltype)
  ct_prop <- prop.table(ct, margin = 1)

  pdf(file.path(OUTPUT_DIR, sprintf("heatmap_res%.2f_vs_celltype.pdf", res)),
      width = 10, height = 6)
  heatmap(ct_prop, scale = "none", margins = c(12, 8),
          main = sprintf("Banksy res=%.2f vs cell types (row-normalized)", res),
          col = colorRampPalette(c("white", "navy"))(50))
  dev.off()
}
cat("Cell type composition heatmaps saved.\n")

# --- 7. Cluster proportions by condition for each resolution ---------------------
plots_cond <- lapply(all_res, function(res) {
  col_name <- paste0("banksy_res_", res)
  df <- as.data.frame(prop.table(table(seu@meta.data[[col_name]], seu$condition),
                                 margin = 2))
  colnames(df) <- c("Cluster", "Condition", "Proportion")
  ggplot(df, aes(x = Condition, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
    scale_x_discrete(labels = levels(seu$condition)) +
    ggtitle(sprintf("res=%.2f", res)) +
    labs(x = NULL, y = "Proportion") +
    theme_publication()
})

p_cond <- wrap_plots(plots_cond, ncol = 2) +
  plot_annotation(title = "Cluster Proportions by Condition",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUTPUT_DIR, "barplot_resolution_sweep_by_condition.pdf"),
       p_cond, width = 12, height = 10)
cat("Condition composition plots saved.\n")

# --- 8. Save updated object ------------------------------------------------------
saveRDS(seu, BANKSY_RDS)
cat("Updated Seurat object saved with new cluster columns.\n")

# Export all cluster assignments
all_cluster_cols <- grep("^banksy_", colnames(seu@meta.data), value = TRUE)
cluster_df <- seu@meta.data[, c("sample", "condition", "Tier1_celltype", all_cluster_cols)]
write.csv(cluster_df, file.path(OUTPUT_DIR, "banksy_all_cluster_assignments.csv"))
cat("All cluster assignments exported.\n")
cat("Done. Outputs saved to:", OUTPUT_DIR, "\n")