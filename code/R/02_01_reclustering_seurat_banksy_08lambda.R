# 02_01_reclustering_seurat_banksy_08lambda.R
# Spatial clustering with Banksy (lambda=0.8, strong spatial context) via
# Seurat interface, Harmony batch correction, and cluster visualization.

library(Seurat)
library(Banksy)
library(harmony)
library(ggplot2)
library(patchwork)
library(SeuratWrappers)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- 1. Load Seurat object from previous step ----------------------------------
seu <- readRDS(SEURAT_RDS)
cat("Loaded Seurat object:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# --- 2. Store spatial coordinates in metadata (needed for RunBanksy) -------------
# Coordinates are in the "spatial" DimReduction; RunBanksy needs them as metadata cols
spatial_coords <- Embeddings(seu, "spatial")
seu$spatial_x <- spatial_coords[, 1]
seu$spatial_y <- spatial_coords[, 2]

# --- 3. Identify expression assay and find variable features --------------------
expr_assay <- DefaultAssay(seu)
cat("Expression assay:", expr_assay, "\n")

# Normalize only if not already done
if (is.null(seu[[expr_assay]]@layers[["data"]]) ||
    identical(GetAssayData(seu, assay = expr_assay, layer = "data"),
              GetAssayData(seu, assay = expr_assay, layer = "counts"))) {
  seu <- NormalizeData(seu, assay = expr_assay)
  cat("Normalized expression data.\n")
} else {
  cat("Data already normalized, skipping.\n")
}

seu <- FindVariableFeatures(seu, assay = expr_assay, nfeatures = 2000)

# --- 4. Run Banksy (Seurat interface) -------------------------------------------
# Uses dimx/dimy to pass coordinates directly — no FOV/image object needed.
# lambda=0.2: moderate spatial context; lambda=0.8: strong spatial context
# k_geom=15: local neighbors; k_geom=30: broader neighborhood
lambdas <- c(BANKSY_LAMBDA)
k_geoms <- c(BANKSY_K_GEOM)

seu <- RunBanksy(seu,
                 lambda = lambdas,
                 k_geom = k_geoms,
                 dimx = "spatial_x",
                 dimy = "spatial_y",
                 assay = expr_assay,
                 slot = "data",
                 features = "variable",
                 group = "sample",
                 split.scale = TRUE)

cat("Banksy features computed.\n")
cat("Available assays:", paste(Assays(seu), collapse = ", "), "\n")

# --- 5. PCA on Banksy-augmented assay ------------------------------------------
seu <- RunPCA(seu, assay = "BANKSY", features = rownames(seu[["BANKSY"]]), npcs = BANKSY_NPCS)

# --- 6. Harmony batch correction across samples ---------------------------------
seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca")
cat("Harmony correction complete.\n")

# --- 7. Clustering sweep over nPCs and resolutions ------------------------------
npc_values <- c(BANKSY_NPC_USE)
resolutions <- c(BANKSY_RESOLUTION)

for (npc in npc_values) {
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:npc)
  seu <- FindClusters(seu,
                      resolution = resolutions,
                      cluster.name = paste0("banksy_clust_pc", npc))
}

# List all cluster columns created
cluster_cols <- grep("^banksy_clust_pc", colnames(seu@meta.data), value = TRUE)
cat("Cluster columns created:\n")
print(cluster_cols)

# --- 8. UMAP for visualization --------------------------------------------------
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:BANKSY_NPC_USE,
               reduction.name = "umap_banksy")

# --- 9. Visualization -----------------------------------------------------------

# Pick a primary clustering for plots
primary_clust <- paste0("banksy_clust_pc", BANKSY_NPC_USE, "_res.", BANKSY_RESOLUTION)
if (!primary_clust %in% colnames(seu@meta.data)) {
  primary_clust <- cluster_cols[1]
}
cat("Primary clustering for plots:", primary_clust, "\n")

# UMAP: Banksy clusters vs cell type
p_clust <- DimPlot(seu, reduction = "umap_banksy", group.by = primary_clust,
                   pt.size = 1, shuffle = TRUE,
                   raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Banksy Clusters") +
  NoAxes()

p_ct <- DimPlot(seu, reduction = "umap_banksy", group.by = "Tier1_celltype",
                pt.size = 1, cols = TIER1_COLORS, shuffle = TRUE,
                raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Cell Type") +
  NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

ggsave(file.path(OUTPUT_DIR, "umap_banksy_clusters.pdf"), p_clust, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "umap_banksy_celltype.pdf"), p_ct, width = 9, height = 6)
ggsave(file.path(OUTPUT_DIR, "umap_banksy_clusters_vs_celltype.pdf"),
       p_clust + p_ct, width = 16, height = 6)

# UMAP: condition and sample
p_cond <- DimPlot(seu, reduction = "umap_banksy", group.by = "condition",
                  pt.size = 1, cols = CONDITION_COLORS, shuffle = TRUE,
                  raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Condition") +
  NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_samp <- DimPlot(seu, reduction = "umap_banksy", group.by = "sample",
                  pt.size = 1, cols = SAMPLE_COLORS, shuffle = TRUE,
                  raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Sample") +
  NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(file.path(OUTPUT_DIR, "umap_banksy_condition_sample.pdf"),
       p_cond + p_samp, width = 13, height = 5)

# Spatial plots per sample colored by Banksy cluster
sample_ids <- levels(seu$sample)
plots_spatial <- lapply(sample_ids, function(s) {
  sub <- subset(seu, subset = sample == s)
  coords_s <- Embeddings(sub, "spatial")
  df <- data.frame(x = coords_s[, 1], y = coords_s[, 2],
                   cluster = sub@meta.data[[primary_clust]])
  ggplot(df, aes(x, y, color = cluster)) +
    geom_point(size = 0.3) +
    coord_fixed() +
    ggtitle(s) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
})

p_spatial <- wrap_plots(plots_spatial, ncol = 4) +
  plot_annotation(title = "Banksy Spatial Domains",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUTPUT_DIR, "spatial_banksy_clusters.pdf"), p_spatial, width = 16, height = 10)

# Cluster composition heatmap: Banksy clusters vs Tier1 cell types
ct <- table(seu@meta.data[[primary_clust]], seu$Tier1_celltype)
ct_prop <- prop.table(ct, margin = 1)

pdf(file.path(OUTPUT_DIR, "heatmap_cluster_vs_celltype.pdf"), width = 10, height = 8)
heatmap(ct_prop, scale = "none", margins = c(12, 8),
        main = "Banksy clusters vs cell types (row-normalized)",
        col = colorRampPalette(c("white", "navy"))(50))
dev.off()

# Cluster proportions by condition (ggplot)
ct_cond <- as.data.frame(prop.table(table(seu@meta.data[[primary_clust]], seu$condition), margin = 2))
colnames(ct_cond) <- c("cluster", "condition", "proportion")

p_bar <- ggplot(ct_cond, aes(x = condition, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.2) +
  labs(x = NULL, y = "Proportion", fill = "Cluster") +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "barplot_cluster_by_condition.pdf"), p_bar, width = 5, height = 6)

# --- 10. Save results -----------------------------------------------------------
BANKSY_RDS <- file.path(PROJECT_DIR, "output", "R", "DKO_banksy_seurat.rds")
saveRDS(seu, BANKSY_RDS)
cat("\nSeurat object with Banksy results saved to:", BANKSY_RDS, "\n")

# Export cluster assignments
cluster_df <- seu@meta.data[, c("sample", "condition", "Tier1_celltype", cluster_cols)]
write.csv(cluster_df, file.path(OUTPUT_DIR, "banksy_cluster_assignments.csv"))
cat("Cluster assignments exported to CSV.\n")
cat("Done. Plots saved to:", OUTPUT_DIR, "\n")
