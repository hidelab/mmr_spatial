# 02_03_merge_domains_08lambda.R
# Hierarchical merging of Banksy spatial domains.
# Computes cluster centroids in Harmony embedding space, builds a dendrogram,
# visualizes cell-type composition per domain, and merges domains by cutting
# the tree at a user-defined k.

library(Seurat)
library(ggplot2)
library(patchwork)
library(pheatmap)

PROJECT_DIR <- "/data/work/Pourya/mmr_spatial"
OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "merge_domains")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# USER PARAMETERS
# =============================================================================
DOMAIN_COL <- "banksy_res_0.8"
MERGE_K <- 4

# Domain color palette
DOMAIN_COLORS <- c(
  "1" = "#F9A825",   # Mammary Glands — yellow
  "2" = "#B71C1C",   # Tumor Core — deep red
  "3" = "#4A148C",   # Immune Engulfing — deep purple
  "4" = "#BDBDBD"    # Stroma — light gray
)

# Domain names (assigned after inspecting composition)
DOMAIN_NAMES <- c(
  "1" = "Mammary Glands",
  "2" = "Tumor Core",
  "3" = "Immune Engulfing",
  "4" = "Stroma"
)

# Display order for domains: Tumor Core, Immune Engulfing, Stroma, Mammary Glands
DOMAIN_ORDER <- c("2", "3", "4", "1")

# --- 1. Load Banksy-clustered Seurat object -----------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_seurat.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", "DKO", seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

domains <- sort(unique(seu@meta.data[[DOMAIN_COL]]))
cat(sprintf("Original domains: %d (labels: %s)\n", length(domains),
            paste(domains, collapse = ", ")))

# --- 2. Compute cluster centroids in Harmony space ----------------------------
harmony_emb <- Embeddings(seu, "harmony")
domain_labels <- seu@meta.data[[DOMAIN_COL]]

centroids <- t(sapply(domains, function(d) {
  colMeans(harmony_emb[domain_labels == d, , drop = FALSE])
}))
rownames(centroids) <- domains
cat(sprintf("Centroids matrix: %d domains x %d dims\n", nrow(centroids), ncol(centroids)))

# --- 3. Hierarchical clustering on centroids ----------------------------------
dist_mat <- dist(centroids, method = "euclidean")
hc <- hclust(dist_mat, method = "ward.D2")

# --- 4. Plot dendrogram -------------------------------------------------------
pdf(file.path(OUTPUT_DIR, "domain_dendrogram.pdf"), width = 10, height = 6)
plot(hc, main = "Hierarchical clustering of Banksy domains (Harmony centroids)",
     xlab = "Domain", ylab = "Ward distance", sub = "")
cut_height <- sort(hc$height, decreasing = TRUE)[MERGE_K - 1]
abline(h = cut_height, col = "red", lty = 2)
text(x = length(domains) * 0.8, y = cut_height * 1.05,
     labels = paste("k =", MERGE_K), col = "red", cex = 0.9)
dev.off()
cat("Dendrogram saved.\n")

# --- 5. Cell-type composition heatmap per domain ------------------------------
ct_table <- table(domain_labels, seu$Tier1_celltype)
ct_prop <- prop.table(ct_table, margin = 1)

# Ordered by dendrogram
ct_prop_ordered <- ct_prop[hc$labels[hc$order], ]

pdf(file.path(OUTPUT_DIR, "heatmap_domain_celltype_composition.pdf"), width = 12, height = 8)
pheatmap(ct_prop_ordered,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(c("#FFFFFF", "#C8DCF0", "#6BAED6", "#2171B5"))(100),
         main = "Cell-type composition per Banksy domain (dendrogram order)",
         fontsize_row = 9, fontsize_col = 9)
dev.off()

pdf(file.path(OUTPUT_DIR, "heatmap_domain_celltype_clustered.pdf"), width = 12, height = 8)
pheatmap(ct_prop,
         clustering_method = "ward.D2",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",
         color = colorRampPalette(c("#FFFFFF", "#C8DCF0", "#6BAED6", "#2171B5"))(100),
         main = "Cell-type composition per Banksy domain (clustered)",
         fontsize_row = 9, fontsize_col = 9)
dev.off()
cat("Cell-type composition heatmaps saved.\n")

# --- 6. Cut tree and assign merged domain labels ------------------------------
merged_assignments <- cutree(hc, k = MERGE_K)

merged_col <- "banksy_domain_merged"
seu@meta.data[[merged_col]] <- as.character(merged_assignments[as.character(domain_labels)])

cat(sprintf("\nMerged %d domains into %d groups:\n", length(domains), MERGE_K))
merge_map <- data.frame(
  merged_domain = merged_assignments,
  original_domain = names(merged_assignments)
)
merge_map <- merge_map[order(merge_map$merged_domain, merge_map$original_domain), ]
for (m in sort(unique(merged_assignments))) {
  orig <- merge_map$original_domain[merge_map$merged_domain == m]
  cat(sprintf("  Merged %d <- original: %s\n", m, paste(orig, collapse = ", ")))
}

# --- 7. Condition composition per merged domain --------------------------------
ct_cond <- as.data.frame(prop.table(table(seu@meta.data[[merged_col]], seu$condition), margin = 1))
colnames(ct_cond) <- c("Domain", "Condition", "Proportion")
ct_cond$Domain <- factor(ct_cond$Domain, levels = as.character(1:MERGE_K))

p_cond <- ggplot(ct_cond, aes(x = Domain, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = "Merged Domain", y = "Proportion", title = "Condition composition per domain") +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "barplot_domain_condition.pdf"), p_cond, width = 6, height = 5)

# --- 8. Cell-type composition per merged domain --------------------------------
ct_merged <- table(seu@meta.data[[merged_col]], seu$Tier1_celltype)
ct_merged_prop <- prop.table(ct_merged, margin = 1)
ct_merged_prop <- ct_merged_prop[DOMAIN_ORDER, ]

pdf(file.path(OUTPUT_DIR, "heatmap_merged_domain_celltype.pdf"), width = 8, height = 6)
pheatmap(ct_merged_prop,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         clustering_method = "ward.D2",
         scale = "none",
         color = colorRampPalette(c("white", "steelblue", "navy"))(100),
         main = paste("Cell-type composition per merged domain (k =", MERGE_K, ")"),
         fontsize_row = 11, fontsize_col = 9,
         labels_row = DOMAIN_NAMES[DOMAIN_ORDER])
dev.off()

write.csv(as.data.frame.matrix(ct_merged_prop),
          file.path(OUTPUT_DIR, "merged_domain_celltype_proportions.csv"))

# --- 8b. Cell-type distribution across merged domains (inverse perspective) ----
# Shows: what % of each cell type is in each domain
ct_merged_inv <- table(seu$Tier1_celltype, seu@meta.data[[merged_col]])
ct_merged_inv_prop <- prop.table(ct_merged_inv, margin = 1)
ct_merged_inv_prop <- ct_merged_inv_prop[, DOMAIN_ORDER]

pdf(file.path(OUTPUT_DIR, "heatmap_celltype_distribution_across_domains.pdf"), width = 8, height = 8)
pheatmap(ct_merged_inv_prop,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         clustering_method = "ward.D2",
         scale = "none",
         color = colorRampPalette(c("white", "steelblue", "navy"))(100),
         main = paste("Distribution of cell types across domains (k =", MERGE_K, ")"),
         fontsize_row = 11, fontsize_col = 10,
         labels_col = DOMAIN_NAMES[DOMAIN_ORDER],
         display_numbers = TRUE,
         number_format = "%.2f",
         number_color = "black",
         fontsize_number = 9)
dev.off()

write.csv(as.data.frame.matrix(ct_merged_inv_prop),
          file.path(OUTPUT_DIR, "celltype_distribution_across_domains.csv"))

# --- 8c. Cell-type distribution by condition -----------------------------------
# Get row ordering from the full dataset heatmap
temp_heatmap <- pheatmap(ct_merged_inv_prop,
                         cluster_rows = TRUE, cluster_cols = FALSE,
                         clustering_method = "ward.D2", silent = TRUE)
row_order <- rownames(ct_merged_inv_prop)[temp_heatmap$tree_row$order]

all_celltypes <- rownames(ct_merged_inv_prop)
all_domains <- DOMAIN_ORDER

for (cond in levels(seu$condition)) {
  idx_cond <- which(seu$condition == cond)
  celltype_cond <- seu$Tier1_celltype[idx_cond]
  domain_cond <- seu@meta.data[[merged_col]][idx_cond]

  ct_cond_inv <- table(celltype_cond, domain_cond)
  ct_cond_inv_prop <- prop.table(ct_cond_inv, margin = 1)

  # Pad with zeros to include all cell types and domains
  ct_cond_full <- matrix(0, nrow = length(all_celltypes), ncol = length(all_domains),
                         dimnames = list(all_celltypes, all_domains))
  ct_cond_full[rownames(ct_cond_inv_prop), colnames(ct_cond_inv_prop)] <- ct_cond_inv_prop
  ct_cond_full <- ct_cond_full[row_order, DOMAIN_ORDER]

  pdf(file.path(OUTPUT_DIR, paste0("heatmap_celltype_distribution_", cond, ".pdf")),
      width = 8, height = 8)
  pheatmap(ct_cond_full,
           cluster_rows = FALSE, cluster_cols = FALSE,
           scale = "none",
           color = colorRampPalette(c("white", "steelblue", "navy"))(100),
           main = sprintf("Cell type distribution across domains — %s (k = %d)", cond, MERGE_K),
           fontsize_row = 11, fontsize_col = 10,
           labels_col = DOMAIN_NAMES[DOMAIN_ORDER],
           display_numbers = TRUE,
           number_format = "%.2f",
           number_color = "black",
           fontsize_number = 9)
  dev.off()

  write.csv(as.data.frame.matrix(ct_cond_full),
            file.path(OUTPUT_DIR, paste0("celltype_distribution_", cond, ".csv")))
  cat(sprintf("Saved condition-specific plots for %s\n", cond))
}

# --- 9. Spatial plots of merged domains per sample ----------------------------
seu@meta.data[[merged_col]] <- factor(seu@meta.data[[merged_col]],
                                       levels = as.character(1:MERGE_K))
sample_ids <- levels(seu$sample)

plots_spatial <- lapply(sample_ids, function(s) {
  sub <- subset(seu, subset = sample == s)
  DimPlot(sub, reduction = "spatial", group.by = merged_col, pt.size = 1,
          cols = DOMAIN_COLORS, shuffle = TRUE,
          raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle(s) +
    NoAxes() +
    NoLegend() +
    coord_fixed()
})

p_spatial <- wrap_plots(plots_spatial, nrow = 1) +
  plot_annotation(title = paste("Merged Banksy Domains (k =", MERGE_K, ")"),
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUTPUT_DIR, "spatial_merged_domains.pdf"),
       p_spatial, width = 4 * length(sample_ids), height = 5)

# UMAP view of merged domains vs cell type
p_dom_umap <- DimPlot(seu, reduction = "umap_banksy", group.by = merged_col,
                      pt.size = 1, cols = DOMAIN_COLORS, shuffle = TRUE,
                      raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle(paste("Merged Domains (k =", MERGE_K, ")")) +
  NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3)))

p_ct_umap <- DimPlot(seu, reduction = "umap_banksy", group.by = "Tier1_celltype",
                     pt.size = 1, cols = TIER1_COLORS, shuffle = TRUE,
                     raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Cell Type") +
  NoAxes() +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))

ggsave(file.path(OUTPUT_DIR, "umap_merged_domains_vs_celltype.pdf"),
       p_dom_umap + p_ct_umap, width = 16, height = 6)
cat("Spatial and UMAP plots saved.\n")

# --- 10. Assign domain names ---------------------------------------------------
seu@meta.data[["banksy_domain_name"]] <- DOMAIN_NAMES[as.character(seu@meta.data[[merged_col]])]
cat("Domain name assignments:\n")
for (d in names(DOMAIN_NAMES)) {
  cat(sprintf("  Domain %s -> %s\n", d, DOMAIN_NAMES[d]))
}

# --- 11. Export merge mapping and save updated object --------------------------
write.csv(merge_map, file.path(OUTPUT_DIR, "domain_merge_mapping.csv"), row.names = FALSE)

saveRDS(seu@meta.data, file.path(OUTPUT_DIR, "DKO_banksy_merged_metadata.rds"))
saveRDS(seu, file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))

cat("\nSeurat object saved to:", file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"), "\n")
cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")