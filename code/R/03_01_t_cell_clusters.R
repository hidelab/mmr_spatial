# 03_01_t_cell_clusters.R
# Sub-clustering of CD8+ T cells on the expression-only (single-cell) Harmony
# embedding. The "harmony" reduction made in 02_01 is built on the BANKSY
# (spatially augmented) PCA and is NOT used here. The scanpy-derived
# "X_pca_harmony" (RNA assay) is used when present; otherwise PCA + Harmony are
# recomputed on the RNA assay of the CD8+ T subset.
#
# Read-only on existing RDS files: the merged object is loaded but never saved.
# New outputs go to output/R/03_t_cell_clusters/ and existing RDS files there
# are never overwritten.

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette and helpers
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

# Script-specific parameters
TCELL_TYPE        <- "CD8+ T"
EXPR_ASSAY        <- "RNA"
SC_HARMONY_RED    <- "X_pca_harmony"   # expression-only Harmony from the h5ad
TCELL_NPC_USE     <- 30
TCELL_RESOLUTIONS <- c(0.2, 0.4,0.5, 0.6, 0.8, 1.0)
TCELL_RESOLUTION  <- 0.5               # primary resolution for markers/plots
TCELL_SEED        <- 42

# Curated markers, plotted in this exact order
TCELL_MARKERS <- c(
  "Ptprc", "Thy1", "Cd3e", "Cd8a", "Cd8b1", "Tnfrsf9", "Lck", "Nkg7", "Gzmb",
  "Prf1", "Ifng", "Tnfsf10", "Fasl", "Mki67", "Nfatc1",
  "Lag3", "Tigit", "Ctla4", "Cxcr6", "Cxcr3", "Ccr5",
  "Il2ra", "Il2rb", "Il2rg","Klrk1",
  "Cd226", "Cd244a"
)

# Bubble plots: size legend levels (% expressing) and color = mean expression
# min-max scaled per gene across clusters (0 = lowest, 1 = highest), on a
# single-hue scale as in Seurat's DotPlot
BUBBLE_SIZE_BREAKS <- c(15, 30, 45, 60, 75)

bubble_data <- function(obj, genes, group_col) {
  d <- DotPlot(obj, features = genes, group.by = group_col)$data
  d$expr_scaled <- ave(d$avg.exp, d$features.plot, FUN = function(x) {
    r <- range(x)
    if (diff(r) == 0) x * 0 else (x - r[1]) / diff(r)
  })
  d
}

bubble_scales <- function(max_size = 6) list(
  scale_size_area(max_size = max_size, limits = c(0, 100), breaks = BUBBLE_SIZE_BREAKS,
                  name = "% expressing"),
  scale_color_gradient(low = "lightgrey", high = "blue", limits = c(0, 1),
                       name = "Scaled\nmean expr.")
)

# Genes on x (in the given order), clusters on y
bubble_plot_x <- function(obj, genes, group_col) {
  bubble_plot_df(bubble_data(obj, genes, group_col), genes)
}

bubble_plot_df <- function(d, genes) {
  d$features.plot <- factor(d$features.plot, levels = genes)
  ggplot(d, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = expr_scaled)) +
    bubble_scales() +
    labs(x = NULL, y = "Cluster") + theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major = element_line(color = "grey92"))
}

# Save an RDS only if the target does not exist yet (never overwrite RDS data)
save_rds_new <- function(object, path) {
  if (file.exists(path)) {
    warning("RDS already exists, not overwritten: ", path)
    return(invisible(FALSE))
  }
  saveRDS(object, path)
  cat("Saved:", path, "\n")
  invisible(TRUE)
}

# --- 1. Load merged Banksy object (read-only) ---------------------------------
BANKSY_MERGED_RDS <- file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds")
seu <- readRDS(BANKSY_MERGED_RDS)
cat("Loaded:", ncol(seu), "cells\n")
cat("Reductions:", paste(Reductions(seu), collapse = ", "), "\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))
seu$banksy_domain_name <- factor(seu$banksy_domain_name, levels = DOMAIN_ORDER)

# Background coordinates of all cells for the spatial maps (section 8)
all_spatial <- Embeddings(seu, "spatial")
bg_df <- data.frame(x = all_spatial[, 1], y = all_spatial[, 2], sample = seu$sample)
max_x_span <- max(tapply(all_spatial[, 1], seu$sample, function(x) diff(range(x))))
max_y_span <- max(tapply(all_spatial[, 2], seu$sample, function(x) diff(range(x))))
half_span <- max(max_x_span, max_y_span) / 2 * 1.05

# --- 2. Pick the expression-only Harmony embedding ----------------------------
# Accept a Harmony reduction only if it was computed on the RNA assay, never on
# the BANKSY assay.
harmony_reds <- grep("harmony", Reductions(seu), value = TRUE, ignore.case = TRUE)
sc_harmony_reds <- harmony_reds[sapply(harmony_reds, function(r) DefaultAssay(seu[[r]]) == EXPR_ASSAY)]
if (SC_HARMONY_RED %in% sc_harmony_reds) {
  harmony_red <- SC_HARMONY_RED
} else if (length(sc_harmony_reds) > 0) {
  harmony_red <- sc_harmony_reds[1]
} else {
  harmony_red <- NULL
}
cat("Single-cell Harmony reduction:",
    if (is.null(harmony_red)) "none found, will recompute on the subset" else harmony_red, "\n")

# --- 3. Subset CD8+ T cells with cell-level QC --------------------------------
keep <- seu$Tier1_celltype == TCELL_TYPE &
        seu$total_counts >= MIN_TOTAL_COUNTS &
        seu$n_genes_by_counts >= MIN_GENES_DETECTED
cat(sprintf("%s cells: %d total, %d after QC (total_counts >= %d, n_genes >= %d)\n",
            TCELL_TYPE, sum(seu$Tier1_celltype == TCELL_TYPE), sum(keep),
            MIN_TOTAL_COUNTS, MIN_GENES_DETECTED))

DefaultAssay(seu) <- EXPR_ASSAY
keep_reds <- unique(c("spatial", harmony_red))
tcell <- DietSeurat(subset(seu, cells = colnames(seu)[keep]),
                    assays = EXPR_ASSAY, dimreducs = keep_reds)
rm(seu); invisible(gc())

print(table(tcell$sample))

# Normalize only if not already done
if (identical(GetAssayData(tcell, assay = EXPR_ASSAY, layer = "data"),
              GetAssayData(tcell, assay = EXPR_ASSAY, layer = "counts"))) {
  tcell <- NormalizeData(tcell, assay = EXPR_ASSAY)
  cat("Normalized expression data.\n")
} else {
  cat("Data already normalized, skipping.\n")
}

# Fallback: RNA-only PCA + Harmony on the CD8+ T subset (all panel genes)
if (is.null(harmony_red)) {
  tcell <- ScaleData(tcell, features = rownames(tcell))
  tcell <- RunPCA(tcell, features = rownames(tcell), npcs = 50,
                  reduction.name = "pca_rna", verbose = FALSE)
  set.seed(TCELL_SEED)
  tcell <- RunHarmony(tcell, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca_rna",
                      reduction.save = "harmony_rna")
  harmony_red <- "harmony_rna"
}

npc_use <- min(TCELL_NPC_USE, ncol(Embeddings(tcell, harmony_red)))

# --- 4. Clustering over a resolution sweep ------------------------------------
set.seed(TCELL_SEED)
tcell <- FindNeighbors(tcell, reduction = harmony_red, dims = 1:npc_use,
                       graph.name = c("tcell_nn", "tcell_snn"))
for (res in TCELL_RESOLUTIONS) {
  tcell <- FindClusters(tcell, graph.name = "tcell_snn", resolution = res,
                        random.seed = TCELL_SEED,
                        cluster.name = paste0("tcell_clust_res.", res))
}
clust_cols <- paste0("tcell_clust_res.", TCELL_RESOLUTIONS)
primary_col <- paste0("tcell_clust_res.", TCELL_RESOLUTION)
Idents(tcell) <- primary_col

cat("Clusters per resolution:\n")
print(sapply(clust_cols, function(cc) length(unique(tcell@meta.data[[cc]]))))

tcell <- RunUMAP(tcell, reduction = harmony_red, dims = 1:npc_use,
                 reduction.name = "umap_tcell", seed.use = TCELL_SEED, verbose = FALSE)

n_clust <- length(levels(tcell@meta.data[[primary_col]]))
CLUSTER_COLORS <- tcell_subtype_colors(levels(tcell@meta.data[[primary_col]]))

# --- 5. UMAPs -----------------------------------------------------------------
p_res <- lapply(clust_cols, function(cc) {
  DimPlot(tcell, reduction = "umap_tcell", group.by = cc, label = TRUE,
          pt.size = 0.5, shuffle = TRUE, raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle(sub("tcell_clust_", "", cc)) + NoAxes() + NoLegend()
})
save_pdf_png(wrap_plots(p_res, nrow = 1), file.path(OUTPUT_DIR, "umap_resolution_sweep"),
             width = 4 * length(p_res), height = 4.5)

p_clust <- DimPlot(tcell, reduction = "umap_tcell", group.by = primary_col, label = TRUE,
                   cols = CLUSTER_COLORS, pt.size = 0.5, shuffle = TRUE,
                   raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle(paste0(TCELL_TYPE, " clusters (res ", TCELL_RESOLUTION, ")")) + NoAxes()
p_cond <- DimPlot(tcell, reduction = "umap_tcell", group.by = "condition",
                  cols = CONDITION_COLORS, pt.size = 0.5, shuffle = TRUE,
                  raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Condition") + NoAxes()
p_samp <- DimPlot(tcell, reduction = "umap_tcell", group.by = "sample",
                  cols = SAMPLE_COLORS, pt.size = 0.5, shuffle = TRUE,
                  raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Sample") + NoAxes()
p_dom <- DimPlot(tcell, reduction = "umap_tcell", group.by = "banksy_domain_name",
                 cols = DOMAIN_COLORS, pt.size = 0.5, shuffle = TRUE,
                 raster = TRUE, raster.dpi = c(300, 300)) +
  ggtitle("Spatial domain") + NoAxes()
save_pdf_png((p_clust | p_cond) / (p_samp | p_dom), file.path(OUTPUT_DIR, "umap_clusters_overview"),
             width = 12, height = 10)

p_feat <- FeaturePlot(tcell, reduction = "umap_tcell",
                      features = intersect(TCELL_MARKERS, rownames(tcell)),
                      order = TRUE, pt.size = 0.3, raster = TRUE, ncol = 6) &
  NoAxes() & theme(plot.title = element_text(size = 10))
save_pdf_png(p_feat, file.path(OUTPUT_DIR, "umap_marker_features"), width = 18, height = 22)

# --- 6. Cluster markers -------------------------------------------------------
markers <- FindAllMarkers(tcell, assay = EXPR_ASSAY, only.pos = TRUE,
                          min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
write.csv(markers, file.path(OUTPUT_DIR, paste0("cluster_markers_", primary_col, ".csv")),
          row.names = FALSE)

markers_present <- intersect(TCELL_MARKERS, rownames(tcell))
# Curated markers (1) per cluster on all cells and (2) per cluster x condition,
# rows "CD8T Cx|WT" and "CD8T Cx|DKO" next to each other. In (2) the per-gene
# color scale is min-max over all cluster x condition rows.
label_clusters <- scale_y_discrete(labels = function(x) paste0("CD8T C", x))
curated_w <- 0.3 * length(markers_present) + 3

p_dot <- bubble_plot_x(tcell, markers_present, primary_col) + label_clusters +
  labs(y = NULL, title = paste0(TCELL_TYPE, ": all cells"))
save_pdf_png(p_dot, file.path(OUTPUT_DIR, "dotplot_curated_markers"),
             width = curated_w, height = 0.35 * n_clust + 2.5)

clust_levels <- levels(tcell@meta.data[[primary_col]])
cond_levels  <- levels(tcell$condition)
tcell$tmp_clust_cond <- factor(
  paste(tcell@meta.data[[primary_col]], tcell$condition, sep = "|"),
  levels = as.vector(t(outer(clust_levels, cond_levels, paste, sep = "|"))))
p_dot_cond <- bubble_plot_x(tcell, markers_present, "tmp_clust_cond") + label_clusters +
  labs(y = NULL, title = paste0(TCELL_TYPE, ": clusters by condition"))
tcell$tmp_clust_cond <- NULL
save_pdf_png(p_dot_cond, file.path(OUTPUT_DIR, "dotplot_curated_markers_by_condition"),
             width = curated_w, height = 0.35 * n_clust * length(cond_levels) + 2.5)

top_markers <- do.call(rbind, lapply(split(markers, markers$cluster), function(df) {
  df <- df[df$p_val_adj < 0.05, ]
  head(df[order(-df$avg_log2FC), ], 5)
}))
top_genes <- unique(top_markers$gene)
if (length(top_genes) > 0) {
  p_top <- bubble_plot_x(tcell, top_genes, primary_col)
  save_pdf_png(p_top, file.path(OUTPUT_DIR, "dotplot_top5_markers"),
               width = 0.3 * length(top_genes) + 3, height = 0.35 * n_clust + 2.5)
}

# Bubble plot of all significant markers (FDR < 0.05). Each gene is grouped
# under the cluster where it has the largest log2FC; within a group, genes are
# ordered by log2FC (strongest at the top).
sig_markers <- markers[markers$p_val_adj < 0.05, ]
if (nrow(sig_markers) > 0) {
  gene_group <- sig_markers[order(-sig_markers$avg_log2FC), ]
  gene_group <- gene_group[!duplicated(gene_group$gene), ]
  gene_group$cluster <- factor(gene_group$cluster, levels = levels(tcell@meta.data[[primary_col]]))
  gene_group <- gene_group[order(gene_group$cluster, -gene_group$avg_log2FC), ]
  write.csv(gene_group[, c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")],
            file.path(OUTPUT_DIR, "bubble_all_sig_markers_genes.csv"), row.names = FALSE)

  dot_df <- bubble_data(tcell, unique(gene_group$gene), primary_col)
  dot_df$group <- paste("Cluster", gene_group$cluster[match(dot_df$features.plot, gene_group$gene)])
  dot_df$features.plot <- factor(dot_df$features.plot, levels = rev(gene_group$gene))

  p_bubble <- ggplot(dot_df, aes(x = id, y = features.plot)) +
    geom_point(aes(size = pct.exp, color = expr_scaled)) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    bubble_scales(max_size = 5) +
    labs(x = "Cluster", y = NULL,
         title = paste0(TCELL_TYPE, " cluster markers (FDR < 0.05)"),
         subtitle = sprintf("%d genes; grouped by cluster with largest log2FC", nrow(gene_group))) +
    theme_publication(base_size = 9) +
    theme(strip.text.y = element_text(angle = 0), panel.grid.major = element_line(color = "grey92"))
  save_pdf_png(p_bubble, file.path(OUTPUT_DIR, "bubble_all_sig_markers"),
               width = 5.5, height = 0.13 * nrow(gene_group) + 2)
}

# Mean scaled expression of curated markers per cluster
avg <- AverageExpression(tcell, assays = EXPR_ASSAY, features = markers_present,
                         group.by = primary_col, layer = "data")[[EXPR_ASSAY]]
avg_z <- t(scale(t(as.matrix(avg))))
colnames(avg_z) <- levels(tcell@meta.data[[primary_col]])
pdf(file.path(OUTPUT_DIR, "heatmap_curated_markers.pdf"), width = 0.4 * n_clust + 3, height = 9)
draw(Heatmap(avg_z, name = "Row z-score",
             col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#FFFFFF", "#B2182B")),
             cluster_rows = FALSE, cluster_columns = TRUE,   # rows keep TCELL_MARKERS order
             clustering_method_columns = "ward.D2",
             row_names_side = "left", border = TRUE,
             column_title = paste0(TCELL_TYPE, " clusters: mean expression (res ", TCELL_RESOLUTION, ")")))
dev.off()

# --- 7. Cluster composition by sample, condition and domain -------------------
meta <- tcell@meta.data
comp_sample <- as.data.frame(prop.table(table(meta[[primary_col]], meta$sample), margin = 2))
colnames(comp_sample) <- c("cluster", "sample", "proportion")
comp_sample$condition <- ifelse(comp_sample$sample %in% WT_SAMPLES, REFERENCE_CONDITION, TEST_CONDITION)
write.csv(comp_sample, file.path(OUTPUT_DIR, "cluster_proportions_by_sample.csv"), row.names = FALSE)

p_bar <- ggplot(comp_sample, aes(x = sample, y = proportion, fill = cluster)) +
  geom_col(color = "white", linewidth = 0.2) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  labs(x = NULL, y = paste("Proportion of", TCELL_TYPE), fill = "Cluster") +
  theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_box <- ggplot(comp_sample, aes(x = condition, y = proportion, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5) +
  facet_wrap(~ cluster, scales = "free_y", nrow = 1) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = paste("Proportion of", TCELL_TYPE)) +
  theme_publication() + theme(legend.position = "none")
save_pdf_png(p_bar, file.path(OUTPUT_DIR, "barplot_cluster_by_sample"), width = 6, height = 5)
save_pdf_png(p_box, file.path(OUTPUT_DIR, "boxplot_cluster_by_condition"),
             width = 1.6 * n_clust + 1, height = 3.5)

comp_domain <- prop.table(table(meta[[primary_col]], meta$banksy_domain_name), margin = 1)
write.csv(as.data.frame.matrix(comp_domain),
          file.path(OUTPUT_DIR, "cluster_distribution_across_domains.csv"))
pdf(file.path(OUTPUT_DIR, "heatmap_cluster_by_domain.pdf"), width = 5, height = 0.35 * n_clust + 2)
draw(Heatmap(as.matrix(unclass(comp_domain)), name = "Proportion",
             col = colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#6BAED6", "#08519C")),
             cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,
             row_names_side = "left", column_names_rot = 45,
             cell_fun = function(j, i, x, y, w, h, fill) {
               grid.text(sprintf("%.2f", comp_domain[i, j]), x, y, gp = gpar(fontsize = 8))
             },
             column_title = "Domain distribution per cluster"))
dev.off()

# --- 8. Spatial maps of clusters per sample -----------------------------------
tc_spatial <- Embeddings(tcell, "spatial")
tc_df <- data.frame(x = tc_spatial[, 1], y = tc_spatial[, 2], sample = tcell$sample,
                    cluster = tcell@meta.data[[primary_col]])
plots_spatial <- lapply(levels(tcell$sample), function(s) {
  p <- ggplot() +
    geom_point(data = bg_df[bg_df$sample == s, ], aes(x, y), color = "grey90", size = 0.05) +
    geom_point(data = tc_df[tc_df$sample == s, ], aes(x, y, color = cluster), size = 0.4) +
    scale_color_manual(values = CLUSTER_COLORS, drop = FALSE) +
    ggtitle(s) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (CENTER_COORDS) {
    p + coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span))
  } else {
    p + coord_fixed()
  }
})
p_spatial <- wrap_plots(plots_spatial, nrow = 1, guides = "collect") &
  guides(color = guide_legend(override.aes = list(size = 3)))
ggsave(file.path(OUTPUT_DIR, "spatial_tcell_clusters.png"), p_spatial,
       width = 4 * length(plots_spatial) + 1, height = 5, dpi = 300)

# --- 9. Export ----------------------------------------------------------------
assign_df <- data.frame(cell = colnames(tcell),
                        meta[, c("sample", "condition", "banksy_domain_name", clust_cols)],
                        Embeddings(tcell, "umap_tcell"), check.names = FALSE)
write.csv(assign_df, file.path(OUTPUT_DIR, "tcell_cluster_assignments.csv"), row.names = FALSE)

save_rds_new(tcell, file.path(OUTPUT_DIR, "DKO_CD8T_clusters.rds"))

cat(sprintf("\nEmbedding: %s, dims 1:%d; primary clustering: %s (%d clusters)\n",
            harmony_red, npc_use, primary_col, n_clust))
cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")
