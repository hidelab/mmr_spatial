# 03_11_clustering_redo.R
# Re-cluster all cells from scratch (expression only) to audit the original
# Tier1 cell-type labels. The original labels are a 1:1 mapping of the scanpy
# "leiden" clusters, computed on all cells including very low-count ones.
#
# Pipeline (RNA assay, starting from raw counts; no spatial information):
#   QC filter (MIN_TOTAL_COUNTS, MIN_GENES_DETECTED) -> LogNormalize (scale
#   factor = median counts per cell) -> scale all panel genes -> PCA ->
#   Harmony by sample -> kNN/SNN graph -> Leiden clustering over a resolution
#   sweep -> UMAP.
#
# Label matching: each new cluster gets the original Tier1 label that holds the
# majority of its cells. Clusters whose majority label holds < MIN_PURITY of
# their cells are labelled "Mixed (<A>/<B>)". Agreement with the original labels
# is reported per resolution (ARI, accuracy) and per cell type (precision,
# recall, F1), plus a confusion matrix and the fate of original CD8+ T cells.
#
# Reproducibility: every stochastic step has a fixed seed (PCA, Harmony, Leiden,
# UMAP, marker subsampling) and BLAS/OpenMP threads are capped; parameters and
# sessionInfo() are written to the output folder.
#
# Read-only on existing data: the merged object is loaded but never saved. All
# outputs go to output/R/03_11_clustering_redo/, and existing RDS files there
# are never overwritten.

# --- Thread limits (set BEFORE BLAS/library load) ----------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS,
           NUMEXPR_NUM_THREADS = N_THREADS)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(N_THREADS)
  RhpcBLASctl::omp_set_num_threads(N_THREADS)
}

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(mclust)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "03_11_clustering_redo")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette and helpers
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

# =============================================================================
# Parameters
# =============================================================================
SEED               <- 42
EXPR_ASSAY         <- "RNA"
CLUSTER_QC_ONLY    <- TRUE      # cluster only cells passing MIN_TOTAL_COUNTS / MIN_GENES_DETECTED
REDO_NPCS          <- 50        # PCs computed
REDO_NPC_USE       <- 30        # Harmony dims used for the graph / UMAP
REDO_K_PARAM       <- 20        # kNN size for FindNeighbors
REDO_ALGORITHM     <- 4         # 4 = Leiden (as the original scanpy labels), 1 = Louvain
REDO_RESOLUTIONS   <- c(0.5, 1.0, 1.5, 2.0)
REDO_RESOLUTION    <- 1.0       # primary resolution for labels / markers / plots
MIN_PURITY         <- 0.5       # min majority-label fraction to accept a match
MARKER_MAX_CELLS   <- 1000      # cells sampled per cluster for FindAllMarkers
ORIG_LABEL_COL     <- "Tier1_celltype"
ORIG_CLUSTER_COL   <- "leiden"

CANONICAL_MARKERS <- c(
  "Krt5", "Pgr", "Mki67",                                   # epithelial / tumor
  "Ptprc", "Cd3e", "Cd8a", "Gzmb", "Foxp3", "Nkg7", "Klrb1c", "Ncr1",   # lymphoid
  "Csf1r", "Mrc1", "Arg1", "Cxcl16", "Itgax",                # macrophage
  "Xcr1", "Clec9a", "Siglech", "Ccr7",                       # DC
  "Ly6g", "S100a8", "Cpa3", "Mcpt4", "Kit",                  # granulocyte / mast
  "Pdgfra", "Acta2", "Pecam1", "Cdh5", "Kdr", "Adipoq",      # stroma
  "Cd19", "Cd79a"                                            # B cells
)

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

params <- list(SEED = SEED, EXPR_ASSAY = EXPR_ASSAY, CLUSTER_QC_ONLY = CLUSTER_QC_ONLY,
               MIN_TOTAL_COUNTS = MIN_TOTAL_COUNTS, MIN_GENES_DETECTED = MIN_GENES_DETECTED,
               REDO_NPCS = REDO_NPCS, REDO_NPC_USE = REDO_NPC_USE, REDO_K_PARAM = REDO_K_PARAM,
               REDO_ALGORITHM = REDO_ALGORITHM,
               REDO_RESOLUTIONS = paste(REDO_RESOLUTIONS, collapse = ", "),
               REDO_RESOLUTION = REDO_RESOLUTION, MIN_PURITY = MIN_PURITY,
               HARMONY_GROUP_VAR = HARMONY_GROUP_VAR, N_THREADS = N_THREADS)
writeLines(c(sprintf("03_11 clustering redo, run %s", format(Sys.time())),
             sprintf("%-20s = %s", names(params), unlist(params))),
           file.path(OUTPUT_DIR, "run_parameters.txt"))

# =============================================================================
# 1. Load merged object (read-only) and keep raw counts only
# =============================================================================
seu_full <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu_full), "cells\n")

seu_full$condition <- gsub("^KO$", TEST_CONDITION, seu_full$condition)
seu_full$condition <- factor(seu_full$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu_full$sample <- factor(seu_full$sample,
                          levels = intersect(names(SAMPLE_COLORS), unique(seu_full$sample)))
seu_full$Tier1_celltype <- factor(seu_full$Tier1_celltype,
                                  levels = intersect(TIER1_ORDER, unique(seu_full$Tier1_celltype)))

qc_pass <- seu_full$total_counts >= MIN_TOTAL_COUNTS &
           seu_full$n_genes_by_counts >= MIN_GENES_DETECTED
keep <- if (CLUSTER_QC_ONLY) qc_pass else rep(TRUE, ncol(seu_full))
cat(sprintf("Cells clustered: %d of %d (QC pass: %d)\n", sum(keep), length(keep), sum(qc_pass)))

# Fresh object from raw counts; carry over metadata and spatial coordinates only
meta_keep <- c("sample", "condition", "total_counts", "n_genes_by_counts",
               ORIG_LABEL_COL, "Tier1_group", ORIG_CLUSTER_COL, "banksy_domain_name")
seu <- CreateSeuratObject(
  counts = GetAssayData(seu_full, assay = EXPR_ASSAY, layer = "counts")[, keep],
  assay = EXPR_ASSAY,
  meta.data = seu_full@meta.data[keep, meta_keep]
)
spatial_xy <- Embeddings(seu_full, "spatial")[keep, ]
rm(seu_full); invisible(gc())
seu$Tier1_celltype <- droplevels(seu$Tier1_celltype)
cat("Fresh object:", ncol(seu), "cells x", nrow(seu), "genes\n")

# =============================================================================
# 2. Normalize, scale, PCA, Harmony
# =============================================================================
scale_factor <- median(seu$nCount_RNA)
cat(sprintf("LogNormalize with scale factor = median counts per cell = %g\n", scale_factor))
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = scale_factor)

# Targeted panel: use every gene rather than selecting variable features
VariableFeatures(seu) <- rownames(seu)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = rownames(seu), npcs = REDO_NPCS, seed.use = SEED, verbose = FALSE)

p_elbow <- ElbowPlot(seu, ndims = REDO_NPCS) +
  geom_vline(xintercept = REDO_NPC_USE, linetype = "dashed", color = "red") +
  ggtitle("PCA elbow (red = dims used)") + theme_publication()
ggsave(file.path(OUTPUT_DIR, "pca_elbow.pdf"), p_elbow, width = 6, height = 4)

set.seed(SEED)
seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca",
                  dims.use = 1:REDO_NPCS, reduction.save = "harmony", verbose = FALSE)
cat("Harmony done.\n")

# =============================================================================
# 3. Graph, Leiden clustering sweep, UMAP
# =============================================================================
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:REDO_NPC_USE, k.param = REDO_K_PARAM,
                     graph.name = c("redo_nn", "redo_snn"), verbose = FALSE)

clust_cols <- paste0("redo_clust_res.", REDO_RESOLUTIONS)
for (i in seq_along(REDO_RESOLUTIONS)) {
  t0 <- Sys.time()
  seu <- FindClusters(seu, graph.name = "redo_snn", resolution = REDO_RESOLUTIONS[i],
                      algorithm = REDO_ALGORITHM, random.seed = SEED,
                      cluster.name = clust_cols[i], verbose = FALSE)
  cat(sprintf("  res %.2f: %d clusters (%.1f min)\n", REDO_RESOLUTIONS[i],
              length(unique(seu@meta.data[[clust_cols[i]]])),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
primary_col <- paste0("redo_clust_res.", REDO_RESOLUTION)

seu <- RunUMAP(seu, reduction = "harmony", dims = 1:REDO_NPC_USE, seed.use = SEED,
               reduction.name = "umap_redo", verbose = FALSE)

# =============================================================================
# 4. Match new clusters to the original labels
# =============================================================================
orig <- as.character(seu@meta.data[[ORIG_LABEL_COL]])

match_clusters <- function(clusters, orig) {
  tab <- table(cluster = clusters, label = orig)
  frac <- prop.table(tab, 1)
  top <- apply(frac, 1, function(r) names(sort(r, decreasing = TRUE))[1:2])
  purity <- apply(frac, 1, max)
  data.frame(
    cluster      = rownames(frac),
    n_cells      = as.integer(rowSums(tab)),
    top_label    = top[1, ],
    purity       = as.numeric(purity),
    second_label = top[2, ],
    second_frac  = apply(frac, 1, function(r) sort(r, decreasing = TRUE)[2]),
    matched_label = ifelse(purity >= MIN_PURITY, top[1, ],
                           paste0("Mixed (", top[1, ], "/", top[2, ], ")")),
    stringsAsFactors = FALSE
  )
}

per_type_agreement <- function(new_label, orig) {
  types <- sort(unique(orig))
  data.frame(
    label     = types,
    n_orig    = sapply(types, function(t) sum(orig == t)),
    n_new     = sapply(types, function(t) sum(new_label == t)),
    precision = sapply(types, function(t) mean(orig[new_label == t] == t)),
    recall    = sapply(types, function(t) mean(new_label[orig == t] == t)),
    row.names = NULL
  ) %>% mutate(F1 = 2 * precision * recall / (precision + recall))
}

match_list <- list(); summary_res <- data.frame()
for (i in seq_along(REDO_RESOLUTIONS)) {
  cc <- clust_cols[i]
  cl <- as.character(seu@meta.data[[cc]])
  mm <- match_clusters(cl, orig)
  label_col <- sub("clust", "label", cc)
  seu@meta.data[[label_col]] <- mm$matched_label[match(cl, mm$cluster)]
  match_list[[cc]] <- mm
  write.csv(mm, file.path(OUTPUT_DIR, paste0("cluster_label_match_", cc, ".csv")), row.names = FALSE)

  new_lab <- seu@meta.data[[label_col]]
  summary_res <- rbind(summary_res, data.frame(
    resolution        = REDO_RESOLUTIONS[i],
    n_clusters        = nrow(mm),
    n_mixed_clusters  = sum(mm$purity < MIN_PURITY),
    frac_cells_mixed  = mean(grepl("^Mixed", new_lab)),
    ARI_vs_Tier1      = adjustedRandIndex(cl, orig),
    ARI_vs_leiden     = adjustedRandIndex(cl, seu@meta.data[[ORIG_CLUSTER_COL]]),
    accuracy_vs_Tier1 = mean(new_lab == orig),
    types_without_cluster = paste(setdiff(unique(orig), mm$matched_label), collapse = "; ")
  ))
}
write.csv(summary_res, file.path(OUTPUT_DIR, "agreement_by_resolution.csv"), row.names = FALSE)
cat("\n== Agreement with original labels by resolution ==\n")
print(summary_res, digits = 3)

primary_label_col <- sub("clust", "label", primary_col)
seu$redo_label <- seu@meta.data[[primary_label_col]]
seu$label_changed <- seu$redo_label != orig

agree_type <- per_type_agreement(seu$redo_label, orig)
write.csv(agree_type, file.path(OUTPUT_DIR, paste0("agreement_per_celltype_", primary_col, ".csv")),
          row.names = FALSE)
cat(sprintf("\n== Per cell type agreement (%s) ==\n", primary_col))
print(agree_type, digits = 3)

# Fate of original CD8+ T cells
cd8_fate <- as.data.frame(table(new_label = seu$redo_label[orig == "CD8+ T"])) %>%
  mutate(fraction = Freq / sum(Freq)) %>% arrange(desc(Freq))
write.csv(cd8_fate, file.path(OUTPUT_DIR, paste0("CD8T_original_cells_new_labels_", primary_col, ".csv")),
          row.names = FALSE)
cat("\n== Original CD8+ T cells -> new labels ==\n")
print(head(cd8_fate, 10), digits = 3)

# Relabelled fraction per sample x original type
relabel_tab <- as.data.frame(table(sample = seu$sample, orig = orig, changed = seu$label_changed)) %>%
  group_by(sample, orig) %>% mutate(frac = Freq / sum(Freq)) %>%
  filter(changed == "TRUE") %>% ungroup()
write.csv(relabel_tab, file.path(OUTPUT_DIR, "relabelled_fraction_sample_x_type.csv"), row.names = FALSE)

# =============================================================================
# 5. Markers of the new clusters
# =============================================================================
Idents(seu) <- primary_col
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25,
                          max.cells.per.ident = MARKER_MAX_CELLS, random.seed = SEED,
                          verbose = FALSE)
mm_primary <- match_list[[primary_col]]
markers$matched_label <- mm_primary$matched_label[match(markers$cluster, mm_primary$cluster)]
write.csv(markers, file.path(OUTPUT_DIR, paste0("cluster_markers_", primary_col, ".csv")), row.names = FALSE)

# =============================================================================
# 6. Plots
# =============================================================================
# Cluster order: by matched label (Tier1 order), then by size
label_rank <- match(sub("^Mixed \\((.*)/.*$", "\\1", mm_primary$matched_label), TIER1_ORDER)
mm_primary <- mm_primary[order(label_rank, grepl("^Mixed", mm_primary$matched_label),
                               -mm_primary$n_cells), ]
cluster_order <- mm_primary$cluster
row_labels <- sprintf("%s: %s (%.0f%%)", mm_primary$cluster, mm_primary$matched_label,
                      100 * mm_primary$purity)

# 6a. Confusion heatmap: new clusters (rows) x original labels (columns)
conf <- prop.table(table(seu@meta.data[[primary_col]], orig), 1)[cluster_order, , drop = FALSE]
conf <- conf[, intersect(TIER1_ORDER, colnames(conf)), drop = FALSE]
pdf(file.path(OUTPUT_DIR, paste0("heatmap_cluster_vs_original_", primary_col, ".pdf")),
    width = 11, height = 0.22 * nrow(conf) + 3)
draw(Heatmap(unclass(conf), name = "Fraction of\ncluster",
             col = colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#6BAED6", "#08519C")),
             cluster_rows = FALSE, cluster_columns = FALSE, row_labels = row_labels,
             row_names_gp = gpar(fontsize = 8), column_names_rot = 45, border = TRUE,
             column_title = sprintf("New clusters (%s) vs original Tier1 labels", primary_col)))
dev.off()

# 6b. Cell-level confusion: original label (rows) x new matched label (columns)
conf_lab <- prop.table(table(orig, seu$redo_label), 1)
col_order <- c(intersect(TIER1_ORDER, colnames(conf_lab)),
               sort(setdiff(colnames(conf_lab), TIER1_ORDER)))
conf_lab <- conf_lab[intersect(TIER1_ORDER, rownames(conf_lab)), col_order, drop = FALSE]
pdf(file.path(OUTPUT_DIR, paste0("heatmap_original_vs_new_label_", primary_col, ".pdf")),
    width = 0.4 * ncol(conf_lab) + 4, height = 7)
draw(Heatmap(unclass(conf_lab), name = "Fraction of\noriginal type",
             col = colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#FDAE6B", "#A63603")),
             cluster_rows = FALSE, cluster_columns = FALSE, border = TRUE,
             column_names_rot = 45, row_names_side = "left",
             cell_fun = function(j, i, x, y, w, h, fill) {
               v <- conf_lab[i, j]
               if (v >= 0.05) grid.text(sprintf("%.2f", v), x, y, gp = gpar(fontsize = 7))
             },
             row_title = "Original label", column_title = "New (matched) label"))
dev.off()

# 6c. UMAPs
label_colors <- TIER1_COLORS[intersect(names(TIER1_COLORS), unique(c(orig, seu$redo_label)))]
mixed_labels <- sort(unique(grep("^Mixed", seu$redo_label, value = TRUE)))
label_colors <- c(label_colors, setNames(rep("grey30", length(mixed_labels)), mixed_labels))
umap_args <- list(reduction = "umap_redo", pt.size = 0.3, shuffle = TRUE, seed = SEED,
                  raster = TRUE, raster.dpi = c(400, 400))

p_new <- do.call(DimPlot, c(list(seu, group.by = primary_col, label = TRUE, label.size = 3), umap_args)) +
  NoLegend() + NoAxes() + ggtitle(sprintf("New clusters (%s)", primary_col))
p_lab <- do.call(DimPlot, c(list(seu, group.by = "redo_label", cols = label_colors), umap_args)) +
  NoAxes() + ggtitle("New matched label") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
p_orig <- do.call(DimPlot, c(list(seu, group.by = ORIG_LABEL_COL, cols = label_colors), umap_args)) +
  NoAxes() + ggtitle("Original Tier1 label") +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 1))
p_samp <- do.call(DimPlot, c(list(seu, group.by = "sample", cols = SAMPLE_COLORS), umap_args)) +
  NoAxes() + ggtitle("Sample (Harmony-corrected)")
p_chg <- do.call(DimPlot, c(list(seu, group.by = "label_changed",
                                 cols = c("FALSE" = "grey85", "TRUE" = "#D73027"),
                                 order = "TRUE"), umap_args)) +
  NoAxes() + ggtitle("Label changed")
save_pdf_png((p_new | p_lab) / (p_orig | p_samp), file.path(OUTPUT_DIR, "umap_redo_overview"),
             width = 16, height = 13)
save_pdf_png(p_chg, file.path(OUTPUT_DIR, "umap_label_changed"), width = 7, height = 6)

# 6d. Dot plots: canonical markers and top markers by new cluster
canon <- intersect(CANONICAL_MARKERS, rownames(seu))
seu@meta.data[[primary_col]] <- factor(seu@meta.data[[primary_col]], levels = rev(cluster_order))
p_canon <- DotPlot(seu, features = canon, group.by = primary_col) + RotatedAxis() +
  scale_y_discrete(labels = setNames(row_labels, cluster_order)) +
  labs(x = NULL, y = NULL, title = "Canonical markers by new cluster") +
  theme_publication(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf_png(p_canon, file.path(OUTPUT_DIR, paste0("dotplot_canonical_", primary_col)),
             width = 0.28 * length(canon) + 5, height = 0.22 * length(cluster_order) + 2.5)

top3 <- markers %>% filter(p_val_adj < 0.05) %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n = 3, with_ties = FALSE) %>% ungroup()
top_genes <- unique(top3$gene[order(match(top3$cluster, cluster_order))])
p_top <- DotPlot(seu, features = top_genes, group.by = primary_col) + RotatedAxis() +
  scale_y_discrete(labels = setNames(row_labels, cluster_order)) +
  labs(x = NULL, y = NULL, title = "Top 3 markers per new cluster") +
  theme_publication(base_size = 9) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
save_pdf_png(p_top, file.path(OUTPUT_DIR, paste0("dotplot_top3_", primary_col)),
             width = 0.2 * length(top_genes) + 5, height = 0.22 * length(cluster_order) + 2.5)

# 6e. Relabelled fraction per sample x original type
p_relab <- ggplot(relabel_tab, aes(x = sample, y = factor(orig, levels = rev(TIER1_ORDER)), fill = frac)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", frac)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#D73027", limits = c(0, 1)) +
  labs(title = "Fraction of cells whose label changed", x = NULL, y = "Original label",
       fill = "Fraction") +
  theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf_png(p_relab, file.path(OUTPUT_DIR, "heatmap_relabelled_fraction"), width = 7, height = 6)

# =============================================================================
# 7. Export (new files only)
# =============================================================================
label_cols <- sub("clust", "label", clust_cols)
assign_df <- data.frame(cell = colnames(seu),
                        seu@meta.data[, c("sample", "condition", "total_counts", "n_genes_by_counts",
                                          ORIG_CLUSTER_COL, ORIG_LABEL_COL, clust_cols, label_cols,
                                          "redo_label", "label_changed")],
                        spatial_xy, Embeddings(seu, "umap_redo"), check.names = FALSE)
assign_df[[primary_col]] <- as.character(assign_df[[primary_col]])
write.csv(assign_df, file.path(OUTPUT_DIR, "redo_cluster_assignments.csv"), row.names = FALSE)

# Lightweight object: counts + normalized data, PCA/Harmony/UMAP, no scale.data
seu[["spatial"]] <- CreateDimReducObject(embeddings = spatial_xy, key = "spatial_", assay = EXPR_ASSAY)
seu_out <- DietSeurat(seu, layers = c("counts", "data"),
                      dimreducs = c("pca", "harmony", "umap_redo", "spatial"), graphs = NULL)
save_rds_new(seu_out, file.path(OUTPUT_DIR, "DKO_reclustered_harmony.rds"))

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))

cat(sprintf("\nPrimary: %s -> %d clusters; accuracy vs original = %.3f; ARI = %.3f\n",
            primary_col, nrow(mm_primary),
            summary_res$accuracy_vs_Tier1[summary_res$resolution == REDO_RESOLUTION],
            summary_res$ARI_vs_Tier1[summary_res$resolution == REDO_RESOLUTION]))
cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")
