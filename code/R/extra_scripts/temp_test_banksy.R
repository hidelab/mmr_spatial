# temp_test_banksy.R
# Reproducibility check for 02_01: re-run the full 02_01 pipeline (variable
# features -> RunBanksy lambda=0.8 -> PCA -> Harmony (seed 42) -> neighbours ->
# Louvain res 0.8 (seed 42)) and confirm the Harmony embedding is identical to
# the stored one. Every intermediate (variable features, BANKSY matrix, PCA,
# Harmony, clusters) is compared so any divergence is localised to its step.
# Harmony is also run a second time in-session to check same-seed determinism.
#
# Input: SEURAT_RDS (the true 02_01 input) if it exists; otherwise the 02_01
# input is rebuilt from DKO_banksy_merged.rds by stripping everything 02_01 and
# later steps added (BANKSY assay, pca/harmony, BANKSY graphs, cluster columns).
#
# Read-only with respect to existing results: the Seurat object is never saved
# and all outputs go to a folder used only by this script.

library(Seurat)
library(Banksy)
library(harmony)
library(SeuratWrappers)
library(dplyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "temp_test_banksy")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

N_CELL_SUBSET <- 5000   # cells used to compare the (large, dense) BANKSY matrices

# --- Helpers -------------------------------------------------------------------
ari <- function(a, b) mclust::adjustedRandIndex(a, b)

# Matrix comparison: exact identity, max abs difference, and per-column |r|
# (|r| is sign-invariant, so a PCA sign flip still shows as ~1).
compare_mat <- function(a, b, step) {
  stopifnot(identical(dim(a), dim(b)))
  a <- as.matrix(a); b <- as.matrix(b)
  col_r <- vapply(seq_len(ncol(a)), function(j) abs(cor(a[, j], b[, j])), numeric(1))
  data.frame(step = step,
             identical = identical(unname(a), unname(b)),
             all_equal_1e8 = isTRUE(all.equal(unname(a), unname(b), tolerance = 1e-8)),
             max_abs_diff = max(abs(a - b)),
             min_abs_col_cor = min(col_r),
             median_abs_col_cor = median(col_r))
}

# --- 1. Reference results from the stored object --------------------------------
merged_rds <- file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds")
ref <- readRDS(merged_rds)
cat("Loaded reference:", ncol(ref), "cells\n")
stopifnot(all(c("pca", "harmony") %in% Reductions(ref)), "BANKSY" %in% Assays(ref))

ref_cells    <- colnames(ref)
ref_varfeat  <- VariableFeatures(ref[[DefaultAssay(ref)]])
if (DefaultAssay(ref) == "BANKSY") ref_varfeat <- VariableFeatures(ref[["RNA"]])
ref_pca      <- Embeddings(ref, "pca")
ref_harmony  <- Embeddings(ref, "harmony")
ref_clusters <- ref@meta.data[[DOMAIN_SOURCE_COL]]

set.seed(7)
cell_subset <- sort(sample(ncol(ref), N_CELL_SUBSET))
ref_banksy_layers <- intersect(c("data", "scale.data"), Layers(ref[["BANKSY"]]))
ref_banksy <- lapply(setNames(ref_banksy_layers, ref_banksy_layers), function(l) {
  m <- LayerData(ref, assay = "BANKSY", layer = l)
  as.matrix(m[, cell_subset, drop = FALSE])
})
ref_banksy_features <- rownames(ref[["BANKSY"]])

# --- 2. Build the 02_01 input -----------------------------------------------------
if (file.exists(SEURAT_RDS)) {
  cat("Using the true 02_01 input:", SEURAT_RDS, "\n")
  rm(ref); invisible(gc())
  seu <- readRDS(SEURAT_RDS)
  input_source <- SEURAT_RDS
} else {
  cat("SEURAT_RDS not found; rebuilding the 02_01 input from", merged_rds, "\n")
  seu <- ref; rm(ref); invisible(gc())
  DefaultAssay(seu) <- "RNA"
  seu[["BANKSY"]] <- NULL
  for (r in intersect(c("pca", "harmony", "umap_banksy"), Reductions(seu))) seu[[r]] <- NULL
  for (g in grep("^BANKSY", Graphs(seu), value = TRUE)) seu[[g]] <- NULL
  added_cols <- grep("banksy|BANKSY|seurat_clusters|domain|spatial_x|spatial_y",
                     colnames(seu@meta.data), value = TRUE, ignore.case = TRUE)
  seu@meta.data[, added_cols] <- NULL
  seu@commands <- list()
  input_source <- paste(merged_rds, "(stripped)")
}
invisible(gc())
stopifnot(identical(colnames(seu), ref_cells))
cat("02_01 input:", ncol(seu), "cells,", nrow(seu), "genes\n")

# ===== 3. 02_01 steps, verbatim ===================================================
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

spatial_coords <- Embeddings(seu, "spatial")
seu$spatial_x <- spatial_coords[, 1]
seu$spatial_y <- spatial_coords[, 2]

expr_assay <- DefaultAssay(seu)
cat("Expression assay:", expr_assay, "\n")

if (is.null(seu[[expr_assay]]@layers[["data"]]) ||
    identical(GetAssayData(seu, assay = expr_assay, layer = "data"),
              GetAssayData(seu, assay = expr_assay, layer = "counts"))) {
  seu <- NormalizeData(seu, assay = expr_assay)
  cat("Normalized expression data.\n")
} else {
  cat("Data already normalized, skipping.\n")
}

seu <- FindVariableFeatures(seu, assay = expr_assay, nfeatures = 2000)

seu <- RunBanksy(seu,
                 lambda = c(BANKSY_LAMBDA),
                 k_geom = c(BANKSY_K_GEOM),
                 dimx = "spatial_x",
                 dimy = "spatial_y",
                 assay = expr_assay,
                 slot = "data",
                 features = "variable",
                 group = "sample",
                 split.scale = TRUE)

seu <- RunPCA(seu, assay = "BANKSY", features = rownames(seu[["BANKSY"]]), npcs = BANKSY_NPCS)

set.seed(42)
seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca")

set.seed(42)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:BANKSY_NPC_USE)
seu <- FindClusters(seu, resolution = c(BANKSY_RESOLUTION),
                    cluster.name = paste0("banksy_clust_pc", BANKSY_NPC_USE))

primary_clust <- paste0("banksy_clust_pc", BANKSY_NPC_USE, "_res.", BANKSY_RESOLUTION)
if (!primary_clust %in% colnames(seu@meta.data)) {
  primary_clust <- grep("^banksy_clust_pc", colnames(seu@meta.data), value = TRUE)[1]
}
new_clusters <- seu@meta.data[[primary_clust]]
# ===== end of 02_01 steps =========================================================

# --- 4. Same-seed Harmony repeat in this session ---------------------------------
set.seed(42)
seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca",
                  reduction.save = "harmony_repeat")

# --- 5. Compare every step against the stored results -----------------------------
new_varfeat <- VariableFeatures(seu[[expr_assay]])
feature_check <- data.frame(
  step = c("variable_features (set)", "variable_features (order)", "BANKSY features (order)"),
  identical = c(setequal(new_varfeat, ref_varfeat),
                identical(new_varfeat, ref_varfeat),
                identical(rownames(seu[["BANKSY"]]), ref_banksy_features)))
print(feature_check)
write.csv(feature_check, file.path(OUTPUT_DIR, "feature_check.csv"), row.names = FALSE)

cmp <- list()
if (identical(rownames(seu[["BANKSY"]]), ref_banksy_features)) {
  for (l in ref_banksy_layers) {
    new_m <- as.matrix(LayerData(seu, assay = "BANKSY", layer = l)[, cell_subset, drop = FALSE])
    cmp[[paste0("banksy_", l)]] <- compare_mat(new_m, ref_banksy[[l]],
                                               sprintf("BANKSY %s (%d-cell subset)", l, N_CELL_SUBSET))
  }
}
cmp$pca             <- compare_mat(Embeddings(seu, "pca"), ref_pca, "PCA vs stored")
cmp$harmony         <- compare_mat(Embeddings(seu, "harmony"), ref_harmony, "Harmony vs stored")
cmp$harmony_repeat  <- compare_mat(Embeddings(seu, "harmony_repeat"), Embeddings(seu, "harmony"),
                                   "Harmony repeat vs Harmony (same session, seed 42)")
cmp_df <- bind_rows(cmp)
print(cmp_df, digits = 4)
write.csv(cmp_df, file.path(OUTPUT_DIR, "embedding_comparison.csv"), row.names = FALSE)

clust_tab <- table(rerun = new_clusters, stored = ref_clusters)
cluster_check <- data.frame(
  n_clusters_stored = length(unique(ref_clusters)),
  n_clusters_rerun  = length(unique(new_clusters)),
  labels_identical  = identical(as.character(new_clusters), as.character(ref_clusters)),
  ari               = ari(new_clusters, ref_clusters))
print(cluster_check)
write.csv(cluster_check, file.path(OUTPUT_DIR, "cluster_check.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(clust_tab), file.path(OUTPUT_DIR, "cluster_crosstab_rerun_vs_stored.csv"))

pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
run_info <- data.frame(
  item  = c("input", "R", "Seurat", "SeuratObject", "Banksy", "harmony", "SeuratWrappers",
            "RcppAnnoy", "BLAS", "n_cells"),
  value = c(input_source, R.version.string, pkg_ver("Seurat"), pkg_ver("SeuratObject"),
            pkg_ver("Banksy"), pkg_ver("harmony"), pkg_ver("SeuratWrappers"),
            pkg_ver("RcppAnnoy"), extSoftVersion()[["BLAS"]], ncol(seu)))
write.csv(run_info, file.path(OUTPUT_DIR, "run_info.csv"), row.names = FALSE)

# --- 6. Verdict -------------------------------------------------------------------
cat("\n================ VERDICT ================\n")
cat("Harmony re-run identical to stored:      ", cmp$harmony$identical, "\n")
cat("Harmony re-run all.equal (tol 1e-8):     ", cmp$harmony$all_equal_1e8, "\n")
cat("Harmony same-seed repeat identical:      ", cmp$harmony_repeat$identical, "\n")
cat("Clusters identical / ARI vs stored:      ", cluster_check$labels_identical, "/",
    round(cluster_check$ari, 4), "\n")
cat("Outputs saved to:", OUTPUT_DIR, "\n")
