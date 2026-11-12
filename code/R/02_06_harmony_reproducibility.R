# 02_06_harmony_reproducibility.R
# Which step of the 02_01 -> 02_03 pipeline makes the clusters and merged
# domains non-reproducible? 02_05 found that re-running Harmony + clustering with
# identical settings recovers the stored clusters with ARI 0.35 only. Starting
# from the stored PCA (Banksy is not re-run), each step is isolated in turn:
#   A. Louvain only          re-cluster the STORED SNN graph
#   B. Neighbours + Louvain  rebuild the SNN graph from the STORED Harmony embedding
#   C. Harmony determinism   two Harmony runs with the same seed (identical?), and
#                            how close they are to the stored embedding
#   D. Harmony seed effect   Harmony with seeds 42, 1, 2 -> cluster -> merge
# Every clustering is merged into MERGE_K domains exactly as in 02_03 (Ward.D2 on
# cluster centroids in the Harmony space it came from), and all partitions are
# compared pairwise (ARI; label-invariant).
#
# Read-only with respect to existing results: the Seurat object is never saved
# and all outputs go to a folder used only by this script.

library(Seurat)
library(harmony)
library(Matrix)
library(FNN)
library(dplyr)
library(ggplot2)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "harmony_reproducibility")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
HARMONY_SEEDS <- c(42, 1, 2)   # 42 = the seed used in 02_01
N_KNN_QUERY   <- 3000          # cells used to compare neighbourhoods between embeddings
KNN_COMPARE   <- 30

# --- Helpers -------------------------------------------------------------------
ari <- function(a, b) mclust::adjustedRandIndex(a, b)

# Merge clusters into MERGE_K domains as in 02_03 (cluster order = factor levels)
merge_domains <- function(emb, clusters) {
  clusters <- droplevels(as.factor(clusters))
  ids <- levels(clusters)
  cen <- rowsum(emb, clusters)[ids, , drop = FALSE] / as.numeric(table(clusters)[ids])
  hc  <- hclust(dist(cen, method = MERGE_DIST_METHOD), method = MERGE_LINK_METHOD)
  unname(cutree(hc, k = MERGE_K)[as.character(clusters)])
}

# Mean Jaccard overlap of the k-nearest-neighbour sets of the same query cells
# in two embeddings (rotation- and sign-invariant similarity of the embeddings).
knn_overlap <- function(emb_a, emb_b, query, k = KNN_COMPARE) {
  na <- get.knnx(emb_a, emb_a[query, , drop = FALSE], k = k + 1)$nn.index[, -1]
  nb <- get.knnx(emb_b, emb_b[query, , drop = FALSE], k = k + 1)$nn.index[, -1]
  mean(vapply(seq_along(query), function(i) {
    length(intersect(na[i, ], nb[i, ])) / length(union(na[i, ], nb[i, ]))
  }, numeric(1)))
}

# Jaccard overlap of the edge sets of two sparse graphs
edge_jaccard <- function(g1, g2) {
  n  <- nrow(g1)
  e1 <- summary(as(g1, "CsparseMatrix")); e2 <- summary(as(g2, "CsparseMatrix"))
  k1 <- e1$i + (e1$j - 1) * n; k2 <- e2$i + (e2$j - 1) * n   # numeric edge keys
  n_shared <- sum(k1 %in% k2)
  n_shared / (length(k1) + length(k2) - n_shared)
}

cluster_embedding <- function(seu, reduction, name) {
  set.seed(42)
  seu <- FindNeighbors(seu, reduction = reduction, dims = 1:BANKSY_NPC_USE,
                       graph.name = paste0(name, c("_nn", "_snn")), verbose = FALSE)
  seu <- FindClusters(seu, graph.name = paste0(name, "_snn"), resolution = BANKSY_RESOLUTION,
                      cluster.name = paste0(name, "_clusters"), verbose = FALSE)
  seu
}

# --- 1. Load object (read-only) and record versions ----------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
stopifnot(all(c("pca", "harmony") %in% Reductions(seu)))

stored_clusters <- seu@meta.data[[DOMAIN_SOURCE_COL]]
stored_domains  <- as.character(seu@meta.data[[COL_DOMAIN]])
stored_harmony  <- Embeddings(seu, "harmony")

pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
versions <- data.frame(
  item  = c("R", "Seurat", "SeuratObject", "harmony", "RcppAnnoy", "BLAS",
            "stored_object_commands", "stored_graphs"),
  value = c(R.version.string, pkg_ver("Seurat"), pkg_ver("SeuratObject"), pkg_ver("harmony"),
            pkg_ver("RcppAnnoy"), extSoftVersion()[["BLAS"]],
            paste(names(seu@commands), collapse = "; "), paste(Graphs(seu), collapse = "; ")))
write.csv(versions, file.path(OUTPUT_DIR, "versions.csv"), row.names = FALSE)
print(versions)

partitions <- list(stored = list(clusters = stored_clusters, domains = stored_domains))
checks <- list()

# --- 2. A: Louvain only, on the stored SNN graph --------------------------------
stored_snn <- grep("snn", Graphs(seu), value = TRUE)[1]
if (!is.na(stored_snn)) {
  cat("A. Re-clustering the stored graph", stored_snn, "...\n")
  seu <- FindClusters(seu, graph.name = stored_snn, resolution = BANKSY_RESOLUTION,
                      cluster.name = "A_stored_graph_clusters", verbose = FALSE)
  cl <- seu$A_stored_graph_clusters
  partitions$A_louvain_stored_graph <- list(clusters = cl, domains = merge_domains(stored_harmony, cl))
} else {
  cat("A. No stored SNN graph in the object; step A skipped.\n")
}

# --- 3. B: neighbours + Louvain on the stored Harmony embedding -----------------
cat("B. Rebuilding the SNN graph from the stored Harmony embedding...\n")
seu <- cluster_embedding(seu, "harmony", "B")
cl  <- seu$B_clusters
partitions$B_nn_louvain_stored_harmony <- list(clusters = cl, domains = merge_domains(stored_harmony, cl))
if (!is.na(stored_snn)) {
  checks$B_snn_edge_jaccard_vs_stored <- edge_jaccard(seu[[stored_snn]], seu[["B_snn"]])
}
seu[["B_nn"]] <- NULL; seu[["B_snn"]] <- NULL

# --- 4. C + D: Harmony runs -----------------------------------------------------
set.seed(7)
query <- sample(ncol(seu), N_KNN_QUERY)
emb_runs <- list()
run_harmony <- function(seed, name) {
  set.seed(seed)
  seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca",
                    dims.use = seq_len(BANKSY_NPCS), reduction.save = name, verbose = FALSE)
  Embeddings(seu, name)
}
for (s in HARMONY_SEEDS) {
  cat(sprintf("C/D. Harmony with seed %d...\n", s))
  emb_runs[[paste0("seed", s)]] <- run_harmony(s, paste0("harmony_seed", s))
}
cat("C. Harmony with seed 42, second run...\n")
emb_runs$seed42_repeat <- run_harmony(42, "harmony_seed42_repeat")

embedding_cmp <- bind_rows(lapply(names(emb_runs), function(r) data.frame(
  run = r,
  max_abs_diff_vs_stored = max(abs(emb_runs[[r]] - stored_harmony)),
  knn_jaccard_vs_stored  = knn_overlap(stored_harmony, emb_runs[[r]], query),
  max_abs_diff_vs_seed42 = max(abs(emb_runs[[r]] - emb_runs$seed42)),
  knn_jaccard_vs_seed42  = knn_overlap(emb_runs$seed42, emb_runs[[r]], query))))
write.csv(embedding_cmp, file.path(OUTPUT_DIR, "embedding_comparison.csv"), row.names = FALSE)
print(embedding_cmp, digits = 3)

# Cluster + merge each Harmony run (the repeat only if it differs from seed 42)
for (r in names(emb_runs)) {
  if (r == "seed42_repeat" && max(abs(emb_runs$seed42_repeat - emb_runs$seed42)) == 0) next
  cat(sprintf("D. Clustering Harmony run %s...\n", r))
  seu[[paste0("H_", r)]] <- CreateDimReducObject(embeddings = emb_runs[[r]], key = paste0("H", gsub("[^A-Za-z0-9]", "", r), "_"),
                                                 assay = DefaultAssay(seu))
  seu <- cluster_embedding(seu, paste0("H_", r), r)
  cl  <- seu@meta.data[[paste0(r, "_clusters")]]
  partitions[[paste0("D_harmony_", r)]] <- list(clusters = cl, domains = merge_domains(emb_runs[[r]], cl))
  seu[[paste0(r, "_nn")]] <- NULL; seu[[paste0(r, "_snn")]] <- NULL
}
rm(seu); invisible(gc())

# --- 5. Pairwise agreement -------------------------------------------------------
pnames <- names(partitions)
pairwise <- bind_rows(lapply(seq_along(pnames), function(i) bind_rows(lapply(seq_along(pnames), function(j) {
  if (j <= i) return(NULL)
  a <- partitions[[pnames[i]]]; b <- partitions[[pnames[j]]]
  data.frame(run_a = pnames[i], run_b = pnames[j],
             n_clusters_a = length(unique(a$clusters)), n_clusters_b = length(unique(b$clusters)),
             cluster_ari = ari(a$clusters, b$clusters),
             domain_ari  = ari(a$domains, b$domains))
}))))
write.csv(pairwise, file.path(OUTPUT_DIR, "pairwise_partition_agreement.csv"), row.names = FALSE)
print(pairwise, digits = 3)

if (length(checks)) {
  checks_df <- data.frame(check = names(checks), value = unlist(checks))
  write.csv(checks_df, file.path(OUTPUT_DIR, "graph_checks.csv"), row.names = FALSE)
  print(checks_df)
}

# Heatmap of pairwise domain / cluster ARI
p_df <- bind_rows(pairwise %>% transmute(a = run_a, b = run_b, level = "Clusters (res. 0.8)", ari = cluster_ari),
                  pairwise %>% transmute(a = run_a, b = run_b, level = sprintf("Domains (k = %d)", MERGE_K),
                                         ari = domain_ari))
p_df <- bind_rows(p_df, p_df %>% rename(a = b, b = a)) %>%
  mutate(a = factor(a, levels = pnames), b = factor(b, levels = pnames))
p <- ggplot(p_df, aes(a, b, fill = ari)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", ari)), size = 3) +
  scale_fill_gradient(low = "white", high = "#08519C", limits = c(0, 1), name = "ARI") +
  facet_wrap(~ level) +
  labs(x = NULL, y = NULL, title = "Reproducibility of clusters and merged domains") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "pairwise_partition_agreement.pdf"), p, width = 13, height = 6)

cat("Outputs saved to:", OUTPUT_DIR, "\n")
