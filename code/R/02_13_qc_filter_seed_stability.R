# 02_13_qc_filter_seed_stability.R
# Are the domains more reproducible across Harmony seeds once low-quality cells
# are removed?
#
# Cells are kept when total_counts >= MIN_TOTAL_COUNTS and n_genes_by_counts >=
# MIN_GENES_DETECTED (the 04_01 DEG filters). On the kept cells the PCA is
# recomputed from the stored Banksy scale.data (each feature re-centred on the
# kept cells; Banksy features themselves are not recomputed, so neighbourhood
# averages still include removed cells), then for every seed in HARMONY_SEEDS:
# Harmony (02_01 settings), SNN on Harmony 1:BANKSY_NPC_USE, Louvain at
# RESOLUTION (random.seed 0, as 02_10) and both domain recipes of 02_10.
# Agreement across seeds is compared with the all-cell results of 02_10 at the
# same resolution.
#
# Read-only: nothing is saved except CSV tables and plots, all in a folder used
# only by this script (no RDS output).

library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "qc_filter_seed_stability")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
HARMONY_SEEDS    <- c(42, 1, 2, 3, 4)   # 42 = the seed used in 02_01
RESOLUTION       <- 0.5
MIN_CLUSTER_FRAC <- 0.005               # minsize: clusters below this fraction of cells are dissolved
RECIPES  <- c(current_k4 = "unweighted", minsize_k5 = "minsize")
RECIPE_K <- c(current_k4 = MERGE_K, minsize_k5 = 5)
ALL_CELL_SUMMARY <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08",
                              "resolution_seed_stability", "stability_summary.csv")

# --- Helpers -------------------------------------------------------------------
ari <- function(a, b) mclust::adjustedRandIndex(a, b)

# Unweighted Ward.D2 tree on cluster centroids (as in 02_03)
centroid_tree <- function(emb, clusters) {
  clusters <- as.character(clusters)
  ids <- sort(unique(clusters))
  cen <- rowsum(emb, clusters)[ids, , drop = FALSE] / as.numeric(table(clusters)[ids])
  hclust(dist(cen, method = MERGE_DIST_METHOD), method = MERGE_LINK_METHOD)
}

# Dissolve clusters with fewer than min_size cells: each of their cells takes the
# label most strongly connected to it in the SNN graph among the kept clusters
# (as in 02_09).
reassign_small <- function(clusters, graph, emb, min_size) {
  cl    <- as.character(clusters)
  sizes <- table(cl)
  keep  <- sort(names(sizes)[sizes >= min_size])
  todo  <- !(cl %in% keep)
  while (any(todo)) {
    done  <- which(!todo)
    lab   <- sparseMatrix(i = done, j = match(cl[done], keep), x = 1,
                          dims = c(length(cl), length(keep)))
    votes <- graph[todo, , drop = FALSE] %*% lab
    has   <- rowSums(votes) > 0
    if (!any(has)) break
    idx <- which(todo)[has]
    cl[idx]   <- keep[max.col(as.matrix(votes[has, , drop = FALSE]), ties.method = "first")]
    todo[idx] <- FALSE
  }
  if (any(todo)) {
    cen <- rowsum(emb[!todo, , drop = FALSE], cl[!todo])[keep, , drop = FALSE] /
           as.numeric(table(cl[!todo])[keep])
    d2  <- outer(rowSums(emb[todo, , drop = FALSE]^2), rowSums(cen^2), "+") -
           2 * emb[todo, , drop = FALSE] %*% t(cen)
    cl[todo] <- keep[max.col(-d2, ties.method = "first")]
  }
  cl
}

# Cell-level domains for one recipe; NA when there are fewer clusters than k
merge_domains <- function(method, emb, clusters, graph, k) {
  labels <- as.character(clusters)
  if (method == "minsize") {
    labels <- reassign_small(labels, graph, emb, MIN_CLUSTER_FRAC * length(labels))
  }
  if (length(unique(labels)) < k) return(rep(NA_integer_, length(labels)))
  unname(cutree(centroid_tree(emb, labels), k = k)[labels])
}

# Relabel `x` so that each label takes the `ref` label it shares most cells with
# under a one-to-one (Hungarian) matching.
match_to_ref <- function(x, ref) {
  tab <- table(x, ref)
  asg <- clue::solve_LSAP(unclass(tab), maximum = TRUE)
  unname(setNames(colnames(tab)[as.integer(asg)], rownames(tab))[as.character(x)])
}

# --- 1. Load object, filter cells, recompute PCA -------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
keep <- seu$total_counts >= MIN_TOTAL_COUNTS & seu$n_genes_by_counts >= MIN_GENES_DETECTED
cat(sprintf("QC filter (total_counts >= %d, n_genes_by_counts >= %d): %d -> %d cells (%.1f%% removed)\n",
            MIN_TOTAL_COUNTS, MIN_GENES_DETECTED, length(keep), sum(keep), 100 * mean(!keep)))

samples   <- as.character(seu$sample)[keep]
samp_fac  <- factor(samples, levels = intersect(names(SAMPLE_COLORS), unique(samples)))
is_tumour <- (seu$Tier1_celltype == "Tumor")[keep]
removed_by_sample <- data.frame(sample = names(table(seu$sample)),
                                n_cells = as.integer(table(seu$sample)),
                                n_kept = as.integer(table(seu$sample[keep])[names(table(seu$sample))]))
removed_by_sample$frac_removed <- 1 - removed_by_sample$n_kept / removed_by_sample$n_cells
write.csv(removed_by_sample, file.path(OUTPUT_DIR, "qc_filter_by_sample.csv"), row.names = FALSE)
print(removed_by_sample, digits = 3)

banksy_scaled <- LayerData(seu, assay = "BANKSY", layer = "scale.data")[, keep]
rm(seu); invisible(gc())
banksy_scaled <- banksy_scaled - rowMeans(banksy_scaled)
cat("Recomputing PCA on kept cells...\n")
pca_emb <- Embeddings(RunPCA(banksy_scaled, npcs = BANKSY_NPCS, assay = "BANKSY", verbose = FALSE))
rm(banksy_scaled); invisible(gc())

# --- 2. Harmony per seed -> SNN -> Louvain -> domains --------------------------
runs <- list()
for (s in HARMONY_SEEDS) {
  name <- paste0("seed", s)
  cat(sprintf("%s: Harmony + SNN + Louvain...\n", name))
  set.seed(s)
  emb <- RunHarmony(pca_emb, meta_data = data.frame(sample = samples), vars_use = HARMONY_GROUP_VAR,
                    verbose = FALSE)
  set.seed(42)
  g   <- FindNeighbors(emb[, 1:BANKSY_NPC_USE], verbose = FALSE)$snn
  cls <- as.character(FindClusters(g, resolution = RESOLUTION, random.seed = 0, verbose = FALSE)[[1]])
  g_dgc <- as(g, "dgCMatrix")
  runs[[name]] <- list(clusters = cls,
                       domains = lapply(setNames(names(RECIPES), names(RECIPES)), function(rc)
                         merge_domains(RECIPES[[rc]], emb[, 1:BANKSY_NPC_USE], cls, g_dgc, RECIPE_K[[rc]])))
  cat(sprintf("  %d clusters\n", length(unique(cls))))
  rm(emb, g, g_dgc); invisible(gc())
}

# --- 3. Agreement across seeds -------------------------------------------------
pairwise <- bind_rows(lapply(combn(names(runs), 2, simplify = FALSE), function(p) {
  a <- runs[[p[1]]]; b <- runs[[p[2]]]
  out <- data.frame(seed_a = p[1], seed_b = p[2],
                    n_clusters_a = length(unique(a$clusters)), n_clusters_b = length(unique(b$clusters)),
                    cluster_ari = ari(a$clusters, b$clusters))
  for (rc in names(RECIPES)) {
    out[[paste0("domain_ari_", rc)]] <- if (anyNA(a$domains[[rc]]) || anyNA(b$domains[[rc]])) NA
                                        else ari(a$domains[[rc]], b$domains[[rc]])
  }
  out
}))
write.csv(pairwise, file.path(OUTPUT_DIR, "pairwise_seed_agreement.csv"), row.names = FALSE)

unanimity <- sapply(names(RECIPES), function(rc) {
  doms <- lapply(runs, function(r) r$domains[[rc]])
  if (any(sapply(doms, anyNA))) return(NA)
  matched <- sapply(doms, match_to_ref, ref = doms[[1]])
  mean(rowSums(matched == matched[, 1]) == ncol(matched))
})

tumor_core <- bind_rows(lapply(names(runs), function(sn) bind_rows(lapply(names(RECIPES), function(rc) {
  d <- runs[[sn]]$domains[[rc]]
  if (anyNA(d)) return(NULL)
  tc <- names(which.max(tapply(is_tumour, d, mean)))
  data.frame(seed = sn, recipe = rc, sample = levels(samp_fac),
             core_frac_of_sample = as.numeric(tapply(d == tc, samp_fac, mean)))
}))))
write.csv(tumor_core, file.path(OUTPUT_DIR, "tumor_core_by_sample.csv"), row.names = FALSE)

# --- 4. Summary vs all cells (02_10) -------------------------------------------
summary_df <- data.frame(
  cells = sprintf("QC-filtered (%d)", length(samples)), resolution = RESOLUTION,
  n_clusters = paste(range(c(pairwise$n_clusters_a, pairwise$n_clusters_b)), collapse = "-"),
  cluster_ari_mean = mean(pairwise$cluster_ari),
  domain_ari_current_k4_mean = mean(pairwise$domain_ari_current_k4),
  domain_ari_current_k4_min  = min(pairwise$domain_ari_current_k4),
  domain_ari_minsize_k5_mean = mean(pairwise$domain_ari_minsize_k5),
  domain_ari_minsize_k5_min  = min(pairwise$domain_ari_minsize_k5),
  unanimous_current_k4 = unanimity[["current_k4"]],
  unanimous_minsize_k5 = unanimity[["minsize_k5"]])
if (file.exists(ALL_CELL_SUMMARY)) {
  summary_df <- bind_rows(summary_df,
    read.csv(ALL_CELL_SUMMARY) %>% filter(resolution == RESOLUTION) %>%
      transmute(cells = "All cells (02_10)", resolution, n_clusters, cluster_ari_mean,
                domain_ari_current_k4_mean, domain_ari_current_k4_min,
                domain_ari_minsize_k5_mean, domain_ari_minsize_k5_min,
                unanimous_current_k4, unanimous_minsize_k5))
}
write.csv(summary_df, file.path(OUTPUT_DIR, "stability_summary.csv"), row.names = FALSE)
cat("\nStability across Harmony seeds:\n")
print(summary_df, digits = 3)

cat("\nTumor Core share of each sample (median over seeds [min-max]):\n")
print(as.data.frame(tumor_core %>%
  group_by(recipe, sample) %>%
  summarise(v = sprintf("%.2f [%.2f-%.2f]", median(core_frac_of_sample),
                        min(core_frac_of_sample), max(core_frac_of_sample)), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = v)))

# --- 5. Plot -------------------------------------------------------------------
p_ari <- pairwise %>%
  pivot_longer(c(cluster_ari, starts_with("domain_ari_")), names_to = "level", values_to = "ari") %>%
  mutate(level = factor(c(cluster_ari = "Clusters", domain_ari_current_k4 = sprintf("Domains, unweighted k = %d", MERGE_K),
                          domain_ari_minsize_k5 = "Domains, minsize k = 5")[level],
                        levels = c("Clusters", sprintf("Domains, unweighted k = %d", MERGE_K), "Domains, minsize k = 5"))) %>%
  ggplot(aes(x = level, y = ari)) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 1.5, color = "grey40") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#D95F02", linewidth = 0.4) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "ARI between Harmony seeds (pairs)",
       title = sprintf("QC-filtered cells, resolution %.1f", RESOLUTION)) +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "qc_filter_seed_stability.pdf"), p_ari, width = 7, height = 5)

cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")
