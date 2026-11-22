# 02_04_merge_robustness.R
# Robustness of the hierarchical merge of Banksy clusters into tissue domains
# (02_03). Asks whether cutting the Ward.D2 dendrogram of Harmony-space cluster
# centroids at k = MERGE_K is supported by the data, by comparing k = 2..K_MAX on:
#   (1) Tree structure   — lifetime of each k-cluster solution along the
#                          dendrogram height axis; cophenetic correlation.
#   (2) Internal indices — silhouette width (centroid- and cell-level),
#                          Calinski-Harabasz index, gap statistic.
#   (3) Stability        — stratified tumour bootstrap, leave-one-sample-out and
#                          Harmony-component subsampling: ARI against the
#                          full-data partition, per-domain Jaccard, consensus PAC.
#   (4) Sensitivity      — k = MERGE_K partition under alternative linkages,
#                          distances and numbers of Harmony components.
# Partition agreement is cell-weighted: every Banksy cluster counts by its number
# of cells, so ARI / Jaccard are the cell-level values.

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(cluster)
library(ComplexHeatmap)
library(circlize)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "merge_robustness")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
K_RANGE         <- 2:10    # candidate numbers of merged domains
N_BOOT          <- 1000    # tumour-bootstrap replicates
N_FEAT_SUB      <- 1000    # Harmony-component subsampling replicates
FEAT_SUB_FRAC   <- 0.8     # fraction of Harmony components kept per subsample
N_GAP_REF       <- 500     # reference data sets for the gap statistic
N_CELLS_SIL     <- 10000   # cells subsampled for the cell-level silhouette
PAC_BOUNDS      <- c(0.1, 0.9)  # consensus interval counted as ambiguous (PAC)
SENS_DIMS       <- BANKSY_NPC_USE  # reduced dimension set for the sensitivity analysis
set.seed(42)

DOMAIN_COL <- DOMAIN_SOURCE_COL

# --- Helpers -------------------------------------------------------------------
# Adjusted Rand index between partitions a and b of the same items, with item
# weights w (cell counts -> cell-level ARI, since all cells of a Banksy cluster
# share its labels).
weighted_ari <- function(a, b, w = rep(1, length(a))) {
  tab <- tapply(w, list(a, b), sum)
  tab[is.na(tab)] <- 0
  pairs <- function(x) sum(x * (x - 1) / 2)
  n_pairs  <- pairs(sum(tab))
  s_ab     <- pairs(tab)
  s_a      <- pairs(rowSums(tab))
  s_b      <- pairs(colSums(tab))
  expected <- s_a * s_b / n_pairs
  max_idx  <- (s_a + s_b) / 2
  if (max_idx == expected) return(NA_real_)
  (s_ab - expected) / (max_idx - expected)
}

# For each reference group, the best cell-weighted Jaccard overlap with any group
# of the comparison partition (clusterboot-style per-cluster stability). Groups
# with no members among the compared items get NA.
best_jaccard <- function(ref, cmp, w, groups) {
  sapply(groups, function(g) {
    in_ref <- ref == g
    if (!any(in_ref)) return(NA_real_)
    max(sapply(unique(cmp), function(h) {
      in_cmp <- cmp == h
      sum(w[in_ref & in_cmp]) / sum(w[in_ref | in_cmp])
    }))
  })
}

# Ward.D2 (or MERGE_LINK_METHOD) tree on centroid rows.
centroid_tree <- function(x, method = MERGE_LINK_METHOD, dist_method = MERGE_DIST_METHOD) {
  hclust(dist(x, method = dist_method), method = method)
}

# Calinski-Harabasz index of a partition of the rows of x.
calinski_harabasz <- function(x, cl) {
  k <- length(unique(cl)); n <- nrow(x)
  grand <- colMeans(x)
  within <- 0; between <- 0
  for (g in unique(cl)) {
    xg <- x[cl == g, , drop = FALSE]
    mg <- colMeans(xg)
    within  <- within + sum(sweep(xg, 2, mg)^2)
    between <- between + nrow(xg) * sum((mg - grand)^2)
  }
  (between / (k - 1)) / (within / (n - k))
}

# Proportion of ambiguous clustering (Senbabaoglu et al. 2014): fraction of item
# pairs whose consensus co-clustering lies strictly inside PAC_BOUNDS.
pac_score <- function(consensus) {
  v <- consensus[upper.tri(consensus)]
  v <- v[!is.na(v)]
  mean(v > PAC_BOUNDS[1] & v < PAC_BOUNDS[2])
}

# --- 1. Load Harmony embedding and cluster labels ------------------------------
BANKSY_MERGED_RDS <- file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds")
seu <- readRDS(BANKSY_MERGED_RDS)
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
# Same cluster order as 02_03 (factor-level order), so cutree() numbers the
# merged domains identically and DOMAIN_NAMES applies.
domains       <- as.character(sort(unique(seu@meta.data[[DOMAIN_COL]])))
harmony_emb   <- Embeddings(seu, "harmony")
domain_labels <- as.character(seu@meta.data[[DOMAIN_COL]])
sample_labels <- as.character(seu$sample)
cond_labels   <- as.character(seu$condition)
merged_saved  <- as.character(seu@meta.data[[COL_DOMAIN]])
rm(seu); invisible(gc())

n_dims  <- ncol(harmony_emb)
cat(sprintf("%d Banksy clusters, %d Harmony components\n", length(domains), n_dims))

# Per sample x cluster embedding sums and cell counts: every resampled centroid
# below is a weighted combination of these, so no pass over cells is needed.
samples   <- sort(unique(sample_labels))
grp       <- paste(sample_labels, domain_labels, sep = "__")
grp_sum   <- rowsum(harmony_emb, grp)
grp_n     <- table(grp)
sum_arr   <- array(0, dim = c(length(samples), length(domains), n_dims),
                   dimnames = list(samples, domains, NULL))
n_mat     <- matrix(0, length(samples), length(domains), dimnames = list(samples, domains))
for (g in rownames(grp_sum)) {
  parts <- strsplit(g, "__")[[1]]
  sum_arr[parts[1], parts[2], ] <- grp_sum[g, ]
  n_mat[parts[1], parts[2]]     <- grp_n[[g]]
}
sample_cond <- tapply(cond_labels, sample_labels, `[`, 1)[samples]

# Centroids from sample weights (a sample drawn twice counts twice). Clusters
# with no cells under the weights are dropped (returned rows = present clusters).
weighted_centroids <- function(sample_w) {
  n_w   <- colSums(n_mat * sample_w)
  keep  <- n_w > 0
  sums  <- apply(sum_arr * sample_w, c(2, 3), sum)
  sums[keep, , drop = FALSE] / n_w[keep]
}

cluster_n <- colSums(n_mat)   # cells per Banksy cluster (ARI / Jaccard weights)

# --- 2. Reference tree (reproduces 02_03) --------------------------------------
centroids <- weighted_centroids(setNames(rep(1, length(samples)), samples))
hc_ref    <- centroid_tree(centroids)
ref_cuts  <- sapply(K_RANGE, function(k) cutree(hc_ref, k = k))
colnames(ref_cuts) <- paste0("k", K_RANGE)

# Check that the reference k = MERGE_K partition equals the saved merged domains
ref_k <- ref_cuts[, paste0("k", MERGE_K)]
saved_map <- as.vector(tapply(merged_saved, domain_labels, `[`, 1)[names(ref_k)])
stopifnot(identical(as.character(unname(ref_k)), saved_map))
cat(sprintf("Reference partition at k = %d reproduces 02_03 merged domains.\n", MERGE_K))

# --- 3. Tree structure: lifetime of each k and cophenetic correlation ----------
# With merge heights h_1 <= ... <= h_{n-1}, a cut yields k clusters for any
# height in [h_{n-k}, h_{n-k+1}); the width of that interval is the k-solution's
# lifetime (larger = more clearly separated groups).
heights  <- sort(hc_ref$height)
n_items  <- nrow(centroids)
lifetime <- sapply(K_RANGE, function(k) heights[n_items - k + 1] - heights[n_items - k])
coph_cor <- cor(cophenetic(hc_ref), dist(centroids, method = MERGE_DIST_METHOD))
cat(sprintf("Cophenetic correlation: %.3f\n", coph_cor))

# --- 4. Internal validity indices ----------------------------------------------
d_centroids <- dist(centroids, method = MERGE_DIST_METHOD)
sil_centroid <- sapply(K_RANGE, function(k)
  mean(silhouette(ref_cuts[, paste0("k", k)], d_centroids)[, "sil_width"]))
ch_index <- sapply(K_RANGE, function(k) calinski_harabasz(centroids, ref_cuts[, paste0("k", k)]))

# Cell-level silhouette on a random subsample of cells in Harmony space, each cell
# labelled by the merged domain of its Banksy cluster.
sil_cells_idx <- sample(nrow(harmony_emb), min(N_CELLS_SIL, nrow(harmony_emb)))
d_cells <- dist(harmony_emb[sil_cells_idx, ])
sil_cell <- sapply(K_RANGE, function(k) {
  lab <- ref_cuts[domain_labels[sil_cells_idx], paste0("k", k)]
  mean(silhouette(lab, d_cells)[, "sil_width"])
})
rm(d_cells); invisible(gc())

# Gap statistic (Tibshirani et al. 2001): reference data uniform over the
# PCA-aligned bounding box of the centroids, clustered with the same tree.
gap <- clusGap(centroids,
               FUNcluster = function(x, k) list(cluster = cutree(centroid_tree(x), k = k)),
               K.max = max(K_RANGE), B = N_GAP_REF, spaceH0 = "scaledPCA", verbose = FALSE)
gap_tab <- as.data.frame(gap$Tab)
gap_k_first  <- maxSE(gap_tab$gap, gap_tab$SE.sim, method = "firstSEmax")
gap_k_tibs   <- maxSE(gap_tab$gap, gap_tab$SE.sim, method = "Tibs2001SEmax")

# --- 5. Resampling stability ---------------------------------------------------
# Each replicate re-estimates centroids, rebuilds the tree, cuts it at every k and
# is compared with the full-data partition at the same k over the clusters
# present in both. Consensus matrices (fraction of replicates in which two Banksy
# clusters share a domain) are accumulated for PAC.
run_replicates <- function(centroid_fun, n_rep, label) {
  ari  <- matrix(NA_real_, n_rep, length(K_RANGE), dimnames = list(NULL, paste0("k", K_RANGE)))
  jac  <- matrix(NA_real_, n_rep, MERGE_K)
  co   <- lapply(K_RANGE, function(k) matrix(0, length(domains), length(domains),
                                              dimnames = list(domains, domains)))
  both <- matrix(0, length(domains), length(domains), dimnames = list(domains, domains))
  names(co) <- paste0("k", K_RANGE)
  for (b in seq_len(n_rep)) {
    cb   <- centroid_fun(b)
    hc_b <- centroid_tree(cb)
    ids  <- rownames(cb)
    both[ids, ids] <- both[ids, ids] + 1
    for (k in K_RANGE) {
      kk  <- paste0("k", k)
      cut <- cutree(hc_b, k = k)
      ari[b, kk] <- weighted_ari(ref_cuts[ids, kk], cut[ids], cluster_n[ids])
      co[[kk]][ids, ids] <- co[[kk]][ids, ids] + outer(cut[ids], cut[ids], "==")
      if (k == MERGE_K) jac[b, ] <- best_jaccard(ref_cuts[ids, kk], cut[ids], cluster_n[ids],
                                                 groups = seq_len(MERGE_K))
    }
  }
  consensus <- lapply(co, function(m) { m <- m / both; m[both == 0] <- NA; m })
  cat(sprintf("  %s: %d replicates done\n", label, n_rep))
  list(ari = ari, jaccard = jac, consensus = consensus)
}

cat("Resampling stability...\n")

# (a) Stratified tumour bootstrap: resample WT and DKO tumours with replacement
# within genotype, so every replicate keeps the 3 WT + 4 DKO design.
boot_weights <- function(b) {
  w <- setNames(rep(0, length(samples)), samples)
  for (cond in unique(sample_cond)) {
    pool <- samples[sample_cond == cond]
    draw <- table(factor(sample(pool, length(pool), replace = TRUE), levels = pool))
    w[pool] <- as.numeric(draw)
  }
  w
}
res_boot <- run_replicates(function(b) weighted_centroids(boot_weights(b)), N_BOOT,
                           "tumour bootstrap")

# (b) Leave-one-sample-out
res_loso <- run_replicates(function(b) {
  w <- setNames(rep(1, length(samples)), samples)
  w[samples[b]] <- 0
  weighted_centroids(w)
}, length(samples), "leave-one-sample-out")
rownames(res_loso$ari) <- samples

# (c) Harmony-component subsampling: keep a random FEAT_SUB_FRAC of components
n_keep_dims <- round(FEAT_SUB_FRAC * n_dims)
res_feat <- run_replicates(function(b) {
  centroids[, sort(sample(n_dims, n_keep_dims)), drop = FALSE]
}, N_FEAT_SUB, "component subsampling")

# --- 6. Sensitivity of the k = MERGE_K partition --------------------------------
sens_specs <- list(
  list(linkage = "ward.D2",  distance = "euclidean",   dims = n_dims),
  list(linkage = "ward.D2",  distance = "euclidean",   dims = SENS_DIMS),
  list(linkage = "average",  distance = "euclidean",   dims = n_dims),
  list(linkage = "complete", distance = "euclidean",   dims = n_dims),
  list(linkage = "ward.D2",  distance = "correlation", dims = n_dims),
  list(linkage = "average",  distance = "correlation", dims = n_dims)
)
sensitivity <- bind_rows(lapply(sens_specs, function(s) {
  x  <- centroids[, seq_len(s$dims), drop = FALSE]
  d  <- if (s$distance == "correlation") as.dist(1 - cor(t(x))) else dist(x, method = s$distance)
  cl <- cutree(hclust(d, method = s$linkage), k = MERGE_K)
  data.frame(linkage = s$linkage, distance = s$distance, n_components = s$dims,
             ari_cells = weighted_ari(ref_k, cl, cluster_n),
             ari_clusters = weighted_ari(ref_k, cl),
             n_singleton_domains = sum(table(cl) == 1))
}))
cat("Sensitivity of the k =", MERGE_K, "partition:\n"); print(sensitivity)

# --- 7. Assemble tables --------------------------------------------------------
summ <- function(m) list(mean = colMeans(m, na.rm = TRUE),
                         lo = apply(m, 2, quantile, 0.025, na.rm = TRUE))
boot_s <- summ(res_boot$ari); feat_s <- summ(res_feat$ari)

k_metrics <- data.frame(
  k                      = K_RANGE,
  lifetime               = lifetime,
  lifetime_rel           = lifetime / max(heights),
  silhouette_centroids   = sil_centroid,
  silhouette_cells       = sil_cell,
  calinski_harabasz      = ch_index,
  gap                    = gap_tab$gap[K_RANGE],
  gap_se                 = gap_tab$SE.sim[K_RANGE],
  ari_boot_mean          = boot_s$mean,
  ari_boot_q025          = boot_s$lo,
  ari_loso_mean          = colMeans(res_loso$ari, na.rm = TRUE),
  ari_loso_min           = apply(res_loso$ari, 2, min, na.rm = TRUE),
  ari_featsub_mean       = feat_s$mean,
  ari_featsub_q025       = feat_s$lo,
  pac_boot               = sapply(res_boot$consensus, pac_score),
  pac_featsub            = sapply(res_feat$consensus, pac_score),
  row.names = NULL
)
write.csv(k_metrics, file.path(OUTPUT_DIR, "k_selection_metrics.csv"), row.names = FALSE)
write.csv(gap_tab %>% mutate(k = seq_len(n())) %>% relocate(k),
          file.path(OUTPUT_DIR, "gap_statistic.csv"), row.names = FALSE)

loso_long <- as.data.frame(res_loso$ari) %>%
  mutate(left_out_sample = rownames(res_loso$ari),
         condition = unname(sample_cond[left_out_sample])) %>%
  pivot_longer(starts_with("k"), names_to = "k", values_to = "ari") %>%
  mutate(k = as.integer(sub("k", "", k)))
write.csv(loso_long, file.path(OUTPUT_DIR, "loso_ari.csv"), row.names = FALSE)

domain_ids <- sort(unique(ref_k))
domain_stability <- data.frame(
  domain           = domain_ids,
  domain_name      = unname(DOMAIN_NAMES[as.character(domain_ids)]),
  n_banksy_clusters = as.integer(table(ref_k)[as.character(domain_ids)]),
  n_cells          = as.numeric(tapply(cluster_n[names(ref_k)], ref_k, sum)[as.character(domain_ids)]),
  jaccard_boot_mean    = colMeans(res_boot$jaccard, na.rm = TRUE),
  jaccard_boot_q025    = apply(res_boot$jaccard, 2, quantile, 0.025, na.rm = TRUE),
  frac_boot_ge_0.75    = colMeans(res_boot$jaccard >= 0.75, na.rm = TRUE),
  jaccard_loso_min     = apply(res_loso$jaccard, 2, min, na.rm = TRUE),
  jaccard_featsub_mean = colMeans(res_feat$jaccard, na.rm = TRUE),
  jaccard_featsub_q025 = apply(res_feat$jaccard, 2, quantile, 0.025, na.rm = TRUE),
  row.names = NULL
)
write.csv(domain_stability, file.path(OUTPUT_DIR, paste0("domain_stability_k", MERGE_K, ".csv")),
          row.names = FALSE)
write.csv(sensitivity, file.path(OUTPUT_DIR, paste0("sensitivity_k", MERGE_K, ".csv")),
          row.names = FALSE)

# --- 8. Plots ------------------------------------------------------------------
metric_labels <- c(
  lifetime_rel         = "Dendrogram lifetime\n(fraction of max height)",
  silhouette_centroids = "Silhouette\n(centroids)",
  silhouette_cells     = "Silhouette\n(cells)",
  calinski_harabasz    = "Calinski-Harabasz",
  gap                  = "Gap statistic",
  ari_boot_mean        = "ARI, tumour bootstrap",
  ari_loso_min         = "ARI, leave-one-out (min)",
  ari_featsub_mean     = "ARI, component subsampling",
  pac_boot             = "PAC, tumour bootstrap\n(lower = more stable)"
)
plot_df <- k_metrics %>%
  select(k, all_of(names(metric_labels))) %>%
  pivot_longer(-k, names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels))
gap_err <- k_metrics %>%
  transmute(k, metric = factor(metric_labels["gap"], levels = metric_labels),
            lo = gap - gap_se, hi = gap + gap_se)

p_metrics <- ggplot(plot_df, aes(x = k, y = value)) +
  geom_vline(xintercept = MERGE_K, linetype = "dashed", color = "#B71C1C") +
  geom_errorbar(data = gap_err, aes(x = k, ymin = lo, ymax = hi), inherit.aes = FALSE,
                width = 0.2, color = "grey40") +
  geom_line(color = "grey30") +
  geom_point(size = 1.8) +
  facet_wrap(~ metric, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "Number of merged domains (k)", y = NULL,
       title = "Selection of the number of merged domains",
       subtitle = sprintf("Dashed line: k = %d used in 02_03", MERGE_K)) +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "k_selection_metrics.pdf"), p_metrics, width = 11, height = 9)

# Consensus heatmap at k = MERGE_K (tumour bootstrap), dendrogram order
cons <- res_boot$consensus[[paste0("k", MERGE_K)]]
ord  <- hc_ref$labels[hc_ref$order]
DOMAIN_COLORS_NUM <- setNames(DOMAIN_COLORS[DOMAIN_NAMES[as.character(domain_ids)]],
                              DOMAIN_NAMES[as.character(domain_ids)])
dom_anno <- unname(DOMAIN_NAMES[as.character(ref_k[ord])])
ha_top  <- HeatmapAnnotation(Domain = dom_anno, col = list(Domain = DOMAIN_COLORS_NUM),
                             show_annotation_name = FALSE)
ha_left <- rowAnnotation(Domain = dom_anno, col = list(Domain = DOMAIN_COLORS_NUM),
                         show_annotation_name = FALSE, show_legend = FALSE)
pdf(file.path(OUTPUT_DIR, paste0("consensus_heatmap_k", MERGE_K, ".pdf")), width = 9, height = 8)
draw(Heatmap(cons[ord, ord], name = "Co-assignment\nfrequency",
             col = colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#6BAED6", "#08519C")),
             cluster_rows = FALSE, cluster_columns = FALSE,
             top_annotation = ha_top, left_annotation = ha_left,
             row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5),
             column_title = sprintf("Banksy cluster co-assignment at k = %d (%d tumour-bootstrap replicates)",
                                    MERGE_K, N_BOOT)))
dev.off()

# Dendrogram with the lifetime of the k = MERGE_K solution shaded
pdf(file.path(OUTPUT_DIR, "dendrogram_k_lifetime.pdf"), width = 10, height = 6)
plot(hc_ref, labels = FALSE, hang = -1, xlab = "Banksy cluster", sub = "",
     ylab = "Ward distance", main = "Domain dendrogram: height interval yielding each k")
k_band <- heights[n_items - MERGE_K + c(0, 1)]
rect(0, k_band[1], n_items + 1, k_band[2], col = adjustcolor("#B71C1C", 0.15), border = NA)
text(n_items * 0.95, mean(k_band), paste("k =", MERGE_K), col = "#B71C1C", adj = 1)
dev.off()

# --- 9. Summary ----------------------------------------------------------------
best <- function(x, fun = which.max) K_RANGE[fun(x)]
summary_lines <- c(
  strrep("=", 78),
  "02_04 robustness of the domain merge (Ward.D2 tree of Harmony centroids)",
  strrep("=", 78),
  sprintf("Banksy clusters: %d   Harmony components: %d   k evaluated: %d-%d",
          length(domains), n_dims, min(K_RANGE), max(K_RANGE)),
  sprintf("Cophenetic correlation of the tree: %.3f", coph_cor),
  "",
  "k favoured by each criterion:",
  sprintf("  Largest dendrogram lifetime:        k = %d", best(lifetime)),
  sprintf("  Max silhouette (centroids):         k = %d", best(sil_centroid)),
  sprintf("  Max silhouette (cells):             k = %d", best(sil_cell)),
  sprintf("  Max Calinski-Harabasz:              k = %d", best(ch_index)),
  sprintf("  Gap statistic (firstSEmax / Tibs2001SEmax): k = %d / %d", gap_k_first, gap_k_tibs),
  sprintf("  Min PAC (tumour bootstrap):         k = %d", best(k_metrics$pac_boot, which.min)),
  "",
  sprintf("At k = %d:", MERGE_K),
  sprintf("  ARI tumour bootstrap  mean %.3f (2.5%% quantile %.3f)",
          k_metrics$ari_boot_mean[k_metrics$k == MERGE_K], k_metrics$ari_boot_q025[k_metrics$k == MERGE_K]),
  sprintf("  ARI leave-one-out     mean %.3f (min %.3f)",
          k_metrics$ari_loso_mean[k_metrics$k == MERGE_K], k_metrics$ari_loso_min[k_metrics$k == MERGE_K]),
  sprintf("  ARI comp. subsampling mean %.3f (2.5%% quantile %.3f)",
          k_metrics$ari_featsub_mean[k_metrics$k == MERGE_K], k_metrics$ari_featsub_q025[k_metrics$k == MERGE_K]),
  sprintf("  PAC tumour bootstrap  %.3f", k_metrics$pac_boot[k_metrics$k == MERGE_K]),
  "  Per-domain Jaccard (tumour bootstrap mean / leave-one-out min):",
  sprintf("    %-18s %.3f / %.3f", domain_stability$domain_name,
          domain_stability$jaccard_boot_mean, domain_stability$jaccard_loso_min),
  "  (Hennig 2007: mean Jaccard >= 0.85 highly stable, >= 0.75 stable, < 0.6 dissolved)",
  "",
  strrep("=", 78)
)
writeLines(summary_lines, file.path(OUTPUT_DIR, "SUMMARY.txt"))
cat(summary_lines, sep = "\n")
cat("Outputs saved to:", OUTPUT_DIR, "\n")