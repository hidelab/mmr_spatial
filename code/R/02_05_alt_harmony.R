# 02_05_alt_harmony.R
# Would the spatial clusters and merged domains change if Harmony were run on
# PCs 1-30 (the dimensions used for neighbour detection) instead of all 50 PCs?
#
# 02_01 computes 50 PCs on the Banksy matrix, runs Harmony on all 50, builds the
# SNN graph on Harmony 1-30 and clusters at resolution 0.8; 02_03 then merges the
# clusters on centroids in all 50 Harmony dimensions. Harmony is not separable by
# dimension (its soft clustering and corrections use every input PC), so
# Harmony 1-30 from a 50-PC run differs from a 30-PC run.
#
# Starting from the stored PCA (Banksy is not re-run), this script compares:
#   current        clusters from 02_01, merged on Harmony(50 PCs) 1-50  [= 02_03]
#   current_dims30 clusters from 02_01, merged on Harmony(50 PCs) 1-30
#   rerun_pc50     Harmony on PCs 1-50 re-run + re-clustered, merged on 1-50
#                  (control: run-to-run reproducibility of the 02_01 pipeline)
#   alt_pc30       Harmony on PCs 1-30 + re-clustered, merged on 1-30
# For each: cluster agreement with 02_01, k = MERGE_K domain agreement with
# 02_03, and the 02_04 tree/stability metrics for k = 2..10.
#
# Read-only with respect to existing results: the Seurat object is never saved
# and all outputs go to a folder used only by this script.

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(cluster)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "alt_harmony_pc30")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
ALT_NPCS   <- BANKSY_NPC_USE   # PCs given to Harmony in the alternative run (30)
K_RANGE    <- 2:10
N_BOOT     <- 1000             # stratified tumour-bootstrap replicates
PAC_BOUNDS <- c(0.1, 0.9)

# --- Helpers -------------------------------------------------------------------
# Adjusted Rand index from a (weighted) contingency table of two partitions.
ari_from_table <- function(tab) {
  pairs <- function(x) sum(x * (x - 1) / 2)
  n_pairs  <- pairs(sum(tab))
  s_a      <- pairs(rowSums(tab))
  s_b      <- pairs(colSums(tab))
  expected <- s_a * s_b / n_pairs
  max_idx  <- (s_a + s_b) / 2
  if (max_idx == expected) return(NA_real_)
  (pairs(tab) - expected) / (max_idx - expected)
}
ari <- function(a, b) ari_from_table(table(a, b))
weighted_ari <- function(a, b, w) {
  tab <- tapply(w, list(a, b), sum)
  tab[is.na(tab)] <- 0
  ari_from_table(tab)
}

# Normalized mutual information (arithmetic-mean normalization) of two partitions.
nmi <- function(a, b) {
  p  <- table(a, b) / length(a)
  pa <- rowSums(p); pb <- colSums(p)
  nz <- p > 0
  mi <- sum(p[nz] * log(p[nz] / outer(pa, pb)[nz]))
  h  <- function(q) -sum(q[q > 0] * log(q[q > 0]))
  2 * mi / (h(pa) + h(pb))
}

centroid_tree <- function(x) {
  hclust(dist(x, method = MERGE_DIST_METHOD), method = MERGE_LINK_METHOD)
}

# Merge clusters into domains and score the merge, as in 02_04: cluster
# centroids in `emb`, Ward.D2 tree, cut at each k; tree lifetime, centroid
# silhouette, stratified tumour bootstrap (cell-weighted ARI vs the full-data cut,
# PAC) and leave-one-sample-out ARI. `clusters` is a factor whose level order
# fixes the cutree() numbering.
merge_and_score <- function(emb, clusters, samples, conditions, label) {
  clusters <- droplevels(clusters)
  ids      <- levels(clusters)
  samp     <- sort(unique(samples))
  grp      <- paste(samples, clusters, sep = "__")
  grp_sum  <- rowsum(emb, grp)
  grp_n    <- table(grp)
  sum_arr  <- array(0, dim = c(length(samp), length(ids), ncol(emb)), dimnames = list(samp, ids, NULL))
  n_mat    <- matrix(0, length(samp), length(ids), dimnames = list(samp, ids))
  for (g in rownames(grp_sum)) {
    parts <- strsplit(g, "__")[[1]]
    sum_arr[parts[1], parts[2], ] <- grp_sum[g, ]
    n_mat[parts[1], parts[2]]     <- grp_n[[g]]
  }
  samp_cond <- tapply(as.character(conditions), samples, `[`, 1)[samp]
  centroids_w <- function(w) {
    n_w  <- colSums(n_mat * w); keep <- n_w > 0
    apply(sum_arr * w, c(2, 3), sum)[keep, , drop = FALSE] / n_w[keep]
  }
  cluster_n <- colSums(n_mat)

  cen  <- centroids_w(setNames(rep(1, length(samp)), samp))
  hc   <- centroid_tree(cen)
  cuts <- sapply(K_RANGE, function(k) cutree(hc, k = k))
  colnames(cuts) <- paste0("k", K_RANGE)

  heights  <- sort(hc$height)
  n_items  <- nrow(cen)
  lifetime <- sapply(K_RANGE, function(k) heights[n_items - k + 1] - heights[n_items - k])
  d_cen    <- dist(cen, method = MERGE_DIST_METHOD)
  sil      <- sapply(K_RANGE, function(k) mean(silhouette(cuts[, paste0("k", k)], d_cen)[, "sil_width"]))

  resample <- function(w) {
    cb <- centroids_w(w); hb <- centroid_tree(cb); keep <- rownames(cb)
    lapply(K_RANGE, function(k) cutree(hb, k = k)[keep])
  }
  boot_ari <- matrix(NA_real_, N_BOOT, length(K_RANGE))
  co   <- lapply(K_RANGE, function(k) matrix(0, length(ids), length(ids), dimnames = list(ids, ids)))
  both <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (b in seq_len(N_BOOT)) {
    w <- setNames(rep(0, length(samp)), samp)
    for (cond in unique(samp_cond)) {
      pool <- samp[samp_cond == cond]
      w[pool] <- as.numeric(table(factor(sample(pool, length(pool), replace = TRUE), levels = pool)))
    }
    cuts_b <- resample(w)
    keep   <- names(cuts_b[[1]])
    both[keep, keep] <- both[keep, keep] + 1
    for (i in seq_along(K_RANGE)) {
      boot_ari[b, i] <- weighted_ari(cuts[keep, i], cuts_b[[i]], cluster_n[keep])
      co[[i]][keep, keep] <- co[[i]][keep, keep] + outer(cuts_b[[i]], cuts_b[[i]], "==")
    }
  }
  pac <- sapply(co, function(m) {
    v <- (m / both)[upper.tri(m)]; v <- v[is.finite(v)]
    mean(v > PAC_BOUNDS[1] & v < PAC_BOUNDS[2])
  })
  loso_ari <- sapply(samp, function(s) {
    w <- setNames(rep(1, length(samp)), samp); w[s] <- 0
    cuts_s <- resample(w); keep <- names(cuts_s[[1]])
    sapply(seq_along(K_RANGE), function(i) weighted_ari(cuts[keep, i], cuts_s[[i]], cluster_n[keep]))
  })

  metrics <- data.frame(
    variant = label, k = K_RANGE,
    lifetime_rel = lifetime / max(heights),
    silhouette_centroids = sil,
    ari_boot_mean = colMeans(boot_ari, na.rm = TRUE),
    ari_loso_mean = rowMeans(loso_ari, na.rm = TRUE),
    ari_loso_min  = apply(loso_ari, 1, min, na.rm = TRUE),
    pac_boot = pac)
  cat(sprintf("  %s: %d clusters merged and scored\n", label, length(ids)))
  # Cell-level domain labels at k = MERGE_K
  domain_cells <- cuts[as.character(clusters), paste0("k", MERGE_K)]
  list(metrics = metrics, domain_cells = unname(domain_cells), n_clusters = length(ids))
}

# --- 1. Load object (read-only) ------------------------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
stopifnot(all(c("pca", "harmony") %in% Reductions(seu)))
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
cat(sprintf("Stored PCA: %d components; stored Harmony: %d components\n",
            ncol(Embeddings(seu, "pca")), ncol(Embeddings(seu, "harmony"))))

clusters_current <- seu@meta.data[[DOMAIN_SOURCE_COL]]
domains_current  <- as.character(seu@meta.data[[COL_DOMAIN]])
samples_all      <- as.character(seu$sample)
cond_all         <- as.character(seu$condition)

# --- 2. Re-run Harmony and clustering ------------------------------------------
# Same settings as 02_01 (seed, batch variable, SNN on 30 components, resolution);
# only the PCs given to Harmony differ between the two runs.
rerun_pipeline <- function(npcs, name) {
  set.seed(42)
  seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca",
                    dims.use = seq_len(npcs), reduction.save = paste0("harmony_", name))
  set.seed(42)
  seu <- FindNeighbors(seu, reduction = paste0("harmony_", name), dims = 1:BANKSY_NPC_USE,
                       graph.name = paste0(name, c("_nn", "_snn")))
  seu <- FindClusters(seu, graph.name = paste0(name, "_snn"), resolution = BANKSY_RESOLUTION,
                      cluster.name = paste0(name, "_clusters"))
  list(emb = Embeddings(seu, paste0("harmony_", name)),
       clusters = seu@meta.data[[paste0(name, "_clusters")]])
}
cat("Re-running Harmony on PCs 1-50 (control)...\n")
run50 <- rerun_pipeline(BANKSY_NPCS, "pc50")
cat(sprintf("Re-running Harmony on PCs 1-%d...\n", ALT_NPCS))
run30 <- rerun_pipeline(ALT_NPCS, "pc30")

harmony_current <- Embeddings(seu, "harmony")
rm(seu); invisible(gc())

# --- 3. Merge and score each variant -------------------------------------------
cat("Merging and scoring...\n")
set.seed(42)
variants <- list(
  current        = merge_and_score(harmony_current, clusters_current, samples_all, cond_all, "current"),
  current_dims30 = merge_and_score(harmony_current[, 1:BANKSY_NPC_USE], clusters_current,
                                   samples_all, cond_all, "current_dims30"),
  rerun_pc50     = merge_and_score(run50$emb, run50$clusters, samples_all, cond_all, "rerun_pc50"),
  alt_pc30       = merge_and_score(run30$emb, run30$clusters, samples_all, cond_all, "alt_pc30")
)
cluster_sets <- list(current = clusters_current, current_dims30 = clusters_current,
                     rerun_pc50 = run50$clusters, alt_pc30 = run30$clusters)

# Reference check: the "current" variant must reproduce the 02_03 domains
stopifnot(ari(variants$current$domain_cells, domains_current) == 1)

# --- 4. Agreement with the current clusters and domains ------------------------
comparison <- bind_rows(lapply(names(variants), function(v) {
  m4 <- variants[[v]]$metrics[variants[[v]]$metrics$k == MERGE_K, ]
  data.frame(
    variant                = v,
    harmony_input_pcs      = c(current = 50, current_dims30 = 50, rerun_pc50 = 50, alt_pc30 = ALT_NPCS)[[v]],
    merge_dims             = c(current = 50, current_dims30 = BANKSY_NPC_USE, rerun_pc50 = 50,
                               alt_pc30 = ALT_NPCS)[[v]],
    n_clusters             = variants[[v]]$n_clusters,
    cluster_ari_vs_current = ari(cluster_sets[[v]], clusters_current),
    cluster_nmi_vs_current = nmi(cluster_sets[[v]], clusters_current),
    domain_ari_vs_current  = ari(variants[[v]]$domain_cells, domains_current),
    domain_nmi_vs_current  = nmi(variants[[v]]$domain_cells, domains_current),
    lifetime_rel_k         = m4$lifetime_rel,
    max_lifetime_k_ge3     = { m <- variants[[v]]$metrics; m$k[m$k >= 3][which.max(m$lifetime_rel[m$k >= 3])] },
    silhouette_k           = m4$silhouette_centroids,
    ari_boot_mean_k        = m4$ari_boot_mean,
    ari_loso_min_k         = m4$ari_loso_min,
    pac_boot_k             = m4$pac_boot)
}))
write.csv(comparison, file.path(OUTPUT_DIR, "variant_comparison.csv"), row.names = FALSE)
print(comparison, digits = 3)

metrics_all <- bind_rows(lapply(variants, `[[`, "metrics"))
write.csv(metrics_all, file.path(OUTPUT_DIR, "merge_metrics_by_variant.csv"), row.names = FALSE)

# Domain cross-tabulations (cells): each variant's k = MERGE_K domains vs 02_03
current_names <- unname(DOMAIN_NAMES[domains_current])
crosstab <- bind_rows(lapply(setdiff(names(variants), "current"), function(v) {
  as.data.frame(table(current_domain = current_names, variant_domain = variants[[v]]$domain_cells),
                responseName = "n_cells") %>% mutate(variant = v, .before = 1)
}))
write.csv(crosstab, file.path(OUTPUT_DIR, "domain_crosstab_vs_current.csv"), row.names = FALSE)

# Cell-level assignments for the re-run variants
assign_df <- data.frame(sample = samples_all, condition = cond_all,
                        current_cluster = as.character(clusters_current), current_domain = current_names,
                        rerun_pc50_cluster = as.character(run50$clusters),
                        rerun_pc50_domain = variants$rerun_pc50$domain_cells,
                        alt_pc30_cluster = as.character(run30$clusters),
                        alt_pc30_domain = variants$alt_pc30$domain_cells)
saveRDS(assign_df, file.path(OUTPUT_DIR, "cell_assignments_alt_harmony.rds"))

# --- 5. Plots ------------------------------------------------------------------
metric_labels <- c(lifetime_rel = "Dendrogram lifetime (fraction of max height)",
                   silhouette_centroids = "Silhouette (centroids)",
                   ari_boot_mean = "ARI, tumour bootstrap",
                   pac_boot = "PAC, tumour bootstrap (lower = more stable)")
p_metrics <- metrics_all %>%
  pivot_longer(all_of(names(metric_labels)), names_to = "metric") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels)) %>%
  ggplot(aes(x = k, y = value, color = variant)) +
  geom_vline(xintercept = MERGE_K, linetype = "dashed", color = "grey50") +
  geom_line() + geom_point(size = 1.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "Number of merged domains (k)", y = NULL, color = "Variant",
       title = "Domain merge: Harmony on 50 vs 30 PCs") +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "merge_metrics_by_variant.pdf"), p_metrics, width = 10, height = 7)

# Domain cross-tab heatmap (row-normalized: where each current domain's cells go)
p_cross <- crosstab %>%
  group_by(variant, current_domain) %>% mutate(frac = n_cells / sum(n_cells)) %>% ungroup() %>%
  mutate(current_domain = factor(current_domain, levels = rev(DOMAIN_ORDER))) %>%
  ggplot(aes(x = variant_domain, y = current_domain, fill = frac)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", frac)), size = 3) +
  scale_fill_gradient(low = "white", high = "#08519C", limits = c(0, 1), name = "Fraction of\ncurrent domain") +
  facet_wrap(~ variant, nrow = 1) +
  labs(x = sprintf("Variant domain (k = %d)", MERGE_K), y = "Current domain (02_03)") +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "domain_crosstab_vs_current.pdf"), p_cross, width = 12, height = 4)

cat("Outputs saved to:", OUTPUT_DIR, "\n")