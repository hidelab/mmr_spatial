# 02_07_no_harmony.R
# Can Harmony be dropped from the spatial-domain pipeline? 02_06 showed that
# Harmony is the only non-deterministic step and that different Harmony seeds
# change the merged domains substantially (domain ARI 0.25-0.44). Samples are
# also nested within genotype, so sample-level correction can absorb genotype
# effects. Here clusters and domains are built directly from the stored Banksy
# PCA with otherwise identical settings (SNN on PCs 1-30, Louvain resolution 0.8,
# Ward.D2 merge of cluster centroids at k = MERGE_K) and compared with the current
# Harmony-based domains on:
#   (1) Mixing          per-cell iLISI over tumours, sample composition of every
#                       domain, tumour dominance of fine clusters
#   (2) Agreement       ARI / best one-to-one relabelling vs the current domains;
#                       cell-type composition and domain-name matching
#   (3) Merge robustness  02_04 metrics (tree lifetime, silhouette, tumour
#                       bootstrap, leave-one-sample-out, per-domain Jaccard, PAC)
#   (4) Determinism     Louvain random seeds 0/1/2 (the only randomness left)
#
# Read-only with respect to existing results: the Seurat object is never saved
# and all outputs go to a folder used only by this script.

library(Seurat)
library(FNN)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(cluster)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "no_harmony")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
K_RANGE       <- 2:10
N_BOOT        <- 1000
PAC_BOUNDS    <- c(0.1, 0.9)
LOUVAIN_SEEDS <- c(0, 1, 2)    # 0 = Seurat default, used in 02_01
N_LISI_QUERY  <- 10000         # cells at which iLISI is evaluated
LISI_K        <- 90            # neighbours for iLISI (~3 x perplexity 30)
LISI_DIMS     <- BANKSY_NPC_USE
set.seed(42)

# --- Helpers -------------------------------------------------------------------
ari <- function(a, b) mclust::adjustedRandIndex(a, b)

ari_from_table <- function(tab) {
  pairs <- function(x) sum(x * (x - 1) / 2)
  n_pairs <- pairs(sum(tab)); s_a <- pairs(rowSums(tab)); s_b <- pairs(colSums(tab))
  expected <- s_a * s_b / n_pairs; max_idx <- (s_a + s_b) / 2
  if (max_idx == expected) return(NA_real_)
  (pairs(tab) - expected) / (max_idx - expected)
}
weighted_ari <- function(a, b, w) {
  tab <- tapply(w, list(a, b), sum); tab[is.na(tab)] <- 0
  ari_from_table(tab)
}

# Fraction of cells agreeing under the best one-to-one relabelling of two
# partitions with the same number of groups (all permutations; k is small).
best_relabel_agreement <- function(a, b) {
  tab <- table(a, b)
  if (nrow(tab) != ncol(tab)) return(NA_real_)
  perms <- function(v) if (length(v) <= 1) list(v) else
    do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
  max(sapply(perms(seq_len(ncol(tab))), function(p) sum(tab[cbind(seq_len(nrow(tab)), p)]))) / sum(tab)
}

# iLISI (Korsunsky et al. 2019), approximated as the inverse Simpson index of
# tumour labels among each query cell's LISI_K nearest neighbours: 1 = all
# neighbours from one tumour, up to 7 = evenly mixed across all tumours.
ilisi <- function(emb, samples, query) {
  nn <- get.knnx(emb, emb[query, , drop = FALSE], k = LISI_K)$nn.index
  apply(nn, 1, function(idx) { p <- table(samples[idx]) / LISI_K; 1 / sum(p^2) })
}

centroid_tree <- function(x) hclust(dist(x, method = MERGE_DIST_METHOD), method = MERGE_LINK_METHOD)

# Merge clusters into domains and score the merge as in 02_04 (cell-weighted
# stability against tumour resampling). `clusters` level order fixes cutree()
# numbering.
merge_and_score <- function(emb, clusters, samples, conditions, label) {
  clusters <- droplevels(as.factor(clusters))
  ids  <- levels(clusters)
  samp <- sort(unique(samples))
  grp  <- paste(samples, clusters, sep = "__")
  grp_sum <- rowsum(emb, grp); grp_n <- table(grp)
  sum_arr <- array(0, dim = c(length(samp), length(ids), ncol(emb)), dimnames = list(samp, ids, NULL))
  n_mat   <- matrix(0, length(samp), length(ids), dimnames = list(samp, ids))
  for (g in rownames(grp_sum)) {
    parts <- strsplit(g, "__")[[1]]
    sum_arr[parts[1], parts[2], ] <- grp_sum[g, ]
    n_mat[parts[1], parts[2]]     <- grp_n[[g]]
  }
  samp_cond <- tapply(as.character(conditions), samples, `[`, 1)[samp]
  centroids_w <- function(w) {
    n_w <- colSums(n_mat * w); keep <- n_w > 0
    apply(sum_arr * w, c(2, 3), sum)[keep, , drop = FALSE] / n_w[keep]
  }
  cluster_n <- colSums(n_mat)

  cen  <- centroids_w(setNames(rep(1, length(samp)), samp))
  hc   <- centroid_tree(cen)
  cuts <- sapply(K_RANGE, function(k) cutree(hc, k = k)); colnames(cuts) <- paste0("k", K_RANGE)
  heights  <- sort(hc$height); n_items <- nrow(cen)
  lifetime <- sapply(K_RANGE, function(k) heights[n_items - k + 1] - heights[n_items - k])
  d_cen    <- dist(cen, method = MERGE_DIST_METHOD)
  sil      <- sapply(K_RANGE, function(k) mean(silhouette(cuts[, paste0("k", k)], d_cen)[, "sil_width"]))

  best_jac <- function(ref, cmp, w) sapply(seq_len(MERGE_K), function(g) {
    in_ref <- ref == g
    if (!any(in_ref)) return(NA_real_)
    max(sapply(unique(cmp), function(h) sum(w[in_ref & cmp == h]) / sum(w[in_ref | cmp == h])))
  })
  resample <- function(w) {
    cb <- centroids_w(w); hb <- centroid_tree(cb); keep <- rownames(cb)
    lapply(K_RANGE, function(k) cutree(hb, k = k)[keep])
  }
  ki <- which(K_RANGE == MERGE_K)
  boot_ari <- matrix(NA_real_, N_BOOT, length(K_RANGE))
  boot_jac <- matrix(NA_real_, N_BOOT, MERGE_K)
  co   <- lapply(K_RANGE, function(k) matrix(0, length(ids), length(ids), dimnames = list(ids, ids)))
  both <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  for (b in seq_len(N_BOOT)) {
    w <- setNames(rep(0, length(samp)), samp)
    for (cond in unique(samp_cond)) {
      pool <- samp[samp_cond == cond]
      w[pool] <- as.numeric(table(factor(sample(pool, length(pool), replace = TRUE), levels = pool)))
    }
    cb <- resample(w); keep <- names(cb[[1]])
    both[keep, keep] <- both[keep, keep] + 1
    for (i in seq_along(K_RANGE)) {
      boot_ari[b, i] <- weighted_ari(cuts[keep, i], cb[[i]], cluster_n[keep])
      co[[i]][keep, keep] <- co[[i]][keep, keep] + outer(cb[[i]], cb[[i]], "==")
    }
    boot_jac[b, ] <- best_jac(cuts[keep, ki], cb[[ki]], cluster_n[keep])
  }
  pac <- sapply(co, function(m) {
    v <- (m / both)[upper.tri(m)]; v <- v[is.finite(v)]
    mean(v > PAC_BOUNDS[1] & v < PAC_BOUNDS[2])
  })
  loso <- lapply(samp, function(s) {
    w <- setNames(rep(1, length(samp)), samp); w[s] <- 0
    cb <- resample(w); keep <- names(cb[[1]])
    list(ari = sapply(seq_along(K_RANGE), function(i) weighted_ari(cuts[keep, i], cb[[i]], cluster_n[keep])),
         jac = best_jac(cuts[keep, ki], cb[[ki]], cluster_n[keep]))
  })
  loso_ari <- sapply(loso, `[[`, "ari"); loso_jac <- t(sapply(loso, `[[`, "jac"))

  metrics <- data.frame(variant = label, k = K_RANGE,
                        lifetime_rel = lifetime / max(heights),
                        silhouette_centroids = sil,
                        ari_boot_mean = colMeans(boot_ari, na.rm = TRUE),
                        ari_loso_mean = rowMeans(loso_ari, na.rm = TRUE),
                        ari_loso_min  = apply(loso_ari, 1, min, na.rm = TRUE),
                        pac_boot = pac)
  jaccard <- data.frame(variant = label, domain = seq_len(MERGE_K),
                        jaccard_boot_mean = colMeans(boot_jac, na.rm = TRUE),
                        jaccard_loso_min  = apply(loso_jac, 2, min, na.rm = TRUE))
  cat(sprintf("  %s: %d clusters merged and scored\n", label, length(ids)))
  list(metrics = metrics, jaccard = jaccard, n_clusters = length(ids),
       domain_cells = unname(cuts[as.character(clusters), paste0("k", MERGE_K)]))
}

# --- 1. Load object (read-only) ------------------------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
stopifnot(all(c("pca", "harmony", "spatial") %in% Reductions(seu)))
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)

samples_all   <- as.character(seu$sample)
cond_all      <- as.character(seu$condition)
celltype_all  <- as.character(seu$Tier1_celltype)
coords        <- Embeddings(seu, "spatial")
pca_emb       <- Embeddings(seu, "pca")
harmony_emb   <- Embeddings(seu, "harmony")
clusters_harm <- seu@meta.data[[DOMAIN_SOURCE_COL]]
domains_harm  <- as.character(seu@meta.data[[COL_DOMAIN]])

# --- 2. Cluster without Harmony ------------------------------------------------
# Same neighbour / clustering settings as 02_01, on the uncorrected PCA.
cat("Clustering on PCA (no Harmony)...\n")
seu <- FindNeighbors(seu, reduction = "pca", dims = 1:BANKSY_NPC_USE,
                     graph.name = c("noharm_nn", "noharm_snn"), verbose = FALSE)
louvain_runs <- lapply(LOUVAIN_SEEDS, function(s) {
  seu <- FindClusters(seu, graph.name = "noharm_snn", resolution = BANKSY_RESOLUTION,
                      random.seed = s, cluster.name = "tmp_clusters", verbose = FALSE)
  seu$tmp_clusters
})
names(louvain_runs) <- paste0("seed", LOUVAIN_SEEDS)
clusters_pca <- louvain_runs$seed0
rm(seu); invisible(gc())

# --- 3. Merge and score --------------------------------------------------------
cat("Merging and scoring...\n")
variants <- list(
  harmony_current = merge_and_score(harmony_emb, clusters_harm, samples_all, cond_all, "harmony_current"),
  no_harmony      = merge_and_score(pca_emb, clusters_pca, samples_all, cond_all, "no_harmony"),
  no_harmony_dims30 = merge_and_score(pca_emb[, 1:BANKSY_NPC_USE], clusters_pca, samples_all, cond_all,
                                      "no_harmony_dims30")
)
stopifnot(ari(variants$harmony_current$domain_cells, domains_harm) == 1)
cluster_sets <- list(harmony_current = clusters_harm, no_harmony = clusters_pca,
                     no_harmony_dims30 = clusters_pca)
emb_sets <- list(harmony_current = harmony_emb, no_harmony = pca_emb, no_harmony_dims30 = pca_emb)

# Louvain-seed determinism of the no-Harmony pipeline (domains merged on all PCs)
louvain_check <- bind_rows(lapply(names(louvain_runs)[-1], function(r) {
  dom_r <- merge_and_score(pca_emb, louvain_runs[[r]], samples_all, cond_all, paste0("louvain_", r))
  data.frame(louvain_seed = r, n_clusters = length(unique(louvain_runs[[r]])),
             cluster_ari_vs_seed0 = ari(louvain_runs[[r]], clusters_pca),
             domain_ari_vs_seed0  = ari(dom_r$domain_cells, variants$no_harmony$domain_cells))
}))
write.csv(louvain_check, file.path(OUTPUT_DIR, "louvain_seed_check.csv"), row.names = FALSE)
print(louvain_check, digits = 3)

# --- 4. Mixing across tumours --------------------------------------------------
lisi_query <- sample(length(samples_all), N_LISI_QUERY)
lisi <- list(harmony = ilisi(harmony_emb[, 1:LISI_DIMS], samples_all, lisi_query),
             pca     = ilisi(pca_emb[, 1:LISI_DIMS], samples_all, lisi_query))

domain_named <- function(v) {
  d <- variants[[v]]$domain_cells
  if (v == "harmony_current") unname(DOMAIN_NAMES[as.character(d)]) else paste0("D", d)
}

# Cell-type composition of each domain, and the best-matching current domain by
# composition (cosine similarity) and by cell overlap
ct_comp <- function(dom) prop.table(table(dom, celltype_all), 1)
comp_current <- ct_comp(domain_named("harmony_current"))
name_match <- bind_rows(lapply(setdiff(names(variants), "harmony_current"), function(v) {
  comp <- ct_comp(domain_named(v))
  cosine <- (comp %*% t(comp_current)) /
    outer(sqrt(rowSums(comp^2)), sqrt(rowSums(comp_current^2)))
  overlap <- prop.table(table(domain_named(v), domain_named("harmony_current")), 1)
  data.frame(variant = v, domain = rownames(comp),
             n_cells = as.integer(table(domain_named(v))[rownames(comp)]),
             best_composition_match = colnames(cosine)[apply(cosine, 1, which.max)],
             composition_cosine = apply(cosine, 1, max),
             best_overlap_match = colnames(overlap)[apply(overlap[rownames(comp), , drop = FALSE], 1, which.max)],
             overlap_fraction = apply(overlap[rownames(comp), , drop = FALSE], 1, max))
}))
write.csv(name_match, file.path(OUTPUT_DIR, "domain_name_matching.csv"), row.names = FALSE)
print(name_match, digits = 3)

celltype_comp <- bind_rows(lapply(names(variants), function(v) {
  as.data.frame(ct_comp(domain_named(v)), responseName = "fraction") %>%
    rename(domain = dom, celltype = celltype_all) %>% mutate(variant = v, .before = 1)
}))
write.csv(celltype_comp, file.path(OUTPUT_DIR, "domain_celltype_composition.csv"), row.names = FALSE)

sample_comp <- bind_rows(lapply(names(variants), function(v) {
  as.data.frame(table(domain = domain_named(v), sample = samples_all), responseName = "n_cells") %>%
    group_by(domain) %>% mutate(frac_of_domain = n_cells / sum(n_cells)) %>%
    group_by(sample) %>% mutate(frac_of_sample = n_cells / sum(n_cells)) %>%
    ungroup() %>% mutate(variant = v, .before = 1)
}))
write.csv(sample_comp, file.path(OUTPUT_DIR, "domain_sample_composition.csv"), row.names = FALSE)

lisi_by_domain <- bind_rows(lapply(names(variants), function(v) {
  emb_name <- if (v == "harmony_current") "harmony" else "pca"
  data.frame(variant = v, domain = domain_named(v)[lisi_query], ilisi = lisi[[emb_name]])
})) %>% group_by(variant, domain) %>%
  summarise(n_query_cells = n(), ilisi_median = median(ilisi), ilisi_mean = mean(ilisi), .groups = "drop")
write.csv(lisi_by_domain, file.path(OUTPUT_DIR, "ilisi_by_domain.csv"), row.names = FALSE)

# --- 5. Summary table ----------------------------------------------------------
summary_df <- bind_rows(lapply(names(variants), function(v) {
  m4  <- variants[[v]]$metrics[variants[[v]]$metrics$k == MERGE_K, ]
  m3  <- variants[[v]]$metrics[variants[[v]]$metrics$k >= 3, ]
  cl  <- cluster_sets[[v]]
  dom_share <- apply(prop.table(table(cl, samples_all), 1), 1, max)
  sc  <- sample_comp[sample_comp$variant == v, ]
  data.frame(
    variant = v,
    n_clusters = variants[[v]]$n_clusters,
    cluster_ari_vs_current = ari(cl, clusters_harm),
    domain_ari_vs_current  = ari(variants[[v]]$domain_cells, domains_harm),
    domain_best_relabel_agreement = best_relabel_agreement(variants[[v]]$domain_cells, domains_harm),
    ilisi_median = median(lisi[[if (v == "harmony_current") "harmony" else "pca"]]),
    clusters_single_tumour_gt80 = sum(dom_share > 0.8),
    clusters_single_tumour_gt50 = sum(dom_share > 0.5),
    min_tumour_share_of_any_domain = min(sc$frac_of_sample),
    max_single_tumour_share_of_domain = max(sc$frac_of_domain),
    lifetime_rel_k = m4$lifetime_rel,
    max_lifetime_k_ge3 = m3$k[which.max(m3$lifetime_rel)],
    silhouette_k = m4$silhouette_centroids,
    ari_boot_mean_k = m4$ari_boot_mean,
    ari_loso_min_k = m4$ari_loso_min,
    pac_boot_k = m4$pac_boot)
}))
write.csv(summary_df, file.path(OUTPUT_DIR, "variant_summary.csv"), row.names = FALSE)
print(summary_df, digits = 3)

write.csv(bind_rows(lapply(variants, `[[`, "metrics")),
          file.path(OUTPUT_DIR, "merge_metrics_by_variant.csv"), row.names = FALSE)
jac <- bind_rows(lapply(names(variants), function(v)
  variants[[v]]$jaccard %>% mutate(domain = if (v == "harmony_current")
    unname(DOMAIN_NAMES[as.character(domain)]) else paste0("D", domain))))
write.csv(jac, file.path(OUTPUT_DIR, paste0("domain_stability_k", MERGE_K, ".csv")), row.names = FALSE)
print(jac, digits = 3)

saveRDS(data.frame(sample = samples_all, condition = cond_all, celltype = celltype_all,
                   harmony_cluster = as.character(clusters_harm), harmony_domain = domain_named("harmony_current"),
                   noharm_cluster = as.character(clusters_pca), noharm_domain = domain_named("no_harmony"),
                   noharm_dims30_domain = domain_named("no_harmony_dims30")),
        file.path(OUTPUT_DIR, "cell_assignments_no_harmony.rds"))

# --- 6. Plots ------------------------------------------------------------------
# Spatial maps: current domains vs no-Harmony domains. No-Harmony domains are
# coloured by their best composition match when that mapping is one-to-one.
nh <- name_match[name_match$variant == "no_harmony", ]
nh_colors <- if (!any(duplicated(nh$best_composition_match))) {
  setNames(DOMAIN_COLORS[nh$best_composition_match], nh$domain)
} else setNames(c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")[seq_len(nrow(nh))], nh$domain)
nh_labels <- setNames(sprintf("%s (~%s)", nh$domain, nh$best_composition_match), nh$domain)

spatial_df <- data.frame(x = coords[, 1], y = coords[, 2], sample = factor(samples_all, levels = names(SAMPLE_COLORS)),
                         current = domain_named("harmony_current"), no_harmony = domain_named("no_harmony"))
p_cur <- ggplot(spatial_df, aes(x, y, color = current)) +
  geom_point(size = 0.05, stroke = 0) + scale_color_manual(values = DOMAIN_COLORS, name = "Current (Harmony)") +
  facet_wrap(~ sample, nrow = 1) + coord_fixed() + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3)))
p_nh <- ggplot(spatial_df, aes(x, y, color = no_harmony)) +
  geom_point(size = 0.05, stroke = 0) + scale_color_manual(values = nh_colors, labels = nh_labels, name = "No Harmony") +
  facet_wrap(~ sample, nrow = 1) + coord_fixed() + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3)))
ggsave(file.path(OUTPUT_DIR, "spatial_domains_harmony_vs_noharmony.png"), p_cur / p_nh,
       width = 24, height = 8, dpi = 200)

# iLISI distributions
p_lisi <- data.frame(embedding = rep(c("Harmony", "PCA (no Harmony)"), each = N_LISI_QUERY),
                     ilisi = c(lisi$harmony, lisi$pca)) %>%
  ggplot(aes(x = ilisi, fill = embedding)) +
  geom_histogram(bins = 60, alpha = 0.6, position = "identity") +
  labs(x = sprintf("iLISI over tumours (%d neighbours; 1 = single tumour, 7 = fully mixed)", LISI_K),
       y = "Cells", fill = NULL) +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "ilisi_harmony_vs_noharmony.pdf"), p_lisi, width = 7, height = 4)

# Cross-tab: where each current domain's cells go without Harmony
p_cross <- as.data.frame(prop.table(table(current = domain_named("harmony_current"),
                                          no_harmony = domain_named("no_harmony")), 1),
                         responseName = "frac") %>%
  ggplot(aes(no_harmony, current, fill = frac)) +
  geom_tile(color = "white") + geom_text(aes(label = sprintf("%.2f", frac)), size = 3) +
  scale_fill_gradient(low = "white", high = "#08519C", limits = c(0, 1), name = "Fraction of\ncurrent domain") +
  scale_x_discrete(labels = nh_labels) +
  labs(x = "No-Harmony domain (best composition match)", y = "Current domain") +
  theme_publication() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "domain_crosstab_harmony_vs_noharmony.pdf"), p_cross, width = 7, height = 4.5)

cat("Outputs saved to:", OUTPUT_DIR, "\n")
