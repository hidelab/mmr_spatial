# 02_08_NicheAssay.R
# Do cell-type neighbourhood niches recover the merged Banksy domains?
# 02_05-02_07 showed that the k = MERGE_K Banksy domains change with the Harmony
# seed (domain ARI 0.25-0.44) and that without Harmony the domains are dominated
# by single tumours. Niches built from the Tier 1 cell-type labels need neither
# an expression embedding nor batch correction. As in Seurat's BuildNicheAssay(),
# each cell is described by the cell-type counts among its NICHE_K nearest cells
# (itself included; neighbours searched within its own section), the counts are
# z-scored per cell type (clipped at 10, as ScaleData) and clustered by k-means.
# BuildNicheAssay() needs an FOV, which this object lacks, so the steps are
# re-implemented here. k-means uses the MacQueen algorithm, which reached the
# same objective as the default Hartigan-Wong ~3x faster and without its
# convergence warnings at this size, and 150 random starts instead of
# BuildNicheAssay's 30 (see KMEANS_NSTART).
#
# For k = 2..10:
#   (1) Determinism       pairwise ARI between k-means seeds
#   (2) Tumour stability  leave-one-tumour-out refits, ARI on the remaining cells
#   (3) Fit               between-niche share of variance, silhouette (subsample)
#   (4) Tumour balance    largest share of any niche coming from one tumour
#   (5) Agreement         ARI / NMI vs the current domains (02_03)
# At k = MERGE_K: niche-to-domain matching, cell-type composition, cross-tab with
# the current domains and sensitivity to the neighbourhood size. Spatial maps for
# k in SPATIAL_K.
#
# Read-only with respect to existing results: the Seurat object is never saved
# and all outputs go to a folder used only by this script.

library(Seurat)
library(FNN)
library(Matrix)
library(parallel)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(cluster)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "niche_assay")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
NICHE_K       <- 20               # neighbourhood size (BuildNicheAssay default)
NICHE_K_SENS  <- c(10, 20, 50)    # neighbourhood sizes compared at k = MERGE_K
K_RANGE       <- 2:10
KMEANS_SEEDS  <- 1:5
KMEANS_NSTART <- 150              # BuildNicheAssay uses 30, which left some seeds in a
                                  # worse local optimum at k = 5-6
KMEANS_ITER   <- 100
N_SIL         <- 5000             # cells used for the silhouette
SPATIAL_K     <- 3:6              # niche numbers shown on the spatial maps
N_CORES       <- min(6, detectCores())
FALLBACK_COLORS <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                     "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#B2DF8A")
set.seed(42)

# --- Helpers -------------------------------------------------------------------
ari <- function(a, b) mclust::adjustedRandIndex(a, b)

# Normalized mutual information (arithmetic-mean normalization) of two partitions.
nmi <- function(a, b) {
  p  <- table(a, b) / length(a)
  pa <- rowSums(p); pb <- colSums(p)
  nz <- p > 0
  mi <- sum(p[nz] * log(p[nz] / outer(pa, pb)[nz]))
  h  <- function(q) -sum(q[q > 0] * log(q[q > 0]))
  2 * mi / (h(pa) + h(pb))
}

# Cell-type counts among each cell's `k` nearest cells (itself included),
# searched within its own section; z-scored per cell type and clipped at 10.
niche_features <- function(k) {
  counts <- matrix(0, length(celltype_all), nlevels(celltype_all),
                   dimnames = list(NULL, levels(celltype_all)))
  for (s in unique(samples_all)) {
    idx <- which(samples_all == s)
    nn  <- cbind(seq_along(idx), get.knn(coords[idx, ], k = k - 1)$nn.index)
    counts[idx, ] <- as.matrix(sparseMatrix(i = rep(seq_along(idx), k),
                                            j = as.integer(celltype_all[idx])[nn], x = 1,
                                            dims = c(length(idx), nlevels(celltype_all))))
  }
  x <- scale(counts)
  x[x > 10] <- 10
  x
}

# k-means with a fixed seed; niches renumbered by decreasing size so that the
# labels do not depend on k-means' arbitrary ordering.
fit_niches <- function(x, k, seed) {
  set.seed(seed)
  fit <- kmeans(x, centers = k, nstart = KMEANS_NSTART, iter.max = KMEANS_ITER,
                algorithm = "MacQueen")
  list(cluster = match(fit$cluster, order(-tabulate(fit$cluster, k))),
       tot_withinss = fit$tot.withinss, betweenss = fit$betweenss, totss = fit$totss)
}

# Name each niche after the current domain holding most of its cells; colour it
# with that domain's colour when the matching is one-to-one.
niche_palette <- function(niche) {
  tab   <- table(niche, current_names)
  match <- setNames(colnames(tab)[apply(tab, 1, which.max)], rownames(tab))
  cols  <- if (!any(duplicated(match))) setNames(DOMAIN_COLORS[match], names(match))
           else setNames(FALLBACK_COLORS[seq_along(match)], names(match))
  list(colors = cols, labels = setNames(sprintf("%s (~%s)", names(match), match), names(match)))
}

# --- 1. Load object (read-only) ------------------------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
stopifnot("spatial" %in% Reductions(seu), !anyNA(seu$Tier1_celltype))
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)

samples_all   <- as.character(seu$sample)
cond_all      <- as.character(seu$condition)
celltype_all  <- factor(seu$Tier1_celltype)
coords        <- Embeddings(seu, "spatial")
domains_curr  <- as.character(seu@meta.data[[COL_DOMAIN]])
current_names <- unname(DOMAIN_NAMES[domains_curr])
rm(seu); invisible(gc())
samp <- intersect(names(SAMPLE_COLORS), unique(samples_all))

# --- 2. Niche features and k-means ---------------------------------------------
cat(sprintf("Building niche features (%d nearest cells)...\n", NICHE_K))
feat <- niche_features(NICHE_K)

cat(sprintf("k-means: %d values of k x %d seeds...\n", length(K_RANGE), length(KMEANS_SEEDS)))
grid <- expand.grid(seed = KMEANS_SEEDS, k = K_RANGE)
fits <- mclapply(seq_len(nrow(grid)), function(i) fit_niches(feat, grid$k[i], grid$seed[i]),
                 mc.cores = N_CORES)

# Reference partition for each k: the seed with the lowest within-niche sum of squares
ref <- lapply(K_RANGE, function(k) {
  f <- fits[grid$k == k]
  f[[which.min(sapply(f, `[[`, "tot_withinss"))]]
})
names(ref) <- paste0("k", K_RANGE)

cat("Leave-one-tumour-out refits...\n")
loso_grid <- expand.grid(sample = samp, k = K_RANGE, stringsAsFactors = FALSE)
loso_ari  <- unlist(mclapply(seq_len(nrow(loso_grid)), function(i) {
  keep <- samples_all != loso_grid$sample[i]
  fit  <- fit_niches(feat[keep, ], loso_grid$k[i], KMEANS_SEEDS[1])
  ari(fit$cluster, ref[[paste0("k", loso_grid$k[i])]]$cluster[keep])
}, mc.cores = N_CORES))
loso_grid$ari <- loso_ari

# --- 3. Metrics by k -----------------------------------------------------------
sil_idx  <- sample(nrow(feat), N_SIL)
sil_dist <- dist(feat[sil_idx, ])
metrics <- bind_rows(lapply(K_RANGE, function(k) {
  f    <- fits[grid$k == k]
  r    <- ref[[paste0("k", k)]]
  pair <- combn(length(f), 2, function(p) ari(f[[p[1]]]$cluster, f[[p[2]]]$cluster))
  lo   <- loso_grid$ari[loso_grid$k == k]
  tumour_share <- prop.table(table(r$cluster, samples_all), 1)
  data.frame(
    k = k,
    var_explained         = r$betweenss / r$totss,
    silhouette            = mean(silhouette(r$cluster[sil_idx], sil_dist)[, "sil_width"]),
    ari_seed_mean         = mean(pair),
    ari_seed_min          = min(pair),
    ari_loso_mean         = mean(lo),
    ari_loso_min          = min(lo),
    max_single_tumour_share = max(tumour_share),
    min_niche_size        = min(tabulate(r$cluster, k)),
    ari_vs_current        = ari(r$cluster, domains_curr),
    nmi_vs_current        = nmi(r$cluster, domains_curr))
}))
write.csv(metrics, file.path(OUTPUT_DIR, "niche_metrics_by_k.csv"), row.names = FALSE)
write.csv(loso_grid, file.path(OUTPUT_DIR, "niche_loso_ari.csv"), row.names = FALSE)
print(metrics, digits = 3)

# --- 4. Niches at k = MERGE_K vs the current domains ---------------------------
niche_k <- paste0("N", ref[[paste0("k", MERGE_K)]]$cluster)

# Best-matching current domain by cell-type composition (cosine) and by overlap
ct_comp <- function(dom) prop.table(table(dom, celltype_all), 1)
comp_current <- ct_comp(current_names)
comp_niche   <- ct_comp(niche_k)
cosine  <- (comp_niche %*% t(comp_current)) /
  outer(sqrt(rowSums(comp_niche^2)), sqrt(rowSums(comp_current^2)))
overlap <- prop.table(table(niche_k, current_names), 1)[rownames(comp_niche), , drop = FALSE]
name_match <- data.frame(
  niche = rownames(comp_niche),
  n_cells = as.integer(table(niche_k)[rownames(comp_niche)]),
  best_composition_match = colnames(cosine)[apply(cosine, 1, which.max)],
  composition_cosine     = apply(cosine, 1, max),
  best_overlap_match     = colnames(overlap)[apply(overlap, 1, which.max)],
  overlap_fraction       = apply(overlap, 1, max))
write.csv(name_match, file.path(OUTPUT_DIR, sprintf("niche_name_matching_k%d.csv", MERGE_K)), row.names = FALSE)
print(name_match, digits = 3)

celltype_comp <- bind_rows(
  as.data.frame(comp_current, responseName = "fraction") %>%
    rename(group = dom, celltype = celltype_all) %>% mutate(partition = "Current domains", .before = 1),
  as.data.frame(comp_niche, responseName = "fraction") %>%
    rename(group = dom, celltype = celltype_all) %>% mutate(partition = sprintf("Niches (k = %d)", MERGE_K), .before = 1))
write.csv(celltype_comp, file.path(OUTPUT_DIR, sprintf("niche_celltype_composition_k%d.csv", MERGE_K)),
          row.names = FALSE)

sample_comp <- as.data.frame(table(niche = niche_k, sample = samples_all), responseName = "n_cells") %>%
  group_by(niche) %>% mutate(frac_of_niche = n_cells / sum(n_cells)) %>%
  group_by(sample) %>% mutate(frac_of_sample = n_cells / sum(n_cells)) %>% ungroup()
write.csv(sample_comp, file.path(OUTPUT_DIR, sprintf("niche_sample_composition_k%d.csv", MERGE_K)),
          row.names = FALSE)

crosstab <- as.data.frame(table(current_domain = current_names, niche = niche_k), responseName = "n_cells")
write.csv(crosstab, file.path(OUTPUT_DIR, sprintf("niche_crosstab_vs_current_k%d.csv", MERGE_K)),
          row.names = FALSE)

# Sensitivity to the neighbourhood size (seed KMEANS_SEEDS[1] for each size)
cat("Neighbourhood-size sensitivity...\n")
sens_fits <- lapply(NICHE_K_SENS, function(nk) {
  if (nk == NICHE_K) return(ref[[paste0("k", MERGE_K)]])
  fit_niches(niche_features(nk), MERGE_K, KMEANS_SEEDS[1])
})
sensitivity <- data.frame(
  niche_k        = NICHE_K_SENS,
  ari_vs_niche_k = sapply(sens_fits, function(f) ari(f$cluster, ref[[paste0("k", MERGE_K)]]$cluster)),
  ari_vs_current = sapply(sens_fits, function(f) ari(f$cluster, domains_curr)),
  nmi_vs_current = sapply(sens_fits, function(f) nmi(f$cluster, domains_curr)))
names(sensitivity)[2] <- sprintf("ari_vs_niche_k%d", NICHE_K)
write.csv(sensitivity, file.path(OUTPUT_DIR, sprintf("niche_neighbourhood_size_k%d.csv", MERGE_K)),
          row.names = FALSE)
print(sensitivity, digits = 3)

# Cell-level niche assignments for every k
assign_df <- data.frame(sample = samples_all, condition = cond_all, celltype = as.character(celltype_all),
                        current_domain = current_names,
                        sapply(ref, function(r) paste0("N", r$cluster)))
names(assign_df)[-(1:4)] <- paste0("niche_", names(ref))
saveRDS(assign_df, file.path(OUTPUT_DIR, "cell_assignments_niche.rds"))

# --- 5. Plots ------------------------------------------------------------------
metric_labels <- c(ari_seed_min = "ARI between k-means seeds (min)",
                   ari_loso_min = "ARI, leave-one-tumour-out (min)",
                   silhouette = "Silhouette (subsample)",
                   var_explained = "Between-niche share of variance",
                   max_single_tumour_share = "Largest single-tumour share of a niche",
                   ari_vs_current = "ARI vs current domains")
p_metrics <- metrics %>%
  pivot_longer(all_of(names(metric_labels)), names_to = "metric") %>%
  mutate(metric = factor(metric_labels[metric], levels = metric_labels)) %>%
  ggplot(aes(x = k, y = value)) +
  geom_vline(xintercept = MERGE_K, linetype = "dashed", color = "grey50") +
  geom_line() + geom_point(size = 1.5) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(x = "Number of niches (k)", y = NULL,
       title = sprintf("Cell-type niches (%d nearest cells)", NICHE_K)) +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "niche_metrics_by_k.pdf"), p_metrics, width = 10, height = 8)

# Cell-type composition of the current domains and the niches
p_comp <- celltype_comp %>%
  mutate(group = factor(group, levels = rev(c(DOMAIN_ORDER, sort(unique(niche_k)))))) %>%
  ggplot(aes(x = celltype, y = group, fill = fraction)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", fraction)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#08519C", name = "Fraction\nof cells") +
  facet_wrap(~ partition, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, sprintf("niche_celltype_composition_k%d.pdf", MERGE_K)), p_comp,
       width = 12, height = 6)

# Cross-tab (row-normalized: where each current domain's cells go)
p_cross <- crosstab %>%
  group_by(current_domain) %>% mutate(frac = n_cells / sum(n_cells)) %>% ungroup() %>%
  mutate(current_domain = factor(current_domain, levels = rev(DOMAIN_ORDER))) %>%
  ggplot(aes(x = niche, y = current_domain, fill = frac)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", frac)), size = 3) +
  scale_fill_gradient(low = "white", high = "#08519C", limits = c(0, 1), name = "Fraction of\ncurrent domain") +
  labs(x = sprintf("Niche (k = %d)", MERGE_K), y = "Current domain (02_03)") +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, sprintf("niche_crosstab_vs_current_k%d.pdf", MERGE_K)), p_cross,
       width = 6, height = 4)

# Spatial maps: current domains, then niches for each k in SPATIAL_K
spatial_row <- function(labels, cols, legend_labels, title) {
  data.frame(x = coords[, 1], y = coords[, 2], group = labels,
             sample = factor(samples_all, levels = samp)) %>%
    ggplot(aes(x, y, color = group)) +
    geom_point(size = 0.2, stroke = 0) +
    scale_color_manual(values = cols, labels = legend_labels, name = title) +
    facet_wrap(~ sample, nrow = 1) + coord_fixed() + theme_void() +
    guides(color = guide_legend(override.aes = list(size = 3)))
}
rows <- c(list(spatial_row(current_names, DOMAIN_COLORS, setNames(names(DOMAIN_COLORS), names(DOMAIN_COLORS)),
                           "Current domains")),
          lapply(SPATIAL_K, function(k) {
            lab <- paste0("N", ref[[paste0("k", k)]]$cluster)
            pal <- niche_palette(lab)
            spatial_row(lab, pal$colors, pal$labels, sprintf("Niches, k = %d", k))
          }))
p_spatial <- wrap_plots(rows, ncol = 1) +
  plot_annotation(title = sprintf("Current domains vs cell-type niches (%d nearest cells)", NICHE_K))
ggsave(file.path(OUTPUT_DIR, "spatial_niches_vs_current.png"), p_spatial,
       width = 24, height = 4 * length(rows), dpi = 200)

cat("Outputs saved to:", OUTPUT_DIR, "\n")
