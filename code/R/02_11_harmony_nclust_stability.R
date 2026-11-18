# 02_11_harmony_nclust_stability.R
# Does a smaller number of Harmony clusters (nclust) make the Banksy domains
# reproducible across Harmony seeds, and does it over-correct?
#
# 02_09/02_10 showed that the domains change with the Harmony seed at every
# clustering resolution. Harmony's default nclust is min(N / 30, 100) = 100 here,
# and its k-means initialisation with 100 centres on ~330k cells stops at the
# Quick-TRANSfer step limit; fewer centres may initialise (and converge) more
# reproducibly, but fewer, broader Harmony clusters also risk over-correcting
# real between-tumour differences. For every nclust in NCLUST_VALUES Harmony is
# run with each seed in HARMONY_SEEDS (other settings as in 02_01), then SNN on
# Harmony 1:BANKSY_NPC_USE, Louvain at RESOLUTION, and both domain recipes of
# 02_10. Reported per nclust:
#   Convergence      Harmony iterations, k-means warnings, run time
#   Reproducibility  kNN overlap of the embeddings, cluster and domain ARI
#                    between seeds, cells given the same domain by all seeds
#   Correction       iLISI over tumours (as 02_07) and cell-type purity of the
#                    neighbourhoods, with the uncorrected PCA as reference;
#                    Tumor Core share of every sample
#   Effect of nclust domain ARI between nclust settings (all seed pairs)
#
# Read-only: nothing is saved except CSV tables and plots, all in a folder used
# only by this script (no RDS output).

library(Seurat)
library(harmony)
library(Matrix)
library(FNN)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "harmony_nclust_stability")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
NCLUST_VALUES    <- c(8, 30, NA)        # NA = Harmony default, min(N / 30, 100)
HARMONY_SEEDS    <- c(42, 1, 2, 3, 4)   # 42 = the seed used in 02_01
RESOLUTION       <- 0.5                 # best worst-case domain ARI in 02_10
MIN_CLUSTER_FRAC <- 0.005               # minsize: clusters below this fraction of cells are dissolved
RECIPES  <- c(current_k4 = "unweighted", minsize_k5 = "minsize")
RECIPE_K <- c(current_k4 = MERGE_K, minsize_k5 = 5)
RECIPE_LABELS <- c(current_k4 = sprintf("Unweighted Ward.D2, k = %d (02_03)", MERGE_K),
                   minsize_k5 = sprintf("Min-size %.1f%% + SNN, Ward.D2, k = 5", 100 * MIN_CLUSTER_FRAC))
N_QUERY <- 5000                         # cells at which iLISI, purity and kNN overlap are evaluated
LISI_K  <- 90                           # neighbours for iLISI and purity (as 02_07)
KNN_K   <- 30                           # neighbours for the kNN overlap between embeddings
set.seed(7)

# --- Helpers -------------------------------------------------------------------
ari <- function(a, b) mclust::adjustedRandIndex(a, b)
nclust_label <- function(n) ifelse(is.na(n), "default (100)", as.character(n))

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

# iLISI (Korsunsky et al. 2019) approximated as in 02_07: inverse Simpson index
# of tumour labels among each query cell's LISI_K nearest neighbours (1 = single
# tumour, 7 = fully mixed). Purity: share of those neighbours with the query
# cell's Tier 1 cell type (lower than PCA = cell types blurred together).
neighbourhood_metrics <- function(emb, query) {
  nn <- get.knnx(emb, emb[query, , drop = FALSE], k = LISI_K)$nn.index
  data.frame(
    ilisi  = apply(nn, 1, function(idx) { p <- table(samples[idx]) / LISI_K; 1 / sum(p^2) }),
    purity = rowMeans(matrix(celltype[nn], nrow = length(query)) == celltype[query]))
}

# Mean Jaccard overlap of the KNN_K-nearest-neighbour sets of the query cells in
# two embeddings (as 02_06).
knn_overlap <- function(nn_a, nn_b) {
  mean(vapply(seq_len(nrow(nn_a)), function(i) {
    length(intersect(nn_a[i, ], nn_b[i, ])) / length(union(nn_a[i, ], nn_b[i, ]))
  }, numeric(1)))
}

# --- 1. Load object (read-only) ------------------------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
samples   <- as.character(seu$sample)
samp_fac  <- factor(samples, levels = intersect(names(SAMPLE_COLORS), unique(samples)))
celltype  <- as.character(seu$Tier1_celltype)
is_tumour <- celltype == "Tumor"
query     <- sample(ncol(seu), N_QUERY)

# Reference: uncorrected PCA
pca_nm <- neighbourhood_metrics(Embeddings(seu, "pca")[, 1:BANKSY_NPC_USE], query)

# --- 2. Harmony (per nclust and seed) -> SNN -> Louvain -> domains -------------
runs <- list()
for (nc in NCLUST_VALUES) {
  for (s in HARMONY_SEEDS) {
    id <- sprintf("nclust_%s__seed%d", nclust_label(nc), s)
    cat(sprintf("%s: Harmony...\n", id))
    msgs <- character(0); warns <- character(0)
    harmony_args <- list(group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca",
                         dims.use = seq_len(BANKSY_NPCS), reduction.save = "harmony_run", verbose = TRUE)
    if (!is.na(nc)) harmony_args$nclust <- nc
    set.seed(s)
    t_harmony <- system.time(
      seu <- withCallingHandlers(do.call(RunHarmony, c(list(seu), harmony_args)),
        message = function(m) { msgs  <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
        warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    )[["elapsed"]]
    conv <- regmatches(msgs, regexpr("converged after [0-9]+ iteration", msgs))

    set.seed(42)
    seu <- FindNeighbors(seu, reduction = "harmony_run", dims = 1:BANKSY_NPC_USE,
                         graph.name = c("run_nn", "run_snn"), verbose = FALSE)
    seu <- FindClusters(seu, graph.name = "run_snn", resolution = RESOLUTION,
                        cluster.name = "run_clusters", verbose = FALSE)
    emb <- Embeddings(seu, "harmony_run")[, 1:BANKSY_NPC_USE]
    g   <- as(seu[["run_snn"]], "dgCMatrix")
    cls <- as.character(seu$run_clusters)

    runs[[id]] <- list(
      nclust = nclust_label(nc), seed = s, clusters = cls,
      domains = lapply(setNames(names(RECIPES), names(RECIPES)), function(rc)
        merge_domains(RECIPES[[rc]], emb, cls, g, RECIPE_K[[rc]])),
      harmony_seconds = t_harmony,
      converged_after = if (length(conv)) as.integer(sub("\\D+(\\d+).*", "\\1", conv[1])) else NA_integer_,
      n_kmeans_warnings = sum(grepl("Quick-TRANSfer", warns)),
      other_warnings = paste(unique(warns[!grepl("Quick-TRANSfer", warns)]), collapse = " | "),
      knn = get.knnx(emb, emb[query, , drop = FALSE], k = KNN_K + 1)$nn.index[, -1],
      nm  = neighbourhood_metrics(emb, query))
    cat(sprintf("  %.0f s, converged after %s iterations, %d k-means warnings, %d clusters\n",
                t_harmony, runs[[id]]$converged_after, runs[[id]]$n_kmeans_warnings,
                length(unique(cls))))

    seu[["harmony_run"]] <- NULL
    seu[["run_nn"]]  <- NULL
    seu[["run_snn"]] <- NULL
    rm(emb, g); invisible(gc())
  }
}
run_info <- bind_rows(lapply(names(runs), function(id) with(runs[[id]], data.frame(
  run = id, nclust = nclust, seed = seed, n_clusters = length(unique(clusters)),
  harmony_seconds = harmony_seconds, converged_after = converged_after,
  n_kmeans_warnings = n_kmeans_warnings, other_warnings = other_warnings,
  ilisi_median = median(nm$ilisi), ilisi_median_tumour_cells = median(nm$ilisi[is_tumour[query]]),
  celltype_purity_mean = mean(nm$purity)))))
write.csv(run_info, file.path(OUTPUT_DIR, "run_info.csv"), row.names = FALSE)
nclust_levels <- unique(run_info$nclust)

# --- 3. Reproducibility across seeds (within each nclust) ----------------------
pairwise <- bind_rows(lapply(nclust_levels, function(nc) {
  ids <- names(runs)[sapply(runs, `[[`, "nclust") == nc]
  bind_rows(lapply(combn(ids, 2, simplify = FALSE), function(p) {
    a <- runs[[p[1]]]; b <- runs[[p[2]]]
    out <- data.frame(nclust = nc, seed_a = a$seed, seed_b = b$seed,
                      knn_jaccard = knn_overlap(a$knn, b$knn),
                      cluster_ari = ari(a$clusters, b$clusters))
    for (rc in names(RECIPES)) {
      out[[paste0("domain_ari_", rc)]] <- if (anyNA(a$domains[[rc]]) || anyNA(b$domains[[rc]])) NA
                                          else ari(a$domains[[rc]], b$domains[[rc]])
    }
    out
  }))
}))
write.csv(pairwise, file.path(OUTPUT_DIR, "pairwise_seed_agreement.csv"), row.names = FALSE)

# Cells given the same (matched) domain by every seed
unanimity <- bind_rows(lapply(nclust_levels, function(nc) {
  ids <- names(runs)[sapply(runs, `[[`, "nclust") == nc]
  bind_rows(lapply(names(RECIPES), function(rc) {
    doms <- lapply(ids, function(id) runs[[id]]$domains[[rc]])
    if (any(sapply(doms, anyNA))) return(data.frame(nclust = nc, recipe = rc, frac_cells_unanimous = NA))
    matched <- sapply(doms, match_to_ref, ref = doms[[1]])
    data.frame(nclust = nc, recipe = rc, frac_cells_unanimous = mean(rowSums(matched == matched[, 1]) == ncol(matched)))
  }))
}))
write.csv(unanimity, file.path(OUTPUT_DIR, "per_cell_unanimity.csv"), row.names = FALSE)

# --- 4. Does nclust change the domains? (all seed pairs across settings) -------
between <- bind_rows(lapply(combn(nclust_levels, 2, simplify = FALSE), function(p) {
  ia <- names(runs)[sapply(runs, `[[`, "nclust") == p[1]]
  ib <- names(runs)[sapply(runs, `[[`, "nclust") == p[2]]
  bind_rows(lapply(ia, function(a) bind_rows(lapply(ib, function(b) {
    out <- data.frame(nclust_a = p[1], nclust_b = p[2], seed_a = runs[[a]]$seed, seed_b = runs[[b]]$seed)
    for (rc in names(RECIPES)) {
      da <- runs[[a]]$domains[[rc]]; db <- runs[[b]]$domains[[rc]]
      out[[paste0("domain_ari_", rc)]] <- if (anyNA(da) || anyNA(db)) NA else ari(da, db)
    }
    out
  }))))
}))
write.csv(between, file.path(OUTPUT_DIR, "domain_agreement_between_nclust.csv"), row.names = FALSE)

# --- 5. Tumor Core share per sample --------------------------------------------
# Tumor Core = the domain with the highest tumour-cell fraction in that run
tumor_core <- bind_rows(lapply(names(runs), function(id) bind_rows(lapply(names(RECIPES), function(rc) {
  d <- runs[[id]]$domains[[rc]]
  if (anyNA(d)) return(NULL)
  tc <- names(which.max(tapply(is_tumour, d, mean)))
  data.frame(nclust = runs[[id]]$nclust, seed = runs[[id]]$seed, recipe = rc, sample = levels(samp_fac),
             tumour_frac_of_core = mean(is_tumour[d == tc]),
             core_frac_of_sample = as.numeric(tapply(d == tc, samp_fac, mean)))
}))))
write.csv(tumor_core, file.path(OUTPUT_DIR, "tumor_core_by_sample.csv"), row.names = FALSE)

# --- 6. Summary ----------------------------------------------------------------
summary_df <- run_info %>%
  group_by(nclust) %>%
  summarise(harmony_seconds_median = median(harmony_seconds),
            converged_after = paste(sort(unique(converged_after)), collapse = ","),
            kmeans_warnings_median = median(n_kmeans_warnings),
            n_clusters = paste(range(n_clusters), collapse = "-"),
            ilisi_median = median(ilisi_median), ilisi_tumour = median(ilisi_median_tumour_cells),
            celltype_purity = median(celltype_purity_mean), .groups = "drop") %>%
  left_join(pairwise %>% group_by(nclust) %>%
              summarise(knn_jaccard = mean(knn_jaccard), cluster_ari = mean(cluster_ari),
                        across(starts_with("domain_ari_"), list(mean = ~mean(.x), min = ~min(.x)),
                               .names = "{.col}_{.fn}"), .groups = "drop"), by = "nclust") %>%
  left_join(unanimity %>% pivot_wider(names_from = recipe, values_from = frac_cells_unanimous,
                                      names_prefix = "unanimous_"), by = "nclust") %>%
  mutate(nclust = factor(nclust, levels = nclust_levels)) %>% arrange(nclust)
write.csv(summary_df, file.path(OUTPUT_DIR, "nclust_summary.csv"), row.names = FALSE)

cat(sprintf("\nReference, uncorrected PCA: iLISI median %.2f (tumour cells %.2f), cell-type purity %.3f\n",
            median(pca_nm$ilisi), median(pca_nm$ilisi[is_tumour[query]]), mean(pca_nm$purity)))
cat("\nSummary by nclust:\n")
print(as.data.frame(summary_df), digits = 3)

cat("\nDomain ARI between nclust settings (mean over seed pairs):\n")
print(as.data.frame(between %>% group_by(nclust_a, nclust_b) %>%
                      summarise(across(starts_with("domain_ari_"), mean), .groups = "drop")), digits = 3)

cat("\nTumor Core share of each sample (median over seeds [min-max]):\n")
print(as.data.frame(tumor_core %>%
  group_by(recipe, nclust, sample) %>%
  summarise(v = sprintf("%.2f [%.2f-%.2f]", median(core_frac_of_sample),
                        min(core_frac_of_sample), max(core_frac_of_sample)), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = v) %>%
  arrange(recipe, factor(nclust, levels = nclust_levels))))

# --- 7. Plots ------------------------------------------------------------------
nclust_fac <- function(x) factor(x, levels = nclust_levels)

p_ari <- pairwise %>%
  pivot_longer(c(knn_jaccard, cluster_ari, starts_with("domain_ari_")), names_to = "level", values_to = "value") %>%
  mutate(level = factor(c(knn_jaccard = "Embedding kNN overlap (Jaccard)", cluster_ari = "Clusters (ARI)",
                          setNames(paste(RECIPE_LABELS, "(ARI)"), paste0("domain_ari_", names(RECIPE_LABELS))))[level],
                        levels = c("Embedding kNN overlap (Jaccard)", "Clusters (ARI)", paste(RECIPE_LABELS, "(ARI)"))),
         nclust = nclust_fac(nclust)) %>%
  ggplot(aes(x = nclust, y = value)) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 1.5, color = "grey40") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#D95F02", linewidth = 0.4) +
  facet_wrap(~ level, nrow = 1, labeller = labeller(level = label_wrap_gen(25))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Harmony nclust", y = "Agreement between seeds (pairs)",
       title = sprintf("Reproducibility across Harmony seeds (resolution %.1f)", RESOLUTION)) +
  theme_publication()

nm_df <- bind_rows(
  data.frame(nclust = "PCA (no Harmony)", ilisi = pca_nm$ilisi, purity = pca_nm$purity),
  bind_rows(lapply(runs, function(r) data.frame(nclust = r$nclust, ilisi = r$nm$ilisi, purity = r$nm$purity)))) %>%
  mutate(nclust = factor(nclust, levels = c("PCA (no Harmony)", nclust_levels)))
p_lisi <- ggplot(nm_df, aes(x = nclust, y = ilisi)) +
  geom_boxplot(outlier.shape = NA, fill = "#A8D8E8") +
  labs(x = NULL, y = sprintf("iLISI over tumours (%d neighbours)", LISI_K), title = "Tumour mixing") +
  theme_publication() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
p_pur <- nm_df %>% group_by(nclust) %>% summarise(purity = mean(purity), .groups = "drop") %>%
  ggplot(aes(x = nclust, y = purity)) +
  geom_col(fill = "#4A9FBC") +
  labs(x = NULL, y = "Mean share of neighbours with the same cell type", title = "Cell-type purity") +
  theme_publication() + theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_tc <- tumor_core %>%
  mutate(nclust = nclust_fac(nclust), recipe = factor(RECIPE_LABELS[recipe], levels = RECIPE_LABELS),
         sample = factor(sample, levels = levels(samp_fac))) %>%
  ggplot(aes(x = sample, y = core_frac_of_sample, color = sample)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.5) +
  scale_color_manual(values = SAMPLE_COLORS, guide = "none") +
  facet_grid(recipe ~ nclust, labeller = labeller(recipe = label_wrap_gen(25),
                                                  nclust = function(x) paste("nclust", x))) +
  labs(x = NULL, y = "Tumor Core share of sample (one point per seed)",
       title = "Tumor Core per sample across Harmony seeds") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "reproducibility_by_nclust.pdf"), p_ari, width = 14, height = 5)
ggsave(file.path(OUTPUT_DIR, "correction_by_nclust.pdf"), p_lisi + p_pur, width = 10, height = 5)
ggsave(file.path(OUTPUT_DIR, "tumor_core_by_sample.pdf"), p_tc, width = 12, height = 6.5)

cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")
