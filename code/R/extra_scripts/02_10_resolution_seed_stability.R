# 02_10_resolution_seed_stability.R
# Does a lower clustering resolution make the clusters and merged domains
# reproducible across Harmony seeds?
#
# 02_09 showed that with resolution BANKSY_RESOLUTION (0.8) both the current
# merge (unweighted Ward.D2, k = MERGE_K) and the minsize k = 5 merge change
# substantially with the Harmony seed. Here Harmony is run once per seed (02_01
# settings), the SNN graph is built once per seed, and Louvain clustering is
# repeated on that graph at each resolution in RESOLUTIONS. Every clustering is
# merged with both recipes and compared across seeds:
#   - pairwise ARI of clusters and of domains between seeds
#   - per-cell agreement: domains of every seed matched one-to-one (Hungarian)
#     to the first seed; fraction of cells given the same domain by all seeds
#   - Tumor Core (the domain with the highest tumour-cell fraction): share of
#     each sample's cells, per seed
#
# Read-only: nothing is saved except CSV tables and plots, all in a folder used
# only by this script (no RDS output).

library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "resolution_seed_stability")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
RESOLUTIONS      <- c(0.8, 0.5, 0.3, 0.1)
HARMONY_SEEDS    <- c(42, 1, 2, 3, 4)   # 42 = the seed used in 02_01
MIN_CLUSTER_FRAC <- 0.005               # minsize: clusters below this fraction of cells are dissolved
RECIPES  <- c(current_k4 = "unweighted", minsize_k5 = "minsize")
RECIPE_K <- c(current_k4 = MERGE_K, minsize_k5 = 5)
RECIPE_LABELS <- c(current_k4 = sprintf("Unweighted Ward.D2, k = %d (02_03)", MERGE_K),
                   minsize_k5 = sprintf("Min-size %.1f%% + SNN, Ward.D2, k = 5", 100 * MIN_CLUSTER_FRAC))

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

# --- 1. Load object (read-only) ------------------------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
samples   <- factor(seu$sample, levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))
is_tumour <- seu$Tier1_celltype == "Tumor"

# --- 2. Harmony per seed; Louvain at every resolution on that seed's graph -----
results <- list()
for (s in HARMONY_SEEDS) {
  name <- paste0("seed", s)
  cat(sprintf("Seed %d: Harmony + SNN graph...\n", s))
  set.seed(s)
  seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca",
                    dims.use = seq_len(BANKSY_NPCS), reduction.save = "harmony_seed", verbose = FALSE)
  set.seed(42)
  seu <- FindNeighbors(seu, reduction = "harmony_seed", dims = 1:BANKSY_NPC_USE,
                       graph.name = c("seed_nn", "seed_snn"), verbose = FALSE)
  emb <- Embeddings(seu, "harmony_seed")
  g   <- as(seu[["seed_snn"]], "dgCMatrix")

  for (r in RESOLUTIONS) {
    seu <- FindClusters(seu, graph.name = "seed_snn", resolution = r,
                        cluster.name = "seed_clusters", verbose = FALSE)
    cls <- as.character(seu$seed_clusters)
    results[[name]][[as.character(r)]] <- list(
      clusters = cls,
      domains  = lapply(setNames(names(RECIPES), names(RECIPES)), function(rc)
        merge_domains(RECIPES[[rc]], emb, cls, g, RECIPE_K[[rc]])))
    cat(sprintf("  resolution %.2f: %d clusters\n", r, length(unique(cls))))
  }
  seu[["harmony_seed"]] <- NULL
  seu[["seed_nn"]]  <- NULL
  seu[["seed_snn"]] <- NULL
  rm(emb, g); invisible(gc())
}
seed_names <- names(results)

# --- 3. Pairwise agreement between seeds ---------------------------------------
pairs <- combn(seed_names, 2, simplify = FALSE)
pairwise <- bind_rows(lapply(RESOLUTIONS, function(r) {
  rr <- as.character(r)
  bind_rows(lapply(pairs, function(p) {
    a <- results[[p[1]]][[rr]]; b <- results[[p[2]]][[rr]]
    out <- data.frame(resolution = r, seed_a = p[1], seed_b = p[2],
                      n_clusters_a = length(unique(a$clusters)), n_clusters_b = length(unique(b$clusters)),
                      cluster_ari = ari(a$clusters, b$clusters))
    for (rc in names(RECIPES)) {
      out[[paste0("domain_ari_", rc)]] <- if (anyNA(a$domains[[rc]]) || anyNA(b$domains[[rc]])) NA
                                          else ari(a$domains[[rc]], b$domains[[rc]])
    }
    out
  }))
}))
write.csv(pairwise, file.path(OUTPUT_DIR, "pairwise_seed_agreement.csv"), row.names = FALSE)

# --- 4. Per-cell agreement across all seeds ------------------------------------
# Domains of every seed are matched to the first seed; a cell is unanimous when
# all seeds give it the same (matched) domain.
unanimity <- bind_rows(lapply(RESOLUTIONS, function(r) {
  rr <- as.character(r)
  bind_rows(lapply(names(RECIPES), function(rc) {
    doms <- lapply(seed_names, function(sn) results[[sn]][[rr]]$domains[[rc]])
    if (any(sapply(doms, anyNA))) {
      return(data.frame(resolution = r, recipe = rc, frac_cells_unanimous = NA, mean_agreement_with_ref = NA))
    }
    ref     <- doms[[1]]
    matched <- sapply(doms, match_to_ref, ref = ref)
    agree   <- matched == matched[, 1]
    data.frame(resolution = r, recipe = rc,
               frac_cells_unanimous    = mean(rowSums(agree) == ncol(agree)),
               mean_agreement_with_ref = mean(agree[, -1]))
  }))
}))
write.csv(unanimity, file.path(OUTPUT_DIR, "per_cell_unanimity.csv"), row.names = FALSE)

# --- 5. Tumor Core share per sample --------------------------------------------
# Tumor Core = the domain with the highest tumour-cell fraction in that run
tumor_core <- bind_rows(lapply(RESOLUTIONS, function(r) {
  rr <- as.character(r)
  bind_rows(lapply(names(RECIPES), function(rc) bind_rows(lapply(seed_names, function(sn) {
    d <- results[[sn]][[rr]]$domains[[rc]]
    if (anyNA(d)) return(NULL)
    tc <- names(which.max(tapply(is_tumour, d, mean)))
    data.frame(resolution = r, recipe = rc, seed = sn, sample = levels(samples),
               tumour_frac_of_core = mean(is_tumour[d == tc]),
               core_cells = as.integer(tapply(d == tc, samples, sum)),
               core_frac_of_sample = as.numeric(tapply(d == tc, samples, mean)))
  }))))
}))
write.csv(tumor_core, file.path(OUTPUT_DIR, "tumor_core_by_sample.csv"), row.names = FALSE)

# --- 6. Summary ----------------------------------------------------------------
summary_df <- pairwise %>%
  group_by(resolution) %>%
  summarise(n_clusters = paste(range(c(n_clusters_a, n_clusters_b)), collapse = "-"),
            cluster_ari_mean = mean(cluster_ari),
            across(starts_with("domain_ari_"), list(mean = ~mean(.x), min = ~min(.x)), .names = "{.col}_{.fn}"),
            .groups = "drop") %>%
  left_join(unanimity %>% select(resolution, recipe, frac_cells_unanimous) %>%
              pivot_wider(names_from = recipe, values_from = frac_cells_unanimous,
                          names_prefix = "unanimous_"), by = "resolution")
write.csv(summary_df, file.path(OUTPUT_DIR, "stability_summary.csv"), row.names = FALSE)
cat("\nStability across Harmony seeds by resolution:\n")
print(as.data.frame(summary_df), digits = 3)

cat("\nTumor Core share of each sample (median over seeds [min-max]):\n")
print(as.data.frame(tumor_core %>%
  group_by(recipe, resolution, sample) %>%
  summarise(v = sprintf("%.2f [%.2f-%.2f]", median(core_frac_of_sample),
                        min(core_frac_of_sample), max(core_frac_of_sample)), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = v) %>%
  arrange(recipe, desc(resolution))))

# --- 7. Plots ------------------------------------------------------------------
p_ari <- pairwise %>%
  pivot_longer(c(cluster_ari, starts_with("domain_ari_")), names_to = "level", values_to = "ari") %>%
  mutate(level = factor(c(cluster_ari = "Clusters",
                          setNames(RECIPE_LABELS, paste0("domain_ari_", names(RECIPE_LABELS))))[level],
                        levels = c("Clusters", RECIPE_LABELS))) %>%
  ggplot(aes(x = factor(resolution), y = ari)) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 1.5, color = "grey40") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#D95F02", linewidth = 0.4) +
  facet_wrap(~ level, nrow = 1, labeller = labeller(level = label_wrap_gen(30))) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Louvain resolution", y = "ARI between Harmony seeds (pairs)",
       title = "Agreement across Harmony seeds") +
  theme_publication()

p_unan <- unanimity %>%
  mutate(recipe = factor(RECIPE_LABELS[recipe], levels = RECIPE_LABELS)) %>%
  ggplot(aes(x = factor(resolution), y = frac_cells_unanimous, fill = recipe)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = setNames(c("#7570B3", "#1B9E77"), RECIPE_LABELS), name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Louvain resolution", y = "Cells with the same domain in all seeds",
       title = "Per-cell agreement") +
  theme_publication() +
  theme(legend.position = "bottom", legend.direction = "vertical")

p_tc <- tumor_core %>%
  mutate(recipe = factor(RECIPE_LABELS[recipe], levels = RECIPE_LABELS)) %>%
  ggplot(aes(x = sample, y = core_frac_of_sample, color = sample)) +
  geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.5) +
  scale_color_manual(values = SAMPLE_COLORS, guide = "none") +
  facet_grid(recipe ~ resolution, labeller = labeller(recipe = label_wrap_gen(25),
                                                      resolution = function(x) paste("resolution", x))) +
  labs(x = NULL, y = "Tumor Core share of sample (one point per seed)",
       title = "Tumor Core per sample across Harmony seeds") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "stability_by_resolution.pdf"),
       (p_ari | p_unan) + plot_layout(widths = c(3, 1.2)), width = 15, height = 5.5)
ggsave(file.path(OUTPUT_DIR, "tumor_core_by_sample.pdf"), p_tc, width = 13, height = 6.5)

cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")
