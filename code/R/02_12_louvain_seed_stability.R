# 02_12_louvain_seed_stability.R
# Is the clustering step (FindClusters) itself a source of the domain
# instability, independent of Harmony?
#
# 02_10/02_11 varied the Harmony seed with the Louvain seed fixed (random.seed =
# 0); clusters agreed at ARI ~0.45-0.57 and domains at ~0.4-0.5 even when the
# embeddings were 95% identical. Here the stored Harmony embedding and the stored
# SNN graph (the ones behind the current domains) are held fixed and only the
# FindClusters seed changes. Configurations, each with seeds LOUVAIN_SEEDS:
#   louvain_nstart10   Louvain, n.start = 10 (Seurat default, as in 02_01)
#   louvain_nstart100  Louvain, n.start = 100
#   leiden             Leiden (algorithm = 4)
# at every resolution in RESOLUTIONS. Each clustering is merged with both domain
# recipes of 02_10 and compared across seeds (pairwise ARI, cells given the same
# domain by all seeds, Tumor Core share per sample). The Harmony-seed results of
# 02_10 are joined for comparison.
#
# Read-only: nothing is saved except CSV tables and plots, all in a folder used
# only by this script (no RDS output). Harmony is not re-run.

library(Seurat)
library(Matrix)
library(parallel)
library(ggplot2)
library(dplyr)
library(tidyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "louvain_seed_stability")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
LOUVAIN_SEEDS    <- 1:5                 # Leiden needs seeds > 0
RESOLUTIONS      <- c(0.8, 0.5)
CONFIGS <- list(louvain_nstart10  = list(algorithm = 1, n.start = 10),
                louvain_nstart100 = list(algorithm = 1, n.start = 100),
                leiden            = list(algorithm = 4, n.start = 10))
MIN_CLUSTER_FRAC <- 0.005               # minsize: clusters below this fraction of cells are dissolved
RECIPES  <- c(current_k4 = "unweighted", minsize_k5 = "minsize")
RECIPE_K <- c(current_k4 = MERGE_K, minsize_k5 = 5)
N_CORES  <- min(5, detectCores())
HARMONY_SEED_SUMMARY <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08",
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
merge_domains <- function(method, clusters, k) {
  labels <- as.character(clusters)
  if (method == "minsize") {
    labels <- reassign_small(labels, snn_dgc, harmony_emb, MIN_CLUSTER_FRAC * length(labels))
  }
  if (length(unique(labels)) < k) return(rep(NA_integer_, length(labels)))
  unname(cutree(centroid_tree(harmony_emb, labels), k = k)[labels])
}

# Relabel `x` so that each label takes the `ref` label it shares most cells with
# under a one-to-one (Hungarian) matching.
match_to_ref <- function(x, ref) {
  tab <- table(x, ref)
  asg <- clue::solve_LSAP(unclass(tab), maximum = TRUE)
  unname(setNames(colnames(tab)[as.integer(asg)], rownames(tab))[as.character(x)])
}

# --- 1. Load stored embedding and graph (read-only) ----------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")
snn_name        <- grep("snn", Graphs(seu), value = TRUE)[1]
snn_graph       <- seu[[snn_name]]
snn_dgc         <- as(snn_graph, "dgCMatrix")
harmony_emb     <- Embeddings(seu, "harmony")
stored_clusters <- as.character(seu@meta.data[[DOMAIN_SOURCE_COL]])
samp_fac        <- factor(seu$sample, levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))
is_tumour       <- seu$Tier1_celltype == "Tumor"
rm(seu); invisible(gc())
cat(sprintf("Fixed inputs: stored Harmony embedding, SNN graph '%s'\n", snn_name))

# --- 2. FindClusters with each config x resolution x seed ----------------------
grid <- expand.grid(config = names(CONFIGS), resolution = RESOLUTIONS, seed = LOUVAIN_SEEDS,
                    stringsAsFactors = FALSE)
cat(sprintf("Running %d clusterings on %d cores...\n", nrow(grid), N_CORES))
runs <- mclapply(seq_len(nrow(grid)), function(i) {
  cfg <- CONFIGS[[grid$config[i]]]
  t   <- system.time(
    cl <- FindClusters(snn_graph, resolution = grid$resolution[i], algorithm = cfg$algorithm,
                       n.start = cfg$n.start, random.seed = grid$seed[i], verbose = FALSE)[[1]]
  )[["elapsed"]]
  cls <- as.character(cl)
  list(clusters = cls, seconds = t,
       domains = lapply(setNames(names(RECIPES), names(RECIPES)), function(rc)
         merge_domains(RECIPES[[rc]], cls, RECIPE_K[[rc]])))
}, mc.cores = N_CORES, mc.preschedule = FALSE)
stopifnot(!any(sapply(runs, inherits, "try-error")))
grid$n_clusters <- sapply(runs, function(r) length(unique(r$clusters)))
grid$seconds    <- sapply(runs, `[[`, "seconds")
write.csv(grid, file.path(OUTPUT_DIR, "run_info.csv"), row.names = FALSE)

# Reference check: Louvain seed 0 on the stored graph reproduces the stored clusters
ref0 <- as.character(FindClusters(snn_graph, resolution = BANKSY_RESOLUTION, random.seed = 0,
                                  verbose = FALSE)[[1]])
cat(sprintf("Louvain seed 0, resolution %.1f vs stored clusters: ARI = %.3f\n",
            BANKSY_RESOLUTION, ari(ref0, stored_clusters)))

# --- 3. Agreement across seeds -------------------------------------------------
settings <- unique(grid[, c("config", "resolution")])
pairwise <- bind_rows(lapply(seq_len(nrow(settings)), function(j) {
  idx <- which(grid$config == settings$config[j] & grid$resolution == settings$resolution[j])
  bind_rows(lapply(combn(idx, 2, simplify = FALSE), function(p) {
    a <- runs[[p[1]]]; b <- runs[[p[2]]]
    out <- data.frame(config = settings$config[j], resolution = settings$resolution[j],
                      seed_a = grid$seed[p[1]], seed_b = grid$seed[p[2]],
                      cluster_ari = ari(a$clusters, b$clusters))
    for (rc in names(RECIPES)) {
      out[[paste0("domain_ari_", rc)]] <- if (anyNA(a$domains[[rc]]) || anyNA(b$domains[[rc]])) NA
                                          else ari(a$domains[[rc]], b$domains[[rc]])
    }
    out
  }))
}))
write.csv(pairwise, file.path(OUTPUT_DIR, "pairwise_seed_agreement.csv"), row.names = FALSE)

unanimity <- bind_rows(lapply(seq_len(nrow(settings)), function(j) {
  idx <- which(grid$config == settings$config[j] & grid$resolution == settings$resolution[j])
  bind_rows(lapply(names(RECIPES), function(rc) {
    doms <- lapply(idx, function(i) runs[[i]]$domains[[rc]])
    if (any(sapply(doms, anyNA))) {
      return(data.frame(config = settings$config[j], resolution = settings$resolution[j],
                        recipe = rc, frac_cells_unanimous = NA))
    }
    matched <- sapply(doms, match_to_ref, ref = doms[[1]])
    data.frame(config = settings$config[j], resolution = settings$resolution[j], recipe = rc,
               frac_cells_unanimous = mean(rowSums(matched == matched[, 1]) == ncol(matched)))
  }))
}))
write.csv(unanimity, file.path(OUTPUT_DIR, "per_cell_unanimity.csv"), row.names = FALSE)

# Tumor Core = the domain with the highest tumour-cell fraction in that run
tumor_core <- bind_rows(lapply(seq_len(nrow(grid)), function(i) bind_rows(lapply(names(RECIPES), function(rc) {
  d <- runs[[i]]$domains[[rc]]
  if (anyNA(d)) return(NULL)
  tc <- names(which.max(tapply(is_tumour, d, mean)))
  data.frame(config = grid$config[i], resolution = grid$resolution[i], seed = grid$seed[i], recipe = rc,
             sample = levels(samp_fac), core_frac_of_sample = as.numeric(tapply(d == tc, samp_fac, mean)))
}))))
write.csv(tumor_core, file.path(OUTPUT_DIR, "tumor_core_by_sample.csv"), row.names = FALSE)

# --- 4. Summary, with the Harmony-seed results of 02_10 for comparison ---------
summary_df <- pairwise %>%
  group_by(config, resolution) %>%
  summarise(cluster_ari_mean = mean(cluster_ari),
            across(starts_with("domain_ari_"), list(mean = ~mean(.x), min = ~min(.x)), .names = "{.col}_{.fn}"),
            .groups = "drop") %>%
  left_join(grid %>% group_by(config, resolution) %>%
              summarise(n_clusters = paste(range(n_clusters), collapse = "-"),
                        seconds_median = median(seconds), .groups = "drop"),
            by = c("config", "resolution")) %>%
  left_join(unanimity %>% pivot_wider(names_from = recipe, values_from = frac_cells_unanimous,
                                      names_prefix = "unanimous_"), by = c("config", "resolution")) %>%
  mutate(varied = "FindClusters seed (fixed embedding + graph)", .before = 1)

if (file.exists(HARMONY_SEED_SUMMARY)) {
  harmony_ref <- read.csv(HARMONY_SEED_SUMMARY) %>%
    filter(resolution %in% RESOLUTIONS) %>%
    transmute(varied = "Harmony seed (02_10, Louvain seed fixed)", config = "louvain_nstart10", resolution,
              cluster_ari_mean, domain_ari_current_k4_mean, domain_ari_current_k4_min,
              domain_ari_minsize_k5_mean, domain_ari_minsize_k5_min, n_clusters,
              unanimous_current_k4, unanimous_minsize_k5)
  summary_df <- bind_rows(summary_df, harmony_ref)
}
summary_df <- summary_df %>% arrange(desc(resolution), config, varied)
write.csv(summary_df, file.path(OUTPUT_DIR, "stability_summary.csv"), row.names = FALSE)
cat("\nStability across seeds:\n")
print(as.data.frame(summary_df %>% select(-any_of("seconds_median"))), digits = 3)

cat("\nTumor Core share of each sample (median over seeds [min-max]):\n")
print(as.data.frame(tumor_core %>%
  group_by(recipe, config, resolution, sample) %>%
  summarise(v = sprintf("%.2f [%.2f-%.2f]", median(core_frac_of_sample),
                        min(core_frac_of_sample), max(core_frac_of_sample)), .groups = "drop") %>%
  pivot_wider(names_from = sample, values_from = v) %>%
  arrange(recipe, desc(resolution), config)))

# --- 5. Plot -------------------------------------------------------------------
p_ari <- pairwise %>%
  pivot_longer(c(cluster_ari, starts_with("domain_ari_")), names_to = "level", values_to = "ari") %>%
  mutate(level = factor(c(cluster_ari = "Clusters", domain_ari_current_k4 = sprintf("Domains, unweighted k = %d", MERGE_K),
                          domain_ari_minsize_k5 = "Domains, minsize k = 5")[level],
                        levels = c("Clusters", sprintf("Domains, unweighted k = %d", MERGE_K), "Domains, minsize k = 5")),
         config = factor(config, levels = names(CONFIGS))) %>%
  ggplot(aes(x = config, y = ari)) +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 1), size = 1.5, color = "grey40") +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "#D95F02", linewidth = 0.4) +
  facet_grid(paste("resolution", resolution) ~ level) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "ARI between FindClusters seeds (pairs)",
       title = "Stored Harmony embedding and SNN graph, only the clustering seed varied") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "louvain_seed_stability.pdf"), p_ari, width = 11, height = 7)

cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")
