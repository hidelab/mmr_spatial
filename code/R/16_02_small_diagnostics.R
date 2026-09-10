# 16_02_small_diagnostics.R
# ============================================================================
# Three small sanity checks on the REAL cohort that decide whether `p_shape`
# from 15_04 is interpretable as written. All three ask the same question from
# different angles: is the tumor-proximity axis (and the graph it is built on)
# CONFOUNDED with condition? If any metric separates WT from DKO, then a
# WT-vs-DKO difference in the trajectory could be an artifact of the axis
# construction rather than a biological remodeling of the gradient.
#
#   1. Graph bandwidth (sigma) by condition
#      sigma = median k-NN distance per sample = the Gaussian kernel width used
#      by score_diffusion(). It is a direct proxy for local cell density: denser
#      tissue -> smaller NN distances. If WT and DKO pack cells at different
#      densities, the diffusion length scale differs systematically by condition.
#      -> confirms whether the k-NN graph "follows actual cell density".
#
#   2. RAW (pre-normalization) tumor field by condition
#      compute_tumor_score(..., normalize = FALSE). The per-sample [0,1] min-max
#      in 15_04 hides absolute differences in field magnitude/spread; here we
#      look at the raw diffused field before that rescaling.
#
#   3. Tumor fraction by condition
#      tapply(is_tumor, samples, mean). If DKO simply has more/less tumor, the
#      source term of the diffusion differs and the axis is not comparable.
#
# Uses the SAME scoring parameters as 15_04 (diffusion/heat, k=30, heat_steps=20)
# so the numbers describe the axis actually used in the trajectory analysis.
#
# Outputs (-> output/R/16_small_diagnostics/):
#   * diagnostics_by_sample.csv     — one row per sample, all three metrics
#   * diagnostics_by_condition.csv  — WT vs DKO summary + underpowered test
#   * diagnostics_by_condition.png  — per-metric dotplot, colored by condition
# ============================================================================

# --- Thread limits (set BEFORE BLAS/library load) ---------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS,
           NUMEXPR_NUM_THREADS = N_THREADS)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(N_THREADS)
  RhpcBLASctl::omp_set_num_threads(N_THREADS)
}

library(Seurat)
library(ggplot2)
library(patchwork)
library(FNN)
library(Matrix)
library(dplyr)

# Centralized parameters, palette, and the diffusion scorers.
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_small_diagnostics")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters — MUST match the axis used in 15_04 for the checks to be relevant
# ============================================================================
KNN_K                <- 30        # k for the spatial k-NN graph (= 15_04 KNN_K)
DIFFUSION_VARIANT    <- "heat"    # score_diffusion variant used in 15_04
DIFFUSION_HEAT_STEPS <- 20        # backward-Euler steps used in 15_04
TARGET_CELLTYPE      <- "Tumor"

# ============================================================================
# Load data (same recoding as 15_04)
# ============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

all_coords <- Embeddings(seu, "spatial")
samples    <- seu$sample
is_tumor   <- seu$Tier1_celltype == TARGET_CELLTYPE

sample_ids <- levels(seu$sample)
# Condition label per sample (from the first cell of each sample).
samp_cond  <- sapply(sample_ids, function(s)
  as.character(seu$condition[which(samples == s)[1]]))

# ============================================================================
# Diagnostic 1 — graph bandwidth (sigma) by sample
#   sigma = median k-NN distance = Gaussian kernel width in score_diffusion().
#   Smaller = denser tissue. Does the k-NN scale follow actual cell density,
#   and does that density differ by condition?
# ============================================================================
cat("\n[1] Graph bandwidth (median k-NN distance, k =", KNN_K, ") per sample...\n")
sigma_by_sample <- sapply(sample_ids, function(s) {
  i <- which(samples == s)
  median(FNN::get.knn(all_coords[i, , drop = FALSE], k = KNN_K)$nn.dist)
})

# ============================================================================
# Diagnostic 2 — RAW (pre-normalization) tumor field by sample
#   Same diffusion as 15_04 but normalize = FALSE, so absolute field magnitude
#   and spread are visible (the per-sample min-max rescale would erase them).
# ============================================================================
cat("[2] RAW diffusion field (normalize = FALSE) per sample...\n")
raw_field <- compute_tumor_score(
  all_coords, is_tumor, samples, method = "diffusion",
  diffusion_method = DIFFUSION_VARIANT, k = KNN_K,
  heat_steps = DIFFUSION_HEAT_STEPS,
  normalize = FALSE, verbose = FALSE)

raw_mean   <- tapply(raw_field, samples, mean,   na.rm = TRUE)[sample_ids]
raw_median <- tapply(raw_field, samples, median, na.rm = TRUE)[sample_ids]
raw_max    <- tapply(raw_field, samples, max,    na.rm = TRUE)[sample_ids]
# Spread of the raw field (dynamic range available before rescaling).
raw_range  <- tapply(raw_field, samples,
                     function(v) diff(range(v, na.rm = TRUE)))[sample_ids]

# ============================================================================
# Diagnostic 3 — tumor fraction by sample
# ============================================================================
cat("[3] Tumor fraction per sample...\n")
tumor_frac <- tapply(as.numeric(is_tumor), samples, mean)[sample_ids]
n_cells    <- as.numeric(table(samples)[sample_ids])

# ============================================================================
# Assemble per-sample table
# ============================================================================
by_sample <- data.frame(
  sample        = sample_ids,
  condition     = samp_cond,
  n_cells       = n_cells,
  sigma_knn     = as.numeric(sigma_by_sample),
  raw_field_mean   = as.numeric(raw_mean),
  raw_field_median = as.numeric(raw_median),
  raw_field_max    = as.numeric(raw_max),
  raw_field_range  = as.numeric(raw_range),
  tumor_frac    = as.numeric(tumor_frac),
  row.names = NULL,
  stringsAsFactors = FALSE
)
by_sample$condition <- factor(by_sample$condition,
                              levels = c(REFERENCE_CONDITION, TEST_CONDITION))
by_sample <- by_sample[order(by_sample$condition, by_sample$sample), ]

cat("\n================ Per-sample diagnostics ================\n")
print(by_sample, row.names = FALSE, digits = 4)

write.csv(by_sample,
          file.path(OUTPUT_DIR, "diagnostics_by_sample.csv"),
          row.names = FALSE)

# ============================================================================
# Condition-level summary + (underpowered) rank-sum test
#   n = 3 WT vs 4 DKO, so p-values are illustrative only. The eyeball test
#   (do the two groups overlap?) matters more than the p-value here.
# ============================================================================
metrics <- c("sigma_knn", "raw_field_mean", "raw_field_median",
             "raw_field_max", "raw_field_range", "tumor_frac")

summarise_metric <- function(m) {
  wt  <- by_sample[[m]][by_sample$condition == REFERENCE_CONDITION]
  dko <- by_sample[[m]][by_sample$condition == TEST_CONDITION]
  p <- tryCatch(
    suppressWarnings(wilcox.test(wt, dko)$p.value),
    error = function(e) NA_real_)
  data.frame(
    metric      = m,
    WT_mean     = mean(wt),
    DKO_mean    = mean(dko),
    WT_range    = paste0(signif(min(wt), 3),  "–", signif(max(wt), 3)),
    DKO_range   = paste0(signif(min(dko), 3), "–", signif(max(dko), 3)),
    ratio_DKO_WT = mean(dko) / mean(wt),
    wilcox_p    = p,
    stringsAsFactors = FALSE)
}
by_condition <- do.call(rbind, lapply(metrics, summarise_metric))

cat("\n================ WT vs DKO summary (n = 3 vs 4; wilcox_p is underpowered) ================\n")
print(by_condition, row.names = FALSE, digits = 4)

write.csv(by_condition,
          file.path(OUTPUT_DIR, "diagnostics_by_condition.csv"),
          row.names = FALSE)

# ============================================================================
# Visual: per-metric dot per sample, colored by condition
# ============================================================================
plot_metric <- function(m, ylab) {
  ggplot(by_sample, aes(x = condition, y = .data[[m]],
                        color = condition, label = sample)) +
    geom_point(size = 3, position = position_jitter(width = 0.08, height = 0)) +
    ggrepel::geom_text_repel(size = 3, show.legend = FALSE, seed = 1) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(x = NULL, y = ylab, title = m) +
    theme_publication() +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold"))
}

labs_map <- c(
  sigma_knn        = "median k-NN dist (sigma)",
  raw_field_mean   = "raw field mean",
  raw_field_median = "raw field median",
  raw_field_max    = "raw field max",
  raw_field_range  = "raw field range",
  tumor_frac       = "tumor fraction")

have_repel <- requireNamespace("ggrepel", quietly = TRUE)
if (!have_repel) {
  plot_metric <- function(m, ylab) {
    ggplot(by_sample, aes(x = condition, y = .data[[m]], color = condition)) +
      geom_point(size = 3, position = position_jitter(width = 0.08, height = 0)) +
      scale_color_manual(values = CONDITION_COLORS) +
      labs(x = NULL, y = ylab, title = m) +
      theme_publication() +
      theme(legend.position = "none",
            plot.title = element_text(size = 10, face = "bold"))
  }
}

panels <- lapply(metrics, function(m) plot_metric(m, labs_map[[m]]))
p_all  <- wrap_plots(panels, ncol = 3) +
  plot_annotation(
    title = "15_04 tumor-axis confounding checks: does any metric separate WT vs DKO?",
    theme = theme(plot.title = element_text(size = 12, face = "bold")))

ggsave(file.path(OUTPUT_DIR, "diagnostics_by_condition.png"),
       p_all, width = 10, height = 6.5, dpi = 150)

cat("\nDone. Tables + plot saved to:", OUTPUT_DIR, "\n")