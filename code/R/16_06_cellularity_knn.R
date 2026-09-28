# 16_06_cellularity_knn.R
# ==============================================================================
# Standalone sanity check for reviewer point #1 on 15_04.
# ==============================================================================
#
# THE CONCERN (reviewer #1)
# -------------------------
# 15_04's tumor-proximity covariate is scored by compute_tumor_score() ->
# score_diffusion() (00_utils.R), which sets the Gaussian kernel bandwidth to
#     sigma <- median(knn_res$nn.dist)        (00_utils.R:169)
# computed INSIDE EACH SAMPLE (the dispatcher loops per sample). Median k-NN
# distance is a cell-density proxy (denser tissue -> smaller distances). If
# cellularity tracks tumor burden, and burden tracks condition (DKO ~0.47x WT),
# then the physical length scale of the covariate differs WT vs DKO, and the
# condition x score difference smooth (p_shape) could be reading out "the ruler
# is stretched differently in each arm" rather than biology.
#
# WHAT THIS SCRIPT DOES (does NOT modify 15_04 / 00_utils.R / compute_tumor_score)
# --------------------------------------------------------------------------------
#   Part A. Is the sigma-driving quantity (median k-NN distance = cellularity)
#           actually different by condition? Per-sample table + WT-vs-DKO test.
#           NOTE: median k-NN distance and cell density are the SAME axis
#           (inversely related), not two independent things.
#   Part B. THE DECISIVE TEST. Re-score every cell two ways -- with each sample's
#           OWN sigma (= exactly 15_04's behavior) and with a single GLOBAL sigma
#           -- and correlate the results. If they are ~identical (rank-corr ~1),
#           the per-sample sigma choice is IMMATERIAL for this dataset, which
#           exonerates 15_04 regardless of any sigma difference. This is the
#           reviewer-proof argument: a Wilcoxon at n=3v4 (floor p~0.029, no power)
#           cannot honestly prove the null on its own.
#   Part C. Plain-language verdict -> SUMMARY.txt.
#
# Outputs -> output/R/16_06_cellularity_knn/
# ==============================================================================

# --- Thread limits (set BEFORE BLAS/library load) ----------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS,
           NUMEXPR_NUM_THREADS = N_THREADS)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(N_THREADS)
  RhpcBLASctl::omp_set_num_threads(N_THREADS)
}

library(Seurat)
library(FNN)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)

# Centralized parameters, palette, and scoring utilities. 00_utils.R is sourced
# ONLY so we can run the real compute_tumor_score() as a validation reference in
# Part B -- we do not modify it.
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_06_cellularity_knn")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Parameters: fixed to 15_04's diffusion settings -------------------------
KNN_K           <- 30    # 15_04 KNN_K
HEAT_STEPS      <- 20    # 15_04 DIFFUSION_HEAT_STEPS
ALPHA           <- 2     # 15_04 DIFFUSION_HEAT_ALPHA
TARGET_CELLTYPE <- "Tumor"
TRAJ_CELLTYPE   <- "CD8+ T"   # the cell type whose covariate feeds the GAM

# =============================================================================
# Load data (same object + recoding as 15_04)
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample    <- factor(seu$sample,
                        levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

all_coords <- Embeddings(seu, "spatial")   # native platform units (no um scale exists)
celltypes  <- as.character(seu$Tier1_celltype)
samples    <- seu$sample
is_tumor   <- celltypes == TARGET_CELLTYPE
sample_ids <- levels(samples)

# Sample -> condition lookup
samp_cond <- tapply(as.character(seu$condition), as.character(samples),
                    function(x) x[1])

# =============================================================================
# PART A — cellularity / k-NN distance by condition
# =============================================================================
cat("\n=== Part A: per-sample cellularity / k-NN distance ===\n")

per_sample <- lapply(sample_ids, function(s) {
  idx      <- which(as.character(samples) == s)
  coords_s <- all_coords[idx, , drop = FALSE]
  nn       <- FNN::get.knn(coords_s, k = KNN_K)$nn.dist
  # sigma_knn = median NN distance = EXACTLY the sigma score_diffusion() uses.
  # density_bbox is a crude cells-per-bounding-box-area proxy (tissue is not
  # rectangular, so it is approximate and only for a second density read).
  bbox_area <- diff(range(coords_s[, 1])) * diff(range(coords_s[, 2]))
  data.frame(
    sample       = s,
    condition    = unname(samp_cond[s]),
    n_cells      = length(idx),
    sigma_knn    = median(nn),
    mean_nn      = mean(nn),
    density_bbox = length(idx) / bbox_area,
    tumor_frac   = mean(is_tumor[idx]),
    stringsAsFactors = FALSE)
})
per_sample <- bind_rows(per_sample)
per_sample$condition <- factor(per_sample$condition,
                               levels = c(REFERENCE_CONDITION, TEST_CONDITION))
write.csv(per_sample, file.path(OUTPUT_DIR, "cellularity_per_sample.csv"),
          row.names = FALSE)
cat("Per-sample table:\n"); print(per_sample, row.names = FALSE)

# WT-vs-DKO tabulation (mirror 16_02 summarise_metric style).
# CAVEAT: with 3 vs 4 samples the Wilcoxon floor is 1/choose(7,3) = 0.029 and the
# test is underpowered -- read the fold-ratio and range overlap, not the p.
summarise_metric <- function(df, col) {
  g <- split(df[[col]], df$condition)
  wt <- g[[REFERENCE_CONDITION]]; dko <- g[[TEST_CONDITION]]
  pw <- tryCatch(suppressWarnings(wilcox.test(wt, dko)$p.value), error = function(e) NA_real_)
  data.frame(
    metric        = col,
    WT_mean       = mean(wt),  DKO_mean      = mean(dko),
    WT_range      = sprintf("%.4g-%.4g", min(wt),  max(wt)),
    DKO_range     = sprintf("%.4g-%.4g", min(dko), max(dko)),
    fold_DKO_WT   = mean(dko) / mean(wt),
    ranges_overlap = !(max(wt) < min(dko) || max(dko) < min(wt)),
    wilcox_p      = pw,
    stringsAsFactors = FALSE)
}
by_condition <- bind_rows(lapply(c("sigma_knn", "density_bbox", "tumor_frac", "n_cells"),
                                 function(m) summarise_metric(per_sample, m)))
write.csv(by_condition, file.path(OUTPUT_DIR, "cellularity_by_condition.csv"),
          row.names = FALSE)
cat(sprintf("\nWT-vs-DKO (n=%d vs %d; Wilcoxon floor = %.3f, UNDERPOWERED):\n",
            sum(per_sample$condition == REFERENCE_CONDITION),
            sum(per_sample$condition == TEST_CONDITION),
            1 / choose(nrow(per_sample), sum(per_sample$condition == REFERENCE_CONDITION))))
print(by_condition, row.names = FALSE)

# --- Plots -------------------------------------------------------------------
p_sig <- ggplot(per_sample, aes(x = condition, y = sigma_knn, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = NA) +
  geom_jitter(width = 0.12, size = 3) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS, guide = "none") +
  labs(x = NULL, y = "sigma = median k-NN distance (native units)",
       title = "Per-sample kernel bandwidth (cellularity) by condition") +
  theme_publication()

p_dens <- ggplot(per_sample, aes(x = density_bbox, y = sigma_knn, color = condition)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = "cells / bbox area (approx. density)", y = "sigma = median k-NN distance",
       title = "sigma tracks cell density (sanity: they are the same axis)") +
  theme_publication()

p_tum <- ggplot(per_sample, aes(x = tumor_frac, y = sigma_knn, color = condition)) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = sample), size = 3, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = "tumor fraction", y = "sigma = median k-NN distance",
       title = "sigma vs tumor burden") +
  theme_publication()

# ggrepel is optional; fall back to plain labels if absent.
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  p_sig  <- p_sig  + geom_text(aes(label = sample), size = 3, vjust = -1, show.legend = FALSE)
  p_dens <- p_dens + geom_text(aes(label = sample), size = 3, vjust = -1, show.legend = FALSE)
  p_tum  <- p_tum  + geom_text(aes(label = sample), size = 3, vjust = -1, show.legend = FALSE)
}
ggsave(file.path(OUTPUT_DIR, "sigma_by_condition.pdf"), p_sig, width = 5, height = 4)
ggsave(file.path(OUTPUT_DIR, "sigma_vs_density.pdf"),   p_dens, width = 5.5, height = 4)
ggsave(file.path(OUTPUT_DIR, "sigma_vs_tumorfrac.pdf"), p_tum, width = 5.5, height = 4)

# =============================================================================
# PART B — score-invariance: per-sample sigma vs global sigma
# =============================================================================
cat("\n=== Part B: score-invariance (per-sample vs global sigma) ===\n")

# Faithful copy of score_diffusion() + .diffuse_field(method="heat") from
# 00_utils.R, but with sigma PASSED IN so it is the only thing that varies.
# L_sym heat, backward-Euler, num/den local fraction, no min-max. Matches 15_04's
# diffusion path exactly when sigma = median(nn.dist) of the sample.
diffuse_heat_sym <- function(W, s_raw, alpha, heat_steps) {
  n <- nrow(W)
  d_vec <- Matrix::rowSums(W); d_vec[d_vec == 0] <- 1
  D_inv_sqrt <- Matrix::Diagonal(n, x = 1 / sqrt(d_vec))
  L_sym <- Matrix::Diagonal(n) - D_inv_sqrt %*% W %*% D_inv_sqrt
  A  <- Matrix::Diagonal(n) + alpha * L_sym
  ch <- Matrix::Cholesky(as(A, "symmetricMatrix"), LDL = FALSE, perm = TRUE)
  s <- s_raw
  for (step in seq_len(heat_steps)) s <- as.numeric(Matrix::solve(ch, s, system = "A"))
  s
}

score_diff_sigma <- function(coords_s, is_tumor_s, sigma, k = KNN_K,
                             heat_steps = HEAT_STEPS, alpha = ALPHA) {
  n <- nrow(coords_s)
  if (n < k + 1) return(rep(NA_real_, n))
  knn_res <- FNN::get.knn(coords_s, k = k)
  rows    <- rep(seq_len(n), each = k)
  cols    <- as.vector(t(knn_res$nn.index))
  dists   <- as.vector(t(knn_res$nn.dist))
  weights <- exp(-dists^2 / (2 * sigma^2))
  W <- Matrix::sparseMatrix(i = rows, j = cols, x = weights, dims = c(n, n))
  W <- (W + Matrix::t(W)) / 2
  num <- diffuse_heat_sym(W, as.numeric(as.integer(is_tumor_s)), alpha, heat_steps)
  den <- diffuse_heat_sym(W, rep(1, n), alpha, heat_steps)
  den <- pmax(den, .Machine$double.eps)
  num / den
}

# One global bandwidth: pooled median NN distance over ALL cells.
sigma_global <- median(unlist(lapply(sample_ids, function(s) {
  idx <- which(as.character(samples) == s)
  FNN::get.knn(all_coords[idx, , drop = FALSE], k = KNN_K)$nn.dist
})))
cat(sprintf("Global sigma (pooled median NN dist) = %.4f\n", sigma_global))
cat(sprintf("Per-sample sigma range: %.4f - %.4f (fold %.2fx)\n",
            min(per_sample$sigma_knn), max(per_sample$sigma_knn),
            max(per_sample$sigma_knn) / min(per_sample$sigma_knn)))

score_persample <- rep(NA_real_, ncol(seu))
score_global    <- rep(NA_real_, ncol(seu))
for (s in sample_ids) {
  idx        <- which(as.character(samples) == s)
  coords_s   <- all_coords[idx, , drop = FALSE]
  is_tumor_s <- is_tumor[idx]
  sig_s      <- per_sample$sigma_knn[per_sample$sample == s]
  score_persample[idx] <- score_diff_sigma(coords_s, is_tumor_s, sigma = sig_s)
  score_global[idx]    <- score_diff_sigma(coords_s, is_tumor_s, sigma = sigma_global)
  cat(sprintf("  scored %s (n=%d, sigma_persample=%.3f)\n", s, length(idx), sig_s))
}

# --- Validation: does our per-sample reimplementation reproduce 15_04's real
#     compute_tumor_score()? Spearman must be ~1.000 and max|delta| tiny. --------
cat("\n--- Validation vs real compute_tumor_score() (00_utils.R) ---\n")
score_real <- compute_tumor_score(all_coords, is_tumor, samples, method = "diffusion",
                                  diffusion_method = "heat", k = KNN_K,
                                  heat_steps = HEAT_STEPS, alpha = ALPHA,
                                  local_fraction = TRUE, normalize = FALSE,
                                  verbose = FALSE)
ok_v <- is.finite(score_real) & is.finite(score_persample)
val_spearman <- suppressWarnings(cor(score_real[ok_v], score_persample[ok_v], method = "spearman"))
val_maxabs   <- max(abs(score_real[ok_v] - score_persample[ok_v]))
cat(sprintf("  reimpl vs real: spearman = %.6f, max|delta| = %.3e (want ~1.0 / ~0)\n",
            val_spearman, val_maxabs))

# --- Compare per-sample-sigma vs global-sigma scores -------------------------
ok <- is.finite(score_persample) & is.finite(score_global)
overall <- data.frame(
  set          = "all_cells",
  n            = sum(ok),
  spearman     = suppressWarnings(cor(score_persample[ok], score_global[ok], method = "spearman")),
  pearson      = suppressWarnings(cor(score_persample[ok], score_global[ok], method = "pearson")),
  max_abs_diff = max(abs(score_persample[ok] - score_global[ok])))

cd8 <- celltypes == TRAJ_CELLTYPE
ok8 <- ok & cd8
overall <- rbind(overall, data.frame(
  set          = "CD8pos_T",
  n            = sum(ok8),
  spearman     = suppressWarnings(cor(score_persample[ok8], score_global[ok8], method = "spearman")),
  pearson      = suppressWarnings(cor(score_persample[ok8], score_global[ok8], method = "pearson")),
  max_abs_diff = max(abs(score_persample[ok8] - score_global[ok8]))))
write.csv(overall, file.path(OUTPUT_DIR, "score_invariance_overall.csv"), row.names = FALSE)
cat("Overall / CD8 score invariance (per-sample sigma vs global sigma):\n")
print(overall, row.names = FALSE)

per_samp_inv <- bind_rows(lapply(sample_ids, function(s) {
  idx  <- which(as.character(samples) == s)
  oks  <- idx[is.finite(score_persample[idx]) & is.finite(score_global[idx])]
  ok8s <- oks[celltypes[oks] == TRAJ_CELLTYPE]
  data.frame(
    sample        = s,
    condition     = unname(samp_cond[s]),
    spearman_all  = suppressWarnings(cor(score_persample[oks], score_global[oks], method = "spearman")),
    spearman_cd8  = if (length(ok8s) > 10)
                      suppressWarnings(cor(score_persample[ok8s], score_global[ok8s], method = "spearman"))
                    else NA_real_,
    n_cd8         = length(ok8s),
    stringsAsFactors = FALSE)
}))
write.csv(per_samp_inv, file.path(OUTPUT_DIR, "score_invariance_per_sample.csv"),
          row.names = FALSE)
cat("Per-sample score invariance:\n"); print(per_samp_inv, row.names = FALSE)

# --- Scatter: per-sample-sigma vs global-sigma score (sampled cells) ---------
set.seed(1)
samp_cells <- sample(which(ok), min(40000, sum(ok)))
df_sc <- data.frame(persample = score_persample[samp_cells],
                    global    = score_global[samp_cells],
                    sample    = samples[samp_cells])
p_inv <- ggplot(df_sc, aes(x = persample, y = global, color = sample)) +
  geom_point(size = 0.2, alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  scale_color_manual(values = SAMPLE_COLORS) +
  coord_equal() +
  labs(x = "Tumor score (per-sample sigma = 15_04)", y = "Tumor score (global sigma)",
       title = "Score invariance to sigma choice",
       subtitle = sprintf("all cells rho = %.4f | CD8+ T rho = %.4f",
                          overall$spearman[overall$set == "all_cells"],
                          overall$spearman[overall$set == "CD8pos_T"])) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "score_persample_vs_global.pdf"), p_inv, width = 6, height = 5.5)
ggsave(file.path(OUTPUT_DIR, "score_persample_vs_global.png"), p_inv, width = 6, height = 5.5, dpi = 200)

# =============================================================================
# PART C — verdict summary
# =============================================================================
sig_fold  <- max(per_sample$sigma_knn) / min(per_sample$sigma_knn)
sig_by    <- by_condition[by_condition$metric == "sigma_knn", ]
rho_cd8   <- overall$spearman[overall$set == "CD8pos_T"]
rho_all   <- overall$spearman[overall$set == "all_cells"]
min_rho_ps <- min(per_samp_inv$spearman_cd8, per_samp_inv$spearman_all, na.rm = TRUE)

verdict <- if (rho_cd8 >= 0.99 && rho_all >= 0.99) {
  "IMMATERIAL: the global-sigma and per-sample-sigma scores are ~identical. The per-sample sigma choice does not change the covariate for this dataset, so 15_04's scoring is defensible as-is on this axis. Cite this script."
} else if (rho_cd8 >= 0.95) {
  "MINOR: scores are highly but not perfectly concordant. Report the correlation; consider the global-sigma fix if a reviewer presses."
} else {
  "MATERIAL: the sigma choice measurably changes the score. The per-sample bandwidth is confounding the covariate; switch to a global sigma (or 15_06) and re-run."
}

con <- file(file.path(OUTPUT_DIR, "SUMMARY.txt"), open = "wt")
writeLines(c(
  strrep("=", 78),
  "16_06 cellularity / k-NN bandwidth sanity check -- SUMMARY",
  strrep("=", 78),
  "",
  "REVIEWER #1 CONCERN: score_diffusion() sets sigma = median(k-NN distance) per",
  "sample; if cellularity tracks condition, the covariate's length scale is",
  "condition-confounded.",
  "",
  "-- Part A: cellularity (sigma = median k-NN distance) by condition --",
  sprintf("  WT mean sigma  = %.4f   (range %s)", sig_by$WT_mean, sig_by$WT_range),
  sprintf("  DKO mean sigma = %.4f   (range %s)", sig_by$DKO_mean, sig_by$DKO_range),
  sprintf("  DKO/WT fold    = %.3f", sig_by$fold_DKO_WT),
  sprintf("  ranges overlap = %s", sig_by$ranges_overlap),
  sprintf("  Wilcoxon p     = %.4f  (n=3v4 floor 0.029; UNDERPOWERED, read the fold/overlap)",
          sig_by$wilcox_p),
  sprintf("  per-sample sigma spread = %.2fx (min %.4f, max %.4f)",
          sig_fold, min(per_sample$sigma_knn), max(per_sample$sigma_knn)),
  "",
  "-- Validation: reimplementation vs real compute_tumor_score() --",
  sprintf("  spearman = %.6f, max|delta| = %.3e  (must be ~1.0 / ~0)", val_spearman, val_maxabs),
  "",
  "-- Part B: score invariance to sigma choice (per-sample vs global) --",
  sprintf("  all cells: spearman = %.5f, pearson = %.5f, max|diff| = %.4f",
          rho_all, overall$pearson[overall$set == "all_cells"],
          overall$max_abs_diff[overall$set == "all_cells"]),
  sprintf("  CD8+ T   : spearman = %.5f, pearson = %.5f, max|diff| = %.4f",
          rho_cd8, overall$pearson[overall$set == "CD8pos_T"],
          overall$max_abs_diff[overall$set == "CD8pos_T"]),
  sprintf("  worst per-sample spearman = %.5f", min_rho_ps),
  "",
  "-- VERDICT --",
  strwrap(verdict, width = 76, prefix = "  "),
  "",
  strrep("=", 78)), con)
close(con)

cat("\n", readLines(file.path(OUTPUT_DIR, "SUMMARY.txt")), sep = "\n")
cat(sprintf("\nAll outputs -> %s\n", OUTPUT_DIR))
