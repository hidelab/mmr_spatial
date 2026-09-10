# 16_03_small_permutation_diagnostic.R
# ============================================================================
# Sample-level permutation calibration for the 15_04 tumor-trajectory shape test
# (reviewer point #4).
#
# ---------------------------------------------------------------------------
# WHY THIS EXISTS (the rationale, in full)
# ---------------------------------------------------------------------------
# 15_04 tests whether the WT->DKO expression trajectory along the tumor-proximity
# axis differs in SHAPE (the `p_shape` from the mgcv ordered-factor difference
# smooth, fit on sample x bin pseudobulk with a per-sample factor smooth `fs`).
# Two facts about that test are BOTH true and must be held together:
#
#   (A) The MODEL is sound. The `fs` term s(score, sample, bs="fs") gives each
#       mouse its own penalized curve, so between-mouse trajectory-shape
#       heterogeneity goes into the error term (not attributed to condition).
#       A variance-component check (Diagnostic 1 below) confirms the per-sample
#       smooths are NOT collapsed (edf well above 0), so the difference smooth is
#       genuinely tested against between-mouse variation.
#
#   (B) The mgcv ASYMPTOTIC p-values are NOT trustworthy in absolute magnitude at
#       n = 3 WT + 4 DKO. They are asymptotic and conditional on the estimated
#       smoothing parameters (including a between-sample variance component that is
#       barely identified from 7 mice). Values like Ptprc's p ~ 1e-16 claim a
#       resolution the design cannot support.
#
# The honest resolution of the design is a SAMPLE-LEVEL PERMUTATION. There are
# choose(7,3) = 35 ways to assign the condition label across the 7 mice. A
# permutation test therefore has a HARD FLOOR of 1/35 ~= 0.029 on any achievable
# p-value, and BH across ~340 genes makes per-gene FDR < 0.1 unreachable
# (0.029 x 340 >> 0.1). So per-gene FDR bars are not a defensible headline; the
# design supports (i) a GLOBAL statement, and (ii) a per-gene RANKING.
#
# This script provides exactly those two, using the SAME pseudobulk_fs model and
# score as 15_04:
#
#   Diagnostic 1 — variance component: report the fitted `fs` edf per reported
#     gene, to show the per-sample smooths are not degenerate (model is sound).
#
#   Diagnostic 2a — per-gene permutation: for each reported gene, refit the fs
#     model under all 35 label assignments, and rank the OBSERVED WT/DKO split's
#     shape statistic among them. permutation p = (# assignments with statistic >=
#     observed) / 35. This calibrates each gene's asymptotic p against the honest
#     sample-level null (expect big genes to sit at rank 1, perm p = 0.029).
#
#   Diagnostic 2b — GLOBAL test: for each of the 35 assignments, run the whole
#     shape test over all tested genes, BH-adjust within that assignment, and count
#     genes at asymptotic FDR < SIG_FDR. If the OBSERVED assignment yields more such
#     genes than every permuted assignment, that is a global "the DKO trajectory
#     differs" statement at p = 1/35 = 0.029 — the strongest the n supports.
#
# NOTE: permutation does not "fix" the model (the model is fine); it calibrates the
# p-values and quantifies the design's true resolution. Present per-gene results as
# a ranking (cf. tradeSeq), lead inference with the global test.
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
library(MASS)
library(Matrix)
library(mgcv)
library(dplyr)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_permutation_diagnostic")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters — MUST match 15_04 so this calibrates the SAME model / score
# ============================================================================
CELLTYPES        <- c("CD8+ T")     # cell type(s) to calibrate
TARGET_CELLTYPE  <- "Tumor"

# Score (local tumor fraction, no per-sample min-max — matches 15_04 #6 default).
# Keep DIFFUSION_HEAT_STEPS / DIFFUSION_HEAT_ALPHA / GAM_K in sync with 15_04 so this
# calibrates the SAME model.
KNN_K            <- 30
DIFFUSION_VARIANT<- "heat"
DIFFUSION_HEAT_STEPS <- 20
DIFFUSION_HEAT_ALPHA <- 2
LOCAL_FRACTION   <- TRUE

# Pseudobulk / GAM (match 15_04)
N_BINS            <- 10
MIN_CELLS_PER_BIN <- 5
GAM_K             <- 5
FS_K              <- 3    # per-sample fs basis dim; small -> order-stable (see notes.md)
SIG_FDR           <- 0.1            # asymptotic-FDR cutoff used for the GLOBAL count
# Match 15_04's common-support + score-axis choices. IMPORTANT: common support
# depends on the condition labels, so it is recomputed INSIDE each permuted
# labeling (not fixed to the observed split) — otherwise the permutation test is
# not a deterministic function of the labels and is invalid.
COMMON_SUPPORT               <- TRUE
MIN_SAMPLES_PER_COND_PER_BIN <- 1
PB_SCORE_AXIS                <- "score"   # "score" (real mean local fraction) | "decile"

# Per-gene permutation report: how many top genes (ranked by observed asymptotic
# p_shape) to profile in full (Diagnostics 1 + 2a). The GLOBAL test (2b) always
# uses ALL tested genes.
PERGENE_REPORT_N  <- 20
# Cap the number of genes entered into the (expensive) global permutation loop.
# NA = all tested genes. Lower it for a quick smoke test.
GLOBAL_MAX_GENES  <- NA

# ============================================================================
# Load + score (identical recoding to 15_04)
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

cat("Computing local-tumor-fraction score...\n")
tumor_score <- compute_tumor_score(
  all_coords, is_tumor, samples, method = "diffusion",
  diffusion_method = DIFFUSION_VARIANT, k = KNN_K, heat_steps = DIFFUSION_HEAT_STEPS,
  alpha = DIFFUSION_HEAT_ALPHA,
  local_fraction = LOCAL_FRACTION, normalize = !LOCAL_FRACTION, verbose = FALSE)

# ============================================================================
# Core fit: the SAME pseudobulk_fs model as 15_04, returning the difference-smooth
# statistic + p-value + fs edf. `wt_set` is the set of samples labelled WT.
# ============================================================================
fit_shape <- function(y_by_group, g_score, g_sample, g_bin, g_n, wt_set) {
  cond <- ifelse(g_sample %in% wt_set, REFERENCE_CONDITION, TEST_CONDITION)
  keep <- rep(TRUE, length(cond))
  # Common support recomputed for THIS labeling: keep only bins with
  # >= MIN_SAMPLES_PER_COND_PER_BIN distinct samples in each condition.
  if (COMMON_SUPPORT) {
    bins <- unique(g_bin)
    ok <- vapply(bins, function(b) {
      m <- g_bin == b
      length(unique(g_sample[m & cond == REFERENCE_CONDITION])) >= MIN_SAMPLES_PER_COND_PER_BIN &&
      length(unique(g_sample[m & cond == TEST_CONDITION]))      >= MIN_SAMPLES_PER_COND_PER_BIN
    }, logical(1))
    keep <- g_bin %in% bins[ok]
  }
  if (sum(keep) < 4 || length(unique(cond[keep])) < 2)
    return(c(F = NA, p = NA, edf_fs = NA))
  d <- data.frame(expr = y_by_group[keep], score = g_score[keep],
                  sample = factor(g_sample[keep]),
                  cond_ord = ordered(cond[keep], levels = c(REFERENCE_CONDITION, TEST_CONDITION)),
                  n = g_n[keep])
  f <- tryCatch(suppressWarnings(
         gam(expr ~ s(score, k = GAM_K) + s(score, by = cond_ord, k = GAM_K) +
                    cond_ord + s(score, sample, bs = "fs", k = FS_K),
             data = d, weights = d$n, method = "REML")),
       error = function(e) NULL)
  if (is.null(f)) return(c(F = NA, p = NA, edf_fs = NA))
  sm <- summary(f); rn <- rownames(sm$s.table)
  dr <- grep(":cond_ord", rn); fr <- grep("sample", rn)
  c(F = if (length(dr)) sm$s.table[dr[1], 3] else NA,
    p = if (length(dr)) sm$s.table[dr[1], 4] else NA,
    edf_fs = if (length(fr)) sm$s.table[fr[1], "edf"] else NA)
}

# ============================================================================
# Per cell type
# ============================================================================
for (ct in CELLTYPES) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
  ct_dir  <- file.path(OUTPUT_DIR, ct_safe)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n=== %s ===\n", ct))

  ct_idx <- which(seu$Tier1_celltype == ct & !is.na(tumor_score))
  if (length(ct_idx) < 50) { cat("  too few cells; skip\n"); next }

  ct_score <- tumor_score[ct_idx]
  ct_sample<- droplevels(seu$sample[ct_idx])
  expr     <- GetAssayData(seu, assay = "RNA", layer = "data")[, ct_idx]
  counts   <- GetAssayData(seu, assay = "RNA", layer = "counts")[, ct_idx]

  # --- test genes: SAME baseline-expression filter as 15_04 / DE scripts ------
  det_frac  <- Matrix::rowSums(counts > 0) / length(ct_idx)
  samp_sum  <- vapply(levels(ct_sample), function(s)
                 Matrix::rowSums(counts[, ct_sample == s, drop = FALSE]),
                 numeric(nrow(counts)))
  n_samp_det <- rowSums(samp_sum >= MIN_COUNTS_PER_SAMPLE)
  test_genes <- rownames(counts)[det_frac >= MIN_GENE_EXPR_FRAC &
                                 n_samp_det >= MIN_SAMPLES_EXPRESSED]
  test_genes <- intersect(test_genes, rownames(expr))
  cat(sprintf("  %d tested genes (%d cells)\n", length(test_genes), length(ct_idx)))

  # --- build pseudobulk group structure ONCE (sample x bin) -------------------
  # Only the condition labels change under permutation, so groups + per-gene group
  # means are computed once and reused across all 35 assignments.
  br <- quantile(ct_score, seq(0, 1, length.out = N_BINS + 1), na.rm = TRUE)
  br[1] <- br[1] - 1e-6
  bin <- as.integer(cut(ct_score, breaks = br, labels = seq_len(N_BINS)))
  grp <- paste(ct_sample, bin, sep = "__")
  gi  <- split(seq_along(grp), grp)
  gi  <- gi[vapply(gi, length, 0L) >= MIN_CELLS_PER_BIN]
  g_sample <- vapply(strsplit(names(gi), "__"), `[`, "", 1)
  g_bin    <- as.integer(vapply(strsplit(names(gi), "__"), `[`, "", 2))
  g_n      <- vapply(gi, length, 0L)
  # x-axis: real mean local tumor fraction per group (matches 15_04 PB_SCORE_AXIS),
  # or evenly-spaced decile index. Fixed across permutations (only labels change).
  g_score  <- if (PB_SCORE_AXIS == "score")
                vapply(gi, function(ix) mean(ct_score[ix]), 0.0)
              else (g_bin - 0.5) / N_BINS
  cat(sprintf("  %d pseudobulk points across %d samples\n",
              length(gi), length(unique(g_sample))))

  # Per-gene group-mean expression matrix (genes x groups) — computed once.
  gene_group_mean <- function(g) vapply(gi, function(ix) mean(expr[g, ix]), 0.0)

  # --- enumerate all label assignments (exhaustive; not random) ---------------
  samps  <- levels(ct_sample)
  n_wt   <- length(intersect(WT_SAMPLES, samps))
  combos <- combn(samps, n_wt)                       # each column = the WT set
  n_perm <- ncol(combos)
  obs_j  <- which(apply(combos, 2, function(x) setequal(x, intersect(WT_SAMPLES, samps))))
  cat(sprintf("  %d label assignments (choose(%d,%d)); permutation floor p = %.3f\n",
              n_perm, length(samps), n_wt, 1 / n_perm))

  # ==========================================================================
  # OBSERVED asymptotic pass (true labels) over ALL test genes
  # ==========================================================================
  obs <- t(vapply(test_genes, function(g)
            fit_shape(gene_group_mean(g), g_score, g_sample, g_bin, g_n, WT_SAMPLES),
            numeric(3)))
  obs <- as.data.frame(obs); obs$gene <- test_genes
  obs$fdr_shape <- p.adjust(obs$p, method = "BH")
  obs <- obs[order(obs$p), ]
  obs_count <- sum(obs$fdr_shape < SIG_FDR, na.rm = TRUE)
  cat(sprintf("  observed: %d genes at asymptotic FDR < %.2f\n", obs_count, SIG_FDR))

  # ==========================================================================
  # Diagnostic 1 + 2a — per-gene permutation for the top PERGENE_REPORT_N genes
  # ==========================================================================
  report_genes <- head(obs$gene, PERGENE_REPORT_N)
  pergene <- lapply(report_genes, function(g) {
    y <- gene_group_mean(g)
    Fperm <- vapply(seq_len(n_perm),
                    function(j) fit_shape(y, g_score, g_sample, g_bin, g_n, combos[, j])["F"],
                    0.0)
    Fobs <- Fperm[obs_j]
    rank <- sum(Fperm >= Fobs, na.rm = TRUE)       # 1 = most extreme
    data.frame(gene = g,
               asym_p_shape = obs$p[obs$gene == g],
               fs_edf       = obs$edf_fs[obs$gene == g],
               F_obs        = Fobs,
               perm_p       = rank / n_perm,
               perm_rank    = sprintf("%d/%d", rank, n_perm),
               stringsAsFactors = FALSE)
  })
  pergene <- bind_rows(pergene)
  write.csv(pergene, file.path(ct_dir, "pergene_permutation_ranks.csv"), row.names = FALSE)
  cat("\n  Per-gene permutation (top genes):\n")
  print(pergene, row.names = FALSE, digits = 3)

  # ==========================================================================
  # Diagnostic 2b — GLOBAL permutation test over ALL tested genes
  #   statistic per assignment = # genes at asymptotic FDR < SIG_FDR
  # ==========================================================================
  glob_genes <- if (is.na(GLOBAL_MAX_GENES)) test_genes
                else head(obs$gene, GLOBAL_MAX_GENES)
  Y <- t(vapply(glob_genes, gene_group_mean, numeric(length(gi))))  # genes x groups

  perm_counts <- integer(n_perm)
  for (j in seq_len(n_perm)) {
    pj <- apply(Y, 1, function(y)
                fit_shape(y, g_score, g_sample, g_bin, g_n, combos[, j])["p"])
    perm_counts[j] <- sum(p.adjust(pj, method = "BH") < SIG_FDR, na.rm = TRUE)
    if (j %% 5 == 0 || j == n_perm)
      cat(sprintf("    global perm %d/%d (count=%d)\n", j, n_perm, perm_counts[j]))
  }
  obs_global <- perm_counts[obs_j]
  global_p   <- sum(perm_counts >= obs_global) / n_perm

  glob_df <- data.frame(assignment = seq_len(n_perm),
                        wt_set = apply(combos, 2, paste, collapse = "+"),
                        n_genes_fdr_sig = perm_counts,
                        is_observed = seq_len(n_perm) == obs_j)
  write.csv(glob_df, file.path(ct_dir, "global_permutation_counts.csv"), row.names = FALSE)

  summ <- data.frame(
    cell_type = ct, n_perm = n_perm, perm_floor_p = 1 / n_perm,
    observed_n_genes_fdr_sig = obs_global,
    max_permuted_count = max(perm_counts),
    global_p = global_p)
  write.csv(summ, file.path(ct_dir, "global_permutation_summary.csv"), row.names = FALSE)
  cat(sprintf("\n  GLOBAL TEST: observed %d genes vs permuted max %d  ->  global p = %.3f\n",
              obs_global, max(perm_counts), global_p))

  # ==========================================================================
  # Plots
  # ==========================================================================
  # (1) Global null: histogram of permuted counts + observed marker
  p_glob <- ggplot(glob_df, aes(x = n_genes_fdr_sig)) +
    geom_histogram(binwidth = 1, fill = "grey70", color = "white") +
    geom_vline(xintercept = obs_global, color = CONDITION_COLORS[[TEST_CONDITION]],
               linewidth = 1.2) +
    annotate("text", x = obs_global, y = Inf,
             label = sprintf("observed = %d\nglobal p = %.3f", obs_global, global_p),
             hjust = 1.05, vjust = 1.5, size = 3.2) +
    labs(x = sprintf("# genes at asymptotic FDR < %.2f", SIG_FDR),
         y = sprintf("# of %d label assignments", n_perm),
         title = sprintf("%s — global sample-level permutation null", ct)) +
    theme_publication()
  ggsave(file.path(ct_dir, "global_permutation_null.pdf"), p_glob, width = 6, height = 4)

  # (2) Calibration: asymptotic vs permutation p for the reported genes
  p_cal <- ggplot(pergene, aes(x = -log10(pmax(asym_p_shape, .Machine$double.eps)),
                               y = perm_p, label = gene)) +
    geom_hline(yintercept = 1 / n_perm, linetype = "dashed", color = "grey50") +
    geom_point(size = 2, color = CONDITION_COLORS[[TEST_CONDITION]]) +
    { if (requireNamespace("ggrepel", quietly = TRUE))
        ggrepel::geom_text_repel(size = 3, seed = 1) else geom_text(size = 3, vjust = -0.6) } +
    labs(x = "-log10 asymptotic p_shape (mgcv)", y = "permutation p (sample-level)",
         title = sprintf("%s — asymptotic vs permutation calibration", ct),
         subtitle = sprintf("dashed = permutation floor 1/%d = %.3f", n_perm, 1 / n_perm)) +
    theme_publication()
  ggsave(file.path(ct_dir, "asymptotic_vs_permutation_p.pdf"), p_cal, width = 6, height = 4.5)

  cat(sprintf("  Done: %s -> %s\n", ct, ct_dir))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")