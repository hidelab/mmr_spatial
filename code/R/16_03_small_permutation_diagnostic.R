# 16_03_small_permutation_diagnostic.R
# ============================================================================
# Sample-level permutation calibration for the 15_04 tumor-trajectory condition
# test (reviewer point #4), for BOTH the SHAPE and the LEVEL question.
#
# ---------------------------------------------------------------------------
# WHY THIS EXISTS (the rationale, in full)
# ---------------------------------------------------------------------------
# 15_04 tests whether the WT->DKO expression trajectory along the tumor-proximity
# axis differs in SHAPE (`p_shape`, the mgcv ordered-factor difference smooth) and
# in overall LEVEL (`p_level`, the parametric condition offset). This script
# calibrates BOTH against the honest sample-level null, using the SAME spec 15_04
# now designates PRIMARY: `pseudobulk_nb_fs` — a NEGATIVE-BINOMIAL GAM on SUMMED
# RAW COUNTS per sample x bin + a log(library-size) OFFSET + a per-sample factor
# smooth s(score, sample, bs="fs"). (Earlier versions of this script calibrated the
# Gaussian `pseudobulk_fs`; the NB count likelihood is the current primary.)
#
# Two facts are BOTH true and must be held together:
#   (A) The MODEL is sound. The `fs` term gives each mouse its own penalized curve,
#       so between-mouse trajectory-SHAPE heterogeneity goes into the error term
#       (not attributed to condition). A variance-component check (the fs edf,
#       reported per gene) confirms the per-sample smooths are not collapsed.
#   (B) The mgcv ASYMPTOTIC p-values are NOT trustworthy in absolute magnitude at
#       n = 3 WT + 4 DKO. They are conditional on a between-sample variance
#       component barely identified from 7 mice; values like ~1e-16 claim a
#       resolution the design cannot support.
#
# The honest resolution is a SAMPLE-LEVEL PERMUTATION: there are choose(7,3) = 35
# ways to assign the condition label across the 7 mice, so any permutation p has a
# HARD FLOOR of 1/35 ~= 0.029, and BH across ~hundreds of genes makes per-gene
# FDR < 0.1 unreachable. The design supports (i) a GLOBAL statement and (ii) a
# per-gene RANKING — computed here for shape AND level:
#
#   Diagnostic 1 — variance component: fs edf per reported gene (model soundness).
#   Diagnostic 2a — per-gene permutation: refit under all 35 assignments; rank the
#     OBSERVED split's shape statistic (F of the difference smooth) and level
#     statistic (|z| of the condition offset) among them -> perm_p_shape/perm_p_level.
#   Diagnostic 2b — GLOBAL test: for each assignment, run the whole test over all
#     genes, BH within-assignment, count genes at FDR < SIG_FDR, SEPARATELY for shape
#     and level. If the observed assignment tops all 35, that is a global statement
#     at p = 1/35 = 0.029.
#
# COMMON SUPPORT: 15_04 requires MIN_SAMPLES_PER_COND_PER_BIN = 2. This script uses
# 1 (see the parameter note) and recomputes support INSIDE each permuted labeling
# (support depends on the labels, so it must be a deterministic function of them or
# the permutation test is invalid).
# ============================================================================

# --- Thread limits (set BEFORE BLAS/library load) ---------------------------
N_THREADS <- 4   # lower than usual: the global loop parallelizes over assignments
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
library(parallel)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_permutation_diagnostic")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters — synced to 15_04 (latest), except common support (see note)
# ============================================================================
CELLTYPES        <- c("CD8+ T")
TARGET_CELLTYPE  <- "Tumor"
MODEL_LABEL      <- "pseudobulk_nb_fs"   # documented: the spec this calibrates

# Score (local tumor fraction, no per-sample min-max — matches 15_04)
KNN_K            <- 30
DIFFUSION_VARIANT<- "heat"
DIFFUSION_HEAT_STEPS <- 20
DIFFUSION_HEAT_ALPHA <- 2
LOCAL_FRACTION   <- TRUE

# Pseudobulk / GAM (match 15_04 latest: fixed 5 bins, min5, GAM_K=5, FS_K=3, FS_M=2)
N_BINS            <- 5
MIN_CELLS_PER_BIN <- 5
GAM_K             <- 5
FS_K              <- 3
FS_M              <- 2
SIG_FDR           <- 0.1
PB_SCORE_AXIS     <- "score"   # "score" (real mean local fraction) | "decile"
APPLY_QC          <- TRUE      # 15_04 latest applies QC to trajectory cells

# COMMON SUPPORT: 15_04 uses 2; here 1 ON PURPOSE. This is the permutation analog
# of the LOSO argument (16_04) — under some of the 35 relabelings a condition arm
# is reduced, and requiring 2 samples/condition/bin would collapse the retained
# support for those labelings, breaking the test. 1 keeps every labeling fittable.
COMMON_SUPPORT               <- TRUE
MIN_SAMPLES_PER_COND_PER_BIN <- 1

# Per-gene permutation report: how many top genes (by observed asymptotic p) to
# profile in full, for EACH of shape and level (union of the two top lists).
PERGENE_REPORT_N  <- 20
# Cap genes entered into the (expensive) global loop. NA = all tested genes.
GLOBAL_MAX_GENES  <- NA
# Parallelism for the global permutation loop (over the 35 assignments).
N_PERM_CORES      <- max(1, min(12, detectCores() - 2))

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

# QC mask (per cell); all-TRUE when APPLY_QC = FALSE. Applied only to trajectory
# cell selection, NOT to the score field (matches 15_04).
qc_pass <- if (APPLY_QC) {
  seu$total_counts >= MIN_TOTAL_COUNTS & seu$n_genes_by_counts >= MIN_GENES_DETECTED
} else rep(TRUE, ncol(seu))
cat(sprintf("QC (APPLY_QC=%s): %d / %d cells pass\n", APPLY_QC, sum(qc_pass), length(qc_pass)))

cat("Computing local-tumor-fraction score...\n")
tumor_score <- compute_tumor_score(
  all_coords, is_tumor, samples, method = "diffusion",
  diffusion_method = DIFFUSION_VARIANT, k = KNN_K, heat_steps = DIFFUSION_HEAT_STEPS,
  alpha = DIFFUSION_HEAT_ALPHA,
  local_fraction = LOCAL_FRACTION, normalize = !LOCAL_FRACTION, verbose = FALSE)

# ============================================================================
# Core fit: the pseudobulk_nb_fs model (NB on summed counts + log-libsize offset +
# fs random smooth). Returns shape + level statistics/p-values + fs edf.
#   y_counts : per-group SUMMED raw counts for one gene
#   libsize  : per-group total counts (all genes) -> log offset
#   wt_set   : the set of samples labelled WT (permuted or observed)
# ============================================================================
fit_nbfs <- function(y_counts, libsize, g_score, g_sample, g_bin, wt_set) {
  cond <- ifelse(g_sample %in% wt_set, REFERENCE_CONDITION, TEST_CONDITION)
  keep <- rep(TRUE, length(cond))
  if (COMMON_SUPPORT) {
    bins <- unique(g_bin)
    ok <- vapply(bins, function(b) {
      m <- g_bin == b
      length(unique(g_sample[m & cond == REFERENCE_CONDITION])) >= MIN_SAMPLES_PER_COND_PER_BIN &&
      length(unique(g_sample[m & cond == TEST_CONDITION]))      >= MIN_SAMPLES_PER_COND_PER_BIN
    }, logical(1))
    keep <- g_bin %in% bins[ok]
  }
  na5 <- c(F_shape = NA, p_shape = NA, z_level = NA, p_level = NA, edf_fs = NA)
  if (sum(keep) < 4 || length(unique(cond[keep])) < 2) return(na5)
  d <- data.frame(expr = y_counts[keep],
                  log_libsize = log(pmax(libsize[keep], 1)),
                  score = g_score[keep],
                  sample = factor(g_sample[keep]),
                  cond_ord = ordered(cond[keep], levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
  f <- tryCatch(suppressWarnings(
         gam(expr ~ s(score, k = GAM_K) + s(score, by = cond_ord, k = GAM_K) +
                    cond_ord + s(score, sample, bs = "fs", k = FS_K, m = FS_M) +
                    offset(log_libsize),
             data = d, family = nb(), method = "REML")),
       error = function(e) NULL)
  if (is.null(f)) return(na5)
  sm <- summary(f)
  # Shape: the difference-smooth s.table row (":cond_ord"). col 3 = test stat, col 4 = p.
  rn <- rownames(sm$s.table); dr <- grep(":cond_ord", rn); fr <- grep("sample", rn)
  # Level: the parametric cond_ord row. Under nb() the p column is "Pr(>|z|)" and
  # the statistic column ends in "value" ("z value") -- grep both by position.
  pr <- rownames(sm$p.table); lr <- grep("cond_ord", pr)
  pcol <- grep("^Pr", colnames(sm$p.table)); zcol <- grep("value$", colnames(sm$p.table))
  c(F_shape = if (length(dr)) sm$s.table[dr[1], 3] else NA,
    p_shape = if (length(dr)) sm$s.table[dr[1], 4] else NA,
    z_level = if (length(lr) && length(zcol)) abs(sm$p.table[lr[1], zcol[1]]) else NA,
    p_level = if (length(lr) && length(pcol)) sm$p.table[lr[1], pcol[1]] else NA,
    edf_fs  = if (length(fr)) sm$s.table[fr[1], "edf"] else NA)
}

# ============================================================================
# Per cell type
# ============================================================================
for (ct in CELLTYPES) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
  ct_dir  <- file.path(OUTPUT_DIR, ct_safe)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n=== %s  [%s] ===\n", ct, MODEL_LABEL))

  ct_idx <- which(seu$Tier1_celltype == ct & !is.na(tumor_score) & qc_pass)
  if (length(ct_idx) < 50) { cat("  too few cells; skip\n"); next }

  ct_score <- tumor_score[ct_idx]
  ct_sample<- droplevels(seu$sample[ct_idx])
  counts   <- GetAssayData(seu, assay = "RNA", layer = "counts")[, ct_idx]

  # --- test genes: SAME baseline-expression filter as 15_04 / DE scripts ------
  det_frac  <- Matrix::rowSums(counts > 0) / length(ct_idx)
  samp_sum  <- vapply(levels(ct_sample), function(s)
                 Matrix::rowSums(counts[, ct_sample == s, drop = FALSE]),
                 numeric(nrow(counts)))
  n_samp_det <- rowSums(samp_sum >= MIN_COUNTS_PER_SAMPLE)
  test_genes <- rownames(counts)[det_frac >= MIN_GENE_EXPR_FRAC &
                                 n_samp_det >= MIN_SAMPLES_EXPRESSED]
  cat(sprintf("  %d tested genes (%d cells)\n", length(test_genes), length(ct_idx)))

  # --- build pseudobulk group structure ONCE (sample x bin) -------------------
  br <- quantile(ct_score, seq(0, 1, length.out = N_BINS + 1), na.rm = TRUE)
  br[1] <- br[1] - 1e-6
  bin <- as.integer(cut(ct_score, breaks = br, labels = seq_len(N_BINS)))
  grp <- paste(ct_sample, bin, sep = "__")
  gi  <- split(seq_along(grp), grp)
  gi  <- gi[vapply(gi, length, 0L) >= MIN_CELLS_PER_BIN]
  g_sample <- vapply(strsplit(names(gi), "__"), `[`, "", 1)
  g_bin    <- as.integer(vapply(strsplit(names(gi), "__"), `[`, "", 2))
  g_n      <- vapply(gi, length, 0L)
  g_score  <- if (PB_SCORE_AXIS == "score")
                vapply(gi, function(ix) mean(ct_score[ix]), 0.0)
              else (g_bin - 0.5) / N_BINS
  # NB pieces: per-group total library size (all genes), and summed counts / gene.
  g_libsize <- vapply(gi, function(ix) sum(counts[, ix]), 0.0)
  gene_group_sum <- function(g) vapply(gi, function(ix) sum(counts[g, ix]), 0.0)
  cat(sprintf("  %d pseudobulk points across %d samples (fixed %d bins, min %d cells)\n",
              length(gi), length(unique(g_sample)), N_BINS, MIN_CELLS_PER_BIN))

  # --- enumerate all label assignments (exhaustive) ---------------------------
  samps  <- levels(ct_sample)
  wt_here<- intersect(WT_SAMPLES, samps)
  n_wt   <- length(wt_here)
  combos <- combn(samps, n_wt)                       # each column = a WT set
  n_perm <- ncol(combos)
  obs_j  <- which(apply(combos, 2, function(x) setequal(x, wt_here)))
  cat(sprintf("  %d label assignments (choose(%d,%d)); permutation floor p = %.3f\n",
              n_perm, length(samps), n_wt, 1 / n_perm))

  # Summed-count matrix (genes x groups) — computed once; labels never change it.
  Y <- t(vapply(test_genes, gene_group_sum, numeric(length(gi))))
  rownames(Y) <- test_genes

  # helper: fit all genes under one WT set -> data.frame with both tests
  fit_all <- function(wt_set) {
    m <- t(apply(Y, 1, function(y)
             fit_nbfs(y, g_libsize, g_score, g_sample, g_bin, wt_set)))
    data.frame(gene = test_genes, F_shape = m[, "F_shape"], p_shape = m[, "p_shape"],
               z_level = m[, "z_level"], p_level = m[, "p_level"],
               edf_fs = m[, "edf_fs"], stringsAsFactors = FALSE)
  }

  # ==========================================================================
  # OBSERVED asymptotic pass (true labels), both tests
  # ==========================================================================
  obs <- fit_all(wt_here)
  obs$fdr_shape <- p.adjust(obs$p_shape, method = "BH")
  obs$fdr_level <- p.adjust(obs$p_level, method = "BH")
  obs <- obs[order(obs$p_shape), ]
  obs_count_shape <- sum(obs$fdr_shape < SIG_FDR, na.rm = TRUE)
  obs_count_level <- sum(obs$fdr_level < SIG_FDR, na.rm = TRUE)
  write.csv(obs, file.path(ct_dir, "observed_asymptotic.csv"), row.names = FALSE)
  cat(sprintf("  observed asymptotic FDR<%.2f: shape %d, level %d genes\n",
              SIG_FDR, obs_count_shape, obs_count_level))

  # ==========================================================================
  # Diagnostic 2b — GLOBAL permutation test (parallel over the 35 assignments)
  #   statistic per assignment = # genes at asymptotic FDR < SIG_FDR (shape/level)
  # ==========================================================================
  glob_genes <- if (is.na(GLOBAL_MAX_GENES)) test_genes else head(obs$gene, GLOBAL_MAX_GENES)
  Yg <- Y[glob_genes, , drop = FALSE]

  cat(sprintf("  global permutation over %d assignments x %d genes (%d cores)...\n",
              n_perm, length(glob_genes), N_PERM_CORES))
  perm_res <- mclapply(seq_len(n_perm), function(j) {
    wt_set <- combos[, j]
    m <- t(apply(Yg, 1, function(y)
             fit_nbfs(y, g_libsize, g_score, g_sample, g_bin, wt_set)))
    ps <- m[, "p_shape"]; pl <- m[, "p_level"]
    c(shape = sum(p.adjust(ps, "BH") < SIG_FDR, na.rm = TRUE),
      level = sum(p.adjust(pl, "BH") < SIG_FDR, na.rm = TRUE))
  }, mc.cores = N_PERM_CORES)
  perm_mat  <- do.call(rbind, perm_res)          # n_perm x c(shape, level)
  cnt_shape <- perm_mat[, "shape"]; cnt_level <- perm_mat[, "level"]

  obs_shape <- cnt_shape[obs_j]; obs_level <- cnt_level[obs_j]
  gp_shape  <- sum(cnt_shape >= obs_shape) / n_perm
  gp_level  <- sum(cnt_level >= obs_level) / n_perm

  glob_df <- data.frame(assignment = seq_len(n_perm),
                        wt_set = apply(combos, 2, paste, collapse = "+"),
                        n_shape_fdr_sig = cnt_shape,
                        n_level_fdr_sig = cnt_level,
                        is_observed = seq_len(n_perm) == obs_j)
  write.csv(glob_df, file.path(ct_dir, "global_permutation_counts.csv"), row.names = FALSE)

  summ <- data.frame(
    cell_type = ct, model = MODEL_LABEL, n_perm = n_perm, perm_floor_p = 1 / n_perm,
    observed_shape = obs_shape, max_permuted_shape = max(cnt_shape), global_p_shape = gp_shape,
    observed_level = obs_level, max_permuted_level = max(cnt_level), global_p_level = gp_level)
  write.csv(summ, file.path(ct_dir, "global_permutation_summary.csv"), row.names = FALSE)
  cat(sprintf("  GLOBAL shape: observed %d vs permuted max %d -> p = %.3f\n",
              obs_shape, max(cnt_shape), gp_shape))
  cat(sprintf("  GLOBAL level: observed %d vs permuted max %d -> p = %.3f\n",
              obs_level, max(cnt_level), gp_level))

  # ==========================================================================
  # Diagnostic 1 + 2a — per-gene permutation for the union of top-N shape & level
  # ==========================================================================
  report_genes <- union(head(obs$gene[order(obs$p_shape)], PERGENE_REPORT_N),
                        head(obs$gene[order(obs$p_level)], PERGENE_REPORT_N))
  # Per gene, refit under all 35 assignments; rank observed shape/level statistic.
  pergene <- bind_rows(mclapply(report_genes, function(g) {
    y <- Y[g, ]
    stats <- vapply(seq_len(n_perm), function(j)
               fit_nbfs(y, g_libsize, g_score, g_sample, g_bin, combos[, j])[c("F_shape", "z_level")],
               numeric(2))
    Fs <- stats["F_shape", ]; Zl <- stats["z_level", ]
    rank_shape <- sum(Fs >= Fs[obs_j], na.rm = TRUE)   # 1 = most extreme
    rank_level <- sum(Zl >= Zl[obs_j], na.rm = TRUE)
    data.frame(gene = g,
               asym_p_shape = obs$p_shape[obs$gene == g], fs_edf = obs$edf_fs[obs$gene == g],
               F_shape_obs = Fs[obs_j], perm_p_shape = rank_shape / n_perm,
               perm_rank_shape = sprintf("%d/%d", rank_shape, n_perm),
               asym_p_level = obs$p_level[obs$gene == g],
               z_level_obs = Zl[obs_j], perm_p_level = rank_level / n_perm,
               perm_rank_level = sprintf("%d/%d", rank_level, n_perm),
               stringsAsFactors = FALSE)
  }, mc.cores = N_PERM_CORES))
  pergene <- pergene[order(pergene$asym_p_shape), ]
  write.csv(pergene, file.path(ct_dir, "pergene_permutation_ranks.csv"), row.names = FALSE)
  cat("\n  Per-gene permutation (union of top shape & level genes):\n")
  print(head(pergene[, c("gene", "asym_p_shape", "perm_p_shape", "asym_p_level", "perm_p_level")], 15),
        row.names = FALSE, digits = 3)

  # ==========================================================================
  # Plots — global null (shape + level) and asymptotic-vs-permutation calibration
  # ==========================================================================
  mk_global_plot <- function(counts_vec, obs_val, gp, lab) {
    ggplot(data.frame(n = counts_vec), aes(x = n)) +
      geom_histogram(binwidth = 1, fill = "grey70", color = "white") +
      geom_vline(xintercept = obs_val, color = CONDITION_COLORS[[TEST_CONDITION]], linewidth = 1.2) +
      annotate("text", x = obs_val, y = Inf,
               label = sprintf("observed = %d\nglobal p = %.3f", obs_val, gp),
               hjust = 1.05, vjust = 1.5, size = 3.2) +
      labs(x = sprintf("# %s genes at asymptotic FDR < %.2f", lab, SIG_FDR),
           y = sprintf("# of %d assignments", n_perm),
           title = sprintf("%s — global %s permutation null [%s]", ct, lab, MODEL_LABEL)) +
      theme_publication()
  }
  ggsave(file.path(ct_dir, "global_permutation_null_shape.pdf"),
         mk_global_plot(cnt_shape, obs_shape, gp_shape, "shape"), width = 6, height = 4)
  ggsave(file.path(ct_dir, "global_permutation_null_level.pdf"),
         mk_global_plot(cnt_level, obs_level, gp_level, "level"), width = 6, height = 4)

  mk_cal_plot <- function(asym, perm, lab) {
    ggplot(data.frame(x = -log10(pmax(asym, .Machine$double.eps)), y = perm,
                      gene = pergene$gene),
           aes(x = x, y = y, label = gene)) +
      geom_hline(yintercept = 1 / n_perm, linetype = "dashed", color = "grey50") +
      geom_point(size = 2, color = CONDITION_COLORS[[TEST_CONDITION]]) +
      { if (requireNamespace("ggrepel", quietly = TRUE))
          ggrepel::geom_text_repel(size = 3, seed = 1, max.overlaps = 20)
        else geom_text(size = 3, vjust = -0.6) } +
      labs(x = sprintf("-log10 asymptotic p_%s (mgcv)", lab), y = "permutation p (sample-level)",
           title = sprintf("%s — %s: asymptotic vs permutation [%s]", ct, lab, MODEL_LABEL),
           subtitle = sprintf("dashed = floor 1/%d = %.3f", n_perm, 1 / n_perm)) +
      theme_publication()
  }
  ggsave(file.path(ct_dir, "asymptotic_vs_permutation_shape.pdf"),
         mk_cal_plot(pergene$asym_p_shape, pergene$perm_p_shape, "shape"), width = 6, height = 4.5)
  ggsave(file.path(ct_dir, "asymptotic_vs_permutation_level.pdf"),
         mk_cal_plot(pergene$asym_p_level, pergene$perm_p_level, "level"), width = 6, height = 4.5)

  cat(sprintf("  Done: %s -> %s\n", ct, ct_dir))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")
