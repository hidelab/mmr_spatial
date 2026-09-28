# 16_04_loso_influence.R
# ============================================================================
# Leave-one-sample-out (LOSO) influence analysis for the 15_04 tumor-trajectory
# condition test, for BOTH the SHAPE and the LEVEL question.
#
# ---------------------------------------------------------------------------
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The sample-level permutation (16_03) answers "is the observed WT/DKO split
# special vs arbitrary regroupings of these 7 mice?" — the correct POPULATION
# inference, but bounded by a 1/choose(7,3)=1/35~=0.029 floor, and with only 35
# regroupings a single leverage mouse can dominate the statistic.
#
# The DIRECT question "does one sample drive the findings?" is about INFLUENCE, not
# the null. LOSO answers it: drop each mouse in turn (keeping the TRUE WT/DKO
# labels), refit the whole pipeline on the remaining 6, and see which shape- and
# level-significant genes survive. A result hinging on a single mouse is not robust
# regardless of any p-value.
#
# This uses the SAME spec 15_04 designates PRIMARY and that 16_03 calibrates:
# `pseudobulk_nb_fs` (NB on summed raw counts + log(library-size) offset + a
# per-sample factor smooth). The full-fit baseline reproduces the reported gene set.
#
# COMMON SUPPORT: the FULL-fit baseline uses MIN_SAMPLES_PER_COND_PER_BIN + 1 = 2,
# matching 15_04 exactly, so the baseline IS the reported gene set — we test the
# robustness of the actual claim, not of a looser baseline we never reported.
#
# LOSO BINNING (LOSO_INHERIT_FULL_BINS = TRUE, default): each drop FREEZES the full
# model's binning — it inherits BOTH (a) the full model's quantile cutoffs (so a cell
# lands in the same bin regardless of which mouse is dropped) AND (b) the full model's
# SET of common-support-passing bins (the min=2 baseline support region). A drop then
# only removes the dropped mouse's cells from those frozen bins and refits. A bin can
# therefore only LEAVE a drop's fit — because removing that mouse pushed it below
# MIN_CELLS_PER_BIN (genuine influence) — and can NEVER ENTER one.
#
# WHY FREEZE THE BIN SET rather than re-run the rule at a looser min_samp: consider a
# tumor-core bin the min=2 baseline EXCLUDED because only DKO2 had cells there among
# the DKO mice (1 DKO sample < 2). If a drop re-ran common support at min=1, dropping
# any OTHER mouse (say DKO1) would suddenly let that bin QUALIFY (1 DKO now suffices)
# and ENTER the drop's fit. That is doubly wrong: (i) the drop's support region no
# longer matches the baseline it is meant to perturb, and (ii) the added bin's DKO arm
# rests on a SINGLE mouse (DKO2) — exactly the pseudoreplication the >=2-samples rule
# and the fs random effect exist to prevent, and tied to the leverage mouse we probe.
# Freezing the passing-bin set eliminates this. It also removes the reason drops could
# not simply reuse min=2 (re-applying min=2 to the reduced data would collapse support)
# — we never re-apply the rule to a drop, we only remove cells. (Set the toggle FALSE
# to recompute edges AND re-run common support per drop, the old anticonservative mode.)
#
# READING THE OUTPUT: the asymptotic mgcv FDRs are miscalibrated in MAGNITUDE at
# this sample size (see 16_03) — LOSO reads the STABILITY of gene membership across
# the 7 drops, not calibrated significance.
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
library(tidyr)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_loso_influence")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters — synced to 15_04 / 16_03, except common support (see note above)
# ============================================================================
CELLTYPES        <- c("CD8+ T")
TARGET_CELLTYPE  <- "Tumor"
MODEL_LABEL      <- "pseudobulk_nb_fs"

# Score (local tumor fraction, no per-sample min-max)
KNN_K            <- 30
DIFFUSION_VARIANT<- "heat"
DIFFUSION_HEAT_STEPS <- 20
DIFFUSION_HEAT_ALPHA <- 2
LOCAL_FRACTION   <- TRUE

# Pseudobulk / GAM (match 15_04 latest)
N_BINS            <- 10
MIN_CELLS_PER_BIN <- 5
GAM_K             <- 5
FS_K              <- 3
FS_M              <- 2
SIG_FDR           <- 0.1
PB_SCORE_AXIS     <- "score"
APPLY_QC          <- TRUE
COMMON_SUPPORT               <- TRUE
MIN_SAMPLES_PER_COND_PER_BIN <- 1     # only the fallback when LOSO_INHERIT_FULL_BINS=FALSE
# FULL-fit baseline uses +1 (=2) to reproduce 15_04 exactly and to DEFINE the frozen
# baseline support region that the LOSO drops inherit (see header).
MIN_SAMPLES_FULL             <- MIN_SAMPLES_PER_COND_PER_BIN + 1L

# How many top genes (by full-fit p) to show in each influence heatmap.
HEATMAP_TOP_N     <- 50

# LOSO drops INHERIT the FULL model's binning (TRUE): BOTH the bin cutoffs AND the set
# of common-support-passing bins are frozen from the full (min=2) fit, so a drop only
# removes the dropped mouse's cells from those fixed bins — a bin can leave (falls below
# MIN_CELLS_PER_BIN) but never enter. This is the clean influence test; see header.
# FALSE = recompute equal-count quantile edges AND re-run common support per drop.
LOSO_INHERIT_FULL_BINS <- TRUE

# ============================================================================
# Load + score (identical recoding to 15_04 / 16_03)
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
# Core fit: pseudobulk_nb_fs (NB summed counts + log-libsize offset + fs). Returns
# shape + level statistics/p-values + fs edf. Labels are the TRUE WT/DKO labels
# (LOSO does not permute); `wt_set` is the retained samples' WT labels.
# ============================================================================
fit_nbfs <- function(y_counts, libsize, g_score, g_sample, g_bin, wt_set,
                     min_samp = MIN_SAMPLES_PER_COND_PER_BIN, allowed_bins = NULL) {
  cond <- ifelse(g_sample %in% wt_set, REFERENCE_CONDITION, TEST_CONDITION)
  keep <- rep(TRUE, length(cond))
  if (!is.null(allowed_bins)) {
    # LOSO drops: restrict to the FULL model's frozen common-support-passing bins.
    # A bin can only LEAVE a drop (its cells fell below MIN_CELLS_PER_BIN -> it is
    # absent from g_bin) and can NEVER ENTER. Re-running common support at a looser
    # min_samp would let a single-sample-dominated bin the min=2 baseline excluded
    # (e.g. a DKO2-only tumor-core bin) qualify and enter a drop of some OTHER mouse
    # — changing the support region AND reintroducing single-mouse pseudoreplication.
    keep <- g_bin %in% allowed_bins
  } else if (COMMON_SUPPORT) {
    bins <- unique(g_bin)
    ok <- vapply(bins, function(b) {
      m <- g_bin == b
      length(unique(g_sample[m & cond == REFERENCE_CONDITION])) >= min_samp &&
      length(unique(g_sample[m & cond == TEST_CONDITION]))      >= min_samp
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
  rn <- rownames(sm$s.table); dr <- grep(":cond_ord", rn); fr <- grep("sample", rn)
  pr <- rownames(sm$p.table); lr <- grep("cond_ord", pr)
  pcol <- grep("^Pr", colnames(sm$p.table)); zcol <- grep("value$", colnames(sm$p.table))
  c(F_shape = if (length(dr)) sm$s.table[dr[1], 3] else NA,
    p_shape = if (length(dr)) sm$s.table[dr[1], 4] else NA,
    z_level = if (length(lr) && length(zcol)) abs(sm$p.table[lr[1], zcol[1]]) else NA,
    p_level = if (length(lr) && length(pcol)) sm$p.table[lr[1], pcol[1]] else NA,
    edf_fs  = if (length(fr)) sm$s.table[fr[1], "edf"] else NA)
}

# ----------------------------------------------------------------------------
# Build the sample x bin pseudobulk structure from a chosen set of cells (NB:
# per-group summed counts + total library size). Bin EDGES are inherited when `breaks`
# is supplied (LOSO drops) or computed as equal-count quantiles otherwise (full model);
# MIN_CELLS_PER_BIN is always applied on the RETAINED cells. Common support is NOT
# decided here — it is handled in fit_nbfs (full: min_samp rule; drops: frozen bin set).
# Returns group indices (into the ct-cell vector) plus the bin edges used.
# ----------------------------------------------------------------------------
build_pb <- function(ct_score, ct_sample, keep_cells, counts, breaks = NULL) {
  sc  <- ct_score[keep_cells]
  sm  <- droplevels(ct_sample[keep_cells])
  if (is.null(breaks)) {
    # master edges: equal-count quantile bins on the supplied cells (the full model)
    br  <- quantile(sc, seq(0, 1, length.out = N_BINS + 1), na.rm = TRUE)
    br[1] <- br[1] - 1e-6
    if (any(duplicated(br))) br <- unique(br)
  } else {
    # inherit fixed edges: a LOSO drop reuses the full model's cutoffs so every cell
    # lands in the same bin regardless of which mouse is dropped. Retained scores are
    # a subset of the full range, so all fall within these edges (no out-of-range NA).
    br <- breaks
  }
  bn  <- as.integer(cut(sc, breaks = br, labels = seq_len(length(br) - 1)))
  grp <- paste(as.character(sm), bn, sep = "__")
  gi  <- split(keep_cells, grp)
  gi  <- gi[vapply(gi, length, 0L) >= MIN_CELLS_PER_BIN]
  if (!length(gi)) return(NULL)
  g_sample  <- vapply(strsplit(names(gi), "__"), `[`, "", 1)
  g_bin     <- as.integer(vapply(strsplit(names(gi), "__"), `[`, "", 2))
  g_n       <- vapply(gi, length, 0L)
  g_score   <- if (PB_SCORE_AXIS == "score")
                 vapply(gi, function(ix) mean(ct_score[ix]), 0.0)
               else (g_bin - 0.5) / N_BINS
  g_libsize <- vapply(gi, function(ix) sum(counts[, ix]), 0.0)
  list(gi = gi, g_sample = g_sample, g_bin = g_bin, g_n = g_n,
       g_score = g_score, g_libsize = g_libsize, breaks = br)
}

# Fit all genes for one pseudobulk structure -> df(gene, F_shape, p_shape, z_level,
# p_level, edf_fs, fdr_shape, fdr_level). Uses SUMMED raw counts per group.
# Set of bins passing common support (>= min_samp distinct samples in BOTH conditions),
# gene-independent. Used to freeze the full model's baseline support region for LOSO.
passing_bins <- function(g_bin, g_sample, wt_set, min_samp) {
  cond <- ifelse(g_sample %in% wt_set, REFERENCE_CONDITION, TEST_CONDITION)
  bins <- sort(unique(g_bin))
  bins[vapply(bins, function(b) {
    m <- g_bin == b
    length(unique(g_sample[m & cond == REFERENCE_CONDITION])) >= min_samp &&
    length(unique(g_sample[m & cond == TEST_CONDITION]))      >= min_samp
  }, logical(1))]
}

fit_all_genes <- function(pb, counts, test_genes, wt_set,
                          min_samp = MIN_SAMPLES_PER_COND_PER_BIN, allowed_bins = NULL) {
  gsum <- function(g) vapply(pb$gi, function(ix) sum(counts[g, ix]), 0.0)
  m <- t(vapply(test_genes, function(g)
           fit_nbfs(gsum(g), pb$g_libsize, pb$g_score, pb$g_sample, pb$g_bin, wt_set,
                    min_samp, allowed_bins),
           numeric(5)))
  df <- data.frame(gene = test_genes, F_shape = m[, 1], p_shape = m[, 2],
                   z_level = m[, 3], p_level = m[, 4], edf_fs = m[, 5],
                   stringsAsFactors = FALSE)
  df$fdr_shape <- p.adjust(df$p_shape, method = "BH")
  df$fdr_level <- p.adjust(df$p_level, method = "BH")
  df
}

# ----------------------------------------------------------------------------
# Per-model bin report: one row per pseudobulk group (sample x bin) for a given
# model (full or a LOSO drop). Records the cells analyzed, cells-per-bin, bin
# edges, library size, condition label, and whether the bin passed common support
# (gene-independent, so computable once per model). Per-model totals are repeated
# on every row for easy filtering/aggregation. All models dumped to one CSV.
# ----------------------------------------------------------------------------
pb_report <- function(pb, model_label, dropped, wt_set, min_samp, allowed_bins = NULL) {
  cond_grp <- ifelse(pb$g_sample %in% wt_set, REFERENCE_CONDITION, TEST_CONDITION)
  if (!is.null(allowed_bins)) {
    # drops: a bin "enters the GAM" iff it is in the frozen baseline set (and present)
    pass_vec <- pb$g_bin %in% allowed_bins
  } else {
    # full model: compute common support at this min_samp (this defines the baseline)
    bins <- sort(unique(pb$g_bin))
    cs_pass <- vapply(bins, function(b) {
      m <- pb$g_bin == b
      length(unique(pb$g_sample[m & cond_grp == REFERENCE_CONDITION])) >= min_samp &&
      length(unique(pb$g_sample[m & cond_grp == TEST_CONDITION]))      >= min_samp
    }, logical(1))
    names(cs_pass) <- as.character(bins)
    pass_vec <- unname(cs_pass[as.character(pb$g_bin)])
  }
  data.frame(
    model                = model_label,
    dropped_sample       = dropped,
    min_samp             = min_samp,
    sample               = pb$g_sample,
    condition            = cond_grp,
    bin                  = pb$g_bin,
    bin_lo               = pb$breaks[pb$g_bin],
    bin_hi               = pb$breaks[pb$g_bin + 1],
    bin_score_mean       = pb$g_score,
    n_cells              = pb$g_n,               # cells in this sample x bin group
    libsize              = pb$g_libsize,         # total counts (all genes) in group
    bin_common_support   = as.logical(pass_vec), # did this bin enter the GAM?
    # per-model totals (repeated per row)
    n_samples_in_model   = length(unique(pb$g_sample)),
    n_bins_defined       = length(pb$breaks) - 1L,
    n_pb_points_in_model = length(pb$gi),
    n_pb_points_analyzed = sum(pass_vec),
    n_cells_in_model     = sum(pb$g_n),
    n_cells_analyzed     = sum(pb$g_n[pass_vec]),
    stringsAsFactors     = FALSE
  )
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

  # --- test genes: SAME baseline-expression filter as 15_04 / 16_03 ----------
  det_frac  <- Matrix::rowSums(counts > 0) / length(ct_idx)
  samp_sum  <- vapply(levels(ct_sample), function(s)
                 Matrix::rowSums(counts[, ct_sample == s, drop = FALSE]),
                 numeric(nrow(counts)))
  n_samp_det <- rowSums(samp_sum >= MIN_COUNTS_PER_SAMPLE)
  test_genes <- rownames(counts)[det_frac >= MIN_GENE_EXPR_FRAC &
                                 n_samp_det >= MIN_SAMPLES_EXPRESSED]
  cat(sprintf("  %d tested genes (%d cells across %d samples)\n",
              length(test_genes), length(ct_idx), nlevels(ct_sample)))

  samps    <- levels(ct_sample)
  ct_local <- seq_along(ct_idx)
  wt_here  <- intersect(WT_SAMPLES, samps)

  # ==========================================================================
  # FULL fit (all samples) — the baseline gene set
  # ==========================================================================
  pb_full   <- build_pb(ct_score, ct_sample, ct_local, counts)
  full      <- fit_all_genes(pb_full, counts, test_genes, wt_here, min_samp = MIN_SAMPLES_FULL)
  full      <- full[order(full$p_shape), ]
  sig_shape <- full$gene[which(full$fdr_shape < SIG_FDR)]
  sig_level <- full$gene[which(full$fdr_level < SIG_FDR)]
  cat(sprintf("  FULL (min_samp=%d, matches 15_04): %d pb points; asymptotic FDR<%.2f -> shape %d, level %d genes\n",
              MIN_SAMPLES_FULL, length(pb_full$gi), SIG_FDR, length(sig_shape), length(sig_level)))

  # Frozen baseline support region: the set of bins passing common support in the FULL
  # (min=2) fit. When LOSO_INHERIT_FULL_BINS, every drop is restricted to THESE bins —
  # a drop can lose a bin (cells < MIN_CELLS_PER_BIN) but never gain one. NULL when the
  # toggle is off -> drops re-run common support at MIN_SAMPLES_PER_COND_PER_BIN.
  baseline_bins <- if (LOSO_INHERIT_FULL_BINS)
    passing_bins(pb_full$g_bin, pb_full$g_sample, wt_here, MIN_SAMPLES_FULL) else NULL
  cat(sprintf("  frozen baseline bins (min_samp=%d): {%s} of %d defined\n",
              MIN_SAMPLES_FULL,
              if (is.null(baseline_bins)) "per-drop re-eval" else paste(sort(baseline_bins), collapse = ","),
              length(pb_full$breaks) - 1L))

  # accumulate per-model bin reports (full first, then each drop)
  report_rows <- list(full = pb_report(pb_full, "full", "(none)", wt_here, MIN_SAMPLES_FULL))

  # ==========================================================================
  # LOSO: drop each sample, refit everything on the remaining 6
  # ==========================================================================
  loso <- list()
  cnt_shape <- integer(length(samps)); names(cnt_shape) <- samps
  cnt_level <- integer(length(samps)); names(cnt_level) <- samps
  for (s in samps) {
    keep_cells <- ct_local[ct_sample != s]
    pb_s <- build_pb(ct_score, ct_sample, keep_cells, counts,
                     breaks = if (LOSO_INHERIT_FULL_BINS) pb_full$breaks else NULL)
    if (is.null(pb_s)) { cat(sprintf("  drop %s: no pb points; skip\n", s)); next }
    wt_s  <- setdiff(wt_here, s)
    fit_s <- fit_all_genes(pb_s, counts, test_genes, wt_s,
                           allowed_bins = baseline_bins)
    rownames(fit_s) <- fit_s$gene
    loso[[s]] <- fit_s
    report_rows[[s]] <- pb_report(pb_s, paste0("drop_", s), s, wt_s,
                                  MIN_SAMPLES_PER_COND_PER_BIN, allowed_bins = baseline_bins)
    cnt_shape[s] <- sum(fit_s$fdr_shape < SIG_FDR, na.rm = TRUE)
    cnt_level[s] <- sum(fit_s$fdr_level < SIG_FDR, na.rm = TRUE)
    cat(sprintf("  drop %-5s (%s): %d pb points, shape %d, level %d genes FDR<%.2f\n",
                s, if (s %in% wt_here) "WT" else "DKO",
                length(pb_s$gi), cnt_shape[s], cnt_level[s], SIG_FDR))
  }

  # ==========================================================================
  # Per-model bin report — one file, all models (full + each drop)
  # ==========================================================================
  model_order  <- c("full", paste0("drop_", samps))
  model_report <- bind_rows(report_rows)
  model_report <- model_report[order(match(model_report$model, model_order),
                                     model_report$bin, model_report$sample), ]
  write.csv(model_report, file.path(ct_dir, "model_bin_report.csv"), row.names = FALSE)
  # compact per-model summary (one row per model) alongside the granular table
  model_summary <- model_report %>%
    group_by(model, dropped_sample, min_samp) %>%
    summarise(n_samples = dplyr::first(n_samples_in_model),
              n_bins_defined = dplyr::first(n_bins_defined),
              n_pb_points = dplyr::first(n_pb_points_in_model),
              n_pb_points_analyzed = dplyr::first(n_pb_points_analyzed),
              n_cells = dplyr::first(n_cells_in_model),
              n_cells_analyzed = dplyr::first(n_cells_analyzed),
              min_cells_per_bin = min(n_cells), max_cells_per_bin = max(n_cells),
              .groups = "drop") %>%
    arrange(match(model, model_order))
  write.csv(model_summary, file.path(ct_dir, "model_bin_summary.csv"), row.names = FALSE)
  cat(sprintf("  wrote bin report: model_bin_report.csv (%d rows), model_bin_summary.csv (%d models)\n",
              nrow(model_report), nrow(model_summary)))
  print(as.data.frame(model_summary), row.names = FALSE)

  # ==========================================================================
  # Raw per-gene statistics for EVERY model — one CSV per model in a subdirectory.
  # Each file: gene, F_shape, p_shape, z_level, p_level, edf_fs, fdr_shape, fdr_level
  # (the untouched fit output; sorted by p_shape). Plus one combined long table.
  # ==========================================================================
  pval_dir <- file.path(ct_dir, "per_model_pvalues")
  dir.create(pval_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(full[order(full$p_shape), ],
            file.path(pval_dir, "pvalues_full.csv"), row.names = FALSE)
  for (s in names(loso)) {
    fs <- loso[[s]][order(loso[[s]]$p_shape), ]
    write.csv(fs, file.path(pval_dir, sprintf("pvalues_drop_%s.csv", s)), row.names = FALSE)
  }
  # combined long file (model column) for convenience
  pval_long <- bind_rows(
    cbind(model = "full", dropped_sample = "(none)", full),
    bind_rows(lapply(names(loso), function(s)
      cbind(model = paste0("drop_", s), dropped_sample = s, loso[[s]]))))
  write.csv(pval_long, file.path(pval_dir, "pvalues_all_models_long.csv"), row.names = FALSE)
  cat(sprintf("  wrote raw per-gene stats for %d models -> %s\n",
              length(loso) + 1L, pval_dir))

  # ==========================================================================
  # Rank-correlation matrices across models (full + each drop): how stable is
  # the whole-panel gene RANKING when each mouse is dropped? Spearman rho on the
  # per-gene p-values, computed separately for p_shape and p_level. rho ~ 1 =>
  # dropping that mouse barely reorders the panel; a low row/column flags a mouse
  # whose removal reshuffles the ranking (leverage on the ORDERING — complementary
  # to the membership-based influence tables below, which read FDR<cutoff in/out).
  # Spearman is monotone-invariant, so raw p vs -log10(p) give identical matrices.
  # ==========================================================================
  rank_corr <- function(pcol, lab) {
    genes  <- full$gene
    mat <- cbind(full = full[[pcol]][match(genes, full$gene)],
                 vapply(names(loso), function(s) loso[[s]][genes, pcol],
                        numeric(length(genes))))
    colnames(mat) <- c("full", names(loso))
    cmat <- cor(mat, method = "spearman", use = "pairwise.complete.obs")
    write.csv(cmat, file.path(ct_dir, sprintf("loso_rank_corr_%s.csv", lab)))
    # order models as full, then each dropped sample in SAMPLE_COLORS order
    ord  <- c("full", intersect(samps, colnames(cmat)))
    long <- as.data.frame(as.table(cmat), stringsAsFactors = FALSE)
    names(long) <- c("m1", "m2", "rho")
    long$m1 <- factor(long$m1, levels = ord)
    long$m2 <- factor(long$m2, levels = rev(ord))
    lo <- min(cmat, na.rm = TRUE)
    p <- ggplot(long, aes(x = m1, y = m2, fill = rho)) +
      geom_tile(color = "white") +
      geom_text(aes(label = sprintf("%.2f", rho)), size = 2.6) +
      scale_fill_gradient(low = "grey92", high = CONDITION_COLORS[[TEST_CONDITION]],
                          limits = c(lo, 1), name = "Spearman\nrho") +
      coord_equal() +
      labs(x = NULL, y = NULL,
           title = sprintf("%s — LOSO rank-correlation of %s p-values [%s]", ct, lab, MODEL_LABEL),
           subtitle = "Spearman rho of the per-gene p across models; low row/col => that mouse reorders the panel") +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(ct_dir, sprintf("loso_rank_corr_%s.pdf", lab)), p,
           width = 6.5, height = 6)
    cat(sprintf("  wrote rank-corr matrix (%s): min rho = %.3f -> loso_rank_corr_%s.{csv,pdf}\n",
                lab, lo, lab))
  }
  rank_corr("p_shape", "shape")
  rank_corr("p_level", "level")

  # ==========================================================================
  # Influence tables (one per test): per gene, FDR in full + each drop + robustness
  # ==========================================================================
  gene_order <- full$gene
  make_influence <- function(fdr_col) {
    tab <- data.frame(gene = gene_order,
                      p_full   = full[[sub("fdr", "p", fdr_col)]][match(gene_order, full$gene)],
                      fdr_full = full[[fdr_col]][match(gene_order, full$gene)],
                      stringsAsFactors = FALSE)
    for (s in names(loso)) tab[[paste0("fdr_drop_", s)]] <- loso[[s]][gene_order, fdr_col]
    drop_cols <- grep("^fdr_drop_", names(tab), value = TRUE)
    tab$n_subsets_sig      <- rowSums(tab[, drop_cols, drop = FALSE] < SIG_FDR, na.rm = TRUE)
    tab$n_subsets          <- length(drop_cols)
    tab$max_fdr_over_drops <- apply(tab[, drop_cols, drop = FALSE], 1, max, na.rm = TRUE)
    tab$sig_full           <- tab$fdr_full < SIG_FDR
    # robust = significant in the full fit AND in every LOSO drop
    tab$robust             <- tab$sig_full & (tab$n_subsets_sig == tab$n_subsets)
    tab
  }
  tab_shape <- make_influence("fdr_shape")
  tab_level <- make_influence("fdr_level")
  write.csv(tab_shape, file.path(ct_dir, "loso_influence_table_shape.csv"), row.names = FALSE)
  write.csv(tab_level, file.path(ct_dir, "loso_influence_table_level.csv"), row.names = FALSE)

  # ==========================================================================
  # Summary (both tests)
  # ==========================================================================
  kills <- function(tab) vapply(names(loso), function(s)
    sum(tab$sig_full & !(tab[[paste0("fdr_drop_", s)]] < SIG_FDR), na.rm = TRUE), integer(1))
  k_shape <- kills(tab_shape); k_level <- kills(tab_level)
  summ <- data.frame(
    cell_type = ct, model = MODEL_LABEL, n_tested = length(test_genes),
    test = c("shape", "level"),
    n_sig_full = c(length(sig_shape), length(sig_level)),
    n_robust_all_drops = c(sum(tab_shape$robust, na.rm = TRUE), sum(tab_level$robust, na.rm = TRUE)),
    most_influential_sample = c(names(k_shape)[which.max(k_shape)], names(k_level)[which.max(k_level)]),
    genes_killed_by_most_influential = c(max(k_shape), max(k_level)),
    min_count_over_drops = c(min(cnt_shape[names(loso)]), min(cnt_level[names(loso)])),
    max_count_over_drops = c(max(cnt_shape[names(loso)]), max(cnt_level[names(loso)])),
    stringsAsFactors = FALSE)
  summ$frac_robust <- round(summ$n_robust_all_drops / pmax(1, summ$n_sig_full), 3)
  write.csv(summ, file.path(ct_dir, "loso_summary.csv"), row.names = FALSE)

  # per-drop count table (long: one row per dropped sample x test)
  cnt_df <- bind_rows(
    data.frame(dropped_sample = c("(none: full)", names(cnt_shape)),
               condition = c("-", ifelse(names(cnt_shape) %in% wt_here, "WT", "DKO")),
               test = "shape",
               n_genes_fdr_sig = c(length(sig_shape), as.integer(cnt_shape)),
               genes_killed_vs_full = c(0L, as.integer(k_shape[names(cnt_shape)]))),
    data.frame(dropped_sample = c("(none: full)", names(cnt_level)),
               condition = c("-", ifelse(names(cnt_level) %in% wt_here, "WT", "DKO")),
               test = "level",
               n_genes_fdr_sig = c(length(sig_level), as.integer(cnt_level)),
               genes_killed_vs_full = c(0L, as.integer(k_level[names(cnt_level)]))))
  write.csv(cnt_df, file.path(ct_dir, "loso_counts_by_dropped_sample.csv"), row.names = FALSE)

  cat(sprintf("\n  SUMMARY shape: %d/%d survive dropping ANY mouse; most influential = %s (kills %d)\n",
              sum(tab_shape$robust, na.rm = TRUE), length(sig_shape),
              names(k_shape)[which.max(k_shape)], max(k_shape)))
  cat(sprintf("  SUMMARY level: %d/%d survive dropping ANY mouse; most influential = %s (kills %d)\n",
              sum(tab_level$robust, na.rm = TRUE), length(sig_level),
              names(k_level)[which.max(k_level)], max(k_level)))
  print(cnt_df, row.names = FALSE)

  # ==========================================================================
  # Plots
  # ==========================================================================
  # (1) sig-gene count per dropped sample, faceted by test
  lvl <- c("(none: full)", names(cnt_shape))
  cnt_df$dropped_sample <- factor(cnt_df$dropped_sample, levels = lvl)
  full_counts <- c(shape = length(sig_shape), level = length(sig_level))
  hlines <- data.frame(test = names(full_counts), yint = as.integer(full_counts))
  p_cnt <- ggplot(cnt_df, aes(x = dropped_sample, y = n_genes_fdr_sig, fill = condition)) +
    geom_col() +
    geom_hline(data = hlines, aes(yintercept = yint), linetype = "dashed", color = "grey40") +
    geom_text(aes(label = n_genes_fdr_sig), vjust = -0.4, size = 2.8) +
    facet_wrap(~ test, nrow = 1, scales = "free_y") +
    labs(x = "dropped sample", y = sprintf("# genes at asymptotic FDR < %.2f", SIG_FDR),
         title = sprintf("%s — LOSO significant-gene count [%s]", ct, MODEL_LABEL),
         subtitle = "dashed = full-fit count; a large drop => that mouse drives the result") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(ct_dir, "loso_gene_counts.pdf"), p_cnt, width = 9, height = 4.5)

  # (2) influence heatmaps: top genes x {full, each drop}, cell = -log10 FDR, x = not sig
  mk_heat <- function(tab, sig_genes, lab) {
    top <- head(sig_genes, HEATMAP_TOP_N)
    if (length(top) < 2) { cat(sprintf("  [%s] <2 sig genes; heatmap skipped\n", lab)); return(invisible()) }
    drop_cols <- grep("^fdr_drop_", names(tab), value = TRUE)
    M <- as.matrix(tab[match(top, tab$gene), c("fdr_full", drop_cols)])
    rownames(M) <- top; colnames(M) <- c("full", sub("^fdr_drop_", "-", drop_cols))
    long <- expand.grid(gene = rownames(M), drop = colnames(M), stringsAsFactors = FALSE)
    long$fdr    <- M[cbind(match(long$gene, rownames(M)), match(long$drop, colnames(M)))]
    long$neglog <- pmin(-log10(pmax(long$fdr, .Machine$double.eps)), 6)
    long$sig    <- long$fdr < SIG_FDR
    long$gene   <- factor(long$gene, levels = rev(top))
    long$drop   <- factor(long$drop, levels = colnames(M))
    p <- ggplot(long, aes(x = drop, y = gene, fill = neglog)) +
      geom_tile(color = "white") +
      geom_text(data = subset(long, !sig), aes(label = "x"), size = 2.4, color = "grey30") +
      scale_fill_gradient(low = "grey92", high = CONDITION_COLORS[[TEST_CONDITION]],
                          name = "-log10 FDR\n(cap 6)") +
      labs(x = "leave-one-out (full, then drop each sample)", y = NULL,
           title = sprintf("%s — LOSO influence on %s genes [%s]", ct, lab, MODEL_LABEL),
           subtitle = "x = not significant when that mouse is dropped") +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 6))
    ggsave(file.path(ct_dir, sprintf("loso_influence_heatmap_%s.pdf", lab)), p,
           width = 7, height = max(4, 0.16 * length(top) + 1.5), limitsize = FALSE)
  }
  mk_heat(tab_shape, sig_shape, "shape")
  mk_heat(tab_level, sig_level, "level")

  cat(sprintf("  Done: %s -> %s\n", ct, ct_dir))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")
