# 15_04_tumor_trajectory_working.R
# Working tumor-proximity trajectory script.
#
# Improvements over 15_01/15_02/15_03/16_01:
#   * Tumor score comes from the shared compute_tumor_score() in 00_utils.R, so
#     the scoring method (KDE / k-NN / diffusion) is a one-line switch.
#   * Sample-aware GAMs address pseudoreplication (3 WT + 4 DKO): the condition
#     test no longer treats individual cells as independent replicates.
#   * Multiple GAM specifications can be run and compared in one pass.
#   * Binning-quality diagnostics for the pseudobulk approach.
#   * A normalized-expression-along-gradient ("wave") plot.
#
# Starts with the KDE score. Outputs -> output/R/15_tumor_trajectory_working.

# --- Thread limits (set BEFORE BLAS/library load) ----------------------------
# gam() fits pull in multithreaded BLAS; cap it so we don't hog the shared box.
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
library(dplyr)
library(tidyr)
library(mgcv)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# Source centralized parameters, palette, and scoring utilities
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "15_tumor_trajectory_working")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters (all switchable)
# =============================================================================
# --- Tumor score -------------------------------------------------------------
SCORE_METHOD  <- "diffusion"  # "kde" | "knn" | "diffusion"  (see 00_utils.R)
KDE_BW_SCALE  <- 1.0          # KDE bandwidth multiplier (>1 = smoother field)
KDE_GRID_N    <- 200          # KDE grid resolution
KNN_K         <- 30           # k for "knn" / "diffusion"
DIFFUSION_VARIANT <- "heat"   # "implicit" | "heat" | "ppr" (when SCORE_METHOD="diffusion")
DIFFUSION_HEAT_STEPS <- 20    # backward-Euler steps for the "heat" variant
DIFFUSION_HEAT_ALPHA <- 2    # per-step timestep alpha (total diffusion time = steps*alpha)
# Local-tumor-fraction score (reviewer point #6): compute (K x)/(K 1), a
# kernel-weighted local tumor fraction that is natively in [0,1] and comparable
# ACROSS samples, then SKIP the per-sample min-max normalization (which otherwise
# stretches a tumor-poor mouse's ceiling to 1). TRUE is the corrected default.
DIFFUSION_LOCAL_FRACTION <- TRUE

TARGET_CELLTYPE <- "Tumor"

# --- Binning / genes ---------------------------------------------------------
N_BINS            <- 5
MIN_CELLS_PER_BIN <- 5
# --- Adaptive binning (optional) ---------------------------------------------
# When ADAPTIVE_N_BINS is TRUE, N_BINS is chosen PER CELL TYPE instead of using
# the fixed value above: candidate bin counts in [ADAPTIVE_N_BINS_MIN,
# ADAPTIVE_N_BINS_MAX] are scanned and the LARGEST n is picked for which at least
# (n - ADAPTIVE_BIN_TOL_K) bins are "valid" — i.e. at most k candidate bins fail.
# A bin is VALID if it would survive the common-support rule: it has
# >= MIN_SAMPLES_PER_COND_PER_BIN distinct samples in BOTH conditions, each
# contributing >= MIN_CELLS_PER_BIN cells (exactly the bins that feed the
# pseudobulk GAM). This keeps the binning as fine as the data support while
# tolerating k sparse tail bins (default k=2). The fixed N_BINS above is used when
# ADAPTIVE_N_BINS is FALSE (and remains the report's nominal value). The chosen n
# and the full scan are logged per cell type to <celltype>/adaptive_nbins_scan.csv.
ADAPTIVE_N_BINS     <- FALSE
ADAPTIVE_N_BINS_MIN <- 5
ADAPTIVE_N_BINS_MAX <- 20
ADAPTIVE_BIN_TOL_K  <- 0
# --- Common support (reviewer follow-up #1) ----------------------------------
# The fs difference smooth must be estimated from DATA in both conditions, not
# extrapolated by shrinkage where one condition has no cells. Because DKO tumor
# burden is ~0.47x WT and the tumor-poor DKO samples top out at low score, the
# upper score bins are WT-dominated / WT-only. COMMON_SUPPORT keeps only bins where
# BOTH conditions contribute >= MIN_SAMPLES_PER_COND_PER_BIN distinct samples (each
# already >= MIN_CELLS_PER_BIN cells), truncating the trajectory to the score range
# with genuine 2-vs-2+ replicated support. The retained score range is reported.
COMMON_SUPPORT               <- TRUE
MIN_SAMPLES_PER_COND_PER_BIN <- 2
# Pseudobulk GAM x-axis (reviewer #2). "score" = REAL mean local tumor fraction of
# each group -> fitted curves in interpretable score units, keeps within-bin
# variation (corrected default). "decile" = the previous evenly-spaced bin-index
# coordinate ((bin-0.5)/N_BINS) -> equal plot width per decile, but throws away the
# comparable score. Both are monotone maps common to both conditions, so neither
# biases the condition test; the choice is interpretability vs even spacing.
PB_SCORE_AXIS                <- "score"   # "score" | "decile"
GENE_SELECTION    <- "binned_variance"  # "binned_variance" | "spearman" | "svg"
TOP_N_GENES       <- 80
# --- Test set vs plot set (reviewer point #5) --------------------------------
# The GAM condition test runs on a BASELINE-EXPRESSION filtered gene set, using the
# SAME gene-level filters as the DESeq2 DE scripts (00_parameters.R): a gene is
# tested if it is detected in >= MIN_GENE_EXPR_FRAC of this cell type's cells AND has
# >= MIN_COUNTS_PER_SAMPLE counts in >= MIN_SAMPLES_EXPRESSED samples. This filter is
# independent of the condition x score interaction under the null, so BH over the
# whole passing set uses the correct m (not confounded with the test). The
# score-dependent GENE_SELECTION (binned_variance top TOP_N_GENES) is used ONLY to
# pick genes shown in the reference/wave heatmaps, never to define the tested set.
# SVG_GENES: placeholder for a future top-200 spatially-variable-gene set.
# When GENE_SELECTION == "svg", genes are taken from this vector (intersected
# with the cell type's expressed genes). Populate later from an SVG analysis.
SVG_GENES         <- NULL

# --- GAM models --------------------------------------------------------------
# Any subset of:
#   "pseudobulk"    — sample x bin weighted GAM, random INTERCEPT s(sample,bs="re").
#   "pseudobulk_fs" — same pseudobulk data, but a factor-SMOOTH random effect
#                     s(score,sample,bs="fs") that gives each sample its own
#                     penalized curve (random intercept AND shape). This is the
#                     statistically correct spec for the p_shape (trajectory-shape)
#                     test: between-mouse shape heterogeneity goes into the error
#                     term instead of being attributed to condition. The fs term
#                     subsumes the random intercept, so no separate re term. Cheap
#                     on the ~140 pseudobulk points. Expect FAR fewer significant
#                     genes than "pseudobulk"/"ranef" — that drop is the correct
#                     answer given the true n = 3 WT + 4 DKO, not a failure.
#   "pseudobulk_nb" — pseudobulk_rs random structure (random INTERCEPT + random
#                     LINEAR SLOPE, select=TRUE + gamma) but fit as a proper
#                     NEGATIVE-BINOMIAL GAM on SUMMED RAW COUNTS with a
#                     log(library-size) OFFSET (reviewer point #3-B). Removes the
#                     mean-of-log-normalized bias of the Gaussian pseudobulk specs:
#                     counts are summed per sample x bin, family=nb() models the
#                     count mean-variance relation, and the offset absorbs
#                     sequencing depth (incl. the systematically higher DKO library
#                     size). Same ordered-factor difference-smooth condition test
#                     (p_shape/p_level, now on the log link). NO n_cells weights —
#                     the offset carries depth (reviewer #3-D). This is the
#                     pseudobulk analog of tradeSeq's cell-level NB-GAM (15_07), but
#                     at the sample x bin level with a between-mouse random structure.
#   "pseudobulk_nb_fs" — COUNT-BASED version of the defensible shape spec: NB-GAM on
#                     summed raw counts + log(library-size) offset (like
#                     pseudobulk_nb) BUT with the per-sample factor SMOOTH random
#                     effect s(score,sample,bs="fs",k=FS_K) of pseudobulk_fs instead
#                     of the rs random slope. Pairs the correct count likelihood with
#                     the correct null for p_shape — between-mouse curve-SHAPE
#                     heterogeneity goes into the error term, not condition. fs
#                     subsumes the per-sample intercept (no re term); FS_K kept small
#                     for order-stability; no n_cells weights (offset carries depth).
#                     The most defensible single spec when a count likelihood is
#                     wanted for the shape test — still bounded by the n=3v4
#                     permutation floor, so report as a ranking, not per-gene FDR.
#   "pseudobulk_rs" — same pseudobulk data, but the per-sample fs CURVE is
#                     replaced by a random INTERCEPT s(sample,bs="re") + random
#                     LINEAR SLOPE s(score,sample,bs="re") — two variance
#                     components instead of a penalized curve family. Far better
#                     conditioned at 7 samples / ~9 points per mouse (no basis
#                     saturation), so it needs NO FS_K clamp and is order-stable
#                     by construction (reviewer: fs is over-parameterized here).
#                     Still puts between-mouse LEVEL + linear-TREND heterogeneity
#                     in the error term; only NON-linear per-mouse shape (which is
#                     unidentifiable at 7 mice anyway) is given up. Fit with
#                     select=TRUE + gamma>1 (see RS_SELECT / RS_GAMMA) so the
#                     difference smooth can shrink to exactly zero. This is the
#                     more standard, more defensible shape spec at this n.
#   "ranef"         — single-cell GAM, random INTERCEPT s(sample,bs="re").
#   "fs"            — single-cell factor-smooth (slow; the pseudobulk_fs variant
#                     is preferred for the shape test).
# Each is fit + saved separately. "pseudobulk" MUST be listed first when combined
# with non-FULL_GENE_MODELS, because its p_shape ranking picks genes for them.
# NOTE on interpretation: pseudobulk/ranef (random intercept only) are
# ANTICONSERVATIVE for p_shape (the by-condition smooth is tested against bin-level
# residuals, not the 7 mice) — treat them as exploratory rankings. pseudobulk_fs is
# the defensible shape test; p_level (level offset) is the trustworthy
# between-mouse level test in all specs.
GAM_MODELS <- c("pseudobulk", "pseudobulk_fs", "pseudobulk_rs", "pseudobulk_nb", "pseudobulk_nb_fs")
# Specs fitted on the sample x bin PSEUDOBULK data (vs the single-cell "fs"/"ranef").
# Used throughout to decide pseudobulk-data build + axis labels.
PSEUDOBULK_MODELS <- c("pseudobulk", "pseudobulk_fs", "pseudobulk_rs", "pseudobulk_nb", "pseudobulk_nb_fs")
# Subset of PSEUDOBULK_MODELS whose pseudobulk is built by SUMMING RAW COUNTS (not
# averaging the log-normalized data layer) and fit with family=nb() + a
# log(library-size) offset. These need the counts matrix + a per-group library size,
# so they build a different pseudobulk object (agg="sum") than the Gaussian specs.
PSEUDOBULK_NB_MODELS <- c("pseudobulk_nb", "pseudobulk_nb_fs")
GAM_K      <- 5
# Basis dimension for the per-sample factor-smooth s(score, sample, bs="fs"). MUST be
# kept small: with only ~9 pseudobulk points per sample, a high-k fs smooth
# near-saturates the model, the REML smoothing-parameter surface becomes multimodal,
# and the fitted p_shape becomes ROW-ORDER DEPENDENT (see notes.md, 2026-09-10). k=3
# is well-determined here and makes the fit exactly order-invariant. The trajectory
# smooths s(score)/by-smooth stay at GAM_K.
FS_K       <- 3
# Penalty ORDER of the per-sample factor-smooth s(score, sample, bs="fs", m=FS_M).
# The other lever (besides FS_K) for controlling per-sample-curve complexity, and the
# reviewer's suggested alternative to clamping FS_K: m sets which derivative the
# penalty acts on. m=2 (default) penalizes the 2nd derivative (curvature) -> pulls each
# per-sample curve toward a straight line (unpenalized null space {1, x}). m=1
# penalizes the 1st derivative (slope) -> pulls toward a flat constant (null space {1}),
# shrinking the per-mouse effective parameters harder. Both notes.md-tested stable at
# n=3v4: the current default is FS_K=3, m=2 (fewer basis fns, curvature penalty); the
# reviewer's alternative FS_K=5, FS_M=1 keeps fine per-sample resolution with a stiffer
# lower-order penalty. Stable configs give DIFFERENT hit lists — a legitimate modelling
# choice, not instability (notes.md 2026-09-10). Applies to every fs term:
# pseudobulk_fs, pseudobulk_nb_fs, and single-cell fs. For a 1-D marginal smooth m>=1.
FS_M       <- 2
# --- pseudobulk_rs (random intercept + random slope) fit settings ------------
# rs uses NO FS_K (it has no per-sample basis to saturate). select=TRUE adds a
# shrinkage penalty on each smooth's null space so the ordered-factor DIFFERENCE
# smooth can shrink to exactly zero (removes the unpenalized-linear boundary that
# makes the by-smooth p-value anticonservative; Marra & Wood 2011). gamma>1
# multiplies the edf cost in REML, biasing toward smoother fits (Wood ~1.4).
RS_SELECT <- TRUE
RS_GAMMA  <- 1.4
# Models fit on the FULL test_genes set (all expression-filtered genes). Any GAM
# model in GAM_MODELS but NOT here is restricted to the top N_GENES_OTHER_MODELS
# pseudobulk genes (cheap cross-check).
FULL_GENE_MODELS     <- c("pseudobulk", "ranef", "pseudobulk_fs", "pseudobulk_rs", "pseudobulk_nb", "pseudobulk_nb_fs")
N_GENES_OTHER_MODELS <- 10   # gene budget for non-full models (e.g. "fs")

# --- Condition-test visualization -------------------------------------------
# Per-gene GAM scatter+spline PDFs. A vector of model names for which to emit the
# per-gene gam_<gene>.pdf splines; character(0) disables them entirely. To keep the
# folder readable, PDFs are emitted ONLY for genes passing BH-FDR < SIG_FDR (shape
# or level) for that model — the same set shown in the significant-trajectory
# heatmap — not for every tested gene. Must name models that are in GAM_MODELS.
# Default: every pseudobulk spec (PSEUDOBULK_MODELS), so each gam_<model>/ folder
# gets its own per-gene splines for side-by-side comparison. Set to c("pseudobulk_fs")
# for just the primary spec, or character(0) to disable entirely. Each gam_<model>/
# also gets a per_sample_curves/ subfolder with one gam_<gene>_per_sample.pdf per
# significant gene: the WT/DKO population curves overlaid with EACH sample's own
# fitted trajectory (random effect included), showing the between-mouse spread.
PLOT_INDIVIDUAL_GAMS <- PSEUDOBULK_MODELS
# BH-FDR cutoff for a gene to enter the significant-trajectory heatmap. A gene is
# kept if fdr_shape < SIG_FDR OR fdr_level < SIG_FDR. Set to 0.1 (was 0.05) to
# admit borderline trajectory-shape genes (e.g. Klrk1, Cxcr6) for visualization;
# the per-gene FDR bars in the heatmap still show exact significance.
SIG_FDR              <- 0.1
HEATMAP_ROW_ORDER    <- "cluster"  # "cluster" (hclust on fitted curves) | "shape" (by p_shape)

# Restricted to CD8+ T cells for now (was IMMUNE_CELLTYPES_DE).
CELLTYPES_TRAJECTORY <- c("CD8+ T")

# --- Quality control ---------------------------------------------------------
# Drop low-quality cells from the *expression trajectory* (same thresholds as
# 08_01 / 15_05: total_counts >= MIN_TOTAL_COUNTS & n_genes >= MIN_GENES_DETECTED,
# from 00_parameters.R). Sparse cells add noise to the per-gene GAM fits and can
# distort the decile x-axis. Default TRUE; set FALSE to use all cells.
# NOTE: QC is applied only when selecting cells for the trajectory, NOT when
# computing the tumor-proximity score — the tumor field is a property of tissue
# architecture, so every physically present cell (QC-pass or not) should shape it.
# Default FALSE: reviewed and the QC filter is intentionally OFF for this analysis.
APPLY_QC <- TRUE

stopifnot(SCORE_METHOD %in% c("kde", "knn", "diffusion"))
stopifnot(all(GAM_MODELS %in% c("fs", "ranef", "pseudobulk", "pseudobulk_fs", "pseudobulk_rs", "pseudobulk_nb", "pseudobulk_nb_fs")))
stopifnot(GENE_SELECTION %in% c("binned_variance", "spearman", "svg"))
stopifnot(HEATMAP_ROW_ORDER %in% c("cluster", "shape"))
stopifnot(PB_SCORE_AXIS %in% c("score", "decile"))
# fs penalty order: 1-D marginal smooth requires 2*m > d = 1, so m >= 1; k must
# exceed the null-space dimension (m + 1) for the penalized part to exist.
stopifnot(FS_M >= 1, FS_K > FS_M)
stopifnot(is.logical(ADAPTIVE_N_BINS), ADAPTIVE_N_BINS_MIN >= 2,
          ADAPTIVE_N_BINS_MAX >= ADAPTIVE_N_BINS_MIN, ADAPTIVE_BIN_TOL_K >= 0)

# =============================================================================
# Parameter report (written FIRST, before any heavy computation)
# =============================================================================
# Snapshot every run parameter into OUTPUT_DIR so each result folder is
# self-documenting and reproducible. Emitted before the readRDS / scoring / GAM
# chunks so that even a crashed or interrupted run leaves a record of exactly
# which settings produced (or attempted) the outputs. Writes both a machine-
# readable run_parameters.csv (group, parameter, value) and a human-readable
# run_parameters.txt.
write_param_report <- function(dir) {
  fmt  <- function(x) paste(as.character(x), collapse = ", ")
  grp  <- function(group, ...) {
    kv <- list(...)
    data.frame(group = group, parameter = names(kv),
               value = vapply(kv, fmt, character(1)),
               row.names = NULL, stringsAsFactors = FALSE)
  }
  pkgver <- function(p) if (requireNamespace(p, quietly = TRUE))
                          as.character(utils::packageVersion(p)) else "NA"
  rep_df <- rbind(
    grp("run",
        script     = "code/R/15_04_tumor_trajectory_working.R",
        run_time   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        host       = as.character(Sys.info()[["nodename"]]),
        R_version  = R.version.string,
        mgcv       = pkgver("mgcv"),
        Seurat     = pkgver("Seurat"),
        ComplexHeatmap = pkgver("ComplexHeatmap"),
        N_THREADS  = N_THREADS,
        input_rds  = "output/R/DKO_banksy_merged.rds"),
    grp("tumor_score",
        SCORE_METHOD = SCORE_METHOD, DIFFUSION_VARIANT = DIFFUSION_VARIANT,
        DIFFUSION_HEAT_STEPS = DIFFUSION_HEAT_STEPS,
        DIFFUSION_HEAT_ALPHA = DIFFUSION_HEAT_ALPHA,
        DIFFUSION_LOCAL_FRACTION = DIFFUSION_LOCAL_FRACTION,
        KNN_K = KNN_K, KDE_BW_SCALE = KDE_BW_SCALE, KDE_GRID_N = KDE_GRID_N,
        TARGET_CELLTYPE = TARGET_CELLTYPE),
    grp("binning_support",
        N_BINS = if (ADAPTIVE_N_BINS) sprintf("%d (nominal; adaptive per cell type)", N_BINS)
                 else N_BINS,
        MIN_CELLS_PER_BIN = MIN_CELLS_PER_BIN,
        ADAPTIVE_N_BINS = ADAPTIVE_N_BINS,
        ADAPTIVE_N_BINS_MIN = ADAPTIVE_N_BINS_MIN,
        ADAPTIVE_N_BINS_MAX = ADAPTIVE_N_BINS_MAX,
        ADAPTIVE_BIN_TOL_K = ADAPTIVE_BIN_TOL_K,
        COMMON_SUPPORT = COMMON_SUPPORT,
        MIN_SAMPLES_PER_COND_PER_BIN = MIN_SAMPLES_PER_COND_PER_BIN,
        PB_SCORE_AXIS = PB_SCORE_AXIS),
    grp("genes",
        GENE_SELECTION = GENE_SELECTION, TOP_N_GENES = TOP_N_GENES,
        MIN_GENE_EXPR_FRAC = MIN_GENE_EXPR_FRAC,
        MIN_COUNTS_PER_SAMPLE = MIN_COUNTS_PER_SAMPLE,
        MIN_SAMPLES_EXPRESSED = MIN_SAMPLES_EXPRESSED),
    grp("gam_models",
        GAM_MODELS = GAM_MODELS, GAM_K = GAM_K, FS_K = FS_K, FS_M = FS_M,
        RS_SELECT = RS_SELECT, RS_GAMMA = RS_GAMMA,
        FULL_GENE_MODELS = FULL_GENE_MODELS,
        N_GENES_OTHER_MODELS = N_GENES_OTHER_MODELS),
    grp("condition_test_viz",
        PLOT_INDIVIDUAL_GAMS = if (length(PLOT_INDIVIDUAL_GAMS))
                                 PLOT_INDIVIDUAL_GAMS else "(none)",
        SIG_FDR = SIG_FDR, HEATMAP_ROW_ORDER = HEATMAP_ROW_ORDER),
    grp("qc",
        APPLY_QC = APPLY_QC, MIN_TOTAL_COUNTS = MIN_TOTAL_COUNTS,
        MIN_GENES_DETECTED = MIN_GENES_DETECTED),
    grp("design",
        REFERENCE_CONDITION = REFERENCE_CONDITION,
        TEST_CONDITION = TEST_CONDITION, WT_SAMPLES = WT_SAMPLES,
        CELLTYPES_TRAJECTORY = CELLTYPES_TRAJECTORY)
  )
  write.csv(rep_df, file.path(dir, "run_parameters.csv"), row.names = FALSE)
  con <- file(file.path(dir, "run_parameters.txt"), open = "wt")
  writeLines(c(strrep("=", 78),
               "15_04 tumor-trajectory working run parameters",
               strrep("=", 78)), con)
  for (g in unique(rep_df$group)) {
    writeLines(sprintf("\n-- %s --", g), con)
    sub <- rep_df[rep_df$group == g, ]
    writeLines(sprintf("  %-30s = %s", sub$parameter, sub$value), con)
  }
  writeLines(c("", strrep("=", 78)), con)
  close(con)
  invisible(rep_df)
}
param_report <- write_param_report(OUTPUT_DIR)
cat(sprintf("Wrote parameter report -> %s\n",
            file.path(OUTPUT_DIR, "run_parameters.{txt,csv}")))
print(param_report, row.names = FALSE)

# =============================================================================
# Load data
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

DefaultAssay(seu) <- "RNA"

# --- QC mask (per cell) ------------------------------------------------------
# Logical vector over ALL cells. Applied later at per-cell-type selection so it
# filters the trajectory expression set but not the tumor-score field. When
# APPLY_QC is FALSE, every cell passes.
if (APPLY_QC) {
  qc_pass <- seu$total_counts >= MIN_TOTAL_COUNTS &
             seu$n_genes_by_counts >= MIN_GENES_DETECTED
  cat(sprintf("QC (total_counts >= %d, n_genes >= %d): %d / %d cells pass (%.1f%%)\n",
              MIN_TOTAL_COUNTS, MIN_GENES_DETECTED, sum(qc_pass), length(qc_pass),
              100 * mean(qc_pass)))
} else {
  qc_pass <- rep(TRUE, ncol(seu))
  cat("QC filter disabled (APPLY_QC = FALSE): using all cells\n")
}

# =============================================================================
# Compute tumor score (per sample, normalized to [0,1]) via 00_utils.R
# =============================================================================
all_coords <- Embeddings(seu, "spatial")
celltypes  <- seu$Tier1_celltype
samples    <- seu$sample
is_tumor   <- celltypes == TARGET_CELLTYPE

cat(sprintf("Computing tumor score (method = %s)...\n", SCORE_METHOD))
tumor_score <- switch(SCORE_METHOD,
  kde = compute_tumor_score(all_coords, is_tumor, samples, method = "kde",
                            bw_scale = KDE_BW_SCALE, grid_n = KDE_GRID_N),
  knn = compute_tumor_score(all_coords, is_tumor, samples, method = "knn",
                            k = KNN_K),
  diffusion = compute_tumor_score(all_coords, is_tumor, samples, method = "diffusion",
                                  diffusion_method = DIFFUSION_VARIANT, k = KNN_K,
                                  heat_steps = DIFFUSION_HEAT_STEPS,
                                  local_fraction = DIFFUSION_LOCAL_FRACTION,
                                  # Local-fraction score is already comparable [0,1];
                                  # skip per-sample min-max (reviewer point #6).
                                  normalize = !DIFFUSION_LOCAL_FRACTION,
                                  alpha = DIFFUSION_HEAT_ALPHA)
)
seu$tumor_score <- tumor_score
cat("Tumor score computed.\n")

# X-axis labels. The raw score (SCORE_LAB) is used where the actual (skewed)
# distribution matters (histogram, cell-density plot). The trajectory plots (GAM
# curves, normalized-expression heatmap) instead use an evenly-spaced *decile*
# coordinate (TRAJ_X_LAB) so the heavily zero-concentrated diffusion score is not
# squished against the left edge — equal plot width per decile of cells.
SCORE_LAB   <- sprintf("Tumor Proximity Score (%s)", toupper(SCORE_METHOD))
TRAJ_X_LAB  <- sprintf("Tumor proximity (decile of %s score)", toupper(SCORE_METHOD))
# X-label for the pseudobulk GAM curves, matching PB_SCORE_AXIS.
GAM_X_LAB   <- if (PB_SCORE_AXIS == "score") {
                 sprintf("Tumor proximity (%s local tumor fraction)", toupper(SCORE_METHOD))
               } else TRAJ_X_LAB

# =============================================================================
# Helpers
# =============================================================================

# --- Binning ----------------------------------------------------------------
# Decile bins of a density/score vector; returns bin assignment + midpoints.
make_bins <- function(score_vec, n_bins = N_BINS) {
  bin_breaks <- quantile(score_vec, probs = seq(0, 1, length.out = n_bins + 1),
                         na.rm = TRUE)
  bin_breaks[1] <- bin_breaks[1] - 1e-6
  if (length(unique(bin_breaks)) < length(bin_breaks)) {
    bin_breaks <- seq(min(score_vec, na.rm = TRUE) - 1e-6,
                      max(score_vec, na.rm = TRUE),
                      length.out = n_bins + 1)
  }
  bins <- cut(score_vec, breaks = bin_breaks, labels = seq_len(n_bins))
  midpoints <- (bin_breaks[-1] + bin_breaks[-(n_bins + 1)]) / 2
  list(bins = bins, breaks = bin_breaks, midpoints = midpoints)
}

# --- Adaptive N_BINS selection ----------------------------------------------
# For a candidate bin count, how many bins would be RETAINED under the common-
# support rule: >= min_samp_per_cond distinct samples in BOTH conditions, each
# with >= min_cells cells. These are exactly the bins that feed the pseudobulk GAM
# (see make_pseudobulk()'s common_support filter), so this is the meaningful
# "valid bin" count whether or not COMMON_SUPPORT is toggled on for the fit.
count_valid_bins <- function(score_vec, sample_vec, condition_vec, n_bins,
                             min_cells = MIN_CELLS_PER_BIN,
                             min_samp_per_cond = MIN_SAMPLES_PER_COND_PER_BIN) {
  b   <- make_bins(score_vec, n_bins)
  df  <- data.frame(bin = as.integer(b$bins),
                    sample = as.character(sample_vec),
                    condition = as.character(condition_vec))
  df  <- df[!is.na(df$bin), ]
  grp <- df %>% dplyr::count(bin, sample, condition, name = "n_cells") %>%
                dplyr::filter(n_cells >= min_cells)
  if (nrow(grp) == 0) return(0L)
  per_bin <- grp %>% dplyr::group_by(bin) %>%
    dplyr::summarise(
      n_wt  = dplyr::n_distinct(sample[condition == REFERENCE_CONDITION]),
      n_dko = dplyr::n_distinct(sample[condition == TEST_CONDITION]),
      .groups = "drop")
  as.integer(sum(per_bin$n_wt >= min_samp_per_cond &
                 per_bin$n_dko >= min_samp_per_cond))
}

# Pick the LARGEST candidate n in [n_min, n_max] with #valid >= n - k. Scans from
# n_max downward so the first satisfying n is the maximum; falls back to n_min
# (with a warning) if none qualify. The full scan table is attached as attr "scan".
choose_adaptive_nbins <- function(score_vec, sample_vec, condition_vec,
                                  k = ADAPTIVE_BIN_TOL_K,
                                  n_min = ADAPTIVE_N_BINS_MIN,
                                  n_max = ADAPTIVE_N_BINS_MAX) {
  cand <- seq.int(n_max, n_min)
  vals <- vapply(cand, function(n)
    count_valid_bins(score_vec, sample_vec, condition_vec, n), integer(1))
  pass <- vals >= (cand - k)
  chosen <- as.integer(if (any(pass)) cand[which(pass)[1]] else n_min)
  if (!any(pass))
    warning(sprintf("adaptive N_BINS: no n in [%d,%d] has >= n-%d valid bins; using %d",
                    n_min, n_max, k, n_min))
  attr(chosen, "scan") <- data.frame(n_bins = cand, n_valid = vals,
                                     required = cand - k, pass = pass)
  chosen
}

# Mean expression per bin (genes x bins).
binned_expression <- function(expr_mat, bins, n_bins = N_BINS) {
  out <- matrix(NA, nrow = nrow(expr_mat), ncol = n_bins,
                dimnames = list(rownames(expr_mat), paste0("bin_", seq_len(n_bins))))
  for (b in seq_len(n_bins)) {
    cells <- which(as.integer(bins) == b)
    if (length(cells) > 1) {
      out[, b] <- rowMeans(expr_mat[, cells, drop = FALSE])
    } else if (length(cells) == 1) {
      out[, b] <- expr_mat[, cells]
    }
  }
  out
}

# --- Gene selection ---------------------------------------------------------
select_genes <- function(expr_mat, score_vec, method = GENE_SELECTION,
                         top_n = TOP_N_GENES) {
  if (method == "svg") {
    if (is.null(SVG_GENES))
      stop("GENE_SELECTION='svg' but SVG_GENES is NULL - populate SVG_GENES first.")
    g <- intersect(SVG_GENES, rownames(expr_mat))
    gene_means <- rowMeans(expr_mat[g, , drop = FALSE])
    return(g[gene_means > 0.05][seq_len(min(top_n, sum(gene_means > 0.05)))])
  }

  if (method == "spearman") {
    rho <- apply(expr_mat, 1, function(g) if (sd(g) == 0) 0 else
                 cor(g, score_vec, method = "spearman"))
    rho <- rho[!is.na(rho)]
    g <- names(sort(abs(rho), decreasing = TRUE))
    gene_means <- rowMeans(expr_mat[g, , drop = FALSE])
    g <- g[gene_means > 0.05]
    return(g[seq_len(min(top_n, length(g)))])
  }

  # binned_variance (default): variance of bin-mean expression
  bins <- make_bins(score_vec)$bins
  be <- binned_expression(expr_mat, bins)
  row_vars  <- apply(be, 1, var, na.rm = TRUE)
  row_means <- rowMeans(be, na.rm = TRUE)
  keep <- row_means > 0.05 & !is.na(row_vars)
  g <- names(sort(row_vars[keep], decreasing = TRUE))
  g[seq_len(min(top_n, length(g)))]
}

# --- Pseudobulk (sample x bin) ----------------------------------------------
make_pseudobulk <- function(expr_mat, score_vec, sample_vec, condition_vec,
                            n_bins = N_BINS, min_cells = MIN_CELLS_PER_BIN,
                            common_support = COMMON_SUPPORT,
                            min_samp_per_cond = MIN_SAMPLES_PER_COND_PER_BIN,
                            agg = c("mean", "sum")) {
  # agg = "mean" -> per-gene mean of the log-normalized data layer (Gaussian specs);
  #       "sum"  -> per-gene SUM of raw counts (NB-GAM specs), and libsize (group
  #                 total counts) is recorded for the log-offset.
  agg <- match.arg(agg)
  b <- make_bins(score_vec, n_bins)
  bins <- as.integer(b$bins)
  groups <- paste(sample_vec, bins, sep = "__")

  pb_list <- list(); pb_meta <- list()
  for (g in unique(groups)) {
    cell_idx <- which(groups == g)
    if (length(cell_idx) < min_cells) next
    parts <- strsplit(g, "__")[[1]]
    bin_id <- as.integer(parts[2])
    if (is.na(bin_id)) next
    pb_list[[g]] <- if (agg == "sum") rowSums(expr_mat[, cell_idx, drop = FALSE])
                    else             rowMeans(expr_mat[, cell_idx, drop = FALSE])
    pb_meta[[g]] <- data.frame(
      group = g, sample = parts[1],
      condition = as.character(condition_vec[cell_idx[1]]),
      bin = bin_id,
      # REAL score of the group = mean local tumor fraction of its cells (reviewer
      # #2), NOT the bin index — so fitted curves are in interpretable score units
      # and within-bin score variation is retained.
      score = mean(score_vec[cell_idx]),
      # NB log-offset source: group total counts across ALL input rows. NA for the
      # mean-of-log build (never read there).
      libsize = if (agg == "sum") sum(expr_mat[, cell_idx]) else NA_real_,
      n_cells = length(cell_idx), stringsAsFactors = FALSE)
  }
  if (length(pb_list) == 0) return(NULL)
  meta <- bind_rows(pb_meta)
  meta$condition <- factor(meta$condition,
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION))

  # --- Common-support filter (reviewer #1) -----------------------------------
  # Keep only bins where BOTH conditions contribute >= min_samp_per_cond distinct
  # samples (each group already has >= min_cells cells). This confines the fs
  # difference smooth to the score range with real 2-vs-2+ data in both conditions
  # instead of letting it extrapolate over WT-only upper bins.
  support <- NULL
  if (common_support) {
    per_bin <- meta %>%
      group_by(bin) %>%
      summarise(n_wt  = dplyr::n_distinct(sample[condition == REFERENCE_CONDITION]),
                n_dko = dplyr::n_distinct(sample[condition == TEST_CONDITION]),
                .groups = "drop")
    keep_bins <- sort(per_bin$bin[per_bin$n_wt  >= min_samp_per_cond &
                                  per_bin$n_dko >= min_samp_per_cond])
    support <- list(per_bin = per_bin, keep_bins = keep_bins)
    meta <- meta[meta$bin %in% keep_bins, ]
    if (nrow(meta) < 4) return(NULL)     # too little common support to fit
    pb_list <- pb_list[meta$group]       # realign expr columns to filtered meta
  }

  # Canonical row order (sample, then bin) so the fit is deterministic and matches
  # any other script building the same pseudobulk (belt-and-suspenders on top of the
  # order-stable fs_k; see notes.md 2026-09-10).
  ord <- order(meta$sample, meta$bin)
  meta <- meta[ord, ]
  pb_list <- pb_list[meta$group]

  list(expr = do.call(cbind, pb_list), meta = meta, support = support)
}

# --- GAM fitting (model-specific) -------------------------------------------
# Returns list(pred, p_shape, p_level, ok). `pred` is the population-level curve
# per condition (random effects excluded).
#
# Condition test uses the mgcv-canonical *ordered-factor difference smooth*
# (Wood 2017, "difference smooths") rather than a two-model deviance LRT. With
# condition coded as an ordered factor (WT = reference), the model is
#   expr ~ s(score) + s(score, by = condition_ord) + condition_ord + <rand eff>
# where s(score) is the WT trajectory, s(score, by = condition_ord) is the DKO
# *deviation* smooth, and condition_ord is a constant DKO offset. This cleanly
# separates two questions and, unlike anova(m0, m1) on penalized smooths, never
# returns NA:
#   * p_shape = does the trajectory SHAPE differ (the deviation smooth != 0) —
#               the trajectory question; the novel result.
#   * p_level = does the OVERALL level differ (parametric offset) — already
#               captured by standard bulk DKO-vs-WT DE.
fit_trajectory_gam <- function(df, model, gam_k = GAM_K, fs_k = FS_K, fs_m = FS_M) {
  df$condition <- droplevels(factor(df$condition,
                                    levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
  df$sample <- droplevels(factor(df$sample))
  if (nlevels(df$condition) < 2) return(list(ok = FALSE))

  # Ordered factor (WT < DKO) -> by-smooth is the DKO-vs-WT *difference* smooth.
  df$cond_ord <- ordered(as.character(df$condition),
                         levels = c(REFERENCE_CONDITION, TEST_CONDITION))

  # Model-specific formula (reference + difference smooth + offset + rand eff)
  # and the random-effect term to exclude when predicting population curves.
  # sel / gam_gamma default to gam()'s defaults (no double penalty, gamma 1); fam
  # defaults to gaussian(). pseudobulk_rs overrides sel/gam_gamma; pseudobulk_nb
  # overrides all three (NB family) and adds a log-library-size offset.
  sel <- FALSE; gam_gamma <- 1; fam <- gaussian()
  if (model == "pseudobulk") {
    if (nlevels(df$sample) < 2) return(list(ok = FALSE))
    # Sample random intercept so the ~n pseudobulk points (from only 7 samples)
    # are not treated as independent replicates -> guards against residual
    # pseudoreplication at the sample level (LEVEL only, not shape).
    form <- expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
                   cond_ord + s(sample, bs = "re")
    excl <- "s(sample)"
    wcol <- df$weight
  } else if (model == "pseudobulk_fs") {
    if (nlevels(df$sample) < 2) return(list(ok = FALSE))
    # Pseudobulk data + factor-SMOOTH random effect: each sample gets its own
    # penalized curve, so between-mouse trajectory-SHAPE heterogeneity is in the
    # error term. This is the correct null for the p_shape test. The fs term
    # subsumes the per-sample intercept, so no separate s(sample, bs="re").
    # Weighted by n_cells like the standard pseudobulk fit. fs term uses the smaller
    # fs_k (not gam_k) to keep the per-sample smooth well-determined -> order-stable.
    form <- expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
                   cond_ord + s(score, sample, bs = "fs", k = fs_k, m = fs_m)
    excl <- "s(score,sample)"
    wcol <- df$weight
  } else if (model == "pseudobulk_rs") {
    if (nlevels(df$sample) < 2) return(list(ok = FALSE))
    # Pseudobulk data + random INTERCEPT s(sample,bs="re") and random LINEAR SLOPE
    # s(score,sample,bs="re"): two variance components in place of the per-sample fs
    # curve. Well-conditioned at 7 samples (no basis saturation) -> order-stable with
    # NO FS_K clamp. Between-mouse level + linear-trend heterogeneity stays in the
    # error term. select=TRUE + gamma>1 (RS_SELECT / RS_GAMMA) make the by-smooth able
    # to shrink to zero and bias toward smoother fits. BOTH random terms are excluded
    # from the population-curve prediction.
    form <- expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
                   cond_ord + s(sample, bs = "re") + s(score, sample, bs = "re")
    excl <- c("s(sample)", "s(score,sample)")
    wcol <- df$weight
    sel  <- RS_SELECT; gam_gamma <- RS_GAMMA
  } else if (model == "pseudobulk_nb") {
    if (nlevels(df$sample) < 2) return(list(ok = FALSE))
    # NB-GAM on SUMMED RAW COUNTS + log(library-size) offset (reviewer #3-B), with the
    # pseudobulk_rs between-mouse random structure: random INTERCEPT s(sample,bs="re")
    # + random LINEAR SLOPE s(score,sample,bs="re"). family=nb() estimates the count
    # dispersion; the offset absorbs sequencing depth (incl. the DKO>WT library-size
    # gap), so NO n_cells weights (reviewer #3-D). select + gamma let the difference
    # smooth shrink to zero, as in rs. Both random terms + the offset are handled in
    # the population-curve prediction (offset held at its mean; see below).
    form <- expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
                   cond_ord + s(sample, bs = "re") + s(score, sample, bs = "re") +
                   offset(log_libsize)
    excl <- c("s(sample)", "s(score,sample)")
    wcol <- NULL
    sel  <- RS_SELECT; gam_gamma <- RS_GAMMA; fam <- nb()
  } else if (model == "pseudobulk_nb_fs") {
    if (nlevels(df$sample) < 2) return(list(ok = FALSE))
    # NB-GAM on SUMMED RAW COUNTS + log(library-size) offset (reviewer #3-B) with the
    # pseudobulk_fs per-sample factor SMOOTH random effect: each sample gets its own
    # penalized curve, so between-mouse trajectory-SHAPE heterogeneity is in the error
    # term (the correct null for p_shape) AND the count likelihood is correct. The fs
    # term subsumes the per-sample intercept, so no separate s(sample,bs="re"). fs_k
    # kept small (order-stable, see notes.md 2026-09-10). No n_cells weights (offset
    # carries depth). select/gamma left at gam() defaults, matching pseudobulk_fs.
    form <- expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
                   cond_ord + s(score, sample, bs = "fs", k = fs_k, m = fs_m) +
                   offset(log_libsize)
    excl <- "s(score,sample)"
    wcol <- NULL
    fam  <- nb()
  } else if (model == "fs") {
    if (nlevels(df$sample) < 2) return(list(ok = FALSE))
    form <- expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
                   cond_ord + s(score, sample, bs = "fs", k = fs_k, m = fs_m)
    excl <- "s(score,sample)"
    wcol <- NULL
  } else if (model == "ranef") {
    if (nlevels(df$sample) < 2) return(list(ok = FALSE))
    form <- expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
                   cond_ord + s(sample, bs = "re")
    excl <- "s(sample)"
    wcol <- NULL
  }

  # Single REML fit serves both the plotted curve and the significance tests.
  # select / gamma are gam() defaults (FALSE / 1) except for pseudobulk_rs.
  fit <- tryCatch(gam(form, data = df, weights = wcol, family = fam, method = "REML",
                      select = sel, gamma = gam_gamma),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE))

  sm <- summary(fit)
  # Difference smooth row: the only s.table term carrying the by-factor (":co").
  s_rows  <- rownames(sm$s.table)
  diff_row <- s_rows[grepl(":cond_ord", s_rows, fixed = TRUE)]
  p_shape  <- if (length(diff_row)) sm$s.table[diff_row[1], "p-value"] else NA_real_
  # Parametric offset row: ordered factor -> "cond_ord.L" (linear contrast).
  # Grep the p-value column rather than hard-coding "Pr(>|t|)": under a non-Gaussian
  # family (e.g. nb()) mgcv reports "Pr(>|z|)", so a fixed name would silently
  # return NA. The p-value is always the last p.table column matching "^Pr(".
  p_rows   <- rownames(sm$p.table)
  lev_row  <- p_rows[grepl("cond_ord", p_rows, fixed = TRUE)]
  pcol     <- grep("^Pr\\(", colnames(sm$p.table))
  p_level  <- if (length(lev_row) && length(pcol))
                sm$p.table[lev_row[1], pcol[1]] else NA_real_

  # mgcv can underflow an extreme smooth/parametric p-value to exactly 0. Floor at
  # machine precision so the CSV / BH-FDR report a representable value (~2.2e-16)
  # rather than a literal 0 (which is not a valid p-value).
  if (!is.na(p_shape)) p_shape <- max(p_shape, .Machine$double.eps)
  if (!is.na(p_level)) p_level <- max(p_level, .Machine$double.eps)

  # Population-level prediction grid (random effects excluded)
  newd <- expand.grid(
    score = seq(min(df$score), max(df$score), length.out = 100),
    cond_ord = ordered(c(REFERENCE_CONDITION, TEST_CONDITION),
                       levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
  newd$sample <- df$sample[1]   # dummy; excluded from prediction
  # NB model carries a log(library-size) offset; predict at the MEAN depth so the
  # population curve is on the log-E[count] scale at a representative library size
  # (a constant offset shift, irrelevant to the z-scored heatmap / difference test).
  if (model %in% PSEUDOBULK_NB_MODELS) newd$log_libsize <- mean(df$log_libsize)
  pred <- predict(fit, newdata = newd, se.fit = TRUE,
                  exclude = excl, newdata.guaranteed = TRUE)
  newd$condition <- factor(as.character(newd$cond_ord),
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  newd$fit   <- pred$fit
  newd$se    <- pred$se.fit
  newd$ci_lo <- pred$fit - 1.96 * pred$se.fit
  newd$ci_hi <- pred$fit + 1.96 * pred$se.fit

  # --- Per-sample fitted curves (random effects INCLUDED) --------------------
  # The population curves above EXCLUDE the sample random effect; here we predict
  # EACH sample's own trajectory WITH its random deviation (fs curve / re
  # intercept+slope). Predicted over that sample's observed score range only, so
  # the per-mouse smooth is interpolated, not extrapolated. Each sample's
  # condition is fixed (a mouse is either WT or DKO). NB specs predict at the mean
  # log(library-size) so the per-sample curves share the population curve's scale.
  # Used for the per-gene per-sample overlay plots; NA-safe if a fit fails.
  samp_levels <- levels(df$sample)
  pred_sample <- bind_rows(lapply(samp_levels, function(s) {
    di <- df[df$sample == s, , drop = FALSE]
    if (nrow(di) < 2) return(NULL)
    grid <- seq(min(di$score), max(di$score), length.out = 60)
    nd <- data.frame(
      score    = grid,
      cond_ord = ordered(as.character(di$cond_ord[1]),
                         levels = c(REFERENCE_CONDITION, TEST_CONDITION)),
      sample   = factor(s, levels = samp_levels))
    if (model %in% PSEUDOBULK_NB_MODELS) nd$log_libsize <- mean(df$log_libsize)
    pr <- tryCatch(as.numeric(predict(fit, newdata = nd, newdata.guaranteed = TRUE)),
                   error = function(e) rep(NA_real_, nrow(nd)))
    data.frame(score = grid, sample = s,
               condition = factor(as.character(di$cond_ord[1]),
                                  levels = c(REFERENCE_CONDITION, TEST_CONDITION)),
               fit = pr, stringsAsFactors = FALSE)
  }))

  list(pred = newd, pred_sample = pred_sample,
       p_shape = p_shape, p_level = p_level, ok = TRUE)
}

# =============================================================================
# Global: spatial score field + per-sample histogram
# =============================================================================
sample_ids <- levels(seu$sample)
max_x_span <- max(tapply(all_coords[, 1], samples, function(x) diff(range(x))))
max_y_span <- max(tapply(all_coords[, 2], samples, function(x) diff(range(x))))
half_span  <- max(max_x_span, max_y_span) / 2 * 1.05

plots_field <- lapply(sample_ids, function(s) {
  idx <- which(samples == s)
  df_plot <- data.frame(x = all_coords[idx, 1], y = all_coords[idx, 2],
                        score = tumor_score[idx])
  ggplot(df_plot, aes(x = x, y = y, color = score)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_viridis_c(option = "magma", limits = c(0, 1), name = "Tumor\nScore") +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    ggtitle(s) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
})
p_field <- wrap_plots(plots_field, nrow = 1) +
  plot_annotation(title = sprintf("Tumor Score Field (%s)", toupper(SCORE_METHOD)),
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_score_field.pdf"),
       p_field, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_score_field.png"),
       p_field, width = 4 * length(sample_ids), height = 5, dpi = 300)

df_hist <- data.frame(score = tumor_score[!is.na(tumor_score)],
                      sample = samples[!is.na(tumor_score)])
p_hist <- ggplot(df_hist, aes(x = score, fill = sample)) +
  geom_histogram(bins = 60, color = NA) +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = SCORE_LAB, y = "Cell count",
       title = "Tumor Score Distribution per Sample") +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "tumor_score_histogram_per_sample.pdf"),
       p_hist, width = 4 * ceiling(length(sample_ids) / 2), height = 6)

# =============================================================================
# Per cell type
# =============================================================================
for (ct in CELLTYPES_TRAJECTORY) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
  ct_dir  <- file.path(OUTPUT_DIR, ct_safe)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n=== %s ===\n", ct))

  # Select this cell type's cells that also have a valid score and pass QC
  # (qc_pass is all-TRUE when APPLY_QC = FALSE).
  ct_idx <- which(celltypes == ct & !is.na(tumor_score) & qc_pass)
  if (length(ct_idx) < 50) {
    cat(sprintf("  Skipping %s: too few cells (%d)\n", ct, length(ct_idx)))
    next
  }

  ct_score     <- tumor_score[ct_idx]
  # Evenly-spaced decile coordinate: each cell's empirical quantile of the score
  # within this cell type, in [0,1]. Deciles land at ~0.05,0.15,...,0.95. Used
  # for the trajectory x-axis so the skewed raw score is not squished left.
  ct_score_dec <- ecdf(ct_score)(ct_score)
  ct_condition <- droplevels(seu$condition[ct_idx])
  ct_sample    <- droplevels(seu$sample[ct_idx])
  expr_mat     <- GetAssayData(seu, assay = "RNA", layer = "data")[, ct_idx]

  cat(sprintf("  Total cells: %d (WT: %d, DKO: %d)\n", length(ct_idx),
              sum(ct_condition == REFERENCE_CONDITION),
              sum(ct_condition == TEST_CONDITION)))

  # --- Adaptive N_BINS (optional) --------------------------------------------
  # Choose N_BINS for THIS cell type as the finest binning the data support (at
  # most ADAPTIVE_BIN_TOL_K invalid bins; see the ADAPTIVE_N_BINS docs). Reassign
  # the module-level N_BINS so every downstream make_bins()/make_pseudobulk()/
  # diagnostic call in this iteration picks it up via their `n_bins = N_BINS`
  # default. Recomputed at the top of each cell-type iteration.
  if (ADAPTIVE_N_BINS) {
    n_sel <- choose_adaptive_nbins(ct_score, ct_sample, ct_condition)
    write.csv(attr(n_sel, "scan"),
              file.path(ct_dir, "adaptive_nbins_scan.csv"), row.names = FALSE)
    cat(sprintf("  Adaptive N_BINS: selected %d bins (tol k=%d, scanned %d-%d, %d valid)\n",
                n_sel, ADAPTIVE_BIN_TOL_K, ADAPTIVE_N_BINS_MIN, ADAPTIVE_N_BINS_MAX,
                attr(n_sel, "scan")$n_valid[match(n_sel, attr(n_sel, "scan")$n_bins)]))
    N_BINS <- as.integer(n_sel)
  }

  # --- Gene selection --------------------------------------------------------
  # TWO distinct sets (reviewer point #5):
  #   test_genes — baseline-expression filter, IDENTICAL to the DESeq2 DE scripts
  #                (00_parameters.R): detected in >= MIN_GENE_EXPR_FRAC of this cell
  #                type's cells AND >= MIN_COUNTS_PER_SAMPLE counts in
  #                >= MIN_SAMPLES_EXPRESSED samples. Independent of the condition x
  #                score interaction, so BH over this whole set uses the correct m.
  #   plot_genes — score-dependent GENE_SELECTION (binned_variance top-N). Used ONLY
  #                for the reference/wave heatmaps, never to define the tested set.
  counts_ct <- GetAssayData(seu, assay = "RNA", layer = "counts")[, ct_idx]
  det_frac  <- Matrix::rowSums(counts_ct > 0) / length(ct_idx)
  samp_sum  <- vapply(levels(ct_sample), function(s)
                 Matrix::rowSums(counts_ct[, ct_sample == s, drop = FALSE]),
                 numeric(nrow(counts_ct)))               # genes x samples (summed counts)
  n_samp_det <- rowSums(samp_sum >= MIN_COUNTS_PER_SAMPLE)
  test_genes <- rownames(counts_ct)[det_frac >= MIN_GENE_EXPR_FRAC &
                                    n_samp_det >= MIN_SAMPLES_EXPRESSED]
  test_genes <- intersect(test_genes, rownames(expr_mat))
  plot_genes <- select_genes(expr_mat, ct_score)
  if (length(test_genes) < 5 || length(plot_genes) < 5) {
    cat(sprintf("  Skipping %s: too few genes (test %d, plot %d)\n",
                ct, length(test_genes), length(plot_genes))); next
  }
  cat(sprintf("  Genes: %d tested (>=%.0f%% cells & >=%d counts in >=%d samples), %d shown in heatmaps\n",
              length(test_genes), MIN_GENE_EXPR_FRAC * 100, MIN_COUNTS_PER_SAMPLE,
              MIN_SAMPLES_EXPRESSED, length(plot_genes)))

  # ---------------------------------------------------------------------------
  # Binning-quality diagnostics
  # ---------------------------------------------------------------------------
  diag_dir <- file.path(ct_dir, "diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  bininfo <- make_bins(ct_score)
  df_score <- data.frame(score = ct_score, sample = ct_sample, condition = ct_condition,
                         bin = as.integer(bininfo$bins))

  # (a) Score histogram with decile bin breaks
  p_bh <- ggplot(df_score, aes(x = score)) +
    geom_histogram(bins = 60, fill = "#5E35B1", color = NA, alpha = 0.8) +
    geom_vline(xintercept = bininfo$breaks, color = "grey30",
               linetype = "dashed", linewidth = 0.3) +
    labs(x = SCORE_LAB, y = "Cell count",
         title = sprintf("%s score distribution + decile breaks", ct)) +
    theme_publication()
  ggsave(file.path(diag_dir, "score_histogram_with_bins.pdf"), p_bh, width = 6, height = 4)

  # (b) n_cells per sample x bin (tile) + threshold flag
  counts <- df_score %>% count(sample, bin, name = "n_cells") %>%
    tidyr::complete(sample, bin = seq_len(N_BINS), fill = list(n_cells = 0)) %>%
    mutate(pass = n_cells >= MIN_CELLS_PER_BIN)
  p_tile <- ggplot(counts, aes(x = bin, y = sample, fill = n_cells)) +
    geom_tile(color = "white") +
    geom_text(aes(label = n_cells, color = pass), size = 2.6, show.legend = FALSE) +
    scale_fill_viridis_c(option = "mako", name = "n cells") +
    scale_color_manual(values = c(`TRUE` = "white", `FALSE` = "grey60")) +
    scale_x_continuous(breaks = seq_len(N_BINS)) +
    labs(x = "Score bin (low -> high)", y = NULL,
         title = sprintf("%s — Cells per sample x bin (>= %d passes)", ct, MIN_CELLS_PER_BIN)) +
    theme_publication()
  ggsave(file.path(diag_dir, "cells_per_sample_bin.pdf"), p_tile, width = 7, height = 4)

  # (c) # bins per sample meeting min-cells threshold
  bins_pass <- counts %>% group_by(sample) %>%
    summarise(n_bins_pass = sum(pass), .groups = "drop")
  p_pass <- ggplot(bins_pass, aes(x = sample, y = n_bins_pass, fill = sample)) +
    geom_col() +
    scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
    geom_hline(yintercept = N_BINS, linetype = "dashed", color = "grey40") +
    labs(x = NULL, y = "# bins with enough cells",
         title = sprintf("%s usable bins per sample (max %d)", ct, N_BINS)) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(diag_dir, "bins_passing_per_sample.pdf"), p_pass, width = 5, height = 4)

  write.csv(counts, file.path(diag_dir, "cells_per_sample_bin.csv"), row.names = FALSE)

  # ---------------------------------------------------------------------------
  # Binned heatmaps (combined + per-condition) — reference view
  # ---------------------------------------------------------------------------
  be_all <- binned_expression(expr_mat[plot_genes, , drop = FALSE], bininfo$bins)
  hm_scaled <- t(scale(t(be_all)))
  hm_df <- as.data.frame(hm_scaled) %>%
    mutate(gene = rownames(hm_scaled)) %>%
    pivot_longer(-gene, names_to = "bin", values_to = "zscore") %>%
    mutate(bin_num = as.integer(gsub("bin_", "", bin)),
           gene = factor(gene, levels = rev(plot_genes)))
  p_hm <- ggplot(hm_df, aes(x = bin_num, y = gene, fill = zscore)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno", name = "Z-score") +
    scale_x_continuous(breaks = c(1, N_BINS), labels = c("Low", "High")) +
    labs(x = TRAJ_X_LAB, y = NULL,
         title = sprintf("%s Expression along tumor score", ct)) +
    theme_publication() + theme(axis.text.y = element_text(size = 7))
  # limitsize = FALSE: height scales with gene count and can exceed ggplot's
  # 50-inch guard at large TOP_N_GENES (e.g. 200 -> 50 in).
  ggsave(file.path(ct_dir, "heatmap_trajectory_combined.pdf"),
         p_hm, width = 6, height = max(4, length(plot_genes) * 0.25),
         limitsize = FALSE)

  # ---------------------------------------------------------------------------
  # Normalized-expression-along-gradient heatmap
  # Per-gene expression is smoothed against the CONTINUOUS score (NO binning)
  # with a univariate GAM, evaluated on a common grid, then min-max normalized
  # per gene. Shown as a heatmap; genes ordered by the score at which the
  # combined profile peaks (forms a cascade). Combined + WT + DKO facets.
  # ---------------------------------------------------------------------------
  norm01 <- function(v) { r <- range(v, na.rm = TRUE)
                          if (diff(r) > 0) (v - r[1]) / diff(r) else rep(0.5, length(v)) }

  # Smooth against the evenly-spaced decile coordinate (not the raw score) so the
  # grid spans deciles uniformly. Grid is fixed [0,1] and shared across facets.
  score_grid <- seq(0, 1, length.out = 100)

  # Smooth one gene vs the decile coordinate; return fitted values on score_grid.
  smooth_gene <- function(y, x, k = 10) {
    fit <- tryCatch(gam(expr ~ s(score, bs = "cr", k = k),
                        data = data.frame(expr = y, score = x), method = "REML"),
                    error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, length(score_grid)))
    as.numeric(predict(fit, newdata = data.frame(score = score_grid)))
  }

  build_profiles <- function(cells) {
    m <- t(vapply(plot_genes,
                  function(g) norm01(smooth_gene(as.numeric(expr_mat[g, cells]),
                                                 ct_score_dec[cells])),
                  numeric(length(score_grid))))
    rownames(m) <- plot_genes
    m
  }

  prof_to_df <- function(m, facet) {
    colnames(m) <- score_grid
    as.data.frame(m) %>%
      mutate(gene = rownames(m)) %>%
      pivot_longer(-gene, names_to = "score", values_to = "norm_expr") %>%
      mutate(score = as.numeric(score), facet = facet)
  }

  prof_all   <- build_profiles(seq_along(ct_idx))
  # Cascade order: sort genes by the score at which their combined profile peaks
  gene_order <- rownames(prof_all)[order(apply(prof_all, 1, which.max))]

  hm_list <- list(prof_to_df(prof_all, "Combined"))
  for (cond in levels(ct_condition)) {
    cc <- which(ct_condition == cond)
    if (length(cc) >= 30) hm_list[[cond]] <- prof_to_df(build_profiles(cc), cond)
  }
  prof_df <- bind_rows(hm_list)
  prof_df$facet <- factor(prof_df$facet,
                          levels = c("Combined", REFERENCE_CONDITION, TEST_CONDITION))
  prof_df$gene  <- factor(prof_df$gene, levels = gene_order)

  p_wave <- ggplot(prof_df, aes(x = score, y = gene, fill = norm_expr)) +
    geom_tile() +
    facet_wrap(~ facet, nrow = 1) +
    scale_fill_viridis_c(option = "inferno", name = "Norm.\nexpr", limits = c(0, 1)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(x = TRAJ_X_LAB, y = NULL,
         title = sprintf("%s Normalized expression along tumor gradient", ct)) +
    theme_publication() +
    theme(axis.text.y = element_text(size = 7))
  ggsave(file.path(ct_dir, "normalized_expression_heatmap.pdf"),
         p_wave, width = 4 * nlevels(prof_df$facet),
         height = max(4, length(plot_genes) * 0.22 + 1), limitsize = FALSE)

  # ---------------------------------------------------------------------------
  # GAM fits, loop over selected models
  # ---------------------------------------------------------------------------
  # Genes ranked by the pseudobulk shape test; populated when "pseudobulk" runs
  # and used to restrict the single-cell (fs/ranef) fits to the top hits.
  pseudobulk_top_genes <- NULL
  for (model in GAM_MODELS) {
    gam_dir <- file.path(ct_dir, paste0("gam_", model))
    dir.create(gam_dir, recursive = TRUE, showWarnings = FALSE)

    # Gene set for this model. Models in FULL_GENE_MODELS use the full
    # expression-filtered test_genes set (correct BH denominator, reviewer #5);
    # any other model (e.g. fs) re-fits only the top N_GENES_OTHER_MODELS
    # pseudobulk hits so the expensive fits stay focused.
    if (model %in% FULL_GENE_MODELS) {
      genes_this_model <- test_genes
    } else if (!is.null(pseudobulk_top_genes)) {
      genes_this_model <- head(pseudobulk_top_genes, N_GENES_OTHER_MODELS)
    } else {
      # Pseudobulk not run this session, fall back to the first N test_genes.
      genes_this_model <- head(test_genes, N_GENES_OTHER_MODELS)
    }
    cat(sprintf("  [%s] fitting %d genes\n", model, length(genes_this_model)))

    # Build the data the model needs
    if (model %in% PSEUDOBULK_MODELS) {
      # NB specs sum RAW COUNTS (+ carry a library-size offset); the Gaussian specs
      # average the log-normalized data layer.
      pb <- if (model %in% PSEUDOBULK_NB_MODELS)
              make_pseudobulk(counts_ct, ct_score, ct_sample, ct_condition, agg = "sum")
            else make_pseudobulk(expr_mat, ct_score, ct_sample, ct_condition, agg = "mean")
      if (is.null(pb) || nrow(pb$meta) < 10) {
        cat(sprintf("  [%s] too few pseudobulk groups (common support too small?), skipping\n", model)); next
      }
      # Report the common-support truncation: which bins / score range survived.
      if (!is.null(pb$support)) {
        kb <- pb$support$keep_bins
        cat(sprintf("  [%s] common support: %d/%d bins retained (bins %s); score range [%.3f, %.3f]\n",
                    model, length(kb), N_BINS, paste(range(kb), collapse = "-"),
                    min(pb$meta$score), max(pb$meta$score)))
      }
    }

    gam_summary <- list()
    fit_curves  <- list()   # per-gene population-level fitted curves (for heatmap)
    plot_store  <- list()   # per-gene curve+points, kept only to plot SIGNIFICANT genes
    for (gene in genes_this_model) {
      if (model %in% PSEUDOBULK_MODELS) {
        # x-axis per PB_SCORE_AXIS (reviewer #2): "score" = REAL mean local tumor
        # fraction of the group (interpretable units, keeps within-bin variation);
        # "decile" = evenly-spaced bin index (equal plot width per decile).
        pb_x <- if (PB_SCORE_AXIS == "score") pb$meta$score
                else (pb$meta$bin - 0.5) / N_BINS
        df_gam <- data.frame(expr = as.numeric(pb$expr[gene, ]),
                             score = pb_x,
                             condition = pb$meta$condition,
                             sample = pb$meta$sample, weight = pb$meta$n_cells)
        # NB-GAM log-offset = log of the group's total counts (floored at 1 to guard
        # against log(0)). expr here is summed raw counts, not mean log-norm.
        if (model %in% PSEUDOBULK_NB_MODELS)
          df_gam$log_libsize <- log(pmax(pb$meta$libsize, 1))
      } else {
        df_gam <- data.frame(expr = as.numeric(expr_mat[gene, ]),
                             score = ct_score_dec, condition = ct_condition,
                             sample = ct_sample)
      }

      res <- fit_trajectory_gam(df_gam, model)
      if (!isTRUE(res$ok)) next

      gam_summary[[gene]] <- data.frame(gene = gene, p_shape = res$p_shape,
                                        p_level = res$p_level)
      fit_curves[[gene]]  <- res$pred[, c("score", "condition", "fit")]

      # Stash the full fitted curve (with CI) + underlying points so the per-gene
      # spline PDFs can be emitted AFTER FDR is known, for the significant genes only.
      if (model %in% PLOT_INDIVIDUAL_GAMS) {
        plot_store[[gene]] <- list(pred = res$pred, pred_sample = res$pred_sample,
                                   df = df_gam,
                                   p_shape = res$p_shape, p_level = res$p_level)
      }
    }

    if (length(gam_summary) > 0) {
      gs <- bind_rows(gam_summary)
      # BH-FDR correction across all tested genes, separately for shape & level.
      gs$fdr_shape <- p.adjust(gs$p_shape, method = "BH")
      gs$fdr_level <- p.adjust(gs$p_level, method = "BH")
      gs <- gs[order(gs$p_shape), ]   # rank by trajectory-shape difference
      write.csv(gs, file.path(ct_dir, sprintf("gam_condition_test_%s.csv", model)),
                row.names = FALSE)
      cat(sprintf("  [%s] GAM results: %d genes\n", model, nrow(gs)))

      # Record the pseudobulk shape ranking so fs/ranef (later in the loop) can
      # re-fit its top hits.
      if (model == "pseudobulk") pseudobulk_top_genes <- gs$gene

      # -----------------------------------------------------------------------
      # Per-gene spline PDFs for the SIGNIFICANT genes only (BH-FDR < SIG_FDR on
      # shape or level) — the readable subset, emitted now that FDR is known.
      # -----------------------------------------------------------------------
      if (model %in% PLOT_INDIVIDUAL_GAMS) {
        sig_genes <- gs$gene[which(gs$fdr_shape < SIG_FDR | gs$fdr_level < SIG_FDR)]
        sig_genes <- intersect(sig_genes, names(plot_store))
        cat(sprintf("  [%s] per-gene spline PDFs for %d significant genes (FDR < %.2f)\n",
                    model, length(sig_genes), SIG_FDR))
        ylab <- if (model %in% PSEUDOBULK_NB_MODELS)
                  "log E[counts] (NB-GAM, depth-adjusted)"
                else if (model %in% PSEUDOBULK_MODELS)
                  "Mean expression (pseudobulk)" else "Expression"
        xlab_gam <- if (model %in% PSEUDOBULK_MODELS) GAM_X_LAB else TRAJ_X_LAB
        # Per-sample fitted-curve overlays go in their own subfolder (one file per
        # gene): WT/DKO population curves + each sample's own random-effect curve.
        ps_dir <- file.path(gam_dir, "per_sample_curves")
        dir.create(ps_dir, recursive = TRUE, showWarnings = FALSE)
        # Combined palette so the single color scale carries BOTH the per-sample
        # lines (SAMPLE_COLORS) and the two thick population lines (condition colors).
        pop_keys <- setNames(
          c(CONDITION_COLORS[[REFERENCE_CONDITION]], CONDITION_COLORS[[TEST_CONDITION]]),
          c(paste0(REFERENCE_CONDITION, " (fitted mean)"),
            paste0(TEST_CONDITION, " (fitted mean)")))
        ps_pal <- c(SAMPLE_COLORS, pop_keys)
        for (gene in sig_genes) {
          st <- plot_store[[gene]]
          subt <- paste(c(
            if (!is.na(st$p_shape)) sprintf("shape p = %.3g", st$p_shape),
            if (!is.na(st$p_level)) sprintf("level p = %.3g", st$p_level)
          ), collapse = "   ")
          # For NB the fitted curve is on the log-E[count] scale at the MEAN library
          # size; put the observed summed counts on that same depth-adjusted log scale
          # (log(count+0.5) - log_libsize + mean log_libsize) so they overlay the curve.
          pts_df <- st$df
          if (model %in% PSEUDOBULK_NB_MODELS)
            pts_df$expr <- log(pts_df$expr + 0.5) - pts_df$log_libsize +
                           mean(pts_df$log_libsize)
          p_gam <- ggplot(st$pred, aes(x = score, y = fit, color = condition, fill = condition)) +
            geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
            geom_line(linewidth = 1) +
            #{ if (model %in% PSEUDOBULK_MODELS)
            #  geom_point(data = pts_df, aes(x = score, y = expr, color = condition, size = weight),
            #             alpha = 0.5, inherit.aes = FALSE) 
            #               } +
            scale_color_manual(values = CONDITION_COLORS) +
            scale_fill_manual(values = CONDITION_COLORS) +
            scale_size_continuous(range = c(1, 4), guide = "none") +
            labs(x = xlab_gam, y = ylab,
                 title = sprintf("%s in %s [%s]", gene, ct, model),
                 subtitle = subt) +
            theme_publication()
          ggsave(file.path(gam_dir, sprintf("gam_%s.pdf", gene)), p_gam, width = 5, height = 4)

          # --- Per-sample fitted-curve overlay -------------------------------
          # Same fitted model, but showing EACH sample's own trajectory (random
          # effect included) behind the WT/DKO population means, so the
          # between-mouse spread the condition test accounts for is visible.
          ps_df <- st$pred_sample
          if (!is.null(ps_df) && nrow(ps_df) > 0 && any(is.finite(ps_df$fit))) {
            pop_df <- st$pred
            pop_df$key <- factor(paste0(as.character(pop_df$condition), " (fitted mean)"),
                                 levels = names(pop_keys))
            p_ps <- ggplot() +
              geom_line(data = ps_df, aes(x = score, y = fit, color = sample, group = sample),
                        linewidth = 0.6, alpha = 0.85) +
              geom_line(data = pop_df, aes(x = score, y = fit, color = key, group = key),
                        linewidth = 1.7) +
              scale_color_manual(values = ps_pal, name = NULL) +
              labs(x = xlab_gam, y = ylab,
                   title = sprintf("%s in %s [%s] — per-sample fits", gene, ct, model),
                   subtitle = subt) +
              theme_publication()
            ggsave(file.path(ps_dir, sprintf("gam_%s_per_sample.pdf", gene)),
                   p_ps, width = 5.5, height = 4)
          }
        }
      }

      # -----------------------------------------------------------------------
      # Significant-trajectory heatmap: WT | DKO GAM-fitted curves, only genes
      # passing BH-FDR < SIG_FDR on shape OR level. Rows are per-gene z-scored
      # across BOTH conditions on a shared scale, so a *level* shift shows as one
      # panel uniformly brighter and a *shape* change as the two panels diverging
      # along the axis. This is the readable replacement for the per-gene curves.
      # -----------------------------------------------------------------------
      sig <- gs[which(gs$fdr_shape < SIG_FDR | gs$fdr_level < SIG_FDR), ]
      if (nrow(sig) >= 3) {
        genes_keep <- sig$gene
        ng <- sum(fit_curves[[genes_keep[1]]]$condition == REFERENCE_CONDITION)
        get_curve <- function(g, cond) {
          d <- fit_curves[[g]]; d <- d[d$condition == cond, ]
          d$fit[order(d$score)]
        }
        mat_wt  <- t(vapply(genes_keep, get_curve, numeric(ng), cond = REFERENCE_CONDITION))
        mat_dko <- t(vapply(genes_keep, get_curve, numeric(ng), cond = TEST_CONDITION))
        rownames(mat_wt) <- rownames(mat_dko) <- genes_keep

        # Shared per-gene z-score across both conditions (preserves level + shape)
        comb <- cbind(mat_wt, mat_dko)
        rmn  <- rowMeans(comb, na.rm = TRUE)
        rsd  <- apply(comb, 1, sd, na.rm = TRUE); rsd[rsd == 0 | is.na(rsd)] <- 1
        z_wt  <- (mat_wt  - rmn) / rsd
        z_dko <- (mat_dko - rmn) / rsd

        # Row ordering: hierarchical clustering on the fitted WT+DKO curves
        # (default, groups co-varying trajectories) or ranked by p_shape.
        if (HEATMAP_ROW_ORDER == "shape") {
          ord   <- order(sig$p_shape)
          z_wt  <- z_wt[ord, , drop = FALSE]; z_dko <- z_dko[ord, , drop = FALSE]
          sig   <- sig[ord, ]
          cluster_rows_arg <- FALSE
        } else {  # "cluster"
          cluster_rows_arg <- as.dendrogram(hclust(dist(cbind(z_wt, z_dko))))
        }

        category <- ifelse(sig$fdr_shape < SIG_FDR & sig$fdr_level < SIG_FDR, "both",
                    ifelse(sig$fdr_shape < SIG_FDR, "shape", "level"))
        cat_cols <- c(shape = "#D55E00", level = "#0072B2", both = "#009E73")
        col_fun  <- colorRamp2(c(-2, -1, 0, 1, 2), viridis::inferno(5))

        # Display -log10 FDR, capped at 10 so a single near-zero p (e.g. Bcl2)
        # does not flatten every other bar. Exact values are in the CSV.
        nlog_cap <- function(p) pmin(-log10(pmax(p, 1e-300)), 10)
        row_ha <- rowAnnotation(
          class = category,
          `-log10 FDRshape` = anno_barplot(nlog_cap(sig$fdr_shape),
                                           gp = gpar(fill = "#D55E00"), border = FALSE),
          `-log10 FDRlevel` = anno_barplot(nlog_cap(sig$fdr_level),
                                           gp = gpar(fill = "#0072B2"), border = FALSE),
          col = list(class = cat_cols),
          annotation_name_gp = gpar(fontsize = 7),
          annotation_legend_param = list(class = list(title = "Sig. in")),
          width = unit(3, "cm"))

        col_note <- HeatmapAnnotation(
          axis = anno_text(c("periphery", rep("", ng - 2), "tumor core"),
                           gp = gpar(fontsize = 7), just = "center"),
          show_annotation_name = FALSE)

        ht_wt <- Heatmap(z_wt, cluster_rows = cluster_rows_arg, cluster_columns = FALSE,
                         col = col_fun, name = "Fitted\nz-score",
                         column_title = REFERENCE_CONDITION,
                         show_column_names = FALSE, row_names_side = "left",
                         row_dend_side = "left", show_row_dend = TRUE,
                         row_dend_width = unit(0.8, "cm"),
                         row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                         bottom_annotation = col_note,
                         width = unit(6.5, "cm"), border = TRUE)
        ht_dko <- Heatmap(z_dko, cluster_rows = FALSE, cluster_columns = FALSE,
                          col = col_fun, name = " ", column_title = TEST_CONDITION,
                          show_column_names = FALSE, show_row_names = FALSE,
                          bottom_annotation = col_note,
                          width = unit(6.5, "cm"), border = TRUE,
                          show_heatmap_legend = FALSE)

        pdf(file.path(ct_dir, sprintf("trajectory_heatmap_significant_%s.pdf", model)),
            width = 13, height = max(3.5, nrow(z_wt) * 0.13 + 1.5))
        draw(ht_wt + ht_dko + row_ha,
             column_title = sprintf("%s significant tumor-proximity trajectories (BH-FDR < %.2f, %d genes)",
                                    ct, SIG_FDR, nrow(z_wt)),
             column_title_gp = gpar(fontsize = 12, fontface = "bold"))
        dev.off()
        cat(sprintf("  [%s] significant-trajectory heatmap: %d genes (%d shape, %d level, %d both)\n",
                    model, nrow(z_wt), sum(category == "shape"),
                    sum(category == "level"), sum(category == "both")))
      } else {
        cat(sprintf("  [%s] < 3 genes pass BH-FDR < %.2f; heatmap skipped\n", model, SIG_FDR))
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Density of cells along the gradient (WT vs DKO)
  # ---------------------------------------------------------------------------
  df_dist <- data.frame(score = ct_score, condition = ct_condition)
  p_dist <- ggplot(df_dist, aes(x = score, fill = condition, color = condition)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(x = SCORE_LAB, y = "Density of cells",
         title = sprintf("%s position along tumor gradient", ct)) +
    theme_publication()
  ggsave(file.path(ct_dir, "density_distribution_along_gradient.pdf"),
         p_dist, width = 5, height = 4)

  cat(sprintf("  Done: %s\n", ct))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")
