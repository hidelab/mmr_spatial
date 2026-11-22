# ==============================================================================
# 15_06_tumor_trajectory_working_revised.R
# Revision of 15_04_tumor_trajectory_working.R
# ==============================================================================
#
# SUMMARY OF EDITS
# ==============================================================================
#
# --- A. THE COVARIATE (see 00_utils2.R for the scoring changes) --------------
#
# A1. Scoring switched to compute_tumor_score2(..., laplacian = "rw").
#     Still heat diffusion / backward Euler, but on the random-walk Laplacian
#     rather than L_sym, so the diffusion operator is row-stochastic and the
#     score is a true weighted local tumor fraction in [0,1] with no degree or
#     tissue-boundary bias. (The old L_sym operator was NOT the problem I first
#     suspected -- a permuted-label control showed it tracks tumor labels
#     correctly -- but it does not preserve constants, so it is biased upward
#     in dense regions and downward at tissue edges.)
#
# A2. GLOBAL kernel bandwidth. Previously sigma <- median(nn.dist) was computed
#     per sample inside score_diffusion(). Measured: holding tumor geometry
#     fixed and changing only cellularity 3x moved sigma 12.68 -> 7.26, i.e. the
#     covariate's physical length scale was a per-sample quantity, and therefore
#     condition-confounded if WT and DKO differ in cellularity (likely, given
#     the readout is CD8 infiltration). Now one bandwidth for the whole cohort.
#
# A3. MIN-MAX NORMALIZATION REMOVED. Previously every sample was rescaled to
#     [0,1]. Measured: a big tumor (frac 0.50, raw IQR 0.925) and a small tumor
#     (frac 0.13, raw IQR 0.176) both ended up spanning exactly [0,1], so tumor
#     burden was normalized away and both endpoints were set by single extreme
#     cells. The rw score needs no rescaling.
#
# A4. POOLED-ECDF X-AXIS REMOVED. Previously ct_score_dec <- ecdf(ct_score) was
#     computed over WT and DKO pooled, so "decile 5" was a different absolute
#     tumor proximity in each condition. Any redistribution of CD8s along the
#     gradient -- a plausible primary DKO phenotype -- would have been read out
#     as a trajectory-SHAPE difference. The decile coordinate is now used for
#     PLOTTING ONLY and never enters a model.
#
# A5. BIN INDEX X-AXIS REMOVED. The pseudobulk model previously used
#     score = (bin - 0.5)/N_BINS, i.e. evenly spaced bin indices. It now uses
#     the actual MEAN SCORE of the cells in each sample x bin group.
#
# A6. Diagnostics are now mandatory and printed/plotted before any modelling:
#     per-sample tumor fraction, nn-distance median, raw score range/IQR,
#     tested against condition. Plus a permuted-tumor-label null check on the
#     score field itself. If any diagnostic separates by condition, the
#     condition x trajectory interaction is not interpretable and the script
#     says so.
#
# --- B. THE MODEL ------------------------------------------------------------
#
# B1. PRIMARY SPEC IS NOW "fs" (factor-smooth), not "pseudobulk"/"ranef" with a
#     random intercept. s(sample, bs="re") absorbs sample-level LEVEL shifts
#     only; it does nothing for sample-specific curve SHAPE, so the difference
#     smooth s(score, by=cond_ord) was being tested against variation across
#     pseudobulk bins rather than across the 7 mice. s(score, sample, bs="fs")
#     gives each mouse its own shrunk smooth and subsumes the random intercept
#     (the redundant s(sample, bs="re") term is dropped when fs is present).
#     "ranef" is retained as a deliberately anticonservative diagnostic so the
#     gap between the two is visible.
#
# B2. NEGATIVE BINOMIAL + LIBRARY SIZE OFFSET. Previously gam() defaulted to
#     gaussian() on log-normalized values. Pseudobulk now SUMS RAW COUNTS per
#     sample x bin (mean-of-log-normalized is biased, and the bias depends on
#     per-cell library sizes, which trend with tumor proximity) and fits
#     family = nb() with offset(log(total_counts)).
#
# B3. p.table / s.table are indexed BY POSITION, not by the column name
#     "Pr(>|t|)", which becomes "Pr(>|z|)" under nb() and would have silently
#     returned NA.
#
# B4. GAM_K raised from 5 to 10 and k.check() is run and logged. k=5 caps the
#     fit at ~4 edf, which smears exactly the sharp margin transition that is
#     the interesting feature. The descriptive "wave" plot previously used
#     k=10 while the test used k=5; these are now harmonized.
#
# B5. BIN COVERAGE ASSERTION. MIN_CELLS_PER_BIN dropped sample x bin cells
#     non-randomly, so the difference smooth at the tumor-core end could rest
#     on one or two mice. Bins are now retained only where BOTH conditions
#     contribute >= MIN_SAMPLES_PER_BIN samples, and the score range is
#     truncated accordingly. Dropped bins are reported.
#
# --- C. INFERENCE ------------------------------------------------------------
#
# C1. GENE SELECTION DECOUPLED FROM TESTING. Previously the 400 genes with the
#     highest binned variance ALONG THE SCORE were selected, then BH-corrected
#     within those 400 -- selection and the interaction test were dependent and
#     the FDR denominator was wrong. Genes are now filtered on EXPRESSION ONLY
#     and BH runs over the full tested set. binned_variance is retained solely
#     for ordering genes in plots.
#
# C2. SAMPLE-LEVEL PERMUTATION TEST. With 3 WT + 4 DKO there are choose(7,3)=35
#     distinct condition assignments, so the smallest achievable permutation
#     p-value is 1/35 = 0.029 and NO gene can reach a permutation-based
#     FDR < 0.1 across hundreds of genes. That is a property of the design.
#     The script therefore reports:
#       (a) a GLOBAL test -- number of genes passing the nominal cutoff in the
#           real assignment vs in all 35 assignments. This supports a statement
#           about the cohort ("the DKO trajectory differs", p >= 1/35) which is
#           the strongest claim n=3v4 licenses.
#       (b) EMPIRICAL per-gene p-values, calibrated against the pooled null of
#           shape statistics from the 34 non-identity assignments, instead of
#           mgcv's asymptotic reference (which is conditional on variance
#           components estimated from 7 samples and has no such floor).
#
# C3. Per-gene results are labelled as a RANKING, and the heatmap annotation is
#     relabelled accordingly, unless the global permutation test is significant.
#
# --- D. SMALLER FIXES --------------------------------------------------------
#
# D1. make_bins() previously fell back SILENTLY from quantile to equal-width
#     breaks when the score had ties. Now warns and logs the breaks.
# D2. A single primary model is declared up front (PRIMARY_MODEL) so two
#     full-gene specs are not both written out and chosen after the fact.
# D3. fit_trajectory_gam() gains an else-branch for unrecognized model names.
# D4. alpha (diffusion step size) is exposed in the parameter block; previously
#     it was hard-coded at 2.0 inside score_diffusion() so the effective
#     diffusion time (alpha * steps = 40) appeared nowhere in the config.
# D5. Single-cell models use bam(discrete=TRUE) rather than gam().
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
library(parallel)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))
#source(file.path(PROJECT_DIR, "code", "R", "00_utils2.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "15_06_tumor_trajectory_revised")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters
# =============================================================================
# --- Tumor score (A1-A3) ------------------------------------------------------
SCORE_METHOD    <- "diffusion"   # "diffusion" | "kernel_fraction" | "knn"
SCORE_GRAPH     <- "knn"         # "knn" | "radius"  <- "radius" recommended
SCORE_LAPLACIAN <- "rw"          # "rw" (correct) | "sym" (legacy, ratio-corrected)
KNN_K           <- 30
DIFF_ALPHA      <- 2.0           # D4: exposed. Diffusion time = ALPHA * STEPS
DIFF_STEPS      <- 20
SCORE_SIGMA     <- NULL          # NULL = one global bandwidth pooled over samples
SCORE_RADIUS    <- NULL          # NULL = 3 * sigma (only used if graph="radius")
TARGET_CELLTYPE <- "Tumor"

# --- Binning (B5) -------------------------------------------------------------
N_BINS              <- 20
MIN_CELLS_PER_BIN   <- 5
MIN_SAMPLES_PER_BIN <- 2   # per condition; bins failing this are dropped

# --- Genes (C1) ---------------------------------------------------------------
MIN_MEAN_EXPR     <- 0.05  # expression-only filter; NOT score-dependent
MIN_PCT_DETECTED  <- 0.05
TOP_N_GENES_PLOT  <- 40   # for plot ordering only, never for the FDR denominator

# --- Model (B1, B2, B4) -------------------------------------------------------
PRIMARY_MODEL   <- "fs"                 # D2: declared up front
GAM_MODELS      <- c("fs", "ranef")     # "ranef" is the anticonservative contrast
GAM_K           <- 10                   # B4
GAM_FAMILY      <- "nb"                 # "nb" | "gaussian" (legacy)

# --- Inference (C2) -----------------------------------------------------------
RUN_PERMUTATION <- TRUE
N_PERM_CORES    <- 4
SIG_FDR         <- 0.1
HEATMAP_ROW_ORDER <- "cluster"

CELLTYPES_TRAJECTORY <- c("CD8+ T")
APPLY_QC <- FALSE

stopifnot(PRIMARY_MODEL %in% GAM_MODELS,
          all(GAM_MODELS %in% c("fs", "ranef")),
          GAM_FAMILY %in% c("nb", "gaussian"))

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

qc_pass <- if (APPLY_QC) {
  seu$total_counts >= MIN_TOTAL_COUNTS & seu$n_genes_by_counts >= MIN_GENES_DETECTED
} else rep(TRUE, ncol(seu))
cat(sprintf("QC: %d / %d cells pass\n", sum(qc_pass), length(qc_pass)))

# =============================================================================
# Tumor score (A1-A3) + mandatory diagnostics (A6)
# =============================================================================
all_coords <- Embeddings(seu, "spatial")
celltypes  <- seu$Tier1_celltype
samples    <- seu$sample
is_tumor   <- celltypes == TARGET_CELLTYPE

tumor_score <- compute_tumor_score2(
  all_coords, is_tumor, samples,
  method = SCORE_METHOD, graph = SCORE_GRAPH, laplacian = SCORE_LAPLACIAN,
  sigma = SCORE_SIGMA, k = KNN_K, radius = SCORE_RADIUS,
  alpha = DIFF_ALPHA, steps = DIFF_STEPS)
seu$tumor_score <- tumor_score

score_diag <- attr(tumor_score, "diagnostics")
score_diag$condition <- as.character(
  seu$condition[match(score_diag$sample, as.character(seu$sample))])
write.csv(score_diag, file.path(OUTPUT_DIR, "score_diagnostics_per_sample.csv"),
          row.names = FALSE)

# ---- A6: does the COVARIATE itself separate by condition? -------------------
# If it does, s(score, by = cond_ord) cannot be read as a biological
# interaction, because the x-axis means something different in each arm.
cat("\n=== Covariate confounding check (A6) ===\n")
print(score_diag)
confound_flag <- FALSE
for (v in c("tumor_frac", "nn_dist_med", "score_median", "score_iqr")) {
  g <- split(score_diag[[v]], score_diag$condition)
  if (length(g) == 2 && all(lengths(g) >= 2)) {
    pw <- suppressWarnings(wilcox.test(g[[1]], g[[2]])$p.value)
    cat(sprintf("  %-14s  %s = %.3f | %s = %.3f   (Wilcoxon p = %.3f)\n",
                v, names(g)[1], mean(g[[1]]), names(g)[2], mean(g[[2]]), pw))
    if (is.finite(pw) && pw < 0.1) confound_flag <- TRUE
  }
}
if (confound_flag) {
  cat("\n  *** WARNING: a covariate property differs by condition at p < 0.1.\n",
      "  *** The score axis is not comparable across arms. Treat any\n",
      "  *** condition x trajectory result as confounded until resolved.\n", sep = "")
}

# ---- A6: permuted-tumor-label null check on the field -----------------------
cat("\n=== Score-field null check (permuted tumor labels) ===\n")
null_chk <- tumor_score_null_check(
  all_coords, is_tumor, samples,
  method = SCORE_METHOD, graph = SCORE_GRAPH, laplacian = SCORE_LAPLACIAN,
  sigma = SCORE_SIGMA, k = KNN_K, radius = SCORE_RADIUS,
  alpha = DIFF_ALPHA, steps = DIFF_STEPS)
saveRDS(null_chk[c("cor_real_null", "real_iqr", "null_iqr")],
        file.path(OUTPUT_DIR, "score_null_check.rds"))

SCORE_LAB  <- "Tumor proximity (weighted local tumor fraction)"
TRAJ_X_LAB <- SCORE_LAB   # A4: model x-axis and plot x-axis are the same thing now

# ---- Field + distribution plots ---------------------------------------------
sample_ids <- levels(seu$sample)
half_span <- max(max(tapply(all_coords[, 1], samples, function(x) diff(range(x)))),
                 max(tapply(all_coords[, 2], samples, function(x) diff(range(x))))) / 2 * 1.05
p_field <- wrap_plots(lapply(sample_ids, function(s) {
  idx <- which(samples == s)
  ggplot(data.frame(x = all_coords[idx, 1], y = all_coords[idx, 2],
                    score = tumor_score[idx]),
         aes(x, y, color = score)) +
    geom_point(size = 0.2, alpha = 0.5) +
    scale_color_viridis_c(option = "magma", limits = c(0, 1), name = "Tumor\nscore") +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    ggtitle(s) + theme_void() + theme(plot.title = element_text(hjust = .5, size = 10))
}), nrow = 1)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_score_field.pdf"),
       p_field, width = 4 * length(sample_ids), height = 5)

# A3: with no rescaling these histograms are directly comparable across mice.
p_hist <- ggplot(data.frame(score = tumor_score, sample = samples,
                            condition = seu$condition),
                 aes(x = score, fill = condition)) +
  geom_histogram(bins = 60, color = NA) +
  scale_fill_manual(values = CONDITION_COLORS) +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = SCORE_LAB, y = "Cell count",
       title = "Tumor score per sample (no per-sample rescaling)") +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "tumor_score_histogram_per_sample.pdf"),
       p_hist, width = 4 * ceiling(length(sample_ids) / 2), height = 6)

# =============================================================================
# Helpers
# =============================================================================

# D1: warn instead of silently switching to equal-width breaks.
make_bins <- function(score_vec, n_bins = N_BINS) {
  br <- quantile(score_vec, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
  br[1] <- br[1] - 1e-9
  if (length(unique(br)) < length(br)) {
    warning(sprintf(paste("make_bins: %d/%d quantile breaks are tied (heavily",
                          "zero-inflated score) -- falling back to EQUAL-WIDTH",
                          "bins. Bin midpoints are no longer quantiles."),
                    sum(duplicated(br)), length(br)), call. = FALSE)
    br <- seq(min(score_vec, na.rm = TRUE) - 1e-9, max(score_vec, na.rm = TRUE),
              length.out = n_bins + 1)
  }
  cat("    bin breaks:", paste(sprintf("%.4f", br), collapse = " "), "\n")
  list(bins = cut(score_vec, breaks = br, labels = seq_len(n_bins)), breaks = br)
}

# --- Pseudobulk: SUMMED RAW COUNTS + library size + REAL mean score (A5, B2) --
make_pseudobulk <- function(count_mat, score_vec, sample_vec, condition_vec,
                            n_bins = N_BINS, min_cells = MIN_CELLS_PER_BIN) {
  b <- make_bins(score_vec, n_bins)
  bins <- as.integer(b$bins)
  groups <- paste(sample_vec, bins, sep = "__")

  keep <- names(which(table(groups) >= min_cells))
  keep <- keep[!grepl("NA$", keep)]
  if (length(keep) == 0) return(NULL)

  expr <- vapply(keep, function(g) {
    Matrix::rowSums(count_mat[, groups == g, drop = FALSE])   # B2: SUM counts
  }, numeric(nrow(count_mat)))
  rownames(expr) <- rownames(count_mat)

  meta <- bind_rows(lapply(keep, function(g) {
    ci <- which(groups == g)
    data.frame(group = g,
               sample = strsplit(g, "__")[[1]][1],
               condition = as.character(condition_vec[ci[1]]),
               bin = as.integer(strsplit(g, "__")[[1]][2]),
               score = mean(score_vec[ci]),               # A5: real mean score
               n_cells = length(ci),
               libsize = sum(count_mat[, ci]),            # B2: offset
               stringsAsFactors = FALSE)
  }))
  meta$condition <- factor(meta$condition,
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION))

  # --- B5: require BOTH conditions to contribute >= MIN_SAMPLES_PER_BIN -------
  cov <- meta %>% group_by(bin, condition) %>%
    summarise(n_samp = n_distinct(sample), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = condition, values_from = n_samp, values_fill = 0)
  ok_bins <- cov$bin[apply(cov[, -1, drop = FALSE], 1,
                           function(r) all(r >= MIN_SAMPLES_PER_BIN))]
  dropped <- setdiff(sort(unique(meta$bin)), ok_bins)
  if (length(dropped)) {
    cat(sprintf("    B5: dropping bins %s (a condition has < %d samples there)\n",
                paste(dropped, collapse = ","), MIN_SAMPLES_PER_BIN))
  }
  if (length(ok_bins) < 4) return(NULL)
  sel <- meta$bin %in% ok_bins
  list(expr = expr[, sel, drop = FALSE], meta = meta[sel, ],
       bin_coverage = cov, dropped_bins = dropped)
}

# --- GAM fit (B1, B2, B3, D3) -------------------------------------------------
# stat_only = TRUE skips prediction; used inside the permutation loop.
fit_trajectory_gam <- function(df, model, gam_k = GAM_K, stat_only = FALSE) {
  df$sample <- droplevels(factor(df$sample))
  df$cond_ord <- ordered(as.character(df$condition),
                         levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  if (nlevels(droplevels(factor(df$condition))) < 2 || nlevels(df$sample) < 2)
    return(list(ok = FALSE))

  if (model == "fs") {
    # B1: per-sample random SMOOTH. fs carries its own intercept, so the
    # redundant s(sample, bs="re") of the old code is deliberately absent.
    rhs  <- ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
              cond_ord + s(score, sample, bs = "fs", k = gam_k)
    excl <- "s(score,sample)"
  } else if (model == "ranef") {
    # Retained ONLY as the anticonservative contrast: intercept-only.
    rhs  <- ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
              cond_ord + s(sample, bs = "re")
    excl <- "s(sample)"
  } else {
    stop("fit_trajectory_gam: unrecognized model '", model, "'")   # D3
  }

  fam     <- if (GAM_FAMILY == "nb") nb() else gaussian()
  rhs_txt <- paste(deparse(rhs[[2]]), collapse = " ")
  off_txt <- if (GAM_FAMILY == "nb") " + offset(log(libsize))" else ""
  form    <- as.formula(paste0("expr ~ ", rhs_txt, off_txt))

  fit <- tryCatch(bam(form, data = df, family = fam, method = "fREML",   # D5
                      discrete = TRUE, nthreads = 1),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE))

  sm <- summary(fit)
  # B3: index BY POSITION. Under nb() the columns are Chi.sq / Pr(>|z|),
  # under gaussian() they are F / Pr(>|t|); name-based lookup returns NA.
  s_rows <- rownames(sm$s.table)
  dr <- s_rows[grepl(":cond_ord", s_rows, fixed = TRUE)]
  p_shape <- if (length(dr)) sm$s.table[dr[1], ncol(sm$s.table)] else NA_real_
  p_rows <- rownames(sm$p.table)
  lr <- p_rows[grepl("cond_ord", p_rows, fixed = TRUE)]
  p_level <- if (length(lr)) sm$p.table[lr[1], ncol(sm$p.table)] else NA_real_

  if (stat_only) return(list(ok = TRUE, p_shape = p_shape, p_level = p_level))

  # B4: report basis adequacy so k=10 can be justified or raised.
  kc <- tryCatch(k.check(fit), error = function(e) NULL)

  newd <- expand.grid(score = seq(min(df$score), max(df$score), length.out = 100),
                      cond_ord = ordered(c(REFERENCE_CONDITION, TEST_CONDITION),
                                         levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
  newd$sample  <- df$sample[1]
  newd$libsize <- exp(mean(log(df$libsize)))   # predict at geometric-mean depth
  pr <- predict(fit, newdata = newd, se.fit = TRUE, exclude = excl,
                newdata.guaranteed = TRUE)
  newd$condition <- factor(as.character(newd$cond_ord),
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  newd$fit <- pr$fit; newd$se <- pr$se.fit
  newd$ci_lo <- pr$fit - 1.96 * pr$se.fit
  newd$ci_hi <- pr$fit + 1.96 * pr$se.fit

  list(pred = newd, p_shape = p_shape, p_level = p_level, k_check = kc, ok = TRUE)
}

# =============================================================================
# Per cell type
# =============================================================================
for (ct in CELLTYPES_TRAJECTORY) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
  ct_dir  <- file.path(OUTPUT_DIR, ct_safe)
  diag_dir <- file.path(ct_dir, "diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n=== %s ===\n", ct))

  ct_idx <- which(celltypes == ct & is.finite(tumor_score) & qc_pass)
  if (length(ct_idx) < 50) { cat("  too few cells; skipping\n"); next }

  ct_score     <- tumor_score[ct_idx]
  ct_condition <- droplevels(seu$condition[ct_idx])
  ct_sample    <- droplevels(seu$sample[ct_idx])
  count_mat    <- GetAssayData(seu, assay = "RNA", layer = "counts")[, ct_idx]
  norm_mat     <- GetAssayData(seu, assay = "RNA", layer = "data")[, ct_idx]
  ct_libsize   <- Matrix::colSums(count_mat)
  # A4: decile coordinate is computed for PLOTTING ONLY and never modelled.
  ct_score_dec <- ecdf(ct_score)(ct_score)

  cat(sprintf("  cells: %d (WT %d, DKO %d) across %d samples\n", length(ct_idx),
              sum(ct_condition == REFERENCE_CONDITION),
              sum(ct_condition == TEST_CONDITION), nlevels(ct_sample)))

  # --- C1: expression-only gene filter ---------------------------------------
  pct_det   <- Matrix::rowMeans(count_mat > 0)
  mean_expr <- Matrix::rowMeans(norm_mat)
  test_genes <- rownames(count_mat)[mean_expr > MIN_MEAN_EXPR &
                                    pct_det   > MIN_PCT_DETECTED]
  cat(sprintf("  C1: testing %d genes (expression filter only; FDR denominator = %d)\n",
              length(test_genes), length(test_genes)))
  if (length(test_genes) < 10) { cat("  too few genes; skipping\n"); next }

  # --- Condition-wise cell density along the gradient (read this FIRST) -------
  # A4: if these differ, a compositional shift can masquerade as a shape change.
  ks <- suppressWarnings(ks.test(ct_score[ct_condition == REFERENCE_CONDITION],
                                 ct_score[ct_condition == TEST_CONDITION]))
  cat(sprintf("  A4: KS test of CD8 score distribution WT vs DKO: D=%.3f p=%.3g\n",
              ks$statistic, ks$p.value))
  p_dist <- ggplot(data.frame(score = ct_score, condition = ct_condition,
                              sample = ct_sample),
                   aes(x = score, color = condition, fill = condition)) +
    geom_density(alpha = 0.3) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(x = SCORE_LAB, y = "Density of cells",
         title = sprintf("%s — position along tumor gradient", ct),
         subtitle = sprintf("KS D = %.3f, p = %.3g — a shift here can masquerade as a shape change",
                            ks$statistic, ks$p.value)) +
    theme_publication()
  ggsave(file.path(ct_dir, "density_distribution_along_gradient.pdf"),
         p_dist, width = 6, height = 4)

  # --- Pseudobulk ------------------------------------------------------------
  pb <- make_pseudobulk(count_mat[test_genes, , drop = FALSE], ct_score,
                        ct_sample, ct_condition)
  if (is.null(pb)) { cat("  pseudobulk failed coverage requirements; skipping\n"); next }
  cat(sprintf("  pseudobulk: %d sample x bin groups, score range %.3f-%.3f\n",
              nrow(pb$meta), min(pb$meta$score), max(pb$meta$score)))
  write.csv(pb$bin_coverage, file.path(diag_dir, "bin_coverage.csv"), row.names = FALSE)

  p_tile <- ggplot(pb$meta, aes(x = factor(bin), y = sample, fill = n_cells)) +
    geom_tile(color = "white") + geom_text(aes(label = n_cells), size = 2.6, color = "white") +
    scale_fill_viridis_c(option = "mako", name = "n cells") +
    labs(x = "Score bin (retained only)", y = NULL,
         title = sprintf("%s — cells per sample x bin after B5 coverage filter", ct)) +
    theme_publication()
  ggsave(file.path(diag_dir, "cells_per_sample_bin.pdf"), p_tile, width = 7, height = 4)

  # ===========================================================================
  # GAM fits
  # ===========================================================================
  results <- list()
  for (model in GAM_MODELS) {
    cat(sprintf("\n  --- model: %s%s ---\n", model,
                if (model == PRIMARY_MODEL) "  [PRIMARY]" else "  [diagnostic contrast]"))
    gam_dir <- file.path(ct_dir, paste0("gam_", model))
    dir.create(gam_dir, recursive = TRUE, showWarnings = FALSE)

    build_df <- function(gene, cond_map = NULL) {
      cond <- pb$meta$condition
      if (!is.null(cond_map)) cond <- factor(cond_map[pb$meta$sample],
                                             levels = c(REFERENCE_CONDITION, TEST_CONDITION))
      data.frame(expr = as.numeric(pb$expr[gene, ]),
                 score = pb$meta$score,          # A5
                 condition = cond,
                 sample = pb$meta$sample,
                 libsize = pb$meta$libsize)      # B2
    }

    fits <- mclapply(test_genes, function(g) {
      r <- fit_trajectory_gam(build_df(g), model)
      if (!isTRUE(r$ok)) return(NULL)
      list(gene = g, p_shape = r$p_shape, p_level = r$p_level,
           pred = r$pred[, c("score", "condition", "fit")], k_check = r$k_check)
    }, mc.cores = N_PERM_CORES)
    names(fits) <- test_genes
    fits <- Filter(Negate(is.null), fits)
    if (length(fits) < 10) { cat("    too few converged fits\n"); next }

    gs <- bind_rows(lapply(fits, function(f) data.frame(gene = f$gene,
                                                        p_shape = f$p_shape,
                                                        p_level = f$p_level)))
    gs$fdr_shape <- p.adjust(gs$p_shape, method = "BH")   # C1: over the FULL set
    gs$fdr_level <- p.adjust(gs$p_level, method = "BH")

    # B4: basis adequacy across genes
    edf_glob <- sapply(fits, function(f) if (is.null(f$k_check)) NA else f$k_check[1, "edf"])
    cat(sprintf("    B4: global smooth edf — median %.2f, 95th pct %.2f (k = %d, ceiling %d)\n",
                median(edf_glob, na.rm = TRUE), quantile(edf_glob, .95, na.rm = TRUE),
                GAM_K, GAM_K - 1))
    if (quantile(edf_glob, .95, na.rm = TRUE) > GAM_K - 1.5)
      cat("    *** edf is saturating k — raise GAM_K.\n")

    # =========================================================================
    # C2: sample-level permutation
    # =========================================================================
    if (RUN_PERMUTATION && model == PRIMARY_MODEL) {
      samp_cond <- pb$meta %>% distinct(sample, condition)
      samps <- samp_cond$sample
      n_ref <- sum(samp_cond$condition == REFERENCE_CONDITION)
      assigns <- combn(length(samps), n_ref, simplify = FALSE)
      cat(sprintf("\n    C2: %d distinct sample-label assignments (choose(%d,%d)).\n",
                  length(assigns), length(samps), n_ref))
      cat(sprintf("        Smallest achievable permutation p = %.4f.\n", 1 / length(assigns)))
      cat(sprintf("        Smallest achievable BH-FDR over %d genes = %.2f",
                  nrow(gs), min(1, nrow(gs) / length(assigns))))
      cat(" -> per-gene FDR control is NOT attainable at this n.\n")

      true_set <- which(samp_cond$condition == REFERENCE_CONDITION)
      perm_stats <- vector("list", length(assigns))
      obs_idx <- NA_integer_
      for (a in seq_along(assigns)) {
        if (setequal(assigns[[a]], true_set)) obs_idx <- a
        cm <- setNames(rep(TEST_CONDITION, length(samps)), as.character(samps))
        cm[assigns[[a]]] <- REFERENCE_CONDITION
        ps <- unlist(mclapply(names(fits), function(g) {
          r <- fit_trajectory_gam(build_df(g, cm), model, stat_only = TRUE)
          if (isTRUE(r$ok)) r$p_shape else NA_real_
        }, mc.cores = N_PERM_CORES))
        perm_stats[[a]] <- -log10(pmax(ps, 1e-300))
        cat(sprintf("      assignment %2d/%d: %d genes at BH-FDR < %.2f%s\n",
                    a, length(assigns),
                    sum(p.adjust(ps, "BH") < SIG_FDR, na.rm = TRUE), SIG_FDR,
                    if (a == obs_idx) "   <- OBSERVED" else ""))
      }

      # (a) GLOBAL test
      counts <- sapply(perm_stats, function(v) sum(p.adjust(10^(-v), "BH") < SIG_FDR, na.rm = TRUE))
      p_global <- mean(counts >= counts[obs_idx])
      cat(sprintf("\n    C2(a) GLOBAL permutation p = %.4f  (observed %d genes vs null median %d)\n",
                  p_global, counts[obs_idx], median(counts[-obs_idx])))

      # (b) EMPIRICAL per-gene p against the pooled null (identity excluded)
      null_pool <- unlist(perm_stats[-obs_idx])
      null_pool <- null_pool[is.finite(null_pool)]
      obs_stat  <- -log10(pmax(gs$p_shape, 1e-300))
      gs$p_shape_emp <- (1 + sapply(obs_stat, function(o) sum(null_pool >= o))) /
                        (1 + length(null_pool))
      gs$fdr_shape_emp <- p.adjust(gs$p_shape_emp, method = "BH")
      cat(sprintf("    C2(b) empirical p: %d genes < 0.05, %d at empirical BH-FDR < %.2f\n",
                  sum(gs$p_shape_emp < 0.05), sum(gs$fdr_shape_emp < SIG_FDR), SIG_FDR))

      saveRDS(list(counts = counts, obs_idx = obs_idx, p_global = p_global,
                   null_pool = null_pool),
              file.path(ct_dir, "permutation_null.rds"))

      p_perm <- ggplot(data.frame(n = counts), aes(x = n)) +
        geom_histogram(bins = 20, fill = "grey70") +
        geom_vline(xintercept = counts[obs_idx], color = "#D55E00", linewidth = 1) +
        labs(x = sprintf("# genes at BH-FDR < %.2f", SIG_FDR), y = "Assignments",
             title = sprintf("%s — sample-label permutation null", ct),
             subtitle = sprintf("observed = %d (red), global p = %.4f",
                                counts[obs_idx], p_global)) +
        theme_publication()
      ggsave(file.path(ct_dir, "permutation_global_test.pdf"), p_perm, width = 6, height = 4)

      attr(gs, "p_global") <- p_global
    }

    gs <- gs[order(gs$p_shape), ]
    write.csv(gs, file.path(ct_dir, sprintf("gam_condition_test_%s.csv", model)),
              row.names = FALSE)
    results[[model]] <- list(gs = gs, fits = fits)

    # =========================================================================
    # Significant-trajectory heatmap (C3: labelled as ranking where appropriate)
    # =========================================================================
    use_emp <- "fdr_shape_emp" %in% names(gs)
    sig_col <- if (use_emp) "fdr_shape_emp" else "fdr_shape"
    sig <- gs[which(gs[[sig_col]] < SIG_FDR | gs$fdr_level < SIG_FDR), ]
    if (nrow(sig) < 3) {
      cat(sprintf("    < 3 genes pass; showing top 30 by shape rank instead\n"))
      sig <- head(gs, 30)
      ranking_only <- TRUE
    } else ranking_only <- !use_emp

    genes_keep <- sig$gene
    ng <- sum(fits[[genes_keep[1]]]$pred$condition == REFERENCE_CONDITION)
    getc <- function(g, cond) { d <- fits[[g]]$pred; d <- d[d$condition == cond, ]
                                d$fit[order(d$score)] }
    mat_wt  <- t(vapply(genes_keep, getc, numeric(ng), cond = REFERENCE_CONDITION))
    mat_dko <- t(vapply(genes_keep, getc, numeric(ng), cond = TEST_CONDITION))
    rownames(mat_wt) <- rownames(mat_dko) <- genes_keep

    comb <- cbind(mat_wt, mat_dko)
    rmn <- rowMeans(comb, na.rm = TRUE)
    rsd <- apply(comb, 1, sd, na.rm = TRUE); rsd[rsd == 0 | is.na(rsd)] <- 1
    z_wt <- (mat_wt - rmn) / rsd; z_dko <- (mat_dko - rmn) / rsd

    cl_rows <- if (HEATMAP_ROW_ORDER == "cluster")
      as.dendrogram(hclust(dist(cbind(z_wt, z_dko)))) else FALSE

    nlog <- function(p) pmin(-log10(pmax(p, 1e-300)), 10)
    ann_lab <- if (use_emp) "-log10 empirical FDR" else "-log10 FDR (RANK ONLY)"
    row_ha <- rowAnnotation(
      shape = anno_barplot(nlog(sig[[sig_col]]), gp = gpar(fill = "#D55E00"), border = FALSE),
      level = anno_barplot(nlog(sig$fdr_level),  gp = gpar(fill = "#0072B2"), border = FALSE),
      annotation_label = c(ann_lab, "-log10 FDR level"),
      annotation_name_gp = gpar(fontsize = 7), width = unit(3, "cm"))
    col_note <- HeatmapAnnotation(
      axis = anno_text(c("low tumor", rep("", ng - 2), "high tumor"),
                       gp = gpar(fontsize = 7), just = "center"),
      show_annotation_name = FALSE)
    cf <- colorRamp2(c(-2, -1, 0, 1, 2), viridis::inferno(5))

    subtitle <- if (!is.null(attr(gs, "p_global")))
      sprintf("global permutation p = %.4f", attr(gs, "p_global")) else
      "no permutation run"
    if (ranking_only) subtitle <- paste(subtitle,
      "| per-gene values are a RANKING, not calibrated FDR")

    ht1 <- Heatmap(z_wt, cluster_rows = cl_rows, cluster_columns = FALSE, col = cf,
                   name = "Fitted\nz-score", column_title = REFERENCE_CONDITION,
                   show_column_names = FALSE, row_names_side = "left",
                   row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                   bottom_annotation = col_note, width = unit(6.5, "cm"), border = TRUE)
    ht2 <- Heatmap(z_dko, cluster_rows = FALSE, cluster_columns = FALSE, col = cf,
                   name = " ", column_title = TEST_CONDITION,
                   show_column_names = FALSE, show_row_names = FALSE,
                   bottom_annotation = col_note, width = unit(6.5, "cm"),
                   border = TRUE, show_heatmap_legend = FALSE)
    pdf(file.path(ct_dir, sprintf("trajectory_heatmap_%s.pdf", model)),
        width = 13, height = max(3.5, nrow(z_wt) * 0.13 + 2))
    draw(ht1 + ht2 + row_ha,
         column_title = sprintf("%s — tumor-proximity trajectories [%s]\n%s",
                                ct, model, subtitle),
         column_title_gp = gpar(fontsize = 11, fontface = "bold"))
    dev.off()
  }

  # --- B1: how much did the random-intercept spec overstate significance? -----
  if (all(c("fs", "ranef") %in% names(results))) {
    m <- merge(results$fs$gs[, c("gene", "p_shape")],
               results$ranef$gs[, c("gene", "p_shape")],
               by = "gene", suffixes = c("_fs", "_ranef"))
    cat(sprintf("\n  B1 contrast: BH-FDR < %.2f — fs: %d genes, ranef: %d genes\n",
                SIG_FDR,
                sum(results$fs$gs$fdr_shape < SIG_FDR),
                sum(results$ranef$gs$fdr_shape < SIG_FDR)))
    p_cmp <- ggplot(m, aes(x = -log10(p_shape_ranef), y = -log10(p_shape_fs))) +
      geom_point(alpha = 0.3, size = 0.8) + geom_abline(slope = 1, linetype = "dashed") +
      labs(x = "-log10 p  (random intercept)", y = "-log10 p  (random smooth)",
           title = sprintf("%s — points below the line are pseudoreplication", ct)) +
      theme_publication()
    ggsave(file.path(diag_dir, "model_spec_comparison.pdf"), p_cmp, width = 5, height = 5)
  }

  saveRDS(results, file.path(ct_dir, "trajectory_results.rds"))
  cat(sprintf("  Done: %s\n", ct))
}

cat("\nAll done. Outputs ->", OUTPUT_DIR, "\n")