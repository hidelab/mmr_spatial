# 15_07_tumor_trajectory_tradeseq.R
# tradeSeq NB-GAM cross-check of the tumor-proximity trajectory (CD8+ T).
#
# WHY THIS SCRIPT EXISTS
# ----------------------
# 15_04 hand-rolls an mgcv difference-smooth test of whether CD8+ T expression
# along the tumor-proximity gradient differs between WT and DKO. tradeSeq
# (Van den Berge et al. 2020, Nat Commun) is the STANDARD, CITABLE implementation
# of exactly that core: per-gene negative-binomial GAMs along a continuous
# pseudotime axis, with
#   * associationTest() — does a gene vary along the axis?
#   * conditionTest()   — do the WT vs DKO smoothers differ? (the p_shape question)
# This script runs tradeSeq on the *identical* diffusion score used by 15_04, so
# the two analyses share one x-axis: same axis, different statistical engine =
# clean orthogonal cross-check (see trajectory.md "orthogonal published method").
# It also uses a proper NB likelihood + library-size offset on RAW COUNTS, which
# is the more defensible likelihood than 15_04's Gaussian-on-log-data (reviewer
# point #3-A in notes.md).
#
# ============================================================================
# READ BEFORE INTERPRETING ANY P-VALUE HERE
# ============================================================================
# tradeSeq's conditionTest is a CELL-LEVEL test with NO biological-replicate
# (mouse / sample) random effect, and tradeSeq does NO common-support filtering.
# With only 3 WT + 4 DKO mice this is precisely the ANTICONSERVATIVE
# pseudoreplication regime that notes.md diagnosed and rejected:
#   * the honest sample-level permutation floor is ~1/35 = 0.029;
#   * the defensible spec is 15_04's pseudobulk_fs (per-sample factor smooth);
#   * the global sample-level permutation p for a shape difference was ~0.54.
# Therefore the association/condition p-values below are EXPLORATORY — a per-gene
# RANKING, NOT honest inference. Expect FAR MORE "significant" genes here than
# pseudobulk_fs reports; that inflation is the cell-level artifact, not biology.
# Read these results ALONGSIDE the pseudobulk_fs conclusion in notes.md, never as
# a replacement for it. Condition smoothers in WT-only high-score regions are
# EXTRAPOLATED (no common support), so trust the periphery->mid range only.
#
# Outputs -> output/R/15_tumor_trajectory_tradeseq.

# --- Thread limits (set BEFORE BLAS/library load) ----------------------------
# tradeSeq parallelizes over GENES with BiocParallel; each per-gene mgcv fit is
# tiny, so we want MANY workers with FEW BLAS threads each (avoid nested
# oversubscription). N_WORKERS * N_BLAS_THREADS ~= half the 48-core box.
N_WORKERS      <- 12
N_BLAS_THREADS <- 2
Sys.setenv(OMP_NUM_THREADS = N_BLAS_THREADS, MKL_NUM_THREADS = N_BLAS_THREADS,
           OPENBLAS_NUM_THREADS = N_BLAS_THREADS,
           VECLIB_MAXIMUM_THREADS = N_BLAS_THREADS,
           NUMEXPR_NUM_THREADS = N_BLAS_THREADS)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(N_BLAS_THREADS)
  RhpcBLASctl::omp_set_num_threads(N_BLAS_THREADS)
}

# --- Dependency guard --------------------------------------------------------
# tradeSeq is NOT installed by default in this env; install it yourself, then run.
if (!requireNamespace("tradeSeq", quietly = TRUE)) {
  stop("tradeSeq is not installed. Install it (Bioconductor), then re-run:\n",
       "  if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')\n",
       "  BiocManager::install('tradeSeq')\n",
       call. = FALSE)
}

library(tradeSeq)
library(SingleCellExperiment)
library(BiocParallel)
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

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "15_tumor_trajectory_tradeseq")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

BPPARAM <- BiocParallel::MulticoreParam(workers = N_WORKERS, progressbar = TRUE)

# =============================================================================
# Parameters
# =============================================================================
# --- Tumor score (IDENTICAL to 15_04's diffusion settings, for comparability) -
SCORE_METHOD         <- "diffusion"   # keep "diffusion" to match 15_04
DIFFUSION_VARIANT    <- "heat"
DIFFUSION_HEAT_STEPS <- 20
DIFFUSION_HEAT_ALPHA <- 2
DIFFUSION_LOCAL_FRACTION <- TRUE       # kernel-weighted local tumor fraction, [0,1]
KNN_K                <- 30
KDE_BW_SCALE         <- 1.0
KDE_GRID_N           <- 200
TARGET_CELLTYPE      <- "Tumor"

# --- tradeSeq ----------------------------------------------------------------
NKNOTS        <- 6         # spline knots per smoother (tradeSeq default)
REUSE_SCE     <- TRUE      # TRUE -> reload a cached tradeseq_sce.rds if present
                           #         (skip the expensive refit; iterate on tests/plots)
RUN_EVALUATEK <- FALSE     # TRUE -> emit an evaluateK() knot-selection diagnostic
EVALK_KRANGE  <- 3:10      # candidate knot counts for evaluateK
EVALK_NGENES  <- 100       # genes sampled by evaluateK

# --- Genes / cell types ------------------------------------------------------
# Test-gene set uses the SAME expression filter as 15_04 / the DESeq2 DE scripts
# (00_parameters.R): detected in >= MIN_GENE_EXPR_FRAC of this cell type's cells
# AND >= MIN_COUNTS_PER_SAMPLE counts in >= MIN_SAMPLES_EXPRESSED samples.
CELLTYPES_TRAJECTORY <- c("CD8+ T")

# --- Visualization -----------------------------------------------------------
TOP_N_PLOT <- 25    # per-gene smoother PDFs + max rows in the heatmap
SIG_FDR    <- 0.1   # BH-FDR cutoff for the condition heatmap gene set
HEATMAP_ROW_ORDER <- "cluster"   # "cluster" | "shape" (by condition p-value)

# --- Quality control ---------------------------------------------------------
# Match 15_04: QC filter is intentionally OFF (the tumor field is a property of
# tissue architecture; every physical cell shapes it). Set TRUE to drop sparse
# cells from the trajectory only (never from the score field).
APPLY_QC <- FALSE

stopifnot(SCORE_METHOD %in% c("kde", "knn", "diffusion"))
stopifnot(HEATMAP_ROW_ORDER %in% c("cluster", "shape"))

# =============================================================================
# Parameter report (written FIRST, before any heavy computation)
# =============================================================================
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
        script     = "code/R/15_07_tumor_trajectory_tradeseq.R",
        run_time   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        host       = as.character(Sys.info()[["nodename"]]),
        R_version  = R.version.string,
        tradeSeq   = pkgver("tradeSeq"),
        mgcv       = pkgver("mgcv"),
        Seurat     = pkgver("Seurat"),
        ComplexHeatmap = pkgver("ComplexHeatmap"),
        N_WORKERS  = N_WORKERS, N_BLAS_THREADS = N_BLAS_THREADS,
        input_rds  = "output/R/DKO_banksy_merged.rds"),
    grp("tumor_score",
        SCORE_METHOD = SCORE_METHOD, DIFFUSION_VARIANT = DIFFUSION_VARIANT,
        DIFFUSION_HEAT_STEPS = DIFFUSION_HEAT_STEPS,
        DIFFUSION_HEAT_ALPHA = DIFFUSION_HEAT_ALPHA,
        DIFFUSION_LOCAL_FRACTION = DIFFUSION_LOCAL_FRACTION,
        KNN_K = KNN_K, TARGET_CELLTYPE = TARGET_CELLTYPE),
    grp("tradeseq",
        NKNOTS = NKNOTS, RUN_EVALUATEK = RUN_EVALUATEK,
        EVALK_KRANGE = paste(range(EVALK_KRANGE), collapse = "-"),
        EVALK_NGENES = EVALK_NGENES),
    grp("genes",
        MIN_GENE_EXPR_FRAC = MIN_GENE_EXPR_FRAC,
        MIN_COUNTS_PER_SAMPLE = MIN_COUNTS_PER_SAMPLE,
        MIN_SAMPLES_EXPRESSED = MIN_SAMPLES_EXPRESSED),
    grp("viz",
        TOP_N_PLOT = TOP_N_PLOT, SIG_FDR = SIG_FDR,
        HEATMAP_ROW_ORDER = HEATMAP_ROW_ORDER),
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
               "15_07 tumor-trajectory tradeSeq — run parameters",
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
# Compute tumor score (per sample) via 00_utils.R — SAME axis as 15_04
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
                                  normalize = !DIFFUSION_LOCAL_FRACTION,
                                  alpha = DIFFUSION_HEAT_ALPHA)
)
seu$tumor_score <- tumor_score
cat("Tumor score computed.\n")

SCORE_LAB <- sprintf("Tumor Proximity Score (%s)", toupper(SCORE_METHOD))

# =============================================================================
# Global: spatial score field + per-sample histogram (confirms same axis as 15_04)
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
# Helper: build a WT|DKO fitted-smoother matrix from predictSmooth() output
# =============================================================================
# predictSmooth(tidy=TRUE) returns a long data frame; column names have varied a
# little across tradeSeq versions, so resolve them defensively.
extract_smooth_matrices <- function(sce, genes, npoints = 100) {
  ps <- tradeSeq::predictSmooth(sce, gene = genes, nPoints = npoints, tidy = TRUE)
  ps <- as.data.frame(ps)
  cn <- tolower(colnames(ps))
  gcol <- colnames(ps)[which(cn %in% c("gene", "genes"))[1]]
  tcol <- colnames(ps)[which(cn %in% c("time", "pseudotime", "t"))[1]]
  ycol <- colnames(ps)[which(cn %in% c("yhat", "y", "fitted"))[1]]
  ccol <- colnames(ps)[which(cn == "condition")[1]]
  stopifnot(!is.na(gcol), !is.na(tcol), !is.na(ycol), !is.na(ccol))

  build_one <- function(cond_label) {
    sub <- ps[as.character(ps[[ccol]]) == cond_label, ]
    w <- sub %>%
      dplyr::select(dplyr::all_of(c(gcol, tcol, ycol))) %>%
      tidyr::pivot_wider(names_from = dplyr::all_of(tcol),
                         values_from = dplyr::all_of(ycol))
    m <- as.matrix(w[, -1, drop = FALSE])
    rownames(m) <- w[[gcol]]
    # order columns by numeric time
    m[, order(as.numeric(colnames(m))), drop = FALSE][genes, , drop = FALSE]
  }
  # tradeSeq labels conditions by their factor levels
  list(wt  = build_one(REFERENCE_CONDITION),
       dko = build_one(TEST_CONDITION))
}

# =============================================================================
# Per cell type
# =============================================================================
for (ct in CELLTYPES_TRAJECTORY) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
  ct_dir  <- file.path(OUTPUT_DIR, ct_safe)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n=== %s ===\n", ct))

  ct_idx <- which(celltypes == ct & !is.na(tumor_score) & qc_pass)
  if (length(ct_idx) < 50) {
    cat(sprintf("  Skipping %s: too few cells (%d)\n", ct, length(ct_idx)))
    next
  }

  ct_score     <- tumor_score[ct_idx]
  ct_condition <- droplevels(seu$condition[ct_idx])
  ct_sample    <- droplevels(seu$sample[ct_idx])

  cat(sprintf("  Total cells: %d (WT: %d, DKO: %d)\n", length(ct_idx),
              sum(ct_condition == REFERENCE_CONDITION),
              sum(ct_condition == TEST_CONDITION)))
  if (nlevels(ct_condition) < 2) {
    cat("  Skipping: need both conditions present.\n"); next
  }

  # --- Test-gene set: same expression filter as 15_04 / DESeq2 ---------------
  counts_ct <- GetAssayData(seu, assay = "RNA", layer = "counts")[, ct_idx]
  det_frac  <- Matrix::rowSums(counts_ct > 0) / length(ct_idx)
  samp_sum  <- vapply(levels(ct_sample), function(s)
                 Matrix::rowSums(counts_ct[, ct_sample == s, drop = FALSE]),
                 numeric(nrow(counts_ct)))
  n_samp_det <- rowSums(samp_sum >= MIN_COUNTS_PER_SAMPLE)
  test_genes <- rownames(counts_ct)[det_frac >= MIN_GENE_EXPR_FRAC &
                                    n_samp_det >= MIN_SAMPLES_EXPRESSED]
  if (length(test_genes) < 5) {
    cat(sprintf("  Skipping %s: too few genes (%d)\n", ct, length(test_genes))); next
  }
  cat(sprintf("  Genes tested: %d (>=%.0f%% cells & >=%d counts in >=%d samples)\n",
              length(test_genes), MIN_GENE_EXPR_FRAC * 100, MIN_COUNTS_PER_SAMPLE,
              MIN_SAMPLES_EXPRESSED))

  # tradeSeq wants a genes x cells RAW COUNT matrix (integer), NB likelihood.
  counts_mat <- as.matrix(counts_ct[test_genes, , drop = FALSE])

  # --- tradeSeq inputs: single spatial "lineage" = the tumor-proximity score --
  pseudotime  <- matrix(ct_score, ncol = 1,
                        dimnames = list(colnames(counts_mat), "Lineage1"))
  cellWeights <- matrix(1, nrow = length(ct_idx), ncol = 1,
                        dimnames = list(colnames(counts_mat), "Lineage1"))
  conditions  <- ct_condition   # factor: WT (ref) / DKO

  # --- Optional knot-selection diagnostic ------------------------------------
  if (RUN_EVALUATEK) {
    cat("  Running evaluateK() knot diagnostic...\n")
    set.seed(1)
    pdf(file.path(ct_dir, "evaluateK_diagnostic.pdf"), width = 10, height = 6)
    ek <- tradeSeq::evaluateK(counts = counts_mat, pseudotime = pseudotime,
                              cellWeights = cellWeights, conditions = conditions,
                              k = EVALK_KRANGE, nGenes = EVALK_NGENES,
                              verbose = TRUE, parallel = TRUE, BPPARAM = BPPARAM)
    dev.off()
  }

  # --- Fit NB-GAMs (condition-specific smoothers on the shared knots) --------
  # REUSE_SCE reloads a cached fit so tests/plots can be iterated without the
  # (expensive) refit; the cache is only reused when its gene set still matches.
  sce_path <- file.path(ct_dir, "tradeseq_sce.rds")
  sce <- NULL
  if (REUSE_SCE && file.exists(sce_path)) {
    cat("  REUSE_SCE: loading cached fit ->", sce_path, "\n")
    sce_cached <- readRDS(sce_path)
    if (setequal(rownames(sce_cached), rownames(counts_mat))) {
      sce <- sce_cached
    } else {
      cat("  cached gene set differs from current test_genes; refitting.\n")
    }
  }
  if (is.null(sce)) {
    cat(sprintf("  Fitting NB-GAMs with tradeSeq (nknots = %d, %d workers)...\n",
                NKNOTS, N_WORKERS))
    set.seed(1)
    sce <- tradeSeq::fitGAM(counts = counts_mat, pseudotime = pseudotime,
                            cellWeights = cellWeights, conditions = conditions,
                            nknots = NKNOTS, parallel = TRUE, BPPARAM = BPPARAM,
                            verbose = TRUE)
    saveRDS(sce, sce_path)   # cache immediately so a later crash keeps the fit
  }
  conv <- SummarizedExperiment::rowData(sce)$tradeSeq$converged
  cat(sprintf("  fitGAM: %d/%d genes converged\n",
              sum(conv, na.rm = TRUE), length(conv)))

  # --- associationTest: does expression vary along the tumor gradient? -------
  assoc <- tradeSeq::associationTest(sce, global = TRUE, lineages = FALSE)
  assoc$gene <- rownames(assoc)
  assoc$fdr  <- p.adjust(assoc$pvalue, method = "BH")
  # Rank by p, then by DESCENDING waldStat: many p-values underflow to exactly 0,
  # so the Wald statistic is what actually orders the top of the (exploratory)
  # ranking — otherwise ties fall back to alphabetical gene order.
  assoc <- assoc[order(assoc$pvalue, -assoc$waldStat),
                 c("gene", setdiff(colnames(assoc), "gene"))]
  write.csv(assoc, file.path(ct_dir, "tradeseq_association_test.csv"), row.names = FALSE)
  cat(sprintf("  associationTest: %d genes, %d at FDR < %.2f (EXPLORATORY, cell-level)\n",
              nrow(assoc), sum(assoc$fdr < SIG_FDR, na.rm = TRUE), SIG_FDR))

  # --- conditionTest: do the WT vs DKO smoothers differ? (p_shape analog) ----
  cond <- tradeSeq::conditionTest(sce, global = TRUE, pairwise = FALSE)
  cond$gene <- rownames(cond)
  cond$fdr  <- p.adjust(cond$pvalue, method = "BH")
  # Same tie-break as associationTest: p, then descending waldStat (many p == 0),
  # so the top-ranked genes / smoother PDFs show the strongest effects, not A-Z.
  cond <- cond[order(cond$pvalue, -cond$waldStat),
               c("gene", setdiff(colnames(cond), "gene"))]
  write.csv(cond, file.path(ct_dir, "tradeseq_condition_test.csv"), row.names = FALSE)
  cat(sprintf("  conditionTest: %d genes, %d at FDR < %.2f (EXPLORATORY — anticonservative at n=3v4)\n",
              nrow(cond), sum(cond$fdr < SIG_FDR, na.rm = TRUE), SIG_FDR))

  # ---------------------------------------------------------------------------
  # Per-gene smoother PDFs for the top conditionTest genes
  # ---------------------------------------------------------------------------
  smooth_dir <- file.path(ct_dir, "smoothers")
  # Clear stale PDFs first: the top-N set changes between runs (e.g. after a
  # re-rank), so a plain overwrite would leave orphaned genes behind.
  unlink(smooth_dir, recursive = TRUE)
  dir.create(smooth_dir, recursive = TRUE, showWarnings = FALSE)
  cond_ok    <- cond[!is.na(cond$pvalue), ]
  top_genes  <- head(cond_ok$gene, TOP_N_PLOT)
  cat(sprintf("  Emitting %d per-gene smoother PDFs (top conditionTest genes)\n",
              length(top_genes)))
  for (gene in top_genes) {
    p_sm <- tryCatch(
      tradeSeq::plotSmoothers(sce, counts = counts_mat, gene = gene,
                              alpha = 0.4, border = TRUE) +
        labs(title = sprintf("%s — %s [tradeSeq NB-GAM]", gene, ct),
             subtitle = sprintf("condition p = %.3g (FDR %.3g) — EXPLORATORY, cell-level",
                                cond_ok$pvalue[cond_ok$gene == gene],
                                cond_ok$fdr[cond_ok$gene == gene])),
      error = function(e) NULL)
    if (!is.null(p_sm))
      ggsave(file.path(smooth_dir, sprintf("gam_%s.pdf", gene)), p_sm,
             width = 5, height = 4)
  }

  # ---------------------------------------------------------------------------
  # WT | DKO fitted-smoother heatmap (matches 15_04 aesthetic)
  # ---------------------------------------------------------------------------
  sig <- cond_ok[which(cond_ok$fdr < SIG_FDR), ]
  if (nrow(sig) < 3) {
    # fall back to the top genes so the figure still renders
    sig <- head(cond_ok, min(TOP_N_PLOT, nrow(cond_ok)))
    cat(sprintf("  < 3 genes at FDR < %.2f; heatmap shows top %d conditionTest genes\n",
                SIG_FDR, nrow(sig)))
  } else if (nrow(sig) > TOP_N_PLOT) {
    sig <- head(sig, TOP_N_PLOT)   # cap rows for readability
  }
  genes_keep <- sig$gene

  sm_mats <- tryCatch(extract_smooth_matrices(sce, genes_keep, npoints = 100),
                      error = function(e) { cat("  predictSmooth failed:",
                                                conditionMessage(e), "\n"); NULL })
  if (!is.null(sm_mats) && nrow(sm_mats$wt) >= 3) {
    mat_wt  <- sm_mats$wt
    mat_dko <- sm_mats$dko

    # Shared per-gene z-score across BOTH conditions (preserves level + shape)
    comb <- cbind(mat_wt, mat_dko)
    rmn  <- rowMeans(comb, na.rm = TRUE)
    rsd  <- apply(comb, 1, sd, na.rm = TRUE); rsd[rsd == 0 | is.na(rsd)] <- 1
    z_wt  <- (mat_wt  - rmn) / rsd
    z_dko <- (mat_dko - rmn) / rsd

    if (HEATMAP_ROW_ORDER == "shape") {
      ord   <- order(sig$pvalue)
      z_wt  <- z_wt[ord, , drop = FALSE]; z_dko <- z_dko[ord, , drop = FALSE]
      sig   <- sig[ord, ]
      cluster_rows_arg <- FALSE
    } else {
      cluster_rows_arg <- as.dendrogram(hclust(dist(cbind(z_wt, z_dko))))
    }

    col_fun  <- colorRamp2(c(-2, -1, 0, 1, 2), viridis::inferno(5))
    nlog_cap <- function(p) pmin(-log10(pmax(p, 1e-300)), 10)
    ng <- ncol(z_wt)
    row_ha <- rowAnnotation(
      `-log10 FDRcond` = anno_barplot(nlog_cap(sig$fdr),
                                      gp = gpar(fill = "#D55E00"), border = FALSE),
      annotation_name_gp = gpar(fontsize = 7),
      width = unit(2.5, "cm"))
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

    pdf(file.path(ct_dir, "tradeseq_condition_heatmap.pdf"),
        width = 12, height = max(3.5, nrow(z_wt) * 0.13 + 1.5))
    draw(ht_wt + ht_dko + row_ha,
         column_title = sprintf("%s — tradeSeq NB-GAM smoothers (EXPLORATORY, cell-level; %d genes)",
                                ct, nrow(z_wt)),
         column_title_gp = gpar(fontsize = 12, fontface = "bold"))
    dev.off()
    cat(sprintf("  condition heatmap: %d genes\n", nrow(z_wt)))
  } else {
    cat("  condition heatmap skipped (too few genes / predictSmooth failed)\n")
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
         title = sprintf("%s — Position along tumor gradient", ct)) +
    theme_publication()
  ggsave(file.path(ct_dir, "density_distribution_along_gradient.pdf"),
         p_dist, width = 5, height = 4)

  # (SCE already cached to tradeseq_sce.rds right after the fit.)
  cat(sprintf("  Done: %s\n", ct))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")
