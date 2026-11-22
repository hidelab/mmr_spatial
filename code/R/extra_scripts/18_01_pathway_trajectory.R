# 18_01_pathway_trajectory.R
# Pathway activity along the tumor-proximity trajectory.
#
# Same approach as 15_04_tumor_trajectory_working.R (tumor score, QC, decile
# bins, common support, sample x bin pseudobulk, ordered-factor difference-smooth
# condition test), but the response is PATHWAY ACTIVITY instead of gene expression:
#   * Each sample x bin pseudobulk (summed raw counts over the FULL gene panel, no
#     expression filter) is scored with PanomiR::pathwaySummary, using the same
#     recipe as 12_pathway_activity.R / PanomiR::differentialPathwayAnalysis:
#     log counts -> within-sample gene ranks, squared ("x2"), trimmed mean per
#     pathway (trim 0.025), then z-normalized across pathways within each sample.
#   * Pathway scores are treated as Gaussian -> a single Gaussian GAM with a
#     per-sample factor-smooth random effect (the 15_04 "pseudobulk_fs" spec).
#     No NB / count models.
#
# Outputs -> output/R/18_pathway_trajectory.

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
library(PanomiR)
library(babelgene)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(tidyr)
library(mgcv)
library(viridis)
library(ComplexHeatmap)
library(circlize)

# Source centralized parameters, palette, and scoring / orthology utilities
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "18_pathway_trajectory")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters (all switchable; defaults match 15_04)
# =============================================================================
# --- Tumor score -------------------------------------------------------------
SCORE_METHOD  <- "diffusion"  # "kde" | "knn" | "diffusion"  (see 00_utils.R)
KDE_BW_SCALE  <- 1.0
KDE_GRID_N    <- 200
KNN_K         <- 30
DIFFUSION_VARIANT <- "heat"
DIFFUSION_HEAT_STEPS <- 20
DIFFUSION_HEAT_ALPHA <- 2
DIFFUSION_LOCAL_FRACTION <- TRUE

TARGET_CELLTYPE <- "Tumor"

# --- Binning / common support ------------------------------------------------
N_BINS            <- 5
MIN_CELLS_PER_BIN <- 5
COMMON_SUPPORT               <- TRUE
MIN_SAMPLES_PER_COND_PER_BIN <- 2
PB_SCORE_AXIS                <- "score"   # "score" | "decile"

# --- Pathway scoring (PanomiR) -----------------------------------------------
# Same settings as PanomiR::differentialPathwayAnalysis defaults used in 12.
# MIN_PATH_SIZE (00_parameters.R) is applied to pathway genes present in the
# full panel; it is a pathway-size filter, not an expression filter.
PANOMIR_METHOD <- "x2"
PANOMIR_TRIM   <- 0.025
# Number of pathways (highest variance of bin-level activity) shown in the
# reference binned heatmap. Display only; never defines the tested set.
TOP_N_PATHWAYS <- 80

# --- GAM ---------------------------------------------------------------------
# Gaussian specs, both weighted by n_cells with the per-sample factor-smooth
# random effect s(score, sample, bs = "fs", k = FS_K, m = FS_M):
#   "pseudobulk_fs"       — 15_04's spec: shared s(score) + DKO difference smooth
#                           s(score, by = cond_ord) + cond_ord offset. p_shape and
#                           p_level; p_level is the average DKO-WT offset ADJUSTED
#                           for any shape difference.
#   "pseudobulk_fs_level" — level-only: drops the difference smooth, so WT and DKO
#                           share one trajectory shape and differ by a constant
#                           offset. p_level only (p_shape = NA). More power for the
#                           offset when shapes really are parallel; misleading when
#                           they are not (compare against pseudobulk_fs p_shape).
GAM_MODELS <- c("pseudobulk_fs", "pseudobulk_fs_level")
GAM_K     <- 5
FS_K      <- 3
FS_M      <- 2

# --- Condition-test visualization -------------------------------------------
PLOT_INDIVIDUAL_GAMS <- TRUE
SIG_FDR              <- 0.1
HEATMAP_ROW_ORDER    <- "cluster"  # "cluster" | "shape"

CELLTYPES_TRAJECTORY <- c("CD8+ T")

# --- Quality control ---------------------------------------------------------
APPLY_QC <- TRUE

stopifnot(SCORE_METHOD %in% c("kde", "knn", "diffusion"))
stopifnot(HEATMAP_ROW_ORDER %in% c("cluster", "shape"))
stopifnot(PB_SCORE_AXIS %in% c("score", "decile"))
stopifnot(FS_M >= 1, FS_K > FS_M)
stopifnot(all(GAM_MODELS %in% c("pseudobulk_fs", "pseudobulk_fs_level")))

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
        script     = "code/R/18_01_pathway_trajectory.R",
        run_time   = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        host       = as.character(Sys.info()[["nodename"]]),
        R_version  = R.version.string,
        mgcv       = pkgver("mgcv"),
        PanomiR    = pkgver("PanomiR"),
        Seurat     = pkgver("Seurat"),
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
        N_BINS = N_BINS, MIN_CELLS_PER_BIN = MIN_CELLS_PER_BIN,
        COMMON_SUPPORT = COMMON_SUPPORT,
        MIN_SAMPLES_PER_COND_PER_BIN = MIN_SAMPLES_PER_COND_PER_BIN,
        PB_SCORE_AXIS = PB_SCORE_AXIS),
    grp("pathways",
        MSIGDB_PATHWAYS = MSIGDB_PATHWAYS, MIN_PATH_SIZE = MIN_PATH_SIZE,
        PANOMIR_METHOD = PANOMIR_METHOD, PANOMIR_TRIM = PANOMIR_TRIM,
        gene_filter = "none (full panel)", TOP_N_PATHWAYS = TOP_N_PATHWAYS),
    grp("gam",
        GAM_MODELS = GAM_MODELS, family = "gaussian",
        GAM_K = GAM_K, FS_K = FS_K, FS_M = FS_M),
    grp("condition_test_viz",
        PLOT_INDIVIDUAL_GAMS = PLOT_INDIVIDUAL_GAMS,
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
               "18_01 pathway-trajectory run parameters",
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
# Pathway reference over the FULL panel
# =============================================================================
# Mouse panel symbols -> human ENSEMBL by orthology (as in 12); unmapped genes
# cannot enter any pathway. No expression-level filter is applied.
panel_genes       <- rownames(seu[["RNA"]])
symbol_to_ensembl <- convert_mouse_to_human_ensembl(panel_genes)
cat(sprintf("Panel: %d genes, %d mapped to human ENSEMBL\n",
            length(panel_genes), length(symbol_to_ensembl)))

# Keep pathways with >= MIN_PATH_SIZE genes in the panel (PanomiR's own
# pathway cleaning step in differentialPathwayAnalysis).
pathways2    <- as.data.frame(readRDS(MSIGDB_PATHWAYS))
displayNames <- read.csv(MSIGDB_NAMES)[, c("standard_name", "display_name")]
path_sizes   <- pathways2 %>%
  filter(ENSEMBL %in% symbol_to_ensembl) %>%
  count(Pathway, name = "n_panel_genes")
pathway_ref  <- pathways2 %>%
  filter(Pathway %in% path_sizes$Pathway[path_sizes$n_panel_genes >= MIN_PATH_SIZE])
cat(sprintf("Pathways with >= %d panel genes: %d\n",
            MIN_PATH_SIZE, length(unique(pathway_ref$Pathway))))

pathway_display <- function(p) {
  dn <- displayNames$display_name[match(p, displayNames$standard_name)]
  truncate_label(pathway_label(p, dn), PATHWAY_NAME_TRUNC)
}

# Pathway activity of each column of a (mouse symbol) summed-count matrix.
# Mirrors PanomiR::differentialPathwayAnalysis: log counts -> pathwaySummary
# (within-column gene ranks over the full panel, "x2", trimmed mean) -> drop
# NA rows -> z-normalize each column across pathways. Every step is per column,
# so each pseudobulk's score does not depend on which other groups are scored.
score_pathways <- function(count_mat) {
  m <- as.matrix(count_mat[rownames(count_mat) %in% names(symbol_to_ensembl), , drop = FALSE])
  rownames(m) <- symbol_to_ensembl[rownames(m)]
  ps <- PanomiR::pathwaySummary(log(m), pathway_ref, id = "ENSEMBL",
                                zNormalize = FALSE, method = PANOMIR_METHOD,
                                trim = PANOMIR_TRIM)
  ps <- ps[rowSums(is.na(ps)) == 0, , drop = FALSE]
  ps <- apply(ps, 2, function(x) (x - mean(x)) / stats::sd(x))
  as.matrix(ps)
}

# =============================================================================
# Compute tumor score (per sample) via 00_utils.R
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

SCORE_LAB   <- sprintf("Tumor Proximity Score (%s)", toupper(SCORE_METHOD))
TRAJ_X_LAB  <- sprintf("Tumor proximity (decile of %s score)", toupper(SCORE_METHOD))
GAM_X_LAB   <- if (PB_SCORE_AXIS == "score") {
                 sprintf("Tumor proximity (%s local tumor fraction)", toupper(SCORE_METHOD))
               } else TRAJ_X_LAB
ACT_LAB     <- "Pathway activity (PanomiR z-score)"

# =============================================================================
# Helpers
# =============================================================================

# Decile bins of a score vector; returns bin assignment + midpoints (as 15_04).
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

# --- Pseudobulk (sample x bin) ----------------------------------------------
# Same grouping, min-cells rule, real-score x coordinate, common-support filter
# and canonical row order as 15_04's make_pseudobulk(), but always SUMS raw
# counts (PanomiR input, as in 12's pseudobulk counts).
make_pseudobulk <- function(count_mat, score_vec, sample_vec, condition_vec,
                            n_bins = N_BINS, min_cells = MIN_CELLS_PER_BIN,
                            common_support = COMMON_SUPPORT,
                            min_samp_per_cond = MIN_SAMPLES_PER_COND_PER_BIN) {
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
    pb_list[[g]] <- Matrix::rowSums(count_mat[, cell_idx, drop = FALSE])
    pb_meta[[g]] <- data.frame(
      group = g, sample = parts[1],
      condition = as.character(condition_vec[cell_idx[1]]),
      bin = bin_id,
      score = mean(score_vec[cell_idx]),
      n_cells = length(cell_idx), stringsAsFactors = FALSE)
  }
  if (length(pb_list) == 0) return(NULL)
  meta <- bind_rows(pb_meta)
  meta$condition <- factor(meta$condition,
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION))

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
    if (nrow(meta) < 4) return(NULL)
  }

  ord <- order(meta$sample, meta$bin)
  meta <- meta[ord, ]
  pb_list <- pb_list[meta$group]

  list(counts = do.call(cbind, pb_list), meta = meta, support = support)
}

# --- GAM fitting -------------------------------------------------------------
# Gaussian pseudobulk_fs spec from 15_04 with the ordered-factor difference
# smooth: p_shape = DKO-vs-WT trajectory-shape difference (by-smooth),
# p_level = constant DKO offset (cond_ord.L). "pseudobulk_fs_level" drops the
# by-smooth (p_shape = NA). Returns population curves (random effect excluded)
# and per-sample curves (random effect included).
fit_trajectory_gam <- function(df, model, gam_k = GAM_K, fs_k = FS_K, fs_m = FS_M) {
  df$condition <- droplevels(factor(df$condition,
                                    levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
  df$sample <- droplevels(factor(df$sample))
  if (nlevels(df$condition) < 2 || nlevels(df$sample) < 2) return(list(ok = FALSE))

  df$cond_ord <- ordered(as.character(df$condition),
                         levels = c(REFERENCE_CONDITION, TEST_CONDITION))

  form <- if (model == "pseudobulk_fs_level") {
    expr ~ s(score, k = gam_k) + cond_ord +
           s(score, sample, bs = "fs", k = fs_k, m = fs_m)
  } else {
    expr ~ s(score, k = gam_k) + s(score, by = cond_ord, k = gam_k) +
           cond_ord + s(score, sample, bs = "fs", k = fs_k, m = fs_m)
  }
  excl <- "s(score,sample)"

  fit <- tryCatch(gam(form, data = df, weights = df$weight, family = gaussian(),
                      method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit)) return(list(ok = FALSE))

  sm <- summary(fit)
  s_rows   <- rownames(sm$s.table)
  diff_row <- s_rows[grepl(":cond_ord", s_rows, fixed = TRUE)]
  p_shape  <- if (length(diff_row)) sm$s.table[diff_row[1], "p-value"] else NA_real_
  p_rows   <- rownames(sm$p.table)
  lev_row  <- p_rows[grepl("cond_ord", p_rows, fixed = TRUE)]
  pcol     <- grep("^Pr\\(", colnames(sm$p.table))
  p_level  <- if (length(lev_row) && length(pcol))
                sm$p.table[lev_row[1], pcol[1]] else NA_real_
  # A NEGATIVE p-value is mgcv's numerical failure on a near-singular fit (e.g. a
  # saturated fs term), not underflow; set it to NA instead of flooring it to a
  # spurious top hit. Exact 0 (underflow) is floored at machine precision.
  if (!is.na(p_shape) && p_shape < 0) p_shape <- NA_real_
  if (!is.na(p_level) && p_level < 0) p_level <- NA_real_
  if (!is.na(p_shape)) p_shape <- max(p_shape, .Machine$double.eps)
  if (!is.na(p_level)) p_level <- max(p_level, .Machine$double.eps)

  newd <- expand.grid(
    score = seq(min(df$score), max(df$score), length.out = 100),
    cond_ord = ordered(c(REFERENCE_CONDITION, TEST_CONDITION),
                       levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
  newd$sample <- df$sample[1]   # dummy; excluded from prediction
  pred <- predict(fit, newdata = newd, se.fit = TRUE,
                  exclude = excl, newdata.guaranteed = TRUE)
  newd$condition <- factor(as.character(newd$cond_ord),
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  newd$fit   <- pred$fit
  newd$se    <- pred$se.fit
  newd$ci_lo <- pred$fit - 1.96 * pred$se.fit
  newd$ci_hi <- pred$fit + 1.96 * pred$se.fit

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
  counts_ct    <- GetAssayData(seu, assay = "RNA", layer = "counts")[, ct_idx]

  cat(sprintf("  Total cells: %d (WT: %d, DKO: %d)\n", length(ct_idx),
              sum(ct_condition == REFERENCE_CONDITION),
              sum(ct_condition == TEST_CONDITION)))

  # ---------------------------------------------------------------------------
  # Binning diagnostics (cells per sample x bin)
  # ---------------------------------------------------------------------------
  diag_dir <- file.path(ct_dir, "diagnostics")
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

  bininfo  <- make_bins(ct_score)
  df_score <- data.frame(score = ct_score, sample = ct_sample, condition = ct_condition,
                         bin = as.integer(bininfo$bins))
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
  write.csv(counts, file.path(diag_dir, "cells_per_sample_bin.csv"), row.names = FALSE)

  # ---------------------------------------------------------------------------
  # Binned pathway-activity heatmap (all cells per bin) — reference view
  # ---------------------------------------------------------------------------
  bin_counts <- vapply(seq_len(N_BINS), function(b)
                  Matrix::rowSums(counts_ct[, which(as.integer(bininfo$bins) == b), drop = FALSE]),
                  numeric(nrow(counts_ct)))
  rownames(bin_counts) <- rownames(counts_ct)
  colnames(bin_counts) <- paste0("bin_", seq_len(N_BINS))
  bin_act <- score_pathways(bin_counts[, colSums(bin_counts) > 0, drop = FALSE])
  write.csv(bin_act, file.path(ct_dir, "pathway_activity_binned.csv"))

  row_vars  <- apply(bin_act, 1, var, na.rm = TRUE)
  plot_paths <- head(names(sort(row_vars, decreasing = TRUE)), TOP_N_PATHWAYS)
  hm_scaled <- t(scale(t(bin_act[plot_paths, , drop = FALSE])))
  hm_df <- as.data.frame(hm_scaled) %>%
    mutate(pathway = pathway_display(rownames(hm_scaled))) %>%
    pivot_longer(-pathway, names_to = "bin", values_to = "zscore") %>%
    mutate(bin_num = as.integer(gsub("bin_", "", bin)),
           pathway = factor(pathway, levels = rev(unique(pathway_display(plot_paths)))))
  p_hm <- ggplot(hm_df, aes(x = bin_num, y = pathway, fill = zscore)) +
    geom_tile() +
    scale_fill_viridis_c(option = "inferno", name = "Z-score") +
    scale_x_continuous(breaks = c(1, N_BINS), labels = c("Low", "High")) +
    labs(x = TRAJ_X_LAB, y = NULL,
         title = sprintf("%s pathway activity along tumor score", ct)) +
    theme_publication() + theme(axis.text.y = element_text(size = 6))
  ggsave(file.path(ct_dir, "heatmap_pathway_trajectory_combined.pdf"),
         p_hm, width = 9, height = max(4, length(plot_paths) * 0.2),
         limitsize = FALSE)

  # ---------------------------------------------------------------------------
  # Pseudobulk pathway activity (shared by all GAM specs)
  # ---------------------------------------------------------------------------
  pb <- make_pseudobulk(counts_ct, ct_score, ct_sample, ct_condition)
  if (is.null(pb) || nrow(pb$meta) < 10) {
    cat("  too few pseudobulk groups (common support too small?), skipping\n"); next
  }
  if (!is.null(pb$support)) {
    kb <- pb$support$keep_bins
    cat(sprintf("  common support: %d/%d bins retained (bins %s); score range [%.3f, %.3f]\n",
                length(kb), N_BINS, paste(range(kb), collapse = "-"),
                min(pb$meta$score), max(pb$meta$score)))
  }

  pb_act <- score_pathways(pb$counts)
  pb_act <- pb_act[, pb$meta$group, drop = FALSE]
  write.csv(pb_act, file.path(ct_dir, "pathway_activity_pseudobulk.csv"))
  write.csv(pb$meta, file.path(ct_dir, "pseudobulk_meta.csv"), row.names = FALSE)
  test_paths <- rownames(pb_act)

  # ---------------------------------------------------------------------------
  # GAM condition test, loop over specs
  # ---------------------------------------------------------------------------
  gs_all <- list()
  for (model in GAM_MODELS) {
    gam_dir <- file.path(ct_dir, paste0("gam_", model))
    dir.create(gam_dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  [%s] fitting %d pathways on %d pseudobulk groups\n",
                model, length(test_paths), nrow(pb$meta)))

    pb_x <- if (PB_SCORE_AXIS == "score") pb$meta$score
            else (pb$meta$bin - 0.5) / N_BINS

    gam_summary <- list()
    fit_curves  <- list()
    plot_store  <- list()
    for (path in test_paths) {
      df_gam <- data.frame(expr = as.numeric(pb_act[path, ]),
                           score = pb_x,
                           condition = pb$meta$condition,
                           sample = pb$meta$sample, weight = pb$meta$n_cells)

      res <- fit_trajectory_gam(df_gam, model)
      if (!isTRUE(res$ok)) next

      gam_summary[[path]] <- data.frame(pathway = path, p_shape = res$p_shape,
                                        p_level = res$p_level)
      fit_curves[[path]]  <- res$pred[, c("score", "condition", "fit")]
      if (PLOT_INDIVIDUAL_GAMS)
        plot_store[[path]] <- list(pred = res$pred, pred_sample = res$pred_sample,
                                   p_shape = res$p_shape, p_level = res$p_level)
    }

    if (length(gam_summary) == 0) { cat("  no GAM fits succeeded\n"); next }

    gs <- bind_rows(gam_summary)
    gs$fdr_shape <- p.adjust(gs$p_shape, method = "BH")
    gs$fdr_level <- p.adjust(gs$p_level, method = "BH")
    gs <- dplyr::left_join(gs, displayNames, by = c("pathway" = "standard_name"))
    gs <- gs[order(gs$p_shape, gs$p_level), ]   # level-only spec: p_shape all NA
    write.csv(gs, file.path(ct_dir, sprintf("gam_condition_test_%s.csv", model)),
              row.names = FALSE)
    cat(sprintf("  [%s] GAM results: %d pathways\n", model, nrow(gs)))
    gs_all[[model]] <- gs

    # ---------------------------------------------------------------------------
    # Per-pathway spline PDFs for SIGNIFICANT pathways only
    # ---------------------------------------------------------------------------
    if (PLOT_INDIVIDUAL_GAMS) {
      sig_paths <- gs$pathway[which(gs$fdr_shape < SIG_FDR | gs$fdr_level < SIG_FDR)]
      sig_paths <- intersect(sig_paths, names(plot_store))
      cat(sprintf("  [%s] per-pathway spline PDFs for %d significant pathways (FDR < %.2f)\n",
                  model, length(sig_paths), SIG_FDR))
      ps_dir <- file.path(gam_dir, "per_sample_curves")
      dir.create(ps_dir, recursive = TRUE, showWarnings = FALSE)
      pop_keys <- setNames(
        c(CONDITION_COLORS[[REFERENCE_CONDITION]], CONDITION_COLORS[[TEST_CONDITION]]),
        c(paste0(REFERENCE_CONDITION, " (fitted mean)"),
          paste0(TEST_CONDITION, " (fitted mean)")))
      ps_pal <- c(SAMPLE_COLORS, pop_keys)
      for (path in sig_paths) {
        st <- plot_store[[path]]
        subt <- paste(c(
          if (!is.na(st$p_shape)) sprintf("shape p = %.3g", st$p_shape),
          if (!is.na(st$p_level)) sprintf("level p = %.3g", st$p_level)
        ), collapse = "   ")
        ttl <- pathway_display(path)
        p_gam <- ggplot(st$pred, aes(x = score, y = fit, color = condition, fill = condition)) +
          geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
          geom_line(linewidth = 1) +
          scale_color_manual(values = CONDITION_COLORS) +
          scale_fill_manual(values = CONDITION_COLORS) +
          labs(x = GAM_X_LAB, y = ACT_LAB,
               title = sprintf("%s\n%s [%s]", ttl, ct, model),
               subtitle = subt) +
          theme_publication() + theme(plot.title = element_text(size = 9))
        ggsave(file.path(gam_dir, sprintf("gam_%s.pdf", path)), p_gam, width = 6, height = 4.5)

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
            labs(x = GAM_X_LAB, y = ACT_LAB,
                 title = sprintf("%s\n%s [%s] — per-sample fits", ttl, ct, model),
                 subtitle = subt) +
            theme_publication() + theme(plot.title = element_text(size = 9))
          ggsave(file.path(ps_dir, sprintf("gam_%s_per_sample.pdf", path)),
                 p_ps, width = 6.5, height = 4.5)
        }
      }
    }

    # ---------------------------------------------------------------------------
    # Significant-trajectory heatmap: WT | DKO GAM-fitted curves (shared z-score)
    # ---------------------------------------------------------------------------
    sig <- gs[which(gs$fdr_shape < SIG_FDR | gs$fdr_level < SIG_FDR), ]
    if (nrow(sig) >= 3) {
      paths_keep <- sig$pathway
      ng <- sum(fit_curves[[paths_keep[1]]]$condition == REFERENCE_CONDITION)
      get_curve <- function(p, cond) {
        d <- fit_curves[[p]]; d <- d[d$condition == cond, ]
        d$fit[order(d$score)]
      }
      mat_wt  <- t(vapply(paths_keep, get_curve, numeric(ng), cond = REFERENCE_CONDITION))
      mat_dko <- t(vapply(paths_keep, get_curve, numeric(ng), cond = TEST_CONDITION))

      comb <- cbind(mat_wt, mat_dko)
      rmn  <- rowMeans(comb, na.rm = TRUE)
      rsd  <- apply(comb, 1, sd, na.rm = TRUE); rsd[rsd == 0 | is.na(rsd)] <- 1
      z_wt  <- (mat_wt  - rmn) / rsd
      z_dko <- (mat_dko - rmn) / rsd
      rownames(z_wt) <- rownames(z_dko) <- pathway_display(paths_keep)

      if (HEATMAP_ROW_ORDER == "shape") {
        ord   <- order(sig$p_shape)
        z_wt  <- z_wt[ord, , drop = FALSE]; z_dko <- z_dko[ord, , drop = FALSE]
        sig   <- sig[ord, ]
        cluster_rows_arg <- FALSE
      } else {
        cluster_rows_arg <- as.dendrogram(hclust(dist(cbind(z_wt, z_dko))))
      }

      sig_sh   <- !is.na(sig$fdr_shape) & sig$fdr_shape < SIG_FDR
      sig_lv   <- !is.na(sig$fdr_level) & sig$fdr_level < SIG_FDR
      category <- ifelse(sig_sh & sig_lv, "both", ifelse(sig_sh, "shape", "level"))
      cat_cols <- c(shape = "#D55E00", level = "#0072B2", both = "#009E73")
      col_fun  <- colorRamp2(c(-2, -1, 0, 1, 2), viridis::inferno(5))

      nlog_cap <- function(p) pmin(-log10(pmax(p, 1e-300)), 10)
      # Level-only spec has no shape test -> omit the FDRshape bar.
      fdr_bars <- list(`-log10 FDRlevel` = anno_barplot(nlog_cap(sig$fdr_level),
                                         gp = gpar(fill = "#0072B2"), border = FALSE, which = "row"))
      if (any(!is.na(sig$fdr_shape)))
        fdr_bars <- c(list(`-log10 FDRshape` = anno_barplot(nlog_cap(sig$fdr_shape),
                                         gp = gpar(fill = "#D55E00"), border = FALSE, which = "row")),
                      fdr_bars)
      row_ha <- do.call(rowAnnotation, c(
        list(class = category), fdr_bars,
        list(col = list(class = cat_cols),
             annotation_name_gp = gpar(fontsize = 7),
             annotation_legend_param = list(class = list(title = "Sig. in")),
             width = unit(3, "cm"))))

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
                       row_names_gp = gpar(fontsize = 6),
                       row_names_max_width = unit(12, "cm"),
                       bottom_annotation = col_note,
                       width = unit(6.5, "cm"), border = TRUE)
      ht_dko <- Heatmap(z_dko, cluster_rows = FALSE, cluster_columns = FALSE,
                        col = col_fun, name = " ", column_title = TEST_CONDITION,
                        show_column_names = FALSE, show_row_names = FALSE,
                        bottom_annotation = col_note,
                        width = unit(6.5, "cm"), border = TRUE,
                        show_heatmap_legend = FALSE)

      pdf(file.path(ct_dir, sprintf("trajectory_heatmap_significant_%s.pdf", model)),
          width = 20, height = max(3.5, nrow(z_wt) * 0.13 + 1.5))
      draw(ht_wt + ht_dko + row_ha,
           column_title = sprintf("%s significant pathway trajectories (BH-FDR < %.2f, %d pathways)",
                                  ct, SIG_FDR, nrow(z_wt)),
           column_title_gp = gpar(fontsize = 12, fontface = "bold"))
      dev.off()
      cat(sprintf("  [%s] significant-trajectory heatmap: %d pathways (%d shape, %d level, %d both)\n",
                  model, nrow(z_wt), sum(category == "shape"),
                  sum(category == "level"), sum(category == "both")))
    } else {
      cat(sprintf("  [%s] < 3 pathways pass BH-FDR < %.2f; heatmap skipped\n", model, SIG_FDR))
    }
  }

  # ---------------------------------------------------------------------------
  # Level test: full (shape-adjusted) vs level-only spec
  # ---------------------------------------------------------------------------
  if (all(c("pseudobulk_fs", "pseudobulk_fs_level") %in% names(gs_all))) {
    cmp <- dplyr::full_join(
      gs_all[["pseudobulk_fs"]][, c("pathway", "display_name", "p_shape", "fdr_shape",
                                    "p_level", "fdr_level")],
      gs_all[["pseudobulk_fs_level"]][, c("pathway", "p_level", "fdr_level")],
      by = "pathway", suffix = c("_full", "_levelonly"))
    cmp$sig_level <- with(cmp, ifelse(
      fdr_level_full < SIG_FDR & fdr_level_levelonly < SIG_FDR, "both",
      ifelse(fdr_level_full < SIG_FDR, "full only",
      ifelse(fdr_level_levelonly < SIG_FDR, "level-only only", "neither"))))
    cmp <- cmp[order(cmp$p_level_levelonly), ]
    write.csv(cmp, file.path(ct_dir, "level_test_full_vs_levelonly.csv"), row.names = FALSE)
    print(table(cmp$sig_level, useNA = "ifany"))

    # Pathways with a significant shape difference are where the level-only
    # offset is least trustworthy; flag them in the scatter.
    cmp$shape_sig <- !is.na(cmp$fdr_shape) & cmp$fdr_shape < SIG_FDR
    p_cmp <- ggplot(cmp, aes(x = -log10(p_level_full), y = -log10(p_level_levelonly),
                             color = sig_level, shape = shape_sig)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(alpha = 0.7, size = 1.6) +
      scale_color_manual(values = c(both = "#009E73", `full only` = "#D55E00",
                                    `level-only only` = "#0072B2", neither = "grey70"),
                         name = sprintf("BH-FDR < %.2f", SIG_FDR)) +
      scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                         name = sprintf("Shape FDR < %.2f", SIG_FDR)) +
      labs(x = "-log10 p_level, full model (shape-adjusted)",
           y = "-log10 p_level, level-only model",
           title = sprintf("%s pathway level test: full vs level-only", ct)) +
      theme_publication()
    ggsave(file.path(ct_dir, "level_test_full_vs_levelonly.pdf"), p_cmp, width = 7, height = 5.5)
  }

  cat(sprintf("  Done: %s\n", ct))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")
