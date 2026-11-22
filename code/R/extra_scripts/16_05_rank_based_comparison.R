# 16_05_rank_based_comparison.R
# Rank-based (Spearman) concordance of the 15_04 tumor-proximity trajectory results
# ACROSS parameter settings and GAM models.
#
# Motivation: we ran 15_04 under 14 different parameter settings, saved to
#   output/R/15_tumor_trajectory_working_v1 ... _v14 (some suffixed, e.g. _v10_K5).
# Each setting produced one gam_condition_test_<model>.csv per fitted GAM model
# (columns: gene, p_shape, p_level, fdr_shape, fdr_level). Before locking a model +
# parameters for the paper, we want to know how STABLE the per-gene ranking is:
#   * Do different parameter settings rank genes the same way (robust) or not?
#   * Do the GAM models agree with each other?
# Spearman correlation between two setups' p-value vectors = correlation of the gene
# RANKS (Spearman is invariant to any monotone transform of p, so correlating raw
# p_shape == correlating -log10 p == correlating the significance rank). High rho =
# the two setups order genes concordantly; low rho = the choice matters.
#
# We do this SEPARATELY for p_shape (trajectory-shape test — the spatial question)
# and p_level (level offset ~ bulk DE), producing two setup x setup heatmaps.
#
# Handles ragged coverage: a setting may have {pseudobulk, pseudobulk_fs} only while
# another has all five specs, and v1 has no CSVs at all — every (experiment, model)
# combination that has a CSV becomes one column; missing ones are simply absent.
# Correlations use pairwise-complete genes.
#
# Outputs -> output/R/16_rank_based_comparison/

# --- Thread limits (set BEFORE BLAS/library load) ----------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS,
           NUMEXPR_NUM_THREADS = N_THREADS)

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Sourced only for PROJECT_DIR (+ palette conventions).
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Parameters
# =============================================================================
OUT_ROOT   <- file.path(PROJECT_DIR, "output", "R")
# Version result folders to scan. Glob matches _v1 .. _v14 incl. suffixed variants.
VERSION_GLOB  <- "15_tumor_trajectory_working_v*"
CELLTYPE_GLOB <- "CD8pos*"                 # celltype subfolder (CD8pos_T)
METRICS       <- c("p_shape", "p_level")   # one heatmap per metric
COR_METHOD    <- "spearman"                # rank-based
SHOW_CELL_VALUES <- TRUE                    # print rho in each cell
# BH-FDR cutoffs at which to count significant genes per setup per test.
SIG_CUTOFFS   <- c(0.05, 0.1, 0.25)

OUTPUT_DIR <- file.path(OUT_ROOT, "16_rank_based_comparison")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Model color strip (matches the specs defined in 15_04).
MODEL_COLORS <- c(
  pseudobulk       = "#8C8C8C",
  pseudobulk_fs    = "#D55E00",
  pseudobulk_rs    = "#0072B2",
  pseudobulk_nb    = "#009E73",
  pseudobulk_nb_fs = "#CC79A7"
)

# GAM_K (spline basis dimension, s(score, k=GAM_K)) color strip. Read per experiment
# from run_parameters.csv (see below). "unknown" = pre-parameter-report runs (v2-v6).
K_COLORS <- c(`K=3` = "#FDE725", `K=4` = "#35B779", `K=5` = "#31688E",
              unknown = "grey80")

# =============================================================================
# Discover result CSVs and load p-values
# =============================================================================
version_dirs <- sort(Sys.glob(file.path(OUT_ROOT, VERSION_GLOB)))
if (length(version_dirs) == 0)
  stop("No version folders matching ", VERSION_GLOB, " under ", OUT_ROOT)

# exp_id = folder name minus the shared prefix, e.g. "v8", "v10_K5", "v12_20bins".
strip_prefix <- function(p) sub("^15_tumor_trajectory_working_", "", basename(p))

records <- list()   # one row per (exp, model) CSV found
long    <- list()   # long-format gene/setup/p_shape/p_level
for (d in version_dirs) {
  exp_id  <- strip_prefix(d)
  ct_dirs <- Sys.glob(file.path(d, CELLTYPE_GLOB))
  if (length(ct_dirs) == 0) {
    records[[length(records) + 1]] <-
      data.frame(exp = exp_id, model = NA, n_genes = 0, present = FALSE)
    next
  }
  csvs <- Sys.glob(file.path(ct_dirs[1], "gam_condition_test_*.csv"))
  if (length(csvs) == 0) {
    records[[length(records) + 1]] <-
      data.frame(exp = exp_id, model = NA, n_genes = 0, present = FALSE)
    next
  }
  for (f in csvs) {
    model <- sub(".*gam_condition_test_(.*)\\.csv$", "\\1", basename(f))
    df <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(df) || !all(c("gene", "p_shape", "p_level") %in% names(df))) {
      records[[length(records) + 1]] <-
        data.frame(exp = exp_id, model = model, n_genes = 0, present = FALSE)
      next
    }
    p_shape <- suppressWarnings(as.numeric(df$p_shape))
    p_level <- suppressWarnings(as.numeric(df$p_level))
    # Prefer the run's own BH-FDR columns (BH over that setup's tested-gene set, so
    # the denominator m is the setup's own); recompute BH from p only if absent
    # (older CSVs predating the fdr_* columns in 15_04).
    fdr_shape <- if ("fdr_shape" %in% names(df)) suppressWarnings(as.numeric(df$fdr_shape))
                 else p.adjust(p_shape, method = "BH")
    fdr_level <- if ("fdr_level" %in% names(df)) suppressWarnings(as.numeric(df$fdr_level))
                 else p.adjust(p_level, method = "BH")
    setup <- paste(exp_id, model, sep = " | ")
    long[[setup]] <- data.frame(
      gene = df$gene, setup = setup, exp = exp_id, model = model,
      p_shape = p_shape, p_level = p_level,
      fdr_shape = fdr_shape, fdr_level = fdr_level,
      stringsAsFactors = FALSE)
    records[[length(records) + 1]] <-
      data.frame(exp = exp_id, model = model, n_genes = nrow(df), present = TRUE)
  }
}

inventory <- bind_rows(records)
write.csv(inventory, file.path(OUTPUT_DIR, "setup_inventory.csv"), row.names = FALSE)
cat("Setup inventory:\n"); print(inventory, row.names = FALSE)

long_df <- bind_rows(long)
if (nrow(long_df) == 0) stop("No gam_condition_test_*.csv files found to compare.")
n_setups <- dplyr::n_distinct(long_df$setup)
cat(sprintf("\nLoaded %d setups (experiment x model) across %d version folders.\n",
            n_setups, length(version_dirs)))

# --- Per-experiment GAM_K (spline basis dim) from run_parameters.csv ----------
# run_parameters.csv is written only by runs from 2026-09-10 onward (v7+); earlier
# versions (v2-v6) predate the parameter report, so their GAM_K is recorded as
# "unknown". FS_K is constant (3) across all runs, so GAM_K is the informative knot.
read_param <- function(version_dir, pname) {
  f <- file.path(version_dir, "run_parameters.csv")
  if (!file.exists(f)) return(NA_character_)
  p <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(p) || !all(c("parameter", "value") %in% names(p))) return(NA_character_)
  v <- p$value[match(pname, p$parameter)]
  if (length(v) == 0 || is.na(v)) NA_character_ else trimws(as.character(v))
}
exp_gamk <- setNames(
  vapply(version_dirs, read_param, character(1), pname = "GAM_K"),
  strip_prefix(version_dirs))

# Metadata per setup, for annotation + ordering (by model, then experiment).
setup_meta <- long_df %>% distinct(setup, exp, model) %>%
  mutate(gam_k = unname(exp_gamk[exp]),
         gam_k = ifelse(is.na(gam_k) | gam_k == "NA", "unknown", paste0("K=", gam_k)))
cat("\nGAM_K per experiment:\n")
print(setup_meta %>% distinct(exp, gam_k) %>% arrange(exp), row.names = FALSE)

# =============================================================================
# Significant-gene counts per setup per test at FDR cutoffs
# =============================================================================
# For each setup (experiment x model) and each test (shape / level), count the
# genes clearing BH-FDR < each SIG_CUTOFFS threshold. This summarises HOW MANY
# hits each modelling choice produces — the counterpart to the rank-concordance
# above (which asks whether they agree on ORDER). Reported as a long CSV, a wide
# CSV, and a setups x (test x cutoff) count heatmap.
sig_rows <- list()
for (s in unique(long_df$setup)) {
  d <- long_df[long_df$setup == s, ]
  for (test in c("shape", "level")) {
    fdr <- d[[paste0("fdr_", test)]]
    n_tested <- sum(!is.na(fdr))
    for (cut in SIG_CUTOFFS) {
      n_sig <- sum(fdr < cut, na.rm = TRUE)
      sig_rows[[length(sig_rows) + 1]] <- data.frame(
        setup = s, exp = d$exp[1], model = d$model[1],
        test = test, cutoff = cut, n_tested = n_tested, n_sig = n_sig,
        frac_sig = if (n_tested > 0) n_sig / n_tested else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
}
sig_long <- bind_rows(sig_rows)
write.csv(sig_long, file.path(OUTPUT_DIR, "significant_gene_counts_long.csv"),
          row.names = FALSE)

# Wide table: one row per setup, columns = <test>_fdr<cutoff> counts (+ n_tested).
combos     <- expand.grid(test = c("shape", "level"), cutoff = SIG_CUTOFFS,
                          stringsAsFactors = FALSE)
combos     <- combos[order(match(combos$test, c("shape", "level")), combos$cutoff), ]
combo_ids  <- paste0(combos$test, "_fdr", combos$cutoff)
combo_labs <- sprintf("%s FDR<%g", combos$test, combos$cutoff)

sig_wide <- sig_long %>%
  mutate(col = paste0(test, "_fdr", cutoff)) %>%
  select(setup, exp, model, col, n_sig) %>%
  pivot_wider(names_from = col, values_from = n_sig) %>%
  left_join(sig_long %>% distinct(setup, n_tested), by = "setup")
# Canonical column order (setup meta, n_tested, then the test x cutoff counts).
sig_wide <- sig_wide[, c("setup", "exp", "model", "n_tested",
                         intersect(combo_ids, names(sig_wide)))]
sig_wide <- sig_wide %>% arrange(match(model, names(MODEL_COLORS)), exp)
write.csv(sig_wide, file.path(OUTPUT_DIR, "significant_gene_counts_wide.csv"),
          row.names = FALSE)
cat("\nSignificant-gene counts per setup (BH-FDR cutoffs ",
    paste(SIG_CUTOFFS, collapse = ", "), "):\n", sep = "")
print(as.data.frame(sig_wide), row.names = FALSE)

# --- Count heatmap: setups (rows) x test-cutoff (cols) ------------------------
make_sig_count_heatmap <- function() {
  ord_setups <- setup_meta %>%
    arrange(match(model, names(MODEL_COLORS)), exp) %>% pull(setup)
  ord_setups <- intersect(ord_setups, unique(sig_long$setup))

  mat <- matrix(NA_integer_, nrow = length(ord_setups), ncol = length(combo_ids),
                dimnames = list(ord_setups, combo_ids))
  for (r in seq_along(ord_setups)) for (cc in seq_along(combo_ids)) {
    v <- sig_long$n_sig[sig_long$setup == ord_setups[r] &
                        sig_long$test  == combos$test[cc] &
                        sig_long$cutoff == combos$cutoff[cc]]
    if (length(v)) mat[r, cc] <- v[1]
  }

  # Colour by count (shared scale); the exact count is printed in every cell.
  mx      <- max(mat, na.rm = TRUE); if (!is.finite(mx) || mx <= 0) mx <- 1
  col_fun <- colorRamp2(c(0, mx / 2, mx), c("#F7F7F7", "#FEB24C", "#BD0026"))

  ann_models <- setup_meta$model[match(ord_setups, setup_meta$setup)]
  ann_k      <- setup_meta$gam_k[match(ord_setups, setup_meta$setup)]
  left_ann <- rowAnnotation(
    model = ann_models, GAM_K = ann_k,
    col = list(model = MODEL_COLORS, GAM_K = K_COLORS),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(model = list(title = "GAM model"),
                                   GAM_K = list(title = "GAM_K (spline k)")))

  cell_fs <- if (length(ord_setups) <= 30) 8 else 6
  ht <- Heatmap(mat, name = "n genes",
                col = col_fun, na_col = "grey85",
                cluster_rows = FALSE, cluster_columns = FALSE,
                left_annotation = left_ann,
                column_split = factor(combos$test, levels = c("shape", "level")),
                column_title_gp = gpar(fontsize = 11, fontface = "bold"),
                column_labels = sprintf("FDR<%g", combos$cutoff),
                row_names_gp = gpar(fontsize = if (length(ord_setups) <= 30) 8 else 6),
                column_names_gp = gpar(fontsize = 9), column_names_rot = 45,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  v <- mat[i, j]
                  if (!is.na(v)) grid.text(v, x, y,
                    gp = gpar(fontsize = cell_fs,
                              col = if (v > mx * 0.6) "white" else "black"))
                },
                heatmap_legend_param = list(at = c(0, round(mx / 2), mx)))

  pdf(file.path(OUTPUT_DIR, "significant_gene_counts_heatmap.pdf"),
      width = max(7, length(combo_ids) * 0.9 + 5),
      height = max(6, length(ord_setups) * 0.32 + 3))
  draw(ht, column_title = "Significant genes per setup per test (BH-FDR cutoffs; CD8+ T)",
       column_title_gp = gpar(fontsize = 13, fontface = "bold"),
       heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat(sprintf("Wrote significant_gene_counts_heatmap.pdf (%d setups x %d test-cutoff cols)\n",
              length(ord_setups), length(combo_ids)))
}
make_sig_count_heatmap()

# =============================================================================
# Per-metric Spearman correlation + heatmap
# =============================================================================
# Build a gene x setup matrix for one metric.
to_matrix <- function(metric) {
  w <- long_df %>%
    select(gene, setup, value = all_of(metric)) %>%
    # a setup should have one value per gene; if a gene repeats, take the first
    group_by(gene, setup) %>% summarise(value = dplyr::first(value), .groups = "drop") %>%
    pivot_wider(names_from = setup, values_from = value)
  m <- as.matrix(w[, -1, drop = FALSE]); rownames(m) <- w$gene
  m
}

make_heatmap <- function(metric) {
  mat <- to_matrix(metric)
  # Order columns by model then experiment so the raw matrix is tidy; clustering
  # (below) reorders the display, but this keeps CSVs readable.
  ord_setups <- setup_meta %>%
    arrange(match(model, names(MODEL_COLORS)), exp) %>% pull(setup)
  ord_setups <- intersect(ord_setups, colnames(mat))
  mat <- mat[, ord_setups, drop = FALSE]

  # Spearman correlation across setups, pairwise-complete over genes.
  cormat <- cor(mat, method = COR_METHOD, use = "pairwise.complete.obs")
  write.csv(round(cormat, 4),
            file.path(OUTPUT_DIR, sprintf("%s_spearman_matrix.csv", metric)))

  # Report overlap counts (how many genes each pair shared) — provenance.
  n_pair <- crossprod(!is.na(mat))          # setups x setups shared-gene counts
  write.csv(n_pair, file.path(OUTPUT_DIR,
            sprintf("%s_pairwise_gene_overlap.csv", metric)))

  # Console summary: within-model vs cross-model mean concordance.
  mdl <- setup_meta$model[match(colnames(cormat), setup_meta$setup)]
  same <- outer(mdl, mdl, `==`); diag(same) <- NA
  offdiag <- cormat; diag(offdiag) <- NA
  cat(sprintf("\n[%s] mean Spearman rho — within-model: %.3f | cross-model: %.3f | overall off-diag: %.3f\n",
              metric,
              mean(offdiag[same], na.rm = TRUE),
              mean(offdiag[!same], na.rm = TRUE),
              mean(offdiag, na.rm = TRUE)))

  # Clustering on 1 - rho (NA distances -> 1 for clustering only).
  dmat <- 1 - cormat; dmat[is.na(dmat)] <- 1
  hc <- hclust(as.dist(dmat), method = "average")
  dend <- as.dendrogram(hc)

  # Color scale: adapt lower bound to the data (rho is in [-1, 1]; realistically
  # [~0, 1] here). White pinned at the midpoint between the floor and 1.
  lo <- min(0, floor(min(cormat, na.rm = TRUE) * 10) / 10)
  col_fun <- colorRamp2(c(lo, (lo + 1) / 2, 1),
                        c("#3B4CC0", "#F7F7F7", "#B40426"))

  # Model + GAM_K annotation (top + left).
  ann_models <- setup_meta$model[match(colnames(cormat), setup_meta$setup)]
  ann_k      <- setup_meta$gam_k[match(colnames(cormat), setup_meta$setup)]
  top_ann <- HeatmapAnnotation(
    model = ann_models, GAM_K = ann_k,
    col = list(model = MODEL_COLORS, GAM_K = K_COLORS),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_legend_param = list(model = list(title = "GAM model"),
                                   GAM_K = list(title = "GAM_K (spline k)")))
  left_ann <- rowAnnotation(
    model = ann_models, GAM_K = ann_k,
    col = list(model = MODEL_COLORS, GAM_K = K_COLORS),
    show_annotation_name = FALSE, show_legend = FALSE)

  lab_fs  <- if (n_setups <= 30) 8 else 6
  cell_fs <- if (n_setups <= 30) 6 else 4

  cell_fun <- if (SHOW_CELL_VALUES) function(j, i, x, y, w, h, fill) {
    v <- cormat[i, j]
    if (!is.na(v)) grid.text(sprintf("%.2f", v), x, y,
                             gp = gpar(fontsize = cell_fs,
                                       col = if (v > (lo + 1) / 2 + 0.35 ||
                                                 v < (lo + 1) / 2 - 0.35) "white" else "black"))
  } else NULL

  metric_lab <- if (metric == "p_shape") "trajectory-shape test (p_shape)" else
                                          "level test (p_level)"
  ht <- Heatmap(cormat, name = "Spearman\nrho",
                col = col_fun, na_col = "grey85",
                cluster_rows = dend, cluster_columns = dend,
                row_dend_width = unit(1.4, "cm"), column_dend_height = unit(1.4, "cm"),
                top_annotation = top_ann, left_annotation = left_ann,
                row_names_gp = gpar(fontsize = lab_fs),
                column_names_gp = gpar(fontsize = lab_fs),
                column_names_rot = 45,
                cell_fun = cell_fun,
                column_title = sprintf(
                  "Rank concordance across %d setups — %s\n(Spearman of per-gene p-value ranks; CD8+ T)",
                  ncol(cormat), metric_lab),
                column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                heatmap_legend_param = list(at = c(lo, (lo + 1) / 2, 1)))

  sz <- max(9, ncol(cormat) * 0.32 + 4.5)
  pdf(file.path(OUTPUT_DIR, sprintf("%s_spearman_heatmap.pdf", metric)),
      width = sz, height = sz)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat(sprintf("[%s] wrote %s_spearman_heatmap.pdf (%d x %d setups)\n",
              metric, metric, ncol(cormat), ncol(cormat)))
  invisible(cormat)
}

for (m in METRICS) make_heatmap(m)

cat("\nDone. Outputs saved to:", OUTPUT_DIR, "\n")
