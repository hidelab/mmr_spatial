# 16_04_loso_influence.R
# ============================================================================
# Leave-one-sample-out (LOSO) influence analysis for the 15_04 tumor-trajectory
# shape test.
#
# ---------------------------------------------------------------------------
# WHY THIS EXISTS
# ---------------------------------------------------------------------------
# The sample-level permutation (16_03) answers "is the observed WT/DKO split
# special vs. arbitrary regroupings of these 7 mice?" (global p ~ 0.20 for
# CD8+ T). It is the correct POPULATION-level inference, but it has two known
# limitations at n = 3 WT + 4 DKO:
#   - a hard floor of 1/choose(7,3) = 1/35 ~= 0.029 on any p-value, and
#   - with only 35 discrete regroupings, a single leverage mouse (the WT-like
#     DKO2, or the tumor-poor DKO1/DKO3) can dominate both the observed
#     statistic AND many permuted ones.
#
# The DIRECT question "does one sample drive the findings?" is not a question
# about the null distribution — it is a question about INFLUENCE. The honest
# tool for it is leave-one-sample-out: drop each mouse in turn (keeping the true
# WT/DKO labels), rerun the whole pipeline on the remaining 6 mice, and see which
# shape-significant genes survive and which evaporate. A result that hinges on a
# single mouse is not robust regardless of what any p-value says.
#
# This uses the SAME stable pseudobulk_fs model and local-fraction score as
# 15_04 / 16_03, so the full-fit baseline reproduces the reported gene set.
#
# WHY NOT SHUFFLE CELL POSITIONS ALONG THE TRAJECTORY (an alternative null we
# considered and rejected): shuffling each cell's score moves the exchangeable
# unit back to the CELL, which (a) reintroduces the pseudoreplication the whole
# analysis was built to remove — the null's variability then comes from
# within-tissue cell resampling, not from between-mouse trajectory-shape
# heterogeneity, so it is anticonservative; and (b) tests the wrong null ("no
# trajectory in either condition" rather than "the two conditions share the same
# trajectory"). LOSO addresses single-sample leverage without either flaw.
#
# READING THE OUTPUT: the asymptotic mgcv FDRs are miscalibrated in MAGNITUDE at
# this sample size (see 16_03) — LOSO reads the STABILITY of gene membership
# across the 7 drops, not calibrated significance.
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

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_loso_influence")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters — MUST match 16_03 / 15_04 so this refits the SAME model / score
# ============================================================================
CELLTYPES        <- c("CD8+ T")
TARGET_CELLTYPE  <- "Tumor"

# Score (local tumor fraction, no per-sample min-max)
KNN_K            <- 30
DIFFUSION_VARIANT<- "heat"
DIFFUSION_HEAT_STEPS <- 20
DIFFUSION_HEAT_ALPHA <- 2
LOCAL_FRACTION   <- TRUE

# Pseudobulk / GAM
N_BINS            <- 10
MIN_CELLS_PER_BIN <- 5
GAM_K             <- 5
FS_K              <- 3
SIG_FDR           <- 0.1
COMMON_SUPPORT               <- TRUE
MIN_SAMPLES_PER_COND_PER_BIN <- 1
PB_SCORE_AXIS                <- "score"   # "score" | "decile"

# How many top genes (by full-fit p_shape) to show in the influence heatmap.
HEATMAP_TOP_N     <- 50

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

cat("Computing local-tumor-fraction score...\n")
tumor_score <- compute_tumor_score(
  all_coords, is_tumor, samples, method = "diffusion",
  diffusion_method = DIFFUSION_VARIANT, k = KNN_K, heat_steps = DIFFUSION_HEAT_STEPS,
  alpha = DIFFUSION_HEAT_ALPHA,
  local_fraction = LOCAL_FRACTION, normalize = !LOCAL_FRACTION, verbose = FALSE)

# ============================================================================
# Core fit: the SAME pseudobulk_fs difference-smooth model as 15_04 / 16_03.
# Labels are the TRUE WT/DKO labels here (LOSO does not permute); `wt_set` is
# passed so common support is computed against the retained samples' labels.
# ============================================================================
fit_shape <- function(y_by_group, g_score, g_sample, g_bin, g_n, wt_set) {
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

# ----------------------------------------------------------------------------
# Build the sample x bin pseudobulk structure from a chosen set of cells.
# Binning quantiles + common support are recomputed on the RETAINED cells, so
# a LOSO subset behaves exactly as if that mouse were never in the study.
# Returns group indices (into the ORIGINAL ct-cell vector) so a genes x ct-cell
# expression matrix can be aggregated without re-subsetting.
# ----------------------------------------------------------------------------
build_pb <- function(ct_idx_local, ct_score, ct_sample, keep_cells) {
  sc  <- ct_score[keep_cells]
  sm  <- droplevels(ct_sample[keep_cells])
  br  <- quantile(sc, seq(0, 1, length.out = N_BINS + 1), na.rm = TRUE)
  br[1] <- br[1] - 1e-6
  # guard against duplicate breakpoints (possible in a sparse subset)
  if (any(duplicated(br))) br <- unique(br)
  bn  <- as.integer(cut(sc, breaks = br, labels = seq_len(length(br) - 1)))
  grp <- paste(as.character(sm), bn, sep = "__")
  # map each retained cell (original ct-local index) to its sample x bin group
  gi  <- split(keep_cells, grp)
  gi  <- gi[vapply(gi, length, 0L) >= MIN_CELLS_PER_BIN]
  if (!length(gi)) return(NULL)
  g_sample <- vapply(strsplit(names(gi), "__"), `[`, "", 1)
  g_bin    <- as.integer(vapply(strsplit(names(gi), "__"), `[`, "", 2))
  g_n      <- vapply(gi, length, 0L)
  g_score  <- if (PB_SCORE_AXIS == "score")
                vapply(gi, function(ix) mean(ct_score[ix]), 0.0)
              else (g_bin - 0.5) / N_BINS
  list(gi = gi, g_sample = g_sample, g_bin = g_bin, g_n = g_n, g_score = g_score)
}

# Fit all genes for one pseudobulk structure; returns data.frame(gene, F, p, edf, fdr)
fit_all_genes <- function(pb, expr, test_genes, wt_set) {
  gmean <- function(g) vapply(pb$gi, function(ix) mean(expr[g, ix]), 0.0)
  m <- t(vapply(test_genes, function(g)
           fit_shape(gmean(g), pb$g_score, pb$g_sample, pb$g_bin, pb$g_n, wt_set),
           numeric(3)))
  df <- data.frame(gene = test_genes, F = m[, 1], p = m[, 2], edf_fs = m[, 3],
                   stringsAsFactors = FALSE)
  df$fdr <- p.adjust(df$p, method = "BH")
  df
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

  # --- test genes: SAME baseline-expression filter as 15_04 / 16_03 ----------
  det_frac  <- Matrix::rowSums(counts > 0) / length(ct_idx)
  samp_sum  <- vapply(levels(ct_sample), function(s)
                 Matrix::rowSums(counts[, ct_sample == s, drop = FALSE]),
                 numeric(nrow(counts)))
  n_samp_det <- rowSums(samp_sum >= MIN_COUNTS_PER_SAMPLE)
  test_genes <- rownames(counts)[det_frac >= MIN_GENE_EXPR_FRAC &
                                 n_samp_det >= MIN_SAMPLES_EXPRESSED]
  test_genes <- intersect(test_genes, rownames(expr))
  cat(sprintf("  %d tested genes (%d cells across %d samples)\n",
              length(test_genes), length(ct_idx), nlevels(ct_sample)))

  samps   <- levels(ct_sample)
  # ct_idx is index into full seu; build_pb works in ct-cell space (1..length(ct_idx)),
  # gene_group_mean indexes expr (genes x ct cells) with those same local indices.
  ct_local <- seq_along(ct_idx)
  wt_here  <- intersect(WT_SAMPLES, samps)

  # ==========================================================================
  # FULL fit (all samples) — the baseline gene set
  # ==========================================================================
  pb_full  <- build_pb(ct_local, ct_score, ct_sample, ct_local)
  full     <- fit_all_genes(pb_full, expr, test_genes, wt_here)
  full     <- full[order(full$p), ]
  sig_full <- full$gene[which(full$fdr < SIG_FDR)]
  cat(sprintf("  FULL: %d pb points; %d genes at asymptotic FDR < %.2f\n",
              length(pb_full$gi), length(sig_full), SIG_FDR))

  # ==========================================================================
  # LOSO: drop each sample, refit everything on the remaining 6
  # ==========================================================================
  loso <- list()
  drop_counts <- integer(length(samps)); names(drop_counts) <- samps
  for (s in samps) {
    keep_cells <- ct_local[ct_sample != s]
    pb_s <- build_pb(ct_local, ct_score, ct_sample, keep_cells)
    if (is.null(pb_s)) { cat(sprintf("  drop %s: no pb points; skip\n", s)); next }
    wt_s <- setdiff(wt_here, s)               # retained WT labels
    fit_s <- fit_all_genes(pb_s, expr, test_genes, wt_s)
    rownames(fit_s) <- fit_s$gene
    loso[[s]] <- fit_s
    drop_counts[s] <- sum(fit_s$fdr < SIG_FDR, na.rm = TRUE)
    cat(sprintf("  drop %-5s (%s): %d pb points, %d genes FDR<%.2f\n",
                s, if (s %in% wt_here) "WT" else "DKO",
                length(pb_s$gi), drop_counts[s], SIG_FDR))
  }

  # ==========================================================================
  # Influence table: per gene, FDR in full + each drop + robustness count
  # ==========================================================================
  gene_order <- full$gene
  tab <- data.frame(gene = gene_order,
                    p_shape_full = full$p[match(gene_order, full$gene)],
                    fdr_full     = full$fdr[match(gene_order, full$gene)],
                    stringsAsFactors = FALSE)
  for (s in names(loso)) {
    tab[[paste0("fdr_drop_", s)]] <- loso[[s]][gene_order, "fdr"]
  }
  drop_cols <- grep("^fdr_drop_", names(tab), value = TRUE)
  # significant in how many of the LOSO subsets
  tab$n_subsets_sig <- rowSums(tab[, drop_cols, drop = FALSE] < SIG_FDR, na.rm = TRUE)
  tab$n_subsets     <- length(drop_cols)
  tab$max_fdr_over_drops <- apply(tab[, drop_cols, drop = FALSE], 1, max, na.rm = TRUE)
  tab$sig_full <- tab$fdr_full < SIG_FDR
  # a gene is "robust" if it is significant in the full fit AND in every LOSO drop
  tab$robust <- tab$sig_full & (tab$n_subsets_sig == tab$n_subsets)
  # which single dropped sample most often breaks a full-fit-significant gene
  write.csv(tab, file.path(ct_dir, "loso_influence_table.csv"), row.names = FALSE)

  # ==========================================================================
  # Summary
  # ==========================================================================
  n_robust <- sum(tab$robust, na.rm = TRUE)
  # for each dropped sample: how many of the full-fit sig genes it kills
  kills <- vapply(names(loso), function(s) {
    sum(tab$sig_full & !(loso[[s]][gene_order, "fdr"] < SIG_FDR), na.rm = TRUE)
  }, integer(1))
  most_influential <- names(kills)[which.max(kills)]

  summ <- data.frame(
    cell_type = ct,
    n_tested = length(test_genes),
    n_sig_full = length(sig_full),
    n_robust_all_drops = n_robust,
    frac_robust = round(n_robust / max(1, length(sig_full)), 3),
    most_influential_sample = most_influential,
    genes_killed_by_most_influential = max(kills),
    min_count_over_drops = min(drop_counts[names(loso)]),
    max_count_over_drops = max(drop_counts[names(loso)]),
    stringsAsFactors = FALSE)
  write.csv(summ, file.path(ct_dir, "loso_summary.csv"), row.names = FALSE)

  # per-drop count table
  cnt_df <- data.frame(dropped_sample = c("(none: full)", names(drop_counts)),
                       condition = c("-", ifelse(names(drop_counts) %in% wt_here, "WT", "DKO")),
                       n_genes_fdr_sig = c(length(sig_full), as.integer(drop_counts)),
                       genes_killed_vs_full = c(0L, as.integer(kills[names(drop_counts)])))
  write.csv(cnt_df, file.path(ct_dir, "loso_counts_by_dropped_sample.csv"), row.names = FALSE)

  cat(sprintf("\n  SUMMARY: %d/%d full-fit shape genes survive dropping ANY single mouse (%.0f%%)\n",
              n_robust, length(sig_full), 100 * n_robust / max(1, length(sig_full))))
  cat(sprintf("  Most influential mouse: %s (removing it drops %d of the %d full-fit genes)\n",
              most_influential, max(kills), length(sig_full)))
  cat("  Per-drop counts:\n"); print(cnt_df, row.names = FALSE)

  # ==========================================================================
  # Plots
  # ==========================================================================
  # (1) shape-sig gene count per dropped sample
  cnt_df$dropped_sample <- factor(cnt_df$dropped_sample, levels = cnt_df$dropped_sample)
  p_cnt <- ggplot(cnt_df, aes(x = dropped_sample, y = n_genes_fdr_sig, fill = condition)) +
    geom_col() +
    geom_hline(yintercept = length(sig_full), linetype = "dashed", color = "grey40") +
    geom_text(aes(label = n_genes_fdr_sig), vjust = -0.4, size = 3) +
    labs(x = "dropped sample", y = sprintf("# genes at asymptotic FDR < %.2f", SIG_FDR),
         title = sprintf("%s — LOSO: shape-significant gene count", ct),
         subtitle = "dashed = full-fit count; large drop => that mouse drives the result") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(ct_dir, "loso_gene_counts.pdf"), p_cnt, width = 6.5, height = 4.2)

  # (2) influence heatmap: top genes x {full, each drop}, cell = -log10 FDR, X = not sig
  top_genes <- head(full$gene[which(full$fdr < SIG_FDR)], HEATMAP_TOP_N)
  if (length(top_genes) >= 2) {
    hm_cols <- c("fdr_full", drop_cols)
    M <- as.matrix(tab[match(top_genes, tab$gene), hm_cols])
    rownames(M) <- top_genes
    colnames(M) <- c("full", sub("^fdr_drop_", "-", drop_cols))
    long <- expand.grid(gene = rownames(M), drop = colnames(M),
                        stringsAsFactors = FALSE)
    long$fdr <- M[cbind(match(long$gene, rownames(M)), match(long$drop, colnames(M)))]
    long$neglog <- pmin(-log10(pmax(long$fdr, .Machine$double.eps)), 6)  # cap at 6
    long$sig <- long$fdr < SIG_FDR
    long$gene <- factor(long$gene, levels = rev(top_genes))
    long$drop <- factor(long$drop, levels = colnames(M))
    p_hm <- ggplot(long, aes(x = drop, y = gene, fill = neglog)) +
      geom_tile(color = "white") +
      geom_text(data = subset(long, !sig), aes(label = "x"), size = 2.4, color = "grey30") +
      scale_fill_gradient(low = "grey92", high = CONDITION_COLORS[[TEST_CONDITION]],
                          name = "-log10 FDR\n(capped 6)") +
      labs(x = "leave-one-out (full, then drop each sample)", y = NULL,
           title = sprintf("%s — LOSO influence on shape genes", ct),
           subtitle = "x = not significant when that mouse is dropped") +
      theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 6))
    ggsave(file.path(ct_dir, "loso_influence_heatmap.pdf"), p_hm,
           width = 7, height = max(4, 0.16 * length(top_genes) + 1.5), limitsize = FALSE)
  }

  cat(sprintf("  Done: %s -> %s\n", ct, ct_dir))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")
