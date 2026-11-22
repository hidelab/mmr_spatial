# 16_01_tumor_contours.R
# Tumor proximity via graph-Laplacian diffusion: compute a smooth, continuous
# tumor-proximity score for every cell by diffusing a binary tumor indicator
# over the spatial k-NN graph. Then fit GAM trajectories and draw spatial
# contour maps of gene expression along this gradient.

library(Seurat)
library(ggplot2)
library(patchwork)
library(FNN)
library(Matrix)
library(dplyr)
library(tidyr)
library(mgcv)
library(viridis)
library(ggnewscale)
library(akima)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_tumor_contours")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters
# =============================================================================
# --- Diffusion method -------------------------------------------------------
# "implicit" : single backward-Euler step, (I + alpha*L) s = s_raw  (mild low-pass)
# "heat"     : HEAT_STEPS backward-Euler steps -> long total diffusion time
#              (Option 1: genuine long-range spread via heat kernel)
# "ppr"      : personalized PageRank / random-walk propagation
#              s = (1-beta)(I - beta*P)^{-1} s_raw, P = D^{-1} W
#              (Option 2: tunable decay length; beta -> 1 spreads farther)
DIFFUSION_METHOD  <- "heat"   # one of "implicit", "heat", "ppr"

K_NEIGHBORS       <- 30
ALPHA_DIFFUSION   <- 2.0     # diffusion strength / step size on normalized Laplacian (eigenvalues in [0,2])
HEAT_STEPS        <- 5      # number of backward-Euler steps for method = "heat"
PPR_BETA          <- 0.95    # random-walk restart complement for method = "ppr" (0<beta<1; ->1 = more spread)
USE_CG_SOLVER     <- TRUE    # TRUE = conjugate gradient (scalable); FALSE = sparse Cholesky (implicit/heat only)
TARGET_CELLTYPE   <- "Tumor"
TOP_N_GENES       <- 20
GAM_K             <- 5
N_CONTOUR_BINS    <- 50      # grid resolution for contour interpolation
CELLTYPES_TRAJECTORY <- IMMUNE_CELLTYPES_DE

stopifnot(DIFFUSION_METHOD %in% c("implicit", "heat", "ppr"))

# =============================================================================
# Diffusion helpers
# =============================================================================
# Conjugate gradient for symmetric positive-definite A (used by implicit/heat).
cg_solve <- function(A, b, tol = 1e-8, maxit = 200) {
  x <- b
  r <- b - as.numeric(A %*% x)
  p <- r
  rsold <- sum(r * r)
  for (iter in seq_len(maxit)) {
    Ap <- as.numeric(A %*% p)
    a  <- rsold / sum(p * Ap)
    x  <- x + a * p
    r  <- r - a * Ap
    rsnew <- sum(r * r)
    if (sqrt(rsnew) < tol) break
    p <- r + (rsnew / rsold) * p
    rsold <- rsnew
  }
  x
}

# Compute a smoothed field from a raw signal over a weighted adjacency W.
# Returns the (unscaled) diffused signal.
compute_diffusion <- function(W, s_raw, method = "implicit",
                              alpha = 1.0, heat_steps = 20, ppr_beta = 0.95,
                              use_cg = FALSE) {
  n <- nrow(W)
  d_vec <- rowSums(W)
  d_vec[d_vec == 0] <- 1  # isolated nodes

  if (method == "ppr") {
    # Random-walk personalized PageRank (asymmetric -> direct solve only)
    P <- Diagonal(n, x = 1 / d_vec) %*% W          # row-stochastic transition
    A <- Diagonal(n) - ppr_beta * P
    return(as.numeric(solve(A, (1 - ppr_beta) * s_raw)))
  }

  # Normalized symmetric Laplacian: L_sym = I - D^{-1/2} W D^{-1/2}
  D_inv_sqrt <- Diagonal(n, x = 1 / sqrt(d_vec))
  L_sym <- Diagonal(n) - D_inv_sqrt %*% W %*% D_inv_sqrt

  if (method == "heat") {
    # Backward-Euler heat kernel: repeat the implicit step heat_steps times.
    # Each step is unconditionally stable; total diffusion time = heat_steps * alpha.
    A <- Diagonal(n) + alpha * L_sym
    s <- s_raw
    if (use_cg) {
      for (step in seq_len(heat_steps)) s <- cg_solve(A, s)
    } else {
      ch <- Cholesky(as(A, "symmetricMatrix"), LDL = FALSE, perm = TRUE)
      for (step in seq_len(heat_steps)) s <- as.numeric(solve(ch, s, system = "A"))
    }
    return(s)
  }

  # method == "implicit": single backward-Euler step
  A <- Diagonal(n) + alpha * L_sym
  if (use_cg) cg_solve(A, s_raw) else as.numeric(solve(A, s_raw))
}

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

# =============================================================================
# Build spatial k-NN graph and compute diffused tumor score (per sample)
# =============================================================================
all_coords <- Embeddings(seu, "spatial")
celltypes  <- seu$Tier1_celltype
samples    <- seu$sample
is_tumor   <- as.numeric(celltypes == TARGET_CELLTYPE)

cat(sprintf("Computing diffused tumor proximity (method=%s, k=%d, alpha=%.1f, heat_steps=%d, ppr_beta=%.2f)...\n",
            DIFFUSION_METHOD, K_NEIGHBORS, ALPHA_DIFFUSION, HEAT_STEPS, PPR_BETA))

tumor_score <- rep(NA_real_, ncol(seu))

for (s in levels(samples)) {
  idx <- which(samples == s)
  if (length(idx) < K_NEIGHBORS + 1) next

  coords_s <- all_coords[idx, , drop = FALSE]
  n <- length(idx)

  # k-NN graph

  knn_res <- get.knn(coords_s, k = K_NEIGHBORS)

  # Build symmetric adjacency with distance-based weights (Gaussian kernel)
  sigma <- median(knn_res$nn.dist)
  rows <- rep(seq_len(n), each = K_NEIGHBORS)
  cols <- as.vector(t(knn_res$nn.index))
  dists <- as.vector(t(knn_res$nn.dist))
  weights <- exp(-dists^2 / (2 * sigma^2))

  W <- sparseMatrix(i = rows, j = cols, x = weights, dims = c(n, n))
  W <- (W + t(W)) / 2  # symmetrize

  # Diffuse the binary tumor indicator over the graph (method set by DIFFUSION_METHOD)
  s_raw <- is_tumor[idx]
  s_smooth <- compute_diffusion(W, s_raw, method = DIFFUSION_METHOD,
                                alpha = ALPHA_DIFFUSION, heat_steps = HEAT_STEPS,
                                ppr_beta = PPR_BETA, use_cg = USE_CG_SOLVER)

  # Rescale to [0, 1]
  s_range <- range(s_smooth)
  if (diff(s_range) > 0) {
    s_smooth <- (s_smooth - s_range[1]) / diff(s_range)
  }

  tumor_score[idx] <- s_smooth
  cat(sprintf("  %s: %d cells (tumor frac: %.2f)\n",
              s, n, mean(is_tumor[idx])))
}

seu$tumor_score <- tumor_score
cat("Diffused tumor score computed.\n")

# =============================================================================
# Histogram of diffused tumor score per sample (diagnostic)
# =============================================================================
df_hist <- data.frame(
  score  = tumor_score[!is.na(tumor_score)],
  sample = samples[!is.na(tumor_score)]
)

p_hist <- ggplot(df_hist, aes(x = score, fill = sample)) +
  geom_histogram(bins = 60, color = NA) +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = "Diffused Tumor Score", y = "Cell count",
       title = sprintf("Tumor Score Distribution per Sample (k=%d, α=%.1f)",
                       K_NEIGHBORS, ALPHA_DIFFUSION)) +
  theme_publication()

n_samp <- length(levels(samples))
ggsave(file.path(OUTPUT_DIR, "tumor_score_histogram_per_sample.pdf"),
       p_hist, width = 4 * ceiling(n_samp / 2), height = 6)
ggsave(file.path(OUTPUT_DIR, "tumor_score_histogram_per_sample.png"),
       p_hist, width = 4 * ceiling(n_samp / 2), height = 6, dpi = 300)

# =============================================================================
# Spatial contour map of the tumor score field
# =============================================================================
sample_ids <- levels(seu$sample)
max_x_span <- max(tapply(all_coords[, 1], samples, function(x) diff(range(x))))
max_y_span <- max(tapply(all_coords[, 2], samples, function(x) diff(range(x))))
half_span  <- max(max_x_span, max_y_span) / 2 * 1.05

plots_score <- lapply(sample_ids, function(s) {
  idx <- which(samples == s)
  df_plot <- data.frame(x = all_coords[idx, 1], y = all_coords[idx, 2],
                        score = tumor_score[idx])
  ggplot(df_plot, aes(x = x, y = y, color = score)) +
    geom_point(size = 0.15, alpha = 0.6) +
    scale_color_viridis_c(option = "magma", limits = c(0, 1),
                          name = "Tumor\nScore") +
    coord_fixed(xlim = c(-half_span, half_span),
                ylim = c(-half_span, half_span)) +
    ggtitle(s) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
})

p_spatial_score <- wrap_plots(plots_score, nrow = 1) +
  plot_annotation(
    title = sprintf("Diffused Tumor Proximity Score (k=%d, α=%.0f)",
                    K_NEIGHBORS, ALPHA_DIFFUSION),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_tumor_score_field.pdf"),
       p_spatial_score, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_score_field.png"),
       p_spatial_score, width = 4 * length(sample_ids), height = 5, dpi = 300)

# =============================================================================
# Contour plot of diffused tumor score (isolines on tissue coordinates)
# =============================================================================
plots_contour <- lapply(sample_ids, function(s) {
  idx <- which(samples == s)
  df_s <- data.frame(x = all_coords[idx, 1], y = all_coords[idx, 2],
                     score = tumor_score[idx])

  grid_interp <- interp(df_s$x, df_s$y, df_s$score,
                         nx = N_CONTOUR_BINS, ny = N_CONTOUR_BINS,
                         duplicate = "mean")

  df_grid <- expand.grid(x = grid_interp$x, y = grid_interp$y)
  df_grid$z <- as.vector(grid_interp$z)
  df_grid <- df_grid[!is.na(df_grid$z), ]

  ggplot(df_grid, aes(x = x, y = y, z = z)) +
    geom_contour_filled(bins = 12) +
    scale_fill_viridis_d(option = "magma", name = "Tumor\nScore") +
    geom_contour(color = "grey20", linewidth = 0.2, alpha = 0.5, bins = 12) +
    coord_fixed(xlim = c(-half_span, half_span),
                ylim = c(-half_span, half_span)) +
    ggtitle(s) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
})

p_contour_field <- wrap_plots(plots_contour, nrow = 1) +
  plot_annotation(
    title = sprintf("Tumor Proximity Score — Contour Map (k=%d, α=%.0f)",
                    K_NEIGHBORS, ALPHA_DIFFUSION),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_tumor_score_contour.pdf"),
       p_contour_field, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_score_contour.png"),
       p_contour_field, width = 4 * length(sample_ids), height = 5, dpi = 300)

# =============================================================================
# Compare score distribution to raw k-NN purity
# =============================================================================
knn_purity <- rep(NA_real_, ncol(seu))
for (s in levels(samples)) {
  idx <- which(samples == s)
  if (length(idx) < K_NEIGHBORS + 1) next
  coords_s <- all_coords[idx, , drop = FALSE]
  knn_res <- get.knn(coords_s, k = K_NEIGHBORS)
  neighbor_tumor <- matrix(is_tumor[idx][knn_res$nn.index],
                           nrow = length(idx), ncol = K_NEIGHBORS)
  knn_purity[idx] <- rowMeans(neighbor_tumor)
}

df_compare <- data.frame(
  score = c(tumor_score[!is.na(tumor_score)], knn_purity[!is.na(knn_purity)]),
  method = rep(c("Diffusion", "k-NN Purity"),
               c(sum(!is.na(tumor_score)), sum(!is.na(knn_purity))))
)

p_compare <- ggplot(df_compare, aes(x = score, fill = method)) +

geom_histogram(bins = 80, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("Diffusion" = "#5E35B1", "k-NN Purity" = "#BDBDBD")) +
  labs(x = "Score", y = "Count",
       title = "Tumor Score Distribution: Diffusion vs k-NN Purity") +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "score_distribution_comparison.pdf"),
       p_compare, width = 7, height = 4)

# =============================================================================
# Per cell type: GAM trajectory + spatial contour plots
# =============================================================================

for (ct in CELLTYPES_TRAJECTORY) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
  ct_dir <- file.path(OUTPUT_DIR, ct_safe)
  dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf("\n=== %s ===\n", ct))

  ct_idx <- which(celltypes == ct & !is.na(tumor_score))
  if (length(ct_idx) < 50) {
    cat(sprintf("  Skipping %s: too few cells (%d)\n", ct, length(ct_idx)))
    next
  }

  ct_score     <- tumor_score[ct_idx]
  ct_condition <- seu$condition[ct_idx]
  ct_sample    <- seu$sample[ct_idx]
  ct_coords    <- all_coords[ct_idx, ]

  expr_mat <- GetAssayData(seu, assay = "RNA", layer = "data")[, ct_idx]

  cat(sprintf("  Total cells: %d (WT: %d, DKO: %d)\n",
              length(ct_idx),
              sum(ct_condition == REFERENCE_CONDITION),
              sum(ct_condition == TEST_CONDITION)))

  # ---------------------------------------------------------------------------
  # Select top variable genes along the diffused gradient (no binning)
  # ---------------------------------------------------------------------------
  # Rank genes by correlation magnitude with tumor_score
  gene_cors <- apply(expr_mat, 1, function(g) {
    if (sd(g) == 0) return(0)
    cor(g, ct_score, method = "spearman")
  })
  gene_cors <- gene_cors[!is.na(gene_cors)]
  top_genes <- names(sort(abs(gene_cors), decreasing = TRUE))[
    seq_len(min(TOP_N_GENES, length(gene_cors)))]

  # Filter genes with minimal expression

  gene_means <- rowMeans(expr_mat[top_genes, , drop = FALSE])
  top_genes <- top_genes[gene_means > 0.05]

  if (length(top_genes) < 5) {
    cat(sprintf("  Skipping %s: too few correlated genes\n", ct))
    next
  }

  # ---------------------------------------------------------------------------
  # GAM fits: WT vs DKO smooth curves along tumor score
  # ---------------------------------------------------------------------------
  gam_dir <- file.path(ct_dir, "gam_curves")
  dir.create(gam_dir, recursive = TRUE, showWarnings = FALSE)

  gam_results <- list()

  for (gene in top_genes) {
    gene_expr <- as.numeric(expr_mat[gene, ])
    df_gam <- data.frame(
      expr = gene_expr,
      score = ct_score,
      condition = ct_condition,
      sample = ct_sample
    )

    fit_gam <- tryCatch(
      gam(expr ~ s(score, by = condition, k = GAM_K) + condition,
          data = df_gam, method = "REML"),
      error = function(e) NULL
    )
    if (is.null(fit_gam)) next

    newdata_pred <- expand.grid(
      score = seq(min(ct_score), max(ct_score), length.out = 100),
      condition = levels(ct_condition)
    )
    pred <- predict(fit_gam, newdata = newdata_pred, se.fit = TRUE)
    newdata_pred$fit   <- pred$fit
    newdata_pred$se    <- pred$se.fit
    newdata_pred$ci_lo <- pred$fit - 1.96 * pred$se.fit
    newdata_pred$ci_hi <- pred$fit + 1.96 * pred$se.fit

    # Condition effect test
    fit_null <- tryCatch(
      gam(expr ~ s(score, k = GAM_K), data = df_gam, method = "REML"),
      error = function(e) NULL
    )
    anova_p <- NA
    if (!is.null(fit_null)) {
      aov_test <- anova(fit_null, fit_gam, test = "Chisq")
      anova_p <- aov_test[["Pr(>Chi)"]][2]
    }

    gam_results[[gene]] <- data.frame(
      gene = gene,
      spearman_rho = gene_cors[gene],
      p_condition = anova_p
    )

    p_gam <- ggplot(newdata_pred,
                    aes(x = score, y = fit, color = condition, fill = condition)) +
      geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi), alpha = 0.2, color = NA) +
      geom_line(linewidth = 1) +
      scale_color_manual(values = CONDITION_COLORS) +
      scale_fill_manual(values = CONDITION_COLORS) +
      labs(x = "Tumor Proximity Score", y = "Expression",
           title = sprintf("%s — %s", gene, ct),
           subtitle = ifelse(!is.na(anova_p),
                             sprintf("Condition p = %.3g | ρ = %.2f",
                                     anova_p, gene_cors[gene]), "")) +
      theme_publication()

    ggsave(file.path(gam_dir, sprintf("gam_%s.pdf", gene)),
           p_gam, width = 5, height = 4)
  }

  # Save GAM summary
  if (length(gam_results) > 0) {
    gam_summary <- bind_rows(gam_results)
    gam_summary <- gam_summary[order(gam_summary$p_condition), ]
    write.csv(gam_summary, file.path(ct_dir, "gam_condition_test.csv"),
              row.names = FALSE)
    cat(sprintf("  GAM results: %d genes\n", nrow(gam_summary)))
  }

  # ---------------------------------------------------------------------------
  # Spatial contour plots: top 6 genes overlaid on tissue coordinates
  # ---------------------------------------------------------------------------
  contour_genes <- head(top_genes, 6)
  contour_dir <- file.path(ct_dir, "spatial_contours")
  dir.create(contour_dir, recursive = TRUE, showWarnings = FALSE)

  for (gene in contour_genes) {
    gene_expr <- as.numeric(expr_mat[gene, ])

    for (s in levels(ct_sample)) {
      s_idx <- which(ct_sample == s)
      if (length(s_idx) < 20) next

      df_contour <- data.frame(
        x = ct_coords[s_idx, 1],
        y = ct_coords[s_idx, 2],
        expr = gene_expr[s_idx],
        score = ct_score[s_idx]
      )

      # Point plot colored by expression with tumor score contour lines
      p_contour <- ggplot(df_contour, aes(x = x, y = y)) +
        stat_density_2d(aes(z = score, color = after_stat(level)),
                        geom = "contour", bins = 8,
                        linewidth = 0.4, alpha = 0.6,
                        data = df_contour) +
        scale_color_gradient(low = "grey70", high = "grey30",
                             name = "Tumor\nScore\nContour") +
        ggnewscale::new_scale_color() +
        geom_point(aes(color = expr), size = 0.5, alpha = 0.7) +
        scale_color_viridis_c(option = "plasma", name = "Expr") +
        coord_fixed() +
        labs(title = sprintf("%s — %s (%s)", gene, ct, s),
             x = NULL, y = NULL) +
        theme_void() +
        theme(plot.title = element_text(size = 9, hjust = 0.5))

      ggsave(file.path(contour_dir, sprintf("contour_%s_%s.pdf", gene, s)),
             p_contour, width = 5, height = 5)
    }
  }

  # ---------------------------------------------------------------------------
  # Panel: top 4 genes, all samples tiled
  # ---------------------------------------------------------------------------
  panel_genes <- head(top_genes, 4)

  for (gene in panel_genes) {
    gene_expr <- as.numeric(expr_mat[gene, ])

    plot_list <- lapply(levels(ct_sample), function(s) {
      s_idx <- which(ct_sample == s)
      if (length(s_idx) < 10) return(NULL)

      df_p <- data.frame(x = ct_coords[s_idx, 1], y = ct_coords[s_idx, 2],
                         expr = gene_expr[s_idx])
      ggplot(df_p, aes(x = x, y = y, color = expr)) +
        geom_point(size = 0.3, alpha = 0.7) +
        scale_color_viridis_c(option = "plasma", name = "Expr") +
        coord_fixed() +
        ggtitle(s) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 8))
    })
    plot_list <- Filter(Negate(is.null), plot_list)

    if (length(plot_list) > 0) {
      p_panel <- wrap_plots(plot_list, nrow = 1) +
        plot_annotation(title = sprintf("%s in %s", gene, ct),
                        theme = theme(plot.title = element_text(
                          size = 12, face = "bold", hjust = 0.5)))
      ggsave(file.path(ct_dir, sprintf("spatial_panel_%s.pdf", gene)),
             p_panel, width = 4 * length(plot_list), height = 4.5)
    }
  }

  # ---------------------------------------------------------------------------
  # Density distribution along the gradient (WT vs DKO)
  # ---------------------------------------------------------------------------
  df_dist <- data.frame(score = ct_score, condition = ct_condition)

  p_dist <- ggplot(df_dist, aes(x = score, fill = condition, color = condition)) +
    geom_density(alpha = 0.4) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(x = "Tumor Proximity Score (diffusion)", y = "Density of cells",
         title = sprintf("%s — Position along Tumor Gradient", ct)) +
    theme_publication()

  ggsave(file.path(ct_dir, "density_distribution_along_gradient.pdf"),
         p_dist, width = 5, height = 4)

  cat(sprintf("  Done: %s\n", ct))
}

cat("\nAll done. Outputs saved to:", OUTPUT_DIR, "\n")