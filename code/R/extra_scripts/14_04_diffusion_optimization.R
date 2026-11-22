# 14_04_diffusion_optimization.R
# ============================================================================
# Grid search over the heat-diffusion tumor-score hyperparameters
# (alpha x heat_steps) to help choose a diffusion strength.
#
# Based on 15_04_tumor_trajectory_working.R + 00_utils.R (score_diffusion /
# .diffuse_field). For every (alpha, heat_steps) combination we recompute the
# per-sample tumor-proximity score field and plot it, so you can see visually
# how the two parameters spread / concentrate the field.
#
# SCORE: kernel-normalized LOCAL TUMOR FRACTION (K x)/(K 1) — the diffused tumor
# indicator divided by the diffused ones vector. This is natively in [0,1] and
# comparable ACROSS samples with NO per-sample min-max normalization, matching the
# corrected score_diffusion(local_fraction = TRUE) used in 00_utils.R / 15_04.
#
# THEORY NOTE: alpha (timestep) and heat_steps (number of backward-Euler steps)
# set the TOTAL diffusion time  T = heat_steps * alpha, i.e. the spatial length
# scale over which "tumor proximity" is measured. Larger T => smoother, longer-
# range field. The scheme is unconditionally stable for any alpha > 0, so it
# never blows up; but as T -> infinity the normalized-Laplacian heat kernel
# converges to the flat stationary (degree) distribution and the field loses all
# spatial contrast. Pick T for a principled length scale, applied GLOBALLY to
# every sample, and judge it from the field/distribution shape — NOT from any
# downstream p-value count.
#
# Efficiency: the Gaussian-weighted spatial k-NN graph depends only on the
# coordinates + k, so it is built ONCE per sample and reused across the whole
# grid; only the (cheap) diffusion solve is repeated per grid cell.
#
# Outputs (PNG ONLY -> output/R/14_diffusion_optimization/):
#   * diffusion_grid_<sample>.png    — one file PER SAMPLE, alpha (rows) x
#                                      heat_steps (cols) facet grid of the field
#   * diffusion_grid_ALL_stacked.png — every sample's grid stacked vertically
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
library(Matrix)
library(dplyr)
library(viridis)

# Centralized parameters, palette, and the diffusion scorers (.diffuse_field).
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "14_diffusion_optimization")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters (all switchable) — EDIT THE GRID HERE
# ============================================================================
# The two axes of the search. Add/remove values freely.
ALPHA_GRID    <- c(0.5, 1, 2, 4, 8)   # diffusion timestep (per step)
HEATSTEP_GRID <- c(5, 10, 20, 40)     # number of backward-Euler steps

KNN_K           <- 30        # k for the spatial k-NN graph
TARGET_CELLTYPE <- "Tumor"   # cell type whose density defines the field

# Plotting only (the diffusion itself always uses ALL cells in a sample).
N_PLOT_CELLS <- 20000        # per-sample cell subsample for display (NA = all)
PT_SIZE      <- 0.15         # point size in the spatial scatter
PLOT_DPI     <- 150          # PNG resolution
RNG_SEED     <- 1            # reproducible plot subsampling

set.seed(RNG_SEED)

# ============================================================================
# Load data
# ============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

all_coords <- Embeddings(seu, "spatial")
samples    <- seu$sample
is_tumor   <- seu$Tier1_celltype == TARGET_CELLTYPE
sample_ids <- levels(seu$sample)

grid <- expand.grid(alpha = ALPHA_GRID, heat_steps = HEATSTEP_GRID,
                    stringsAsFactors = FALSE)
cat(sprintf("Grid: %d alpha x %d heat_steps = %d combinations over %d samples\n",
            length(ALPHA_GRID), length(HEATSTEP_GRID), nrow(grid), length(sample_ids)))

# ============================================================================
# Helpers
# ============================================================================

# Build the Gaussian-weighted, symmetric spatial k-NN adjacency for one sample.
# (Same construction as score_diffusion() in 00_utils.R, factored out so it is
# computed once per sample and reused across the whole grid.)
build_graph <- function(coords, k = KNN_K) {
  n <- nrow(coords)
  knn_res <- FNN::get.knn(coords, k = k)
  sigma   <- median(knn_res$nn.dist)
  rows    <- rep(seq_len(n), each = k)
  cols    <- as.vector(t(knn_res$nn.index))
  dists   <- as.vector(t(knn_res$nn.dist))
  weights <- exp(-dists^2 / (2 * sigma^2))
  W <- Matrix::sparseMatrix(i = rows, j = cols, x = weights, dims = c(n, n))
  (W + Matrix::t(W)) / 2  # symmetrize
}

# Consistent facet labels
alpha_levels <- sort(unique(ALPHA_GRID))
step_levels  <- sort(unique(HEATSTEP_GRID))
mk_alpha_f <- function(a) factor(a, levels = alpha_levels,
                                 labels = paste0("alpha = ", alpha_levels))
mk_step_f  <- function(s) factor(s, levels = step_levels,
                                 labels = paste0("steps = ", step_levels))

# One sample's field grid as a ggplot (alpha rows x heat_steps cols).
plot_sample_grid <- function(df, title) {
  ggplot(df, aes(x = x, y = y, color = score)) +
    geom_point(size = PT_SIZE, alpha = 0.55) +
    scale_color_viridis_c(option = "magma", limits = c(0, 1),
                          name = "Local tumor\nfraction") +
    facet_grid(alpha_f ~ step_f) +
    coord_fixed() +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          strip.text = element_text(size = 8),
          legend.position = "right",
          panel.spacing = unit(0.2, "lines"))
}

# ============================================================================
# Grid search: per sample build graph once, then diffuse for every grid cell
# ============================================================================
sample_plots <- list()

for (s in sample_ids) {
  idx        <- which(samples == s)
  coords_s   <- all_coords[idx, , drop = FALSE]
  is_tumor_s <- as.numeric(as.integer(is_tumor[idx]))
  n_s        <- length(idx)
  cat(sprintf("\n[%s] %d cells (tumor frac %.2f) — building graph...\n",
              s, n_s, mean(is_tumor_s)))

  W <- build_graph(coords_s, KNN_K)

  # Cell subsample for plotting (diffusion still uses all cells).
  plot_idx <- if (is.na(N_PLOT_CELLS) || n_s <= N_PLOT_CELLS) seq_len(n_s)
              else sort(sample.int(n_s, N_PLOT_CELLS))

  ones <- rep(1, n_s)
  df_list <- vector("list", nrow(grid))
  for (g in seq_len(nrow(grid))) {
    a  <- grid$alpha[g]
    hs <- grid$heat_steps[g]
    # Kernel-normalized LOCAL TUMOR FRACTION (Kx)/(K1): the diffused tumor
    # indicator divided by the diffused ones vector. Natively in [0,1] and
    # comparable across samples with NO per-sample min-max (matches the corrected
    # score_diffusion(local_fraction=TRUE) in 00_utils.R / 15_04).
    num <- .diffuse_field(W, is_tumor_s, method = "heat", alpha = a, heat_steps = hs)
    den <- .diffuse_field(W, ones,       method = "heat", alpha = a, heat_steps = hs)
    sc  <- num / pmax(den, .Machine$double.eps)
    df_list[[g]] <- data.frame(
      x = coords_s[plot_idx, 1], y = coords_s[plot_idx, 2],
      score = sc[plot_idx], alpha = a, heat_steps = hs)
    cat(sprintf("  alpha=%-4g steps=%-3d  (T=%.1f)\n", a, hs, a * hs))
  }

  df_s <- bind_rows(df_list)
  df_s$alpha_f <- mk_alpha_f(df_s$alpha)
  df_s$step_f  <- mk_step_f(df_s$heat_steps)

  p_s <- plot_sample_grid(df_s, s)
  sample_plots[[s]] <- p_s

  # --- Save this sample separately -----------------------------------------
  ggsave(file.path(OUTPUT_DIR, sprintf("diffusion_grid_%s.png", s)),
         p_s,
         width  = max(4, length(step_levels) * 1.9),
         height = max(3.5, length(alpha_levels) * 1.9),
         dpi = PLOT_DPI, limitsize = FALSE)
}

# ============================================================================
# Stacked: every sample's grid stacked vertically in one figure
# ============================================================================
p_stack <- wrap_plots(sample_plots, ncol = 1) +
  plot_annotation(
    title = sprintf("Heat-diffusion tumor-score grid search (k=%d) — all samples",
                    KNN_K),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "diffusion_grid_ALL_stacked.png"),
       p_stack,
       width  = max(4, length(step_levels) * 1.9),
       height = max(3.5, length(alpha_levels) * 1.9) * length(sample_ids),
       dpi = PLOT_DPI, limitsize = FALSE)

cat("\nAll done. PNGs saved to:", OUTPUT_DIR, "\n")