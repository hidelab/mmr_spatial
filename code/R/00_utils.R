# ==============================================================================
# 00_utils.R — Reusable tumor-proximity scoring functions
# ==============================================================================
# Source this file at the top of trajectory scripts:
#   source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))
#
# Provides three ways to define a per-cell "tumor proximity / density" score,
# refactored from 15_01 (k-NN purity), 15_02 (KDE) and 16_01 (graph diffusion),
# plus a single dispatcher `compute_tumor_score()` that loops per sample and
# normalizes. Each low-level scorer operates on ONE sample's coordinates.
#
# The sourcing script is expected to load the required libraries:
#   FNN     (get.knn)        -> score_knn_purity, score_diffusion
#   MASS    (kde2d)          -> score_kde
#   Matrix  (sparseMatrix..) -> score_diffusion
# ==============================================================================

# ------------------------------------------------------------------------------
# Diffusion helpers (moved from 16_01)
# ------------------------------------------------------------------------------

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
# Returns the (unscaled) diffused signal. Internal helper for score_diffusion().
.diffuse_field <- function(W, s_raw, method = "implicit",
                           alpha = 1.0, heat_steps = 20, ppr_beta = 0.95,
                           use_cg = FALSE) {
  n <- nrow(W)
  d_vec <- Matrix::rowSums(W)
  d_vec[d_vec == 0] <- 1  # isolated nodes

  if (method == "ppr") {
    # Random-walk personalized PageRank (asymmetric -> direct solve only)
    P <- Matrix::Diagonal(n, x = 1 / d_vec) %*% W        # row-stochastic transition
    A <- Matrix::Diagonal(n) - ppr_beta * P
    return(as.numeric(Matrix::solve(A, (1 - ppr_beta) * s_raw)))
  }

  # Normalized symmetric Laplacian: L_sym = I - D^{-1/2} W D^{-1/2}
  D_inv_sqrt <- Matrix::Diagonal(n, x = 1 / sqrt(d_vec))
  L_sym <- Matrix::Diagonal(n) - D_inv_sqrt %*% W %*% D_inv_sqrt

  if (method == "heat") {
    # Backward-Euler heat kernel: repeat the implicit step heat_steps times.
    # Each step is unconditionally stable; total diffusion time = heat_steps * alpha.
    A <- Matrix::Diagonal(n) + alpha * L_sym
    s <- s_raw
    if (use_cg) {
      for (step in seq_len(heat_steps)) s <- cg_solve(A, s)
    } else {
      ch <- Matrix::Cholesky(as(A, "symmetricMatrix"), LDL = FALSE, perm = TRUE)
      for (step in seq_len(heat_steps)) s <- as.numeric(Matrix::solve(ch, s, system = "A"))
    }
    return(s)
  }

  # method == "implicit": single backward-Euler step
  A <- Matrix::Diagonal(n) + alpha * L_sym
  if (use_cg) cg_solve(A, s_raw) else as.numeric(Matrix::solve(A, s_raw))
}

# ------------------------------------------------------------------------------
# Low-level scorers — each takes ONE sample's coords + is_tumor, returns a
# raw (unnormalized) per-cell score vector of length nrow(coords).
# ------------------------------------------------------------------------------

# k-NN purity: fraction of a cell's k nearest spatial neighbors that are tumor.
# (from 15_01)
score_knn_purity <- function(coords, is_tumor, k = 30) {
  n <- nrow(coords)
  if (n < k + 1) return(rep(NA_real_, n))
  is_tumor <- as.integer(is_tumor)
  knn_res <- FNN::get.knn(coords, k = k)
  neighbor_is_tumor <- matrix(is_tumor[knn_res$nn.index], nrow = n, ncol = k)
  rowMeans(neighbor_is_tumor)
}

# KDE: 2D kernel density of tumor-cell positions, bilinearly interpolated to
# every cell. Larger bw_scale -> smoother / more spread-out field. (from 15_02)
score_kde <- function(coords, is_tumor, bw_scale = 2.0, grid_n = 200, bw = NULL) {
  n <- nrow(coords)
  tumor_idx <- which(as.logical(is_tumor))
  if (length(tumor_idx) < 10) return(rep(NA_real_, n))

  tumor_coords <- coords[tumor_idx, , drop = FALSE]

  # Estimate bandwidth from tumor cells (or use provided), then scale it.
  if (is.null(bw)) {
    bw <- c(MASS::bandwidth.nrd(tumor_coords[, 1]),
            MASS::bandwidth.nrd(tumor_coords[, 2]))
  }
  bw <- bw * bw_scale

  x_range <- range(coords[, 1])
  y_range <- range(coords[, 2])

  kde_fit <- MASS::kde2d(tumor_coords[, 1], tumor_coords[, 2],
                         h = bw, n = grid_n, lims = c(x_range, y_range))
  x_grid <- kde_fit$x
  y_grid <- kde_fit$y
  z_grid <- kde_fit$z

  # Bilinear interpolation of the KDE grid at each cell's position.
  cell_density <- numeric(n)
  for (i in seq_len(n)) {
    cx <- coords[i, 1]
    cy <- coords[i, 2]
    xi <- findInterval(cx, x_grid, all.inside = TRUE)
    yi <- findInterval(cy, y_grid, all.inside = TRUE)

    x1 <- x_grid[xi]; x2 <- x_grid[xi + 1]
    y1 <- y_grid[yi]; y2 <- y_grid[yi + 1]

    fQ11 <- z_grid[xi, yi]
    fQ21 <- z_grid[xi + 1, yi]
    fQ12 <- z_grid[xi, yi + 1]
    fQ22 <- z_grid[xi + 1, yi + 1]

    denom <- (x2 - x1) * (y2 - y1)
    cell_density[i] <- (fQ11 * (x2 - cx) * (y2 - cy) +
                        fQ21 * (cx - x1) * (y2 - cy) +
                        fQ12 * (x2 - cx) * (cy - y1) +
                        fQ22 * (cx - x1) * (cy - y1)) / denom
  }
  cell_density
}

# Graph diffusion: diffuse a binary tumor indicator over a Gaussian-weighted
# symmetric spatial k-NN graph. method one of "implicit", "heat", "ppr".
# (from 16_01)
score_diffusion <- function(coords, is_tumor, method = "ppr", k = 30,
                            alpha = 2.0, heat_steps = 20, ppr_beta = 0.95,
                            use_cg = FALSE) {
  stopifnot(method %in% c("implicit", "heat", "ppr"))
  n <- nrow(coords)
  if (n < k + 1) return(rep(NA_real_, n))

  knn_res <- FNN::get.knn(coords, k = k)

  # Symmetric adjacency with distance-based (Gaussian kernel) weights.
  sigma   <- median(knn_res$nn.dist)
  rows    <- rep(seq_len(n), each = k)
  cols    <- as.vector(t(knn_res$nn.index))
  dists   <- as.vector(t(knn_res$nn.dist))
  weights <- exp(-dists^2 / (2 * sigma^2))

  W <- Matrix::sparseMatrix(i = rows, j = cols, x = weights, dims = c(n, n))
  W <- (W + Matrix::t(W)) / 2  # symmetrize

  .diffuse_field(W, as.numeric(as.integer(is_tumor)), method = method,
                 alpha = alpha, heat_steps = heat_steps, ppr_beta = ppr_beta,
                 use_cg = use_cg)
}

# ------------------------------------------------------------------------------
# Dispatcher — single entry point. Loops per sample, calls the chosen scorer
# with `...`, min-max normalizes to [0, 1] within each sample, and returns a
# full-length vector aligned to `samples`.
#
#   compute_tumor_score(coords, is_tumor, samples, method = "kde", bw_scale = 2)
#   compute_tumor_score(coords, is_tumor, samples, method = "knn", k = 30)
#   compute_tumor_score(coords, is_tumor, samples, method = "diffusion",
#                       diffusion_method = "ppr", ppr_beta = 0.95)
# ------------------------------------------------------------------------------
compute_tumor_score <- function(coords, is_tumor, samples,
                                method = c("kde", "knn", "diffusion"),
                                normalize = TRUE, verbose = TRUE, ...) {
  method <- match.arg(method)
  n_all  <- nrow(coords)
  score  <- rep(NA_real_, n_all)
  dots   <- list(...)

  # For "diffusion", the internal graph-diffusion variant is passed as
  # `diffusion_method` to avoid clashing with the top-level `method`.
  samp_levels <- if (is.factor(samples)) levels(samples) else unique(samples)

  for (s in samp_levels) {
    idx <- which(samples == s)
    if (length(idx) == 0) next
    coords_s   <- coords[idx, , drop = FALSE]
    is_tumor_s <- is_tumor[idx]

    raw <- switch(method,
      knn = do.call(score_knn_purity, c(list(coords_s, is_tumor_s), dots)),
      kde = do.call(score_kde,        c(list(coords_s, is_tumor_s), dots)),
      diffusion = {
        d <- dots
        if (!is.null(d$diffusion_method)) {
          d$method <- d$diffusion_method
          d$diffusion_method <- NULL
        }
        do.call(score_diffusion, c(list(coords_s, is_tumor_s), d))
      }
    )

    if (normalize && all(is.finite(raw))) {
      rng <- range(raw)
      raw <- if (diff(rng) > 0) (raw - rng[1]) / diff(rng) else rep(0.5, length(raw))
    }

    score[idx] <- raw
    if (verbose) {
      cat(sprintf("  %s: %d cells (tumor frac: %.2f)\n",
                  s, length(idx), mean(as.numeric(as.integer(is_tumor_s)))))
    }
  }
  score
}
