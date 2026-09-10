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
#
# local_fraction = TRUE returns a CROSS-SAMPLE-COMPARABLE score: the diffusion
# operator K applied to the tumor indicator, divided by K applied to a ones
# vector — i.e. (K x)/(K 1), a Nadaraya-Watson kernel-weighted LOCAL TUMOR
# FRACTION. Because K is entrywise nonnegative ((I + alpha*L_sym) and (I - beta*P)
# are M-matrices, so their inverses are >= 0), (K x)/(K 1) lies in [0,1] in EVERY
# sample with NO per-sample rescaling. This removes the min-max distortion where a
# tumor-poor mouse would have its ceiling stretched to 1 (reviewer point #6). Call
# compute_tumor_score(..., normalize = FALSE) when using it, so the comparable
# scale is preserved.
score_diffusion <- function(coords, is_tumor, method = "ppr", k = 30,
                            alpha = 2.0, heat_steps = 20, ppr_beta = 0.95,
                            use_cg = FALSE, local_fraction = FALSE) {
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

  num <- .diffuse_field(W, as.numeric(as.integer(is_tumor)), method = method,
                        alpha = alpha, heat_steps = heat_steps, ppr_beta = ppr_beta,
                        use_cg = use_cg)
  if (!local_fraction) return(num)

  # Denominator: same operator applied to a ones vector -> kernel-weighted local
  # tumor fraction. (For "ppr" the row-stochastic kernel gives K 1 = 1 exactly.)
  den <- .diffuse_field(W, rep(1, n), method = method,
                        alpha = alpha, heat_steps = heat_steps, ppr_beta = ppr_beta,
                        use_cg = use_cg)
  den <- pmax(den, .Machine$double.eps)  # guard isolated / near-zero nodes
  num / den
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



# ==============================================================================
# 00_utils2.R — compute_tumor_score2(): revised tumor-proximity scoring
# ==============================================================================
# Source AFTER 00_utils.R (does not overwrite anything in it):
#   source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))
#   source(file.path(PROJECT_DIR, "code", "R", "00_utils2.R"))
#
# WHAT CHANGED vs compute_tumor_score() / score_diffusion() in 00_utils.R
# ------------------------------------------------------------------------------
# (1) RANDOM-WALK LAPLACIAN instead of the symmetric normalized Laplacian.
#     Still heat diffusion, still backward Euler -- only the operator changes.
#     L_rw = I - D^-1 W  (vs  L_sym = I - D^-1/2 W D^-1/2).
#     Why: (I + alpha*L_rw)^-1 is *exactly* row-stochastic, i.e. a weighted
#     average. Constants are preserved (diffusing an all-ones vector returns
#     all ones), so the diffused binary tumor indicator is a genuine
#     "expected tumor fraction in an alpha*steps neighborhood", bounded in
#     [0,1] with NO rescaling, and with no degree or tissue-boundary bias.
#     L_sym does not preserve constants: it is biased upward in high-degree
#     (dense) regions and downward at tissue edges / next to holes.
#     Implemented without losing the Cholesky speedup via the similarity
#     transform L_rw = D^-1/2 L_sym D^1/2, so each backward-Euler step is
#     still one sparse symmetric solve against a cached factorization.
#
# (2) GLOBAL KERNEL BANDWIDTH. The old code did sigma <- median(nn.dist)
#     *inside* each sample. Measured effect: holding tumor geometry fixed and
#     changing only cellularity 3x moved sigma from 12.68 to 7.26. That makes
#     the physical length scale of the covariate a per-sample quantity, and
#     therefore condition-confounded whenever cellularity differs by
#     condition. sigma is now computed ONCE by pooling nearest-neighbour
#     distances across all samples (or supplied directly in physical units).
#
# (3) NO MIN-MAX NORMALIZATION. The old dispatcher rescaled every sample to
#     [0,1] via (x - min)/(max - min). Measured effect: a big tumor (frac 0.50,
#     raw field IQR 0.925) and a small tumor (frac 0.13, raw field IQR 0.176)
#     both ended up spanning exactly [0,1], so tumor burden -- plausibly the
#     phenotype -- was normalized away, and both endpoints were set by single
#     extreme cells. With the rw Laplacian the score is already in [0,1] and
#     means the same thing in every sample, so no rescaling is applied.
#     (The old min-max was also applied unconditionally to method="knn", whose
#     output is already a local tumor fraction -- that was strictly lossy.)
#
# (4) OPTIONAL FIXED-RADIUS GRAPH. k-NN adjacency is density-adaptive: k
#     neighbours span a shorter physical distance in dense tissue. graph =
#     "radius" builds edges at a fixed physical radius instead, so one graph
#     step means the same distance everywhere. Kept "knn" as the default only
#     for continuity with previous runs; "radius" is the recommended setting.
#
# (5) method = "kernel_fraction": a single fixed-bandwidth Gaussian-weighted
#     local tumor fraction, (G x) / (G 1). No graph walk at all, so it is
#     density-invariant by construction. Intended as the sensitivity analysis
#     that shows a result is not an artefact of the diffusion graph.
#
# (6) FAILS LOUDLY. The old dispatcher did `if (normalize && all(is.finite(raw)))`
#     which silently left a sample on a raw scale, mixed in with [0,1] samples
#     from other mice, whenever any cell was non-finite. Now an error.
#
# (7) RETURNS DIAGNOSTICS. attr(score, "diagnostics") is a per-sample data
#     frame (tumor fraction, sigma, median nn distance, raw range, isolated
#     node count). Tabulate it against condition BEFORE trusting any
#     condition x trajectory test: if any column separates by condition, the
#     x-axis is condition-confounded and the difference smooth is not
#     interpretable.
#
# Requires: FNN, Matrix. Optional: dbscan (fast radius graph; falls back to a
# capped k-NN search if absent).
# ==============================================================================


# ------------------------------------------------------------------------------
# Graph construction
# ------------------------------------------------------------------------------

# Symmetric Gaussian-weighted spatial graph.
#   graph = "knn"    : k nearest neighbours (density-adaptive spacing)
#   graph = "radius" : all neighbours within `radius` physical units (fixed spacing)
# `sigma` is the Gaussian bandwidth in the SAME physical units as `coords`.
# self_loops guards against isolated nodes in radius graphs.
.build_spatial_graph <- function(coords, sigma, graph = c("knn", "radius"),
                                 k = 30, radius = NULL, self_loops = TRUE,
                                 knn_cap = 200) {
  graph <- match.arg(graph)
  n <- nrow(coords)

  if (graph == "knn") {
    if (n < k + 1) return(NULL)
    knn_res <- FNN::get.knn(coords, k = k)
    i <- rep(seq_len(n), each = k)
    j <- as.vector(t(knn_res$nn.index))
    d <- as.vector(t(knn_res$nn.dist))
  } else {
    if (is.null(radius)) radius <- 3 * sigma
    if (requireNamespace("dbscan", quietly = TRUE)) {
      fr <- dbscan::frNN(coords, eps = radius)
      lens <- vapply(fr$id, length, integer(1))
      i <- rep(seq_len(n), lens)
      j <- unlist(fr$id, use.names = FALSE)
      d <- unlist(fr$dist, use.names = FALSE)
    } else {
      # Fallback: capped k-NN search, then drop edges beyond the radius.
      kk <- min(knn_cap, n - 1L)
      knn_res <- FNN::get.knn(coords, k = kk)
      keep <- knn_res$nn.dist <= radius
      if (any(keep[, kk])) {
        warning(sprintf(paste("radius graph: %d cells hit the k-NN cap (%d) --",
                              "neighbourhoods are truncated. Install dbscan",
                              "for an exact radius graph."),
                        sum(keep[, kk]), kk), call. = FALSE)
      }
      i <- rep(seq_len(n), each = kk)[as.vector(t(keep))]
      j <- as.vector(t(knn_res$nn.index))[as.vector(t(keep))]
      d <- as.vector(t(knn_res$nn.dist))[as.vector(t(keep))]
    }
  }

  if (length(i) == 0) return(NULL)
  W <- Matrix::sparseMatrix(i = i, j = j, x = exp(-d^2 / (2 * sigma^2)),
                            dims = c(n, n))
  W <- (W + Matrix::t(W)) / 2                 # symmetrize
  if (self_loops) W <- W + Matrix::Diagonal(n) # avoids isolated nodes; exp(0)=1
  W
}


# ------------------------------------------------------------------------------
# Heat diffusion on the RANDOM-WALK Laplacian
# ------------------------------------------------------------------------------
# Backward Euler: s <- (I + alpha*L_rw)^-1 s, repeated `steps` times.
# Total diffusion time = alpha * steps.
#
# Each application is row-stochastic, so:
#   * the result is a weighted average of s_raw over the diffused neighbourhood
#   * s_raw in [0,1]  =>  result in [0,1]
#   * s_raw == 1 everywhere  =>  result == 1 everywhere (checked below)
# The last property is what removes the degree and boundary bias present in
# the L_sym formulation.
.diffuse_heat_rw <- function(W, s_raw, alpha = 2.0, steps = 20,
                             check_tol = 1e-6) {
  n <- nrow(W)
  d <- Matrix::rowSums(W)
  d[d == 0] <- 1

  Dh   <- sqrt(d)            # D^{1/2} as a plain vector
  Dhi  <- 1 / Dh             # D^{-1/2}

  # L_rw = D^{-1/2} L_sym D^{1/2}  =>  (I + a L_rw)^-1 = D^{-1/2} (I + a L_sym)^-1 D^{1/2}
  # so we can keep using a symmetric Cholesky factorization.
  L_sym <- Matrix::Diagonal(n) -
           Matrix::Diagonal(n, x = Dhi) %*% W %*% Matrix::Diagonal(n, x = Dhi)
  A  <- Matrix::Diagonal(n) + alpha * L_sym
  ch <- Matrix::Cholesky(as(A, "symmetricMatrix"), LDL = FALSE, perm = TRUE)

  one_step <- function(v) Dhi * as.numeric(Matrix::solve(ch, Dh * v, system = "A"))

  run <- function(v) { for (ii in seq_len(steps)) v <- one_step(v); v }

  s <- run(as.numeric(s_raw))

  # Constant-preservation check. For the rw Laplacian this is exact up to
  # solver tolerance; a failure means the graph or the factorization is wrong.
  ones <- run(rep(1, n))
  if (max(abs(ones - 1)) > check_tol) {
    warning(sprintf(paste("diffusion did not preserve constants (max dev %.2e);",
                          "dividing by the diffused all-ones vector to correct."),
                    max(abs(ones - 1))), call. = FALSE)
    s <- s / ones
  }
  s
}


# ------------------------------------------------------------------------------
# Single-step fixed-bandwidth weighted tumor fraction (density-invariant)
# ------------------------------------------------------------------------------
.kernel_fraction <- function(W, s_raw) {
  num <- as.numeric(W %*% s_raw)
  den <- Matrix::rowSums(W)
  den[den == 0] <- NA_real_
  num / den
}


# ==============================================================================
# compute_tumor_score2() — main entry point
# ==============================================================================
#   coords     : n x 2 matrix of spatial coordinates (ALL samples)
#   is_tumor   : logical/0-1 vector, length n
#   samples    : sample id per cell, length n
#   method     : "diffusion" (heat, rw Laplacian) | "kernel_fraction" | "knn"
#   graph      : "knn" | "radius"
#   sigma      : Gaussian bandwidth in physical units. NULL = compute ONE global
#                value by pooling nn distances across all samples.
#   radius     : radius-graph cutoff, physical units. NULL = 3 * sigma.
#   alpha,steps: backward-Euler step size and count. Diffusion time = alpha*steps.
#                Note the old code hard-coded alpha = 2.0 and never exposed it,
#                so the effective diffusion time (40) was invisible in configs.
#   laplacian  : "rw" (default, correct) | "sym" (legacy; ratio-corrected)
#
# Returns a numeric vector aligned to `samples`, with attr "diagnostics".
# ==============================================================================
compute_tumor_score2 <- function(coords, is_tumor, samples,
                                 method    = c("diffusion", "kernel_fraction", "knn"),
                                 graph     = c("knn", "radius"),
                                 sigma     = NULL,
                                 k         = 30,
                                 radius    = NULL,
                                 alpha     = 2.0,
                                 steps     = 20,
                                 laplacian = c("rw", "sym"),
                                 self_loops = TRUE,
                                 verbose   = TRUE) {

  method    <- match.arg(method)
  graph     <- match.arg(graph)
  laplacian <- match.arg(laplacian)

  stopifnot(nrow(coords) == length(is_tumor),
            nrow(coords) == length(samples))

  n_all       <- nrow(coords)
  is_tumor    <- as.numeric(as.integer(as.logical(is_tumor)))
  samp_levels <- if (is.factor(samples)) levels(samples) else unique(samples)
  samp_levels <- samp_levels[samp_levels %in% unique(as.character(samples))]

  # ---- Pass 1: ONE global bandwidth -----------------------------------------
  # This is the fix for the per-sample sigma <- median(knn_res$nn.dist) in
  # 00_utils.R. Per-sample nn-distance medians are also recorded so you can
  # check directly whether cellularity differs by condition.
  nnd_median <- setNames(rep(NA_real_, length(samp_levels)), samp_levels)
  nnd_pool   <- list()
  for (s in samp_levels) {
    idx <- which(as.character(samples) == s)
    if (length(idx) < k + 1) next
    dd <- FNN::get.knn(coords[idx, , drop = FALSE], k = k)$nn.dist
    nnd_median[s] <- median(dd)
    nnd_pool[[s]] <- as.vector(dd)
  }
  sigma_global <- if (!is.null(sigma)) sigma else median(unlist(nnd_pool))
  if (is.null(radius)) radius <- 3 * sigma_global

  if (verbose) {
    cat(sprintf("compute_tumor_score2: method=%s graph=%s laplacian=%s\n",
                method, graph, laplacian))
    cat(sprintf("  global sigma = %.3f (per-sample nn medians: %s)\n",
                sigma_global,
                paste(sprintf("%s=%.2f", names(nnd_median), nnd_median),
                      collapse = ", ")))
    if (graph == "radius") cat(sprintf("  radius = %.3f\n", radius))
    if (method == "diffusion")
      cat(sprintf("  diffusion time = alpha * steps = %.1f * %d = %.1f\n",
                  alpha, steps, alpha * steps))
  }

  # ---- Pass 2: score each sample on the shared bandwidth ---------------------
  score <- rep(NA_real_, n_all)
  diag_rows <- list()

  for (s in samp_levels) {
    idx      <- which(as.character(samples) == s)
    coords_s <- coords[idx, , drop = FALSE]
    tum_s    <- is_tumor[idx]

    if (method == "knn") {
      # Local tumor fraction among k nearest neighbours. Already in [0,1];
      # deliberately NOT rescaled (the old dispatcher min-maxed this).
      raw <- score_knn_purity(coords_s, tum_s, k = k)
      n_iso <- 0L
    } else {
      W <- .build_spatial_graph(coords_s, sigma = sigma_global, graph = graph,
                                k = k, radius = radius, self_loops = self_loops)
      if (is.null(W)) {
        raw <- rep(NA_real_, length(idx)); n_iso <- NA_integer_
      } else {
        n_iso <- sum(Matrix::rowSums(W) == 0)
        raw <- if (method == "kernel_fraction") {
          .kernel_fraction(W, tum_s)
        } else if (laplacian == "rw") {
          .diffuse_heat_rw(W, tum_s, alpha = alpha, steps = steps)
        } else {
          # Legacy L_sym path, ratio-corrected so it at least preserves
          # constants. Kept only for reproducing old runs.
          num <- .diffuse_field(W, tum_s, method = "heat",
                                alpha = alpha, heat_steps = steps)
          den <- .diffuse_field(W, rep(1, length(idx)), method = "heat",
                                alpha = alpha, heat_steps = steps)
          num / den
        }
      }
    }

    # ---- Fail loudly (was a silent skip in 00_utils.R) -----------------------
    if (!all(is.finite(raw))) {
      stop(sprintf("compute_tumor_score2: %d non-finite scores in sample '%s'. ",
                   sum(!is.finite(raw)), s),
           "Refusing to continue -- mixing scales across samples silently ",
           "is what the old code did.")
    }
    if (min(raw) < -1e-8 || max(raw) > 1 + 1e-8) {
      warning(sprintf("sample '%s': score outside [0,1] (range %.4f-%.4f)",
                      s, min(raw), max(raw)), call. = FALSE)
    }

    score[idx] <- raw
    diag_rows[[s]] <- data.frame(
      sample        = s,
      n_cells       = length(idx),
      tumor_frac    = mean(tum_s),
      sigma         = sigma_global,
      nn_dist_med   = unname(nnd_median[s]),
      score_min     = min(raw),
      score_max     = max(raw),
      score_median  = median(raw),
      score_iqr     = unname(diff(quantile(raw, c(0.25, 0.75)))),
      n_isolated    = n_iso,
      stringsAsFactors = FALSE)

    if (verbose) {
      cat(sprintf("  %s: n=%d tumor_frac=%.3f score[%.3f, %.3f] IQR=%.3f\n",
                  s, length(idx), mean(tum_s), min(raw), max(raw),
                  diag_rows[[s]]$score_iqr))
    }
  }

  attr(score, "diagnostics") <- do.call(rbind, diag_rows)
  attr(score, "params") <- list(method = method, graph = graph,
                                laplacian = laplacian, sigma = sigma_global,
                                k = k, radius = radius, alpha = alpha,
                                steps = steps, self_loops = self_loops)
  score
}


# ==============================================================================
# Null check: permute tumor labels within a sample, preserving their count.
# A correctly specified score field should come out flat. Run this once per
# dataset -- it is the cheapest possible sanity check on the operator.
# ==============================================================================
tumor_score_null_check <- function(coords, is_tumor, samples, seed = 1, ...) {
  set.seed(seed)
  perm <- is_tumor
  for (s in unique(as.character(samples))) {
    idx <- which(as.character(samples) == s)
    perm[idx] <- sample(is_tumor[idx])
  }
  real <- compute_tumor_score2(coords, is_tumor, samples, verbose = FALSE, ...)
  null <- compute_tumor_score2(coords, perm,     samples, verbose = FALSE, ...)
  ok <- is.finite(real) & is.finite(null)
  out <- list(
    cor_real_null = suppressWarnings(cor(real[ok], null[ok], method = "spearman")),
    real_iqr = diff(quantile(real[ok], c(0.25, 0.75))),
    null_iqr = diff(quantile(null[ok], c(0.25, 0.75))),
    real_score = real, null_score = null)
  cat(sprintf("Null check: spearman(real, permuted-label) = %.3f\n", out$cor_real_null))
  cat(sprintf("  real field IQR = %.4f   permuted field IQR = %.4f  (ratio %.1fx)\n",
              out$real_iqr, out$null_iqr, out$real_iqr / out$null_iqr))
  cat("  Want: correlation near 0 and a permuted IQR much smaller than the real one.\n")
  invisible(out)
}