# ==============================================================================
# 00_utils.R — Reusable helper functions
# ==============================================================================
# Source this file at the top of an analysis script:
#   source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))
#
# Sections:
#   1. Tumor-proximity scoring (trajectory scripts 14-16): compute_tumor_score(),
#      compute_tumor_score2() and their low-level scorers.
#   2. Pathway analysis (12_pathway_activity.R): orthology mapping, gene-set
#      similarity, EnrichmentMap-style pathway modules and their plots. See the
#      "PATHWAY ANALYSIS HELPERS" header at the end of this file.
#
# --- Section 1: tumor-proximity scoring ---
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



# ==============================================================================
# PATHWAY ANALYSIS HELPERS (used by 12_pathway_activity.R)
# ==============================================================================
# Functions for the differential pathway-activity analysis of domain pseudobulk
# profiles (PanomiR, DKO vs WT) and for grouping the tested MSigDB pathways into
# EnrichmentMap-style "pathway modules": sets of pathways that share most of
# their genes, so that results read as a few biological themes instead of dozens
# of overlapping pathway names. They are written for pathway analysis and assume
# its data structures:
#
#   pathways        MSigDB gene-membership table (columns Pathway, ENSEMBL), as
#                   read from MSIGDB_PATHWAYS
#   M               sparse 0/1 pathway x gene matrix (rows = MSigDB standard
#                   names, columns = human ENSEMBL IDs); build_pathway_gene_matrix()
#   inter, sim      pathway x pathway intersection-size / similarity matrices
#   module          named integer vector, pathway -> module ID (M1 = largest)
#   pathway tables  PanomiR/limma results: standard_name, display_name, logFC,
#                   P.Value, adj.P.Val (+ module columns once modules are built)
#
# Contents:
#   Orthology         convert_mouse_to_human_ensembl, human_ensembl_to_mouse
#   Small utilities   domain_from_pseudobulk_file, pathway_label, significant_pathways,
#                     truncate_label, wrap_word_lines, save_pdf_png
#   Gene-set overlap  build_pathway_gene_matrix, pathway_intersections, pathway_unions,
#                     pathway_min_sizes, pathway_similarity, pathway_similarity_label,
#                     pathway_module_rule
#   Modules           cluster_pathway_modules, module_medoids, module_threshold_sweep,
#                     module_membership_matrix, module_gene_frequency,
#                     module_pathway_table, module_gene_table, module_summary_table,
#                     format_core_genes
#   Plots             plot_similarity_distribution, plot_threshold_sensitivity,
#                     plot_module_sizes, plot_similarity_heatmap_all,
#                     plot_similarity_heatmap_significant, plot_module_network,
#                     plot_pathway_barplot, plot_pathway_gene_heatmap
#
# Requirements: the sourcing script attaches ggplot2 (the plot helpers call it
# unqualified) and sources 00_parameters.R (some defaults read FDR_THRESHOLD,
# PATHWAY_NAME_TRUNC, MODULE_CORE_GENE_FRAC and MODULE_CORE_GENE_MIN_N).
# Matrix, pheatmap, igraph, ggraph, ggrepel, ggtext, cluster, dplyr and
# babelgene only need to be installed; they are called with ::.
# ==============================================================================

# Colours shared by the pathway plots: DE direction of a pathway (DKO vs WT) and
# a single-hue ramp for pathway-similarity heatmaps
.pathway_direction_colors <- c(Up = "#C62828", Down = "#1565C0", NS = "grey92")
.pathway_similarity_ramp  <- grDevices::colorRampPalette(c("white", "#263238"))(100)


# ------------------------------------------------------------------------------
# Pathway analysis: mouse <-> human orthology
# ------------------------------------------------------------------------------
# MSigDB pathways are defined on human ENSEMBL IDs while the panel is mouse, so
# genes are mapped by true ortholog (babelgene / HCOP consensus of 12 databases),
# NOT by naive symbol case conversion. Uppercasing fails for the many genes whose
# human ortholog has a different symbol (e.g. Trp53 -> TP53, H2-K1 -> HLA-A) and
# can spuriously collide with unrelated human genes.

# Mouse symbols -> human ENSEMBL, resolved to a 1:1 mapping (named vector,
# names = mouse symbols) so pseudobulk rows can be renamed for PanomiR.
convert_mouse_to_human_ensembl <- function(mouse_symbols) {
  ortho <- babelgene::orthologs(
    genes   = unique(mouse_symbols),
    species = "mouse",
    human   = FALSE  # input genes are the non-human (mouse) species
  )

  # Keep only rows with a human ENSEMBL ID
  ortho <- ortho[!is.na(ortho$human_ensembl) & ortho$human_ensembl != "", ]

  # Resolve to a 1:1 mouse-symbol <-> human-ENSEMBL mapping, preferring the
  # ortholog with the strongest cross-database support:
  ortho <- ortho[order(-ortho$support_n), ]
  #  (a) one human ENSEMBL per mouse symbol
  ortho <- ortho[!duplicated(ortho$symbol), ]
  #  (b) one mouse symbol per human ENSEMBL (avoids duplicate row names downstream)
  ortho <- ortho[!duplicated(ortho$human_ensembl), ]

  setNames(ortho$human_ensembl, ortho$symbol)
}

# Human ENSEMBL -> mouse symbols (reverse direction), used to show a pathway's
# genes in the mouse data. A human gene may have several mouse orthologs (e.g.
# HLA-A -> H2-K1/H2-Q7/...), so all are kept (columns ENSEMBL, mouse_symbol) and
# later intersected with the genes actually present in the data.
human_ensembl_to_mouse <- function(ensembl_ids) {
  ortho <- babelgene::orthologs(
    genes   = unique(ensembl_ids),
    species = "mouse",
    human   = TRUE  # input genes are human (ENSEMBL IDs)
  )
  ortho <- ortho[!is.na(ortho$symbol) & ortho$symbol != "", ]
  ortho$ENSEMBL <- ortho$human_ensembl
  ortho$mouse_symbol <- ortho$symbol
  ortho
}


# ------------------------------------------------------------------------------
# Pathway analysis: small utilities
# ------------------------------------------------------------------------------

# Domain name from a script-04 pseudobulk file used as pathway-analysis input,
# e.g. ".../pseudobulk_domain_3_Tumor Core.csv" -> "Tumor Core"
domain_from_pseudobulk_file <- function(path) {
  name <- gsub("pseudobulk_domain_", "", basename(path))
  name <- gsub("\\.csv$", "", name)
  gsub("^[0-9]_", "", name)
}

# Readable pathway name: the MSigDB display name, or the standard name if missing
pathway_label <- function(standard_name, display_name) {
  ifelse(is.na(display_name) | display_name == "", standard_name, display_name)
}

# Rows of a pathway table that are significant: FDR-adjusted p-value
# (adj.P.Val, Benjamini-Hochberg from limma/PanomiR) below fdr
significant_pathways <- function(df, fdr = FDR_THRESHOLD) {
  df[!is.na(df$adj.P.Val) & df$adj.P.Val < fdr, ]
}

# Cut long pathway / module names to n characters, marking the cut with "*"
truncate_label <- function(x, n) {
  ifelse(nchar(x) > n, paste0(substr(x, 1, n), "*"), x)
}

# Greedy word wrap for pathway-plot labels: the line number of each word when
# lines are <= width characters
wrap_word_lines <- function(words, width) {
  line <- integer(length(words))
  n <- 1L
  len <- 0L
  for (i in seq_along(words)) {
    w <- nchar(words[i])
    if (len > 0 && len + 1 + w > width) {
      n <- n + 1L
      len <- 0L
    }
    len <- len + (len > 0) + w
    line[i] <- n
  }
  line
}

# Save a pathway-analysis ggplot as <stem>.pdf and <stem>.png
save_pdf_png <- function(p, stem, width, height) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot = p, device = cairo_pdf, width = width, height = height,
                  limitsize = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), plot = p, width = width, height = height, dpi = 300,
                  limitsize = FALSE)
}


# ------------------------------------------------------------------------------
# Pathway analysis: matrix-based gene-set similarity between pathways
# ------------------------------------------------------------------------------
# All pairwise overlaps come from one sparse matrix product on the pathway x gene
# matrix M, rather than looping over pathway pairs.

# Binary pathway x gene incidence matrix (sparse) over the valid genes (those
# PanomiR scored). Keeps pathways with >= min_size valid genes, mirroring
# PanomiR's filter, so the rows are exactly the tested pathways.
build_pathway_gene_matrix <- function(pathways, genes, min_size) {
  pw <- unique(pathways[pathways$ENSEMBL %in% genes, c("Pathway", "ENSEMBL")])
  path_sizes <- table(pw$Pathway)
  keep <- sort(names(path_sizes)[path_sizes >= min_size])
  pw <- pw[pw$Pathway %in% keep, ]
  genes <- sort(unique(genes))
  Matrix::sparseMatrix(
    i = match(pw$Pathway, keep),
    j = match(pw$ENSEMBL, genes),
    x = 1,
    dims = c(length(keep), length(genes)),
    dimnames = list(keep, genes)
  )
}

# |A ∩ B| for every pathway pair: M %*% t(M). The diagonal holds the set sizes |A|.
pathway_intersections <- function(M) {
  as.matrix(Matrix::tcrossprod(M))
}

# |A ∪ B| = |A| + |B| - |A ∩ B|
pathway_unions <- function(inter) {
  set_sizes <- diag(inter)
  outer(set_sizes, set_sizes, "+") - inter
}

# min(|A|, |B|), the denominator of the overlap coefficient
pathway_min_sizes <- function(inter) {
  set_sizes <- diag(inter)
  outer(set_sizes, set_sizes, pmin)
}

# Pairwise pathway similarity from the intersection matrix:
#   jaccard  = |A ∩ B| / |A ∪ B|
#   overlap  = |A ∩ B| / min(|A|, |B|)
#   combined = k * jaccard + (1 - k) * overlap  (EnrichmentMap's combined coefficient)
pathway_similarity <- function(inter, method = c("jaccard", "overlap", "combined"), k = 0.5) {
  method <- match.arg(method)
  jaccard <- inter / pathway_unions(inter)
  overlap <- inter / pathway_min_sizes(inter)
  switch(method,
    jaccard  = jaccard,
    overlap  = overlap,
    combined = k * jaccard + (1 - k) * overlap
  )
}

# Plot / caption names for the similarity measure and the module rule
pathway_similarity_label <- function(method) {
  c(jaccard = "Jaccard index", overlap = "Overlap coefficient",
    combined = "Combined Jaccard/overlap")[[method]]
}

pathway_module_rule <- function(method, threshold, linkage) {
  sprintf("%s ≥ %.2f, %s linkage", pathway_similarity_label(method), threshold, linkage)
}


# ------------------------------------------------------------------------------
# Pathway analysis: pathway modules
# ------------------------------------------------------------------------------
# Modules come from hierarchical clustering on 1 - similarity, cut at
# 1 - threshold. With average linkage, the pathways in a module have an average
# pairwise similarity >= threshold. (Plain connected components of the
# thresholded pathway network chain nearly all pathways into one giant module on
# a targeted gene panel, so they are not used.)

# Cut the hierarchical tree at 1 - threshold. Modules are renumbered by size
# (M1 = largest; ties broken by the alphabetically first member).
cluster_pathway_modules <- function(sim, threshold, linkage) {
  hc <- hclust(as.dist(1 - sim), method = linkage)
  cl <- cutree(hc, h = 1 - threshold)
  cl_sizes <- table(cl)
  first_member <- tapply(names(cl), cl, min)[names(cl_sizes)]
  new_id <- setNames(seq_along(cl_sizes), names(cl_sizes)[order(-as.integer(cl_sizes), first_member)])
  list(hclust = hc, module = setNames(as.integer(new_id[as.character(cl)]), names(cl)))
}

# Medoid of each module (used to name it): the member with the highest mean
# similarity to the other members
module_medoids <- function(sim, module) {
  members <- split(names(module), module)
  vapply(members, function(m) m[which.max(rowMeans(sim[m, m, drop = FALSE]))], character(1))
}

# Module statistics when the same tree is cut at each threshold in `thresholds`
# (the threshold-sensitivity diagnostic)
module_threshold_sweep <- function(sim, hc, thresholds) {
  d <- as.dist(1 - sim)
  do.call(rbind, lapply(thresholds, function(thr) {
    cl <- cutree(hc, h = 1 - thr)
    cl_sizes <- table(cl)
    sil <- if (length(cl_sizes) > 1 && length(cl_sizes) < length(cl)) {
      mean(cluster::silhouette(cl, d)[, "sil_width"])
    } else {
      NA_real_
    }
    data.frame(threshold = thr, n_modules = length(cl_sizes), n_singletons = sum(cl_sizes == 1),
               largest_module = max(cl_sizes), mean_silhouette = sil)
  }))
}

# Module x pathway membership matrix (sparse): 1 if the pathway belongs to the module
module_membership_matrix <- function(module) {
  ids <- sort(unique(module))
  Matrix::sparseMatrix(
    i = match(module, ids),
    j = seq_along(module),
    x = 1,
    dims = c(length(ids), length(module)),
    dimnames = list(ids, names(module))
  )
}

# Fraction of each module's pathways that contain each gene: (A %*% M) / module
# sizes, where A is the module x pathway membership matrix and M the pathway x
# gene matrix. Passing a row subset of M (e.g. significant pathways only) gives
# the fraction among those pathways.
module_gene_frequency <- function(M, module) {
  A <- module_membership_matrix(module[rownames(M)])
  as.matrix(A %*% M) / Matrix::rowSums(A)
}

# Pathway table for modules: every tested pathway with its gene count, module,
# module size, medoid, readable label, module label ("M<id>: <medoid name>") and
# its DE statistics from `dep`; ordered by module, then P.Value.
module_pathway_table <- function(M, inter, module, medoids, dep) {
  df <- data.frame(
    standard_name = rownames(M),
    n_genes       = diag(inter),
    module        = module[rownames(M)],
    row.names     = NULL
  )
  df$module_size <- as.integer(table(module)[as.character(df$module)])
  df$medoid <- medoids[as.character(df$module)]
  df$is_medoid <- df$standard_name == df$medoid
  df <- dplyr::left_join(
    df,
    dep[, c("standard_name", "display_name", "logFC", "P.Value", "adj.P.Val")],
    by = "standard_name"
  )
  df$label <- pathway_label(df$standard_name, df$display_name)
  df$module_label <- paste0("M", df$module, ": ", df$label[match(df$medoid, df$standard_name)])
  df[order(df$module, df$P.Value), ]
}

# Module gene table: how often each gene occurs across all of a module's pathways
# (frac) and across its significant pathways (frac_sig; sig_names = significant
# pathways in this domain). Core genes are in more than core_frac of the module's
# significant pathways and in at least core_min_n of them (modules with >= 2
# significant pathways only). gene_symbols maps human ENSEMBL -> mouse symbol.
module_gene_table <- function(M, module, sig_names, gene_symbols, core_frac, core_min_n) {
  freq <- module_gene_frequency(M, module)
  nz <- which(freq > 0, arr.ind = TRUE)
  genes <- data.frame(
    module  = as.integer(rownames(freq)[nz[, "row"]]),
    ensembl = colnames(freq)[nz[, "col"]],
    frac    = freq[nz]
  )
  genes$gene <- gene_symbols[genes$ensembl]
  genes$module_size <- as.integer(table(module)[as.character(genes$module)])
  genes$n_pathways <- round(genes$frac * genes$module_size)

  n_sig_module <- as.integer(table(module[sig_names])[as.character(genes$module)])
  genes$n_sig_module <- ifelse(is.na(n_sig_module), 0L, n_sig_module)
  genes$frac_sig <- NA_real_
  if (length(sig_names) > 0) {
    # Same matrix product, restricted to the significant pathways (rows of M).
    # Modules without significant pathways have no row, so they stay NA.
    freq_sig <- module_gene_frequency(M[sig_names, , drop = FALSE], module)
    genes$frac_sig <- freq_sig[cbind(
      match(as.character(genes$module), rownames(freq_sig)),
      match(genes$ensembl, colnames(freq_sig))
    )]
  }
  genes$n_sig_pathways <- round(genes$frac_sig * genes$n_sig_module)
  genes$is_core <- genes$n_sig_module >= 2 & !is.na(genes$frac_sig) &
    genes$frac_sig > core_frac & genes$n_sig_pathways >= core_min_n
  genes[order(genes$module, -genes$frac_sig, -genes$frac, genes$gene), ]
}

# Per-module summary: size, significant / up / down counts (adj.P.Val < fdr),
# best adjusted p-value, mean logFC, cohesion (mean pairwise similarity within
# the module) and its core genes (from module_gene_table, most frequent first).
module_summary_table <- function(modules_df, sim, module, module_genes, fdr) {
  cohesion <- vapply(split(names(module), module), function(m) {
    if (length(m) < 2) return(NA_real_)
    s <- sim[m, m]
    mean(s[upper.tri(s)])
  }, numeric(1))
  out <- as.data.frame(dplyr::summarise(
    dplyr::group_by(modules_df, module, module_label, module_size),
    n_significant = sum(adj.P.Val < fdr, na.rm = TRUE),
    n_up          = sum(adj.P.Val < fdr & logFC > 0, na.rm = TRUE),
    n_down        = sum(adj.P.Val < fdr & logFC < 0, na.rm = TRUE),
    min_adj_P     = min(adj.P.Val, na.rm = TRUE),
    mean_logFC    = mean(logFC, na.rm = TRUE),
    .groups = "drop"
  ))
  out$mean_within_similarity <- cohesion[as.character(out$module)]

  core <- module_genes[module_genes$is_core, ]
  core_genes <- tapply(core$gene, core$module, paste, collapse = ";")
  out$n_core_genes <- as.integer(table(factor(core$module, levels = out$module)))
  out$core_genes <- unname(core_genes[as.character(out$module)])
  out
}

# Label text for a module's core genes on pathway bar plots: all of them,
# comma-separated and wrapped onto lines of <= width characters. With gene_lfc
# (named gene-level log2FC vector) it returns ggtext markdown in which each gene
# is red (log2FC > 0) or blue (log2FC < 0); genes without a log2FC keep the
# default colour. Both versions break lines at the same places.
format_core_genes <- function(core_genes, gene_lfc = NULL, width = 70) {
  # as.character: read.csv makes an all-NA column logical (domains with no core genes)
  vapply(strsplit(as.character(core_genes), ";"), function(g) {
    if (all(is.na(g))) return("no gene in most pathways")
    commas <- c(rep(",", length(g) - 1), "")
    words <- paste0(g, commas)
    line <- wrap_word_lines(words, width)
    sep <- "\n"
    if (!is.null(gene_lfc)) {
      lfc <- gene_lfc[g]
      col <- ifelse(!is.na(lfc) & lfc > 0, "#C62828", ifelse(!is.na(lfc) & lfc < 0, "#1565C0", NA))
      g <- ifelse(is.na(col), g, sprintf("<span style='color:%s'>%s</span>", col, g))
      words <- paste0(g, commas)
      sep <- "<br>"
    }
    paste(vapply(split(words, line), paste, character(1), collapse = " "), collapse = sep)
  }, character(1))
}


# ------------------------------------------------------------------------------
# Pathway analysis: module diagnostic plots
# ------------------------------------------------------------------------------
# One call per diagnostic; `stem` / `file` is the output path (without extension
# for ggplot figures saved as PDF + PNG, with ".png" for pheatmap figures).

# Distribution of pairwise pathway similarity (pairs sharing >= 1 gene), with
# the module threshold marked
plot_similarity_distribution <- function(sim, threshold, sim_label, domain_name, stem) {
  pair_sim <- sim[upper.tri(sim)]
  p <- ggplot(data.frame(similarity = pair_sim[pair_sim > 0]), aes(x = similarity)) +
    geom_histogram(binwidth = 0.02, boundary = 0, fill = "grey55", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = threshold, linetype = "dashed", colour = "#C62828") +
    annotate("text", x = threshold, y = Inf, label = " threshold",
             hjust = 0, vjust = 1.5, size = 3.5, colour = "#C62828") +
    scale_y_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000)) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = paste0("Pairwise pathway similarity in ", domain_name),
      subtitle = sprintf("%d pathways; %.0f%% of pairs share no genes (not shown); %.1f%% of pairs ≥ %.2f",
                         nrow(sim), 100 * mean(pair_sim == 0), 100 * mean(pair_sim >= threshold), threshold),
      x = sim_label, y = "Pathway pairs (log scale)"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 14, face = "bold"))
  save_pdf_png(p, stem, 8, 5)
}

# Module count, singletons, largest module and mean silhouette width as the
# threshold varies (sweep from module_threshold_sweep), chosen threshold marked
plot_threshold_sensitivity <- function(sweep, threshold, sim_label, module_rule, domain_name, stem) {
  metric_names <- c(n_modules = "Number of modules", n_singletons = "Singleton modules",
                    largest_module = "Largest module (pathways)", mean_silhouette = "Mean silhouette width")
  sweep_long <- do.call(rbind, lapply(names(metric_names), function(m) {
    data.frame(threshold = sweep$threshold, metric = metric_names[[m]], value = sweep[[m]])
  }))
  sweep_long$metric <- factor(sweep_long$metric, levels = metric_names)
  p <- ggplot(sweep_long, aes(x = threshold, y = value)) +
    geom_vline(xintercept = threshold, linetype = "dashed", colour = "#C62828") +
    geom_line(colour = "grey30", linewidth = 0.6) +
    geom_point(colour = "grey30", size = 1.5) +
    facet_wrap(~metric, ncol = 2, scales = "free_y") +
    labs(
      title = paste0("Module threshold sensitivity in ", domain_name),
      subtitle = paste0("Same tree re-cut at each threshold; dashed line = chosen threshold (", module_rule, ")"),
      x = paste0(sim_label, " threshold"), y = NULL
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), strip.background = element_blank(),
          strip.text = element_text(face = "bold"), plot.title = element_text(size = 14, face = "bold"))
  save_pdf_png(p, stem, 9, 6)
}

# Number of modules of each size
plot_module_sizes <- function(module, module_rule, domain_name, stem) {
  size_counts <- as.data.frame(table(module_size = as.integer(table(module))))
  p <- ggplot(size_counts, aes(x = module_size, y = Freq)) +
    geom_col(fill = "grey45", width = 0.8) +
    labs(
      title = paste0("Module sizes in ", domain_name),
      subtitle = sprintf("%d pathways in %d modules (%d singletons); %s",
                         length(module), length(unique(module)), sum(table(module) == 1), module_rule),
      x = "Pathways per module", y = "Number of modules"
    ) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size = 14, face = "bold"))
  save_pdf_png(p, stem, 8, 5)
}

# Similarity heatmap of all tested pathways, ordered by the module tree (hc),
# annotated with each pathway's DE direction (named vector: Up / Down / NS)
plot_similarity_heatmap_all <- function(sim, hc, direction, sim_label, domain_name, file) {
  pheatmap::pheatmap(
    sim,
    cluster_rows = hc,
    cluster_cols = hc,
    treeheight_row = 40,
    treeheight_col = 0,
    color = .pathway_similarity_ramp,
    breaks = seq(0, 1, length.out = 101),
    border_color = NA,
    annotation_row = data.frame(Direction = direction[rownames(sim)], row.names = rownames(sim)),
    annotation_colors = list(Direction = .pathway_direction_colors),
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = sprintf("%s: all %d tested pathways (%s)", domain_name, nrow(sim), sim_label),
    filename = file,
    width = 8,
    height = 7
  )
}

# Similarity heatmap of the significant pathways (sig_df: rows of the module
# pathway table), split by module. Modules are contiguous in the tree's leaf
# order, so ordering by it groups them.
plot_similarity_heatmap_significant <- function(sim, hc, module, sig_df, direction, sim_label, domain_name, file) {
  leaf_order <- hc$labels[hc$order]
  sig_ord <- sig_df$standard_name[order(match(sig_df$standard_name, leaf_order))]
  sig_mod <- module[sig_ord]
  sig_lab <- sig_df$label[match(sig_ord, sig_df$standard_name)]
  pheatmap::pheatmap(
    sim[sig_ord, sig_ord],
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = which(diff(sig_mod) != 0),
    gaps_col = which(diff(sig_mod) != 0),
    color = .pathway_similarity_ramp,
    breaks = seq(0, 1, length.out = 101),
    border_color = NA,
    annotation_row = data.frame(Direction = direction[sig_ord], row.names = sig_ord),
    annotation_colors = list(Direction = .pathway_direction_colors[c("Up", "Down")]),
    labels_row = paste0("M", sig_mod, "  ", truncate_label(sig_lab, 60)),
    show_colnames = FALSE,
    fontsize_row = 6,
    main = sprintf("%s: %d significant pathways, split by module (%s)",
                   domain_name, length(sig_ord), sim_label),
    filename = file,
    width = 11,
    height = max(5, length(sig_ord) * 0.11 + 2)
  )
}

# EnrichmentMap-style network of the significant pathways (sig_df). Edges join
# pathway pairs with similarity >= threshold. Modules with >= 2 significant
# pathways get a shaded hull and a label at their centroid.
plot_module_network <- function(sim, sig_df, threshold, fdr, sim_label, domain_name, stem) {
  sub_sim <- sim[sig_df$standard_name, sig_df$standard_name]
  g <- igraph::graph_from_adjacency_matrix(
    sub_sim * (sub_sim >= threshold),
    mode = "undirected", weighted = TRUE, diag = FALSE
  )
  igraph::V(g)$logFC <- sig_df$logFC
  igraph::V(g)$n_genes <- sig_df$n_genes
  igraph::V(g)$module <- sig_df$module
  lay <- ggraph::create_layout(g, layout = "stress")

  multi <- sig_df$module %in% sig_df$module[duplicated(sig_df$module)]
  multi_nodes <- as.data.frame(lay)[multi, c("x", "y", "module")]
  # Convex hull per module; a thick rounded outline pads it (and turns 2-node hulls into capsules)
  hulls <- do.call(rbind, lapply(split(multi_nodes, multi_nodes$module), function(d) d[grDevices::chull(d$x, d$y), ]))
  centroids <- aggregate(cbind(x, y) ~ module, data = multi_nodes, FUN = mean)
  centroids$text <- truncate_label(sig_df$module_label[match(centroids$module, sig_df$module)], 45)

  p <- ggraph::ggraph(lay)
  if (nrow(multi_nodes) > 0) {
    p <- p + geom_polygon(data = hulls, aes(x = x, y = y, group = module), inherit.aes = FALSE,
                          fill = "grey92", colour = "grey92", linewidth = 6, linejoin = "round")
  }
  if (igraph::ecount(g) > 0) {
    p <- p +
      ggraph::geom_edge_link(aes(edge_width = weight), colour = "grey65", alpha = 0.6) +
      ggraph::scale_edge_width(range = c(0.2, 1.5), name = sim_label)
  }
  p <- p +
    ggraph::geom_node_point(aes(fill = logFC, size = n_genes), shape = 21, colour = "grey30", stroke = 0.3) +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828", midpoint = 0,
                         name = "Effect\nsize") +
    scale_size_continuous(range = c(1.5, 6), name = "Genes")
  if (nrow(centroids) > 0) {
    p <- p + ggrepel::geom_label_repel(
      data = centroids, aes(x = x, y = y, label = text), inherit.aes = FALSE,
      size = 2.8, fill = "white", label.size = 0.2, box.padding = 0.6, min.segment.length = 0,
      segment.colour = "grey40", segment.size = 0.3, max.overlaps = Inf, seed = 1
    )
  }
  p <- p +
    labs(
      title = paste0("Pathway module network in ", domain_name),
      subtitle = sprintf("%d significant pathways (adj.P < %.2f); edges: %s ≥ %.2f",
                         nrow(sig_df), fdr, sim_label, threshold)
    ) +
    theme_void() +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          plot.margin = margin(10, 10, 10, 10))
  save_pdf_png(p, stem, 10, 8)
}


# ------------------------------------------------------------------------------
# Pathway analysis: result plots
# ------------------------------------------------------------------------------

# Signed -log10(P) bar plot of DE pathways, saved as <out_stem>.pdf/.png.
# Bars are grouped into blocks by pathway module (columns module, module_size,
# module_label from the module pathway table). Block labels sit on the left and
# pathway names on the right. Each block is labelled with its module medoid, or
# with show_core_genes = TRUE with the module's core genes only (column
# core_genes_text); singleton modules are unlabelled. Adding color_genes = TRUE
# colours the genes by gene-level DEG direction (column core_genes_md, ggtext
# markdown). module_rule is the caption text for how modules were built
# (pathway_module_rule()).
plot_pathway_barplot <- function(paths, title, out_stem, module_rule,
                                 show_core_genes = FALSE, color_genes = FALSE,
                                 name_trunc = PATHWAY_NAME_TRUNC,
                                 core_frac = MODULE_CORE_GENE_FRAC,
                                 core_min_n = MODULE_CORE_GENE_MIN_N,
                                 min_module_members = 1) {

  # Readable, truncated pathway names; signed -log10 p-value: sign(logFC) * -log10(P.Value)
  paths$label <- truncate_label(pathway_label(paths$standard_name, paths$display_name), name_trunc)
  paths$signed_neg_log_p <- sign(paths$logFC) * (-log10(paths$P.Value))

  # Order module blocks by the mean signed -log10 p-value of their plotted
  # pathways (top = most up), and pathways within a block by their own value.
  # Facet rows run top to bottom, while y levels within a block run bottom to top.
  module_score <- tapply(paths$signed_neg_log_p, paths$module, mean)
  module_order <- names(module_score)[order(-module_score, as.integer(names(module_score)))]
  paths$module_f <- factor(paths$module, levels = module_order)
  paths <- paths[order(paths$signed_neg_log_p), ]
  paths$label <- factor(paths$label, levels = paths$label)

  # Block labels (blank for singletons): the module's core genes, or else
  # "M<id>: <medoid> (n=<pathways in module>)"
  multi <- paths$module_size > min_module_members
  if (show_core_genes) {
    strip_text <- ifelse(multi, paths$core_genes_text, "")
    caption <- sprintf(paste0("Blocks = pathway modules (%s)\n",
                              "Labels = genes in >%.0f%% (and >= %d) of the module's significant pathways"),
                       module_rule, 100 * core_frac, core_min_n)
  } else {
    strip_text <- ifelse(multi, paste0(truncate_label(paths$module_label, 60), " (n=", paths$module_size, ")"), "")
    caption <- paste0("Blocks = pathway modules (", module_rule, "); label = module medoid, n = pathways in module")
  }
  # Line counts come from the plain text; the coloured (markdown) version has the same line breaks
  n_lines <- lengths(strsplit(strip_text, "\n"))
  if (color_genes) {
    strip_text <- ifelse(multi, paths$core_genes_md, "")
    caption <- paste0(caption, "\nGene colour = DEG log2FC in this domain (red > 0, blue < 0, black = not tested)")
  }
  first_in_module <- !duplicated(paths$module)
  strip_labels <- setNames(strip_text, paths$module)[first_in_module]
  n_lines <- n_lines[first_in_module]

  # Multi-line labels use a smaller font, and bars get tall enough that every
  # block fits its label (a label line is ~1.2 x font size)
  strip_size <- if (show_core_genes) 8 else 9
  bar_height <- 0.22
  if (show_core_genes) {
    n_bars <- as.integer(table(paths$module)[names(strip_labels)])
    bar_height <- max(0.26, n_lines * strip_size * 1.2 / 72 / n_bars)
  }
  # Left-hand block labels are right-aligned against the bars.
  # gridtext (ggtext) spaces lines more tightly than element_text, hence the larger lineheight
  strip_element <- if (color_genes) {
    ggtext::element_markdown(size = strip_size, angle = 0, hjust = 1, halign = 1, lineheight = 1.2,
                             margin = margin(0, 4, 0, 0))
  } else {
    element_text(size = strip_size, angle = 0, hjust = 1, lineheight = 0.95, margin = margin(0, 4, 0, 0))
  }

  # switch = "y" puts the block labels on the left; pathway names go on the right axis
  p <- ggplot(paths, aes(x = signed_neg_log_p, y = label, fill = logFC)) +
    geom_col(width = 0.8, orientation = "y") +
    facet_grid(rows = vars(module_f), scales = "free_y", space = "free_y", switch = "y",
               labeller = as_labeller(strip_labels)) +
    scale_y_discrete(position = "right") +
    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828", midpoint = 0,
                         name = "Effect\nsize") +
    labs(
      title = title,
      x = "Signed -log10(P-value)",
      y = NULL,
      caption = caption
    ) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
      panel.spacing.y = unit(2, "pt"),
      strip.background = element_blank(),
      strip.text.y.left = strip_element,
      strip.clip = "off",
      axis.title.x = element_text(size = 14),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 1),
      axis.text.y.right = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 0),
      legend.text = element_text(size = 11),
      legend.position = "right",
      legend.title = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 9, colour = "grey30", hjust = 0),
      # Title and caption span the whole figure (the bar panel sits between two label columns)
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )

  plot_height <- max(4, nrow(paths) * bar_height + length(module_order) * 0.04 + 0.3)
  save_pdf_png(p, out_stem, width = 12, height = plot_height)
}

# Heatmap of one pathway's genes across samples (log2 CPM, row-scaled), saved as
# <out_dir>/<pathway standard name>.png. cpm has mouse symbols as rows;
# ensembl_to_mouse comes from human_ensembl_to_mouse(). Pathways with fewer than
# two genes in the data are skipped.
plot_pathway_gene_heatmap <- function(pathway_name, title, pathways, ensembl_to_mouse, cpm, out_dir,
                                      condition_colors) {
  pathway_ensembl <- pathways$ENSEMBL[pathways$Pathway == pathway_name]
  genes <- ensembl_to_mouse$mouse_symbol[ensembl_to_mouse$ENSEMBL %in% pathway_ensembl]
  genes <- genes[genes %in% rownames(cpm)]
  if (length(genes) < 2) return(invisible(NULL))

  mat_log <- log2(as.matrix(cpm[genes, , drop = FALSE]) + 1)
  pheatmap::pheatmap(
    mat_log,
    main = title,
    fontsize = 6,
    scale = "row",
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    treeheight_row = 0,
    treeheight_col = 0,
    annotation_col = data.frame(Condition = gsub("[0-9]*", "", colnames(mat_log)), row.names = colnames(mat_log)),
    annotation_colors = list(Condition = condition_colors),
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = max(4, min(8, 120 / length(genes))),
    fontsize_col = 6,
    filename = file.path(out_dir, paste0(substr(gsub("[^A-Za-z0-9_]", "_", pathway_name), 1, 80), ".png")),
    width = 4,
    height = max(4, length(genes) * 0.2 + 2)
  )
}