# 16_07_inspect_bins.R
# ============================================================================
# Spatial inspection of the tumor-proximity SCORE BINS used by the 15_04 / 16_04
# trajectory analysis. Standalone, read-only diagnostic.
#
# ---------------------------------------------------------------------------
# WHAT IT DOES
# ---------------------------------------------------------------------------
# Recomputes the diffusion local-tumor-fraction score and the quantile bins
# EXACTLY as 15_04 / 16_04 do (same compute_tumor_score() settings, same
# make_bins() quantile cut on the trajectory cell type), then draws a grid of
# spatial maps to show WHERE each bin's cells sit in tissue:
#
#   rows    = score bins (Bin 1 = farthest from tumor ... Bin N = tumor core)
#   columns = samples (WT first, then DKO)
#
# In every panel all cells are grey; the cells belonging to THAT row's bin (for
# that sample) are highlighted in red. So row 1 highlights every Bin-1 cell
# across all samples, row 2 highlights Bin-2 cells, etc. This makes the spatial
# meaning of the trajectory bins directly visible and lets you sanity-check that
# low bins are tumor-distal and high bins are tumor-core.
#
# This does NOT fit any model or modify 15_04 / 16_04 — it only visualizes their
# binning. Bins are defined on the trajectory cell type (default CD8+ T), exactly
# as the pseudobulk GAM uses them.
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

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "16_inspect_bins")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Parameters — synced to 15_04 / 16_04 (diffusion score + fixed-quantile bins)
# ============================================================================
CELLTYPE         <- "CD8+ T"       # trajectory cell type whose bins we inspect
TARGET_CELLTYPE  <- "Tumor"

# Score (local tumor fraction, no per-sample min-max) — matches 15_04 / 16_04
KNN_K            <- 30
DIFFUSION_VARIANT<- "heat"
DIFFUSION_HEAT_STEPS <- 20
DIFFUSION_HEAT_ALPHA <- 2
LOCAL_FRACTION   <- TRUE

# Binning — matches 16_04 latest (fixed N_BINS quantile bins on the trajectory cells)
N_BINS           <- 5
APPLY_QC         <- TRUE           # QC applied to trajectory-cell selection (15_04)

# Plot styling
HIGHLIGHT_COLOR       <- "#D62728" # red for the current bin's cells
BG_COLOR              <- "grey85"  # all other cells
PT_SIZE_BG            <- 0.15      # background points (trajectory-cell-type grid)
PT_SIZE_HL            <- 0.60      # highlighted bin cells (trajectory-cell-type grid)
PT_SIZE_BG_ALL        <- 0.06      # background points (all-cells grid — much denser)
PT_SIZE_HL_ALL        <- 0.10      # highlighted bin cells (all-cells grid)
SHOW_TUMOR_BACKGROUND <- FALSE     # TRUE = tint tumor cells to show the tumor mass
TUMOR_BG_COLOR        <- "#9ECAE1"

# ============================================================================
# Load + score (identical recoding to 15_04 / 16_04)
# ============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

all_coords <- Embeddings(seu, "spatial")
samples    <- seu$sample
celltypes  <- seu$Tier1_celltype
is_tumor   <- celltypes == TARGET_CELLTYPE

qc_pass <- if (APPLY_QC) {
  seu$total_counts >= MIN_TOTAL_COUNTS & seu$n_genes_by_counts >= MIN_GENES_DETECTED
} else rep(TRUE, ncol(seu))
cat(sprintf("QC (APPLY_QC=%s): %d / %d cells pass\n", APPLY_QC, sum(qc_pass), length(qc_pass)))

cat("Computing local-tumor-fraction score...\n")
tumor_score <- compute_tumor_score(
  all_coords, is_tumor, samples, method = "diffusion",
  diffusion_method = DIFFUSION_VARIANT, k = KNN_K, heat_steps = DIFFUSION_HEAT_STEPS,
  alpha = DIFFUSION_HEAT_ALPHA,
  local_fraction = LOCAL_FRACTION, normalize = !LOCAL_FRACTION, verbose = FALSE)

# ============================================================================
# Bins — copied verbatim from 15_04 make_bins() for exact consistency
# ============================================================================
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

# ============================================================================
# Coordinates — center each sample about its own midrange (as in 15_05) so all
# samples share the symmetric plotting window; does not affect the score/bins.
# Background = ALL cells (grey), computed once and reused by every grid.
# ============================================================================
sample_ids <- levels(seu$sample)
cx <- ave(all_coords[, 1], samples, FUN = function(x) x - mean(range(x)))
cy <- ave(all_coords[, 2], samples, FUN = function(x) x - mean(range(x)))
max_x_span <- max(tapply(all_coords[, 1], samples, function(x) diff(range(x))))
max_y_span <- max(tapply(all_coords[, 2], samples, function(x) diff(range(x))))
half_span  <- max(max_x_span, max_y_span) / 2 * 1.05

# Background = ALL cells (no `bin` column -> repeated in every facet row).
bg_df    <- data.frame(x = cx, y = cy,
                       sample = factor(as.character(samples), levels = sample_ids))
tumor_df <- bg_df[is_tumor, , drop = FALSE]   # optional tumor tint

# Rasterize point layers when ggrastr is available (the grey background is repeated
# in every bin row, so a 10-row grid draws ~ n_cells x N_BINS points).
gp <- function(data, color, size, alpha) {
  g <- geom_point(data = data, aes(x, y), color = color, size = size, alpha = alpha)
  if (requireNamespace("ggrastr", quietly = TRUE)) ggrastr::rasterise(g, dpi = 200) else g
}

# ============================================================================
# make_grid(): rows = bins, columns = samples; the given cells' bin members are
# highlighted red on the grey all-cell background. Used for BOTH the trajectory
# cell type (CD8+ T) and ALL cells regardless of type.
# ============================================================================
make_grid <- function(fg_idx, title_lab, out_prefix, pt_hl, pt_bg = PT_SIZE_BG) {
  fscore   <- tumor_score[fg_idx]
  fsample  <- factor(as.character(samples[fg_idx]), levels = sample_ids)
  bininfo  <- make_bins(fscore, N_BINS)
  fbin     <- as.integer(bininfo$bins)
  cat(sprintf("%s: %d cells binned into %d quantile bins\n",
              title_lab, length(fg_idx), N_BINS))

  bin_labels <- sprintf("Bin %d\n[%.3f, %.3f]", seq_len(N_BINS),
                        bininfo$breaks[-(N_BINS + 1)], bininfo$breaks[-1])
  fg_df <- data.frame(
    x = cx[fg_idx], y = cy[fg_idx], sample = fsample,
    bin = factor(fbin, levels = seq_len(N_BINS), labels = bin_labels))

  # cells-per-bin-per-sample table
  cell_tab <- as.data.frame(table(sample = fsample, bin = fbin), stringsAsFactors = FALSE)
  cell_tab$condition <- ifelse(cell_tab$sample %in% WT_SAMPLES,
                               REFERENCE_CONDITION, TEST_CONDITION)
  cell_tab <- cell_tab[order(as.integer(cell_tab$bin), cell_tab$sample),
                       c("sample", "condition", "bin", "Freq")]
  names(cell_tab)[names(cell_tab) == "Freq"] <- "n_cells"
  write.csv(cell_tab, file.path(OUTPUT_DIR, sprintf("cells_per_bin_per_sample_%s.csv", out_prefix)),
            row.names = FALSE)

  p <- ggplot()
  if (SHOW_TUMOR_BACKGROUND && nrow(tumor_df) > 0)
    p <- p + gp(tumor_df, TUMOR_BG_COLOR, pt_bg, 0.6)
  p <- p +
    gp(bg_df, BG_COLOR, pt_bg, 0.5) +
    gp(fg_df, HIGHLIGHT_COLOR, pt_hl, 0.9) +
    facet_grid(bin ~ sample, switch = "y") +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    labs(title = sprintf("%s tumor-proximity bins in space (red = cells in that bin)", title_lab),
         subtitle = sprintf("diffusion local-tumor-fraction score, %d quantile bins (Bin 1 = tumor-distal, Bin %d = tumor-core)%s",
                            N_BINS, N_BINS,
                            if (SHOW_TUMOR_BACKGROUND) sprintf("; blue = %s cells", TARGET_CELLTYPE) else "")) +
    theme_void(base_size = 11) +
    theme(strip.text.x = element_text(size = 9, face = "bold"),
          strip.text.y.left = element_text(size = 7, angle = 0),
          panel.spacing = unit(0.15, "lines"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(8, 8, 8, 8))

  w <- 2.4 * length(sample_ids) + 1.8
  h <- 2.0 * N_BINS + 1.2
  ggsave(file.path(OUTPUT_DIR, sprintf("bin_inspection_spatial_%s.pdf", out_prefix)),
         p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(OUTPUT_DIR, sprintf("bin_inspection_spatial_%s.png", out_prefix)),
         p, width = w, height = h, dpi = 200, limitsize = FALSE)
}

# ============================================================================
# (1) Trajectory cell type (default CD8+ T) — same selection as 16_04
# ============================================================================
ct_idx <- which(celltypes == CELLTYPE & !is.na(tumor_score) & qc_pass)
if (length(ct_idx) < 50) stop(sprintf("Too few %s cells (%d)", CELLTYPE, length(ct_idx)))
ct_prefix <- gsub("[^A-Za-z0-9]+", "_", gsub("\\+", "pos", CELLTYPE))
make_grid(ct_idx, CELLTYPE, ct_prefix, pt_hl = PT_SIZE_HL)

# ============================================================================
# (2) ALL cells regardless of cell type — whole-tissue score deciles. Bins are
# recomputed over ALL scored cells, so each bin holds ~1/N_BINS of the tissue;
# points are small so the dense red does not swamp the panel.
# ============================================================================
all_idx <- which(!is.na(tumor_score))
make_grid(all_idx, "All cells", "all_cells", pt_hl = PT_SIZE_HL_ALL, pt_bg = PT_SIZE_BG_ALL)

cat(sprintf("Done. Grids (%d bins x %d samples) -> %s\n",
            N_BINS, length(sample_ids), OUTPUT_DIR))
