# 19_01_tumor_core_domain.R
# Score-defined tumor core domain + WT vs DKO differential expression within it.
#
# The tumor core is defined from a tumor-proximity score computed with the shared
# compute_tumor_score() used in 15_04_tumor_trajectory_working.R. SCORE_METHOD
# switches between:
#   "diffusion" - heat-kernel graph diffusion of the binary tumor indicator,
#                 returned as a kernel-weighted LOCAL TUMOR FRACTION (K x)/(K 1);
#   "knn"       - k-NN purity: fraction of the KNN_K nearest spatial neighbors
#                 that are tumor.
# Both are local tumor fractions natively in [0,1] and comparable across samples,
# so neither is min-max rescaled per sample (unlike 15_04's "knn" run, which
# normalizes) and one CORE_THRESHOLD applies to every slide. A cell belongs to the
# tumor core when its score is >= CORE_THRESHOLD. Optionally (CORE_CLOSING =
# "geodesic") the thresholded core is then morphologically CLOSED with geodesic
# distances on the spatial k-NN graph, which fills holes and seals the frayed
# boundary (mainly needed for the noisier "knn" score). The core is then made a
# CLOSED region (CORE_FILL_HOLES): every non-core area fully enclosed by core is
# filled, whatever its size. Finally small spatially isolated core fragments are
# removed (MIN_CORE_COMPONENT_SIZE). The result
# is a per-slide spatial domain, analogous to the BANKSY merged domains, stored in
# the `tumor_zone` column ("Tumor Core" / "Non-core").
#
# DE then mirrors 04_01_DEG_DESeq_merged_domains.R: QC-filtered cells are summed
# into sample x zone pseudobulks and DESeq2 (~ condition) is run per zone.
#
# Outputs -> output/R/19_tumor_core_domain/<SCORE_METHOD>[_geodesic_r<radius>][_filled]

# --- Thread limits (set BEFORE BLAS/library load) ----------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS,
           NUMEXPR_NUM_THREADS = N_THREADS)

library(Seurat)
library(DESeq2)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(FNN)
library(Matrix)
library(igraph)
library(dplyr)
library(tidyr)

# Source centralized parameters, palette, and scoring utilities
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

# =============================================================================
# Parameters
# =============================================================================
# --- Tumor score (see 15_04 / 00_utils.R) -----------------------------------
SCORE_METHOD             <- "diffusion"  # "diffusion" | "knn"
KNN_K                    <- 30      # k for "knn" purity and the "diffusion" graph
DIFFUSION_VARIANT        <- "heat"  # "implicit" | "heat" | "ppr"
DIFFUSION_HEAT_STEPS     <- 10      # backward-Euler steps
DIFFUSION_HEAT_ALPHA     <- 0.5       # per-step timestep (total time = steps * alpha)
DIFFUSION_LOCAL_FRACTION <- TRUE    # (K x)/(K 1): comparable [0,1], no min-max
TARGET_CELLTYPE          <- "Tumor"

# --- Tumor core definition ---------------------------------------------------
# Because both scores are local tumor fractions on the same [0,1] scale in every
# sample, a single absolute threshold defines the core consistently across slides:
# CORE_THRESHOLD = 0.25 -> >= 25% of the (kernel-weighted / k-NN) neighborhood is
# tumor. The two scores have different spreads, so revisit the threshold when
# switching SCORE_METHOD.
CORE_THRESHOLD <- 0.25
# Spatial cleanup: two core cells are linked if either is among the other's
# CORE_COMPONENT_K nearest spatial neighbors (within the sample);
# connected components smaller than MIN_CORE_COMPONENT_SIZE cells are relabelled
# Non-core. Set MIN_CORE_COMPONENT_SIZE to 0 to disable.
CORE_COMPONENT_K        <- 10
MIN_CORE_COMPONENT_SIZE <- 50
# Geodesic closing (applied after thresholding, before the component cleanup).
# Distances are shortest paths on the CORE_COMPONENT_K-NN graph with Euclidean
# edge lengths, so they follow the tissue and never jump across empty space.
#   dilation: cells within CLOSING_RADIUS of the thresholded core join it;
#   erosion : of those, keep cells > CLOSING_RADIUS from the nearest cell left
#             outside the dilation.
# Closing only ADDS cells: it fills non-core holes / inlets narrower than about
# 2 * CLOSING_RADIUS and leaves the outer boundary where it was. Radius is in
# spatial coordinate units (microns for Xenium); the per-sample median
# nearest-neighbor spacing is printed for calibration.
CORE_CLOSING   <- "none"   # "none" | "geodesic"
CLOSING_RADIUS <- 50
# Hole filling (applied after closing, before the component cleanup). The tumor
# core is a closed region, so no non-core area may sit inside it. Cells are
# rasterized on a HOLE_FILL_PIXEL grid (spatial units); a pixel is core when a
# majority of its cells are core (empty pixels are background). Background is
# flood-filled (4-connectivity) from the padded grid border; background pixels it
# cannot reach are enclosed by core and become core, regardless of size. Each
# cell then takes the label of its pixel, so a lone non-core cell inside the
# tumor is absorbed. Regions open to the tissue or slide edge (reached through
# empty pixels) are never filled. Pick HOLE_FILL_PIXEL ~2x the cell spacing so
# the core raster is solid (no pinholes between cells).
CORE_FILL_HOLES <- TRUE
HOLE_FILL_PIXEL <- 15

ZONE_CORE    <- "Tumor Core"
ZONE_NONCORE <- "Non-core"
ZONE_COLORS  <- setNames(c(DOMAIN_COLORS[["Tumor Core"]], "grey85"),
                         c(ZONE_CORE, ZONE_NONCORE))

# --- Differential expression -------------------------------------------------
# Zones to test (WT vs DKO within each). Add ZONE_NONCORE for the complement.
DE_ZONES <- c(ZONE_CORE)
# Restrict the pseudobulk to these cell types (NULL = all cells in the zone, as
# in 04_01). E.g. c("Tumor") for tumor-intrinsic DE within the core.
DE_CELLTYPES <- NULL

stopifnot(SCORE_METHOD %in% c("diffusion", "knn"),
          DIFFUSION_VARIANT %in% c("implicit", "heat", "ppr"),
          CORE_CLOSING %in% c("none", "geodesic"), CLOSING_RADIUS > 0,
          is.logical(CORE_FILL_HOLES), HOLE_FILL_PIXEL > 0,
          CORE_THRESHOLD > 0, CORE_THRESHOLD < 1,
          all(DE_ZONES %in% c(ZONE_CORE, ZONE_NONCORE)))

SCORE_LAB <- switch(SCORE_METHOD,
                    diffusion = "Diffusion local tumor fraction",
                    knn       = sprintf("k-NN tumor purity (k = %d)", KNN_K))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "19_tumor_core_domain",
                        paste0(SCORE_METHOD,
                               if (CORE_CLOSING == "geodesic")
                                 sprintf("_geodesic_r%g", CLOSING_RADIUS),
                               if (CORE_FILL_HOLES) "_filled"))
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load data
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells,", nrow(seu), "genes\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
DefaultAssay(seu) <- "RNA"

# =============================================================================
# 2. Tumor score (all cells; QC is applied only at the DE step so that
#    every physically present cell shapes the field, as in 15_04)
# =============================================================================
all_coords <- Embeddings(seu, "spatial")
samples    <- seu$sample
is_tumor   <- seu$Tier1_celltype == TARGET_CELLTYPE

cat(sprintf("Computing tumor score (method = %s, k = %d)...\n", SCORE_METHOD, KNN_K))
tumor_score <- switch(SCORE_METHOD,
  # k-NN purity is already a [0,1] local fraction: no per-sample min-max.
  knn = compute_tumor_score(all_coords, is_tumor, samples, method = "knn",
                            k = KNN_K, normalize = FALSE),
  diffusion = compute_tumor_score(all_coords, is_tumor, samples, method = "diffusion",
                                  diffusion_method = DIFFUSION_VARIANT, k = KNN_K,
                                  heat_steps = DIFFUSION_HEAT_STEPS,
                                  local_fraction = DIFFUSION_LOCAL_FRACTION,
                                  normalize = !DIFFUSION_LOCAL_FRACTION,
                                  alpha = DIFFUSION_HEAT_ALPHA)
)
seu$tumor_score <- tumor_score

# =============================================================================
# 3. Define the tumor core
# =============================================================================
# Drop core components smaller than min_size within one sample. Edges join core
# cells that appear in each other's k-NN lists (computed on ALL cells of the
# sample, so non-core cells in between break the connection).
clean_core_components <- function(coords, in_core, k, min_size) {
  if (min_size <= 0 || sum(in_core) == 0) return(list(core = in_core, comp = NA))
  nn    <- FNN::get.knn(coords, k = k)$nn.index
  from  <- rep(seq_len(nrow(coords)), each = k)
  to    <- as.vector(t(nn))
  keep  <- in_core[from] & in_core[to]
  core_idx <- which(in_core)
  g <- igraph::make_empty_graph(n = nrow(coords), directed = FALSE)
  g <- igraph::add_edges(g, rbind(from[keep], to[keep]))
  g <- igraph::induced_subgraph(g, core_idx)
  comp <- igraph::components(g)
  big  <- comp$csize[comp$membership] >= min_size
  out  <- in_core
  out[core_idx[!big]] <- FALSE
  comp_id <- rep(NA_integer_, nrow(coords))
  comp_id[core_idx[big]] <- comp$membership[big]
  list(core = out, comp = comp_id, n_comp = sum(comp$csize >= min_size),
       n_dropped = sum(!big))
}

# Morphological closing of in_core within one sample using geodesic distances
# on the weighted k-NN graph. Multi-source shortest paths are computed with one
# Dijkstra run from a virtual super-source joined to every source cell at zero
# cost. Returns the closed core (a superset of in_core).
geodesic_closing <- function(coords, in_core, k, radius) {
  n   <- nrow(coords)
  knn <- FNN::get.knn(coords, k = k)
  g <- igraph::graph_from_edgelist(cbind(rep(seq_len(n), each = k),
                                         as.vector(t(knn$nn.index))), directed = FALSE)
  igraph::E(g)$weight <- as.vector(t(knn$nn.dist))
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "min"))
  cat(sprintf("    median nearest-neighbor spacing = %.2f\n", median(knn$nn.dist[, 1])))

  dist_to <- function(src) {
    if (!any(src)) return(rep(Inf, n))
    g2 <- igraph::add_vertices(g, 1)
    g2 <- igraph::add_edges(g2, rbind(n + 1, which(src)), weight = 0)
    as.numeric(igraph::distances(g2, v = n + 1, to = seq_len(n),
                                 algorithm = "dijkstra"))
  }
  dilated <- dist_to(in_core) <= radius
  dist_to(!dilated) > radius
}

# Make the core a closed region within one sample: rasterize, flood-fill the
# background from the grid border, fill every unreachable background pixel, and
# relabel each cell by its pixel. Returns the cell-level core after filling, with
# the filled core region (pixel centers, incl. cell-free enclosed pixels) attached
# as attr "tiles" for plotting. fill = FALSE skips the flood fill (majority raster
# only), so the same function draws an unfilled core.
fill_core_holes <- function(coords, in_core, pixel, fill = TRUE) {
  # Pixel indices with a 1-pixel empty border so the outside is one component
  ix <- floor((coords[, 1] - min(coords[, 1])) / pixel) + 2L
  iy <- floor((coords[, 2] - min(coords[, 2])) / pixel) + 2L
  nx <- max(ix) + 1L
  ny <- max(iy) + 1L
  pid <- (iy - 1L) * nx + ix

  core_frac <- tapply(in_core, pid, mean)
  core_px   <- rep(FALSE, nx * ny)
  core_px[as.integer(names(core_frac))] <- core_frac >= 0.5

  # 4-connected background graph: right and up neighbors
  bg  <- !core_px
  id  <- seq_len(nx * ny)
  col <- (id - 1L) %% nx + 1L
  h   <- id[col < nx & bg & c(bg[-1], FALSE)]
  v   <- id[id <= nx * (ny - 1L) & bg & c(bg[-(1:nx)], rep(FALSE, nx))]
  g   <- igraph::make_graph(as.vector(rbind(c(h, v), c(h + 1L, v + nx))),
                            n = nx * ny, directed = FALSE)
  memb    <- igraph::components(g)$membership
  outside <- if (fill) bg & memb == memb[1] else bg  # pixel 1 is on the empty border

  region <- which(!outside)
  out <- !outside[pid]
  attr(out, "tiles") <- data.frame(
    x = min(coords[, 1]) + ((region - 1L) %% nx - 1L + 0.5) * pixel,
    y = min(coords[, 2]) + ((region - 1L) %/% nx - 1L + 0.5) * pixel)
  out
}

in_core_raw <- !is.na(tumor_score) & tumor_score >= CORE_THRESHOLD
in_core     <- in_core_raw
core_comp   <- rep(NA_integer_, ncol(seu))
comp_log    <- data.frame()

for (s in levels(samples)) {
  idx <- which(samples == s)
  core_s <- in_core_raw[idx]
  if (CORE_CLOSING == "geodesic" && any(core_s)) {
    cat(sprintf("  %s: geodesic closing (radius %g)\n", s, CLOSING_RADIUS))
    core_s <- geodesic_closing(all_coords[idx, , drop = FALSE], core_s,
                               k = CORE_COMPONENT_K, radius = CLOSING_RADIUS)
  }
  core_closed <- core_s
  if (CORE_FILL_HOLES && any(core_s)) {
    core_s <- as.vector(fill_core_holes(all_coords[idx, , drop = FALSE], core_s,
                                        pixel = HOLE_FILL_PIXEL))
  }
  cc  <- clean_core_components(all_coords[idx, , drop = FALSE], core_s,
                               k = CORE_COMPONENT_K, min_size = MIN_CORE_COMPONENT_SIZE)
  in_core[idx]   <- cc$core
  core_comp[idx] <- cc$comp
  comp_log <- rbind(comp_log, data.frame(
    sample = s, n_core_raw = sum(in_core_raw[idx]),
    n_closing_added = sum(core_closed) - sum(in_core_raw[idx]),
    n_holefill_added = sum(core_s & !core_closed),
    n_holefill_removed = sum(core_closed & !core_s),
    n_components_kept = if (is.null(cc$n_comp)) NA else cc$n_comp,
    n_cells_dropped = if (is.null(cc$n_dropped)) 0 else cc$n_dropped))
}

seu$tumor_zone <- factor(ifelse(in_core, ZONE_CORE, ZONE_NONCORE),
                         levels = c(ZONE_CORE, ZONE_NONCORE))
seu$tumor_core_component <- core_comp

cat(sprintf("Tumor core (score >= %.2f, closing = %s, fill holes = %s, components >= %d cells): %d / %d cells\n",
            CORE_THRESHOLD, CORE_CLOSING, CORE_FILL_HOLES, MIN_CORE_COMPONENT_SIZE,
            sum(in_core), length(in_core)))
print(comp_log)

# Per-cell assignment, reusable by downstream scripts
write.csv(data.frame(cell = colnames(seu), sample = seu$sample,
                     condition = seu$condition, celltype = seu$Tier1_celltype,
                     x = all_coords[, 1], y = all_coords[, 2],
                     tumor_score = tumor_score, tumor_zone = seu$tumor_zone,
                     core_component = core_comp),
          file.path(OUTPUT_DIR, "tumor_core_assignment.csv"), row.names = FALSE)

# =============================================================================
# 4. Core diagnostics
# =============================================================================
diag_dir <- file.path(OUTPUT_DIR, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# --- 4a. Per-sample core summary ---
core_summary <- data.frame(sample = seu$sample, condition = seu$condition,
                           core = in_core, tumor = is_tumor) %>%
  group_by(sample, condition) %>%
  summarise(n_cells = n(),
            n_core = sum(core),
            pct_cells_in_core = 100 * mean(core),
            pct_tumor_in_core = 100 * sum(core & tumor) / max(sum(core), 1),
            pct_tumor_cells_captured = 100 * sum(core & tumor) / max(sum(tumor), 1),
            .groups = "drop") %>%
  left_join(comp_log, by = "sample")
cat("\nPer-sample tumor core summary:\n")
print(as.data.frame(core_summary))
write.csv(core_summary, file.path(OUTPUT_DIR, "tumor_core_summary.csv"), row.names = FALSE)

# --- 4b. Score histogram with threshold ---
df_hist <- data.frame(score = tumor_score, sample = samples)
p_hist <- ggplot(df_hist, aes(x = score, fill = sample)) +
  geom_histogram(bins = 60, color = NA) +
  geom_vline(xintercept = CORE_THRESHOLD, linetype = "dashed", color = "black") +
  scale_fill_manual(values = SAMPLE_COLORS, guide = "none") +
  facet_wrap(~ sample, scales = "free_y") +
  labs(x = SCORE_LAB, y = "Cell count",
       title = sprintf("Tumor score per sample (core threshold = %.2f)", CORE_THRESHOLD)) +
  theme_publication()
ggsave(file.path(diag_dir, "tumor_score_histogram_threshold.pdf"),
       p_hist, width = 4 * ceiling(length(levels(samples)) / 2), height = 6)

# --- 4c. Spatial maps: score field and core domain ---
max_x_span <- max(tapply(all_coords[, 1], samples, function(x) diff(range(x))))
max_y_span <- max(tapply(all_coords[, 2], samples, function(x) diff(range(x))))
half_span  <- max(max_x_span, max_y_span) / 2 * 1.05

spatial_panel <- function(s, fill_layer) {
  idx <- which(samples == s)
  df  <- data.frame(x = all_coords[idx, 1], y = all_coords[idx, 2],
                    score = tumor_score[idx], zone = seu$tumor_zone[idx])
  df  <- df[order(df$zone, decreasing = TRUE), ]  # draw core on top
  ggplot(df, aes(x = x, y = y)) + fill_layer +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    ggtitle(s) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
}

p_score <- wrap_plots(lapply(levels(samples), spatial_panel,
  fill_layer = list(geom_point(aes(color = score), size = 0.2, alpha = 0.5),
                    scale_color_viridis_c(option = "magma", limits = c(0, 1),
                                          name = "Tumor\nscore"))),
  nrow = 1, guides = "collect")
# Core drawn as a REGION (HOLE_FILL_PIXEL tiles of the final core, including
# cell-free enclosed pixels) over all cells in grey, so cell-free gaps inside the
# tumor are not mistaken for holes.
zone_panel <- function(s) {
  idx   <- which(samples == s)
  tiles <- attr(fill_core_holes(all_coords[idx, , drop = FALSE], in_core[idx],
                                pixel = HOLE_FILL_PIXEL, fill = CORE_FILL_HOLES),
                "tiles")
  tiles$zone <- factor(ZONE_CORE, levels = names(ZONE_COLORS))
  cells <- data.frame(x = all_coords[idx, 1], y = all_coords[idx, 2],
                      zone = factor(ZONE_NONCORE, levels = names(ZONE_COLORS)))
  ggplot(mapping = aes(x = x, y = y)) +
    geom_point(data = cells, aes(color = zone), size = 0.2, alpha = 0.6) +
    geom_tile(data = tiles, aes(fill = zone), width = HOLE_FILL_PIXEL,
              height = HOLE_FILL_PIXEL) +
    scale_color_manual(values = ZONE_COLORS, name = NULL, drop = FALSE,
                       guide = guide_legend(override.aes = list(size = 3))) +
    scale_fill_manual(values = ZONE_COLORS, name = NULL, drop = FALSE) +
    coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span)) +
    ggtitle(s) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 10))
}
p_zone <- wrap_plots(lapply(levels(samples), zone_panel), nrow = 1, guides = "collect")
p_spatial <- (p_score / p_zone) +
  plot_annotation(title = sprintf("%s and score-defined tumor core (>= %.2f)",
                                  SCORE_LAB, CORE_THRESHOLD),
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_core.pdf"),
       p_spatial, width = 4 * length(levels(samples)), height = 9)
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_core.png"),
       p_spatial, width = 4 * length(levels(samples)), height = 9, dpi = 300)

# --- 4d. Cell-type composition of the core per sample ---
comp_df <- data.frame(sample = seu$sample, celltype = seu$Tier1_celltype)[in_core, ] %>%
  count(sample, celltype) %>%
  group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup()
p_comp <- ggplot(comp_df, aes(x = sample, y = prop, fill = celltype)) +
  geom_col() +
  scale_fill_manual(values = TIER1_COLORS) +
  labs(x = NULL, y = "Fraction of core cells", fill = "Cell type",
       title = "Cell-type composition of the tumor core") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(diag_dir, "tumor_core_celltype_composition.pdf"), p_comp, width = 7, height = 5)
write.csv(comp_df, file.path(diag_dir, "tumor_core_celltype_composition.csv"), row.names = FALSE)

# --- 4e. Agreement with the BANKSY merged domains ---
banksy_domain <- DOMAIN_NAMES[as.character(seu@meta.data[[COL_DOMAIN]])]
overlap_tab <- table(zone = seu$tumor_zone, banksy_domain = banksy_domain)
cat("\nScore-defined zone vs BANKSY merged domain:\n")
print(overlap_tab)
write.csv(as.data.frame.matrix(overlap_tab), file.path(diag_dir, "zone_vs_banksy_domain.csv"))

banksy_core <- banksy_domain == "Tumor Core"
jaccard_df <- data.frame(sample = seu$sample, a = in_core, b = banksy_core) %>%
  group_by(sample) %>%
  summarise(n_score_core = sum(a), n_banksy_core = sum(b),
            n_both = sum(a & b), jaccard = sum(a & b) / max(sum(a | b), 1),
            .groups = "drop")
cat("\nJaccard overlap with BANKSY 'Tumor Core' domain:\n")
print(as.data.frame(jaccard_df))
write.csv(jaccard_df, file.path(diag_dir, "jaccard_vs_banksy_tumor_core.csv"), row.names = FALSE)

# =============================================================================
# 5. Pre-filter cells (same QC as 04_01)
# =============================================================================
keep_cells <- seu$total_counts >= MIN_TOTAL_COUNTS &
              seu$n_genes_by_counts >= MIN_GENES_DETECTED
if (!is.null(DE_CELLTYPES)) keep_cells <- keep_cells & seu$Tier1_celltype %in% DE_CELLTYPES

seu_filt <- subset(seu, cells = colnames(seu)[keep_cells])
cat(sprintf("\nCell filtering: %d -> %d cells (total_counts >= %d, n_genes >= %d%s)\n",
            ncol(seu), ncol(seu_filt), MIN_TOTAL_COUNTS, MIN_GENES_DETECTED,
            if (is.null(DE_CELLTYPES)) "" else
              paste0(", cell types: ", paste(DE_CELLTYPES, collapse = ", "))))

filter_report <- data.frame(
  zone   = levels(seu$tumor_zone),
  before = as.numeric(table(seu$tumor_zone)),
  after  = as.numeric(table(seu_filt$tumor_zone))
)
filter_report$removed     <- filter_report$before - filter_report$after
filter_report$pct_removed <- round(100 * filter_report$removed / filter_report$before, 1)
write.csv(filter_report, file.path(OUTPUT_DIR, "cell_filtering_report.csv"), row.names = FALSE)
print(filter_report)

# =============================================================================
# 6. Pseudobulk aggregation (sample x zone)
# =============================================================================
cat("\nAggregating pseudobulk profiles...\n")
counts_mat <- GetAssayData(seu_filt, assay = "RNA", layer = "counts")

seu_filt$pseudobulk_group <- paste(seu_filt$sample, seu_filt$tumor_zone, sep = "__")
group_counts <- table(seu_filt$pseudobulk_group)

pb_list <- list()
for (g in names(group_counts)) {
  if (group_counts[[g]] >= MIN_CELLS_PER_SAMPLE) {
    cells_in_group <- colnames(seu_filt)[seu_filt$pseudobulk_group == g]
    pb_list[[g]] <- Matrix::rowSums(counts_mat[, cells_in_group, drop = FALSE])
  }
}
pb_mat <- do.call(cbind, pb_list)
colnames(pb_mat) <- names(pb_list)

sample_condition <- setNames(as.character(seu$condition), as.character(seu$sample))
pb_meta <- data.frame(
  group   = names(pb_list),
  sample  = sapply(strsplit(names(pb_list), "__"), `[`, 1),
  zone    = sapply(strsplit(names(pb_list), "__"), `[`, 2),
  n_cells = as.numeric(group_counts[names(pb_list)]),
  stringsAsFactors = FALSE
)
pb_meta$condition <- sample_condition[pb_meta$sample]
rownames(pb_meta) <- pb_meta$group
cat(sprintf("Pseudobulk matrix: %d genes x %d groups (min %d cells per group)\n",
            nrow(pb_mat), ncol(pb_mat), MIN_CELLS_PER_SAMPLE))
write.csv(pb_meta, file.path(OUTPUT_DIR, "pseudobulk_metadata.csv"), row.names = FALSE)

pb_dir <- file.path(OUTPUT_DIR, "pseudobulk_profiles")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)
for (z in unique(pb_meta$zone)) {
  z_groups <- pb_meta$group[pb_meta$zone == z]
  z_pb <- pb_mat[, z_groups, drop = FALSE]
  colnames(z_pb) <- pb_meta[z_groups, "sample"]
  write.csv(as.data.frame(z_pb),
            file.path(pb_dir, paste0("pseudobulk_", gsub("[^A-Za-z0-9]+", "_", z), ".csv")))
}

# =============================================================================
# 7. DESeq2 per zone (WT vs DKO)
# =============================================================================
all_results   <- list()
summary_stats <- data.frame()

for (z in DE_ZONES) {
  cat(sprintf("\n--- Zone: %s ---\n", z))
  z_groups <- pb_meta$group[pb_meta$zone == z]
  z_meta   <- pb_meta[z_groups, ]
  z_counts <- pb_mat[, z_groups, drop = FALSE]

  n_wt  <- sum(z_meta$condition == REFERENCE_CONDITION)
  n_dko <- sum(z_meta$condition == TEST_CONDITION)
  cat(sprintf("  Samples: %d %s, %d %s\n", n_wt, REFERENCE_CONDITION, n_dko, TEST_CONDITION))
  if (n_wt < 2 || n_dko < 2) {
    cat("  SKIPPED: need at least 2 samples per condition.\n")
    next
  }

  # Gene filter 1: detected in >= MIN_GENE_EXPR_FRAC of the zone's cells
  z_cells <- colnames(seu_filt)[seu_filt$tumor_zone == z]
  gene_det_frac <- Matrix::rowSums(counts_mat[, z_cells, drop = FALSE] > 0) / length(z_cells)
  z_counts <- z_counts[gene_det_frac[rownames(z_counts)] >= MIN_GENE_EXPR_FRAC, , drop = FALSE]
  cat(sprintf("  Genes expressed in >=%.0f%% of cells: %d\n",
              MIN_GENE_EXPR_FRAC * 100, nrow(z_counts)))

  # Gene filter 2: >= MIN_COUNTS_PER_SAMPLE counts in >= MIN_SAMPLES_EXPRESSED samples
  keep_genes <- rowSums(z_counts >= MIN_COUNTS_PER_SAMPLE) >= MIN_SAMPLES_EXPRESSED
  z_counts <- z_counts[keep_genes, , drop = FALSE]
  cat(sprintf("  Genes after sample-level filter (>=%d counts in >=%d samples): %d\n",
              MIN_COUNTS_PER_SAMPLE, MIN_SAMPLES_EXPRESSED, nrow(z_counts)))
  if (nrow(z_counts) < 10) {
    cat("  SKIPPED: too few genes after filtering.\n")
    next
  }

  z_meta$condition <- factor(z_meta$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  dds <- DESeqDataSetFromMatrix(countData = round(z_counts), colData = z_meta,
                                design = ~ condition)
  dds <- tryCatch(DESeq(dds, quiet = TRUE), error = function(e) {
    cat(sprintf("  ERROR in DESeq: %s\n", e$message))
    NULL
  })
  if (is.null(dds)) next

  res_df <- as.data.frame(results(dds, contrast = c("condition", TEST_CONDITION,
                                                    REFERENCE_CONDITION)))
  res_df$gene   <- rownames(res_df)
  res_df$domain <- z
  res_df <- res_df[order(res_df$pvalue), ]

  res_df$significance <- "NS"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange > LFC_THRESHOLD] <- "Up in DKO"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange < -LFC_THRESHOLD] <- "Down in DKO"
  all_results[[z]] <- res_df

  n_up     <- sum(res_df$significance == "Up in DKO")
  n_down   <- sum(res_df$significance == "Down in DKO")
  n_tested <- sum(!is.na(res_df$padj))
  cat(sprintf("  DEGs (FDR < %.2f, |LFC| > %.2f): %d up, %d down (of %d tested)\n",
              FDR_THRESHOLD, LFC_THRESHOLD, n_up, n_down, n_tested))

  summary_stats <- rbind(summary_stats, data.frame(
    domain = z, n_samples_WT = n_wt, n_samples_DKO = n_dko,
    n_cells = sum(z_meta$n_cells), n_genes_tested = n_tested,
    n_up_DKO = n_up, n_down_DKO = n_down, stringsAsFactors = FALSE))

  write.csv(res_df, file.path(OUTPUT_DIR,
                              paste0("DEG_", gsub("[^A-Za-z0-9]+", "_", z), ".csv")),
            row.names = FALSE)
}

# =============================================================================
# 8. Save results
# =============================================================================
write.csv(bind_rows(all_results), file.path(OUTPUT_DIR, "DEG_all_tumor_zones.csv"),
          row.names = FALSE)
write.csv(summary_stats, file.path(OUTPUT_DIR, "DEG_summary.csv"), row.names = FALSE)

writeLines(c(
  "19_01 tumor core domain run parameters",
  sprintf("run_time                 = %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("SCORE_METHOD             = %s", SCORE_METHOD),
  sprintf("DIFFUSION_VARIANT        = %s", DIFFUSION_VARIANT),
  sprintf("KNN_K                    = %d", KNN_K),
  sprintf("DIFFUSION_HEAT_STEPS     = %d", DIFFUSION_HEAT_STEPS),
  sprintf("DIFFUSION_HEAT_ALPHA     = %g", DIFFUSION_HEAT_ALPHA),
  sprintf("DIFFUSION_LOCAL_FRACTION = %s", DIFFUSION_LOCAL_FRACTION),
  sprintf("CORE_THRESHOLD           = %g", CORE_THRESHOLD),
  sprintf("CORE_COMPONENT_K         = %d", CORE_COMPONENT_K),
  sprintf("MIN_CORE_COMPONENT_SIZE  = %d", MIN_CORE_COMPONENT_SIZE),
  sprintf("CORE_CLOSING             = %s", CORE_CLOSING),
  sprintf("CLOSING_RADIUS           = %g", CLOSING_RADIUS),
  sprintf("CORE_FILL_HOLES          = %s", CORE_FILL_HOLES),
  sprintf("HOLE_FILL_PIXEL          = %g", HOLE_FILL_PIXEL),
  sprintf("DE_ZONES                 = %s", paste(DE_ZONES, collapse = ", ")),
  sprintf("DE_CELLTYPES             = %s",
          if (is.null(DE_CELLTYPES)) "(all)" else paste(DE_CELLTYPES, collapse = ", ")),
  sprintf("MIN_TOTAL_COUNTS         = %d", MIN_TOTAL_COUNTS),
  sprintf("MIN_GENES_DETECTED       = %d", MIN_GENES_DETECTED),
  sprintf("MIN_CELLS_PER_SAMPLE     = %d", MIN_CELLS_PER_SAMPLE),
  sprintf("MIN_GENE_EXPR_FRAC       = %g", MIN_GENE_EXPR_FRAC),
  sprintf("MIN_COUNTS_PER_SAMPLE    = %d", MIN_COUNTS_PER_SAMPLE),
  sprintf("MIN_SAMPLES_EXPRESSED    = %d", MIN_SAMPLES_EXPRESSED),
  sprintf("FDR_THRESHOLD            = %g", FDR_THRESHOLD),
  sprintf("LFC_THRESHOLD            = %g", LFC_THRESHOLD)
), file.path(OUTPUT_DIR, "run_parameters.txt"))

cat("\n== DEG Summary ==\n")
print(summary_stats)
cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR))
