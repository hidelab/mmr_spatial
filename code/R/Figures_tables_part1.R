# Figures_tables_part1.R
# Publication figures and tables, part 1.
#
# Figure 7a: spatial maps of Tier1 cell types for every sample. WT samples on the
#   top row, DKO on the bottom; within each row samples are ordered by their
#   number of Tier1 "Tumor" cells (most first). All panels share one µm-per-inch
#   scale, so the single scale bar applies to every panel. Saved with a white and
#   with a black background.
# Figure 7b: immune-cell proportions per sample, DKO vs WT (re-styled from
#   11_01's immune_proportions_boxplot.pdf). Brackets mark cell types with a
#   Dirichlet-regression FDR-adjusted P below F7B_SIG_ALPHA; statistics are read
#   from 11_01's outputs, not recomputed.
# Figure 7c: UMAPs of QC-pass CD8+ T cells, WT and DKO side by side, coloured by
#   the res 0.5 clustering used downstream (03_01), clusters labelled C0-C3.
# Figures 7d, 7e: volcanoes of DKO vs WT DEGs within CD8+ T C1 and C2 (03_03
#   pseudobulk DESeq2, re-styled from its volcano_C1/C2.pdf).
# Figure 7f: dot plot of the 03_01 curated markers per CD8+ T cluster (as
#   03_01's dotplot_curated_markers, genes on rows, legend underneath).
# Figures 7g, 7h: DKO vs WT colocalization of CD8+ T C1 (g) and C2 (h) with tumor
#   cells, as 20_01's condition_DKOvsWT_colocalization_Tumor: per-sample log2
#   observed/expected fraction of links to tumor (expected from random cells),
#   bar = condition mean, P = cell-level binomial GLMM (enrichment model, LRT).
#   g stacked over h, shared y axis. Values read from 20_01's outputs, not recomputed.
# Figure 7i: DKO2 spatial map on black (as 20_01's spatial_map_DKO2): all cells,
#   CD8+ T C2, CD8+ T C1 and tumor in a 2 x 2 grid, the overlay of the three
#   underneath; the coloured sub-panel titles are the colour key.
# Figure 7j: spatial domains (banksy_domain_name), same layout, sample order,
#   scale and backgrounds as 7a.
# Figure 7k: volcano of DKO vs WT pseudobulk DESeq2 in the Tumor Core domain
#   (re-styled from 04_02's volcano_merged_domain_3_Tumor Core.pdf); every DEG
#   (padj < FDR_THRESHOLD, |log2FC| > LFC_THRESHOLD, as 04_01) is labelled.
#   Values read from 04_01's output, not recomputed.
# Figure 7l: pathways up in DKO in the Tumor Core, grouped into the three 12_*
#   pathway modules with >= 3 such pathways (re-styled from 12's
#   "Tumor Core_module_up_pathways_barplot"); each block has a bold group title
#   above the module's common genes. Values read from 12's outputs.
# Figure 7 layout preview: all panels on one page (rows a|b, c/d-e|f|g/h|i, j/k|l),
#   assembled from the panel PNGs, black versions where there are two.
#
# Input:  output/R/DKO_banksy_merged.rds (Tier1_celltype, sample, condition,
#         "spatial" reduction in µm, centred per sample by 01_01)
#         output/R/11_spatial_dirichlet/immune_proportions_per_sample.csv
#         output/R/11_spatial_dirichlet/dirichlet_DKOvsWT_immune.csv
#         output/R/03_t_cell_clusters/tcell_cluster_assignments.csv (cells,
#         res 0.5 clusters, UMAP; DKO_CD8T_clusters.rds predates res 0.5)
#         output/R/03_t_cell_clusters/DEG_pseudobulk/per_subtype/DEG_C{1,2}.csv
#         output/R/20_cell_cell_colocalization/CD8T_C{1,2}/
#           focal_target_colocalization_tests.csv, condition_DKOvsWT_glmm.csv
#         output/R/04_DEG_merged_domains/per_domain/DEG_merged_domain_3_Tumor Core.csv
#         output/R/12_DE_pathways/Tumor Core/Tumor Core_{diffPathways,pathway_modules,
#           module_summary}.csv
# Output: output/figures/Figure7a_spatial_celltypes_{white,black}.{pdf,png}
#         output/figures/Figure7b_immune_proportions.{pdf,png}
#         output/figures/Figure7c_CD8T_umap.{pdf,png}
#         output/figures/Figure7de_CD8T_C1_C2_DEG_volcanoes.{pdf,png}
#         output/figures/Figure7f_CD8T_curated_markers.{pdf,png}
#         output/figures/Figure7gh_CD8T_tumor_colocalization.{pdf,png}
#         output/figures/Figure7i_DKO2_spatial_map.{pdf,png}
#         output/figures/Figure7j_spatial_domains_{white,black}.{pdf,png}
#         output/figures/Figure7k_tumor_core_DEG_volcano.{pdf,png}
#         output/figures/Figure7l_tumor_core_pathway_modules.{pdf,png}
#         output/figures/Figure7_layout_preview.pdf
#         output/tables/Table7de_CD8T_C1_C2_DEG.csv
#         output/tables/Table7gh_tumor_colocalization_glmm.csv
#         output/tables/Table7k_tumor_core_DEG.csv
#         output/tables/Table7l_tumor_core_pathways.csv
#         output/tables/Table7a_spatial_celltypes_cells.csv
#         output/tables/Table7b_immune_proportions_dirichlet.csv

library(Seurat)
library(ggplot2)
library(patchwork)
library(scattermore)
library(ggrepel)

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

FIG_DIR <- file.path(PROJECT_DIR, "output", "figures")
TAB_DIR <- file.path(PROJECT_DIR, "output", "tables")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

SEED <- 42

# --- Figure-wide style (Nature / Cell) ----------------------------------------
FIG_FONT     <- "Arial"
FIG_TEXT_PT  <- 6     # legend text, scale-bar label
FIG_TITLE_PT <- 7     # panel titles, legend title
FIG_DPI      <- 600   # raster resolution of the point layers in the PDF
PAGE_WIDTH   <- 7.2   # in, 183 mm two-column page

BG_THEMES <- list(
  white = list(bg = "white", fg = "black"),
  black = list(bg = "black", fg = "white")
)

# ==============================================================================
# 1. Data
# ==============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))

cells <- data.frame(
  sample    = as.character(seu$sample),
  condition = as.character(seu$condition),
  celltype  = factor(as.character(seu$Tier1_celltype), levels = TIER1_ORDER),
  domain    = factor(as.character(seu$banksy_domain_name), levels = DOMAIN_ORDER),
  x         = Embeddings(seu, "spatial")[, 1],
  y         = Embeddings(seu, "spatial")[, 2]
)
stopifnot(!anyNA(cells$celltype), !anyNA(cells$domain))

# CD8+ T subset for 7c/7f: the QC-pass cells and res 0.5 clusters from 03_01
TCELL_RES_COL <- "tcell_clust_res.0.5"
tc_assign <- read.csv(file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters",
                                "tcell_cluster_assignments.csv"), check.names = FALSE)
stopifnot(all(tc_assign$cell %in% colnames(seu)))
tc_levels <- sort(unique(tc_assign[[TCELL_RES_COL]]))
tc_assign$cluster <- factor(paste0("C", tc_assign[[TCELL_RES_COL]]),
                            levels = paste0("C", tc_levels))

DefaultAssay(seu) <- "RNA"
tcell <- DietSeurat(subset(seu, cells = tc_assign$cell), assays = "RNA", dimreducs = NULL)
rm(seu); invisible(gc())
# Same normalization rule as 03_01: log-normalize only if "data" is raw counts
if (identical(GetAssayData(tcell, assay = "RNA", layer = "data"),
              GetAssayData(tcell, assay = "RNA", layer = "counts"))) {
  tcell <- NormalizeData(tcell, assay = "RNA")
}
tcell$cluster <- tc_assign$cluster[match(colnames(tcell), tc_assign$cell)]

# ==============================================================================
# 2. Figure 7a — spatial cell-type maps
# ==============================================================================
F7A_NCOL        <- 4      # panels per row (max samples per condition)
F7A_MAPS_W      <- 3.6    # in, 2 x 4 grid of maps (7a and 7j)
F7A_LEGEND_W    <- 1.05   # in, 7a legend column right of the maps
F7A_WIDTH       <- F7A_MAPS_W + F7A_LEGEND_W
F7A_PTSIZE      <- 1.6    # scattermore point size (pixels at FIG_DPI)
F7A_SCALEBAR_UM <- 500
F7A_PAD         <- 1.03   # margin around the largest tissue

# Sample order: by condition (WT row, DKO row), then tumor-cell count, descending
tumor_n <- tapply(cells$celltype == "Tumor", cells$sample, sum)
f7a_order <- lapply(c(REFERENCE_CONDITION, TEST_CONDITION), function(cond) {
  s <- unique(cells$sample[cells$condition == cond])
  s[order(-tumor_n[s])]
})
names(f7a_order) <- c(REFERENCE_CONDITION, TEST_CONDITION)
stopifnot(all(lengths(f7a_order) <= F7A_NCOL),
          sum(F7A_NCOL - lengths(f7a_order)) >= 1)   # need one free slot for the scale bar

cat("Figure 7a sample order (tumor cells):\n")
for (cond in names(f7a_order))
  cat(sprintf("  %-4s %s\n", cond,
              paste(sprintf("%s (%d)", f7a_order[[cond]], tumor_n[f7a_order[[cond]]]),
                    collapse = ", ")))

# Centre each tissue on its bounding box; one common square window for all panels
for (s in unique(cells$sample)) {
  i <- cells$sample == s
  cells$x[i] <- cells$x[i] - mean(range(cells$x[i]))
  cells$y[i] <- cells$y[i] - mean(range(cells$y[i]))
}
half_span <- F7A_PAD * max(abs(c(cells$x, cells$y)))
lims <- c(-half_span, half_span)

panel_in <- F7A_MAPS_W / F7A_NCOL
panel_px <- round(panel_in * FIG_DPI)
F7A_MAPS_H   <- 2 * (panel_in + 0.13)   # two rows of square panels + titles
F7A_LEGEND_H <- 0.62

# Shuffle draw order so no cell type is systematically hidden under another
set.seed(SEED)
cells <- cells[sample(nrow(cells)), ]

# Table 7a: one row per cell, everything needed to redraw Figure 7a.
#   x_um, y_um   coordinates as plotted (µm, centred on each tissue's bounding box)
# Rows are in figure panel order, then in plotting order (later rows drawn on top).
# Every panel uses the same square window, x and y in [-half_span, half_span].
f7a_slots <- unlist(lapply(f7a_order, function(s) c(s, rep(NA, F7A_NCOL - length(s)))))
t7a <- data.frame(
  cell_id   = rownames(cells),
  sample    = cells$sample,
  condition = cells$condition,
  cell_type = as.character(cells$celltype),
  color_hex = unname(TIER1_COLORS[as.character(cells$celltype)]),
  x_um      = round(cells$x, 2),
  y_um      = round(cells$y, 2)
)
stopifnot(!anyNA(t7a))
t7a <- t7a[order(match(t7a$sample, f7a_slots), seq_len(nrow(t7a))), ]
data.table::fwrite(t7a, file.path(TAB_DIR, "Table7a_spatial_celltypes_cells.csv"))
cat(sprintf("Table 7a: %d cells; plot window +/- %.1f µm, scale bar %d µm\n",
            nrow(t7a), half_span, F7A_SCALEBAR_UM))

# One sample panel, points coloured by column `var` of `cells` (7a cell types,
# 7j spatial domains); `colors` is a named palette in legend order
f7a_panel <- function(s, style, scale_bar = FALSE, var = "celltype",
                      colors = TIER1_COLORS[TIER1_ORDER], legend_title = "Cell type") {
  d <- if (is.na(s)) cells[0, ] else cells[cells$sample == s, ]
  d$group <- d[[var]]
  p <- ggplot(d, aes(x, y)) +
    geom_scattermore(aes(colour = group), pointsize = F7A_PTSIZE,
                     pixels = c(panel_px, panel_px), interpolate = FALSE) +
    # Invisible layer that carries the legend (scattermore draws no keys)
    geom_point(data = data.frame(x = NA_real_, y = NA_real_,
                                 group = factor(names(colors), levels = names(colors))),
               aes(colour = group), na.rm = TRUE) +
    scale_colour_manual(values = colors, limits = names(colors), drop = FALSE,
                        name = legend_title) +
    coord_fixed(xlim = lims, ylim = lims, expand = FALSE) +
    labs(title = if (is.na(s)) " " else s) +
    theme_void(base_family = FIG_FONT) +
    theme(plot.title = element_text(size = FIG_TITLE_PT, colour = style$fg, hjust = 0.5,
                                    margin = margin(b = 2)),
          plot.background = element_rect(fill = style$bg, colour = NA),
          panel.background = element_rect(fill = style$bg, colour = NA),
          plot.margin = margin(1, 1, 1, 1))
  if (scale_bar) {
    # Centred in the panel; label above the bar
    p <- p +
      annotate("segment", x = -F7A_SCALEBAR_UM / 2, xend = F7A_SCALEBAR_UM / 2,
               y = 0, yend = 0, colour = style$fg, linewidth = 0.5, lineend = "butt") +
      annotate("text", x = 0, y = 0.07 * diff(lims),
               label = paste0(F7A_SCALEBAR_UM, " µm"), hjust = 0.5, vjust = 0,
               size = FIG_TEXT_PT / .pt, colour = style$fg, family = FIG_FONT)
  }
  p
}

# Legend either in a column right of the maps (legend_w, one column) or in rows
# under them (legend_h, legend_ncol columns)
f7a_figure <- function(style, var = "celltype", colors = TIER1_COLORS[TIER1_ORDER],
                       legend_title = "Cell type", legend_side = "right",
                       legend_w = F7A_LEGEND_W, legend_ncol = 4, legend_h = F7A_LEGEND_H) {
  slots <- unlist(lapply(f7a_order, function(s) c(s, rep(NA, F7A_NCOL - length(s)))))
  # The scale bar goes in the first empty slot (same coordinates as the maps),
  # so it never overlaps tissue
  sb_slot <- which(is.na(slots))[1]
  panels <- lapply(seq_along(slots), function(k)
    f7a_panel(slots[k], style, scale_bar = isTRUE(k == sb_slot), var = var,
              colors = colors, legend_title = legend_title))
  maps <- wrap_plots(panels, ncol = F7A_NCOL)
  layout <- if (legend_side == "right") {
    (maps | guide_area()) + plot_layout(guides = "collect", widths = c(F7A_MAPS_W, legend_w))
  } else {
    (maps / guide_area()) + plot_layout(guides = "collect", heights = c(F7A_MAPS_H, legend_h))
  }
  layout &
    guides(colour = guide_legend(ncol = if (legend_side == "right") 1 else legend_ncol,
                                 byrow = FALSE, override.aes = list(size = 1.8))) &
    theme(legend.position = legend_side,
          legend.title = element_text(size = FIG_TITLE_PT, colour = style$fg,
                                      family = FIG_FONT),
          legend.title.position = "top",
          legend.text = element_text(size = FIG_TEXT_PT, colour = style$fg,
                                     family = FIG_FONT, margin = margin(l = 1)),
          legend.key.size = unit(7, "pt"),
          legend.key.spacing.x = unit(4, "pt"),
          legend.key.spacing.y = unit(0, "pt"),
          legend.background = element_rect(fill = style$bg, colour = NA),
          legend.key = element_rect(fill = style$bg, colour = NA),
          legend.margin = margin(0, 0, 0, 0),
          plot.background = element_rect(fill = style$bg, colour = NA))
}

f7a_height <- F7A_MAPS_H
for (bg in names(BG_THEMES)) {
  p <- f7a_figure(BG_THEMES[[bg]])
  save_pdf_png(p, file.path(FIG_DIR, paste0("Figure7a_spatial_celltypes_", bg)),
               width = F7A_WIDTH, height = f7a_height)
}
cat(sprintf("Figure 7a: %.1f x %.1f in, saved to %s\n", F7A_WIDTH, f7a_height, FIG_DIR))

# ==============================================================================
# 3. Figure 7b — immune-cell proportions, DKO vs WT
# ==============================================================================
F7B_WIDTH     <- PAGE_WIDTH - F7A_WIDTH - 0.05   # in, rest of row 1 next to 7a
F7B_HEIGHT    <- f7a_height
F7B_SIG_ALPHA <- 0.05    # bracket shown when Dirichlet FDR-adjusted P < this
F7B_DODGE     <- 0.75

DIRICHLET_DIR <- file.path(PROJECT_DIR, "output", "R", "11_spatial_dirichlet")
imm_prop <- read.csv(file.path(DIRICHLET_DIR, "immune_proportions_per_sample.csv"),
                     check.names = FALSE)
imm_stat <- read.csv(file.path(DIRICHLET_DIR, "dirichlet_DKOvsWT_immune.csv"),
                     check.names = FALSE)

imm_prop$celltype  <- factor(imm_prop$celltype, levels = IMMUNE_CELLTYPES)
imm_prop$condition <- factor(imm_prop$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
stopifnot(!anyNA(imm_prop$celltype), setequal(imm_stat$celltype, IMMUNE_CELLTYPES))

star_label <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", "*"))

# One bracket per significant cell type, spanning its WT and DKO boxes
y_top <- tapply(imm_prop$proportion, imm_prop$celltype, max)
f7b_sig <- imm_stat[imm_stat$p.adj < F7B_SIG_ALPHA, c("celltype", "p.adj")]
f7b_sig$x     <- match(f7b_sig$celltype, IMMUNE_CELLTYPES)
f7b_sig$y     <- y_top[f7b_sig$celltype] + 0.03 * max(imm_prop$proportion)
f7b_sig$tick  <- 0.015 * max(imm_prop$proportion)
f7b_sig$label <- star_label(f7b_sig$p.adj)
f7b_sig$x0    <- f7b_sig$x - F7B_DODGE / 4
f7b_sig$x1    <- f7b_sig$x + F7B_DODGE / 4

cat("Figure 7b brackets (Dirichlet FDR-adjusted P):\n")
print(f7b_sig[, c("celltype", "p.adj", "label")], row.names = FALSE)

f7b_figure <- function(style) {
  ggplot(imm_prop, aes(celltype, proportion, fill = condition)) +
    geom_boxplot(position = position_dodge(F7B_DODGE), width = 0.6, outliers = FALSE,
                 linewidth = 0.25, colour = style$fg, alpha = 0.85) +
    geom_point(aes(group = condition), position = position_dodge(F7B_DODGE),
               shape = 21, size = 0.9, stroke = 0.2, colour = style$fg,
               show.legend = FALSE) +
    geom_segment(data = f7b_sig, aes(x = x0, xend = x1, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.25, colour = style$fg) +
    geom_segment(data = rbind(transform(f7b_sig, xt = x0), transform(f7b_sig, xt = x1)),
                 aes(x = xt, xend = xt, y = y, yend = y - tick),
                 inherit.aes = FALSE, linewidth = 0.25, colour = style$fg) +
    geom_text(data = f7b_sig, aes(x = x, y = y, label = label), inherit.aes = FALSE,
              vjust = 0.2, size = FIG_TITLE_PT / .pt, colour = style$fg, family = FIG_FONT) +
    scale_fill_manual(values = CONDITION_COLORS, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    labs(x = NULL, y = "Proportion of immune cells") +
    theme_classic(base_size = FIG_TITLE_PT, base_family = FIG_FONT) +
    theme(axis.line = element_line(linewidth = 0.25, colour = style$fg),
          axis.ticks = element_line(linewidth = 0.25, colour = style$fg),
          axis.ticks.length = unit(1.5, "pt"),
          axis.text = element_text(size = FIG_TEXT_PT, colour = style$fg),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size = FIG_TITLE_PT, colour = style$fg),
          legend.position = "inside",
          legend.position.inside = c(1, 1),
          legend.justification = c(1, 1),
          legend.direction = "horizontal",
          legend.text = element_text(size = FIG_TEXT_PT, colour = style$fg),
          legend.key.size = unit(7, "pt"),
          legend.background = element_rect(fill = NA, colour = NA),
          legend.key = element_rect(fill = NA, colour = NA),
          panel.background = element_rect(fill = style$bg, colour = NA),
          plot.background = element_rect(fill = style$bg, colour = NA),
          plot.margin = margin(4, 4, 2, 4))
}

save_pdf_png(f7b_figure(BG_THEMES$white), file.path(FIG_DIR, "Figure7b_immune_proportions"),
             width = F7B_WIDTH, height = F7B_HEIGHT)
cat(sprintf("Figure 7b: %.1f x %.1f in, saved to %s\n", F7B_WIDTH, F7B_HEIGHT, FIG_DIR))

# Table 7b: Dirichlet-regression statistics behind Figure 7b, in figure order.
# 11_01 fits the "common" parametrization, so each coefficient is the DKO-vs-WT
# change in log(alpha_k), not in the proportion; its sign does not give the
# direction of the proportion change (11_01's "direction" column assumes it does).
# Direction here is taken from the observed mean proportions instead.
mean_prop <- tapply(imm_prop$proportion, list(imm_prop$celltype, imm_prop$condition), mean)
t7b <- imm_stat[match(IMMUNE_CELLTYPES, imm_stat$celltype), ]
t7b <- data.frame(
  cell_type        = t7b$celltype,
  mean_prop_WT     = mean_prop[t7b$celltype, REFERENCE_CONDITION],
  mean_prop_DKO    = mean_prop[t7b$celltype, TEST_CONDITION],
  coef_log_alpha   = t7b$estimate,
  P                = t7b$pval,
  P_adj_BH         = t7b$p.adj,
  direction        = ifelse(mean_prop[t7b$celltype, TEST_CONDITION] >
                              mean_prop[t7b$celltype, REFERENCE_CONDITION],
                            paste("Higher in", TEST_CONDITION), paste("Lower in", TEST_CONDITION)),
  significant      = t7b$p.adj < F7B_SIG_ALPHA
)
names(t7b)[2:3] <- paste0("mean_prop_", c(REFERENCE_CONDITION, TEST_CONDITION))
write.csv(t7b, file.path(TAB_DIR, "Table7b_immune_proportions_dirichlet.csv"), row.names = FALSE)

# ==============================================================================
# 4. Figure 7c — CD8+ T cell UMAP; Figure 7f — curated markers dot plot
# ==============================================================================
# Row 2 = (7c over the 7d | 7e volcanoes) | 7f | 7g over 7h | 7i; widths sum to
# the page width
F7_ROW2_H  <- 3.4                    # in, row 2 height
F7C_COL_W  <- 2.95                   # in, column of 7c and the two volcanoes under it
F7C_SIZE   <- 1.2                    # in, width of each of the WT / DKO UMAPs
F7C_WIDTH  <- F7C_COL_W               # two UMAPs + shared cluster key on the right
F7C_HEIGHT <- 1.45
F7DE_HEIGHT <- F7_ROW2_H - F7C_HEIGHT - 0.05   # 7d (C1) | 7e (C2) volcanoes under 7c
F7F_WIDTH  <- 1.3                    # genes down, clusters across; legend underneath
F7F_HEIGHT <- F7_ROW2_H
F7C_PTSIZE <- 3.0                    # scattermore point size (pixels at FIG_DPI)

# Curated markers, as TCELL_MARKERS in 03_01_t_cell_clusters.R
F7F_MARKERS <- c(
  "Ptprc", "Thy1", "Cd3e", "Cd8a", "Cd8b1", "Tnfrsf9", "Lck", "Nkg7", "Gzmb",
  "Prf1", "Ifng", "Tnfsf10", "Fasl", "Mki67", "Nfatc1",
  "Lag3", "Tigit", "Ctla4", "Cxcr6", "Cxcr3", "Ccr5",
  "Il2ra", "Il2rb", "Il2rg", "Klrk1",
  "Cd226", "Cd244a"
)
F7F_MARKERS <- intersect(F7F_MARKERS, rownames(tcell))

# Same colours as the 03_* figures (assigned to cluster levels in order)
F7C_COLORS <- setNames(tcell_subtype_colors(tc_levels), levels(tc_assign$cluster))

theme_fig <- function() {
  theme_classic(base_size = FIG_TITLE_PT, base_family = FIG_FONT) +
    theme(axis.line = element_line(linewidth = 0.25, colour = "black"),
          axis.ticks = element_line(linewidth = 0.25, colour = "black"),
          axis.ticks.length = unit(1.5, "pt"),
          axis.text = element_text(size = FIG_TEXT_PT, colour = "black"),
          axis.title = element_text(size = FIG_TITLE_PT, colour = "black"),
          plot.title = element_text(size = FIG_TITLE_PT, colour = "black", hjust = 0.5,
                                    margin = margin(b = 3)),
          plot.title.position = "plot",
          legend.title = element_text(size = FIG_TEXT_PT, colour = "black"),
          legend.text = element_text(size = FIG_TEXT_PT, colour = "black"),
          plot.margin = margin(3, 3, 3, 3))
}

# Volcano of a DESeq2 table (columns log2FoldChange, padj, gene, significance with
# "Up in DKO" / "Down in DKO" / "NS", as written by 03_03 and 04_01). Every DEG is
# labelled; DEG counts in the top corners. compact = short corner labels and axis
# titles for narrow panels.
volcano_plot <- function(deg, title, compact = FALSE) {
  deg <- deg[!is.na(deg$padj), ]
  deg$significance <- factor(deg$significance, levels = c("Down in DKO", "NS", "Up in DKO"))
  deg_sig <- deg[deg$significance != "NS", ]
  n_dir <- table(deg_sig$significance)
  corner <- function(dir) if (compact) sprintf("%s (%d)", dir, n_dir[[paste(dir, "in DKO")]]) else
    sprintf("%s in %s (%d)", dir, TEST_CONDITION, n_dir[[paste(dir, "in DKO")]])
  ggplot(deg, aes(log2FoldChange, -log10(padj))) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), linetype = "dashed",
               linewidth = 0.25, colour = "grey60") +
    geom_hline(yintercept = -log10(FDR_THRESHOLD), linetype = "dashed",
               linewidth = 0.25, colour = "grey60") +
    geom_point(data = deg[deg$significance == "NS", ], colour = "grey75", size = 0.6,
               stroke = 0) +
    geom_point(data = deg_sig, aes(colour = significance), size = 1.1, stroke = 0) +
    geom_text_repel(data = deg_sig, aes(label = gene), size = (FIG_TEXT_PT - 1) / .pt,
                    family = FIG_FONT, fontface = "italic", segment.size = 0.2,
                    segment.colour = "grey50", min.segment.length = 0, box.padding = 0.15,
                    point.padding = 0.05, max.overlaps = Inf, seed = SEED) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
             size = (FIG_TEXT_PT - compact) / .pt, family = FIG_FONT,
             colour = CONDITION_COLORS[[REFERENCE_CONDITION]], label = corner("Down")) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.2,
             size = (FIG_TEXT_PT - compact) / .pt, family = FIG_FONT,
             colour = CONDITION_COLORS[[TEST_CONDITION]], label = corner("Up")) +
    scale_colour_manual(values = c("Down in DKO" = CONDITION_COLORS[[REFERENCE_CONDITION]],
                                   "Up in DKO" = CONDITION_COLORS[[TEST_CONDITION]]),
                        guide = "none") +
    scale_x_continuous(labels = scales::label_number(style_negative = "minus")) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
    labs(title = title,
         x = if (compact) expression(log[2]~"FC (DKO/WT)") else
           expression(log[2]~"fold change (DKO/WT)"),
         y = if (compact) expression(-log[10]~"adj."~italic(P)) else
           expression(-log[10]~"adjusted"~italic(P))) +
    theme_fig() +
    theme(plot.margin = margin(10, 3, 3, 3))   # room for the panel letter
}

# --- 7c: UMAP -----------------------------------------------------------------
umap_df <- data.frame(u1 = tc_assign$umaptcell_1, u2 = tc_assign$umaptcell_2,
                      cluster = tc_assign$cluster, condition = tc_assign$condition)
set.seed(SEED)
umap_df <- umap_df[sample(nrow(umap_df)), ]

u_rng <- apply(umap_df[, c("u1", "u2")], 2, range)
u_span <- max(u_rng[2, ] - u_rng[1, ])
arrow_len <- 0.18 * u_span
ax0 <- u_rng[1, ] - 0.04 * u_span        # UMAP axis-arrow origin, bottom-left

# One UMAP per condition on the same axes; arrows on the first, one shared key
u_lim <- list(x = c(ax0[1] - 0.06 * u_span, u_rng[2, 1]), y = c(ax0[2] - 0.06 * u_span, u_rng[2, 2]))
umap_px <- round(F7C_SIZE * FIG_DPI * c(1, diff(u_lim$y) / diff(u_lim$x)))

f7c_panel <- function(cond, arrows = FALSE) {
  p <- ggplot(umap_df[umap_df$condition == cond, ], aes(u1, u2)) +
    geom_scattermore(aes(colour = cluster), pointsize = F7C_PTSIZE, pixels = umap_px,
                     interpolate = FALSE) +
    # Invisible layer that carries the legend (scattermore draws no keys). C0 and
    # C2 overlap in UMAP space, so a key is used instead of on-plot labels.
    geom_point(data = data.frame(u1 = NA_real_, u2 = NA_real_,
                                 cluster = factor(names(F7C_COLORS), levels = names(F7C_COLORS))),
               aes(colour = cluster), na.rm = TRUE)
  if (arrows) p <- p +
    annotate("segment", x = ax0[1], y = ax0[2], xend = ax0[1] + arrow_len, yend = ax0[2],
             linewidth = 0.25, arrow = arrow(length = unit(2, "pt"), type = "closed")) +
    annotate("segment", x = ax0[1], y = ax0[2], xend = ax0[1], yend = ax0[2] + arrow_len,
             linewidth = 0.25, arrow = arrow(length = unit(2, "pt"), type = "closed")) +
    annotate("text", x = ax0[1] + arrow_len / 2, y = ax0[2] - 0.02 * u_span, label = "UMAP1",
             vjust = 1, size = (FIG_TEXT_PT - 1) / .pt, family = FIG_FONT) +
    annotate("text", x = ax0[1] - 0.02 * u_span, y = ax0[2] + arrow_len / 2, label = "UMAP2",
             angle = 90, vjust = 0, size = (FIG_TEXT_PT - 1) / .pt, family = FIG_FONT)
  p + scale_colour_manual(values = F7C_COLORS, limits = names(F7C_COLORS), name = NULL) +
    coord_fixed(xlim = u_lim$x, ylim = u_lim$y, expand = FALSE, clip = "off") +
    labs(title = sprintf("%s (n = %s)", cond, format(sum(umap_df$condition == cond), big.mark = ","))) +
    theme_void(base_family = FIG_FONT) +
    theme(plot.title = element_text(size = FIG_TEXT_PT, hjust = 0.5, face = "bold",
                                    colour = CONDITION_COLORS[[cond]], margin = margin(b = 2)),
          plot.margin = margin(2, 3, 2, 3))
}

f7c <- (f7c_panel(REFERENCE_CONDITION, arrows = TRUE) | f7c_panel(TEST_CONDITION)) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "CD8 T cell clusters",
                  theme = theme(plot.title = element_text(size = FIG_TITLE_PT, family = FIG_FONT,
                                                          hjust = 0.5, margin = margin(b = 2)))) &
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 1.8))) &
  theme(legend.position = "right",
        legend.text = element_text(size = FIG_TEXT_PT, family = FIG_FONT, margin = margin(l = 0.5)),
        legend.key.size = unit(6, "pt"),
        legend.key.spacing.y = unit(1, "pt"),
        legend.margin = margin(0, 0, 0, 4))

save_pdf_png(f7c, file.path(FIG_DIR, "Figure7c_CD8T_umap"),
             width = F7C_WIDTH, height = F7C_HEIGHT)

# --- 7f: curated markers dot plot ---------------------------------------------
# Colour = mean expression min-max scaled per gene across clusters, size = %
# expressing (as bubble_data() in 03_01)
dot_df <- DotPlot(tcell, features = F7F_MARKERS, group.by = "cluster", assay = "RNA")$data
dot_df$expr_scaled <- ave(dot_df$avg.exp, dot_df$features.plot, FUN = function(x) {
  r <- range(x)
  if (diff(r) == 0) x * 0 else (x - r[1]) / diff(r)
})
dot_df$features.plot <- factor(dot_df$features.plot, levels = rev(F7F_MARKERS))
dot_df$id <- factor(dot_df$id, levels = rev(levels(tc_assign$cluster)))   # C3, C2, C1, C0

f7f <- ggplot(dot_df, aes(id, features.plot)) +
  geom_point(aes(size = pct.exp, colour = expr_scaled)) +
  scale_size_area(max_size = 3, limits = c(0, 100), breaks = c(25, 50, 75),
                  name = "% expressing") +
  scale_colour_gradient(low = "grey85", high = "blue", limits = c(0, 1), breaks = c(0, 0.5, 1),
                        name = "Scaled mean expr.") +
  labs(x = NULL, y = NULL, title = "Curated markers") +
  theme_fig() +
  theme(axis.text.y = element_text(face = "italic"),
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.25),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.title.position = "top",
        legend.justification = "center",
        legend.key.size = unit(6, "pt"),
        legend.key.height = unit(6, "pt"),
        legend.spacing.y = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt")) +
  guides(colour = guide_colourbar(order = 1, barwidth = unit(40, "pt"),
                                  barheight = unit(4, "pt")),
         size = guide_legend(order = 2, keywidth = unit(7, "pt")))

save_pdf_png(f7f, file.path(FIG_DIR, "Figure7f_CD8T_curated_markers"),
             width = F7F_WIDTH, height = F7F_HEIGHT)

# --- 7d, 7e: DKO vs WT DEGs within CD8+ T C1 and C2 (03_03 pseudobulk DESeq2)
F7DE_SUBTYPES <- c("C1", "C2")
CD8_DEG_DIR <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters", "DEG_pseudobulk",
                         "per_subtype")
cd8_deg <- lapply(setNames(F7DE_SUBTYPES, F7DE_SUBTYPES), function(st)
  read.csv(file.path(CD8_DEG_DIR, paste0("DEG_", st, ".csv"))))
f7de <- wrap_plots(lapply(F7DE_SUBTYPES, function(st)
  volcano_plot(cd8_deg[[st]], sprintf("CD8+ T %s: DKO vs WT", st), compact = TRUE)), nrow = 1) +
  plot_layout(axis_titles = "collect")
save_pdf_png(f7de, file.path(FIG_DIR, "Figure7de_CD8T_C1_C2_DEG_volcanoes"),
             width = F7C_COL_W, height = F7DE_HEIGHT)

# Table 7de: all genes tested within CD8+ T C1 and C2, DEGs first
t7de <- do.call(rbind, lapply(F7DE_SUBTYPES, function(st) {
  d <- cd8_deg[[st]]
  d <- d[order(d$significance == "NS", d$padj, na.last = TRUE), ]
  data.frame(subtype = paste("CD8+ T", st), gene = d$gene, base_mean = d$baseMean,
             log2FC_DKO_vs_WT = d$log2FoldChange, lfcSE = d$lfcSE, Wald_stat = d$stat,
             P = d$pvalue, P_adj_BH = d$padj, call = d$significance)
}))
write.csv(t7de, file.path(TAB_DIR, "Table7de_CD8T_C1_C2_DEG.csv"), row.names = FALSE)
cat("Figures 7d/7e volcanoes, DEGs per subtype:\n")
print(table(t7de$subtype, t7de$call))

cat(sprintf("Figure 7c: %.2f in wide; 7f: %.2f x %.2f in; %d CD8+ T cells, clusters %s\n",
            F7C_WIDTH, F7F_WIDTH, F7F_HEIGHT, nrow(tc_assign),
            paste(names(table(tc_assign$cluster)), table(tc_assign$cluster), sep = "=", collapse = ", ")))

# ==============================================================================
# 5. Figures 7g, 7h — DKO vs WT colocalization of CD8+ T C1 / C2 with tumor
# ==============================================================================
F7GH_SUBTYPES <- c(g = "1", h = "2")   # panel -> CD8+ T cluster (res 0.5)
F7GH_TARGET   <- "Tumor"
F7GH_WIDTH    <- 1.7                   # in; g stacked over h, row 2 height
F7GH_HEIGHT   <- F7_ROW2_H

COLOC_DIR <- file.path(PROJECT_DIR, "output", "R", "20_cell_cell_colocalization")

# P value as a plotmath label: "P = 0.012" or "P = 8.1 x 10^-4" (as in 20_01)
format_p_math <- function(p) {
  if (is.na(p)) return("italic(P)~'not available'")
  if (p >= 0.001) return(sprintf("italic(P)~'='~'%s'", formatC(p, digits = 2, format = "fg")))
  e <- floor(log10(p)); m <- round(p / 10^e, 1)
  if (m >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("italic(P)~'='~%.1f%%*%%10^%d", m, e)
}

coloc <- lapply(F7GH_SUBTYPES, function(cl) {
  d <- file.path(COLOC_DIR, paste0("CD8T_C", cl))
  list(samples = subset(read.csv(file.path(d, "focal_target_colocalization_tests.csv")),
                        target == F7GH_TARGET),
       glmm = subset(read.csv(file.path(d, "condition_DKOvsWT_glmm.csv")),
                     target == F7GH_TARGET & grepl("^enrichment", model)))
})
stopifnot(all(sapply(coloc, function(x) nrow(x$glmm) == 1)))

# Shared y range so e and f are directly comparable; bracket sits above all points
y_all <- c(0, unlist(lapply(coloc, function(x) x$samples$log2_enrich_vs_random)))
y_br  <- max(y_all) + 0.10 * diff(range(y_all))
y_lim <- c(min(y_all) - 0.08 * diff(range(y_all)), y_br + 0.22 * diff(range(y_all)))

f7gh_panel <- function(panel) {
  d <- coloc[[panel]]$samples
  d$condition <- factor(d$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  ggplot(d, aes(condition, log2_enrich_vs_random)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25, colour = "grey60") +
    stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "errorbar", width = 0.45,
                 linewidth = 0.4, colour = "black") +
    geom_point(aes(fill = condition), shape = 21, size = 1.8, stroke = 0.25, colour = "black",
               position = position_jitter(width = 0.08, height = 0, seed = SEED)) +
    annotate("segment", x = 1, xend = 2, y = y_br, yend = y_br, linewidth = 0.25) +
    annotate("text", x = 1.5, y = y_br, label = format_p_math(coloc[[panel]]$glmm$p_lrt),
             parse = TRUE, vjust = -0.5, size = FIG_TEXT_PT / .pt, family = FIG_FONT) +
    scale_fill_manual(values = CONDITION_COLORS, guide = "none") +
    scale_y_continuous(labels = scales::label_number(style_negative = "minus")) +
    coord_cartesian(xlim = c(0.4, 2.6), ylim = y_lim, expand = FALSE) +
    labs(title = sprintf("CD8+ T C%s \u2013 tumor", F7GH_SUBTYPES[[panel]]),
         x = NULL, y = "Sample-specific colocalization\n(log2 observed/expected)") +
    theme_fig()
}

# g stacked over h; the shared y-axis title gets its own column so it cannot run
# into the panel titles
f7gh_ylab <- ggplot() +
  annotate("text", x = 0, y = 0, angle = 90, size = FIG_TITLE_PT / .pt, family = FIG_FONT,
           label = "Sample-specific colocalization\n(log2 observed/expected)") +
  theme_void()
f7gh <- f7gh_ylab + (f7gh_panel("g") / f7gh_panel("h") & labs(y = NULL)) +
  plot_layout(widths = c(0.36, 1))
save_pdf_png(f7gh, file.path(FIG_DIR, "Figure7gh_CD8T_tumor_colocalization"),
             width = F7GH_WIDTH, height = F7GH_HEIGHT)

# Table 7gh: GLMM behind the P values in 7g/7h (enrichment model only)
t7gh <- do.call(rbind, lapply(names(F7GH_SUBTYPES), function(panel) {
  g <- coloc[[panel]]$glmm
  data.frame(panel = paste0("7", panel),
             subtype = paste0("CD8+ T C", F7GH_SUBTYPES[[panel]]),
             target = g$target,
             model = g$model,
             log_odds_DKO_vs_WT = g$log_odds_DKO_vs_WT,
             se = g$se,
             ci95_low = g$log_odds_DKO_vs_WT - qnorm(0.975) * g$se,
             ci95_high = g$log_odds_DKO_vs_WT + qnorm(0.975) * g$se,
             odds_ratio = exp(g$log_odds_DKO_vs_WT),
             P_LRT = g$p_lrt,
             n_samples_WT = sum(coloc[[panel]]$samples$condition == REFERENCE_CONDITION),
             n_samples_DKO = sum(coloc[[panel]]$samples$condition == TEST_CONDITION),
             n_cells_subtype = sum(coloc[[panel]]$samples$n_focal))
}))
write.csv(t7gh, file.path(TAB_DIR, "Table7gh_tumor_colocalization_glmm.csv"), row.names = FALSE)
cat("Figure 7g/7h GLMM (enrichment model):\n")
print(t7gh[, c("panel", "subtype", "log_odds_DKO_vs_WT", "se", "P_LRT", "n_cells_subtype")], row.names = FALSE)

# ==============================================================================
# 6. Figure 7i — DKO2 spatial map: all cells, CD8+ T C2, CD8+ T C1, tumor, overlay
# ==============================================================================
F7I_SAMPLE      <- "DKO2"
F7I_WIDTH       <- PAGE_WIDTH - F7C_COL_W - F7F_WIDTH - F7GH_WIDTH - 0.15   # in, rest of row 2 (3 gaps);
                                     # small panels 1/4 of it, overlay 1/2
F7I_SCALEBAR_UM <- 500
F7I_POPS   <- c("CD8+ T C2", "CD8+ T C1", "Tumor")
F7I_COLORS <- c("CD8+ T C2" = F7C_COLORS[["C2"]], "CD8+ T C1" = F7C_COLORS[["C1"]],
                "Tumor" = "#FF3B30")   # CD8+ T clusters as in 7c
F7I_PTSIZE <- c("CD8+ T C2" = 2.2, "CD8+ T C1" = 2.2, "Tumor" = 1.2)   # pixels at F7I_DPI
F7I_DPI    <- 300   # = PNG resolution, so sparse T-cell points are not thinned by downscaling

# Populations as in 20_01: Tier1 tumor; QC-pass CD8+ T cells by res 0.5 cluster
df_map <- cells[cells$sample == F7I_SAMPLE, c("x", "y", "celltype")]
tc_cl <- as.character(tc_assign$cluster[match(rownames(df_map), tc_assign$cell)])
df_map$pop <- ifelse(df_map$celltype == "Tumor", "Tumor",
                     ifelse(!is.na(tc_cl) & paste("CD8+ T", tc_cl) %in% F7I_POPS,
                            paste("CD8+ T", tc_cl), NA))
cat("Figure 7i cells:\n"); print(table(df_map$pop, useNA = "ifany"))

x_rng <- range(df_map$x); y_rng <- range(df_map$y)
sb <- data.frame(x = x_rng[1] + 0.04 * diff(x_rng), y = y_rng[1] - 0.04 * diff(y_rng))
sb$xend <- sb$x + F7I_SCALEBAR_UM
y_lim_map <- c(sb$y - 0.02 * diff(y_rng), y_rng[2])
map_asp <- diff(y_lim_map) / diff(x_rng)

map_theme <- theme_void(base_family = FIG_FONT) +
  theme(plot.background = element_rect(fill = "black", colour = NA),
        panel.background = element_rect(fill = "black", colour = NA),
        plot.title = element_text(colour = "white", size = FIG_TEXT_PT - 1, hjust = 0.5,
                                  margin = margin(b = 1)),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5))

f7i_panel <- function(title, pops, width_in, bg_colour = "grey30", scale_bar = FALSE,
                      title_colour = "white") {
  map_px <- round(width_in * F7I_DPI * c(1, map_asp))
  p <- ggplot(df_map, aes(x, y)) +
    geom_scattermore(colour = bg_colour, pointsize = 0.8, pixels = map_px)
  for (pp in rev(F7I_POPS)) {                 # tumor first, CD8+ T C2 on top
    if (!pp %in% pops) next
    p <- p + geom_scattermore(data = df_map[df_map$pop %in% pp, ], colour = F7I_COLORS[[pp]],
                              pointsize = F7I_PTSIZE[[pp]], pixels = map_px)
  }
  if (scale_bar) {
    p <- p + geom_segment(data = sb, aes(x = x, xend = xend, y = y, yend = y),
                          colour = "white", linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = sb$x, y = sb$y,
               label = paste(F7I_SCALEBAR_UM, "µm"), colour = "white",
               size = (FIG_TEXT_PT - 1) / .pt, hjust = 0, vjust = -0.5, family = FIG_FONT)
  }
  p + coord_fixed(xlim = x_rng, ylim = y_lim_map, expand = FALSE) +
    labs(title = title) + map_theme +
    theme(plot.title = element_text(colour = title_colour))
}

small <- function(pp) f7i_panel(pp, pp, F7I_WIDTH / 2, title_colour = F7I_COLORS[[pp]])
f7i <- wrap_plots(
  A = f7i_panel("All cells", character(0), F7I_WIDTH / 2, bg_colour = "grey60", scale_bar = TRUE),
  B = small("CD8+ T C2"), C = small("CD8+ T C1"), D = small("Tumor"),
  E = f7i_panel("Overlay", F7I_POPS, F7I_WIDTH, scale_bar = TRUE),
  design = "AB\nCD\nEE", heights = c(1, 1, 2)) +
  plot_annotation(title = F7I_SAMPLE,
                  theme = theme(plot.background = element_rect(fill = "black", colour = NA),
                                plot.title = element_text(colour = "white", size = FIG_TITLE_PT,
                                                          hjust = 0.5, family = FIG_FONT)))
f7i_height <- 2 * F7I_WIDTH * map_asp + 0.45   # 2 small rows + overlay, + titles
save_pdf_png(f7i, file.path(FIG_DIR, "Figure7i_DKO2_spatial_map"),
             width = F7I_WIDTH, height = f7i_height)
cat(sprintf("Figure 7i: %.2f x %.2f in\n", F7I_WIDTH, f7i_height))

# ==============================================================================
# 7. Figure 7j — spatial domain maps (same layout as 7a)
# ==============================================================================
F7J_LEGEND_H <- 0.32
f7j_height <- F7A_MAPS_H + F7J_LEGEND_H
for (bg in names(BG_THEMES)) {
  p <- f7a_figure(BG_THEMES[[bg]], var = "domain", colors = DOMAIN_COLORS[DOMAIN_ORDER],
                  legend_title = "Spatial domain", legend_side = "bottom",
                  legend_ncol = length(DOMAIN_ORDER),
                  legend_h = F7J_LEGEND_H)
  save_pdf_png(p, file.path(FIG_DIR, paste0("Figure7j_spatial_domains_", bg)),
               width = F7A_MAPS_W, height = f7j_height)
}
cat(sprintf("Figure 7j: %.1f x %.1f in\n", F7A_MAPS_W, f7j_height))

# ==============================================================================
# 8. Figure 7k — DKO vs WT DEGs in the Tumor Core domain (volcano)
# ==============================================================================
F7K_DOMAIN <- "3"                      # Tumor Core (DOMAIN_NAMES)
F7K_WIDTH  <- 3.0                      # under 7j in row 3
F7K_HEIGHT <- 1.45

deg <- read.csv(file.path(PROJECT_DIR, "output", "R", "04_DEG_merged_domains", "per_domain",
                          sprintf("DEG_merged_domain_%s_%s.csv", F7K_DOMAIN,
                                  DOMAIN_NAMES[[F7K_DOMAIN]])))
deg <- deg[!is.na(deg$padj), ]
n_dir <- table(factor(deg$significance, levels = c("Down in DKO", "NS", "Up in DKO")))
f7k <- volcano_plot(deg, sprintf("%s: %s vs %s", DOMAIN_NAMES[[F7K_DOMAIN]], TEST_CONDITION,
                                 REFERENCE_CONDITION))
save_pdf_png(f7k, file.path(FIG_DIR, "Figure7k_tumor_core_DEG_volcano"),
             width = F7K_WIDTH, height = F7K_HEIGHT)

# Table 7k: all genes tested in the Tumor Core domain, DEGs first
t7k <- deg[order(deg$significance == "NS", deg$padj),
           c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
             "significance")]
names(t7k) <- c("gene", "base_mean", "log2FC_DKO_vs_WT", "lfcSE", "Wald_stat", "P",
                "P_adj_BH", "call")
write.csv(t7k, file.path(TAB_DIR, "Table7k_tumor_core_DEG.csv"), row.names = FALSE)
cat(sprintf("Figure 7k: %d genes tested, %d up / %d down in DKO\n", nrow(deg),
            n_dir[["Up in DKO"]], n_dir[["Down in DKO"]]))

# ==============================================================================
# 8b. Figure 7l — Tumor Core pathways up in DKO, grouped by pathway module
# ==============================================================================
# Module id (12_* "_pathway_modules.csv") -> group title, in plotting order
F7L_GROUPS <- c("6"  = "Killer cell-mediated cytotoxicity",
                "8"  = "T cell activation",
                "14" = "IL2/IL15 signaling")
F7L_WIDTH      <- PAGE_WIDTH - F7A_MAPS_W - 0.05   # right of 7j/7k in row 3
F7L_HEIGHT     <- NA                    # set below to the 7j + 7k column height
F7L_NAME_CHARS <- 60                    # pathway-name line length (all names fit one line)
F7L_CUT_CHARS  <- 44                    # where the one very long KEGG name is cut
F7L_BOX_W      <- 0.9                   # in, common-genes box column
F7L_BAR_W      <- 0.6                   # in, bar panel

PATH_DIR <- file.path(PROJECT_DIR, "output", "R", "12_DE_pathways", "Tumor Core")
pw  <- read.csv(file.path(PATH_DIR, "Tumor Core_diffPathways.csv"), row.names = 1)
pwm <- read.csv(file.path(PATH_DIR, "Tumor Core_pathway_modules.csv"))
pws <- read.csv(file.path(PATH_DIR, "Tumor Core_module_summary.csv"))
pw$module <- pwm$module[match(pw$standard_name, pwm$standard_name)]

# Same selection as 12_pathway_activity.R: significant (adj.P.Val < FDR_THRESHOLD),
# up in DKO, in modules with >= 3 such pathways
pw_up <- pw[!is.na(pw$adj.P.Val) & pw$adj.P.Val < FDR_THRESHOLD & pw$logFC > 0, ]
n_mod <- table(pw_up$module)
pw_up <- pw_up[pw_up$module %in% names(n_mod)[n_mod >= 3], ]
stopifnot(setequal(unique(pw_up$module), as.integer(names(F7L_GROUPS))))
pw_up$neg_log_p <- -log10(pw_up$P.Value)

# Pathway label: name without its database prefix, on one line (the database is
# kept in Table 7l).
# The one very long KEGG name loses its boilerplate prefix and is cut at a step.
pw_up$source <- sub("^\\(([^)]*)\\).*$", "\\1", pw_up$display_name)
pw_up$name   <- sub("^\\([^)]*\\) ", "", pw_up$display_name)
pw_up$name   <- sub("Pathway Definition from KEGG: ", "", pw_up$name, fixed = TRUE)
long <- nchar(pw_up$name) > F7L_NAME_CHARS
pw_up$name[long] <- paste0(sub(" -> [^>]*$", "", substr(pw_up$name[long], 1, F7L_CUT_CHARS)),
                           " -> \u2026")
wrap_name <- function(x) {
  w <- strsplit(x, " ")[[1]]
  paste(vapply(split(w, wrap_word_lines(w, F7L_NAME_CHARS)), paste, character(1), collapse = " "),
        collapse = "<br>")
}
pw_up$label <- vapply(pw_up$name, wrap_name, character(1))

x_max <- max(pw_up$neg_log_p) * 1.05

# Bars of one group, with the group title in bold above them
f7l_bars <- function(mod, last = FALSE) {
  d <- pw_up[pw_up$module == as.integer(mod), ]
  d <- d[order(d$neg_log_p), ]
  d$label <- factor(d$label, levels = d$label)
  ggplot(d, aes(neg_log_p, label)) +
    geom_col(width = 0.7, fill = CONDITION_COLORS[[TEST_CONDITION]]) +
    scale_x_continuous(limits = c(0, x_max), expand = expansion(mult = c(0, 0))) +
    labs(title = F7L_GROUPS[[mod]], y = NULL,
         x = if (last) expression(-log[10]~italic(P)) else NULL) +
    theme_fig() +
    theme(plot.title = element_text(size = FIG_TEXT_PT, face = "bold", hjust = 0,
                                    margin = margin(b = 1.5)),
          plot.title.position = "plot",
          axis.text.y = ggtext::element_markdown(size = FIG_TEXT_PT - 1, hjust = 1,
                                                 lineheight = 1.05),
          axis.ticks.y = element_blank(),
          axis.text.x = if (last) element_text(size = FIG_TEXT_PT) else element_blank(),
          axis.ticks.x = if (last) element_line(linewidth = 0.25) else element_blank(),
          plot.margin = margin(2, 2, if (last) 1 else 5, 2))
}

# The group's common genes (12_* module core genes) in a box right of its bars
f7l_genes <- function(mod) {
  genes <- strsplit(pws$core_genes[pws$module == as.integer(mod)], ";")[[1]]
  ggplot() +
    ggtext::geom_textbox(aes(x = 0, y = 0.5,
                             label = paste0("<span style='color:grey35'>Common genes</span><br>*",
                                            paste(genes, collapse = ", "), "*")),
                         width = unit(1, "npc"), hjust = 0, vjust = 0.5, halign = 0,
                         size = (FIG_TEXT_PT - 1) / .pt, family = FIG_FONT, lineheight = 1.1,
                         box.colour = "grey60", box.size = 0.25, box.r = unit(1.5, "pt"),
                         box.padding = unit(c(2, 3, 2, 3), "pt"), fill = "grey97") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() + theme(plot.margin = margin(2, 0, 2, 4))
}

f7l_n <- as.integer(table(pw_up$module)[names(F7L_GROUPS)])
f7l_plots <- unlist(lapply(seq_along(F7L_GROUPS), function(k) {
  mod <- names(F7L_GROUPS)[k]
  list(f7l_bars(mod, last = k == length(F7L_GROUPS)), f7l_genes(mod))
}), recursive = FALSE)
f7l <- wrap_plots(f7l_plots, ncol = 2, widths = c(F7L_BAR_W, F7L_BOX_W), heights = f7l_n) +
  plot_annotation(title = sprintf("Pathways up in %s, %s", TEST_CONDITION,
                                  tolower(DOMAIN_NAMES[[F7K_DOMAIN]])),
                  theme = theme(plot.title = element_text(size = FIG_TITLE_PT, family = FIG_FONT,
                                                          hjust = 0.5, margin = margin(b = 2))))
F7L_HEIGHT <- f7j_height + 0.05 + F7K_HEIGHT
save_pdf_png(f7l, file.path(FIG_DIR, "Figure7l_tumor_core_pathway_modules"),
             width = F7L_WIDTH, height = F7L_HEIGHT)

# Table 7l: the pathways in 7l with their group, statistics and module common genes
t7l <- do.call(rbind, lapply(names(F7L_GROUPS), function(mod) {
  d <- pw_up[pw_up$module == as.integer(mod), ]
  d <- d[order(d$P.Value), ]
  data.frame(group = F7L_GROUPS[[mod]], module = as.integer(mod), pathway = d$display_name,
             standard_name = d$standard_name, logFC_DKO_vs_WT = d$logFC, t = d$t,
             P = d$P.Value, P_adj_BH = d$adj.P.Val,
             common_genes = gsub(";", ", ", pws$core_genes[pws$module == as.integer(mod)]))
}))
write.csv(t7l, file.path(TAB_DIR, "Table7l_tumor_core_pathways.csv"), row.names = FALSE)
cat(sprintf("Figure 7l: %d pathways in %d groups, %.2f x %.2f in\n", nrow(t7l),
            length(F7L_GROUPS), F7L_WIDTH, F7L_HEIGHT))

# ==============================================================================
# 9. Figure 7 layout preview — panel PNGs placed on one page at their own sizes
# ==============================================================================
F7_ROWS <- list(   # rows of columns; a column's panels are stacked top to bottom
  list(c(a = "Figure7a_spatial_celltypes_black"), c(b = "Figure7b_immune_proportions")),
  list(c(c = "Figure7c_CD8T_umap", d = "Figure7de_CD8T_C1_C2_DEG_volcanoes"),
       c(f = "Figure7f_CD8T_curated_markers"),
       c(g = "Figure7gh_CD8T_tumor_colocalization"), c(i = "Figure7i_DKO2_spatial_map")),
  list(c(j = "Figure7j_spatial_domains_black", k = "Figure7k_tumor_core_DEG_volcano"),
       c(l = "Figure7l_tumor_core_pathway_modules")))
F7_DARK   <- c("a", "i", "j")              # panels on black: white letters
# Second letters inside two-panel files, at (x, y) fractions of the file's size
F7_EXTRA  <- list(d = list(e = c(0.53, 0)), g = list(h = c(0, 0.5)))
F7_MARGIN <- 0.1                           # in, page margin
F7_GAP    <- 0.05                          # in, between rows and panels

img_in <- function(m) c(w = ncol(m), h = nrow(m)) / 300    # save_pdf_png writes 300 dpi
prev_img <- lapply(F7_ROWS, function(r) lapply(r, function(col)
  lapply(col, function(f) png::readPNG(file.path(FIG_DIR, paste0(f, ".png"))))))
col_w <- function(col) max(sapply(col, function(m) img_in(m)[["w"]]))
col_h <- function(col) sum(sapply(col, function(m) img_in(m)[["h"]])) + F7_GAP * (length(col) - 1)
row_h <- sapply(prev_img, function(r) max(sapply(r, col_h)))
row_w <- sapply(prev_img, function(r) sum(sapply(r, col_w)) + F7_GAP * (length(r) - 1))
page_w <- max(row_w) + 2 * F7_MARGIN
page_h <- sum(row_h) + F7_GAP * (length(row_h) - 1) + 2 * F7_MARGIN

panel_letter <- function(l, x, y, dark) {
  grid::grid.text(l, x = unit(x + 0.02, "in"), y = unit(y - 0.02, "in"), just = c("left", "top"),
                  gp = grid::gpar(fontfamily = FIG_FONT, fontface = "bold", fontsize = 8,
                                  col = if (dark) "white" else "black"))
}
cairo_pdf(file.path(FIG_DIR, "Figure7_layout_preview.pdf"), width = page_w, height = page_h)
grid::grid.newpage()
y_row <- page_h - F7_MARGIN
for (k in seq_along(prev_img)) {
  x <- F7_MARGIN
  for (col in prev_img[[k]]) {
    y <- y_row
    for (p in names(col)) {
      sz <- img_in(col[[p]])
      grid::grid.raster(col[[p]], x = unit(x, "in"), y = unit(y, "in"), width = unit(sz[["w"]], "in"),
                        height = unit(sz[["h"]], "in"), just = c("left", "top"))
      panel_letter(p, x, y, p %in% F7_DARK)
      for (l in names(F7_EXTRA[[p]]))
        panel_letter(l, x + F7_EXTRA[[p]][[l]][1] * sz[["w"]],
                     y - F7_EXTRA[[p]][[l]][2] * sz[["h"]], FALSE)
      y <- y - sz[["h"]] - F7_GAP
    }
    x <- x + col_w(col) + F7_GAP
  }
  y_row <- y_row - row_h[k] - F7_GAP
}
invisible(dev.off())
cat(sprintf("Figure 7 preview: %.2f x %.2f in; row widths %s in, row heights %s in\n", page_w, page_h,
            paste(sprintf("%.2f", row_w), collapse = ", "), paste(sprintf("%.2f", row_h), collapse = ", ")))

writeLines(capture.output(sessionInfo()), file.path(FIG_DIR, "sessionInfo_part1.txt"))
