# 07_02_plot_subjects_merged.R
# Per-sample spatial plots showing each merged domain highlighted individually.
# For each sample, one figure with panels per domain (domain-colored highlight).

library(Seurat)
library(ggplot2)
library(patchwork)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "07_domain_plots", "per_subject")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

DOMAIN_COL <- COL_DOMAIN

# Shared point size for every spatial panel. In raster mode Seurat's on-screen
# dot size depends on each panel's rendered geometry, so we also (a) keep every
# panel legend-free and (b) pin identical coord limits (see build_sample_panels)
# to guarantee the dots render at the same size across all panels.
PT_SIZE <- .3

# Map DOMAIN_COLORS (keyed by name) to numeric IDs for DimPlot
DOMAIN_COLORS_NUM <- setNames(DOMAIN_COLORS[DOMAIN_NAMES[as.character(1:MERGE_K)]],
                              as.character(1:MERGE_K))

# --- Load merged-domain Seurat object ---
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# Order the per-domain panels by the display order defined in 00_parameters.R
# (DOMAIN_ORDER is keyed by name; map to the numeric IDs stored in DOMAIN_COL).
domains_present    <- unique(as.character(seu@meta.data[[DOMAIN_COL]]))
domain_ids_ordered <- names(DOMAIN_NAMES)[match(DOMAIN_ORDER, DOMAIN_NAMES)]
domains            <- domain_ids_ordered[domain_ids_ordered %in% domains_present]
sample_ids <- levels(seu$sample)

cat(sprintf("Plotting %d samples, each with %d merged domain panels\n",
            length(sample_ids), length(domains)))

# ------------------------------------------------------------------------------
# Panel builder: for one sample, return the ordered list of spatial panels:
#   [[1]]  = All Domains (every domain in its original color, leftmost)
#   [[2+]] = one highlight panel per domain in DOMAIN_ORDER
#            (domain cells in their original domain color, all other cells grey)
#
#   show_titles : draw column titles ("All Domains" / domain names) above panels
#   row_label   : if given, print it as a left-side row label on the All Domains
#                 panel (used to label rows in the stacked all-samples figure)
# ------------------------------------------------------------------------------
build_sample_panels <- function(s, show_titles = TRUE, row_label = NULL) {
  # Use `cells =` (not NSE `subset =`) so this works inside a function
  sub <- subset(seu, cells = colnames(seu)[seu$sample == s])

  # Shared coordinate limits for this sample, so every panel is drawn on the
  # exact same canvas under coord_fixed() -> identical on-screen point size.
  emb <- Embeddings(sub, reduction = "spatial")
  xr  <- range(emb[, 1])
  yr  <- range(emb[, 2])

  # Per-domain highlight panels: domain cells keep their original color, rest grey.
  # unname() the domain color so the named vector element does not mangle the
  # "Domain" key into "Domain.<id>" (which would break the color mapping).
  domain_panels <- lapply(domains, function(dom) {
    sub$highlight <- ifelse(as.character(sub@meta.data[[DOMAIN_COL]]) == dom, "Domain", "Other")
    sub$highlight <- factor(sub$highlight, levels = c("Other", "Domain"))

    DimPlot(sub, reduction = "spatial", group.by = "highlight", pt.size = PT_SIZE,
            cols = c("Other" = "grey85", "Domain" = unname(DOMAIN_COLORS_NUM[dom])),
            order = c("Domain")#,
            #raster = TRUE, raster.dpi = c(300, 300)
            ) +
      ggtitle(if (show_titles) DOMAIN_NAMES[dom] else NULL) +
      NoAxes() +
      NoLegend() +
      coord_fixed(xlim = xr, ylim = yr)
  })

  # All-domains panel (leftmost) — every domain in its original color.
  # Kept legend-free (the per-domain panels are titled + colored, so they double
  # as the key) so its drawing area matches the others exactly.
  sub@meta.data[[DOMAIN_COL]] <- factor(sub@meta.data[[DOMAIN_COL]],
                                         levels = names(DOMAIN_COLORS_NUM))
  p_all <- DimPlot(sub, reduction = "spatial", group.by = DOMAIN_COL, pt.size = PT_SIZE,
                   cols = DOMAIN_COLORS_NUM, shuffle = TRUE#,
                   #raster = TRUE, raster.dpi = c(300, 300)
                   ) +
    ggtitle(if (show_titles) "All Domains" else NULL) +
    NoAxes() +
    NoLegend() +
    coord_fixed(xlim = xr, ylim = yr)

  if (!is.null(row_label)) {
    # Keep only the y-axis title as a bold row label for the stacked figure
    p_all <- p_all +
      labs(y = row_label) +
      theme(axis.title.y = element_text(size = 24, face = "bold", angle = 0),
            axis.title.x = element_blank(),
            axis.text    = element_blank(),
            axis.ticks   = element_blank(),
            axis.line    = element_blank())
  }

  c(list(p_all), domain_panels)
}

# --- Per-sample figures (one row of panels per sample) ---
for (s in sample_ids) {
  cat(sprintf("  Sample %s ...\n", s))

  panels <- build_sample_panels(s, show_titles = TRUE, row_label = NULL)
  p <- wrap_plots(panels, nrow = 1) +
    plot_annotation(title = s,
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  s_safe <- gsub("[^A-Za-z0-9_]", "_", s)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_subject_", s_safe, ".pdf")),
         p, width = 4 * (length(domains) + 1), height = 5)
    ggsave(file.path(OUTPUT_DIR, paste0("spatial_subject_", s_safe, ".jpg")),
         p, width = 4 * (length(domains) + 1), height = 5)
}

# --- Combined figure: all samples stacked (one row per sample) ---
cat("  Building stacked all-samples figure ...\n")
rows <- lapply(seq_along(sample_ids), function(i) {
  # Column titles only on the top row; sample name as the left row label
  panels <- build_sample_panels(sample_ids[i],
                                 show_titles = (i == 1),
                                 row_label   = sample_ids[i])
  wrap_plots(panels, nrow = 1)
})
combined <- wrap_plots(rows, ncol = 1)

ggsave(file.path(OUTPUT_DIR, "spatial_subjects_all_stacked.pdf"),
       combined,
       width  = 4 * (length(domains) + 1),
       height = 5 * length(sample_ids),
       limitsize = FALSE)

ggsave(file.path(OUTPUT_DIR, "spatial_subjects_all_stacked.jpg"),
       combined,
       width  = 4 * (length(domains) + 1),
       height = 5 * length(sample_ids),
       limitsize = FALSE)


cat("Done. Plots saved to:", OUTPUT_DIR, "\n")