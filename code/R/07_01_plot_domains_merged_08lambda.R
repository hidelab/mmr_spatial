# 07_01_plot_domains_merged_08lambda.R
# Spatial plots highlighting each merged Banksy domain individually.
# For each merged domain, cells belonging to that domain are colored with
# their domain color, all other cells are gray.

library(Seurat)
library(ggplot2)
library(patchwork)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "07_domain_plots")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

DOMAIN_COL <- COL_DOMAIN

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

domains <- sort(unique(as.character(seu@meta.data[[DOMAIN_COL]])))
sample_ids <- levels(seu$sample)

cat(sprintf("Plotting %d merged domains across %d samples\n", length(domains), length(sample_ids)))

# --- Individual domain highlight plots ---
for (dom in domains) {
  dom_name <- DOMAIN_NAMES[dom]
  dom_color <- DOMAIN_COLORS_NUM[dom]
  cat(sprintf("  Domain %s (%s) ...\n", dom, dom_name))

  plots <- lapply(sample_ids, function(s) {
    sub <- subset(seu, subset = sample == s)
    sub$highlight <- ifelse(as.character(sub@meta.data[[DOMAIN_COL]]) == dom, "Domain", "Other")
    sub$highlight <- factor(sub$highlight, levels = c("Other", "Domain"))

    DimPlot(sub, reduction = "spatial", group.by = "highlight", pt.size = 1,
            cols = c("Other" = "grey85", "Domain" = dom_color),
            order = c("Domain"),
            raster = TRUE, raster.dpi = c(300, 300)) +
      ggtitle(s) +
      NoAxes() +
      NoLegend() +
      coord_fixed()
  })

  p <- wrap_plots(plots, nrow = 1) +
    plot_annotation(title = paste0(dom_name, " (Domain ", dom, ")"),
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  ggsave(file.path(OUTPUT_DIR, paste0("spatial_domain_", dom, "_", gsub(" ", "_", dom_name), ".pdf")),
         p, width = 4 * length(sample_ids), height = 5)
}

# --- Combined: all domains colored together ---
seu@meta.data[[DOMAIN_COL]] <- factor(seu@meta.data[[DOMAIN_COL]],
                                       levels = names(DOMAIN_COLORS_NUM))

plots_all <- lapply(sample_ids, function(s) {
  sub <- subset(seu, subset = sample == s)
  DimPlot(sub, reduction = "spatial", group.by = DOMAIN_COL, pt.size = 1,
          cols = DOMAIN_COLORS_NUM, shuffle = TRUE,
          raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle(s) +
    NoAxes() +
    NoLegend() +
    coord_fixed()
})

p_all <- wrap_plots(plots_all, nrow = 1) +
  plot_annotation(title = "Merged Banksy Domains",
                  theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
ggsave(file.path(OUTPUT_DIR, "spatial_all_domains.pdf"),
       p_all, width = 4 * length(sample_ids), height = 5)

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")