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

cat(sprintf("Plotting %d samples, each with %d merged domain panels\n",
            length(sample_ids), length(domains)))

for (s in sample_ids) {
  cat(sprintf("  Sample %s ...\n", s))

  sub <- subset(seu, subset = sample == s)

  # Per-domain highlight panels
  plots <- lapply(domains, function(dom) {
    sub$highlight <- ifelse(as.character(sub@meta.data[[DOMAIN_COL]]) == dom, "Domain", "Other")
    sub$highlight <- factor(sub$highlight, levels = c("Other", "Domain"))

    DimPlot(sub, reduction = "spatial", group.by = "highlight", pt.size = 1,
            cols = c("Other" = "grey85", "Domain" = DOMAIN_COLORS_NUM[dom]),
            order = c("Domain"),
            raster = TRUE, raster.dpi = c(300, 300)) +
      ggtitle(DOMAIN_NAMES[dom]) +
      NoAxes() +
      NoLegend() +
      coord_fixed()
  })

  # All domains together panel
  sub@meta.data[[DOMAIN_COL]] <- factor(sub@meta.data[[DOMAIN_COL]],
                                         levels = names(DOMAIN_COLORS_NUM))
  p_all <- DimPlot(sub, reduction = "spatial", group.by = DOMAIN_COL, pt.size = 1,
                   cols = DOMAIN_COLORS_NUM, shuffle = TRUE,
                   raster = TRUE, raster.dpi = c(300, 300)) +
    ggtitle("All Domains") +
    NoAxes() +
    coord_fixed() +
    guides(color = guide_legend(override.aes = list(size = 3)))

  p <- wrap_plots(c(plots, list(p_all)), nrow = 1) +
    plot_annotation(title = s,
                    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  s_safe <- gsub("[^A-Za-z0-9_]", "_", s)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_subject_", s_safe, ".pdf")),
         p, width = 4 * (length(domains) + 1), height = 5)
}

cat("Done. Plots saved to:", OUTPUT_DIR, "\n")