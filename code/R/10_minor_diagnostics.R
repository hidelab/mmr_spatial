# 10_minor_diagnostics.R
# Quick diagnostic: check for CD4 T cells (Cd4+) and B cells (Cd19+) in the data.
# Count positive cells, spatial plotting, and co-expression with Cd3d/Ptprc for Cd4+ cells.

library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)
library(dplyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "10_minor_diagnostics")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Load Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

DefaultAssay(seu) <- "RNA"

# =============================================================================
# 1. Count Cd4+ and Cd19+ cells
# =============================================================================
counts_mat <- GetAssayData(seu, assay = "RNA", layer = "counts")

cd4_expr <- counts_mat["Cd4", ]
cd19_expr <- counts_mat["Cd19", ]

cd4_pos <- cd4_expr > 0
cd19_pos <- cd19_expr > 0

cat("\n=== Cell counts ===\n")
cat(sprintf("Cd4+ cells:  %d / %d (%.2f%%)\n", sum(cd4_pos), ncol(seu), 100 * sum(cd4_pos) / ncol(seu)))
cat(sprintf("Cd19+ cells: %d / %d (%.2f%%)\n", sum(cd19_pos), ncol(seu), 100 * sum(cd19_pos) / ncol(seu)))
cat(sprintf("Cd4+ & Cd19+: %d\n", sum(cd4_pos & cd19_pos)))

# Breakdown by cell type
cat("\n--- Cd4+ cells by Tier1_celltype ---\n")
print(sort(table(seu$Tier1_celltype[cd4_pos]), decreasing = TRUE))

cat("\n--- Cd19+ cells by Tier1_celltype ---\n")
print(sort(table(seu$Tier1_celltype[cd19_pos]), decreasing = TRUE))

# Breakdown by condition
cat("\n--- Cd4+ cells by condition ---\n")
print(table(seu$condition[cd4_pos]))

cat("\n--- Cd19+ cells by condition ---\n")
print(table(seu$condition[cd19_pos]))

# Breakdown by sample
cat("\n--- Cd4+ cells by sample ---\n")
print(table(seu$sample[cd4_pos]))

cat("\n--- Cd19+ cells by sample ---\n")
print(table(seu$sample[cd19_pos]))

# =============================================================================
# 2. Spatial plots of Cd4 and Cd19 expression
# =============================================================================
cat("\nGenerating spatial plots...\n")

sample_ids <- levels(seu$sample)
sample_ids <- sample_ids[sample_ids %in% unique(seu$sample)]

for (gene in c("Cd4", "Cd19")) {
  cat(sprintf("  Spatial: %s\n", gene))

  plots <- lapply(sample_ids, function(s) {
    sub <- subset(seu, subset = sample == s)
    FeaturePlot(sub, features = gene, reduction = "spatial", pt.size = 1,
                alpha = 0.5, raster = TRUE, raster.dpi = c(300, 300)) +
      scale_color_viridis(option = "magma", name = gene) +
      ggtitle(s) +
      NoAxes() +
      coord_fixed()
  })

  p <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = sprintf("%s expression — all cells", gene),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", gene, ".pdf")),
         p, width = 4 * length(sample_ids), height = 5)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", gene, ".png")),
         p, width = 4 * length(sample_ids), height = 5, dpi = 300)
}

# =============================================================================
# 3. Cd4+ cells: co-expression with Cd3d and Ptprc
# =============================================================================
cat("\n=== Cd4+ cells: co-expression analysis ===\n")

cd3d_expr <- counts_mat["Cd3d", ]
ptprc_expr <- counts_mat["Ptprc", ]

cd4_cells <- colnames(seu)[cd4_pos]
cat(sprintf("Cd4+ cells: %d\n", length(cd4_cells)))
cat(sprintf("  also Cd3d+:  %d (%.1f%%)\n",
            sum(cd3d_expr[cd4_cells] > 0),
            100 * sum(cd3d_expr[cd4_cells] > 0) / length(cd4_cells)))
cat(sprintf("  also Ptprc+: %d (%.1f%%)\n",
            sum(ptprc_expr[cd4_cells] > 0),
            100 * sum(ptprc_expr[cd4_cells] > 0) / length(cd4_cells)))
cat(sprintf("  Cd3d+ & Ptprc+: %d (%.1f%%)\n",
            sum(cd3d_expr[cd4_cells] > 0 & ptprc_expr[cd4_cells] > 0),
            100 * sum(cd3d_expr[cd4_cells] > 0 & ptprc_expr[cd4_cells] > 0) / length(cd4_cells)))

# Breakdown by celltype for Cd4+Cd3d+ (likely true T cells)
cd4_cd3d_pos <- cd4_pos & cd3d_expr > 0
cat("\n--- Cd4+Cd3d+ cells by Tier1_celltype ---\n")
print(sort(table(seu$Tier1_celltype[cd4_cd3d_pos]), decreasing = TRUE))

# =============================================================================
# 4. Spatial plots of Cd4+ cells colored by Cd3d and Ptprc
# =============================================================================
cat("\nSpatial plots of Cd4+ cells...\n")

seu_cd4 <- subset(seu, cells = cd4_cells)
cat(sprintf("Subsetting to %d Cd4+ cells for spatial plotting\n", ncol(seu_cd4)))

for (gene in c("Cd3d", "Ptprc", "Cd4")) {
  cat(sprintf("  Spatial (Cd4+ cells): %s\n", gene))

  plots <- lapply(sample_ids, function(s) {
    sub <- subset(seu_cd4, subset = sample == s)
    if (ncol(sub) == 0) {
      return(ggplot() + ggtitle(s) + theme_void())
    }
    FeaturePlot(sub, features = gene, reduction = "spatial", pt.size = 2,
                alpha = 0.6, raster = TRUE, raster.dpi = c(300, 300)) +
      scale_color_viridis(option = "magma", name = gene) +
      ggtitle(s) +
      NoAxes() +
      coord_fixed()
  })

  p <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = sprintf("%s in Cd4+ cells", gene),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  ggsave(file.path(OUTPUT_DIR, paste0("spatial_Cd4pos_", gene, ".pdf")),
         p, width = 4 * length(sample_ids), height = 5)
}

# =============================================================================
# 5. Summary table
# =============================================================================
summary_df <- data.frame(
  marker = c("Cd4", "Cd19", "Cd4+Cd3d+", "Cd4+Ptprc+", "Cd4+Cd3d+Ptprc+"),
  n_cells = c(
    sum(cd4_pos),
    sum(cd19_pos),
    sum(cd4_pos & cd3d_expr > 0),
    sum(cd4_pos & ptprc_expr > 0),
    sum(cd4_pos & cd3d_expr > 0 & ptprc_expr > 0)
  ),
  pct_total = round(100 * c(
    sum(cd4_pos),
    sum(cd19_pos),
    sum(cd4_pos & cd3d_expr > 0),
    sum(cd4_pos & ptprc_expr > 0),
    sum(cd4_pos & cd3d_expr > 0 & ptprc_expr > 0)
  ) / ncol(seu), 3)
)
write.csv(summary_df, file.path(OUTPUT_DIR, "marker_positive_summary.csv"), row.names = FALSE)

cat("\n\n=== Summary ===\n")
print(summary_df)
cat("\nDone (sections 1-5). Output saved to:", OUTPUT_DIR, "\n")

# =============================================================================
# 6. Klrk1 expression in DC populations
# =============================================================================
# Biological question: Do dendritic cells express Klrk1 (NKG2D receptor)?
# Where are KLRK1+ DCs located spatially, and does their proportion change
# across Banksy spatial domains?

cat("\n\n##############################################################\n")
cat("# 6. Klrk1 expression in DC populations\n")
cat("##############################################################\n\n")

# DC_TYPES, DOMAIN_NAMES, DOMAIN_COLORS defined in 00_parameters.R
DOMAIN_COL <- COL_DOMAIN

# Subset to DCs
dc_cells <- colnames(seu)[seu$Tier1_celltype %in% DC_TYPES]
cat(sprintf("Total DC cells: %d\n", length(dc_cells)))
cat("  By subtype:\n")
print(table(seu$Tier1_celltype[dc_cells]))

# Klrk1 expression in DCs
klrk1_expr <- counts_mat["Klrk1", ]
klrk1_dc_pos <- klrk1_expr[dc_cells] > 0

cat(sprintf("\nKlrk1+ among all DCs: %d / %d (%.2f%%)\n",
            sum(klrk1_dc_pos), length(dc_cells),
            100 * sum(klrk1_dc_pos) / length(dc_cells)))

# Per DC subtype breakdown
cat("\n--- Klrk1+ cells per DC subtype ---\n")
dc_klrk1_summary <- data.frame(
  dc_type = DC_TYPES,
  n_total = NA_integer_,
  n_klrk1_pos = NA_integer_,
  pct_klrk1_pos = NA_real_
)
for (i in seq_along(DC_TYPES)) {
  cells_i <- colnames(seu)[seu$Tier1_celltype == DC_TYPES[i]]
  n_pos <- sum(klrk1_expr[cells_i] > 0)
  dc_klrk1_summary$n_total[i] <- length(cells_i)
  dc_klrk1_summary$n_klrk1_pos[i] <- n_pos
  dc_klrk1_summary$pct_klrk1_pos[i] <- round(100 * n_pos / length(cells_i), 2)
}
print(dc_klrk1_summary)

# By condition
cat("\n--- Klrk1+ DCs by condition ---\n")
for (dc_type in c("All DCs", DC_TYPES)) {
  if (dc_type == "All DCs") {
    cells_i <- dc_cells
  } else {
    cells_i <- colnames(seu)[seu$Tier1_celltype == dc_type]
  }
  for (cond in c(REFERENCE_CONDITION, TEST_CONDITION)) {
    cells_cond <- cells_i[seu$condition[cells_i] == cond]
    n_pos <- sum(klrk1_expr[cells_cond] > 0)
    cat(sprintf("  %s | %s: %d / %d (%.2f%%)\n",
                dc_type, cond, n_pos, length(cells_cond),
                ifelse(length(cells_cond) > 0, 100 * n_pos / length(cells_cond), 0)))
  }
}

# Save summary table
write.csv(dc_klrk1_summary,
          file.path(OUTPUT_DIR, "klrk1_dc_subtype_summary.csv"), row.names = FALSE)

# =============================================================================
# 7. Spatial plots of Klrk1+ DCs
# =============================================================================
cat("\n--- Spatial plots of Klrk1+ DCs ---\n")

seu_dc <- subset(seu, cells = dc_cells)
seu_dc$klrk1_pos <- ifelse(klrk1_expr[colnames(seu_dc)] > 0, "Klrk1+", "Klrk1-")

# 7a. Spatial plot: all DCs colored by Klrk1 status (per sample)
plots_dc_spatial <- lapply(sample_ids, function(s) {
  sub <- subset(seu_dc, subset = sample == s)
  if (ncol(sub) == 0) return(ggplot() + ggtitle(s) + theme_void())
  # Plot Klrk1- first, Klrk1+ on top
  sub$klrk1_pos <- factor(sub$klrk1_pos, levels = c("Klrk1-", "Klrk1+"))
  ord <- order(sub$klrk1_pos)  # Klrk1+ cells on top
  coords <- Embeddings(sub, "spatial")[ord, ]
  df <- data.frame(x = coords[, 1], y = coords[, 2],
                   status = sub$klrk1_pos[ord])
  ggplot(df, aes(x, y, color = status)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_manual(values = c("Klrk1-" = "grey80", "Klrk1+" = "#D32F2F"),
                       name = NULL) +
    ggtitle(s) +
    coord_fixed() +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
})

p_dc_klrk1_spatial <- wrap_plots(plots_dc_spatial, nrow = 1) +
  plot_annotation(
    title = "Klrk1 expression in all DCs (spatial)",
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

ggsave(file.path(OUTPUT_DIR, "spatial_DC_all_Klrk1.pdf"),
       p_dc_klrk1_spatial, width = 4 * length(sample_ids), height = 5)
ggsave(file.path(OUTPUT_DIR, "spatial_DC_all_Klrk1.png"),
       p_dc_klrk1_spatial, width = 4 * length(sample_ids), height = 5, dpi = 300)

# 7b. Spatial plots per DC subtype
for (dc_type in DC_TYPES) {
  dc_type_safe <- gsub("[^A-Za-z0-9]", "_", dc_type)
  seu_dc_sub <- subset(seu_dc, subset = Tier1_celltype == dc_type)
  if (ncol(seu_dc_sub) == 0) next

  plots_sub <- lapply(sample_ids, function(s) {
    sub <- subset(seu_dc_sub, subset = sample == s)
    if (ncol(sub) == 0) return(ggplot() + ggtitle(s) + theme_void())
    sub$klrk1_pos <- factor(sub$klrk1_pos, levels = c("Klrk1-", "Klrk1+"))
    ord <- order(sub$klrk1_pos)
    coords <- Embeddings(sub, "spatial")[ord, ]
    df <- data.frame(x = coords[, 1], y = coords[, 2],
                     status = sub$klrk1_pos[ord])
    ggplot(df, aes(x, y, color = status)) +
      geom_point(size = 1.5, alpha = 0.7) +
      scale_color_manual(values = c("Klrk1-" = "grey80", "Klrk1+" = "#D32F2F"),
                         name = NULL) +
      ggtitle(s) +
      coord_fixed() +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 11))
  })

  p_sub <- wrap_plots(plots_sub, nrow = 1) +
    plot_annotation(
      title = sprintf("Klrk1 expression in %s (spatial)", dc_type),
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))

  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", dc_type_safe, "_Klrk1.pdf")),
         p_sub, width = 4 * length(sample_ids), height = 5)
  ggsave(file.path(OUTPUT_DIR, paste0("spatial_", dc_type_safe, "_Klrk1.png")),
         p_sub, width = 4 * length(sample_ids), height = 5, dpi = 300)
}

# =============================================================================
# 8. Klrk1+ DC proportion across spatial domains
# =============================================================================
cat("\n--- Klrk1+ DC proportion across spatial domains ---\n")

# Assign domain names to DC subset
seu_dc$domain_name <- DOMAIN_NAMES[as.character(seu_dc[[DOMAIN_COL]][, 1])]
seu_dc$domain_name <- factor(seu_dc$domain_name, levels = DOMAIN_ORDER)

# 8a. Proportion of Klrk1+ among ALL DCs per domain (split by condition)
dc_domain_df <- data.frame(
  cell = colnames(seu_dc),
  dc_type = seu_dc$Tier1_celltype,
  condition = seu_dc$condition,
  sample = seu_dc$sample,
  domain = seu_dc$domain_name,
  klrk1_pos = seu_dc$klrk1_pos == "Klrk1+",
  stringsAsFactors = FALSE
)

# Remove cells with no domain assignment
dc_domain_df <- dc_domain_df[!is.na(dc_domain_df$domain), ]

# Summary: all DCs combined
prop_all_dc <- dc_domain_df %>%
  group_by(domain, condition) %>%
  summarise(n_total = n(), n_pos = sum(klrk1_pos),
            pct_pos = 100 * n_pos / n_total, .groups = "drop")

cat("\nKlrk1+ proportion (all DCs) by domain and condition:\n")
print(as.data.frame(prop_all_dc))

# Summary per DC subtype
prop_per_dc <- dc_domain_df %>%
  group_by(dc_type, domain, condition) %>%
  summarise(n_total = n(), n_pos = sum(klrk1_pos),
            pct_pos = 100 * n_pos / n_total, .groups = "drop")

cat("\nKlrk1+ proportion per DC subtype by domain and condition:\n")
print(as.data.frame(prop_per_dc))

# Save tables
write.csv(prop_all_dc,
          file.path(OUTPUT_DIR, "klrk1_DC_all_proportion_by_domain.csv"), row.names = FALSE)
write.csv(prop_per_dc,
          file.path(OUTPUT_DIR, "klrk1_DC_subtype_proportion_by_domain.csv"), row.names = FALSE)

# 8b. Bar plot: % Klrk1+ DCs per domain, split by condition (all DCs)
p_bar_all <- ggplot(prop_all_dc, aes(x = domain, y = pct_pos, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = CONDITION_COLORS, name = "Condition") +
  labs(x = NULL, y = "% Klrk1+ DCs",
       title = "Klrk1+ proportion among all DCs by spatial domain") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "barplot_klrk1_DC_all_by_domain.pdf"),
       p_bar_all, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "barplot_klrk1_DC_all_by_domain.png"),
       p_bar_all, width = 7, height = 5, dpi = 300)

# 8c. Bar plot per DC subtype: faceted
p_bar_sub <- ggplot(prop_per_dc, aes(x = domain, y = pct_pos, fill = condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = CONDITION_COLORS, name = "Condition") +
  facet_wrap(~dc_type, scales = "free_y") +
  labs(x = NULL, y = "% Klrk1+ DCs",
       title = "Klrk1+ proportion per DC subtype by spatial domain") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "barplot_klrk1_DC_subtype_by_domain.pdf"),
       p_bar_sub, width = 10, height = 5)
ggsave(file.path(OUTPUT_DIR, "barplot_klrk1_DC_subtype_by_domain.png"),
       p_bar_sub, width = 10, height = 5, dpi = 300)

# 8d. Sample-level proportions (for statistical robustness)
prop_sample_level <- dc_domain_df %>%
  group_by(dc_type, domain, condition, sample) %>%
  summarise(n_total = n(), n_pos = sum(klrk1_pos),
            pct_pos = 100 * n_pos / n_total, .groups = "drop")

write.csv(prop_sample_level,
          file.path(OUTPUT_DIR, "klrk1_DC_proportion_by_domain_sample.csv"), row.names = FALSE)

# Dot plot with sample-level points (all DCs combined)
prop_sample_all <- dc_domain_df %>%
  group_by(domain, condition, sample) %>%
  summarise(n_total = n(), n_pos = sum(klrk1_pos),
            pct_pos = 100 * n_pos / n_total, .groups = "drop")

p_dot_all <- ggplot(prop_sample_all, aes(x = domain, y = pct_pos, color = condition)) +
  geom_point(aes(shape = condition), size = 3,
             position = position_dodge(width = 0.5), alpha = 0.8) +
  stat_summary(aes(group = condition), fun = mean, geom = "crossbar",
               width = 0.4, position = position_dodge(width = 0.5),
               linewidth = 0.5) +
  scale_color_manual(values = CONDITION_COLORS, name = "Condition") +
  labs(x = NULL, y = "% Klrk1+ DCs",
       title = "Klrk1+ DC proportion by spatial domain (sample-level)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "dotplot_klrk1_DC_all_by_domain_samples.pdf"),
       p_dot_all, width = 7, height = 5)
ggsave(file.path(OUTPUT_DIR, "dotplot_klrk1_DC_all_by_domain_samples.png"),
       p_dot_all, width = 7, height = 5, dpi = 300)

cat("\nDone (sections 6-8). Klrk1 DC analysis complete.\n")
cat("Output saved to:", OUTPUT_DIR, "\n")
