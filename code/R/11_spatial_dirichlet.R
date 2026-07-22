# 11_spatial_dirichlet.R
# Dirichlet multinomial regression for differential immune cell type proportions
# between DKO and WT conditions.

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)
library(DirichletReg)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "11_spatial_dirichlet")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Parameters (from 00_parameters.R)
# =============================================================================

# Immune cell types to include in composition analysis
CELLTYPES_OF_INTEREST <- IMMUNE_CELLTYPES

# =============================================================================
# 1. Load Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# QC filter
keep_qc <- seu$total_counts >= MIN_TOTAL_COUNTS &
            seu$n_genes_by_counts >= MIN_GENES_DETECTED
seu <- subset(seu, cells = colnames(seu)[keep_qc])
cat(sprintf("After QC: %d cells\n", ncol(seu)))

# Keep only immune cell types of interest
seu <- subset(seu, subset = Tier1_celltype %in% CELLTYPES_OF_INTEREST)
seu$Tier1_celltype <- droplevels(factor(seu$Tier1_celltype, levels = CELLTYPES_OF_INTEREST))
cat(sprintf("After filtering to immune cell types: %d cells\n", ncol(seu)))

# =============================================================================
# 2. Build count table: sample x cell type
# =============================================================================
count_df <- seu@meta.data %>%
  group_by(sample, condition, Tier1_celltype, .drop = TRUE) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Tier1_celltype, values_from = n, values_fill = 0) %>%
  as.data.frame()




ct_cols <- CELLTYPES_OF_INTEREST

cat("\nCell counts per sample:\n")
print(count_df)

# =============================================================================
# 3. Prepare proportions with pseudocount
# =============================================================================
count_df$total <- rowSums(count_df[, ct_cols])

# Add small pseudocount (proportional to sample size) to avoid zeros
prop_df <- count_df
for (ct in ct_cols) {
  prop_df[[ct]] <- count_df[[ct]] + (count_df$total / sum(count_df$total))
}
prop_df$total_adj <- rowSums(prop_df[, ct_cols])
for (ct in ct_cols) {
  prop_df[[ct]] <- prop_df[[ct]] / prop_df$total_adj
}

# =============================================================================
# 4. Dirichlet regression: DKO vs WT
# =============================================================================
cat("\n\n================================================================\n")
cat("=== Dirichlet regression: DKO vs WT (immune cell composition) ===\n")
cat("================================================================\n")

dr_data <- as.data.frame(prop_df)
dr_data$condition <- factor(dr_data$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))

AL <- DR_data(dr_data[, ct_cols])

fit <- DirichReg(AL ~ condition, data = dr_data)
fit_null <- DirichReg(AL ~ 1, data = dr_data)

cat("\nLikelihood ratio test (condition effect):\n")
lr_test <- anova(fit_null, fit)
print(lr_test)

# Extract p-values
u <- summary(fit)
coef_mat <- u$coef.mat
pvals_rows <- grep("Intercept", rownames(coef_mat), invert = TRUE)
pvals <- coef_mat[pvals_rows, 4]
estimates <- coef_mat[pvals_rows, 1]

results_df <- data.frame(
  celltype = u$varnames,
  estimate = estimates,
  pval = pvals,
  stringsAsFactors = FALSE
)
results_df$p.adj <- p.adjust(results_df$pval, method = "fdr")
results_df$direction <- ifelse(results_df$estimate > 0, "Up in DKO", "Up in WT")
results_df$significant <- results_df$p.adj < FDR_THRESHOLD
results_df <- results_df[order(results_df$pval), ]

cat("\nResults (positive estimate = higher proportion in DKO):\n")
print(results_df)

# Fitted proportions per condition
fitted_vals <- fit$fitted.values$mu
fitted_vals <- as.data.frame(fitted_vals)
fitted_vals$condition <- dr_data$condition
fitted_vals$sample <- dr_data$sample

fitted_summary <- fitted_vals %>%
  group_by(condition) %>%
  summarise(across(all_of(ct_cols), mean), .groups = "drop")

cat("\nFitted mean proportions:\n")
print(as.data.frame(fitted_summary))

# Save results
write.csv(results_df, file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_immune.csv"), row.names = FALSE)
write.csv(count_df, file.path(OUTPUT_DIR, "immune_counts_per_sample.csv"), row.names = FALSE)

# Save proportions table
prop_out <- count_df %>%
  pivot_longer(cols = all_of(ct_cols), names_to = "celltype", values_to = "n") %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
write.csv(prop_out, file.path(OUTPUT_DIR, "immune_proportions_per_sample.csv"), row.names = FALSE)

# =============================================================================
# 5. Visualization: Observed proportions stacked bar
# =============================================================================
cat("\nGenerating plots...\n")

obs_prop <- count_df %>%
  pivot_longer(cols = all_of(ct_cols), names_to = "celltype", values_to = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = CELLTYPES_OF_INTEREST))

# Add condition
obs_prop <- obs_prop %>%
  left_join(distinct(count_df[, c("sample", "condition")]), by = c("sample","condition"))

p_stack <- ggplot(obs_prop, aes(x = sample, y = prop, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = TIER1_COLORS, name = "Cell type") +
  labs(title = "Immune cell composition per sample",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "immune_composition_stacked_bar.pdf"), p_stack, width = 10, height = 6)

# =============================================================================
# 6. Visualization: Boxplot of proportions by condition (with significance)
# =============================================================================

# Prepare significance labels from Dirichlet results
sig_labels <- results_df %>%
  mutate(
    label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.1   ~ "^",
      TRUE          ~ ""
    ),
    celltype = factor(celltype, levels = CELLTYPES_OF_INTEREST)
  )

# Get max proportion per celltype for placing asterisks
y_positions <- obs_prop %>%
  group_by(celltype) %>%
  summarise(y_max = max(prop) * 1.05, .groups = "drop")

sig_labels <- sig_labels %>%
  left_join(y_positions, by = "celltype") %>%
  filter(label != "")

p_box <- ggplot(obs_prop, aes(x = celltype, y = prop, fill = condition)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  geom_point(aes(group = condition), position = position_dodge(width = 0.75),
             size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Immune cell proportions: DKO vs WT",
       subtitle = "Dirichlet regression adj. p-value: *** <0.001, ** <0.01, * <0.05, ^ <0.1",
       x = NULL, y = "Proportion of immune cells") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

if (nrow(sig_labels) > 0) {
  p_box <- p_box +
    geom_text(data = sig_labels, aes(x = celltype, y = y_max, label = label),
              inherit.aes = FALSE, size = 5, vjust = 0)
}

ggsave(file.path(OUTPUT_DIR, "immune_proportions_boxplot.pdf"), p_box, width = 10, height = 5)

# =============================================================================
# 6b. Visualization: Collapsed stacked bar (DKO vs WT)
# =============================================================================
collapsed_prop <- count_df %>%
  group_by(condition) %>%
  summarise(across(all_of(ct_cols), sum), .groups = "drop") %>%
  pivot_longer(cols = all_of(ct_cols), names_to = "celltype", values_to = "n") %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = CELLTYPES_OF_INTEREST))

p_collapsed <- ggplot(collapsed_prop, aes(x = condition, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = TIER1_COLORS, name = "Cell type") +
  labs(title = "Immune cell composition: DKO vs WT (all samples pooled)",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "immune_composition_collapsed.pdf"), p_collapsed, width = 7, height = 6)

# =============================================================================
# 7. Visualization: Dot plot of Dirichlet results
# =============================================================================
results_df$celltype <- factor(results_df$celltype, levels = rev(CELLTYPES_OF_INTEREST))

p_dot <- ggplot(results_df, aes(x = estimate, y = celltype,
                                 size = -log10(p.adj),
                                 color = estimate)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient2(low = "#4575B4", mid = "grey70", high = "#D73027",
                        midpoint = 0, name = "Estimate") +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
  labs(title = "Dirichlet regression: DKO vs WT",
       subtitle = "Positive estimate = enriched in DKO",
       x = "Estimate (DKO vs WT)", y = NULL) +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_dotplot.pdf"), p_dot, width = 7, height = 6)

# #############################################################################
# EPITHELIAL & FIBROBLAST: Dirichlet regression DKO vs WT
# #############################################################################
cat("\n\n================================================================\n")
cat("=== Dirichlet regression: DKO vs WT (Epithelial & Fibroblasts) ===\n")
cat("================================================================\n")

EPI_FIB_CELLTYPES <- EPITHELIAL_STROMAL_CELLTYPES

# Reload full object for non-immune cell types
seu_ef <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
seu_ef$condition <- gsub("^KO$", TEST_CONDITION, seu_ef$condition)
seu_ef$condition <- factor(seu_ef$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu_ef$sample <- factor(seu_ef$sample,
                        levels = intersect(names(SAMPLE_COLORS), unique(seu_ef$sample)))

# QC filter
keep_qc_ef <- seu_ef$total_counts >= MIN_TOTAL_COUNTS &
              seu_ef$n_genes_by_counts >= MIN_GENES_DETECTED
seu_ef <- subset(seu_ef, cells = colnames(seu_ef)[keep_qc_ef])

# Keep only cell types present in data
available_ef <- intersect(EPI_FIB_CELLTYPES, unique(seu_ef$Tier1_celltype))
seu_ef <- subset(seu_ef, subset = Tier1_celltype %in% available_ef)
seu_ef$Tier1_celltype <- droplevels(factor(seu_ef$Tier1_celltype, levels = available_ef))
cat(sprintf("Epithelial & fibroblast analysis: %d cells, %d cell types\n", ncol(seu_ef), length(available_ef)))

# Build count table
count_df_ef <- seu_ef@meta.data %>%
  group_by(sample, condition, Tier1_celltype, .drop = TRUE) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Tier1_celltype, values_from = n, values_fill = 0) %>%
  as.data.frame()

ct_cols_ef <- available_ef

cat("\nCell counts per sample (epithelial & fibroblasts):\n")
print(count_df_ef)

# Proportions with pseudocount
count_df_ef$total <- rowSums(count_df_ef[, ct_cols_ef])
prop_df_ef <- count_df_ef
for (ct in ct_cols_ef) {
  prop_df_ef[[ct]] <- count_df_ef[[ct]] + (count_df_ef$total / sum(count_df_ef$total))
}
prop_df_ef$total_adj <- rowSums(prop_df_ef[, ct_cols_ef])
for (ct in ct_cols_ef) {
  prop_df_ef[[ct]] <- prop_df_ef[[ct]] / prop_df_ef$total_adj
}

# Dirichlet regression
dr_data_ef <- as.data.frame(prop_df_ef)
dr_data_ef$condition <- factor(dr_data_ef$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
AL_ef <- DR_data(dr_data_ef[, ct_cols_ef])

fit_ef <- DirichReg(AL_ef ~ condition, data = dr_data_ef)
fit_null_ef <- DirichReg(AL_ef ~ 1, data = dr_data_ef)

cat("\nLikelihood ratio test (condition effect, epithelial & fibroblasts):\n")
lr_test_ef <- anova(fit_null_ef, fit_ef)
print(lr_test_ef)

# Extract results
u_ef <- summary(fit_ef)
coef_mat_ef <- u_ef$coef.mat
pvals_rows_ef <- grep("Intercept", rownames(coef_mat_ef), invert = TRUE)
pvals_ef <- coef_mat_ef[pvals_rows_ef, 4]
estimates_ef <- coef_mat_ef[pvals_rows_ef, 1]

results_df_ef <- data.frame(
  celltype = u_ef$varnames,
  estimate = estimates_ef,
  pval = pvals_ef,
  stringsAsFactors = FALSE
)
results_df_ef$p.adj <- p.adjust(results_df_ef$pval, method = "fdr")
results_df_ef$direction <- ifelse(results_df_ef$estimate > 0, "Up in DKO", "Up in WT")
results_df_ef$significant <- results_df_ef$p.adj < FDR_THRESHOLD
results_df_ef <- results_df_ef[order(results_df_ef$pval), ]

cat("\nResults (epithelial & fibroblasts):\n")
print(results_df_ef)

# Save
write.csv(results_df_ef, file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_epi_fib.csv"), row.names = FALSE)
write.csv(count_df_ef, file.path(OUTPUT_DIR, "epi_fib_counts_per_sample.csv"), row.names = FALSE)

prop_out_ef <- count_df_ef %>%
  pivot_longer(cols = all_of(ct_cols_ef), names_to = "celltype", values_to = "n") %>%
  group_by(sample) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
write.csv(prop_out_ef, file.path(OUTPUT_DIR, "epi_fib_proportions_per_sample.csv"), row.names = FALSE)

# --- Plots for epithelial & fibroblasts ---

# Stacked bar per sample
obs_prop_ef <- count_df_ef %>%
  pivot_longer(cols = all_of(ct_cols_ef), names_to = "celltype", values_to = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = available_ef))

obs_prop_ef <- obs_prop_ef %>%
  left_join(distinct(count_df_ef[, c("sample", "condition")]), by = c("sample", "condition"))

p_stack_ef <- ggplot(obs_prop_ef, aes(x = sample, y = prop, fill = celltype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = TIER1_COLORS, name = "Cell type") +
  labs(title = "Epithelial & fibroblast composition per sample",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "epi_fib_composition_stacked_bar.pdf"), p_stack_ef, width = 10, height = 6)

# Collapsed bar
collapsed_prop_ef <- count_df_ef %>%
  group_by(condition) %>%
  summarise(across(all_of(ct_cols_ef), sum), .groups = "drop") %>%
  pivot_longer(cols = all_of(ct_cols_ef), names_to = "celltype", values_to = "n") %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(celltype = factor(celltype, levels = available_ef))

p_collapsed_ef <- ggplot(collapsed_prop_ef, aes(x = condition, y = prop, fill = celltype)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = TIER1_COLORS, name = "Cell type") +
  labs(title = "Epithelial & fibroblast composition: DKO vs WT (pooled)",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "epi_fib_composition_collapsed.pdf"), p_collapsed_ef, width = 7, height = 6)

# Boxplot with significance
sig_labels_ef <- results_df_ef %>%
  mutate(
    label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      p.adj < 0.1   ~ "^",
      TRUE          ~ ""
    ),
    celltype = factor(celltype, levels = available_ef)
  )

y_positions_ef <- obs_prop_ef %>%
  group_by(celltype) %>%
  summarise(y_max = max(prop) * 1.05, .groups = "drop")

sig_labels_ef <- sig_labels_ef %>%
  left_join(y_positions_ef, by = "celltype") %>%
  filter(label != "")

p_box_ef <- ggplot(obs_prop_ef, aes(x = celltype, y = prop, fill = condition)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  geom_point(aes(group = condition), position = position_dodge(width = 0.75),
             size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Epithelial & fibroblast proportions: DKO vs WT",
       subtitle = "Dirichlet regression adj. p-value: *** <0.001, ** <0.01, * <0.05, ^ <0.1",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")

if (nrow(sig_labels_ef) > 0) {
  p_box_ef <- p_box_ef +
    geom_text(data = sig_labels_ef, aes(x = celltype, y = y_max, label = label),
              inherit.aes = FALSE, size = 5, vjust = 0)
}

ggsave(file.path(OUTPUT_DIR, "epi_fib_proportions_boxplot.pdf"), p_box_ef, width = 10, height = 5)

# Dot plot
results_df_ef$celltype <- factor(results_df_ef$celltype, levels = rev(available_ef))

p_dot_ef <- ggplot(results_df_ef, aes(x = estimate, y = celltype,
                                       size = -log10(p.adj),
                                       color = estimate)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient2(low = "#4575B4", mid = "grey70", high = "#D73027",
                        midpoint = 0, name = "Estimate") +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
  labs(title = "Dirichlet regression: DKO vs WT (Epithelial & Fibroblasts)",
       subtitle = "Positive estimate = enriched in DKO",
       x = "Estimate (DKO vs WT)", y = NULL) +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_epi_fib_dotplot.pdf"), p_dot_ef, width = 7, height = 5)

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("Dirichlet regression (DKO vs WT) complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("\n--- Immune only ---\n"))
cat(sprintf("Samples: %d %s, %d %s\n",
            sum(count_df$condition == REFERENCE_CONDITION), REFERENCE_CONDITION,
            sum(count_df$condition == TEST_CONDITION), TEST_CONDITION))
cat(sprintf("Significant cell types (FDR < %.2f): %d / %d\n", FDR_THRESHOLD,
            sum(results_df$significant), nrow(results_df)))
cat(sprintf("\n--- Epithelial & Fibroblasts ---\n"))
cat(sprintf("Significant cell types (FDR < %.2f): %d / %d\n", FDR_THRESHOLD,
            sum(results_df_ef$significant), nrow(results_df_ef)))
cat("\nFiles generated:\n")
cat("  Immune:\n")
cat("    - dirichlet_DKOvsWT_immune.csv\n")
cat("    - immune_counts_per_sample.csv\n")
cat("    - immune_proportions_per_sample.csv\n")
cat("    - immune_composition_stacked_bar.pdf\n")
cat("    - immune_composition_collapsed.pdf\n")
cat("    - immune_proportions_boxplot.pdf\n")
cat("    - dirichlet_DKOvsWT_dotplot.pdf\n")
cat("  Epithelial & Fibroblasts:\n")
cat("    - dirichlet_DKOvsWT_epi_fib.csv\n")
cat("    - epi_fib_counts_per_sample.csv\n")
cat("    - epi_fib_proportions_per_sample.csv\n")
cat("    - epi_fib_composition_stacked_bar.pdf\n")
cat("    - epi_fib_composition_collapsed.pdf\n")
cat("    - epi_fib_proportions_boxplot.pdf\n")
cat("    - dirichlet_DKOvsWT_epi_fib_dotplot.pdf\n")
