# 03_02_dirichlet_test.R
# Dirichlet multinomial regression for differential composition of CD8+ T cell
# subtypes (clusters from 03_01) between DKO and WT. Same model as
# 11_01_spatial_dirichlet.R: one observation per sample, subtype proportions
# within CD8+ T cells, ~ condition.
#
# Read-only on existing RDS files.

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)
library(DirichletReg)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

TCELL_DIR  <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters")
OUTPUT_DIR <- file.path(TCELL_DIR, "dirichlet")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Parameters
# =============================================================================
TCELL_RDS              <- file.path(TCELL_DIR, "DKO_CD8T_clusters.rds")  # from 03_01
TCELL_CLUSTER_COL      <- "tcell_clust_res.0.4"
TCELL_EXCLUDE_CLUSTERS <- character(0)   # e.g. "C0" to drop a contaminated cluster

# =============================================================================
# 1. Load CD8+ T object (already QC-filtered in 03_01)
# =============================================================================
seu <- readRDS(TCELL_RDS)
cat("Loaded:", ncol(seu), "CD8+ T cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample, levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

clusters <- seu@meta.data[[TCELL_CLUSTER_COL]]
SUBTYPE_ORDER <- paste0("C", levels(clusters))
seu$subtype <- factor(paste0("C", clusters), levels = SUBTYPE_ORDER)
CLUSTER_COLORS <- tcell_subtype_colors(SUBTYPE_ORDER)

if (length(TCELL_EXCLUDE_CLUSTERS) > 0) {
  seu <- subset(seu, cells = colnames(seu)[!seu$subtype %in% TCELL_EXCLUDE_CLUSTERS])
  cat("Excluded subtypes:", paste(TCELL_EXCLUDE_CLUSTERS, collapse = ", "), "\n")
}
seu$subtype <- droplevels(seu$subtype)
subtype_cols <- levels(seu$subtype)
cat(sprintf("Subtype analysis: %d cells, %d subtypes (%s)\n",
            ncol(seu), length(subtype_cols), TCELL_CLUSTER_COL))

# =============================================================================
# 2. Build count table: sample x subtype
# =============================================================================
count_df <- seu@meta.data %>%
  group_by(sample, condition, subtype, .drop = TRUE) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = subtype, values_from = n, values_fill = 0) %>%
  as.data.frame()

cat("\nCell counts per sample:\n")
print(count_df)

# =============================================================================
# 3. Proportions with pseudocount
# =============================================================================
count_df$total <- rowSums(count_df[, subtype_cols])

# Add small pseudocount (proportional to sample size) to avoid zeros
prop_df <- count_df
for (st in subtype_cols) {
  prop_df[[st]] <- count_df[[st]] + (count_df$total / sum(count_df$total))
}
prop_df$total_adj <- rowSums(prop_df[, subtype_cols])
for (st in subtype_cols) {
  prop_df[[st]] <- prop_df[[st]] / prop_df$total_adj
}

# =============================================================================
# 4. Dirichlet regression: DKO vs WT
# =============================================================================
dr_data <- as.data.frame(prop_df)
dr_data$condition <- factor(dr_data$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))

AL <- DR_data(dr_data[, subtype_cols])

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
  subtype = u$varnames,
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
fitted_vals <- as.data.frame(fit$fitted.values$mu)
fitted_vals$condition <- dr_data$condition
fitted_summary <- fitted_vals %>%
  group_by(condition) %>%
  summarise(across(all_of(subtype_cols), mean), .groups = "drop")

cat("\nFitted mean proportions:\n")
print(as.data.frame(fitted_summary))

# Save results
write.csv(results_df, file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_CD8T_subtypes.csv"), row.names = FALSE)
write.csv(count_df, file.path(OUTPUT_DIR, "CD8T_subtype_counts_per_sample.csv"), row.names = FALSE)
write.csv(data.frame(model = c("null", "condition"),
                     loglik = c(fit_null$logLik, fit$logLik),
                     df = c(fit_null$npar, fit$npar),
                     LR_pval = c(NA, lr_test$`Pr(>Chi)`[2])),
          file.path(OUTPUT_DIR, "dirichlet_LRT_condition.csv"), row.names = FALSE)

obs_prop <- count_df %>%
  pivot_longer(cols = all_of(subtype_cols), names_to = "subtype", values_to = "n") %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(subtype = factor(subtype, levels = subtype_cols))
write.csv(obs_prop, file.path(OUTPUT_DIR, "CD8T_subtype_proportions_per_sample.csv"), row.names = FALSE)

# =============================================================================
# 5. Visualization: stacked bars
# =============================================================================
cat("\nGenerating plots...\n")

p_stack <- ggplot(obs_prop, aes(x = sample, y = prop, fill = subtype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CLUSTER_COLORS, name = "Subtype") +
  labs(title = "CD8+ T subtype composition per sample",
       x = NULL, y = "Proportion of CD8+ T cells") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "CD8T_subtype_composition_stacked_bar.pdf"), p_stack, width = 7, height = 5)

collapsed_prop <- count_df %>%
  group_by(condition) %>%
  summarise(across(all_of(subtype_cols), sum), .groups = "drop") %>%
  pivot_longer(cols = all_of(subtype_cols), names_to = "subtype", values_to = "n") %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(subtype = factor(subtype, levels = subtype_cols))

p_collapsed <- ggplot(collapsed_prop, aes(x = condition, y = prop, fill = subtype)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = CLUSTER_COLORS, name = "Subtype") +
  labs(title = "CD8+ T subtype composition: DKO vs WT (pooled)",
       x = NULL, y = "Proportion") +
  theme_publication() +
  theme(legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "CD8T_subtype_composition_collapsed.pdf"), p_collapsed, width = 5, height = 5)

# =============================================================================
# 6. Visualization: boxplot of proportions by condition (with significance)
# =============================================================================
sig_labels <- results_df %>%
  mutate(
    label = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    subtype = factor(subtype, levels = subtype_cols)
  )

y_positions <- obs_prop %>%
  group_by(subtype) %>%
  summarise(y_max = max(prop) * 1.05, .groups = "drop")

sig_labels <- sig_labels %>%
  left_join(y_positions, by = "subtype") %>%
  filter(label != "")

p_box <- ggplot(obs_prop, aes(x = subtype, y = prop, fill = condition)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  geom_point(aes(group = condition), position = position_dodge(width = 0.75),
             size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "CD8+ T subtype proportions: DKO vs WT",
       subtitle = "\nDirichlet regression adj. p-value: *** <0.001, ** <0.01, * <0.05",
       x = NULL, y = "Proportion of CD8+ T cells") +
  theme_publication() +
  theme(legend.position = "top")

if (nrow(sig_labels) > 0) {
  p_box <- p_box +
    geom_text(data = sig_labels, aes(x = subtype, y = y_max, label = label),
              inherit.aes = FALSE, size = 5, vjust = 0)
}

ggsave(file.path(OUTPUT_DIR, "CD8T_subtype_proportions_boxplot.pdf"), p_box, width = 6, height = 5)

# =============================================================================
# 7. Visualization: dot plot of Dirichlet results
# =============================================================================
results_df$subtype <- factor(results_df$subtype, levels = rev(subtype_cols))

p_dot <- ggplot(results_df, aes(x = estimate, y = subtype,
                                size = -log10(p.adj),
                                color = estimate)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient2(low = "#4575B4", mid = "grey70", high = "#D73027",
                        midpoint = 0, name = "Estimate") +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
  labs(title = "Dirichlet regression: DKO vs WT (CD8+ T subtypes)",
       subtitle = "Positive estimate = enriched in DKO",
       x = "Estimate (DKO vs WT)", y = NULL) +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_CD8T_subtypes_dotplot.pdf"), p_dot, width = 7, height = 4)

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("Dirichlet regression (DKO vs WT, CD8+ T subtypes) complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("Samples: %d %s, %d %s\n",
            sum(count_df$condition == REFERENCE_CONDITION), REFERENCE_CONDITION,
            sum(count_df$condition == TEST_CONDITION), TEST_CONDITION))
cat(sprintf("Global LRT p-value (condition): %.3g\n", lr_test$`Pr(>Chi)`[2]))
cat(sprintf("Significant subtypes (FDR < %.2f): %d / %d\n", FDR_THRESHOLD,
            sum(results_df$significant), nrow(results_df)))
