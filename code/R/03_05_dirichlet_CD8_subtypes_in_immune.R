# 03_05_dirichlet_CD8_subtypes_in_immune.R
# Dirichlet regression (DKO vs WT) on immune composition with CD8+ T cells split
# into their 03_01 subtypes. Same model as the immune analysis in
# 11_01_spatial_dirichlet.R (IMMUNE_CELLTYPES, QC filter, pseudocount,
# AL ~ condition), except that "CD8+ T" is replaced by CD8T_C0 ... CD8T_Cn, so
# each subtype is tested as a fraction of all immune cells.
#
# Read-only on existing RDS files.

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DirichletReg)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

TCELL_DIR  <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters")
OUTPUT_DIR <- file.path(TCELL_DIR, "dirichlet_immune")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Parameters
# =============================================================================
TCELL_ASSIGN_CSV  <- file.path(TCELL_DIR, "tcell_cluster_assignments.csv")  # from 03_01
TCELL_CLUSTER_COL <- "tcell_clust_res.0.4"
TCELL_TYPE        <- "CD8+ T"

# =============================================================================
# 1. Load metadata (QC-filtered immune cells) and CD8+ T subtypes
# =============================================================================
meta <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))@meta.data
meta$condition <- gsub("^KO$", TEST_CONDITION, meta$condition)
meta$condition <- factor(meta$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
meta$sample <- factor(meta$sample, levels = intersect(names(SAMPLE_COLORS), unique(meta$sample)))
meta <- meta[meta$total_counts >= MIN_TOTAL_COUNTS &
             meta$n_genes_by_counts >= MIN_GENES_DETECTED &
             meta$Tier1_celltype %in% IMMUNE_CELLTYPES, ]
cat(sprintf("QC-filtered immune cells: %d\n", nrow(meta)))

tc <- read.csv(TCELL_ASSIGN_CSV, stringsAsFactors = FALSE, check.names = FALSE)
subtype_levels <- paste0("CD8T_C", sort(unique(tc[[TCELL_CLUSTER_COL]])))
cd8_cells <- rownames(meta)[meta$Tier1_celltype == TCELL_TYPE]
stopifnot(setequal(cd8_cells, tc$cell))   # 03_01 used the same QC filter

meta$component <- as.character(meta$Tier1_celltype)
meta[tc$cell, "component"] <- paste0("CD8T_C", tc[[TCELL_CLUSTER_COL]])

comp_cols <- unlist(lapply(IMMUNE_CELLTYPES, function(ct)
  if (ct == TCELL_TYPE) subtype_levels else ct))
meta$component <- factor(meta$component, levels = comp_cols)

COMP_COLORS <- c(TIER1_COLORS[setdiff(IMMUNE_CELLTYPES, TCELL_TYPE)],
                 setNames(tcell_subtype_colors(subtype_levels), subtype_levels))

# =============================================================================
# 2. Count table, pseudocount proportions (as in 11_01)
# =============================================================================
count_df <- meta %>%
  group_by(sample, condition, component, .drop = TRUE) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = component, values_from = n, values_fill = 0) %>%
  as.data.frame()
count_df <- count_df[, c("sample", "condition", comp_cols)]
cat("\nCell counts per sample:\n"); print(count_df)

count_df$total <- rowSums(count_df[, comp_cols])
prop_df <- count_df
for (ct in comp_cols) {
  prop_df[[ct]] <- count_df[[ct]] + (count_df$total / sum(count_df$total))
}
prop_df$total_adj <- rowSums(prop_df[, comp_cols])
for (ct in comp_cols) {
  prop_df[[ct]] <- prop_df[[ct]] / prop_df$total_adj
}

# =============================================================================
# 3. Dirichlet regression: DKO vs WT
# =============================================================================
dr_data <- as.data.frame(prop_df)
# DR_data needs syntactic column names; map back afterwards
safe_names <- make.names(comp_cols)
AL <- DR_data(setNames(dr_data[, comp_cols], safe_names))

fit <- DirichReg(AL ~ condition, data = dr_data)
fit_null <- DirichReg(AL ~ 1, data = dr_data)
lr_test <- anova(fit_null, fit)
cat("\nLikelihood ratio test (condition effect):\n"); print(lr_test)

u <- summary(fit)
coef_mat <- u$coef.mat
rows <- grep("Intercept", rownames(coef_mat), invert = TRUE)
results_df <- data.frame(
  component = comp_cols[match(u$varnames, safe_names)],
  estimate = coef_mat[rows, 1],
  pval = coef_mat[rows, 4],
  stringsAsFactors = FALSE
)
results_df$p.adj <- p.adjust(results_df$pval, method = "fdr")
results_df$direction <- ifelse(results_df$estimate > 0, "Up in DKO", "Up in WT")
results_df$significant <- results_df$p.adj < FDR_THRESHOLD
results_df <- results_df[order(results_df$pval), ]
cat("\nResults (positive estimate = higher proportion in DKO):\n"); print(results_df, digits = 3)

write.csv(results_df, file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_immune_CD8subtypes.csv"), row.names = FALSE)
write.csv(count_df, file.path(OUTPUT_DIR, "immune_CD8subtypes_counts_per_sample.csv"), row.names = FALSE)

# =============================================================================
# 4. Plots
# =============================================================================
obs_prop <- count_df %>%
  pivot_longer(all_of(comp_cols), names_to = "component", values_to = "n") %>%
  group_by(sample) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  mutate(component = factor(component, levels = comp_cols))
write.csv(obs_prop, file.path(OUTPUT_DIR, "immune_CD8subtypes_proportions_per_sample.csv"), row.names = FALSE)

sig_labels <- results_df %>%
  mutate(label = case_when(p.adj < 0.001 ~ "***", p.adj < 0.01 ~ "**",
                           p.adj < 0.05 ~ "*", p.adj < FDR_THRESHOLD ~ "^", TRUE ~ ""),
         component = factor(component, levels = comp_cols)) %>%
  left_join(obs_prop %>% group_by(component) %>% summarise(y_max = max(prop) * 1.05),
            by = "component") %>%
  filter(label != "")

p_box <- ggplot(obs_prop, aes(x = component, y = prop, fill = condition)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.8) +
  geom_point(aes(group = condition), position = position_dodge(width = 0.75),
             size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Immune composition with CD8+ T subtypes: DKO vs WT",
       subtitle = "Dirichlet regression adj. p-value: *** <0.001, ** <0.01, * <0.05, ^ <0.1",
       x = NULL, y = "Proportion of immune cells") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
if (nrow(sig_labels) > 0) {
  p_box <- p_box + geom_text(data = sig_labels, aes(x = component, y = y_max, label = label),
                             inherit.aes = FALSE, size = 5, vjust = 0)
}
ggsave(file.path(OUTPUT_DIR, "immune_CD8subtypes_proportions_boxplot.pdf"), p_box, width = 11, height = 5)

p_stack <- ggplot(obs_prop, aes(x = sample, y = prop, fill = component)) +
  geom_col() +
  scale_fill_manual(values = COMP_COLORS, name = NULL) +
  labs(title = "Immune composition per sample (CD8+ T split by subtype)",
       x = NULL, y = "Proportion") +
  theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "immune_CD8subtypes_stacked_bar.pdf"), p_stack, width = 8, height = 6)

results_df$component <- factor(results_df$component, levels = rev(comp_cols))
p_dot <- ggplot(results_df, aes(x = estimate, y = component, size = -log10(p.adj), color = estimate)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_gradient2(low = "#4575B4", mid = "grey70", high = "#D73027", midpoint = 0,
                        name = "Estimate") +
  scale_size_continuous(name = "-log10(FDR)", range = c(2, 8)) +
  labs(title = "Dirichlet regression: DKO vs WT (immune, CD8+ T subtypes)",
       subtitle = "Positive estimate = enriched in DKO", x = "Estimate (DKO vs WT)", y = NULL) +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "dirichlet_DKOvsWT_immune_CD8subtypes_dotplot.pdf"), p_dot, width = 7, height = 6)

cat(sprintf("\nGlobal LRT p-value (condition): %.3g\n", lr_test$`Pr(>Chi)`[2]))
cat("Output directory:", OUTPUT_DIR, "\n")
