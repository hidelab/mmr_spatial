# 09_proportion_analysis.R
# Cell and domain proportion analysis across conditions and samples.
# Uses merged Banksy spatial domains from 02_03_merge_domains_08lambda.R.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# Parameters
# =============================================================================
# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "09_proportion_analysis")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

DOMAIN_COL <- COL_DOMAIN

# =============================================================================
# 1. Load metadata and filter low-quality cells
# =============================================================================
meta <- readRDS(file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08",
                          "merge_domains", "DKO_banksy_merged_metadata.rds"))
cat("Loaded metadata:", nrow(meta), "cells\n")

# Recode condition
meta$condition <- gsub("^KO$", TEST_CONDITION, meta$condition)

# Filter low-quality cells
n_before <- nrow(meta)
meta <- meta[meta$total_counts >= MIN_TOTAL_COUNTS &
             meta$n_genes_by_counts >= MIN_GENES_DETECTED, ]
n_after <- nrow(meta)
cat(sprintf("Cell filtering: %d -> %d cells (removed %d, %.1f%%)\n",
            n_before, n_after, n_before - n_after,
            100 * (n_before - n_after) / n_before))
cat(sprintf("  Filters: total_counts >= %d, n_genes_by_counts >= %d\n",
            MIN_TOTAL_COUNTS, MIN_GENES_DETECTED))

# Assign domain names
meta$domain_name <- DOMAIN_NAMES[as.character(meta[[DOMAIN_COL]])]
meta$domain_name <- factor(meta$domain_name, levels = DOMAIN_ORDER)

cat("Conditions:", paste(unique(meta$condition), collapse = ", "), "\n")
cat("Samples:", paste(sort(unique(meta$sample)), collapse = ", "), "\n")
cat("Domains:", paste(levels(meta$domain_name), collapse = ", "), "\n")
cat("Cell types:", length(unique(meta$Tier1_celltype)), "\n\n")

# =============================================================================
# 2. Domain proportions per sample
# =============================================================================
cat("== Domain proportions per sample ==\n")

prop_domain_sample <- as.data.frame(meta) %>%
  count(sample, condition, domain_name) %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Stacked bar: domain proportions per sample
p1 <- ggplot(prop_domain_sample, aes(x = sample, y = prop, fill = domain_name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = DOMAIN_COLORS) +
  labs(title = "Domain proportions per sample",
       x = NULL, y = "Proportion", fill = "Domain") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "domain_proportions_per_sample.pdf"), p1, width = 8, height = 5)

# =============================================================================
# 3. Domain proportions by condition (aggregated)
# =============================================================================
cat("== Domain proportions by condition ==\n")

prop_domain_cond <- as.data.frame(meta) %>%
  count(condition, domain_name) %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p2 <- ggplot(prop_domain_cond, aes(x = condition, y = prop, fill = domain_name)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = DOMAIN_COLORS) +
  labs(title = "Domain proportions by condition",
       x = NULL, y = "Proportion", fill = "Domain") +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "domain_proportions_by_condition.pdf"), p2, width = 5, height = 5)

# t-test on per-sample proportions (samples are the biological replicates, not cells)
ttest_domain <- prop_domain_sample %>%
  group_by(domain_name) %>%
  summarise(
    p_value = tryCatch(
      t.test(prop ~ condition)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("t-test p = %.3f", p_value))

cat("t-test on per-sample domain proportions (WT vs DKO):\n")
print(as.data.frame(ttest_domain))

# Per-sample dot plot for statistical comparison
p2b <- ggplot(prop_domain_sample, aes(x = condition, y = prop, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.3) +
  geom_jitter(width = 0.1, size = 2.5) +
  geom_text(data = ttest_domain, aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, size = 2.8, vjust = 1.5) +
  scale_color_manual(values = CONDITION_COLORS) +
  facet_wrap(~ domain_name, scales = "free_y") +
  labs(title = "Domain proportions: WT vs DKO (per sample)",
       x = NULL, y = "Proportion of cells in domain") +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "domain_proportions_by_condition_boxplot.pdf"), p2b, width = 8, height = 6)

# =============================================================================
# 4. Cell type proportions per domain
# =============================================================================
cat("== Cell type proportions per domain ==\n")

prop_ct_domain <- as.data.frame(meta) %>%
  count(domain_name, Tier1_celltype) %>%
  group_by(domain_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Stacked bar: cell types within each domain
p3 <- ggplot(prop_ct_domain, aes(x = domain_name, y = prop, fill = Tier1_celltype)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = TIER1_COLORS, breaks = TIER1_ORDER) +
  labs(title = "Cell type composition per domain",
       x = NULL, y = "Proportion", fill = "Cell type") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.text = element_text(size = 7))

ggsave(file.path(OUTPUT_DIR, "celltype_composition_per_domain.pdf"), p3, width = 10, height = 6)

# =============================================================================
# 5. Cell type proportions per domain, split by condition
# =============================================================================
cat("== Cell type proportions per domain by condition ==\n")

prop_ct_domain_cond <- as.data.frame(meta) %>%
  count(domain_name, condition, Tier1_celltype) %>%
  group_by(domain_name, condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p4 <- ggplot(prop_ct_domain_cond, aes(x = condition, y = prop, fill = Tier1_celltype)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = TIER1_COLORS, breaks = TIER1_ORDER) +
  facet_wrap(~ domain_name, nrow = 1) +
  labs(title = "Cell type composition per domain, split by condition",
       x = NULL, y = "Proportion", fill = "Cell type") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 7))

ggsave(file.path(OUTPUT_DIR, "celltype_composition_per_domain_by_condition.pdf"), p4, width = 12, height = 5)

# =============================================================================
# 6. Cell type proportions across domains per sample (heatmap-style)
# =============================================================================
cat("== Cell type x domain x sample proportions ==\n")

prop_ct_sample_domain <- as.data.frame(meta) %>%
  count(sample, condition, domain_name, Tier1_celltype) %>%
  group_by(sample, domain_name) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Focus on major immune/tumor cell types for readability
p5 <- ggplot(prop_ct_sample_domain,
             aes(x = sample, y = Tier1_celltype, fill = prop)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~ domain_name, nrow = 1) +
  scale_fill_viridis_c(option = "plasma", limits = c(0, NA)) +
  labs(title = "Cell type proportions per sample within each domain",
       x = NULL, y = NULL, fill = "Proportion") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        strip.text = element_text(size = 9, face = "bold"))

ggsave(file.path(OUTPUT_DIR, "celltype_proportions_heatmap_sample_domain.pdf"), p5, width = 16, height = 7)

# =============================================================================
# 7. Domain size (absolute cell counts) per sample
# =============================================================================
cat("== Absolute cell counts per domain ==\n")

counts_domain_sample <- as.data.frame(meta) %>%
  count(sample, condition, domain_name)

p6 <- ggplot(counts_domain_sample, aes(x = sample, y = n, fill = domain_name)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = DOMAIN_COLORS) +
  labs(title = "Cell counts per domain per sample",
       x = NULL, y = "Number of cells", fill = "Domain") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "cell_counts_per_domain_sample.pdf"), p6, width = 9, height = 5)

# =============================================================================
# 8. Immune Engulfing / Tumor Core ratio per sample
# =============================================================================
cat("== Immune Engulfing to Tumor Core ratio ==\n")

ratio_df <- as.data.frame(meta) %>%
  count(sample, condition, domain_name) %>%
  filter(domain_name %in% c("Immune Engulfing", "Tumor Core")) %>%
  pivot_wider(names_from = domain_name, values_from = n, values_fill = 0) %>%
  mutate(ratio = `Immune Engulfing` / `Tumor Core`,
         log2_ratio = log2(ratio))

cat("Immune Engulfing / Tumor Core ratio per sample:\n")
print(ratio_df %>% select(sample, condition, `Immune Engulfing`, `Tumor Core`, ratio, log2_ratio))

t_p <- t.test(log2_ratio ~ condition, data = ratio_df)$p.value
cat(sprintf("t-test on log2(ratio) p-value: %.4f\n", t_p))

p7 <- ggplot(ratio_df, aes(x = condition, y = log2_ratio, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.3) +
  geom_jitter(width = 0.1, size = 3) +
  geom_text(aes(label = sample), vjust = -1, size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = CONDITION_COLORS) +
  annotate("text", x = 1.5, y = max(ratio_df$log2_ratio) * 1.1,
           label = sprintf("t-test p = %.3f", t_p), size = 3.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Immune Engulfing / Tumor Core ratio",
       subtitle = "Log2-transformed per-sample cell count ratio between domains",
       x = NULL, y = "log2(Immune Engulfing / Tumor Core)") +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "immune_engulfing_tumor_core_ratio.pdf"), p7, width = 5, height = 5)

write.csv(ratio_df, file.path(OUTPUT_DIR, "immune_engulfing_tumor_core_ratio.csv"), row.names = FALSE)

# =============================================================================
# 9. Immune cell ratio: Stroma (domain 4) vs Immune Engulfing (domain 3)
# =============================================================================
cat("== Immune cell ratio: Stroma vs Immune Engulfing ==\n")

# IMMUNE_CELLTYPES defined in 00_parameters.R


# Count immune cells of interest per sample per domain (only domains 3 & 4)
immune_counts <- as.data.frame(meta) %>%
  filter(Tier1_celltype %in% IMMUNE_CELLTYPES,
         domain_name %in% c("Immune Engulfing", "Stroma")) %>%
  count(sample, condition, domain_name, Tier1_celltype)

# Compute ratio: Stroma / Immune Engulfing per cell type per sample
immune_ratio <- immune_counts %>%
  pivot_wider(names_from = domain_name, values_from = n, values_fill = 0) %>%
  mutate(ratio = (Stroma + 1) / (`Immune Engulfing` + 1),
         log2_ratio = log2(ratio))

cat("Immune cell counts and log2 ratios (Stroma / Immune Engulfing):\n")
print(immune_ratio)

# t-test on log2 ratio per cell type
ttest_results <- immune_ratio %>%
  group_by(Tier1_celltype) %>%
  summarise(
    p_value = tryCatch(
      t.test(log2_ratio ~ condition)$p.value,
      error = function(e) NA_real_
    ),
    mean_WT = mean(log2_ratio[condition == REFERENCE_CONDITION], na.rm = TRUE),
    mean_DKO = mean(log2_ratio[condition == TEST_CONDITION], na.rm = TRUE),
    .groups = "drop"
  )

cat("\nt-test results on log2(Stroma/Immune Engulfing), WT vs DKO:\n")
print(as.data.frame(ttest_results))

# Plot: faceted by cell type
p_label <- immune_ratio %>%
  group_by(Tier1_celltype) %>%
  summarise(y_pos = max(log2_ratio) * 1.15, .groups = "drop") %>%
  left_join(ttest_results %>% select(Tier1_celltype, p_value), by = "Tier1_celltype") %>%
  mutate(label = sprintf("p = %.3f", p_value))

p8 <- ggplot(immune_ratio,
             aes(x = condition, y = log2_ratio, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.3) +
  geom_jitter(width = 0.1, size = 2.5) +
  geom_text(data = p_label, aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = CONDITION_COLORS) +
  facet_wrap(~ Tier1_celltype, scales = "free_y", nrow = 2) +
  labs(title = "Immune cell ratio: Stroma / Immune Engulfing (log2)",
       subtitle = "Per-sample log2 ratio of immune cells in domain 4 vs domain 3",
       x = NULL, y = "log2(Stroma / Immune Engulfing)") +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "immune_cell_ratio_stroma_vs_immune_engulfing.pdf"), p8, width = 9, height = 6)

write.csv(immune_ratio, file.path(OUTPUT_DIR, "immune_cell_ratio_stroma_vs_immune_engulfing.csv"), row.names = FALSE)
write.csv(ttest_results, file.path(OUTPUT_DIR, "immune_cell_ratio_ttest_results.csv"), row.names = FALSE)

# =============================================================================
# 10. Immune cell ratio: Immune Engulfing / Tumor Core per cell type
# =============================================================================
cat("== Immune cell ratio: Immune Engulfing / Tumor Core ==\n")

# Count immune cells per sample, domain (Immune Engulfing & Tumor Core)
immune_counts <- as.data.frame(meta) %>%
  filter(Tier1_celltype %in% IMMUNE_CELLTYPES,
         domain_name %in% c("Immune Engulfing", "Tumor Core")) %>%
  count(sample, condition, domain_name, Tier1_celltype) %>%
  pivot_wider(names_from = domain_name, values_from = n, values_fill = 0)

# Log2 ratio: log2((Immune Engulfing + 1) / (Tumor Core + 1))
immune_counts <- immune_counts %>%
  mutate(ratio = (`Immune Engulfing` + 1) / (`Tumor Core` + 1),
         log2_ratio = log2(ratio))

cat("Immune Engulfing / Tumor Core log2 ratio per cell type per sample:\n")
print(as.data.frame(immune_counts))

# t-test on log2 ratio per cell type
ttest_ratio <- immune_counts %>%
  group_by(Tier1_celltype) %>%
  summarise(
    p_value = tryCatch(
      t.test(log2_ratio ~ condition)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  )

cat("\nt-test on log2(Immune Engulfing / Tumor Core), WT vs DKO:\n")
print(as.data.frame(ttest_ratio))

# Plot labels
p_ratio_label <- immune_counts %>%
  group_by(Tier1_celltype) %>%
  summarise(y_pos = max(log2_ratio) * 1.1, .groups = "drop") %>%
  left_join(ttest_ratio, by = "Tier1_celltype") %>%
  mutate(label = sprintf("p = %.3f", p_value))

p9 <- ggplot(immune_counts, aes(x = condition, y = log2_ratio, color = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.4, alpha = 0.3) +
  geom_jitter(width = 0.1, size = 2.5) +
  geom_text(data = p_ratio_label, aes(x = 1.5, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = CONDITION_COLORS) +
  facet_wrap(~ Tier1_celltype, scales = "free_y", nrow = 2) +
  labs(title = "Immune cell ratio: Immune Engulfing / Tumor Core (log2)",
       subtitle = "Per-sample log2 ratio (pseudocount +1 added to each domain)",
       x = NULL, y = "log2(Immune Engulfing / Tumor Core)") +
  theme_publication() +
  theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "immune_cell_ratio_engulfing_vs_tumor_core.pdf"), p9, width = 9, height = 6)

write.csv(immune_counts, file.path(OUTPUT_DIR, "immune_cell_ratio_engulfing_vs_tumor_core.csv"), row.names = FALSE)
write.csv(ttest_ratio, file.path(OUTPUT_DIR, "immune_cell_ratio_engulfing_vs_tumor_core_ttest.csv"), row.names = FALSE)

# =============================================================================
# 11. Summary table: proportions
# =============================================================================
cat("== Summary tables ==\n")

# Domain proportions per sample
write.csv(prop_domain_sample %>% select(-n) %>% pivot_wider(names_from = domain_name, values_from = prop),
          file.path(OUTPUT_DIR, "domain_proportions_per_sample.csv"), row.names = FALSE)

# Cell type proportions per domain
write.csv(prop_ct_domain %>% select(-n) %>% pivot_wider(names_from = domain_name, values_from = prop),
          file.path(OUTPUT_DIR, "celltype_proportions_per_domain.csv"), row.names = FALSE)

# Cell type proportions per domain x condition
write.csv(prop_ct_domain_cond,
          file.path(OUTPUT_DIR, "celltype_proportions_per_domain_condition.csv"), row.names = FALSE)

# =============================================================================
# Done
# =============================================================================
cat("\n========================================\n")
cat("Proportion analysis complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("\nFiles generated:\n")
cat("  - domain_proportions_per_sample.pdf / .csv\n")
cat("  - domain_proportions_by_condition.pdf\n")
cat("  - domain_proportions_by_condition_boxplot.pdf\n")
cat("  - celltype_composition_per_domain.pdf\n")
cat("  - celltype_composition_per_domain_by_condition.pdf\n")
cat("  - celltype_proportions_heatmap_sample_domain.pdf\n")
cat("  - cell_counts_per_domain_sample.pdf\n")
cat("  - celltype_proportions_per_domain.csv\n")
cat("  - celltype_proportions_per_domain_condition.csv\n")
