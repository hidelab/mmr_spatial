# 04_DEG_DESeq_merged_domains.R
# Differential gene expression between WT and DKO within merged Banksy spatial
# domains using pseudobulk aggregation + DESeq2.

library(Seurat)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)

# =============================================================================
# Parameters
# =============================================================================
# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "04_DEG_merged_domains")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# Merged domain column
DOMAIN_COL <- COL_DOMAIN

# =============================================================================
# 1. Load merged-domain Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded merged-domain Seurat object:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

cat("Conditions:", paste(levels(seu$condition), collapse = ", "), "\n")
cat("Samples:", paste(levels(seu$sample), collapse = ", "), "\n")

# Verify domain column exists
if (!DOMAIN_COL %in% colnames(seu@meta.data)) {
  stop("Domain column '", DOMAIN_COL, "' not found. Available columns:\n",
       paste(grep("banksy", colnames(seu@meta.data), value = TRUE), collapse = ", "))
}

seu$domain <- as.character(seu@meta.data[[DOMAIN_COL]])
cat("Merged spatial domains:", length(unique(seu$domain)), "\n")
cat("Domain sizes:\n")
print(sort(table(seu$domain), decreasing = TRUE))
cat("\n")

# =============================================================================
# 2. Pre-filtering diagnostics
# =============================================================================
cat("== Pre-filtering diagnostics ==\n")
diag_dir <- file.path(OUTPUT_DIR, "diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

qc_df <- data.frame(
  cell = colnames(seu),
  total_counts = seu$total_counts,
  n_genes = seu$n_genes_by_counts,
  sample = seu$sample,
  condition = seu$condition,
  domain = seu$domain,
  stringsAsFactors = FALSE
)

# --- 2a. Per-sample summary statistics ---
sample_summary <- qc_df %>%
  group_by(sample, condition) %>%
  summarise(
    n_cells = n(),
    median_counts = median(total_counts),
    mean_counts = mean(total_counts),
    sd_counts = sd(total_counts),
    q25_counts = quantile(total_counts, 0.25),
    q75_counts = quantile(total_counts, 0.75),
    total_library = sum(total_counts),
    median_genes = median(n_genes),
    mean_genes = mean(n_genes),
    .groups = "drop"
  ) %>%
  arrange(condition, sample)

cat("\nPer-sample read count summary:\n")
print(as.data.frame(sample_summary))
write.csv(sample_summary, file.path(diag_dir, "sample_read_summary.csv"), row.names = FALSE)

# --- 2b. Total library size per sample (bar plot) ---
p_lib <- ggplot(sample_summary, aes(x = reorder(sample, -total_library),
                                    y = total_library / 1e6, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Total library size per sample",
       x = NULL, y = "Total UMI counts (millions)", fill = "Condition") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(diag_dir, "library_size_per_sample.pdf"), p_lib, width = 6, height = 4)

# --- 2c. Violin plot: total counts per domain, colored by condition ---
p_vln_domain <- ggplot(qc_df, aes(x = domain, y = total_counts, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = CONDITION_COLORS) +
  geom_hline(yintercept = MIN_TOTAL_COUNTS, linetype = "dashed", color = "red") +
  scale_y_log10() +
  labs(title = "Total counts per cell by merged domain and condition",
       x = "Merged Domain", y = "Total counts (log10)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(diag_dir, "counts_violin_per_domain.pdf"), p_vln_domain, width = 10, height = 5)

# --- 2d. Domain composition: cell types within each merged domain ---
domain_ct_table <- table(seu$domain, seu$Tier1_celltype)
domain_ct_prop <- prop.table(domain_ct_table, margin = 1)

pdf(file.path(diag_dir, "domain_celltype_composition.pdf"), width = 10, height = 8)
pheatmap(as.matrix(domain_ct_prop),
         color = colorRampPalette(c("white", "navy"))(50),
         cluster_rows = TRUE, cluster_cols = TRUE,
         main = "Cell type composition per merged spatial domain",
         fontsize_row = 9, fontsize_col = 8)
dev.off()

# --- 2e. Median counts per sample x domain ---
sample_domain_summary <- qc_df %>%
  group_by(sample, condition, domain) %>%
  summarise(
    n_cells = n(),
    median_counts = median(total_counts),
    total_counts_sum = sum(total_counts),
    .groups = "drop"
  )

p_median_heatmap <- ggplot(sample_domain_summary,
                           aes(x = sample, y = domain, fill = median_counts)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_cells), size = 2.5, color = "black") +
  scale_fill_viridis_c(option = "plasma", trans = "log10") +
  labs(title = "Median counts per cell (color) & cell count (text)",
       subtitle = "Per sample x merged spatial domain",
       x = NULL, y = "Merged Domain", fill = "Median\ncounts") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(diag_dir, "median_counts_sample_domain.pdf"), p_median_heatmap,
       width = 9, height = 7)

write.csv(sample_domain_summary, file.path(diag_dir, "sample_domain_depth_summary.csv"),
          row.names = FALSE)

cat("\nDiagnostic plots saved to:", diag_dir, "\n\n")

# =============================================================================
# 3. Pre-filter cells
# =============================================================================
n_before <- ncol(seu)

keep_cells <- seu$total_counts >= MIN_TOTAL_COUNTS &
              seu$n_genes_by_counts >= MIN_GENES_DETECTED

seu_filt <- subset(seu, cells = colnames(seu)[keep_cells])
n_after <- ncol(seu_filt)

cat(sprintf("Cell filtering: %d -> %d cells (removed %d, %.1f%%)\n",
            n_before, n_after, n_before - n_after,
            100 * (n_before - n_after) / n_before))
cat(sprintf("  Filters: total_counts >= %d, n_genes_by_counts >= %d\n",
            MIN_TOTAL_COUNTS, MIN_GENES_DETECTED))

# Report cells removed per domain
filter_report <- data.frame(
  domain = names(table(seu$domain)),
  before = as.numeric(table(seu$domain)),
  after = as.numeric(table(seu_filt$domain)[names(table(seu$domain))])
)
filter_report$after[is.na(filter_report$after)] <- 0
filter_report$removed <- filter_report$before - filter_report$after
filter_report$pct_removed <- round(100 * filter_report$removed / filter_report$before, 1)
write.csv(filter_report, file.path(OUTPUT_DIR, "cell_filtering_report.csv"), row.names = FALSE)
cat("\nFiltering report per domain:\n")
print(filter_report)
cat("\n")

# =============================================================================
# 4. Pseudobulk aggregation
# =============================================================================
# Aggregate raw counts by sample + merged spatial domain
cat("Aggregating pseudobulk profiles...\n")

# Get raw counts from the RNA assay (not BANKSY)
counts_mat <- GetAssayData(seu_filt, assay = "RNA", layer = "counts")

# Create grouping variable
seu_filt$pseudobulk_group <- paste(seu_filt$sample, seu_filt$domain, sep = "__")

# Count cells per group
group_counts <- table(seu_filt$pseudobulk_group)
cat(sprintf("Total pseudobulk groups: %d\n", length(group_counts)))

# Aggregate counts by summing across cells in each group
groups <- unique(seu_filt$pseudobulk_group)
pb_list <- list()

for (g in groups) {
  cells_in_group <- colnames(seu_filt)[seu_filt$pseudobulk_group == g]
  if (length(cells_in_group) >= MIN_CELLS_PER_SAMPLE) {
    sub_mat <- counts_mat[, cells_in_group, drop = FALSE]
    pb_list[[g]] <- Matrix::rowSums(sub_mat)
  }
}

# Build pseudobulk count matrix
pb_mat <- do.call(cbind, pb_list)
colnames(pb_mat) <- names(pb_list)

# Build metadata for pseudobulk samples
pb_meta <- data.frame(
  group = names(pb_list),
  sample = sapply(strsplit(names(pb_list), "__"), `[`, 1),
  domain = sapply(strsplit(names(pb_list), "__"), `[`, 2),
  n_cells = as.numeric(group_counts[names(pb_list)]),
  stringsAsFactors = FALSE
)
pb_meta$condition <- ifelse(grepl("^WT", pb_meta$sample), REFERENCE_CONDITION, TEST_CONDITION)
rownames(pb_meta) <- pb_meta$group

cat(sprintf("Pseudobulk matrix: %d genes x %d samples (after min %d cells filter)\n",
            nrow(pb_mat), ncol(pb_mat), MIN_CELLS_PER_SAMPLE))

# --- Save pseudobulk profiles per domain ---
pb_dir <- file.path(OUTPUT_DIR, "pseudobulk_profiles")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)

for (dom in unique(pb_meta$domain)) {
  dom_samples <- pb_meta$group[pb_meta$domain == dom]
  dom_pb <- pb_mat[, dom_samples, drop = FALSE]
  colnames(dom_pb) <- pb_meta[dom_samples, "sample"]
  dom_safe <- gsub("[^A-Za-z0-9_]", "_", dom)
  dom_name <- DOMAIN_NAMES[dom]
  write.csv(as.data.frame(dom_pb), file.path(pb_dir, paste0("pseudobulk_domain_", dom_safe, "_", dom_name, ".csv")))
}
cat(sprintf("Pseudobulk profiles saved to: %s\n", pb_dir))

# =============================================================================
# 5. Run DESeq2 per merged spatial domain
# =============================================================================
domains <- unique(pb_meta$domain)
cat(sprintf("\nRunning DESeq2 for %d merged spatial domains...\n", length(domains)))

all_results <- list()
summary_stats <- data.frame()

for (dom in domains) {
  cat(sprintf("\n--- Merged Domain %s ---\n", dom))

  # Subset to this domain
  dom_samples <- pb_meta$group[pb_meta$domain == dom]
  dom_meta <- pb_meta[dom_samples, ]
  dom_counts <- pb_mat[, dom_samples, drop = FALSE]

  # Check we have replicates in both conditions
  n_wt <- sum(dom_meta$condition == "WT")
  n_dko <- sum(dom_meta$condition == "DKO")
  cat(sprintf("  Samples: %d WT, %d DKO\n", n_wt, n_dko))

  if (n_wt < 2 || n_dko < 2) {
    cat("  SKIPPED: need at least 2 samples per condition.\n")
    next
  }

  # Filter genes by cell-level coverage: keep genes expressed in >= MIN_GENE_EXPR_FRAC of cells
  dom_cells <- colnames(seu_filt)[seu_filt$domain == dom]
  dom_cell_counts <- counts_mat[, dom_cells, drop = FALSE]
  gene_det_frac <- Matrix::rowSums(dom_cell_counts > 0) / length(dom_cells)
  genes_pass_frac <- names(gene_det_frac[gene_det_frac >= MIN_GENE_EXPR_FRAC])
  dom_counts <- dom_counts[intersect(rownames(dom_counts), genes_pass_frac), , drop = FALSE]
  cat(sprintf("  Genes expressed in >=%.0f%% of cells: %d\n",
              MIN_GENE_EXPR_FRAC * 100, nrow(dom_counts)))

  # Filter lowly expressed genes: require detection in a minimum number of samples
  n_samples_detected <- rowSums(dom_counts >= MIN_COUNTS_PER_SAMPLE)
  keep_genes <- n_samples_detected >= MIN_SAMPLES_EXPRESSED
  dom_counts <- dom_counts[keep_genes, , drop = FALSE]
  cat(sprintf("  Genes after sample-level filter (>=%d counts in >=%d samples): %d\n",
              MIN_COUNTS_PER_SAMPLE, MIN_SAMPLES_EXPRESSED, nrow(dom_counts)))

  if (nrow(dom_counts) < 10) {
    cat("  SKIPPED: too few genes after filtering.\n")
    next
  }

  # DESeq2
  dom_meta$condition <- factor(dom_meta$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  dds <- DESeqDataSetFromMatrix(
    countData = round(dom_counts),  # ensure integer
    colData = dom_meta,
    design = ~ condition
  )

  dds <- tryCatch({
    DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  ERROR in DESeq: %s\n", e$message))
    return(NULL)
  })

  if (is.null(dds)) next

  res <- results(dds, contrast = c("condition", TEST_CONDITION, REFERENCE_CONDITION))
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$domain <- dom
  res_df <- res_df[order(res_df$pvalue), ]

  # Classify significance
  res_df$significance <- "NS"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange > LFC_THRESHOLD] <- "Up in DKO"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange < -LFC_THRESHOLD] <- "Down in DKO"

  all_results[[dom]] <- res_df

  # Summary
  n_up <- sum(res_df$significance == "Up in DKO")
  n_down <- sum(res_df$significance == "Down in DKO")
  n_tested <- sum(!is.na(res_df$padj))
  cat(sprintf("  DEGs (FDR < %.2f, |LFC| > %.1f): %d up, %d down (of %d tested)\n",
              FDR_THRESHOLD, LFC_THRESHOLD, n_up, n_down, n_tested))

  summary_stats <- rbind(summary_stats, data.frame(
    domain = dom,
    n_samples_WT = n_wt,
    n_samples_DKO = n_dko,
    n_genes_tested = n_tested,
    n_up_DKO = n_up,
    n_down_DKO = n_down,
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# 6. Save results
# =============================================================================
# Combined results table
combined_res <- bind_rows(all_results)
write.csv(combined_res, file.path(OUTPUT_DIR, "DEG_all_merged_domains.csv"), row.names = FALSE)

# Summary table
write.csv(summary_stats, file.path(OUTPUT_DIR, "DEG_summary.csv"), row.names = FALSE)
cat("\n\n== DEG Summary ==\n")
print(summary_stats)

# Per-domain result files
dom_dir <- file.path(OUTPUT_DIR, "per_domain")
dir.create(dom_dir, recursive = TRUE, showWarnings = FALSE)
for (dom in names(all_results)) {
  dom_safe <- gsub("[^A-Za-z0-9_]", "_", dom)
  dom_name <- DOMAIN_NAMES[dom]
  write.csv(all_results[[dom]],
            file.path(dom_dir, paste0("DEG_merged_domain_", dom_safe, "_", dom_name, ".csv")),
            row.names = FALSE)
}

# =============================================================================
# 7. Volcano plots
# =============================================================================
cat("\nGenerating volcano plots...\n")
volcano_dir <- file.path(OUTPUT_DIR, "volcano_plots")
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

for (dom in names(all_results)) {
  res_df <- all_results[[dom]]
  res_df <- res_df[!is.na(res_df$padj) & is.finite(res_df$padj), ]

  if (nrow(res_df) == 0) next

  dom_safe <- gsub("[^A-Za-z0-9_]", "_", dom)
  dom_name <- DOMAIN_NAMES[dom]

  # Top genes to label
  top_genes <- res_df %>%
    filter(significance != "NS") %>%
    arrange(padj) %>%
    head(VOLCANO_TOP_N)

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Up in DKO" = "#D73027",
                                  "Down in DKO" = "#4575B4",
                                  "NS" = "grey60")) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(FDR_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_genes,
                    aes(label = gene),
                    size = 2.5, max.overlaps = 15,
                    color = "black") +
    labs(title = paste0("DKO vs WT — Domain ", dom, " (", gsub("_", " ", dom_name), ")"),
         x = "Log2 Fold Change (DKO / WT)",
         y = "-log10(adjusted p-value)",
         color = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ggsave(file.path(volcano_dir, paste0("volcano_merged_domain_", dom_safe, "_", dom_name, ".pdf")),
         p, width = 7, height = 6)
}

# =============================================================================
# 8. Summary bar plot: number of DEGs per merged domain
# =============================================================================
if (nrow(summary_stats) > 0) {
  deg_long <- summary_stats %>%
    select(domain, n_up_DKO, n_down_DKO) %>%
    pivot_longer(cols = c(n_up_DKO, n_down_DKO),
                 names_to = "direction", values_to = "n_DEGs") %>%
    mutate(direction = ifelse(direction == "n_up_DKO", "Up in DKO", "Down in DKO"),
           n_DEGs_signed = ifelse(direction == "Up in DKO", n_DEGs, -n_DEGs))

  p_bar <- ggplot(deg_long, aes(x = reorder(domain, -abs(n_DEGs_signed)),
                                y = n_DEGs_signed, fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Up in DKO" = "#D73027", "Down in DKO" = "#4575B4")) +
    geom_hline(yintercept = 0, color = "black") +
    coord_flip() +
    labs(title = sprintf("DEGs per merged domain (FDR < %.2f, |LFC| > %.1f)",
                         FDR_THRESHOLD, LFC_THRESHOLD),
         x = "Merged Domain", y = "Number of DEGs", fill = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "DEG_count_per_merged_domain.pdf"), p_bar, width = 8, height = 6)
}

# =============================================================================
# 9. Heatmap of top DEGs across merged domains (signed -log10 p-value)
# =============================================================================
cat("Generating cross-domain heatmap (signed -log10 p-value)...\n")

# Collect all significant DEGs across domains (no slice_head limit)
sig_genes <- combined_res %>%
  filter(significance != "NS") %>%
  pull(gene) %>%
  unique()

if (length(sig_genes) > 0) {
  # Compute signed -log10(pvalue): sign(LFC) * -log10(pvalue)
  # For non-significant comparisons, set value to 0
  heatmap_df <- combined_res %>%
    filter(gene %in% sig_genes) %>%
    mutate(
      signed_logp = ifelse(
        significance != "NS",
        sign(log2FoldChange) * -log10(pvalue),
        0
      ),
      domain_label = DOMAIN_NAMES[domain]
    ) %>%
    select(gene, domain_label, signed_logp) %>%
    pivot_wider(names_from = domain_label, values_from = signed_logp) %>%
    tibble::column_to_rownames("gene")

  # Replace NA (gene not tested in a domain) with 0
  heatmap_df[is.na(heatmap_df)] <- 0

  # Cap extreme values for visualization at 99th percentile
  cap_val <- quantile(abs(as.matrix(heatmap_df)[as.matrix(heatmap_df) != 0]), 0.99)
  heatmap_df[heatmap_df > cap_val] <- cap_val
  heatmap_df[heatmap_df < -cap_val] <- -cap_val

  # Only plot if we have enough genes
  if (nrow(heatmap_df) >= 3 && ncol(heatmap_df) >= 2) {
    pdf(file.path(OUTPUT_DIR, "DEG_heatmap_top_genes.pdf"),
        width = max(8, ncol(heatmap_df) * 0.5 + 3),
        height = max(8, nrow(heatmap_df) * 0.25 + 2))
    pheatmap(as.matrix(heatmap_df),
             color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
             breaks = seq(-cap_val, cap_val, length.out = 101),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             treeheight_row = 0,
             main = "Significant DEGs — signed -log10(p-value) (DKO vs WT) by Merged Domain",
             fontsize_row = 10,
             fontsize_col = 10)
    dev.off()
  }
}

# =============================================================================
# 10. MA plots
# =============================================================================
cat("Generating MA plots...\n")
ma_dir <- file.path(OUTPUT_DIR, "MA_plots")
dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)

for (dom in names(all_results)) {
  res_df <- all_results[[dom]]
  res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$baseMean), ]

  if (nrow(res_df) == 0) next

  dom_safe <- gsub("[^A-Za-z0-9_]", "_", dom)
  dom_name <- DOMAIN_NAMES[dom]

  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = c("Up in DKO" = "#D73027",
                                  "Down in DKO" = "#4575B4",
                                  "NS" = "grey60")) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    labs(title = paste0("MA plot — DKO vs WT — Domain ", dom, " (", gsub("_", " ", dom_name), ")"),
         x = "log10(Mean Expression + 1)",
         y = "Log2 Fold Change",
         color = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ggsave(file.path(ma_dir, paste0("MA_merged_domain_", dom_safe, "_", dom_name, ".pdf")), p_ma, width = 7, height = 5)
}

# =============================================================================
# 11. Pseudobulk QC: sample-level cell counts and library sizes
# =============================================================================
cat("Generating pseudobulk QC plots...\n")

# Cell count per pseudobulk sample
p_qc <- ggplot(pb_meta, aes(x = reorder(group, -n_cells), y = n_cells, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CONDITION_COLORS) +
  coord_flip() +
  labs(title = "Cells per pseudobulk sample (by merged domain)", x = NULL, y = "Number of cells") +
  theme_publication() +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(OUTPUT_DIR, "pseudobulk_cell_counts.pdf"), p_qc,
       width = 10, height = max(6, nrow(pb_meta) * 0.2 + 2))

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("DEG analysis by merged spatial domain complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("Parameters:\n"))
cat(sprintf("  DOMAIN_COL = %s\n", DOMAIN_COL))
cat(sprintf("  MIN_TOTAL_COUNTS = %d\n", MIN_TOTAL_COUNTS))
cat(sprintf("  MIN_GENES_DETECTED = %d\n", MIN_GENES_DETECTED))
cat(sprintf("  MIN_CELLS_PER_SAMPLE = %d\n", MIN_CELLS_PER_SAMPLE))
cat(sprintf("  FDR_THRESHOLD = %.2f\n", FDR_THRESHOLD))
cat(sprintf("  LFC_THRESHOLD = %.1f\n", LFC_THRESHOLD))
cat(sprintf("Merged domains tested: %d\n", nrow(summary_stats)))
cat(sprintf("Total significant DEGs: %d up, %d down\n",
            sum(summary_stats$n_up_DKO), sum(summary_stats$n_down_DKO)))
cat("\nFiles generated:\n")
cat("  - DEG_all_merged_domains.csv (combined results)\n")
cat("  - DEG_summary.csv (per-domain summary)\n")
cat("  - per_domain/ (individual result CSVs)\n")
cat("  - volcano_plots/ (per-domain volcano PDFs)\n")
cat("  - MA_plots/ (per-domain MA PDFs)\n")
cat("  - DEG_heatmap_top_genes.pdf\n")
cat("  - DEG_count_per_merged_domain.pdf\n")
cat("  - pseudobulk_cell_counts.pdf\n")
cat("  - cell_filtering_report.csv\n")
cat("  - diagnostics/domain_celltype_composition.pdf\n")
