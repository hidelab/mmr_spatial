# 03_DEG_DESeq_tier1.R
# Differential gene expression between WT and DKO within Tier1 cell types
# using pseudobulk aggregation + DESeq2.
# Cells are pre-filtered by minimum total counts before pseudobulking.

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

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "03_DEG_pseudobulk")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# 1. Load Seurat object
# =============================================================================
seu <- readRDS(SEURAT_RDS)
cat("Loaded Seurat object:", ncol(seu), "cells,", nrow(seu), "genes\n")

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
cat("Cell types:", length(levels(seu$Tier1_celltype)), "\n\n")

# =============================================================================
# 2. Pre-filtering diagnostics: read distribution per cell and sample
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
  celltype = seu$Tier1_celltype,
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

# --- 2c. Distribution of total counts per cell, split by sample ---
p_counts_sample <- ggplot(qc_df, aes(x = total_counts, fill = condition)) +
  geom_histogram(bins = 80, alpha = 0.7) +
  facet_wrap(~ sample, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = CONDITION_COLORS) +
  geom_vline(xintercept = MIN_TOTAL_COUNTS, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(title = "Total counts per cell by sample",
       subtitle = paste("Red line = MIN_TOTAL_COUNTS =", MIN_TOTAL_COUNTS),
       x = "Total counts (log10)", y = "Number of cells") +
  theme_publication()

ggsave(file.path(diag_dir, "counts_per_cell_by_sample.pdf"), p_counts_sample, width = 12, height = 6)

# --- 2d. Violin plot: total counts per sample ---
p_vln <- ggplot(qc_df, aes(x = sample, y = total_counts, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = CONDITION_COLORS) +
  geom_hline(yintercept = MIN_TOTAL_COUNTS, linetype = "dashed", color = "red") +
  scale_y_log10() +
  labs(title = "Total counts per cell (violin + box)",
       subtitle = paste("Red line = MIN_TOTAL_COUNTS =", MIN_TOTAL_COUNTS),
       x = NULL, y = "Total counts (log10)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(diag_dir, "counts_violin_per_sample.pdf"), p_vln, width = 7, height = 5)

# --- 2e. Violin plot: total counts per cell type, colored by condition ---
p_vln_ct <- ggplot(qc_df, aes(x = celltype, y = total_counts, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7, position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.15, outlier.size = 0.2, alpha = 0.8,
               position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = CONDITION_COLORS) +
  geom_hline(yintercept = MIN_TOTAL_COUNTS, linetype = "dashed", color = "red") +
  scale_y_log10() +
  labs(title = "Total counts per cell by cell type and condition",
       x = NULL, y = "Total counts (log10)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

ggsave(file.path(diag_dir, "counts_violin_per_celltype.pdf"), p_vln_ct, width = 12, height = 5)

# --- 2f. Genes detected per cell by sample ---
p_genes <- ggplot(qc_df, aes(x = sample, y = n_genes, fill = condition)) +
  geom_violin(scale = "width", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = CONDITION_COLORS) +
  geom_hline(yintercept = MIN_GENES_DETECTED, linetype = "dashed", color = "red") +
  scale_y_log10() +
  labs(title = "Genes detected per cell",
       subtitle = paste("Red line = MIN_GENES_DETECTED =", MIN_GENES_DETECTED),
       x = NULL, y = "Genes detected (log10)") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(diag_dir, "genes_per_cell_by_sample.pdf"), p_genes, width = 7, height = 5)

# --- 2g. Median counts per sample x celltype (potential pseudobulk bias) ---
sample_ct_summary <- qc_df %>%
  group_by(sample, condition, celltype) %>%
  summarise(
    n_cells = n(),
    median_counts = median(total_counts),
    total_counts_sum = sum(total_counts),
    .groups = "drop"
  )

p_median_heatmap <- ggplot(sample_ct_summary,
                           aes(x = sample, y = celltype, fill = median_counts)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_cells), size = 2.2, color = "black") +
  scale_fill_viridis_c(option = "plasma", trans = "log10") +
  labs(title = "Median counts per cell (color) & cell count (text)",
       subtitle = "Check for sample-level depth differences that affect pseudobulk",
       x = NULL, y = NULL, fill = "Median\ncounts") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(diag_dir, "median_counts_sample_celltype.pdf"), p_median_heatmap,
       width = 9, height = 8)

# --- 2h. Total library per pseudobulk group (before filtering) ---
# This shows how much of the pseudobulk signal is driven by library depth
pb_depth <- qc_df %>%
  group_by(sample, condition, celltype) %>%
  summarise(
    total_lib = sum(total_counts),
    n_cells = n(),
    mean_counts_per_cell = mean(total_counts),
    .groups = "drop"
  )

p_pb_depth <- ggplot(pb_depth, aes(x = n_cells, y = total_lib / 1e3, color = condition)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = CONDITION_COLORS) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
  labs(title = "Pseudobulk library size vs cell count",
       subtitle = "Deviation from line suggests per-cell depth differences",
       x = "Number of cells in group",
       y = "Total pseudobulk library (thousands)") +
  theme_publication()

ggsave(file.path(diag_dir, "pseudobulk_depth_vs_ncells.pdf"), p_pb_depth, width = 7, height = 5)

write.csv(sample_ct_summary, file.path(diag_dir, "sample_celltype_depth_summary.csv"),
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

# Report cells removed per cell type
filter_report <- data.frame(
  celltype = names(table(seu$Tier1_celltype)),
  before = as.numeric(table(seu$Tier1_celltype)),
  after = as.numeric(table(seu_filt$Tier1_celltype)[names(table(seu$Tier1_celltype))])
)
filter_report$after[is.na(filter_report$after)] <- 0
filter_report$removed <- filter_report$before - filter_report$after
filter_report$pct_removed <- round(100 * filter_report$removed / filter_report$before, 1)
write.csv(filter_report, file.path(OUTPUT_DIR, "cell_filtering_report.csv"), row.names = FALSE)
cat("\nFiltering report per cell type:\n")
print(filter_report)
cat("\n")

# =============================================================================
# 4. Pseudobulk aggregation
# =============================================================================
# Aggregate raw counts by sample + cell type
cat("Aggregating pseudobulk profiles...\n")

# Get raw counts
counts_mat <- GetAssayData(seu_filt, assay = DefaultAssay(seu_filt), layer = "counts")

# Create grouping variable
seu_filt$pseudobulk_group <- paste(seu_filt$sample, seu_filt$Tier1_celltype, sep = "__")

# Count cells per group
group_counts <- table(seu_filt$pseudobulk_group)
cat(sprintf("Total pseudobulk groups: %d\n", length(group_counts)))

# Aggregate counts by summing across cells in each group
groups <- unique(seu_filt$pseudobulk_group)
pb_list <- list()

for (g in groups) {
  cells_in_group <- colnames(seu_filt)[seu_filt$pseudobulk_group == g]
  if (length(cells_in_group) >= MIN_CELLS_PER_SAMPLE) {
    if (length(cells_in_group) == 1) {
      pb_list[[g]] <- counts_mat[, cells_in_group]
    } else {
      pb_list[[g]] <- Matrix::rowSums(counts_mat[, cells_in_group])
    }
  }
}

# Build pseudobulk count matrix
pb_mat <- do.call(cbind, pb_list)
colnames(pb_mat) <- names(pb_list)

# Build metadata for pseudobulk samples
pb_meta <- data.frame(
  group = names(pb_list),
  sample = sapply(strsplit(names(pb_list), "__"), `[`, 1),
  celltype = sapply(strsplit(names(pb_list), "__"), `[`, 2),
  n_cells = as.numeric(group_counts[names(pb_list)]),
  stringsAsFactors = FALSE
)
pb_meta$condition <- ifelse(pb_meta$sample %in% WT_SAMPLES, REFERENCE_CONDITION, TEST_CONDITION)
rownames(pb_meta) <- pb_meta$group

cat(sprintf("Pseudobulk matrix: %d genes x %d samples (after min %d cells filter)\n",
            nrow(pb_mat), ncol(pb_mat), MIN_CELLS_PER_SAMPLE))

# --- Save pseudobulk profiles per cell type ---
pb_dir <- file.path(OUTPUT_DIR, "pseudobulk_profiles")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)

for (ct in unique(pb_meta$celltype)) {
  ct_samples <- pb_meta$group[pb_meta$celltype == ct]
  ct_pb <- pb_mat[, ct_samples, drop = FALSE]
  # Simplify column names to sample ID
  colnames(ct_pb) <- pb_meta[ct_samples, "sample"]
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  write.csv(as.data.frame(ct_pb), file.path(pb_dir, paste0("pseudobulk_", ct_safe, ".csv")))
}
cat(sprintf("Pseudobulk profiles saved to: %s\n", pb_dir))

# =============================================================================
# 5. Run DESeq2 per cell type
# =============================================================================
cell_types <- unique(pb_meta$celltype)
cat(sprintf("\nRunning DESeq2 for %d cell types...\n", length(cell_types)))

# Significance labels (derived from parameters)
label_up <- paste("Up in", TEST_CONDITION)
label_down <- paste("Down in", TEST_CONDITION)

all_results <- list()
summary_stats <- data.frame()

for (ct in cell_types) {
  cat(sprintf("\n--- %s ---\n", ct))

  # Subset to this cell type
  ct_samples <- pb_meta$group[pb_meta$celltype == ct]
  ct_meta <- pb_meta[ct_samples, ]
  ct_counts <- pb_mat[, ct_samples, drop = FALSE]

  # Check we have replicates in both conditions
  n_wt <- sum(ct_meta$condition == REFERENCE_CONDITION)
  n_dko <- sum(ct_meta$condition == TEST_CONDITION)
  cat(sprintf("  Samples: %d %s, %d %s\n", n_wt, REFERENCE_CONDITION, n_dko, TEST_CONDITION))


  if (n_wt < 2 || n_dko < 2) {
    cat("  SKIPPED: need at least 2 samples per condition.\n")
    next
  }

  # Filter genes by cell-level coverage: keep genes expressed in >= MIN_GENE_EXPR_FRAC of cells
  ct_cells <- colnames(seu_filt)[seu_filt$Tier1_celltype == ct]
  ct_cell_counts <- counts_mat[, ct_cells, drop = FALSE]
  gene_det_frac <- Matrix::rowSums(ct_cell_counts > 0) / length(ct_cells)
  genes_pass_frac <- names(gene_det_frac[gene_det_frac >= MIN_GENE_EXPR_FRAC])
  ct_counts <- ct_counts[intersect(rownames(ct_counts), genes_pass_frac), , drop = FALSE]
  cat(sprintf("  Genes expressed in >=%.0f%% of cells: %d\n",
              MIN_GENE_EXPR_FRAC * 100, nrow(ct_counts)))

  # Filter lowly expressed genes: require detection in a minimum number of samples
  n_samples_detected <- rowSums(ct_counts >= MIN_COUNTS_PER_SAMPLE)
  keep_genes <- n_samples_detected >= MIN_SAMPLES_EXPRESSED
  ct_counts <- ct_counts[keep_genes, , drop = FALSE]
  cat(sprintf("  Genes after sample-level filter (>=%d counts in >=%d samples): %d\n",
              MIN_COUNTS_PER_SAMPLE, MIN_SAMPLES_EXPRESSED, nrow(ct_counts)))

  if (nrow(ct_counts) < 10) {
    cat("  SKIPPED: too few genes after filtering.\n")
    next
  }

  # DESeq2
  ct_meta$condition <- factor(ct_meta$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  dds <- DESeqDataSetFromMatrix(
    countData = round(ct_counts),  # ensure integer
    colData = ct_meta,
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
  res_df$celltype <- ct
  res_df <- res_df[order(res_df$pvalue), ]

  # Classify significance
  res_df$significance <- "NS"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange > LFC_THRESHOLD] <- label_up
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange < -LFC_THRESHOLD] <- label_down

  all_results[[ct]] <- res_df

  # Summary
  n_up <- sum(res_df$significance == label_up)
  n_down <- sum(res_df$significance == label_down)
  n_tested <- sum(!is.na(res_df$padj))
  cat(sprintf("  DEGs (FDR < %.2f, |LFC| > %.1f): %d up, %d down (of %d tested)\n",
              FDR_THRESHOLD, LFC_THRESHOLD, n_up, n_down, n_tested))

  summary_stats <- rbind(summary_stats, data.frame(
    celltype = ct,
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
write.csv(combined_res, file.path(OUTPUT_DIR, "DEG_all_celltypes.csv"), row.names = FALSE)

# Summary table
write.csv(summary_stats, file.path(OUTPUT_DIR, "DEG_summary.csv"), row.names = FALSE)
cat("\n\n== DEG Summary ==\n")
print(summary_stats)

# Per-cell-type result files
ct_dir <- file.path(OUTPUT_DIR, "per_celltype")
dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)
for (ct in names(all_results)) {
  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  write.csv(all_results[[ct]],
            file.path(ct_dir, paste0("DEG_", ct_safe, ".csv")),
            row.names = FALSE)
}

# =============================================================================
# 7. Volcano plots
# =============================================================================
cat("\nGenerating volcano plots...\n")
volcano_dir <- file.path(OUTPUT_DIR, "volcano_plots")
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

for (ct in names(all_results)) {
  res_df <- all_results[[ct]]
  res_df <- res_df[!is.na(res_df$padj) & is.finite(res_df$padj), ]

  if (nrow(res_df) == 0) next

  # Top genes to label
  top_genes <- res_df %>%
    filter(significance != "NS") %>%
    arrange(padj) %>%
    head(VOLCANO_TOP_N)

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = setNames(c("#D73027", "#4575B4", "grey60"),
                                         c(label_up, label_down, "NS"))) +
    geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(FDR_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_genes,
                    aes(label = gene),
                    size = 2.5, max.overlaps = 15,
                    color = "black") +
    labs(title = paste(TEST_CONDITION, "vs", REFERENCE_CONDITION, "—", ct),
         x = paste0("Log2 Fold Change (", TEST_CONDITION, " / ", REFERENCE_CONDITION, ")"),
         y = "-log10(adjusted p-value)",
         color = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  ggsave(file.path(volcano_dir, paste0("volcano_", ct_safe, ".pdf")),
         p, width = 7, height = 6)
}

# =============================================================================
# 8. Summary bar plot: number of DEGs per cell type
# =============================================================================
if (nrow(summary_stats) > 0) {
  deg_long <- summary_stats %>%
    select(celltype, n_up_DKO, n_down_DKO) %>%
    pivot_longer(cols = c(n_up_DKO, n_down_DKO),
                 names_to = "direction", values_to = "n_DEGs") %>%
    mutate(direction = ifelse(direction == "n_up_DKO", label_up, label_down),
           n_DEGs_signed = ifelse(direction == label_up, n_DEGs, -n_DEGs))

  p_bar <- ggplot(deg_long, aes(x = reorder(celltype, -abs(n_DEGs_signed)),
                                y = n_DEGs_signed, fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = setNames(c("#D73027", "#4575B4"), c(label_up, label_down))) +
    geom_hline(yintercept = 0, color = "black") +
    coord_flip() +
    labs(title = sprintf("DEGs per cell type (FDR < %.2f, |LFC| > %.1f)",
                         FDR_THRESHOLD, LFC_THRESHOLD),
         x = NULL, y = "Number of DEGs", fill = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "DEG_count_per_celltype.pdf"), p_bar, width = 8, height = 6)
}

# =============================================================================
# 9. Heatmap of top DEGs across cell types (signed -log10 p-value)
# =============================================================================
cat("Generating cross-cell-type heatmap (signed -log10 p-value)...\n")

# Collect all significant DEGs across cell types (no slice_head limit)
sig_genes <- combined_res %>%
  filter(significance != "NS") %>%
  pull(gene) %>%
  unique()

if (length(sig_genes) > 0) {
  # Compute signed -log10(pvalue): sign(LFC) * -log10(pvalue)
  # For non-significant comparisons, set value to 0
  heatmap_df <- combined_res %>%
    filter(gene %in% sig_genes) %>%
    mutate(signed_logp = ifelse(
      significance != "NS",
      sign(log2FoldChange) * -log10(pvalue),
      0
    )) %>%
    select(gene, celltype, signed_logp) %>%
    pivot_wider(names_from = celltype, values_from = signed_logp) %>%
    tibble::column_to_rownames("gene")

  # Replace NA (gene not tested in a cell type) with 0
  heatmap_df[is.na(heatmap_df)] <- 0

  # Cap extreme values for visualization
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
             main = paste("Significant DEGs — signed -log10(p-value)", paste0("(", TEST_CONDITION, " vs ", REFERENCE_CONDITION, ")")),
             fontsize_row = 10,
             fontsize_col = 10,
            treeheight_row = 0)
    dev.off()
  }
}

# =============================================================================
# 10. MA plots
# =============================================================================
cat("Generating MA plots...\n")
ma_dir <- file.path(OUTPUT_DIR, "MA_plots")
dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)

for (ct in names(all_results)) {
  res_df <- all_results[[ct]]
  res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$baseMean), ]

  if (nrow(res_df) == 0) next

  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = setNames(c("#D73027", "#4575B4", "grey60"),
                                         c(label_up, label_down, "NS"))) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    labs(title = paste("MA plot —", TEST_CONDITION, "vs", REFERENCE_CONDITION, "—", ct),
         x = "log10(Mean Expression + 1)",
         y = "Log2 Fold Change",
         color = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ct_safe <- gsub("[^A-Za-z0-9_]", "_", ct)
  ggsave(file.path(ma_dir, paste0("MA_", ct_safe, ".pdf")), p_ma, width = 7, height = 5)
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
  labs(title = "Cells per pseudobulk sample", x = NULL, y = "Number of cells") +
  theme_publication() +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(OUTPUT_DIR, "pseudobulk_cell_counts.pdf"), p_qc,
       width = 10, height = max(6, nrow(pb_meta) * 0.2 + 2))

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("DEG analysis complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("Parameters:\n"))
cat(sprintf("  MIN_TOTAL_COUNTS = %d\n", MIN_TOTAL_COUNTS))
cat(sprintf("  MIN_GENES_DETECTED = %d\n", MIN_GENES_DETECTED))
cat(sprintf("  MIN_CELLS_PER_SAMPLE = %d\n", MIN_CELLS_PER_SAMPLE))
cat(sprintf("  FDR_THRESHOLD = %.2f\n", FDR_THRESHOLD))
cat(sprintf("  LFC_THRESHOLD = %.1f\n", LFC_THRESHOLD))
cat(sprintf("Cell types tested: %d\n", nrow(summary_stats)))
cat(sprintf("Total significant DEGs: %d up, %d down\n",
            sum(summary_stats$n_up_DKO), sum(summary_stats$n_down_DKO)))
cat("\nFiles generated:\n")
cat("  - DEG_all_celltypes.csv (combined results)\n")
cat("  - DEG_summary.csv (per-cell-type summary)\n")
cat("  - per_celltype/ (individual result CSVs)\n")
cat("  - volcano_plots/ (per-cell-type volcano PDFs)\n")
cat("  - MA_plots/ (per-cell-type MA PDFs)\n")
cat("  - DEG_heatmap_top_genes.pdf\n")
cat("  - DEG_count_per_celltype.pdf\n")
cat("  - pseudobulk_cell_counts.pdf\n")
cat("  - cell_filtering_report.csv\n")

