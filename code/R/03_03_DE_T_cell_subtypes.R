# 03_03_DE_T_cell_subtypes.R
# Differential gene expression between DKO and WT within each CD8+ T cell
# subtype (clusters from 03_01), using pseudobulk aggregation (sample x subtype)
# + DESeq2. Same filters and model as 04_00_DEG_DESeq_tier1.R. Cells were
# already QC-filtered (MIN_TOTAL_COUNTS, MIN_GENES_DETECTED) in 03_01.
#
# Read-only on existing RDS files.

library(Seurat)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# Parameters
# =============================================================================
# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

TCELL_DIR  <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters")
OUTPUT_DIR <- file.path(TCELL_DIR, "DEG_pseudobulk")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

TCELL_RDS              <- file.path(TCELL_DIR, "DKO_CD8T_clusters.rds")  # from 03_01
TCELL_CLUSTER_COL      <- "tcell_clust_res.0.4"
TCELL_EXCLUDE_CLUSTERS <- character(0)   # e.g. "C0" to drop a contaminated cluster
EXPR_ASSAY             <- "RNA"

# =============================================================================
# 1. Load CD8+ T object
# =============================================================================
seu <- readRDS(TCELL_RDS)
cat("Loaded:", ncol(seu), "CD8+ T cells,", nrow(seu), "genes\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample, levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))
seu$subtype <- factor(paste0("C", seu@meta.data[[TCELL_CLUSTER_COL]]),
                      levels = paste0("C", levels(seu@meta.data[[TCELL_CLUSTER_COL]])))

if (length(TCELL_EXCLUDE_CLUSTERS) > 0) {
  seu <- subset(seu, cells = colnames(seu)[!seu$subtype %in% TCELL_EXCLUDE_CLUSTERS])
  cat("Excluded subtypes:", paste(TCELL_EXCLUDE_CLUSTERS, collapse = ", "), "\n")
}
seu$subtype <- droplevels(seu$subtype)

# Defensive re-check of the cell-level QC applied in 03_01
keep_cells <- seu$total_counts >= MIN_TOTAL_COUNTS &
              seu$n_genes_by_counts >= MIN_GENES_DETECTED
seu <- subset(seu, cells = colnames(seu)[keep_cells])

cat("Cells per sample x subtype:\n")
print(table(seu$subtype, seu$sample))

# =============================================================================
# 2. Pseudobulk aggregation (sample x subtype)
# =============================================================================
cat("\nAggregating pseudobulk profiles...\n")
counts_mat <- GetAssayData(seu, assay = EXPR_ASSAY, layer = "counts")

seu$pseudobulk_group <- paste(seu$sample, seu$subtype, sep = "__")
group_counts <- table(seu$pseudobulk_group)

pb_list <- list()
for (g in names(group_counts)) {
  if (group_counts[[g]] >= MIN_CELLS_PER_SAMPLE) {
    cells_in_group <- colnames(seu)[seu$pseudobulk_group == g]
    pb_list[[g]] <- Matrix::rowSums(counts_mat[, cells_in_group, drop = FALSE])
  }
}

pb_mat <- do.call(cbind, pb_list)
colnames(pb_mat) <- names(pb_list)

pb_meta <- data.frame(
  group = names(pb_list),
  sample = sapply(strsplit(names(pb_list), "__"), `[`, 1),
  subtype = sapply(strsplit(names(pb_list), "__"), `[`, 2),
  n_cells = as.numeric(group_counts[names(pb_list)]),
  stringsAsFactors = FALSE
)
pb_meta$condition <- ifelse(pb_meta$sample %in% WT_SAMPLES, REFERENCE_CONDITION, TEST_CONDITION)
rownames(pb_meta) <- pb_meta$group

cat(sprintf("Pseudobulk matrix: %d genes x %d samples (after min %d cells filter)\n",
            nrow(pb_mat), ncol(pb_mat), MIN_CELLS_PER_SAMPLE))
write.csv(pb_meta, file.path(OUTPUT_DIR, "pseudobulk_samples.csv"), row.names = FALSE)

pb_dir <- file.path(OUTPUT_DIR, "pseudobulk_profiles")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)
for (st in unique(pb_meta$subtype)) {
  st_samples <- pb_meta$group[pb_meta$subtype == st]
  st_pb <- pb_mat[, st_samples, drop = FALSE]
  colnames(st_pb) <- pb_meta[st_samples, "sample"]
  write.csv(as.data.frame(st_pb), file.path(pb_dir, paste0("pseudobulk_", st, ".csv")))
}

# =============================================================================
# 3. Run DESeq2 per subtype
# =============================================================================
subtypes <- intersect(levels(seu$subtype), unique(pb_meta$subtype))
cat(sprintf("\nRunning DESeq2 for %d subtypes...\n", length(subtypes)))

label_up <- paste("Up in", TEST_CONDITION)
label_down <- paste("Down in", TEST_CONDITION)

all_results <- list()
summary_stats <- data.frame()

for (st in subtypes) {
  cat(sprintf("\n--- %s ---\n", st))

  st_samples <- pb_meta$group[pb_meta$subtype == st]
  st_meta <- pb_meta[st_samples, ]
  st_counts <- pb_mat[, st_samples, drop = FALSE]

  n_wt <- sum(st_meta$condition == REFERENCE_CONDITION)
  n_dko <- sum(st_meta$condition == TEST_CONDITION)
  cat(sprintf("  Samples: %d %s, %d %s\n", n_wt, REFERENCE_CONDITION, n_dko, TEST_CONDITION))

  if (n_wt < 2 || n_dko < 2) {
    cat("  SKIPPED: need at least 2 samples per condition.\n")
    next
  }

  # Keep genes expressed in >= MIN_GENE_EXPR_FRAC of this subtype's cells
  st_cells <- colnames(seu)[seu$subtype == st]
  gene_det_frac <- Matrix::rowSums(counts_mat[, st_cells, drop = FALSE] > 0) / length(st_cells)
  genes_pass_frac <- names(gene_det_frac[gene_det_frac >= MIN_GENE_EXPR_FRAC])
  st_counts <- st_counts[intersect(rownames(st_counts), genes_pass_frac), , drop = FALSE]
  cat(sprintf("  Genes expressed in >=%.0f%% of cells: %d\n",
              MIN_GENE_EXPR_FRAC * 100, nrow(st_counts)))

  # Require detection in a minimum number of pseudobulk samples
  keep_genes <- rowSums(st_counts >= MIN_COUNTS_PER_SAMPLE) >= MIN_SAMPLES_EXPRESSED
  st_counts <- st_counts[keep_genes, , drop = FALSE]
  cat(sprintf("  Genes after sample-level filter (>=%d counts in >=%d samples): %d\n",
              MIN_COUNTS_PER_SAMPLE, MIN_SAMPLES_EXPRESSED, nrow(st_counts)))

  if (nrow(st_counts) < 10) {
    cat("  SKIPPED: too few genes after filtering.\n")
    next
  }

  st_meta$condition <- factor(st_meta$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  dds <- DESeqDataSetFromMatrix(countData = round(st_counts), colData = st_meta,
                                design = ~ condition)

  dds <- tryCatch({
    DESeq(dds, quiet = TRUE)
  }, error = function(e) {
    cat(sprintf("  ERROR in DESeq: %s\n", e$message))
    return(NULL)
  })
  if (is.null(dds)) next

  # Cook's distance outlier filtering is disabled: with 3-4 samples per group
  # DESeq2 cannot replace outliers, so flagged genes would otherwise get NA.
  res <- results(dds, contrast = c("condition", TEST_CONDITION, REFERENCE_CONDITION),
                 cooksCutoff = FALSE)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  res_df$subtype <- st
  res_df <- res_df[order(res_df$pvalue), ]

  res_df$significance <- "NS"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange > LFC_THRESHOLD] <- label_up
  res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                      res_df$log2FoldChange < -LFC_THRESHOLD] <- label_down

  all_results[[st]] <- res_df

  n_up <- sum(res_df$significance == label_up)
  n_down <- sum(res_df$significance == label_down)
  n_tested <- sum(!is.na(res_df$padj))
  cat(sprintf("  DEGs (FDR < %.2f, |LFC| > %.2f): %d up, %d down (of %d tested)\n",
              FDR_THRESHOLD, LFC_THRESHOLD, n_up, n_down, n_tested))

  summary_stats <- rbind(summary_stats, data.frame(
    subtype = st,
    n_samples_WT = n_wt,
    n_samples_DKO = n_dko,
    n_cells = length(st_cells),
    n_genes_tested = n_tested,
    n_up_DKO = n_up,
    n_down_DKO = n_down,
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# 4. Save results
# =============================================================================
combined_res <- bind_rows(all_results)
write.csv(combined_res, file.path(OUTPUT_DIR, "DEG_all_subtypes.csv"), row.names = FALSE)
write.csv(summary_stats, file.path(OUTPUT_DIR, "DEG_summary.csv"), row.names = FALSE)
cat("\n\n== DEG Summary ==\n")
print(summary_stats)

st_dir <- file.path(OUTPUT_DIR, "per_subtype")
dir.create(st_dir, recursive = TRUE, showWarnings = FALSE)
for (st in names(all_results)) {
  write.csv(all_results[[st]], file.path(st_dir, paste0("DEG_", st, ".csv")), row.names = FALSE)
}

# =============================================================================
# 5. Volcano plots
# =============================================================================
cat("\nGenerating volcano plots...\n")
volcano_dir <- file.path(OUTPUT_DIR, "volcano_plots")
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

plot_list <- list()
for (st in names(all_results)) {
  res_df <- all_results[[st]]
  res_df <- res_df[!is.na(res_df$padj) & is.finite(res_df$padj), ]
  if (nrow(res_df) == 0) next

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
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 15, color = "black") +
    labs(title = paste(TEST_CONDITION, "vs", REFERENCE_CONDITION, "in CD8+ T", st),
         x = paste0("Log2 Fold Change (", TEST_CONDITION, " / ", REFERENCE_CONDITION, ")"),
         y = "-log10(adjusted p-value)",
         color = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  plot_list[[st]] <- p
  ggsave(file.path(volcano_dir, paste0("volcano_", st, ".pdf")), p, width = 4, height = 3.5)
}

if (length(plot_list) > 0) {
  ggsave(file.path(volcano_dir, "volcano_all_subtypes.pdf"),
         wrap_plots(plot_list, ncol = 2), width = 9, height = 4 * ceiling(length(plot_list) / 2))
}

# =============================================================================
# 6. Summary bar plot: number of DEGs per subtype
# =============================================================================
if (nrow(summary_stats) > 0) {
  deg_long <- summary_stats %>%
    select(subtype, n_up_DKO, n_down_DKO) %>%
    pivot_longer(cols = c(n_up_DKO, n_down_DKO),
                 names_to = "direction", values_to = "n_DEGs") %>%
    mutate(direction = ifelse(direction == "n_up_DKO", label_up, label_down),
           n_DEGs_signed = ifelse(direction == label_up, n_DEGs, -n_DEGs))

  p_bar <- ggplot(deg_long, aes(x = subtype, y = n_DEGs_signed, fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = setNames(c("#D73027", "#4575B4"), c(label_up, label_down))) +
    geom_hline(yintercept = 0, color = "black") +
    coord_flip() +
    labs(title = sprintf("DEGs per CD8+ T subtype (FDR < %.2f, |LFC| > %.2f)",
                         FDR_THRESHOLD, LFC_THRESHOLD),
         x = NULL, y = "Number of DEGs", fill = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ggsave(file.path(OUTPUT_DIR, "DEG_count_per_subtype.pdf"), p_bar, width = 7, height = 4)
}

# =============================================================================
# 7. Heatmap of DEGs across subtypes (signed -log10 p-value)
# =============================================================================
sig_genes <- unique(combined_res$gene[combined_res$significance != "NS"])

if (length(sig_genes) >= 3 && length(all_results) >= 2) {
  heatmap_df <- combined_res %>%
    filter(gene %in% sig_genes) %>%
    mutate(signed_logp = ifelse(significance != "NS",
                                sign(log2FoldChange) * -log10(pvalue), 0)) %>%
    select(gene, subtype, signed_logp) %>%
    pivot_wider(names_from = subtype, values_from = signed_logp) %>%
    tibble::column_to_rownames("gene")
  heatmap_df[is.na(heatmap_df)] <- 0

  cap_val <- quantile(abs(as.matrix(heatmap_df)[as.matrix(heatmap_df) != 0]), 0.99)
  heatmap_df[heatmap_df > cap_val] <- cap_val
  heatmap_df[heatmap_df < -cap_val] <- -cap_val

  pdf(file.path(OUTPUT_DIR, "DEG_heatmap_subtypes.pdf"),
      width = max(5, ncol(heatmap_df) * 0.6 + 3),
      height = max(6, nrow(heatmap_df) * 0.25 + 2))
  pheatmap(as.matrix(heatmap_df),
           color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
           breaks = seq(-cap_val, cap_val, length.out = 101),
           cluster_rows = TRUE, cluster_cols = FALSE,
           main = paste0("CD8+ T subtype DEGs, signed -log10(p) (",
                         TEST_CONDITION, " vs ", REFERENCE_CONDITION, ")"),
           fontsize_row = 9, fontsize_col = 10, treeheight_row = 0)
  dev.off()
}

# =============================================================================
# 8. MA plots
# =============================================================================
ma_dir <- file.path(OUTPUT_DIR, "MA_plots")
dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)

for (st in names(all_results)) {
  res_df <- all_results[[st]]
  res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$baseMean), ]
  if (nrow(res_df) == 0) next

  p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significance)) +
    geom_point(alpha = 0.5, size = 0.8) +
    scale_color_manual(values = setNames(c("#D73027", "#4575B4", "grey60"),
                                         c(label_up, label_down, "NS"))) +
    geom_hline(yintercept = 0, color = "black") +
    geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
               linetype = "dashed", color = "grey40") +
    labs(title = paste("MA plot:", TEST_CONDITION, "vs", REFERENCE_CONDITION, "in CD8+ T", st),
         x = "log10(Mean Expression + 1)", y = "Log2 Fold Change", color = NULL) +
    theme_publication() +
    theme(legend.position = "bottom")

  ggsave(file.path(ma_dir, paste0("MA_", st, ".pdf")), p_ma, width = 6, height = 4.5)
}

# =============================================================================
# 9. Pseudobulk QC: cells per pseudobulk sample
# =============================================================================
p_qc <- ggplot(pb_meta, aes(x = reorder(group, -n_cells), y = n_cells, fill = condition)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CONDITION_COLORS) +
  coord_flip() +
  labs(title = "Cells per pseudobulk sample", x = NULL, y = "Number of cells") +
  theme_publication() +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(OUTPUT_DIR, "pseudobulk_cell_counts.pdf"), p_qc,
       width = 7, height = max(4, nrow(pb_meta) * 0.2 + 2))

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("CD8+ T subtype DEG analysis complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("Clustering: %s\n", TCELL_CLUSTER_COL))
cat(sprintf("Subtypes tested: %d\n", nrow(summary_stats)))
cat(sprintf("Total significant DEGs: %d up, %d down\n",
            sum(summary_stats$n_up_DKO), sum(summary_stats$n_down_DKO)))
