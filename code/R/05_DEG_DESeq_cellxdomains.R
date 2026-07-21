# 05_DEG_DESeq_cellxdomains.R
# Differential gene expression across spatial domains using
# pseudobulk aggregation + DESeq2.
# Comparisons: (a) Tumor Core vs Immune Engulfing
#              (b) Immune Engulfing vs Stroma
#              (c) Tumor Core vs Stroma
# Pseudobulk is aggregated per sample (each sample = one replicate per domain).
# Runs for each cell type in CELLTYPES_OF_INTEREST, saving to per-celltype subfolders.

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
PROJECT_DIR <- "/data/work/Pourya/mmr_spatial"
BASE_OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "05_DEG_across_domains")
dir.create(BASE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# Cell-level QC filters
MIN_TOTAL_COUNTS <- 50
MIN_CELLS_PER_SAMPLE <- 10
MIN_GENES_DETECTED <- 10
FDR_THRESHOLD <- 0.1
LFC_THRESHOLD <- 0.25

# Gene-level filtering
MIN_GENE_EXPR_FRAC <- 0.05
MIN_COUNTS_PER_SAMPLE <- 5
MIN_SAMPLES_EXPRESSED <- 3

# --- Cell types of interest ---
CELLTYPES_OF_INTEREST <- c(
  "CD8+ T",
  "Treg",
  "DC (Ccr7+)",
  "pDC",
  "NK",
  "Macrophage (Cxcl16+)",
  "cDC1",
  "Neutrophil"
)

# Domain column and names
DOMAIN_COL <- "banksy_domain_merged"
DOMAIN_NAMES <- c(
  "2" = "Tumor_Core",
  "3" = "Immune_Engulfing",
  "4" = "Stroma"
)

# Pairwise comparisons: test vs reference
# Positive LFC = higher in test domain
COMPARISONS <- list(
  list(test = "Tumor_Core", ref = "Immune_Engulfing", test_id = "2", ref_id = "3"),
  list(test = "Immune_Engulfing", ref = "Stroma", test_id = "3", ref_id = "4"),
  list(test = "Tumor_Core", ref = "Stroma", test_id = "2", ref_id = "4")
)

# =============================================================================
# 1. Load Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", "DKO", seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

seu$domain <- as.character(seu@meta.data[[DOMAIN_COL]])

# =============================================================================
# 2. Loop over cell types of interest
# =============================================================================

# Helper: make a filesystem-safe folder name from cell type
safe_name <- function(ct) {
  gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
}

for (CELLTYPE in CELLTYPES_OF_INTEREST) {

  ct_safe <- safe_name(CELLTYPE)
  OUTPUT_DIR <- file.path(BASE_OUTPUT_DIR, ct_safe)
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

  cat("\n\n################################################################\n")
  cat(sprintf("## Processing: %s\n", CELLTYPE))
  cat(sprintf("## Output: %s\n", OUTPUT_DIR))
  cat("################################################################\n")

  # ===========================================================================
  # 2a. Subset to cell type and filter
  # ===========================================================================
  cat(sprintf("\nSubsetting to %s cells...\n", CELLTYPE))

  seu_ct <- subset(seu, subset = Tier1_celltype == CELLTYPE)
  cat(sprintf("  %s cells: %d\n", CELLTYPE, ncol(seu_ct)))

  # QC filter
  n_before <- ncol(seu_ct)
  keep_cells <- seu_ct$total_counts >= MIN_TOTAL_COUNTS &
                seu_ct$n_genes_by_counts >= MIN_GENES_DETECTED
  seu_ct <- subset(seu_ct, cells = colnames(seu_ct)[keep_cells])
  n_after <- ncol(seu_ct)
  cat(sprintf("  After QC filtering: %d -> %d cells (removed %d, %.1f%%)\n",
              n_before, n_after, n_before - n_after,
              100 * (n_before - n_after) / n_before))

  # Keep only relevant domains (Tumor Core, Immune Engulfing, Stroma)
  seu_ct <- subset(seu_ct, subset = domain %in% names(DOMAIN_NAMES))
  cat(sprintf("  After domain filter (domains 2,3,4): %d cells\n", ncol(seu_ct)))

  if (ncol(seu_ct) < MIN_CELLS_PER_SAMPLE * 2) {
    cat(sprintf("  SKIPPED %s: too few cells after filtering (%d).\n", CELLTYPE, ncol(seu_ct)))
    next
  }

  # Report cells per domain x sample
  cat("\nCells per domain x sample:\n")
  print(table(seu_ct$sample, seu_ct$domain))
  cat("\n")

  # ===========================================================================
  # 3. Pseudobulk aggregation (per sample x domain)
  # ===========================================================================
  cat("Aggregating pseudobulk profiles...\n")

  counts_mat <- GetAssayData(seu_ct, assay = "RNA", layer = "counts")

  # Group by sample + domain
  seu_ct$pseudobulk_group <- paste(seu_ct$sample, seu_ct$domain, sep = "__")
  group_counts <- table(seu_ct$pseudobulk_group)
  cat(sprintf("Total pseudobulk groups: %d\n", length(group_counts)))

  # Aggregate
  groups <- unique(seu_ct$pseudobulk_group)
  pb_list <- list()

  for (g in groups) {
    cells_in_group <- colnames(seu_ct)[seu_ct$pseudobulk_group == g]
    if (length(cells_in_group) >= MIN_CELLS_PER_SAMPLE) {
      sub_mat <- counts_mat[, cells_in_group, drop = FALSE]
      pb_list[[g]] <- Matrix::rowSums(sub_mat)
    } else {
      cat(sprintf("  Skipped %s: only %d cells\n", g, length(cells_in_group)))
    }
  }

  if (length(pb_list) < 4) {
    cat(sprintf("  SKIPPED %s: too few pseudobulk samples (%d).\n", CELLTYPE, length(pb_list)))
    next
  }

  pb_mat <- do.call(cbind, pb_list)
  colnames(pb_mat) <- names(pb_list)

  # Build metadata
  pb_meta <- data.frame(
    group = names(pb_list),
    sample = sapply(strsplit(names(pb_list), "__"), `[`, 1),
    domain = sapply(strsplit(names(pb_list), "__"), `[`, 2),
    n_cells = as.numeric(group_counts[names(pb_list)]),
    stringsAsFactors = FALSE
  )
  pb_meta$domain_name <- DOMAIN_NAMES[pb_meta$domain]
  pb_meta$condition <- ifelse(grepl("^WT", pb_meta$sample), "WT", "DKO")
  pb_meta$condition <- factor(pb_meta$condition, levels = c("WT", "DKO"))
  rownames(pb_meta) <- pb_meta$group

  cat(sprintf("Pseudobulk matrix: %d genes x %d samples (after min %d cells filter)\n",
              nrow(pb_mat), ncol(pb_mat), MIN_CELLS_PER_SAMPLE))
  cat("\nPseudobulk samples per domain:\n")
  print(table(pb_meta$domain_name))
  cat("\n")

  # ===========================================================================
  # 4. Run DESeq2 for each pairwise domain comparison
  # ===========================================================================
  all_results <- list()
  summary_stats <- data.frame()

  for (comp in COMPARISONS) {
    comp_name <- paste0(comp$test, "_vs_", comp$ref)
    cat(sprintf("\n=== %s vs %s ===\n", comp$test, comp$ref))

    # Subset to the two domains
    comp_samples <- pb_meta$group[pb_meta$domain %in% c(comp$test_id, comp$ref_id)]
    comp_meta <- pb_meta[comp_samples, ]
    comp_counts <- pb_mat[, comp_samples, drop = FALSE]

    # Check replicates
    n_test <- sum(comp_meta$domain == comp$test_id)
    n_ref <- sum(comp_meta$domain == comp$ref_id)
    cat(sprintf("  Samples: %d %s, %d %s\n", n_test, comp$test, n_ref, comp$ref))

    if (n_test < 2 || n_ref < 2) {
      cat("  SKIPPED: need at least 2 samples per domain.\n")
      next
    }

    # Gene filtering: expressed in >= MIN_GENE_EXPR_FRAC of cells in these domains
    comp_cells <- colnames(seu_ct)[seu_ct$domain %in% c(comp$test_id, comp$ref_id)]
    comp_cell_counts <- counts_mat[, comp_cells, drop = FALSE]
    gene_det_frac <- Matrix::rowSums(comp_cell_counts > 0) / length(comp_cells)
    genes_pass_frac <- names(gene_det_frac[gene_det_frac >= MIN_GENE_EXPR_FRAC])
    comp_counts <- comp_counts[intersect(rownames(comp_counts), genes_pass_frac), , drop = FALSE]
    cat(sprintf("  Genes expressed in >=%.0f%% of cells: %d\n",
                MIN_GENE_EXPR_FRAC * 100, nrow(comp_counts)))

    # Sample-level filter
    n_samples_detected <- rowSums(comp_counts >= MIN_COUNTS_PER_SAMPLE)
    keep_genes <- n_samples_detected >= MIN_SAMPLES_EXPRESSED
    comp_counts <- comp_counts[keep_genes, , drop = FALSE]
    cat(sprintf("  Genes after sample-level filter (>=%d counts in >=%d samples): %d\n",
                MIN_COUNTS_PER_SAMPLE, MIN_SAMPLES_EXPRESSED, nrow(comp_counts)))

    if (nrow(comp_counts) < 10) {
      cat("  SKIPPED: too few genes after filtering.\n")
      next
    }

    # DESeq2: domain as the factor, adjusting for condition (WT/DKO)
    comp_meta$domain_factor <- factor(comp_meta$domain_name,
                                      levels = c(comp$ref, comp$test))
    dds <- DESeqDataSetFromMatrix(
      countData = round(comp_counts),
      colData = comp_meta,
      design = ~ condition + domain_factor
    )

    dds <- tryCatch({
      DESeq(dds, quiet = TRUE)
    }, error = function(e) {
      cat(sprintf("  ERROR in DESeq: %s\n", e$message))
      return(NULL)
    })

    if (is.null(dds)) next

    res <- results(dds, contrast = c("domain_factor", comp$test, comp$ref))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$comparison <- comp_name
    res_df <- res_df[order(res_df$pvalue), ]

    # Classify significance
    res_df$significance <- "NS"
    res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                        res_df$log2FoldChange > LFC_THRESHOLD] <- paste0("Up in ", comp$test)
    res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                        res_df$log2FoldChange < -LFC_THRESHOLD] <- paste0("Up in ", comp$ref)

    all_results[[comp_name]] <- res_df

    # Summary
    n_up <- sum(grepl("^Up in", res_df$significance) & res_df$log2FoldChange > 0)
    n_down <- sum(grepl("^Up in", res_df$significance) & res_df$log2FoldChange < 0)
    n_tested <- sum(!is.na(res_df$padj))
    cat(sprintf("  DEGs (FDR < %.2f, |LFC| > %.2f): %d up in %s, %d up in %s (of %d tested)\n",
                FDR_THRESHOLD, LFC_THRESHOLD, n_up, comp$test, n_down, comp$ref, n_tested))

    # Cell counts per domain in this comparison
    n_cells_test <- sum(seu_ct$domain == comp$test_id)
    n_cells_ref <- sum(seu_ct$domain == comp$ref_id)

    summary_stats <- rbind(summary_stats, data.frame(
      comparison = comp_name,
      n_cells_test = n_cells_test,
      n_cells_ref = n_cells_ref,
      n_samples_test = n_test,
      n_samples_ref = n_ref,
      n_genes_tested = n_tested,
      n_up_test = n_up,
      n_up_ref = n_down,
      stringsAsFactors = FALSE
    ))
  }

  # ===========================================================================
  # 5. Save results
  # ===========================================================================
  if (length(all_results) == 0) {
    cat(sprintf("\n  No successful comparisons for %s. Skipping outputs.\n", CELLTYPE))
    next
  }

  combined_res <- bind_rows(all_results)
  write.csv(combined_res, file.path(OUTPUT_DIR, paste0("DEG_", ct_safe, "_all_comparisons.csv")),
            row.names = FALSE)

  write.csv(summary_stats, file.path(OUTPUT_DIR, paste0("DEG_", ct_safe, "_summary.csv")),
            row.names = FALSE)
  cat("\n\n== DEG Summary ==\n")
  print(summary_stats)

  # Per-comparison result files
  comp_dir <- file.path(OUTPUT_DIR, "per_comparison")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  for (comp_name in names(all_results)) {
    write.csv(all_results[[comp_name]],
              file.path(comp_dir, paste0("DEG_", ct_safe, "_", comp_name, ".csv")),
              row.names = FALSE)
  }

  # ===========================================================================
  # 6. Volcano plots
  # ===========================================================================
  cat("\nGenerating volcano plots...\n")
  volcano_dir <- file.path(OUTPUT_DIR, "volcano_plots")
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

  for (comp_name in names(all_results)) {
    res_df <- all_results[[comp_name]]
    res_df <- res_df[!is.na(res_df$padj) & is.finite(res_df$padj), ]

    if (nrow(res_df) == 0) next

    # Parse comparison name for labels
    parts <- strsplit(comp_name, "_vs_")[[1]]
    test_name <- gsub("_", " ", parts[1])
    ref_name <- gsub("_", " ", parts[2])

    # Top genes to label
    top_genes <- res_df %>%
      filter(significance != "NS") %>%
      arrange(padj) %>%
      head(20)

    # Color values
    sig_levels <- unique(res_df$significance)
    sig_colors <- c("NS" = "grey60")
    for (s in sig_levels[sig_levels != "NS"]) {
      if (grepl(parts[1], s)) sig_colors[s] <- "#D73027"
      else sig_colors[s] <- "#4575B4"
    }

    p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = sig_colors) +
      geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
                 linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(FDR_THRESHOLD),
                 linetype = "dashed", color = "grey40") +
      geom_text_repel(data = top_genes,
                      aes(label = gene),
                      size = 2.5, max.overlaps = 15,
                      color = "black") +
      labs(title = paste0(CELLTYPE, ": ", test_name, " vs ", ref_name),
           x = paste0("Log2 Fold Change (", test_name, " / ", ref_name, ")"),
           y = "-log10(adjusted p-value)",
           color = NULL) +
      theme_publication() +
      theme(legend.position = "bottom")

    ggsave(file.path(volcano_dir, paste0("volcano_", ct_safe, "_", comp_name, ".pdf")),
           p, width = 7, height = 6)
  }

  # ===========================================================================
  # 7. Summary bar plot
  # ===========================================================================
  if (nrow(summary_stats) > 0) {
    deg_long <- summary_stats %>%
      select(comparison, n_up_test, n_up_ref) %>%
      pivot_longer(cols = c(n_up_test, n_up_ref),
                   names_to = "direction", values_to = "n_DEGs") %>%
      mutate(direction = ifelse(direction == "n_up_test", "Up in test domain", "Up in ref domain"),
             n_DEGs_signed = ifelse(direction == "Up in test domain", n_DEGs, -n_DEGs))

    p_bar <- ggplot(deg_long, aes(x = comparison, y = n_DEGs_signed, fill = direction)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Up in test domain" = "#D73027", "Up in ref domain" = "#4575B4")) +
      geom_hline(yintercept = 0, color = "black") +
      coord_flip() +
      labs(title = sprintf("%s DEGs across domains (FDR < %.2f, |LFC| > %.2f)",
                           CELLTYPE, FDR_THRESHOLD, LFC_THRESHOLD),
           x = NULL, y = "Number of DEGs", fill = NULL) +
      theme_publication() +
      theme(legend.position = "bottom")

    ggsave(file.path(OUTPUT_DIR, "DEG_count_per_comparison.pdf"), p_bar, width = 8, height = 5)
  }

  # ===========================================================================
  # 8. Heatmap of top DEGs across comparisons (signed -log10 p-value)
  # ===========================================================================
  cat("Generating cross-comparison heatmap (signed -log10 p-value)...\n")

  # All significant DEGs (no slice_head limit)
  sig_genes <- combined_res %>%
    filter(significance != "NS") %>%
    pull(gene) %>%
    unique()

  if (length(sig_genes) > 0) {
    heatmap_df <- combined_res %>%
      filter(gene %in% sig_genes) %>%
      mutate(signed_logp = ifelse(
        significance != "NS",
        sign(log2FoldChange) * -log10(pvalue),
        0
      )) %>%
      select(gene, comparison, signed_logp) %>%
      pivot_wider(names_from = comparison, values_from = signed_logp) %>%
      tibble::column_to_rownames("gene")

    heatmap_df[is.na(heatmap_df)] <- 0

    # Cap at 99th percentile
    cap_val <- quantile(abs(as.matrix(heatmap_df)[as.matrix(heatmap_df) != 0]), 0.99)
    heatmap_df[heatmap_df > cap_val] <- cap_val
    heatmap_df[heatmap_df < -cap_val] <- -cap_val

    if (nrow(heatmap_df) >= 3 && ncol(heatmap_df) >= 2) {
      pdf(file.path(OUTPUT_DIR, "DEG_heatmap_top_genes.pdf"),
          width = max(8, ncol(heatmap_df) * 1.5 + 3),
          height = max(8, nrow(heatmap_df) * 0.25 + 2))
      pheatmap(as.matrix(heatmap_df),
               color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
               breaks = seq(-cap_val, cap_val, length.out = 101),
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               treeheight_row = 0,
               main = paste0(CELLTYPE, " — Significant DEGs signed -log10(p) across domain comparisons"),
               fontsize_row = 10,
               fontsize_col = 10)
      dev.off()
    }
  }

  # ===========================================================================
  # 9. MA plots
  # ===========================================================================
  cat("Generating MA plots...\n")
  ma_dir <- file.path(OUTPUT_DIR, "MA_plots")
  dir.create(ma_dir, recursive = TRUE, showWarnings = FALSE)

  for (comp_name in names(all_results)) {
    res_df <- all_results[[comp_name]]
    res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$baseMean), ]

    if (nrow(res_df) == 0) next

    parts <- strsplit(comp_name, "_vs_")[[1]]
    test_name <- gsub("_", " ", parts[1])
    ref_name <- gsub("_", " ", parts[2])

    # Color values
    sig_levels <- unique(res_df$significance)
    sig_colors <- c("NS" = "grey60")
    for (s in sig_levels[sig_levels != "NS"]) {
      if (grepl(parts[1], s)) sig_colors[s] <- "#D73027"
      else sig_colors[s] <- "#4575B4"
    }

    p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significance)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = sig_colors) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
                 linetype = "dashed", color = "grey40") +
      labs(title = paste0(CELLTYPE, " MA plot: ", test_name, " vs ", ref_name),
           x = "log10(Mean Expression + 1)",
           y = "Log2 Fold Change",
           color = NULL) +
      theme_publication() +
      theme(legend.position = "bottom")

    ggsave(file.path(ma_dir, paste0("MA_", ct_safe, "_", comp_name, ".pdf")), p_ma, width = 7, height = 5)
  }

  cat(sprintf("\n  Done: %s -> %s\n", CELLTYPE, OUTPUT_DIR))

} # end cell type loop

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("DEG analysis across spatial domains complete!\n")
cat("========================================\n")
cat(sprintf("Base output directory: %s\n", BASE_OUTPUT_DIR))
cat(sprintf("Cell types processed: %s\n", paste(CELLTYPES_OF_INTEREST, collapse = ", ")))
cat(sprintf("\nParameters:\n"))
cat(sprintf("  MIN_TOTAL_COUNTS = %d\n", MIN_TOTAL_COUNTS))
cat(sprintf("  MIN_GENES_DETECTED = %d\n", MIN_GENES_DETECTED))
cat(sprintf("  MIN_CELLS_PER_SAMPLE = %d\n", MIN_CELLS_PER_SAMPLE))
cat(sprintf("  FDR_THRESHOLD = %.2f\n", FDR_THRESHOLD))
cat(sprintf("  LFC_THRESHOLD = %.2f\n", LFC_THRESHOLD))
cat(sprintf("\nPer-celltype subfolders:\n"))
for (ct in CELLTYPES_OF_INTEREST) {
  cat(sprintf("  - %s/\n", safe_name(ct)))
}