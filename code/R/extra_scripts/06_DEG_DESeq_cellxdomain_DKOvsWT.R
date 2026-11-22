# 06_DEG_DESeq_cellxdomain_DKOvsWT.R
# Differential gene expression between DKO and WT within each spatial domain
# using pseudobulk aggregation + DESeq2.
# Comparisons: (a) DKO vs WT in Tumor Core
#              (b) DKO vs WT in Immune Engulfing
#              (c) DKO vs WT in Stroma
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
# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

BASE_OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "06_DEG_DKOvsWT_per_domain")
dir.create(BASE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# Cell types of interest (immune subtypes for domain DE)
CELLTYPES_OF_INTEREST <- IMMUNE_CELLTYPES_DE

# Domain column
DOMAIN_COL <- COL_DOMAIN

# Domain names for this analysis (exclude Mammary Glands)
DOMAIN_NAMES <- c(
  "2" = "Tumor_Core",
  "3" = "Immune_Engulfing",
  "4" = "Stroma"
)

# =============================================================================
# 1. Load Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
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

  # Keep only relevant domains
  seu_ct <- subset(seu_ct, subset = domain %in% names(DOMAIN_NAMES))
  cat(sprintf("  After domain filter (domains 2,3,4): %d cells\n", ncol(seu_ct)))

  if (ncol(seu_ct) < MIN_CELLS_PER_SAMPLE * 2) {
    cat(sprintf("  SKIPPED %s: too few cells after filtering (%d).\n", CELLTYPE, ncol(seu_ct)))
    next
  }

  # Report cells per domain x sample x condition
  cat("\nCells per domain x sample:\n")
  print(table(seu_ct$sample, seu_ct$domain))
  cat("\nCells per domain x condition:\n")
  print(table(seu_ct$condition, seu_ct$domain))
  cat("\n")

  # ===========================================================================
  # 3. Pseudobulk aggregation (per sample x domain)
  # ===========================================================================
  cat("Aggregating pseudobulk profiles...\n")

  counts_mat <- GetAssayData(seu_ct, assay = "RNA", layer = "counts")

  seu_ct$pseudobulk_group <- paste(seu_ct$sample, seu_ct$domain, sep = "__")
  group_counts <- table(seu_ct$pseudobulk_group)
  cat(sprintf("Total pseudobulk groups: %d\n", length(group_counts)))

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

  pb_meta <- data.frame(
    group = names(pb_list),
    sample = sapply(strsplit(names(pb_list), "__"), `[`, 1),
    domain = sapply(strsplit(names(pb_list), "__"), `[`, 2),
    n_cells = as.numeric(group_counts[names(pb_list)]),
    stringsAsFactors = FALSE
  )
  pb_meta$domain_name <- DOMAIN_NAMES[pb_meta$domain]
  pb_meta$condition <- ifelse(grepl("^WT", pb_meta$sample), REFERENCE_CONDITION, TEST_CONDITION)
  rownames(pb_meta) <- pb_meta$group

  cat(sprintf("Pseudobulk matrix: %d genes x %d samples (after min %d cells filter)\n",
              nrow(pb_mat), ncol(pb_mat), MIN_CELLS_PER_SAMPLE))
  cat("\nPseudobulk samples per domain x condition:\n")
  print(table(pb_meta$domain_name, pb_meta$condition))
  cat("\n")

  # ===========================================================================
  # 4. Run DESeq2: DKO vs WT within each domain
  # ===========================================================================
  all_results <- list()
  summary_stats <- data.frame()

  for (dom_id in names(DOMAIN_NAMES)) {
    dom_name <- DOMAIN_NAMES[dom_id]
    cat(sprintf("\n=== %s (domain %s): %s vs %s ===\n", dom_name, dom_id, TEST_CONDITION, REFERENCE_CONDITION))

    # Subset to this domain
    dom_samples <- pb_meta$group[pb_meta$domain == dom_id]
    dom_meta <- pb_meta[dom_samples, ]
    dom_counts <- pb_mat[, dom_samples, drop = FALSE]

    # Check replicates per condition
    n_test <- sum(dom_meta$condition == TEST_CONDITION)
    n_ref <- sum(dom_meta$condition == REFERENCE_CONDITION)
    cat(sprintf("  Samples: %d %s, %d %s\n", n_ref, REFERENCE_CONDITION, n_test, TEST_CONDITION))

    if (n_test < 2 || n_ref < 2) {
      cat("  SKIPPED: need at least 2 samples per condition.\n")
      next
    }

    # Gene filtering: expressed in >= MIN_GENE_EXPR_FRAC of cells in this domain
    dom_cells <- colnames(seu_ct)[seu_ct$domain == dom_id]
    dom_cell_counts <- counts_mat[, dom_cells, drop = FALSE]
    gene_det_frac <- Matrix::rowSums(dom_cell_counts > 0) / length(dom_cells)
    genes_pass_frac <- names(gene_det_frac[gene_det_frac >= MIN_GENE_EXPR_FRAC])
    dom_counts <- dom_counts[intersect(rownames(dom_counts), genes_pass_frac), , drop = FALSE]
    cat(sprintf("  Genes expressed in >=%.0f%% of cells: %d\n",
                MIN_GENE_EXPR_FRAC * 100, nrow(dom_counts)))

    # Sample-level filter
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
      countData = round(dom_counts),
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
    res_df$domain <- dom_name
    res_df <- res_df[order(res_df$pvalue), ]

    # Classify significance
    res_df$significance <- "NS"
    res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                        res_df$log2FoldChange > LFC_THRESHOLD] <- paste0("Up in ", TEST_CONDITION)
    res_df$significance[!is.na(res_df$padj) & res_df$padj < FDR_THRESHOLD &
                        res_df$log2FoldChange < -LFC_THRESHOLD] <- paste0("Up in ", REFERENCE_CONDITION)

    all_results[[dom_name]] <- res_df

    # Summary
    n_up <- sum(res_df$significance == paste0("Up in ", TEST_CONDITION))
    n_down <- sum(res_df$significance == paste0("Up in ", REFERENCE_CONDITION))
    n_tested <- sum(!is.na(res_df$padj))
    cat(sprintf("  DEGs (FDR < %.2f, |LFC| > %.2f): %d up in %s, %d up in %s (of %d tested)\n",
                FDR_THRESHOLD, LFC_THRESHOLD, n_up, TEST_CONDITION, n_down, REFERENCE_CONDITION, n_tested))

    # Cell counts per condition in this domain
    n_cells_test <- sum(seu_ct$domain == dom_id & seu_ct$condition == TEST_CONDITION)
    n_cells_ref <- sum(seu_ct$domain == dom_id & seu_ct$condition == REFERENCE_CONDITION)

    summary_stats <- rbind(summary_stats, data.frame(
      domain = dom_name,
      n_cells_WT = n_cells_ref,
      n_cells_KO = n_cells_test,
      n_samples_WT = n_ref,
      n_samples_KO = n_test,
      n_genes_tested = n_tested,
      n_up_KO = n_up,
      n_up_WT = n_down,
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
  write.csv(combined_res, file.path(OUTPUT_DIR, paste0("DEG_", ct_safe, "_DKOvsWT_all_domains.csv")),
            row.names = FALSE)

  write.csv(summary_stats, file.path(OUTPUT_DIR, paste0("DEG_", ct_safe, "_DKOvsWT_summary.csv")),
            row.names = FALSE)
  cat("\n\n== DEG Summary ==\n")
  print(summary_stats)

  # Per-domain result files
  dom_dir <- file.path(OUTPUT_DIR, "per_domain")
  dir.create(dom_dir, recursive = TRUE, showWarnings = FALSE)
  for (dom_name in names(all_results)) {
    write.csv(all_results[[dom_name]],
              file.path(dom_dir, paste0("DEG_", ct_safe, "_DKOvsWT_", dom_name, ".csv")),
              row.names = FALSE)
  }

  # ===========================================================================
  # 6. Volcano plots
  # ===========================================================================
  cat("\nGenerating volcano plots...\n")
  volcano_dir <- file.path(OUTPUT_DIR, "volcano_plots")
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

  for (dom_name in names(all_results)) {
    res_df <- all_results[[dom_name]]
    res_df <- res_df[!is.na(res_df$padj) & is.finite(res_df$padj), ]

    if (nrow(res_df) == 0) next

    top_genes <- res_df %>%
      filter(significance != "NS") %>%
      arrange(padj) %>%
      head(VOLCANO_TOP_N)

    sig_colors <- setNames(
      c("#D73027", "#4575B4", "grey60"),
      c(paste0("Up in ", TEST_CONDITION), paste0("Up in ", REFERENCE_CONDITION), "NS")
    )

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
      labs(title = paste0(CELLTYPE, " in ", gsub("_", " ", dom_name), ": ", TEST_CONDITION, " vs ", REFERENCE_CONDITION),
           x = paste0("Log2 Fold Change (", TEST_CONDITION, " / ", REFERENCE_CONDITION, ")"),
           y = "-log10(adjusted p-value)",
           color = NULL) +
      theme_publication() +
      theme(legend.position = "bottom")

    ggsave(file.path(volcano_dir, paste0("volcano_", ct_safe, "_DKOvsWT_", dom_name, ".pdf")),
           p, width = 7, height = 6)
  }

  # ===========================================================================
  # 7. Summary bar plot
  # ===========================================================================
  if (nrow(summary_stats) > 0) {
    deg_long <- summary_stats %>%
      select(domain, n_up_KO, n_up_WT) %>%
      pivot_longer(cols = c(n_up_KO, n_up_WT),
                   names_to = "direction", values_to = "n_DEGs") %>%
      mutate(direction = ifelse(direction == "n_up_KO",
                                paste0("Up in ", TEST_CONDITION),
                                paste0("Up in ", REFERENCE_CONDITION)),
             n_DEGs_signed = ifelse(grepl(TEST_CONDITION, direction), n_DEGs, -n_DEGs))

    p_bar <- ggplot(deg_long, aes(x = domain, y = n_DEGs_signed, fill = direction)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = setNames(c("#D73027", "#4575B4"),
                                          c(paste0("Up in ", TEST_CONDITION),
                                            paste0("Up in ", REFERENCE_CONDITION)))) +
      geom_hline(yintercept = 0, color = "black") +
      coord_flip() +
      labs(title = sprintf("%s DEGs: %s vs %s per domain (FDR < %.2f, |LFC| > %.2f)",
                           CELLTYPE, TEST_CONDITION, REFERENCE_CONDITION, FDR_THRESHOLD, LFC_THRESHOLD),
           x = NULL, y = "Number of DEGs", fill = NULL) +
      theme_publication() +
      theme(legend.position = "bottom")

    ggsave(file.path(OUTPUT_DIR, "DEG_count_per_domain.pdf"), p_bar, width = 8, height = 5)
  }

  # ===========================================================================
  # 8. Heatmap of top DEGs across domains (signed -log10 p-value)
  # ===========================================================================
  cat("Generating cross-domain heatmap (signed -log10 p-value)...\n")

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
      select(gene, domain, signed_logp) %>%
      pivot_wider(names_from = domain, values_from = signed_logp) %>%
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
               main = paste0(CELLTYPE, " — Significant DEGs signed -log10(p) (", TEST_CONDITION, " vs ", REFERENCE_CONDITION, ") per domain"),
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

  for (dom_name in names(all_results)) {
    res_df <- all_results[[dom_name]]
    res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$baseMean), ]

    if (nrow(res_df) == 0) next

    sig_colors <- setNames(
      c("#D73027", "#4575B4", "grey60"),
      c(paste0("Up in ", TEST_CONDITION), paste0("Up in ", REFERENCE_CONDITION), "NS")
    )

    p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = significance)) +
      geom_point(alpha = 0.5, size = 0.8) +
      scale_color_manual(values = sig_colors) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      geom_hline(yintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD),
                 linetype = "dashed", color = "grey40") +
      labs(title = paste0(CELLTYPE, " MA plot: ", TEST_CONDITION, " vs ", REFERENCE_CONDITION, " in ", gsub("_", " ", dom_name)),
           x = "log10(Mean Expression + 1)",
           y = "Log2 Fold Change",
           color = NULL) +
      theme_publication() +
      theme(legend.position = "bottom")

    ggsave(file.path(ma_dir, paste0("MA_", ct_safe, "_DKOvsWT_", dom_name, ".pdf")), p_ma, width = 7, height = 5)
  }

  cat(sprintf("\n  Done: %s -> %s\n", CELLTYPE, OUTPUT_DIR))

} # end cell type loop

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("DEG analysis: DKO vs WT per domain complete!\n")
cat("========================================\n")
cat(sprintf("Base output directory: %s\n", BASE_OUTPUT_DIR))
cat(sprintf("Cell types processed: %s\n", paste(CELLTYPES_OF_INTEREST, collapse = ", ")))
cat(sprintf("\nParameters:\n"))
cat(sprintf("  REFERENCE_CONDITION = %s\n", REFERENCE_CONDITION))
cat(sprintf("  TEST_CONDITION = %s\n", TEST_CONDITION))
cat(sprintf("  MIN_TOTAL_COUNTS = %d\n", MIN_TOTAL_COUNTS))
cat(sprintf("  MIN_GENES_DETECTED = %d\n", MIN_GENES_DETECTED))
cat(sprintf("  MIN_CELLS_PER_SAMPLE = %d\n", MIN_CELLS_PER_SAMPLE))
cat(sprintf("  FDR_THRESHOLD = %.2f\n", FDR_THRESHOLD))
cat(sprintf("  LFC_THRESHOLD = %.2f\n", LFC_THRESHOLD))
cat(sprintf("\nPer-celltype subfolders:\n"))
for (ct in CELLTYPES_OF_INTEREST) {
  cat(sprintf("  - %s/\n", safe_name(ct)))
}