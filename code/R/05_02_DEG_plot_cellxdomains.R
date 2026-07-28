# 05_02_DEG_DESeq_cellxdomains.R
# Plotting DEG results from pairwise domain comparisons (per cell type).
# Reads saved DEG tables from 05_DEG_DESeq_cellxdomains.R.

library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)

# =============================================================================
# Parameters
# =============================================================================
source("./code/R/00_parameters.R")

BASE_OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "05_DEG_across_domains")

source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

CELLTYPES_OF_INTEREST <- IMMUNE_CELLTYPES_DE

# Helper: make a filesystem-safe folder name from cell type
safe_name <- function(ct) {
  gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
}

# =============================================================================
# Loop over cell types
# =============================================================================
for (CELLTYPE in CELLTYPES_OF_INTEREST) {

  ct_safe <- safe_name(CELLTYPE)
  OUTPUT_DIR <- file.path(BASE_OUTPUT_DIR, ct_safe)

  # =========================================================================
  # 1. Load DEG results
  # =========================================================================
  combined_file <- file.path(OUTPUT_DIR, paste0("DEG_", ct_safe, "_all_comparisons.csv"))
  summary_file <- file.path(OUTPUT_DIR, paste0("DEG_", ct_safe, "_summary.csv"))

  if (!file.exists(combined_file)) {
    cat(sprintf("Skipping %s: no results file found.\n", CELLTYPE))
    next
  }

  combined_res <- read.csv(combined_file, stringsAsFactors = FALSE)
  summary_stats <- read.csv(summary_file, stringsAsFactors = FALSE)

  cat(sprintf("\n################################################################\n"))
  cat(sprintf("## Plotting: %s (%d rows, %d comparisons)\n",
              CELLTYPE, nrow(combined_res), length(unique(combined_res$comparison))))
  cat(sprintf("################################################################\n"))

  all_results <- split(combined_res, combined_res$comparison)

  # =========================================================================
  # 2. Volcano plots
  # =========================================================================
  cat("Generating volcano plots...\n")
  volcano_dir <- file.path(OUTPUT_DIR, "volcano_plots")
  dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

  for (comp_name in names(all_results)) {
    res_df <- all_results[[comp_name]]
    res_df <- res_df[!is.na(res_df$padj) & is.finite(res_df$padj), ]

    if (nrow(res_df) == 0) next

    parts <- strsplit(comp_name, "_vs_")[[1]]
    test_name <- gsub("_", " ", parts[1])
    ref_name <- gsub("_", " ", parts[2])

    top_genes <- res_df %>%
      filter(significance != "NS") %>%
      arrange(padj) %>%
      head(VOLCANO_TOP_N)

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

  # =========================================================================
  # 3. Summary bar plot
  # =========================================================================
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

  # =========================================================================
  # 4. Heatmap of top DEGs across comparisons (signed -log10 p-value)
  # =========================================================================
  cat("Generating cross-comparison heatmap (signed -log10 p-value)...\n")

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

    cap_val <- quantile(abs(as.matrix(heatmap_df)[as.matrix(heatmap_df) != 0]), 0.99)
    heatmap_df[heatmap_df > cap_val] <- cap_val
    heatmap_df[heatmap_df < -cap_val] <- -cap_val

    colnames(heatmap_df) <- gsub("_", " ", colnames(heatmap_df))

    if (nrow(heatmap_df) >= 3 && ncol(heatmap_df) >= 2) {
      pdf(file.path(OUTPUT_DIR, "DEG_heatmap_top_genes.pdf"),
          width = max(2, ncol(heatmap_df) * 1.2 + 1),
          height = max(8, nrow(heatmap_df) * 0.25 + 2))
      pheatmap(as.matrix(heatmap_df),
               color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
               breaks = seq(-cap_val, cap_val, length.out = 101),
               cluster_rows = TRUE,
               cluster_cols = FALSE,
               treeheight_row = 0,
               main = paste0(CELLTYPE, " DEGs across domains"),
               fontsize_row = 10,
               fontsize_col = 10)
      dev.off()
    }
  }

  # =========================================================================
  # 5. MA plots
  # =========================================================================
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

  cat(sprintf("\n  Done plotting: %s\n", CELLTYPE))

} # end cell type loop

# =============================================================================
# Done
# =============================================================================
cat("\n\n========================================\n")
cat("DEG plotting across spatial domains complete!\n")
cat("========================================\n")
cat(sprintf("Base output directory: %s\n", BASE_OUTPUT_DIR))
cat(sprintf("Cell types processed: %s\n", paste(CELLTYPES_OF_INTEREST, collapse = ", ")))
cat("\nPlots generated per cell type:\n")
cat("  - volcano_plots/ (per-comparison volcano PDFs)\n")
cat("  - MA_plots/ (per-comparison MA PDFs)\n")
cat("  - DEG_heatmap_top_genes.pdf\n")
cat("  - DEG_count_per_comparison.pdf\n")