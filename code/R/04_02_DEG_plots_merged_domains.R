# 04_02_DEG_DESeq_merged_domains.R
# Plotting DEG results from merged spatial domains (DKO vs WT).
# Reads saved DEG tables from 04_DEG_DESeq_merged_domains.R.

library(ggplot2)
library(ggrepel)
library(pheatmap)
library(dplyr)
library(tidyr)

# =============================================================================
# Parameters
# =============================================================================
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "04_DEG_merged_domains")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# 1. Load DEG results
# =============================================================================
combined_res <- read.csv(file.path(OUTPUT_DIR, "DEG_all_merged_domains.csv"),
                         stringsAsFactors = FALSE)
summary_stats <- read.csv(file.path(OUTPUT_DIR, "DEG_summary.csv"),
                           stringsAsFactors = FALSE)

cat(sprintf("Loaded DEG results: %d rows across %d domains\n",
            nrow(combined_res), length(unique(combined_res$domain))))

# Build per-domain list for iteration
all_results <- split(combined_res, combined_res$domain)

# =============================================================================
# 2. Volcano plots
# =============================================================================
cat("Generating volcano plots...\n")
volcano_dir <- file.path(OUTPUT_DIR, "volcano_plots")
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)

for (dom in names(all_results)) {
  res_df <- all_results[[dom]]
  res_df <- res_df[!is.na(res_df$padj) & is.finite(res_df$padj), ]

  if (nrow(res_df) == 0) next

  dom_safe <- gsub("[^A-Za-z0-9_]", "_", dom)
  dom_name <- DOMAIN_NAMES[dom]

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
    labs(title = paste0("DKO vs WT in Domain ", dom, " (", gsub("_", " ", dom_name), ")"),
         x = "Log2 Fold Change (DKO / WT)",
         y = "-log10(adjusted p-value)",
         color = NULL) +
    theme_publication(base_size = 8) +
    theme(legend.position = "bottom")

  ggsave(file.path(volcano_dir, paste0("volcano_merged_domain_", dom_safe, "_", dom_name, ".pdf")),
         p, width = 4, height = 3)
}

# =============================================================================
# 3. Summary bar plot: number of DEGs per merged domain
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
# 4. Heatmap of top DEGs across merged domains (signed -log10 p-value)
# =============================================================================
cat("Generating cross-domain heatmap (signed -log10 p-value)...\n")

sig_genes <- combined_res %>%
  filter(significance != "NS") %>%
  pull(gene) %>%
  unique()

if (length(sig_genes) > 0) {
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

  heatmap_df[is.na(heatmap_df)] <- 0

  cap_val <- quantile(abs(as.matrix(heatmap_df)[as.matrix(heatmap_df) != 0]), 0.99)
  heatmap_df[heatmap_df > cap_val] <- cap_val
  heatmap_df[heatmap_df < -cap_val] <- -cap_val

  if (nrow(heatmap_df) >= 3 && ncol(heatmap_df) >= 2) {
    pdf(file.path(OUTPUT_DIR, "DEG_heatmap_top_genes.pdf"),
        width = max(4, ncol(heatmap_df) * 0.3 + 3),
        height = max(8, nrow(heatmap_df) * 0.25 + 2))
    pheatmap(as.matrix(heatmap_df),
             color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
             breaks = seq(-cap_val, cap_val, length.out = 101),
             cluster_rows = TRUE,
             cluster_cols = FALSE,
             treeheight_row = 0,
             main = "DEGs (DKO vs WT) by Domain",
             fontsize_row = 10,
             fontsize_col = 10)
    dev.off()
  }
}

# =============================================================================
# 5. MA plots
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
# Done
# =============================================================================
cat("\n========================================\n")
cat("DEG plotting for merged spatial domains complete!\n")
cat("========================================\n")
cat(sprintf("Output directory: %s\n", OUTPUT_DIR))
cat("\nPlots generated:\n")
cat("  - volcano_plots/ (per-domain volcano PDFs)\n")
cat("  - MA_plots/ (per-domain MA PDFs)\n")
cat("  - DEG_heatmap_top_genes.pdf\n")
cat("  - DEG_count_per_merged_domain.pdf\n")