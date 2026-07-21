# 13_DE_pathways_immuneXdomains.R
# Differential pathway activity analysis for immune cell types across spatial domains.
# For each immune cell type x domain combination, run PanomiR differential pathway
# analysis (DKO vs WT), produce bar plots and gene heatmaps for significant pathways.

rm(list = ls())

library(PanomiR)
library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(forcats)
library(pheatmap)
library(stringr)

# =============================================================================
# Parameters
# =============================================================================
PROJECT_DIR <- "/data/work/Pourya/mmr_spatial"
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

BASE_OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "13_DE_pathways_immuneXdomains")
dir.create(BASE_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Pathway gene membership
pathways2 <- readRDS(file.path(PROJECT_DIR, "input/MSigDBPathGeneTab2024.RDS"))
pathway_displayNames <- read.csv(file.path(PROJECT_DIR, "input/MSigDB_display_names_v2024.csv"))
displayNames <- pathway_displayNames[, c("standard_name", "display_name")]

# Cell types of interest (immune)
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

# Domain definitions
DOMAIN_COL <- "banksy_domain_merged"
DOMAIN_NAMES <- c(
    "2" = "Tumor_Core",
    "3" = "Immune_Engulfing",
    "4" = "Stroma"
)

# QC thresholds
MIN_TOTAL_COUNTS <- 50
MIN_GENES_DETECTED <- 10
MIN_CELLS_PER_SAMPLE <- 10
MIN_PATH_SIZE <- 10

# Helper: filesystem-safe name
safe_name <- function(ct) {
    gsub("[^A-Za-z0-9_]", "_", gsub("\\+", "pos", ct))
}

# =============================================================================
# Helper: Convert mouse gene symbols to human ENSEMBL IDs
# =============================================================================
convert_mouse_to_human_ensembl <- function(mouse_symbols) {
    human_symbols <- toupper(mouse_symbols)
    mapping <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = human_symbols,
        keytype = "SYMBOL",
        columns = "ENSEMBL"
    )
    mapping <- mapping[!is.na(mapping$ENSEMBL), ]
    mapping <- mapping[!duplicated(mapping$SYMBOL), ]
    result <- setNames(mapping$ENSEMBL, mapping$SYMBOL)
    return(result)
}

# =============================================================================
# Build ENSEMBL -> mouse symbol mapping for heatmaps
# =============================================================================
ensembl_to_mouse <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(pathways2$ENSEMBL),
    keytype = "ENSEMBL",
    columns = "SYMBOL"
)
ensembl_to_mouse <- ensembl_to_mouse[!is.na(ensembl_to_mouse$SYMBOL), ]
ensembl_to_mouse <- ensembl_to_mouse[!duplicated(ensembl_to_mouse$ENSEMBL), ]
ensembl_to_mouse$mouse_symbol <- str_to_title(ensembl_to_mouse$SYMBOL)

# =============================================================================
# 1. Load Seurat object
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells,", nrow(seu), "genes\n")

# Recode condition
seu$condition <- gsub("^KO$", "DKO", seu$condition)
seu$condition <- factor(seu$condition, levels = c("WT", "DKO"))
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))
seu$domain <- as.character(seu@meta.data[[DOMAIN_COL]])

# =============================================================================
# 2. Loop over immune cell types
# =============================================================================
for (CELLTYPE in CELLTYPES_OF_INTEREST) {

    ct_safe <- safe_name(CELLTYPE)
    OUTPUT_DIR <- file.path(BASE_OUTPUT_DIR, ct_safe)
    dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

    cat("\n\n################################################################\n")
    cat(sprintf("## Processing: %s\n", CELLTYPE))
    cat("################################################################\n")

    # =========================================================================
    # 2a. Subset to cell type and filter
    # =========================================================================
    seu_ct <- subset(seu, subset = Tier1_celltype == CELLTYPE)
    cat(sprintf("  %s cells: %d\n", CELLTYPE, ncol(seu_ct)))

    # QC filter
    keep_cells <- seu_ct$total_counts >= MIN_TOTAL_COUNTS &
                  seu_ct$n_genes_by_counts >= MIN_GENES_DETECTED
    seu_ct <- subset(seu_ct, cells = colnames(seu_ct)[keep_cells])
    cat(sprintf("  After QC: %d cells\n", ncol(seu_ct)))

    # Keep only relevant domains
    seu_ct <- subset(seu_ct, subset = domain %in% names(DOMAIN_NAMES))
    cat(sprintf("  After domain filter: %d cells\n", ncol(seu_ct)))

    if (ncol(seu_ct) < MIN_CELLS_PER_SAMPLE * 2) {
        cat(sprintf("  SKIPPED %s: too few cells (%d).\n", CELLTYPE, ncol(seu_ct)))
        next
    }

    # =========================================================================
    # 2b. Pseudobulk aggregation (per sample x domain)
    # =========================================================================
    cat("  Aggregating pseudobulk profiles...\n")

    counts_mat <- GetAssayData(seu_ct, assay = "RNA", layer = "counts")
    seu_ct$pseudobulk_group <- paste(seu_ct$sample, seu_ct$domain, sep = "__")
    groups <- unique(seu_ct$pseudobulk_group)
    group_counts <- table(seu_ct$pseudobulk_group)

    pb_list <- list()
    for (g in groups) {
        cells_in_group <- colnames(seu_ct)[seu_ct$pseudobulk_group == g]
        if (length(cells_in_group) >= MIN_CELLS_PER_SAMPLE) {
            sub_mat <- counts_mat[, cells_in_group, drop = FALSE]
            pb_list[[g]] <- Matrix::rowSums(sub_mat)
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
    pb_meta$condition <- ifelse(grepl("^WT", pb_meta$sample), "WT", "DKO")
    rownames(pb_meta) <- pb_meta$group

    cat(sprintf("  Pseudobulk matrix: %d genes x %d samples\n", nrow(pb_mat), ncol(pb_mat)))

    # =========================================================================
    # 2c. Run DE pathway analysis per domain
    # =========================================================================
    for (dom_id in names(DOMAIN_NAMES)) {
        dom_name <- DOMAIN_NAMES[dom_id]
        cat(sprintf("\n  --- %s / %s ---\n", CELLTYPE, dom_name))

        # Subset to this domain
        dom_samples <- pb_meta$group[pb_meta$domain == dom_id]
        dom_meta <- pb_meta[dom_samples, ]

        # Check we have both conditions
        n_wt <- sum(dom_meta$condition == "WT")
        n_dko <- sum(dom_meta$condition == "DKO")
        if (n_wt < 2 || n_dko < 2) {
            cat(sprintf("  SKIPPED: not enough replicates (WT=%d, DKO=%d)\n", n_wt, n_dko))
            next
        }

        # Extract counts for this domain
        dom_counts <- as.data.frame(as.matrix(pb_mat[, dom_samples, drop = FALSE]))

        # Filter low-expressed genes
        gene_detected <- rowSums(dom_counts > 0) >= 2
        dom_counts <- dom_counts[gene_detected, , drop = FALSE]

        cat(sprintf("  Genes after filtering: %d\n", nrow(dom_counts)))

        if (nrow(dom_counts) < 50) {
            cat("  SKIPPED: too few genes after filtering.\n")
            next
        }

        # Convert mouse gene symbols to human ENSEMBL IDs
        mouse_symbols <- rownames(dom_counts)
        symbol_to_ensembl <- convert_mouse_to_human_ensembl(mouse_symbols)

        human_symbols_upper <- toupper(mouse_symbols)
        has_mapping <- human_symbols_upper %in% names(symbol_to_ensembl)
        dom_counts_ensembl <- dom_counts[has_mapping, , drop = FALSE]
        rownames(dom_counts_ensembl) <- symbol_to_ensembl[toupper(rownames(dom_counts_ensembl))]

        cat(sprintf("  Genes with ENSEMBL mapping: %d\n", nrow(dom_counts_ensembl)))

        # Build covariates
        covariates <- data.frame(
            sample = colnames(dom_counts_ensembl),
            condition = dom_meta[colnames(dom_counts_ensembl), "condition"],
            row.names = colnames(dom_counts_ensembl)
        )
        covariates$condition <- factor(covariates$condition, levels = c("WT", "DKO"))

        # Output directory for this cell type x domain
        dep_dir <- file.path(OUTPUT_DIR, dom_name)
        dir.create(dep_dir, recursive = TRUE, showWarnings = FALSE)

        # Run differential pathway analysis
        tryCatch({
            de.paths <- PanomiR::differentialPathwayAnalysis(
                geneCounts = dom_counts_ensembl,
                pathways = pathways2,
                covariates = covariates,
                condition = "condition",
                outDir = dep_dir,
                saveOutName = paste0(ct_safe, "_", dom_name, "_DEP_DKO_vs_WT.RDS"),
                minPathSize = MIN_PATH_SIZE,
                contrastConds = "conditionDKO-conditionWT"
            )

            # Save results
            temp <- de.paths$DEP
            temp <- tibble::rownames_to_column(temp, "standard_name")
            temp <- dplyr::left_join(temp, displayNames, by = "standard_name")
            temp <- temp[order(temp$P.Value), ]
            write.csv(temp,
                      file = file.path(dep_dir, paste0(ct_safe, "_", dom_name, "_diffPathways.csv")),
                      row.names = TRUE)

            # Save residuals
            saveRDS(de.paths$PathwayResiduals,
                    file.path(dep_dir, paste0(ct_safe, "_", dom_name, "_Residuals.RDS")))

            n_sig <- sum(temp$adj.P.Val < 0.1, na.rm = TRUE)
            cat(sprintf("  Significant pathways (adj.P < 0.1): %d\n", n_sig))

            # =================================================================
            # Bar plot of top significant pathways
            # =================================================================
            sig_paths <- temp[temp$adj.P.Val < 0.1, ]
            if (nrow(sig_paths) > 0) {
                sig_paths <- head(sig_paths, 30)

                sig_paths$label <- ifelse(
                    is.na(sig_paths$display_name) | sig_paths$display_name == "",
                    sig_paths$standard_name,
                    sig_paths$display_name
                )

                # Truncate long names
                if (any(str_length(sig_paths$label) > 70)) {
                    long_inds <- str_length(sig_paths$label) > 70
                    sig_paths[long_inds, ]$label <- str_sub(sig_paths[long_inds, ]$label, 1, 70)
                    sig_paths[long_inds, ]$label <- paste0(sig_paths[long_inds, ]$label, "*")
                }

                sig_paths$signed_neg_log_p <- sign(sig_paths$logFC) * (-log10(sig_paths$P.Value))
                sig_paths <- sig_paths[order(sig_paths$signed_neg_log_p), ]
                sig_paths$label <- factor(sig_paths$label, levels = sig_paths$label)

                p <- ggplot(sig_paths, aes(x = label, y = signed_neg_log_p, fill = logFC)) +
                    geom_bar(stat = "identity", width = 0.8) +
                    coord_flip() +
                    scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828",
                                         midpoint = 0, name = "log2FC") +
                    labs(
                        title = paste0(CELLTYPE, " - ", gsub("_", " ", dom_name)),
                        y = "Signed -log10(P-value)",
                        x = NULL
                    ) +
                    theme_bw() +
                    theme(
                        panel.grid.major.x = element_blank(),
                        panel.grid.minor.x = element_blank(),
                        panel.border = element_rect(colour = "grey80", fill = NA, linewidth = 0.5),
                        axis.title.x = element_text(size = 14),
                        axis.title.y = element_blank(),
                        axis.text.x = element_text(size = 12),
                        axis.text.y = element_text(size = 10),
                        legend.text = element_text(size = 11),
                        legend.position = "right",
                        legend.title = element_text(size = 11),
                        plot.title = element_text(size = 14, face = "bold")
                    )

                plot_height <- max(4, nrow(sig_paths) * 0.3 + 2)

                ggsave(
                    filename = file.path(dep_dir, paste0(ct_safe, "_", dom_name, "_barplot.pdf")),
                    plot = p, device = cairo_pdf,
                    width = 10, height = plot_height, limitsize = FALSE
                )
                ggsave(
                    filename = file.path(dep_dir, paste0(ct_safe, "_", dom_name, "_barplot.png")),
                    plot = p, width = 10, height = plot_height, dpi = 300, limitsize = FALSE
                )

                # =============================================================
                # Heatmaps for significant pathways
                # =============================================================
                heatmap_dir <- file.path(dep_dir, "heatmaps")
                dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

                # Use original mouse-symbol counts normalized to CPM
                lib_sizes <- colSums(dom_counts)
                cpm <- sweep(dom_counts, 2, lib_sizes, FUN = "/") * 1e6

                for (j in seq_len(nrow(sig_paths))) {
                    pathway_name <- sig_paths$standard_name[j]
                    pathway_label <- sig_paths$label[j]

                    # Get ENSEMBL IDs for this pathway
                    pathway_ensembl <- pathways2$ENSEMBL[pathways2$Pathway == pathway_name]
                    matched <- ensembl_to_mouse[ensembl_to_mouse$ENSEMBL %in% pathway_ensembl, ]
                    mouse_genes <- matched$mouse_symbol

                    # Filter to genes present in our dataset
                    genes_in_data <- mouse_genes[mouse_genes %in% rownames(cpm)]
                    if (length(genes_in_data) < 2) next

                    mat_log <- log2(as.matrix(cpm[genes_in_data, , drop = FALSE]) + 1)

                    col_annotation <- data.frame(
                        Condition = dom_meta[colnames(mat_log), "condition"],
                        row.names = colnames(mat_log)
                    )
                    ann_colors <- list(Condition = CONDITION_COLORS)

                    safe_filename <- substr(gsub("[^A-Za-z0-9_]", "_", pathway_name), 1, 80)
                    plot_title <- as.character(pathway_label)
                    if (nchar(plot_title) > 70) {
                        plot_title <- paste0(substr(plot_title, 1, 70), "*")
                    }

                    pheatmap(
                        mat_log,
                        main = plot_title,
                        fontsize_main = 8,
                        scale = "row",
                        cluster_rows = TRUE,
                        cluster_cols = FALSE,
                        treeheight_row = 0,
                        treeheight_col = 0,
                        annotation_col = col_annotation,
                        annotation_colors = ann_colors,
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        fontsize_row = max(4, min(8, 120 / length(genes_in_data))),
                        fontsize_col = 10,
                        filename = file.path(heatmap_dir, paste0(safe_filename, ".png")),
                        width = 6,
                        height = max(4, length(genes_in_data) * 0.2 + 2)
                    )
                }
                cat(sprintf("  Heatmaps saved to: %s\n", heatmap_dir))
            }
        }, error = function(e) {
            cat(sprintf("  ERROR in %s / %s: %s\n", CELLTYPE, dom_name, e$message))
        })
    }
}

cat("\n\nDone! All results saved to:", BASE_OUTPUT_DIR, "\n")