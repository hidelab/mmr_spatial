rm(list = ls())

library(PanomiR)
library(dplyr)
library(babelgene)  # proper mouse<->human orthology mapping (offline, HCOP-based)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

# Directory to the pathway gene membership
pathways2 <- readRDS(MSIGDB_PATHWAYS)
pathway_displayNames <- read.csv(MSIGDB_NAMES)
displayNames <- pathway_displayNames[, c("standard_name", "display_name")]

# List all pseudobulk domain profiles
pseudobulk_dir <- file.path(PROJECT_DIR, "output/R/04_DEG_merged_domains/pseudobulk_profiles")
pseudobulk_files <- list.files(pseudobulk_dir, pattern = "pseudobulk_domain_.*\\.csv", full.names = TRUE)

# Output directory for DE pathways
out.dir0 <- file.path(PROJECT_DIR, "output/R/12_DE_pathways/")
if (!dir.exists(out.dir0)) dir.create(out.dir0, recursive = TRUE)

# --- Helper: Convert mouse gene symbols to human ENSEMBL IDs via orthology ---
# NOTE: We map by true ortholog (babelgene / HCOP consensus of 12 databases),
# NOT by naive symbol case conversion. Uppercasing fails for the many genes
# whose human ortholog has a different symbol (e.g. Trp53 -> TP53, H2-K1 -> HLA-A)
# and can spuriously collide with unrelated human genes.
convert_mouse_to_human_ensembl <- function(mouse_symbols) {
    ortho <- babelgene::orthologs(
        genes   = unique(mouse_symbols),
        species = "mouse",
        human   = FALSE  # input genes are the non-human (mouse) species
    )

    # Keep only rows with a human ENSEMBL ID
    ortho <- ortho[!is.na(ortho$human_ensembl) & ortho$human_ensembl != "", ]

    # Resolve to a 1:1 mouse-symbol <-> human-ENSEMBL mapping, preferring the
    # ortholog with the strongest cross-database support:
    ortho <- ortho[order(-ortho$support_n), ]
    #  (a) one human ENSEMBL per mouse symbol
    ortho <- ortho[!duplicated(ortho$symbol), ]
    #  (b) one mouse symbol per human ENSEMBL (avoids duplicate row names downstream)
    ortho <- ortho[!duplicated(ortho$human_ensembl), ]

    setNames(ortho$human_ensembl, ortho$symbol)
}

# --- Loop over all domain pseudobulk profiles ---
for (pb_file in pseudobulk_files) {

    # Extract domain name from filename
    domain_name <- gsub("pseudobulk_domain_", "", basename(pb_file))
    domain_name <- gsub("\\.csv$", "", domain_name)
    domain_name <- gsub("^[0-9]_", "", domain_name)
    cat("\n========================================\n")
    cat("Processing domain:", domain_name, "\n")
    cat("========================================\n")

    # Read pseudobulk counts
    genes.counts2 <- read.csv(pb_file, row.names = 1)
    genes.counts2 <- as.data.frame(genes.counts2)

    cat("Original dimensions:", nrow(genes.counts2), "genes x", ncol(genes.counts2), "samples\n")

    # Convert mouse gene symbols (row names) to human ENSEMBL IDs via orthology
    mouse_symbols <- rownames(genes.counts2)
    symbol_to_ensembl <- convert_mouse_to_human_ensembl(mouse_symbols)

    # Filter to genes that have a valid ortholog mapping
    has_mapping <- rownames(genes.counts2) %in% names(symbol_to_ensembl)
    genes.counts2 <- genes.counts2[has_mapping, , drop = FALSE]

    # Replace row names with human ENSEMBL IDs
    rownames(genes.counts2) <- symbol_to_ensembl[rownames(genes.counts2)]

    cat("After ENSEMBL conversion:", nrow(genes.counts2), "genes mapped\n")

    # Build covariates from sample names
    covariates <- data.frame(
        sample = colnames(genes.counts2),
        condition = gsub("[0-9]*", "", colnames(genes.counts2)),
        row.names = colnames(genes.counts2)
    )

    # Ensure condition is a factor with WT as reference
    covariates$condition <- factor(covariates$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))

    cat("Samples:", paste(colnames(genes.counts2), collapse = ", "), "\n")
    cat("Conditions:", paste(levels(covariates$condition), collapse = " vs "), "\n")

    # Run differential pathway analysis
    de.paths <- PanomiR::differentialPathwayAnalysis(
        geneCounts = genes.counts2,
        pathways = pathways2,
        covariates = covariates,
        condition = "condition",
        outDir = out.dir0,
        saveOutName = paste0(domain_name, "_DEP_DKO_vs_WT.RDS"),
        minPathSize = MIN_PATH_SIZE,
        contrastConds = paste0("condition", TEST_CONDITION, "-condition", REFERENCE_CONDITION)
    )

    # Save DE pathway results
    temp <- de.paths$DEP
    temp <- tibble::rownames_to_column(temp, "standard_name")
    temp <- dplyr::left_join(temp, displayNames, by = "standard_name")
    temp <- temp[order(temp$P.Value), ]
    temp %>%
        write.csv(
            file = file.path(out.dir0, paste0(domain_name, "_diffPathways.csv")),
            row.names = TRUE
        )

    # Save pathway residuals
    de.paths$PathwayResiduals %>%
        saveRDS(file.path(out.dir0, paste0(domain_name, "_Residuals.RDS")))

    cat("Results saved for domain:", domain_name, "\n")
    cat("  Significant pathways (P < 0.05):", sum(temp$P.Value < 0.05, na.rm = TRUE), "\n")
}

cat("\n\nDone! All domain DE pathway results saved to:", out.dir0, "\n")

# ============================================================================
# Plot top significant pathways (adj.P.Val < 0.1) per domain
# ============================================================================
library(ggplot2)
library(forcats)
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# Read all DE pathway results and combine
all_dep_results <- list()
for (pb_file in pseudobulk_files) {
    domain_name <- gsub("pseudobulk_domain_", "", basename(pb_file))
    domain_name <- gsub("\\.csv$", "", domain_name)
    domain_name <- gsub("^[0-9]_","",domain_name)

    res_file <- file.path(out.dir0, paste0(domain_name, "_diffPathways.csv"))
    if (file.exists(res_file)) {
        dep <- read.csv(res_file, row.names = 1)
        dep$domain <- domain_name
        all_dep_results[[domain_name]] <- dep
    }
}

# Plot per domain
for (domain_name in names(all_dep_results)) {
    dep <- all_dep_results[[domain_name]]

    # Filter to adj.P.Val < 0.1
    sig_paths <- dep[dep$adj.P.Val < FDR_THRESHOLD, ]

    if (nrow(sig_paths) == 0) {
        cat(sprintf("No significant pathways (adj.P.Val < %.2f) for domain: %s\n", FDR_THRESHOLD, domain_name))
        next
    }

    # Take top pathways by P.Value (already ordered)
    sig_paths <- head(sig_paths, TOP_PATHWAYS_PLOT)

    # Use display_name as label; fall back to standard_name if missing
    sig_paths$label <- ifelse(
        is.na(sig_paths$display_name) | sig_paths$display_name == "",
        sig_paths$standard_name,
        sig_paths$display_name
    )

    # Truncate long pathway names
    if (any(stringr::str_length(sig_paths$label) > PATHWAY_NAME_TRUNC)) {
        long_inds <- stringr::str_length(sig_paths$label) > PATHWAY_NAME_TRUNC
        sig_paths[long_inds, ]$label <- stringr::str_sub(sig_paths[long_inds, ]$label, 1, PATHWAY_NAME_TRUNC)
        sig_paths[long_inds, ]$label <- paste0(sig_paths[long_inds, ]$label, "*")
    }

    # Compute signed -log10 p-value: sign(logFC) * -log10(P.Value)
    sig_paths$signed_neg_log_p <- sign(sig_paths$logFC) * (-log10(sig_paths$P.Value))

    # Order by signed -log10 p-value for the plot
    sig_paths <- sig_paths[order(sig_paths$signed_neg_log_p), ]
    sig_paths$label <- factor(sig_paths$label, levels = sig_paths$label)

    p <- ggplot(sig_paths, aes(x = label, y = signed_neg_log_p, fill = logFC)) +
        geom_bar(stat = "identity", width = 0.8) +
        coord_flip() +
        scale_fill_gradient2(low = "#1565C0", mid = "white", high = "#C62828", midpoint = 0,
                             name = "Effect\nsize") +
        labs(
            title = paste0("DE Pathways in ", gsub("_", " ", domain_name)),
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
            axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(size = 10, angle = 0, vjust = 0.5, hjust = 1),
            legend.text = element_text(size = 11),
            legend.position = "right",
            legend.title = element_text(size = 11),
            plot.title = element_text(size = 14, face = "bold")
        )

    plot_height <- max(4, nrow(sig_paths) * 0.2 + 2)

    ggsave(
        filename = file.path(out.dir0, paste0(domain_name, "_top_pathways_barplot.pdf")),
        plot = p,
        device = cairo_pdf,
        width = 10,
        height = plot_height,
        limitsize = FALSE
    )

    ggsave(
        filename = file.path(out.dir0, paste0(domain_name, "_top_pathways_barplot.png")),
        plot = p,
        width = 10,
        height = plot_height,
        dpi = 300,
        limitsize = FALSE
    )
    cat("Saved bar plot for domain:", domain_name, "\n")
}

# ============================================================================
# Heatmaps: per significant pathway, show gene expression across samples
# ============================================================================
library(pheatmap)

# Build human ENSEMBL -> mouse symbol mapping via orthology (reverse direction).
# Again, this is a true ortholog lookup rather than case conversion; note a human
# gene may have several mouse orthologs (e.g. HLA-A -> H2-K1/H2-Q7/...), so we keep
# all mouse orthologs and later intersect with the genes actually present in the data.
ensembl_to_mouse <- babelgene::orthologs(
    genes   = unique(pathways2$ENSEMBL),
    species = "mouse",
    human   = TRUE  # input genes are human (ENSEMBL IDs)
)
ensembl_to_mouse <- ensembl_to_mouse[!is.na(ensembl_to_mouse$symbol) & ensembl_to_mouse$symbol != "", ]
# Normalize column names used downstream: ENSEMBL (human) + mouse_symbol
ensembl_to_mouse$ENSEMBL <- ensembl_to_mouse$human_ensembl
ensembl_to_mouse$mouse_symbol <- ensembl_to_mouse$symbol

for (pb_file in pseudobulk_files) {

    domain_name <- gsub("pseudobulk_domain_", "", basename(pb_file))
    domain_name <- gsub("\\.csv$", "", domain_name)
    domain_name <- gsub("^[0-9]_","",domain_name)


    # Read original pseudobulk counts (mouse gene symbols as row names)
    raw_counts <- read.csv(pb_file, row.names = 1)
    raw_counts <- as.data.frame(raw_counts)

    # Normalize to CPM
    lib_sizes <- colSums(raw_counts)
    cpm <- sweep(raw_counts, 2, lib_sizes, FUN = "/") * 1e6

    # Read DE pathway results for this domain
    res_file <- file.path(out.dir0, paste0(domain_name, "_diffPathways.csv"))
    if (!file.exists(res_file)) next
    dep <- read.csv(res_file, row.names = 1)

    # Filter to significant pathways
    sig_paths <- dep[dep$adj.P.Val < FDR_THRESHOLD, ]
    if (nrow(sig_paths) == 0) next

    # Create output subdirectory for heatmaps
    heatmap_dir <- file.path(out.dir0, "heatmaps", domain_name)
    if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir, recursive = TRUE)

    cat("\nGenerating heatmaps for domain:", domain_name,
        "(", nrow(sig_paths), "significant pathways)\n")

    for (i in seq_len(nrow(sig_paths))) {
        pathway_name <- sig_paths$standard_name[i]
        pathway_label <- ifelse(
            is.na(sig_paths$display_name[i]) | sig_paths$display_name[i] == "",
            pathway_name,
            sig_paths$display_name[i]
        )

        # Get ENSEMBL IDs for this pathway
        pathway_ensembl <- pathways2$ENSEMBL[pathways2$Pathway == pathway_name]

        # Map to mouse symbols
        matched <- ensembl_to_mouse[ensembl_to_mouse$ENSEMBL %in% pathway_ensembl, ]
        mouse_genes <- matched$mouse_symbol

        # Filter to genes present in our dataset
        genes_in_data <- mouse_genes[mouse_genes %in% rownames(cpm)]

        if (length(genes_in_data) < 2) next

        # Extract CPM matrix for these genes
        mat <- as.matrix(cpm[genes_in_data, , drop = FALSE])

        # Log-transform for better visualization
        mat_log <- log2(mat + 1)

        # Column annotation for condition
        col_annotation <- data.frame(
            Condition = gsub("[0-9]*", "", colnames(mat_log)),
            row.names = colnames(mat_log)
        )

        ann_colors <- list(Condition = CONDITION_COLORS)

        # Truncate title for filename
        safe_filename <- gsub("[^A-Za-z0-9_]", "_", pathway_name)
        safe_filename <- substr(safe_filename, 1, 80)

        # Truncate label for plot title
        plot_title <- pathway_label
        if (nchar(plot_title) > PATHWAY_NAME_TRUNC) {
            plot_title <- paste0(substr(plot_title, 1, PATHWAY_NAME_TRUNC), "*")
        }

        # Generate heatmap
        pheatmap(
            mat_log,
            main = plot_title,
            fontsize = 6,
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
            fontsize_col = 6,
            filename = file.path(heatmap_dir, paste0(safe_filename, ".png")),
            width = 4,
            height = max(4, length(genes_in_data) * 0.2 + 2)
        )
    }
    cat("Heatmaps saved to:", heatmap_dir, "\n")
}
