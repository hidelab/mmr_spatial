# ==============================================================================
# 12_pathway_activity.R — Differential pathway activity (DKO vs WT) per domain
# ==============================================================================
#   1. PanomiR differential pathway analysis on each domain's pseudobulk profile
#   2. EnrichmentMap-style pathway modules (pathways sharing genes) + diagnostics
#   3. Bar plots of significant pathways, grouped and labelled by module
#   4. Per-pathway gene heatmaps
# Helper functions live in the pathway-analysis section of 00_utils.R.
# ==============================================================================
rm(list = ls())

library(PanomiR)
library(dplyr)
library(babelgene)  # proper mouse<->human orthology mapping (offline, HCOP-based)
library(Matrix)
library(ggplot2)
library(ggraph)     # igraph is used via igraph:: (attaching it masks dplyr functions)
library(pheatmap)

# Centralized parameters (defines PROJECT_DIR and all shared constants), shared
# helpers and the project colour palette
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

FDR_THRESHOLD <- 0.1

# Pathway gene membership and display names
pathways2 <- readRDS(MSIGDB_PATHWAYS)
displayNames <- read.csv(MSIGDB_NAMES)[, c("standard_name", "display_name")]

# All pseudobulk domain profiles
pseudobulk_dir <- file.path(PROJECT_DIR, "output/R/04_DEG_merged_domains/pseudobulk_profiles")
pseudobulk_files <- list.files(pseudobulk_dir, pattern = "pseudobulk_domain_.*\\.csv", full.names = TRUE)

# Output directory for DE pathways (results for each domain go in their own subfolder)
out.dir0 <- file.path(PROJECT_DIR, "output/R/12_DE_pathways/")
if (!dir.exists(out.dir0)) dir.create(out.dir0, recursive = TRUE)

# Path of a per-domain output file: <out.dir0>/<domain>/<domain><suffix>
domain_file <- function(domain_name, suffix) file.path(out.dir0, domain_name, paste0(domain_name, suffix))

# ============================================================================
# 1. Differential pathway analysis (PanomiR) per domain
# ============================================================================
# Genes scored by PanomiR in each domain (reused to build pathway modules),
# and their human ENSEMBL -> mouse symbol mapping (for labelling module genes)
domain_genes <- list()
domain_gene_symbols <- list()

for (pb_file in pseudobulk_files) {
    domain_name <- domain_from_pseudobulk_file(pb_file)
    cat("\n========================================\n")
    cat("Processing domain:", domain_name, "\n")
    cat("========================================\n")

    # Per-domain output folder (must exist before PanomiR is called)
    domain_dir <- file.path(out.dir0, domain_name)
    if (!dir.exists(domain_dir)) dir.create(domain_dir, recursive = TRUE)

    # Pseudobulk counts; mouse symbols -> human ENSEMBL by orthology (unmapped genes dropped)
    genes.counts2 <- read.csv(pb_file, row.names = 1)
    cat("Original dimensions:", nrow(genes.counts2), "genes x", ncol(genes.counts2), "samples\n")
    symbol_to_ensembl <- convert_mouse_to_human_ensembl(rownames(genes.counts2))
    genes.counts2 <- genes.counts2[rownames(genes.counts2) %in% names(symbol_to_ensembl), , drop = FALSE]
    rownames(genes.counts2) <- symbol_to_ensembl[rownames(genes.counts2)]
    cat("After ENSEMBL conversion:", nrow(genes.counts2), "genes mapped\n")
    domain_genes[[domain_name]] <- rownames(genes.counts2)
    domain_gene_symbols[[domain_name]] <- setNames(names(symbol_to_ensembl), symbol_to_ensembl)

    # Covariates from sample names (e.g. "DKO2" -> "DKO"), with WT as reference
    covariates <- data.frame(
        sample = colnames(genes.counts2),
        condition = factor(gsub("[0-9]*", "", colnames(genes.counts2)),
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION)),
        row.names = colnames(genes.counts2)
    )
    cat("Samples:", paste(colnames(genes.counts2), collapse = ", "), "\n")
    cat("Conditions:", paste(levels(covariates$condition), collapse = " vs "), "\n")

    de.paths <- PanomiR::differentialPathwayAnalysis(
        geneCounts = genes.counts2,
        pathways = pathways2,
        covariates = covariates,
        condition = "condition",
        outDir = domain_dir,
        saveOutName = paste0(domain_name, "_DEP_DKO_vs_WT.RDS"),
        minPathSize = MIN_PATH_SIZE,
        contrastConds = paste0("condition", TEST_CONDITION, "-condition", REFERENCE_CONDITION)
    )

    # Save DE pathway results (with display names, sorted by P.Value) and residuals
    temp <- dplyr::left_join(tibble::rownames_to_column(de.paths$DEP, "standard_name"), displayNames,
                             by = "standard_name")
    temp <- temp[order(temp$P.Value), ]
    write.csv(temp, file = domain_file(domain_name, "_diffPathways.csv"), row.names = TRUE)
    saveRDS(de.paths$PathwayResiduals, domain_file(domain_name, "_Residuals.RDS"))

    cat("Results saved for domain:", domain_name, "\n")
    cat(sprintf("  Significant pathways (adj.P.Val < %.2f): %d\n", FDR_THRESHOLD,
                sum(temp$adj.P.Val < FDR_THRESHOLD, na.rm = TRUE)))
}

cat("\n\nDone! All domain DE pathway results saved to per-domain folders in:", out.dir0, "\n")

# ============================================================================
# 2. Pathway modules: group pathways that share genes (EnrichmentMap-style)
# ============================================================================
# Many MSigDB pathways overlap heavily on this gene panel (e.g. the Reactome
# DNA replication / repair family), which makes the bar plots hard to read.
# For each domain we build a binary pathway x gene matrix over the genes PanomiR
# scored (same >= MIN_PATH_SIZE filter, so the pathways match the tested set),
# compute pairwise gene-set similarity with matrix products, cluster pathways
# hierarchically and cut the tree at 1 - PATHWAY_SIM_THRESHOLD, and name each
# module after its medoid. Core genes are the genes shared by most of a module's
# significant pathways.
sim_label <- pathway_similarity_label(PATHWAY_SIM_METHOD)
module_rule <- pathway_module_rule(PATHWAY_SIM_METHOD, PATHWAY_SIM_THRESHOLD, PATHWAY_MODULE_LINKAGE)

for (domain_name in names(domain_genes)) {
    diag_dir <- file.path(out.dir0, domain_name, "module_diagnostics")
    if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
    diag_file <- function(suffix) file.path(diag_dir, paste0(domain_name, suffix))

    # Similarity over the genes PanomiR scored in this domain, then modules
    M <- build_pathway_gene_matrix(pathways2, domain_genes[[domain_name]], MIN_PATH_SIZE)
    inter <- pathway_intersections(M)
    sim <- pathway_similarity(inter, PATHWAY_SIM_METHOD, PATHWAY_SIM_COMBINED_K)
    clust <- cluster_pathway_modules(sim, PATHWAY_SIM_THRESHOLD, PATHWAY_MODULE_LINKAGE)
    module <- clust$module

    # Tables: pathways with modules and DE statistics, module genes (core genes), module summary
    dep <- read.csv(domain_file(domain_name, "_diffPathways.csv"), row.names = 1)
    modules_df <- module_pathway_table(M, inter, module, module_medoids(sim, module), dep)
    sig_df <- significant_pathways(modules_df, FDR_THRESHOLD)
    module_genes <- module_gene_table(M, module, sig_df$standard_name, domain_gene_symbols[[domain_name]],
                                      MODULE_CORE_GENE_FRAC, MODULE_CORE_GENE_MIN_N)
    module_summary <- module_summary_table(modules_df, sim, module, module_genes, FDR_THRESHOLD)
    write.csv(modules_df, domain_file(domain_name, "_pathway_modules.csv"), row.names = FALSE)
    write.csv(module_genes, domain_file(domain_name, "_module_genes.csv"), row.names = FALSE)
    write.csv(module_summary, domain_file(domain_name, "_module_summary.csv"), row.names = FALSE)
    saveRDS(
        list(incidence = M, intersection = inter, similarity = sim, hclust = clust$hclust,
             module = module, method = PATHWAY_SIM_METHOD, threshold = PATHWAY_SIM_THRESHOLD,
             linkage = PATHWAY_MODULE_LINKAGE),
        domain_file(domain_name, "_pathway_similarity.RDS")
    )
    cat(sprintf("\n%s: %d pathways -> %d modules (%d singletons)\n",
                domain_name, nrow(M), length(unique(module)), sum(table(module) == 1)))

    # Diagnostics: similarity distribution, threshold sweep, module sizes, similarity heatmap
    plot_similarity_distribution(sim, PATHWAY_SIM_THRESHOLD, sim_label, domain_name,
                                 diag_file("_similarity_distribution"))
    threshold_sweep <- module_threshold_sweep(sim, clust$hclust, seq(0.05, 0.95, by = 0.05))
    write.csv(threshold_sweep, diag_file("_threshold_sweep.csv"), row.names = FALSE)
    plot_threshold_sensitivity(threshold_sweep, PATHWAY_SIM_THRESHOLD, sim_label, module_rule, domain_name,
                               diag_file("_threshold_sensitivity"))
    plot_module_sizes(module, module_rule, domain_name, diag_file("_module_sizes"))
    direction <- setNames(
        ifelse(modules_df$adj.P.Val < FDR_THRESHOLD, ifelse(modules_df$logFC > 0, "Up", "Down"), "NS"),
        modules_df$standard_name
    )
    plot_similarity_heatmap_all(sim, clust$hclust, direction, sim_label, domain_name,
                                diag_file("_similarity_heatmap_all.png"))

    # Significant-pathway heatmap and module network need at least two significant pathways
    cat(sprintf("  %d significant pathways in %d modules\n", nrow(sig_df), length(unique(sig_df$module))))
    if (nrow(sig_df) < 2) next
    plot_similarity_heatmap_significant(sim, clust$hclust, module, sig_df, direction, sim_label, domain_name,
                                        diag_file("_similarity_heatmap_significant.png"))
    plot_module_network(sim, sig_df, PATHWAY_SIM_THRESHOLD, FDR_THRESHOLD, sim_label, domain_name,
                        diag_file("_module_network"))
}

# ============================================================================
# 3. Bar plots of significant pathways (adj.P.Val < FDR_THRESHOLD) per domain
# ============================================================================
# Gene-level DEG log2 fold changes (DKO vs WT) per domain from script 04, used to
# colour module genes. Domains are stored there by number; DOMAIN_NAMES maps them.
deg_all <- read.csv(file.path(PROJECT_DIR, "output", "R", "04_DEG_merged_domains", "DEG_all_merged_domains.csv"))
deg_lfc <- lapply(split(deg_all, DOMAIN_NAMES[as.character(deg_all$domain)]),
                  function(d) setNames(d$log2FoldChange, d$gene))

for (pb_file in pseudobulk_files) {
    domain_name <- domain_from_pseudobulk_file(pb_file)
    domain_label <- gsub("_", " ", domain_name)
    if (!file.exists(domain_file(domain_name, "_diffPathways.csv"))) next

    # DE results with module annotation and core-gene label text (plain and DEG-coloured)
    dep <- read.csv(domain_file(domain_name, "_diffPathways.csv"), row.names = 1)
    modules_df <- read.csv(domain_file(domain_name, "_pathway_modules.csv"))
    dep <- dplyr::left_join(dep, modules_df[, c("standard_name", "module", "module_size", "module_label")],
                            by = "standard_name")
    module_summary <- read.csv(domain_file(domain_name, "_module_summary.csv"))
    core_by_path <- module_summary$core_genes[match(dep$module, module_summary$module)]
    domain_lfc <- deg_lfc[[domain_name]]
    if (is.null(domain_lfc)) domain_lfc <- numeric(0)  # no DEG results: genes stay uncoloured
    dep$core_genes_text <- paste0("Common genes: ",format_core_genes(core_by_path))
    dep$core_genes_md <-  paste0("Common genes: ", format_core_genes(core_by_path, gene_lfc = domain_lfc))

    # Significant = FDR-adjusted p-value (adj.P.Val, BH) < FDR_THRESHOLD
    sig_paths <- significant_pathways(dep, FDR_THRESHOLD)
    if (nrow(sig_paths) == 0) {
        cat(sprintf("No significant pathways (adj.P.Val < %.2f) for domain: %s\n", FDR_THRESHOLD, domain_name))
        next
    }
    # sig_paths is ordered by P.Value, so head() on each direction picks the most significant
    up_sig   <- sig_paths[sig_paths$logFC > 0, ]
    down_sig <- sig_paths[sig_paths$logFC < 0, ]
    if (nrow(up_sig) + nrow(down_sig) == 0) {
        cat(sprintf("No up/down significant pathways for domain: %s\n", domain_name))
        next
    }

    # Top-pathway plots: top 10 up + top 10 down, and top TOP_PATHWAYS_PLOT per direction
    top_plots <- list(
        "_top_pathways_barplot"      = list(rbind(head(up_sig, 10), head(down_sig, 10)), "DE Pathways in "),
        "_top_up_pathways_barplot"   = list(head(up_sig, TOP_PATHWAYS_PLOT), "Upregulated Pathways in "),
        "_top_down_pathways_barplot" = list(head(down_sig, TOP_PATHWAYS_PLOT), "Downregulated Pathways in ")
    )
    for (suffix in names(top_plots)) {
        if (nrow(top_plots[[suffix]][[1]]) == 0) next
        plot_pathway_barplot(top_plots[[suffix]][[1]], paste0(top_plots[[suffix]][[2]], domain_label),
                             domain_file(domain_name, suffix), module_rule)
    }

    # Module plots, one per direction: significant pathways in modules with >= 2
    # significant pathways in that direction, labelled with core genes. Each comes in
    # two versions: plain gene names, and gene names coloured by DEG direction.
    module_dirs <- list(
        up   = list(up_sig,   "Upregulated Modules with 2+ Significant Pathways, "),
        down = list(down_sig, "Downregulated Modules with 2+ Significant Pathways, ")
    )
    for (dir in names(module_dirs)) {
        dir_paths <- module_dirs[[dir]][[1]]
        n_per_module <- table(dir_paths$module)
        dir_paths <- dir_paths[dir_paths$module %in% names(n_per_module)[n_per_module >= 3], ]
        if (nrow(dir_paths) == 0) next
        for (colored in c(FALSE)) {
            plot_pathway_barplot(
                dir_paths, paste0(module_dirs[[dir]][[2]], domain_label),
                domain_file(domain_name, paste0("_module_", dir, "_pathways_barplot", if (colored) "_colored" else "")),
                module_rule, show_core_genes = TRUE, color_genes = colored
            )
        }
    }
    cat("Saved bar plots for domain:", domain_name, "\n")
}

# ============================================================================
# 4. Heatmaps: per significant pathway, gene expression (log2 CPM) across samples
# ============================================================================
ensembl_to_mouse <- human_ensembl_to_mouse(pathways2$ENSEMBL)

for (pb_file in pseudobulk_files) {
    domain_name <- domain_from_pseudobulk_file(pb_file)
    if (!file.exists(domain_file(domain_name, "_diffPathways.csv"))) next
    sig_paths <- significant_pathways(read.csv(domain_file(domain_name, "_diffPathways.csv"), row.names = 1),
                                      FDR_THRESHOLD)
    if (nrow(sig_paths) == 0) next

    # CPM from the original pseudobulk counts (mouse gene symbols as row names)
    raw_counts <- read.csv(pb_file, row.names = 1)
    cpm <- sweep(raw_counts, 2, colSums(raw_counts), FUN = "/") * 1e6

    heatmap_dir <- file.path(out.dir0, domain_name, "heatmaps")
    if (!dir.exists(heatmap_dir)) dir.create(heatmap_dir, recursive = TRUE)
    cat("\nGenerating heatmaps for domain:", domain_name, "(", nrow(sig_paths), "significant pathways)\n")

    titles <- truncate_label(pathway_label(sig_paths$standard_name, sig_paths$display_name), PATHWAY_NAME_TRUNC)
    for (i in seq_len(nrow(sig_paths))) {
        plot_pathway_gene_heatmap(sig_paths$standard_name[i], titles[i], pathways2, ensembl_to_mouse, cpm,
                                  heatmap_dir, CONDITION_COLORS)
    }
    cat("Heatmaps saved to:", heatmap_dir, "\n")
}