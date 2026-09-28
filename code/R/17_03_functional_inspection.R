# 17_03_functional_inspection.R
# GSEA on a tumor-trajectory condition-test gene list (15_04 output).
#
# PURPOSE
#   Take a gam_condition_test_<model>.csv from 15_04 (default: CD8+ T,
#   pseudobulk_nb_fs) — or ANY CSV with a `gene` column plus `p_shape`/`p_level` —
#   rank its genes by significance, and run rank-based GSEA (clusterProfiler /
#   fgsea) against the project MSigDB collection to see which biological functions
#   sit at the TOP of the list.
#
# USAGE
#   Rscript code/R/17_03_functional_inspection.R [input_csv] [label]
#     input_csv : path to a gam_condition_test_*.csv    (default below)
#     label     : output-subfolder name                 (default derived from path)
#   Or edit DEFAULT_CSV and source interactively.
#
# WHAT IT DOES
#   Ranks genes by BOTH p_shape (WT->DKO trajectory-SHAPE difference — the novel
#   spatial result) and p_level (overall level offset ~ bulk DE), each in its own
#   gsea_<metric>/ subfolder, over ALL 3,646 pathways in the project MSigDB DB.
#
# TWO METHODOLOGICAL NOTES (both intrinsic to this input / panel)
#   * UNSIGNED metric. The condition test reports two-sided p-values with no effect
#     direction, so genes are ranked by -log10(p): a ONE-SIDED significance GSEA
#     (enrichment among the most-significant genes). NES > 0 = a function whose
#     members concentrate at the significant top of the list. This is the honest
#     ranking given the input — we cannot orient sets up/down without a signed stat.
#   * UNIVERSE = the tested panel. This is a targeted ~680-gene immune panel, so a
#     genome-wide ORA universe would be invalid. Rank-based GSEA is self-contained:
#     the ranked list IS the universe, and gene-set sizes are counted after
#     intersecting with it — the correct choice here. (GSEA-only by design.)

# --- Thread limits (set BEFORE BLAS/library load) ----------------------------
N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS,
           NUMEXPR_NUM_THREADS = N_THREADS)

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(babelgene)   # offline (HCOP-based) mouse<->human orthology, as in 12/13
})

source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# =============================================================================
# Arguments & parameters
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

# Default list: CD8+ T pseudobulk_nb_fs (the count-based defensible-shape spec).
DEFAULT_CSV <- file.path(PROJECT_DIR, "output", "R", "15_tumor_trajectory_working",
                         "CD8pos_T", "gam_condition_test_pseudobulk_fs.csv")
INPUT_CSV <- if (length(args) >= 1 && nzchar(args[1])) args[1] else DEFAULT_CSV
if (!file.exists(INPUT_CSV))
  stop(sprintf("Input CSV not found: %s\n(Pass a path as the first argument.)", INPUT_CSV))

# Output-subfolder label. Derived as "<celltype>__<model>" from a canonical
# .../<celltype>/gam_condition_test_<model>.csv path; override via arg 2.
derive_label <- function(path) {
  fn    <- sub("\\.csv$", "", basename(path))
  model <- sub("^gam_condition_test_", "", fn)
  ct    <- basename(dirname(path))
  gsub("[^A-Za-z0-9_]+", "_", paste(ct, model, sep = "__"))
}
LABEL <- if (length(args) >= 2 && nzchar(args[2])) args[2] else derive_label(INPUT_CSV)

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "17_functional_inspection", LABEL)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Rank by BOTH shape and level (separate runs). Only columns present in the CSV
# are used; a missing metric is skipped with a message.
RANK_METRICS    <- c("p_shape", "p_level")
# GSEA gene-set size bounds. minGSSize is counted AFTER intersecting with the
# ranked panel genes, so on a ~680-gene panel keep it small (5) or almost nothing
# is testable; MIN_PATH_SIZE (=10, for the genome-scale 12/13 workflow) is too
# strict here.
GSEA_MIN_GSSIZE <- 5
GSEA_MAX_GSSIZE <- 500
GSEA_PADJ_SIG   <- 0.25              # BH-FDR cutoff flagged as "significant" / plotted
TOP_N_PLOT      <- TOP_PATHWAYS_PLOT # 30 (00_parameters.R)
set.seed(42)                         # fgsea permutation reproducibility

# =============================================================================
# Load list + gene sets
# =============================================================================
tab <- read.csv(INPUT_CSV, stringsAsFactors = FALSE)
if (!"gene" %in% colnames(tab)) stop("Input CSV has no `gene` column.")
cat(sprintf("Input : %s\n  %d genes; columns: %s\n  label -> %s\n",
            INPUT_CSV, nrow(tab), paste(colnames(tab), collapse = ", "), LABEL))

# --- Mouse symbol -> human ENSEMBL via true orthology (same as 12_pathway_activity)
# The project MSigDB DB is keyed on human ENSEMBL, so map by ortholog (NOT case
# conversion) to avoid the Trp53->TP53 / H2-K1->HLA-A failures.
convert_mouse_to_human_ensembl <- function(mouse_symbols) {
  ortho <- babelgene::orthologs(genes = unique(mouse_symbols),
                                species = "mouse", human = FALSE)
  ortho <- ortho[!is.na(ortho$human_ensembl) & ortho$human_ensembl != "", ]
  ortho <- ortho[order(-ortho$support_n), ]   # strongest cross-DB support first
  ortho <- ortho[!duplicated(ortho$symbol), ]        # 1 human ENSEMBL per mouse symbol
  ortho <- ortho[!duplicated(ortho$human_ensembl), ] # 1 mouse symbol per human ENSEMBL
  setNames(ortho$human_ensembl, ortho$symbol)
}
sym2ens <- convert_mouse_to_human_ensembl(tab$gene)
ens2sym <- setNames(names(sym2ens), sym2ens)   # reverse map for leading-edge labels
cat(sprintf("Orthology: %d / %d panel genes mapped to a human ENSEMBL ID\n",
            sum(tab$gene %in% names(sym2ens)), nrow(tab)))

# --- MSigDB TERM2GENE (human ENSEMBL) + TERM2NAME (display names) -------------
pathways2 <- readRDS(MSIGDB_PATHWAYS)
term2gene <- unique(pathways2[, c("Pathway", "ENSEMBL")])
disp      <- read.csv(MSIGDB_NAMES, stringsAsFactors = FALSE)
term2name <- unique(disp[, c("standard_name", "display_name")])
term2name$display_name <- ifelse(is.na(term2name$display_name) |
                                   term2name$display_name == "",
                                 term2name$standard_name, term2name$display_name)
cat(sprintf("Gene sets: %d pathways in the project MSigDB collection\n",
            length(unique(term2gene$Pathway))))

# Record the run parameters (self-documenting output folder).
writeLines(c(
  "17_03 functional inspection (GSEA) — run parameters",
  strrep("=", 60),
  sprintf("run_time        = %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("input_csv       = %s", INPUT_CSV),
  sprintf("label           = %s", LABEL),
  sprintf("rank_metrics    = %s", paste(RANK_METRICS, collapse = ", ")),
  sprintf("n_genes_input   = %d", nrow(tab)),
  sprintf("n_genes_mapped  = %d", sum(tab$gene %in% names(sym2ens))),
  sprintf("n_pathways_db   = %d", length(unique(term2gene$Pathway))),
  sprintf("GSEA_MIN_GSSIZE = %d", GSEA_MIN_GSSIZE),
  sprintf("GSEA_MAX_GSSIZE = %d", GSEA_MAX_GSSIZE),
  sprintf("GSEA_PADJ_SIG   = %.2f", GSEA_PADJ_SIG),
  sprintf("metric          = -log10(p)  (unsigned, one-sided GSEA)"),
  sprintf("clusterProfiler = %s", as.character(utils::packageVersion("clusterProfiler")))
), file.path(OUTPUT_DIR, "run_parameters.txt"))

# =============================================================================
# GSEA per ranking metric
# =============================================================================
run_gsea <- function(metric) {
  if (!metric %in% colnames(tab)) {
    cat(sprintf("  [%s] column absent in input; skipping\n", metric)); return(invisible(NULL))
  }
  sub_dir <- file.path(OUTPUT_DIR, paste0("gsea_", metric))
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)

  # Ranked list: -log10(p), larger = more significant = top. Collapse duplicate
  # ENSEMBL (many->one orthology) to the MOST significant member.
  d     <- tab[!is.na(tab[[metric]]) & tab$gene %in% names(sym2ens), , drop = FALSE]
  if (nrow(d) < GSEA_MIN_GSSIZE) {
    cat(sprintf("  [%s] too few mapped genes (%d); skipping\n", metric, nrow(d)))
    return(invisible(NULL))
  }
  score <- -log10(pmax(d[[metric]], .Machine$double.eps))
  agg   <- tapply(score, sym2ens[d$gene], max)
  # tapply() returns a 1-D array; clusterProfiler::GSEA needs a plain named numeric.
  geneList <- sort(setNames(as.numeric(agg), names(agg)), decreasing = TRUE)
  cat(sprintf("  [%s] ranked geneList: %d genes  (-log10 p range %.2f .. %.2f)\n",
              metric, length(geneList), min(geneList), max(geneList)))

  # Save the exact ranked list used (for reproducibility / re-plotting).
  write.csv(data.frame(ensembl = names(geneList),
                       mouse_symbol = ens2sym[names(geneList)],
                       neglog10_p = as.numeric(geneList)),
            file.path(sub_dir, sprintf("ranked_genelist_%s.csv", metric)),
            row.names = FALSE)

  res <- tryCatch(
    clusterProfiler::GSEA(geneList = geneList, TERM2GENE = term2gene,
                          TERM2NAME = term2name, minGSSize = GSEA_MIN_GSSIZE,
                          maxGSSize = GSEA_MAX_GSSIZE, pvalueCutoff = 1,
                          eps = 0, seed = TRUE, verbose = FALSE),
    error = function(e) { cat("  GSEA error:", conditionMessage(e), "\n"); NULL })
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    cat(sprintf("  [%s] no gene sets testable at these size bounds\n", metric))
    return(invisible(NULL))
  }

  rdf <- as.data.frame(res)
  # Translate leading-edge (core_enrichment) ENSEMBL ids back to mouse symbols.
  rdf$core_enrichment_mouse <- vapply(strsplit(rdf$core_enrichment, "/"),
    function(ids) paste(stats::na.omit(ens2sym[ids]), collapse = "/"), character(1))
  write.csv(rdf, file.path(sub_dir, sprintf("gsea_%s_all_pathways.csv", metric)),
            row.names = FALSE)
  saveRDS(res, file.path(sub_dir, sprintf("gsea_%s_result.rds", metric)))
  cat(sprintf("  [%s] %d pathways tested; %d at BH-FDR < %.2f (NES>0: %d)\n",
              metric, nrow(rdf), sum(rdf$p.adjust < GSEA_PADJ_SIG),
              GSEA_PADJ_SIG, sum(rdf$p.adjust < GSEA_PADJ_SIG & rdf$NES > 0)))

  # --- Bar plot: top pathways enriched among the most-significant genes --------
  # One-sided metric -> interpret NES>0 (concentration at the significant top).
  plt <- rdf[rdf$NES > 0, , drop = FALSE]
  plt <- plt[order(plt$p.adjust, plt$pvalue), , drop = FALSE]
  plt <- head(plt, TOP_N_PLOT)
  if (nrow(plt) >= 1) {
    lbl <- ifelse(nchar(plt$Description) > PATHWAY_NAME_TRUNC,
                  paste0(substr(plt$Description, 1, PATHWAY_NAME_TRUNC), "*"),
                  plt$Description)
    plt$label       <- factor(lbl, levels = rev(lbl))
    plt$neglog_padj <- -log10(pmax(plt$p.adjust, 1e-300))
    p <- ggplot(plt, aes(x = label, y = NES, fill = neglog_padj)) +
      geom_col(width = 0.75) + coord_flip() +
      scale_fill_viridis_c(option = "viridis", name = "-log10\nBH-FDR") +
      labs(title = sprintf("GSEA (%s) — %s", metric, LABEL),
           subtitle = sprintf("Top %d pathways enriched among most-significant genes (NES > 0)",
                              nrow(plt)),
           x = NULL, y = "Normalized Enrichment Score (NES)") +
      theme_publication()
    h <- max(4, nrow(plt) * 0.30 + 2)
    ggsave(file.path(sub_dir, sprintf("gsea_%s_top_pathways.pdf", metric)),
           p, width = 11, height = h, limitsize = FALSE)
    ggsave(file.path(sub_dir, sprintf("gsea_%s_top_pathways.png", metric)),
           p, width = 11, height = h, dpi = 300, limitsize = FALSE)
  } else {
    cat(sprintf("  [%s] no NES>0 pathway to plot\n", metric))
  }

  # --- Classic GSEA running-score plot for the single top pathway --------------
  if (requireNamespace("enrichplot", quietly = TRUE) && nrow(plt) >= 1) {
    top_id <- rdf$ID[order(rdf$p.adjust, rdf$pvalue)][1]
    gp <- tryCatch(enrichplot::gseaplot2(res, geneSetID = top_id,
                                         title = plt$Description[1]),
                   error = function(e) NULL)
    if (!is.null(gp))
      ggsave(file.path(sub_dir, sprintf("gsea_%s_top1_runningscore.pdf", metric)),
             gp, width = 7, height = 5)
  }
  invisible(res)
}

for (m in RANK_METRICS) {
  cat(sprintf("\n=== GSEA ranked by %s ===\n", m))
  run_gsea(m)
}

cat(sprintf("\nDone. Outputs -> %s\n", OUTPUT_DIR))
