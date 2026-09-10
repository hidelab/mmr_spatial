# 17_01_verify_DKO2_knockout.R
# Sanity check: did the Mlh1/Msh2 knockout "take" in every DKO sample, or is
# one sample (suspected DKO2) a failed KO that behaves like WT?
#
# Three levels of evidence:
#   (1) DIRECT   — Mlh1/Msh2 (KO'd genes) expression in TUMOR cells, per sample.
#                  Partners Msh6 (dimerizes with Msh2) / Pms2 (with Mlh1) shown too.
#   (2) GLOBAL   — per-sample tumor-cell pseudobulk: correlation heatmap +
#                  hierarchical clustering + PCA. Does DKO2 sit among the WTs?
#   (3) FUNCTIONAL (corroboration) — a small IFN/cytotoxic signature score in
#                  tumor cells per sample (circular w.r.t. project biology, so
#                  treated as supporting only).

N_THREADS <- 8
Sys.setenv(OMP_NUM_THREADS = N_THREADS, MKL_NUM_THREADS = N_THREADS,
           OPENBLAS_NUM_THREADS = N_THREADS, VECLIB_MAXIMUM_THREADS = N_THREADS)

library(Seurat)
library(ggplot2)
library(patchwork)

source(file.path("code", "R", "00_parameters.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "17_verify_DKO2_knockout")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

KO_GENES      <- c("Mlh1", "Msh2")            # the knocked-out genes
PARTNER_GENES <- c("Msh6", "Pms2")            # obligate heterodimer partners
TUMOR_CELLTYPE <- "Tumor"

# =============================================================================
# Load
# =============================================================================
seu <- readRDS(SEURAT_RDS)
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
sample_levels <- c(WT_SAMPLES, DKO_SAMPLES)
sample_levels <- sample_levels[sample_levels %in% unique(as.character(seu$sample))]
seu$sample    <- factor(seu$sample, levels = sample_levels)
DefaultAssay(seu) <- "RNA"

# Ensure a log-normalized data layer exists (this object ships counts only)
if (all(dim(GetAssayData(seu, assay = "RNA", layer = "data")) == 0) ||
    length(GetAssayData(seu, assay = "RNA", layer = "data")@x) == 0) {
  cat("Normalizing RNA counts (data layer was empty)...\n")
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
}

# gene availability
all_probe_genes <- c(KO_GENES, PARTNER_GENES)
present <- all_probe_genes[all_probe_genes %in% rownames(seu)]
absent  <- setdiff(all_probe_genes, rownames(seu))
cat("KO/partner genes present in panel:", paste(present, collapse = ", "), "\n")
if (length(absent) > 0)
  cat("!! NOT in panel:", paste(absent, collapse = ", "), "\n")

# =============================================================================
# (1) DIRECT — KO-gene expression in TUMOR cells, per sample
# =============================================================================
tumor_idx <- which(seu$Tier1_celltype == TUMOR_CELLTYPE)
cat(sprintf("Tumor cells: %d\n", length(tumor_idx)))

expr_data <- GetAssayData(seu, assay = "RNA", layer = "data")   # log-normalized
smp_tum   <- droplevels(seu$sample[tumor_idx])
cond_tum  <- seu$condition[tumor_idx]

direct_tbl <- do.call(rbind, lapply(present, function(g) {
  v <- as.numeric(expr_data[g, tumor_idx])
  data.frame(
    gene       = g,
    sample     = levels(smp_tum),
    mean_expr  = tapply(v, smp_tum, mean),
    pct_expr   = tapply(v, smp_tum, function(x) 100 * mean(x > 0)),
    n_tumor    = as.integer(table(smp_tum)),
    condition  = tapply(as.character(cond_tum), smp_tum,
                        function(x) x[1]),
    row.names  = NULL
  )
}))
direct_tbl$sample    <- factor(direct_tbl$sample, levels = sample_levels)
direct_tbl$condition <- factor(direct_tbl$condition,
                               levels = c(REFERENCE_CONDITION, TEST_CONDITION))
write.csv(direct_tbl,
          file.path(OUTPUT_DIR, "KO_gene_expression_by_sample_TUMOR.csv"),
          row.names = FALSE)
cat("\n--- KO-gene expression in tumor cells (per sample) ---\n")
print(direct_tbl[order(direct_tbl$gene, direct_tbl$sample),
                 c("gene", "sample", "condition", "mean_expr", "pct_expr", "n_tumor")],
      digits = 3, row.names = FALSE)

# bar plots: mean expr and % expressing, faceted by gene, colored by condition
p_mean <- ggplot(direct_tbl, aes(sample, mean_expr, fill = condition)) +
  geom_col() +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "KO-gene mean expression in TUMOR cells", y = "mean log-norm expr", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_pct <- ggplot(direct_tbl, aes(sample, pct_expr, fill = condition)) +
  geom_col() +
  facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "KO-gene % cells expressing in TUMOR cells", y = "% cells > 0", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "direct_KO_gene_expression_TUMOR.pdf"),
       p_mean / p_pct, width = 3.2 * length(present), height = 8)

# =============================================================================
# (1b) WITHIN-SAMPLE CONTROL — tumor vs non-tumor KO-gene ratio, per sample
# -----------------------------------------------------------------------------
# The KO is tumor-cell-specific (Cre driver), so stromal/immune cells in the
# SAME sample retain both alleles and act as an internal reference that cancels
# per-sample technical sensitivity. Successful KO => tumor/non-tumor << 1;
# WT and a FAILED KO => tumor/non-tumor ~ 1.
# =============================================================================
nontumor_idx <- which(seu$Tier1_celltype != TUMOR_CELLTYPE)
smp_all      <- droplevels(seu$sample)
is_tum_all   <- seu$Tier1_celltype == TUMOR_CELLTYPE

ratio_tbl <- do.call(rbind, lapply(present, function(g) {
  v <- as.numeric(expr_data[g, ])
  tum_mean <- tapply(v[is_tum_all],  smp_all[is_tum_all],  mean)
  non_mean <- tapply(v[!is_tum_all], smp_all[!is_tum_all], mean)
  data.frame(
    gene        = g,
    sample      = names(tum_mean),
    tumor_mean  = as.numeric(tum_mean),
    nontumor_mean = as.numeric(non_mean[names(tum_mean)]),
    log2_ratio  = log2((as.numeric(tum_mean) + 1e-6) /
                       (as.numeric(non_mean[names(tum_mean)]) + 1e-6)),
    row.names   = NULL
  )
}))
ratio_tbl$sample    <- factor(ratio_tbl$sample, levels = sample_levels)
ratio_tbl$condition <- factor(ifelse(ratio_tbl$sample %in% WT_SAMPLES,
                                     REFERENCE_CONDITION, TEST_CONDITION),
                              levels = c(REFERENCE_CONDITION, TEST_CONDITION))
write.csv(ratio_tbl,
          file.path(OUTPUT_DIR, "KO_gene_tumor_vs_nontumor_ratio.csv"),
          row.names = FALSE)
cat("\n--- Tumor vs non-tumor KO-gene log2 ratio (per sample; <0 = depleted in tumor) ---\n")
print(ratio_tbl[order(ratio_tbl$gene, ratio_tbl$sample),
                c("gene", "sample", "condition", "tumor_mean", "nontumor_mean", "log2_ratio")],
      digits = 3, row.names = FALSE)

p_ratio <- ggplot(ratio_tbl, aes(sample, log2_ratio, fill = condition)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(~ gene, nrow = 1) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "KO-gene tumor vs non-tumor (internal control)",
       subtitle = "log2(tumor / non-tumor); negative = depleted in tumor (KO took)",
       y = "log2(tumor / non-tumor)", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "within_sample_tumor_vs_nontumor_ratio.pdf"),
       p_ratio, width = 3.2 * length(present), height = 4.5)

# =============================================================================
# (2) GLOBAL — per-sample tumor-cell pseudobulk clustering / PCA
# =============================================================================
counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, tumor_idx]
# sum counts per sample -> CPM-like -> log1p
pb <- sapply(levels(smp_tum), function(s) Matrix::rowSums(counts[, smp_tum == s, drop = FALSE]))
pb_cpm  <- t(t(pb) / colSums(pb) * 1e4)
pb_log  <- log1p(pb_cpm)
# use variable genes for structure
gene_var <- apply(pb_log, 1, var)
top_var  <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(500, sum(gene_var > 0)))]
pb_top   <- pb_log[top_var, ]

samp_cond <- setNames(as.character(cond_tum), as.character(smp_tum))
samp_cond <- samp_cond[!duplicated(names(samp_cond))][levels(smp_tum)]

# correlation heatmap (Spearman) with clustering
cor_mat <- cor(pb_top, method = "spearman")
write.csv(cor_mat, file.path(OUTPUT_DIR, "tumor_pseudobulk_sample_correlation.csv"))

hc <- hclust(as.dist(1 - cor_mat), method = "average")
pdf(file.path(OUTPUT_DIR, "tumor_pseudobulk_dendrogram.pdf"), width = 6, height = 5)
plot(hc, main = "Tumor pseudobulk — sample clustering (1 - Spearman)",
     xlab = "", sub = "")
dev.off()

# PCA of samples
pca <- prcomp(t(pb_top), scale. = TRUE)
var_expl <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pca_df <- data.frame(sample = rownames(pca$x),
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     condition = factor(samp_cond[rownames(pca$x)],
                                        levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 4) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "Tumor pseudobulk PCA (top-500 variable genes)",
       x = sprintf("PC1 (%.1f%%)", var_expl[1]),
       y = sprintf("PC2 (%.1f%%)", var_expl[2])) +
  theme_publication()
ggsave(file.path(OUTPUT_DIR, "tumor_pseudobulk_PCA.pdf"), p_pca, width = 6, height = 5)

cat("\nDone. Outputs ->", OUTPUT_DIR, "\n")
