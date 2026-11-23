# 20_03_ligand_receptor_cell_chat.R
# Ligand-receptor (LR) interactions of one CD8+ T subtype (default: cluster 2
# from 03_01) with its spatial neighbors, per sample, with emphasis on DKO2.
# Split out of 20_01 (which keeps the colocalization tests); same cell labels
# and neighbor graph.
#
# Neighbor graph (per sample, all cells): each cell is linked to its K_NN = 30
# nearest neighbors that lie within RADIUS = 50 um; links are undirected (a link
# exists if either cell is among the other's neighbors).
#
# LR pairs: CellChatDB mouse, Secreted / Cell-Cell Contact / ECM-Receptor pairs
# whose genes are all on the panel; complexes = min of subunits. A pair is
# tested for a partner type only if the focal subtype has >= MIN_LINKS_LR links
# to it and both genes are detected on >= MIN_LR_LINK_FRAC of those links.
# For focal -> partner type t (focal expresses ligand) and partner type t ->
# focal (focal expresses receptor), score = mean over focal-t links of
# ligand(sender) x receptor(receiver) on log-normalized expression.
# Null: subtype labels permuted among QC-pass CD8+ T cells (graph, positions and
# expression fixed), so a significant pair is specific to the focal subtype's
# contacts. One-sided permutation p (score greater than null); FDR within sample.
#
# Read-only on existing RDS files. Outputs -> output/R/20_ligand_receptor

library(Seurat)
library(Matrix)
library(FNN)
library(ggplot2)
library(dplyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "20_ligand_receptor")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette and helpers
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

# =============================================================================
# Parameters
# =============================================================================
SEED              <- 42
K_NN              <- 30
RADIUS            <- 50          # spatial units (microns)
N_PERM_LR         <- 5000        # LR permutations
TCELL_ASSIGN_CSV  <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters",
                               "tcell_cluster_assignments.csv")   # from 03_01
TCELL_CLUSTER_COL <- "tcell_clust_res.0.5"
FOCAL_CLUSTER     <- "2"
FOCUS_SAMPLE      <- "DKO2"
TCELL_TYPE        <- "CD8+ T"
MIN_LINKS_LR      <- 30          # min focal-partner links to test an LR pair
MIN_LR_LINK_FRAC  <- 0.10        # min fraction of focal-partner links with both genes detected
LR_DB_FILE        <- file.path(PROJECT_DIR, "input", "LR", "CellChatDB.mouse.rda")
LR_DB_URL         <- "https://github.com/jinworks/CellChat/raw/main/data/CellChatDB.mouse.rda"
LR_ANNOTATIONS    <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor")

FOCAL <- paste0("CD8T_C", FOCAL_CLUSTER)

params <- list(SEED = SEED, K_NN = K_NN, RADIUS = RADIUS, N_PERM_LR = N_PERM_LR,
               TCELL_CLUSTER_COL = TCELL_CLUSTER_COL, FOCAL = FOCAL, FOCUS_SAMPLE = FOCUS_SAMPLE,
               MIN_LINKS_LR = MIN_LINKS_LR, MIN_LR_LINK_FRAC = MIN_LR_LINK_FRAC,
               LR_ANNOTATIONS = paste(LR_ANNOTATIONS, collapse = "; "))
writeLines(c(sprintf("20_03 ligand-receptor (CellChatDB), run %s", format(Sys.time())),
             sprintf("%-18s = %s", names(params), unlist(params))),
           file.path(OUTPUT_DIR, "run_parameters.txt"))

# =============================================================================
# 1. Load data (read-only)
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
meta <- seu@meta.data
meta$condition <- gsub("^KO$", TEST_CONDITION, meta$condition)
meta$condition <- factor(meta$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
meta$sample <- factor(meta$sample, levels = intersect(names(SAMPLE_COLORS), unique(meta$sample)))
xy <- Embeddings(seu, "spatial")

# LR database (downloaded once into input/LR if missing)
if (!file.exists(LR_DB_FILE)) {
  dir.create(dirname(LR_DB_FILE), recursive = TRUE, showWarnings = FALSE)
  download.file(LR_DB_URL, LR_DB_FILE, mode = "wb")
}
load(LR_DB_FILE)   # CellChatDB.mouse
lr_db <- CellChatDB.mouse
panel <- rownames(seu[["RNA"]])

expand_units <- function(x) {
  if (x %in% rownames(lr_db$complex)) {
    s <- unlist(lr_db$complex[x, ]); s[s != ""]
  } else x
}
ia <- lr_db$interaction[lr_db$interaction$annotation %in% LR_ANNOTATIONS, ]
ia$lig_genes <- lapply(ia$ligand, expand_units)
ia$rec_genes <- lapply(ia$receptor, expand_units)
ia <- ia[sapply(ia$lig_genes, function(g) all(g %in% panel)) &
         sapply(ia$rec_genes, function(g) all(g %in% panel)), ]
ia <- ia[!duplicated(ia$interaction_name), ]
cat(sprintf("LR pairs with all genes on the panel: %d\n", nrow(ia)))

lr_genes <- unique(c(unlist(ia$lig_genes), unlist(ia$rec_genes)))
expr_data <- GetAssayData(seu, assay = "RNA", layer = "data")[lr_genes, , drop = FALSE]
rm(seu); invisible(gc())

# Complex expression = min over subunits; returns cells x pairs
unit_expr <- function(cells, gene_lists) {
  E <- as.matrix(t(expr_data[, cells, drop = FALSE]))
  out <- sapply(gene_lists, function(g) do.call(pmin, lapply(g, function(x) E[, x])))
  matrix(out, nrow = length(cells))
}

# Cell labels: Tier1, with QC-pass CD8+ T split into 03_01 subtypes
tc <- read.csv(TCELL_ASSIGN_CSV, stringsAsFactors = FALSE, check.names = FALSE)
meta$label <- as.character(meta$Tier1_celltype)
meta$label[meta$label == TCELL_TYPE] <- "CD8T_lowQC"
meta[tc$cell, "label"] <- paste0("CD8T_C", tc[[TCELL_CLUSTER_COL]])
subtype_levels <- paste0("CD8T_C", sort(unique(tc[[TCELL_CLUSTER_COL]])))
label_levels <- c(setdiff(TIER1_ORDER, TCELL_TYPE), subtype_levels, "CD8T_lowQC")
meta$label <- factor(meta$label, levels = intersect(label_levels, unique(meta$label)))
stopifnot(FOCAL %in% levels(meta$label))

# =============================================================================
# 2. Helpers
# =============================================================================
build_graph <- function(coords) {
  n <- nrow(coords)
  nn <- get.knn(coords, k = K_NN)
  i <- rep(seq_len(n), K_NN); j <- as.vector(nn$nn.index); d <- as.vector(nn$nn.dist)
  keep <- d <= RADIUS
  a <- pmin(i[keep], j[keep]); b <- pmax(i[keep], j[keep])
  u <- !duplicated(a * (n + 1) + b)
  a <- a[u]; b <- b[u]
  sparseMatrix(i = c(a, b), j = c(b, a), x = 1, dims = c(n, n))
}

# Random 0/1 selection matrix: n_perm rows, each with `size` ones among `pool` columns
perm_matrix <- function(n_perm, pool_size, size) {
  rows <- rep(seq_len(n_perm), each = size)
  cols <- as.vector(replicate(n_perm, sample.int(pool_size, size)))
  sparseMatrix(i = rows, j = cols, x = 1, dims = c(n_perm, pool_size))
}

p_greater <- function(obs, null) (1 + sum(null >= obs)) / (1 + length(null))

# =============================================================================
# 3. Per-sample LR scores on focal links
# =============================================================================
samples <- levels(meta$sample)
lr_res <- list()

for (si in seq_along(samples)) {
  s <- samples[si]
  set.seed(SEED + si)
  idx <- which(meta$sample == s)
  m <- meta[idx, ]
  A <- build_graph(xy[idx, ])

  is_focal <- m$label == FOCAL
  cd8_pool <- which(m$label %in% subtype_levels)
  n_focal <- sum(is_focal)
  cat(sprintf("\n--- %s: %d cells, %d links, focal %d, CD8 %d ---\n",
              s, length(idx), sum(A) / 2, n_focal, length(cd8_pool)))
  if (n_focal < MIN_CELLS_PER_SAMPLE) { cat("  SKIPPED: too few focal cells\n"); next }

  cd8_cells <- idx[cd8_pool]
  L_cd8 <- unit_expr(rownames(meta)[cd8_cells], ia$lig_genes)
  R_cd8 <- unit_expr(rownames(meta)[cd8_cells], ia$rec_genes)
  L_all <- unit_expr(rownames(meta)[idx], ia$lig_genes)
  R_all <- unit_expr(rownames(meta)[idx], ia$rec_genes)
  A_cd8 <- A[cd8_pool, , drop = FALSE]
  P_lr <- perm_matrix(N_PERM_LR, length(cd8_pool), n_focal)
  focal_in_pool <- as.numeric(is_focal[cd8_pool])

  for (pt in levels(m$label)) {
    in_t <- as.numeric(m$label == pt)
    links <- as.vector(A_cd8 %*% in_t)           # links of each CD8 cell to partner type
    n_links_obs <- sum(links * focal_in_pool)
    if (n_links_obs < MIN_LINKS_LR) next
    n_links_null <- as.vector(P_lr %*% links)

    for (dir in c("focal_to_partner", "partner_to_focal")) {
      # per CD8 cell: sum over its links to type pt of ligand(sender) x receptor(receiver),
      # and number of those links on which both genes are detected
      if (dir == "focal_to_partner") {
        S  <- L_cd8 * as.matrix(A_cd8 %*% (R_all * in_t))
        Sn <- (L_cd8 > 0) * as.matrix(A_cd8 %*% ((R_all > 0) * in_t))
      } else {
        S  <- R_cd8 * as.matrix(A_cd8 %*% (L_all * in_t))
        Sn <- (R_cd8 > 0) * as.matrix(A_cd8 %*% ((L_all > 0) * in_t))
      }
      obs <- colSums(S * focal_in_pool) / n_links_obs
      frac_detected <- colSums(Sn * focal_in_pool) / n_links_obs
      keep <- which(frac_detected >= MIN_LR_LINK_FRAC)
      if (length(keep) == 0) next
      null <- as.matrix(P_lr %*% S[, keep, drop = FALSE]) / pmax(n_links_null, 1)
      lr_res[[paste(s, pt, dir)]] <- data.frame(
        sample = s, condition = as.character(m$condition[1]), partner = pt, direction = dir,
        interaction = ia$interaction_name[keep], pathway = ia$pathway_name[keep],
        ligand = ia$ligand[keep], receptor = ia$receptor[keep], annotation = ia$annotation[keep],
        n_links = n_links_obs, frac_links_detected = frac_detected[keep],
        score = obs[keep], null_mean = colMeans(null),
        log2_ratio = log2((obs[keep] + 1e-6) / (colMeans(null) + 1e-6)),
        p = mapply(function(o, k) p_greater(o, null[, k]), obs[keep], seq_along(keep))
      )
    }
  }
}

# =============================================================================
# 4. Assemble results
# =============================================================================
lr_df <- bind_rows(lr_res) %>%
  group_by(sample) %>% mutate(padj = p.adjust(p, "fdr")) %>% ungroup() %>%
  arrange(sample, p)
write.csv(lr_df, file.path(OUTPUT_DIR, "focal_LR_scores_all.csv"), row.names = FALSE)
lr_sig <- lr_df %>% filter(padj < FDR_THRESHOLD, score > 0)
write.csv(lr_sig, file.path(OUTPUT_DIR, "focal_LR_significant.csv"), row.names = FALSE)
cat(sprintf("\nSignificant LR pairs (FDR < %.2f) per sample:\n", FDR_THRESHOLD))
print(table(factor(lr_sig$sample, levels = samples)))

cxcl16_lr <- lr_df %>% filter(interaction == "CXCL16_CXCR6")
write.csv(cxcl16_lr, file.path(OUTPUT_DIR, "focal_LR_CXCL16_CXCR6.csv"), row.names = FALSE)

# =============================================================================
# 5. Plots
# =============================================================================
sig_label <- function(p) case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
                                   p < 0.05 ~ "*", p < FDR_THRESHOLD ~ "^", TRUE ~ "")

# 5a. LR dot plot for the focus sample (significant pairs)
d_lr <- lr_sig %>% filter(sample == FOCUS_SAMPLE)
if (nrow(d_lr) > 0) {
  d_lr <- d_lr %>%
    mutate(pair = paste0(interaction, ifelse(direction == "focal_to_partner", " (focal L)", " (focal R)")),
           partner = factor(partner, levels = label_levels))
  p_lr <- ggplot(d_lr, aes(x = partner, y = reorder(pair, log2_ratio))) +
    geom_point(aes(size = -log10(padj), color = log2_ratio)) +
    scale_color_gradient(low = "#FDD49E", high = "#B30000", name = "log2(score /\nnull mean)") +
    scale_size_continuous(name = "-log10(FDR)", range = c(1.5, 5)) +
    labs(title = sprintf("%s: LR pairs on %s links (FDR < %.2f)", FOCUS_SAMPLE, FOCAL, FDR_THRESHOLD),
         subtitle = "Null: subtype labels shuffled among CD8+ T cells", x = "Partner type", y = NULL) +
    theme_publication(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf_png(p_lr, file.path(OUTPUT_DIR, paste0("LR_significant_", FOCUS_SAMPLE)),
               width = 0.35 * n_distinct(d_lr$partner) + 6, height = 0.16 * n_distinct(d_lr$pair) + 2.5)
}

# 5b. CXCL16 -> CXCR6 (partner ligand, focal receptor) across samples and partners
d_cx <- cxcl16_lr %>% filter(direction == "partner_to_focal") %>%
  mutate(sample = factor(sample, levels = samples), partner = factor(partner, levels = rev(label_levels)))
if (nrow(d_cx) > 0) {
  p_cx <- ggplot(d_cx, aes(x = sample, y = partner)) +
    geom_point(aes(size = score, color = log2_ratio)) +
    geom_text(aes(label = sig_label(padj)), nudge_x = 0.3, size = 3) +
    scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", name = "log2(score /\nnull mean)") +
    labs(title = sprintf("CXCL16 (partner) -> CXCR6 (%s) on graph links", FOCAL),
         x = NULL, y = "Partner type", size = "Mean L x R") +
    theme_publication(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf_png(p_cx, file.path(OUTPUT_DIR, "LR_CXCL16_CXCR6_across_samples"), width = 7, height = 6)
}

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
cat("\nOutputs saved to:", OUTPUT_DIR, "\nDone.\n")
