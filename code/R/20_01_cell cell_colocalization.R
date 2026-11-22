# 20_01_cell cell_colocalization.R
# Spatial colocalization and ligand-receptor (LR) interactions of one CD8+ T
# subtype (default: cluster 2 from 03_01) with other cells, on a per-sample
# spatial neighbor graph, with emphasis on DKO2.
#
# Neighbor graph (per sample, all cells): each cell is linked to its K_NN = 30
# nearest neighbors that lie within RADIUS = 50 um; links are undirected (a link
# exists if either cell is among the other's neighbors).
#
# Colocalization of the focal subtype (FOCAL) with a target set = number of
# graph links between FOCAL cells and target cells. Targets:
#   (a) Tumor cells (Tier1 "Tumor")
#   (b) Cxcl16-expressing cells (raw Cxcl16 count >= CXCL16_MIN_COUNTS = 3), any
#       cell type except CD8+ T. Cxcl16 is detected at 1 count in ~25% of every
#       cell type (background), so 1 count is not specific; thresholds in
#       CXCL16_SENSITIVITY are reported as a sensitivity analysis.
# Statistics (per sample; permutations keep the graph and all positions fixed):
#   Test 1  FOCAL vs other CD8+ T subtypes. Statistic: log odds ratio
#           log[(L(F,T)+.5)/(L(F,notT)+.5)] - log[(L(O,T)+.5)/(L(O,notT)+.5)],
#           L = links, F = focal cells, O = other QC-pass CD8+ T cells.
#           Null: subtype labels permuted among CD8+ T cells. Asks whether the
#           focal subtype contacts the target more than other CD8+ T cells do,
#           controlling for where CD8+ T cells sit in the tissue.
#   Test 2  FOCAL vs random cells. Statistic: fraction of focal links that hit
#           the target. Null: the same number of cells drawn at random from all
#           non-target cells of the sample (z-score, empirical p).
#   Test 3  DKO2 vs other samples. Statistic: logOR(DKO2) - mean(logOR other
#           samples), null from the per-sample Test 1 permutations (independent
#           across samples). Also DKO vs WT on per-sample logOR (Welch t-test).
# Neighborhood profile: neighbor-type counts of focal cells for every cell type,
# z-scored against both nulls above (FDR across cell types within sample).
#
# LR analysis on focal links (CellChatDB mouse, Secreted / Cell-Cell Contact /
# ECM-Receptor pairs whose genes are all on the panel; complexes = min of
# subunits). A pair is tested for a partner type only if both genes are detected
# on >= MIN_LR_LINK_FRAC of the focal-partner links. For focal -> partner type t (focal expresses ligand) and
# partner type t -> focal (focal expresses receptor), score = mean over focal-t
# links of ligand(sender) x receptor(receiver) on log-normalized expression.
# Null: subtype labels permuted among CD8+ T cells (as Test 1), so a significant
# pair is specific to the focal subtype's contacts. FDR within sample.
#
# Read-only on existing RDS files. Outputs -> output/R/20_cell_cell_colocalization

library(Seurat)
library(Matrix)
library(FNN)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(lme4)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "20_cell_cell_colocalization")
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
N_PERM            <- 2000        # colocalization permutations
N_PERM_LR         <- 5000        # LR permutations
TCELL_ASSIGN_CSV  <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters",
                               "tcell_cluster_assignments.csv")   # from 03_01
TCELL_CLUSTER_COL <- "tcell_clust_res.0.5"
FOCAL_CLUSTER     <- "2"
COMPARE_CLUSTERS  <- c("1")     # other subtypes also tested vs random cells (exclusion)
FOCUS_SAMPLE      <- "DKO2"
TCELL_TYPE        <- "CD8+ T"
TUMOR_TYPE        <- "Tumor"
CXCL16_GENE       <- "Cxcl16"
CXCL16_MIN_COUNTS <- 3           # >= 3 counts: Cxcl16 is detected at 1 count in ~25% of every cell type
CXCL16_SENSITIVITY <- c(2, 5)   # extra thresholds reported as sensitivity targets
MIN_LINKS_LR      <- 30          # min focal-partner links to test an LR pair
MIN_LR_LINK_FRAC  <- 0.10        # min fraction of focal-partner links with both genes detected
LR_DB_FILE        <- file.path(PROJECT_DIR, "input", "LR", "CellChatDB.mouse.rda")
LR_DB_URL         <- "https://github.com/jinworks/CellChat/raw/main/data/CellChatDB.mouse.rda"
LR_ANNOTATIONS    <- c("Secreted Signaling", "Cell-Cell Contact", "ECM-Receptor")

FOCAL <- paste0("CD8T_C", FOCAL_CLUSTER)
COMPARE_SUBTYPES <- paste0("CD8T_C", COMPARE_CLUSTERS)

params <- list(SEED = SEED, K_NN = K_NN, RADIUS = RADIUS, N_PERM = N_PERM, N_PERM_LR = N_PERM_LR,
               TCELL_CLUSTER_COL = TCELL_CLUSTER_COL, FOCAL = FOCAL, FOCUS_SAMPLE = FOCUS_SAMPLE,
               CXCL16_MIN_COUNTS = CXCL16_MIN_COUNTS,
               CXCL16_SENSITIVITY = paste(CXCL16_SENSITIVITY, collapse = ", "),
               MIN_LINKS_LR = MIN_LINKS_LR, MIN_LR_LINK_FRAC = MIN_LR_LINK_FRAC,
               LR_ANNOTATIONS = paste(LR_ANNOTATIONS, collapse = "; "))
writeLines(c(sprintf("20_01 cell-cell colocalization, run %s", format(Sys.time())),
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

lr_genes <- unique(c(unlist(ia$lig_genes), unlist(ia$rec_genes), CXCL16_GENE))
expr_data <- GetAssayData(seu, assay = "RNA", layer = "data")[lr_genes, , drop = FALSE]
cxcl16_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[CXCL16_GENE, ]
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

is_cd8_any <- grepl("^CD8T_", meta$label)
meta$is_tumor  <- meta$label == TUMOR_TYPE
meta$is_cxcl16 <- cxcl16_counts[rownames(meta)] >= CXCL16_MIN_COUNTS & !is_cd8_any
TARGETS <- c(Tumor = "is_tumor", Cxcl16_expressing = "is_cxcl16")
for (th in CXCL16_SENSITIVITY) {
  col <- paste0("is_cxcl16_ge", th)
  meta[[col]] <- cxcl16_counts[rownames(meta)] >= th & !is_cd8_any
  TARGETS[paste0("Cxcl16_ge", th, "_sensitivity")] <- col
}

cat("\nCells per sample (focal / CD8 QC-pass / tumor / Cxcl16+):\n")
print(meta %>% group_by(sample) %>%
        summarise(focal = sum(label == FOCAL), cd8 = sum(label %in% subtype_levels),
                  tumor = sum(is_tumor), cxcl16 = sum(is_cxcl16), .groups = "drop"))

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

log_or <- function(fT, fN, oT, oN) log((fT + 0.5) / (fN + 0.5)) - log((oT + 0.5) / (oN + 0.5))
p_greater <- function(obs, null) (1 + sum(null >= obs)) / (1 + length(null))
p_two     <- function(obs, null) min(1, 2 * min(p_greater(obs, null), (1 + sum(null <= obs)) / (1 + length(null))))

# =============================================================================
# 3. Per-sample colocalization and LR analysis
# =============================================================================
samples <- levels(meta$sample)
target_res <- list(); profile_res <- list(); perm_store <- list(); lr_res <- list()
cell_res <- list(); compare_res <- list()
graph_summary <- list()

for (si in seq_along(samples)) {
  s <- samples[si]
  set.seed(SEED + si)
  idx <- which(meta$sample == s)
  m <- meta[idx, ]
  A <- build_graph(xy[idx, ])
  deg <- rowSums(A)
  graph_summary[[s]] <- data.frame(sample = s, n_cells = length(idx), n_links = sum(A) / 2,
                                   mean_degree = mean(deg), frac_isolated = mean(deg == 0))

  is_focal <- m$label == FOCAL
  cd8_pool <- which(m$label %in% subtype_levels)
  n_focal <- sum(is_focal)
  cat(sprintf("\n--- %s: %d cells, %d links, mean degree %.1f, focal %d, CD8 %d ---\n",
              s, length(idx), sum(A) / 2, mean(deg), n_focal, length(cd8_pool)))
  if (n_focal < MIN_CELLS_PER_SAMPLE) { cat("  SKIPPED: too few focal cells\n"); next }

  P_cd8 <- perm_matrix(N_PERM, length(cd8_pool), n_focal)

  # --- 3a. Target tests -------------------------------------------------------
  perm_store[[s]] <- list()
  for (tg in names(TARGETS)) {
    is_t <- m[[TARGETS[[tg]]]]
    nT <- as.vector(A %*% as.numeric(is_t))     # links to target, per cell
    nN <- deg - nT

    fT <- sum(nT[is_focal]); fN <- sum(nN[is_focal])
    oth <- setdiff(cd8_pool, which(is_focal))
    oT <- sum(nT[oth]); oN <- sum(nN[oth])
    lor_obs <- log_or(fT, fN, oT, oN)

    # Test 1: permute subtype labels among CD8+ T cells
    cT_tot <- sum(nT[cd8_pool]); cN_tot <- sum(nN[cd8_pool])
    pfT <- as.vector(P_cd8 %*% nT[cd8_pool]); pfN <- as.vector(P_cd8 %*% nN[cd8_pool])
    lor_null <- log_or(pfT, pfN, cT_tot - pfT, cN_tot - pfN)
    perm_store[[s]][[tg]] <- list(obs = lor_obs, null = lor_null)

    # Test 2: focal vs random non-target cells
    pool <- which(!is_t)
    P_rand <- perm_matrix(N_PERM, length(pool), n_focal)
    frac_obs <- fT / (fT + fN)
    frac_null <- as.vector(P_rand %*% nT[pool]) / pmax(as.vector(P_rand %*% deg[pool]), 1)

    target_res[[paste(s, tg)]] <- data.frame(
      sample = s, condition = as.character(m$condition[1]), target = tg,
      n_focal = n_focal, n_target = sum(is_t),
      focal_links_target = fT, focal_links_total = fT + fN,
      focal_frac_links_target = frac_obs,
      focal_frac_cells_touching = mean(nT[is_focal] > 0),
      otherCD8_frac_links_target = oT / max(oT + oN, 1),
      otherCD8_frac_cells_touching = mean(nT[oth] > 0),
      logOR_vs_otherCD8 = lor_obs,
      p_vs_otherCD8_greater = p_greater(lor_obs, lor_null),
      p_vs_otherCD8_two_sided = p_two(lor_obs, lor_null),
      random_frac_links_target = mean(frac_null),
      log2_enrich_vs_random = log2(frac_obs / mean(frac_null)),
      log2_null_lo = log2(quantile(frac_null, 0.025) / mean(frac_null)),
      log2_null_hi = log2(quantile(frac_null, 0.975) / mean(frac_null)),
      z_vs_random = (frac_obs - mean(frac_null)) / sd(frac_null),
      p_vs_random_greater = p_greater(frac_obs, frac_null),
      p_vs_random_two_sided = p_two(frac_obs, frac_null)
    )

    # Per-cell link counts of QC-pass CD8+ T cells, for the condition GLMMs
    cell_res[[paste(s, tg)]] <- data.frame(
      sample = s, condition = as.character(m$condition[1]), target = tg,
      cell = rownames(m)[cd8_pool], is_focal = is_focal[cd8_pool],
      n_target = nT[cd8_pool], n_other = nN[cd8_pool],
      random_frac = mean(frac_null)
    )
  }

  # --- 3b. Neighborhood profile over all cell types --------------------------
  onehot <- sparse.model.matrix(~ label - 1, data = m)
  colnames(onehot) <- levels(m$label)
  M <- as.matrix(A %*% onehot)                   # cells x neighbor-type counts
  obs_counts <- colSums(M[is_focal, , drop = FALSE])
  null_cd8 <- as.matrix(P_cd8 %*% M[cd8_pool, ])
  P_all <- perm_matrix(N_PERM, length(idx), n_focal)
  null_all <- as.matrix(P_all %*% M)
  profile_res[[s]] <- data.frame(
    sample = s, condition = as.character(m$condition[1]), partner = colnames(M),
    focal_links = obs_counts,
    focal_frac = obs_counts / sum(obs_counts),
    z_vs_otherCD8 = (obs_counts - colMeans(null_cd8)) / apply(null_cd8, 2, sd),
    p_vs_otherCD8 = mapply(function(o, k) p_two(o, null_cd8[, k]), obs_counts, seq_along(obs_counts)),
    z_vs_random = (obs_counts - colMeans(null_all)) / apply(null_all, 2, sd),
    p_vs_random = mapply(function(o, k) p_two(o, null_all[, k]), obs_counts, seq_along(obs_counts)),
    log2_enrich_vs_otherCD8 = log2((obs_counts + 0.5) / (colMeans(null_cd8) + 0.5)),
    log2_enrich_vs_random = log2((obs_counts + 0.5) / (colMeans(null_all) + 0.5))
  ) %>% mutate(padj_vs_otherCD8 = p.adjust(p_vs_otherCD8, "fdr"),
               padj_vs_random = p.adjust(p_vs_random, "fdr"))

  # --- 3c. LR scores on focal links -------------------------------------------
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

  # --- 3d. Comparison subtypes vs random cells (run last, so the focal results
  #         above use the same random draws with or without this block) -------
  for (cl in COMPARE_SUBTYPES) {
    is_c <- m$label == cl
    n_c <- sum(is_c)
    if (n_c < MIN_CELLS_PER_SAMPLE) next
    for (tg in names(TARGETS)) {
      is_t <- m[[TARGETS[[tg]]]]
      nT <- as.vector(A %*% as.numeric(is_t))
      pool <- which(!is_t)
      Pc <- perm_matrix(N_PERM, length(pool), n_c)
      frac_obs <- sum(nT[is_c]) / max(sum(deg[is_c]), 1)
      frac_null <- as.vector(Pc %*% nT[pool]) / pmax(as.vector(Pc %*% deg[pool]), 1)
      compare_res[[paste(s, tg, cl)]] <- data.frame(
        sample = s, condition = as.character(m$condition[1]), target = tg, subtype = cl,
        n_cells = n_c, frac_links_target = frac_obs, random_frac_links_target = mean(frac_null),
        log2_enrich_vs_random = log2(frac_obs / mean(frac_null)),
        log2_null_lo = log2(quantile(frac_null, 0.025) / mean(frac_null)),
        log2_null_hi = log2(quantile(frac_null, 0.975) / mean(frac_null)),
        p_vs_random_two_sided = p_two(frac_obs, frac_null), row.names = NULL
      )
    }
  }
}

# =============================================================================
# 4. Assemble results and cross-sample tests
# =============================================================================
graph_df <- bind_rows(graph_summary)
write.csv(graph_df, file.path(OUTPUT_DIR, "graph_summary_per_sample.csv"), row.names = FALSE)

target_df <- bind_rows(target_res) %>%
  group_by(target) %>%
  mutate(padj_vs_otherCD8 = p.adjust(p_vs_otherCD8_two_sided, "fdr"),
         padj_vs_random = p.adjust(p_vs_random_two_sided, "fdr")) %>%
  ungroup()
write.csv(target_df, file.path(OUTPUT_DIR, "focal_target_colocalization_tests.csv"), row.names = FALSE)
cat("\n== Focal-target colocalization ==\n")
print(as.data.frame(target_df %>% select(sample, target, focal_frac_links_target,
                                         otherCD8_frac_links_target, logOR_vs_otherCD8,
                                         p_vs_otherCD8_two_sided, z_vs_random, p_vs_random_two_sided)),
      digits = 3)

# Test 3: focus sample vs the other samples, and DKO vs WT
cross_df <- bind_rows(lapply(names(TARGETS), function(tg) {
  ss <- names(perm_store)[sapply(perm_store, function(x) !is.null(x[[tg]]))]
  others <- setdiff(ss, FOCUS_SAMPLE)
  obs <- perm_store[[FOCUS_SAMPLE]][[tg]]$obs - mean(sapply(others, function(o) perm_store[[o]][[tg]]$obs))
  null <- perm_store[[FOCUS_SAMPLE]][[tg]]$null -
    rowMeans(sapply(others, function(o) perm_store[[o]][[tg]]$null))
  d <- target_df %>% filter(target == tg)
  tt <- t.test(d$logOR_vs_otherCD8[d$condition == TEST_CONDITION],
               d$logOR_vs_otherCD8[d$condition == REFERENCE_CONDITION])
  data.frame(target = tg,
             focus_minus_others_logOR = obs,
             p_focus_vs_others_greater = p_greater(obs, null),
             p_focus_vs_others_two_sided = p_two(obs, null),
             mean_logOR_DKO = unname(tt$estimate[1]), mean_logOR_WT = unname(tt$estimate[2]),
             p_DKO_vs_WT_welch = tt$p.value)
}))
write.csv(cross_df, file.path(OUTPUT_DIR, "focal_target_cross_sample_tests.csv"), row.names = FALSE)
cat(sprintf("\n== %s vs other samples (logOR focal vs other CD8), and DKO vs WT ==\n", FOCUS_SAMPLE))
print(cross_df, digits = 3)

profile_df <- bind_rows(profile_res)
write.csv(profile_df, file.path(OUTPUT_DIR, "focal_neighborhood_profile.csv"), row.names = FALSE)

# =============================================================================
# 4b. DKO vs WT: does colocalization change between conditions?
# =============================================================================
# Sample-level: one value per sample (3 WT + 4 DKO are the replicates).
#   Welch t-test; exact permutation test over all choose(7, 3) = 35 condition
#   labelings (smallest attainable two-sided p = 2/35 = 0.057); leave-one-
#   sample-out to check whether a single sample (e.g. DKO2) drives the result.
# Metrics per target:
#   raw contact      logit(fraction of focal links reaching the target);
#                    confounded by how much target tissue each slide has
#   vs random        log2(observed / expected) fraction, expected from random
#                    cells of the same slide (corrects for target abundance)
#   vs other CD8     log odds ratio focal vs other CD8+ T (Test 1 statistic)
# Cell-level GLMMs (binomial on per-cell link counts, random intercepts for
# sample and cell [overdispersion]):
#   enrichment  focal cells: links ~ condition + offset(logit expected), the
#               condition effect = change in enrichment vs random cells
#   preference  CD8+ T cells: links ~ is_focal * condition with a random slope
#               for is_focal by sample; the interaction = change in the focal
#               subtype's preference over other CD8+ T
# Likelihood-ratio tests; with 7 samples these p values are approximate.
cond_test <- function(v, cond, smp) {
  ok <- is.finite(v); v <- v[ok]; cond <- cond[ok]; smp <- smp[ok]
  dko <- cond == TEST_CONDITION
  diff_fun <- function(g) mean(v[g]) - mean(v[!g])
  d_obs <- diff_fun(dko)
  combs <- combn(length(v), sum(!dko))
  d_perm <- apply(combs, 2, function(w) { g <- rep(TRUE, length(v)); g[w] <- FALSE; diff_fun(g) })
  loo <- sapply(seq_along(v), function(i) {
    vv <- v[-i]; dd <- dko[-i]
    if (sum(dd) < 2 || sum(!dd) < 2) return(c(NA, NA))
    c(mean(vv[dd]) - mean(vv[!dd]), t.test(vv[dd], vv[!dd])$p.value)
  })
  colnames(loo) <- smp
  data.frame(mean_DKO = mean(v[dko]), mean_WT = mean(v[!dko]), diff_DKO_minus_WT = d_obs,
             p_welch = t.test(v[dko], v[!dko])$p.value,
             p_exact_perm = mean(abs(d_perm) >= abs(d_obs) - 1e-12),
             loo_diff_min = min(loo[1, ], na.rm = TRUE), loo_diff_max = max(loo[1, ], na.rm = TRUE),
             loo_p_welch_max = max(loo[2, ], na.rm = TRUE),
             diff_without_focus = if (FOCUS_SAMPLE %in% smp) loo[1, FOCUS_SAMPLE] else NA,
             p_welch_without_focus = if (FOCUS_SAMPLE %in% smp) loo[2, FOCUS_SAMPLE] else NA,
             row.names = NULL)
}

target_df <- target_df %>%
  mutate(logit_focal_frac = qlogis(pmin(pmax(focal_frac_links_target, 1e-4), 1 - 1e-4)))
COND_METRICS <- c(logit_focal_frac = "raw contact (logit fraction of links)",
                  log2_enrich_vs_random = "enrichment vs random cells (log2 obs/exp)",
                  logOR_vs_otherCD8 = "preference vs other CD8+ T (log OR)")

cond_sample_df <- bind_rows(lapply(names(TARGETS), function(tg) bind_rows(lapply(names(COND_METRICS), function(mt) {
  d <- target_df %>% filter(target == tg)
  data.frame(target = tg, metric = mt, cond_test(d[[mt]], d$condition, d$sample))
}))))

fit_lrt <- function(f1, f0, data) {
  ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  m1 <- tryCatch(glmer(f1, data = data, family = binomial, control = ctrl), error = function(e) NULL)
  m0 <- tryCatch(glmer(f0, data = data, family = binomial, control = ctrl), error = function(e) NULL)
  if (is.null(m1) || is.null(m0)) return(c(est = NA, se = NA, p_lrt = NA))
  term <- grep("condition", names(fixef(m1)), value = TRUE)
  term <- if (any(grepl(":", term))) grep(":", term, value = TRUE) else term
  c(est = unname(fixef(m1)[term]), se = unname(sqrt(diag(vcov(m1)))[match(term, names(fixef(m1)))]),
    p_lrt = anova(m0, m1)$`Pr(>Chisq)`[2])
}

cell_df <- bind_rows(cell_res) %>%
  filter(n_target + n_other > 0) %>%
  mutate(condition = factor(condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION)),
         is_focal = factor(ifelse(is_focal, "focal", "other"), levels = c("other", "focal")),
         offset_rand = qlogis(random_frac))
cond_glmm_df <- bind_rows(lapply(names(TARGETS), function(tg) {
  d <- cell_df %>% filter(target == tg)
  enr <- fit_lrt(cbind(n_target, n_other) ~ condition + offset(offset_rand) + (1 | sample) + (1 | cell),
                 cbind(n_target, n_other) ~ 1 + offset(offset_rand) + (1 | sample) + (1 | cell),
                 d %>% filter(is_focal == "focal"))
  # random slope: the focal-vs-other contrast may differ between samples
  pref <- fit_lrt(cbind(n_target, n_other) ~ is_focal * condition + (1 + is_focal | sample) + (1 | cell),
                  cbind(n_target, n_other) ~ is_focal + condition + (1 + is_focal | sample) + (1 | cell), d)
  data.frame(target = tg,
             model = c("enrichment_vs_random (focal cells)", "preference_vs_otherCD8 (interaction)"),
             log_odds_DKO_vs_WT = c(enr[["est"]], pref[["est"]]),
             se = c(enr[["se"]], pref[["se"]]),
             p_lrt = c(enr[["p_lrt"]], pref[["p_lrt"]]))
}))

write.csv(cond_sample_df, file.path(OUTPUT_DIR, "condition_DKOvsWT_sample_level.csv"), row.names = FALSE)
write.csv(cond_glmm_df, file.path(OUTPUT_DIR, "condition_DKOvsWT_glmm.csv"), row.names = FALSE)
cat("\n== DKO vs WT, sample level (Welch, exact permutation, leave-one-out) ==\n")
print(cond_sample_df %>% select(target, metric, mean_DKO, mean_WT, diff_DKO_minus_WT, p_welch,
                                p_exact_perm, diff_without_focus, p_welch_without_focus), digits = 3)
cat("\n== DKO vs WT, cell-level binomial GLMMs ==\n")
print(cond_glmm_df, digits = 3)

# Neighborhood profile: DKO vs WT per partner type
cond_profile_df <- bind_rows(lapply(c("log2_enrich_vs_random", "log2_enrich_vs_otherCD8"), function(mt) {
  bind_rows(lapply(unique(profile_df$partner), function(pt) {
    d <- profile_df %>% filter(partner == pt)
    if (sum(is.finite(d[[mt]])) < nrow(d) || nrow(d) < length(samples)) return(NULL)
    data.frame(metric = mt, partner = pt, cond_test(d[[mt]], d$condition, d$sample))
  })) %>% mutate(padj_welch = p.adjust(p_welch, "fdr"))
}))
write.csv(cond_profile_df, file.path(OUTPUT_DIR, "condition_DKOvsWT_neighborhood_profile.csv"),
          row.names = FALSE)

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

# 5a. Target tests: logOR focal vs other CD8 with permutation null range, per sample
null_rng <- bind_rows(lapply(names(perm_store), function(s) bind_rows(lapply(names(TARGETS), function(tg) {
  x <- perm_store[[s]][[tg]]
  data.frame(sample = s, target = tg, lo = quantile(x$null, 0.025), hi = quantile(x$null, 0.975))
}))))
d5a <- target_df %>% left_join(null_rng, by = c("sample", "target")) %>%
  mutate(sample = factor(sample, levels = samples), label = sig_label(padj_vs_otherCD8))
p_lor <- ggplot(d5a, aes(x = sample, y = logOR_vs_otherCD8, color = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = lo, ymax = hi), color = "grey70", width = 0.3) +
  geom_point(size = 3) +
  geom_text(aes(label = label), vjust = -1, color = "black", size = 4) +
  facet_wrap(~ target, scales = "free_y") +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = sprintf("%s contacts vs other CD8+ T cells (kNN k=%d, r=%g um)", FOCAL, K_NN, RADIUS),
       subtitle = "Grey bar = 95% permutation null (subtype labels shuffled among CD8+ T); stars = FDR",
       x = NULL, y = "log odds ratio (focal vs other CD8+ T)", color = NULL) +
  theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf_png(p_lor, file.path(OUTPUT_DIR, "target_logOR_vs_otherCD8"), width = 10, height = 4.5)

# 5b. Fraction of links hitting the target: focal vs other CD8 vs random
d5b <- target_df %>%
  select(sample, condition, target, focal = focal_frac_links_target,
         `other CD8+ T` = otherCD8_frac_links_target, random = random_frac_links_target) %>%
  pivot_longer(c(focal, `other CD8+ T`, random), names_to = "group", values_to = "frac") %>%
  mutate(sample = factor(sample, levels = samples),
         group = factor(group, levels = c("focal", "other CD8+ T", "random"),
                        labels = c(FOCAL, "other CD8+ T", "random cells")))
p_frac <- ggplot(d5b, aes(x = sample, y = frac, fill = group)) +
  geom_col(position = position_dodge(0.8), width = 0.75) +
  facet_wrap(~ target, scales = "free_y") +
  scale_fill_manual(values = c(setNames("#CC79A7", FOCAL), "other CD8+ T" = "#7F7F7F",
                               "random cells" = "#D9D9D9")) +
  labs(title = "Fraction of graph links that reach the target", x = NULL,
       y = "Fraction of links", fill = NULL) +
  theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf_png(p_frac, file.path(OUTPUT_DIR, "target_fraction_links"), width = 10, height = 4.5)

# 5c. Neighborhood profile heatmaps (z-scores)
for (nul in c("otherCD8", "random")) {
  d <- profile_df %>%
    mutate(z = .data[[paste0("z_vs_", nul)]], padj = .data[[paste0("padj_vs_", nul)]],
           sample = factor(sample, levels = samples),
           partner = factor(partner, levels = rev(label_levels))) %>%
    filter(!is.na(partner), is.finite(z))
  zcap <- max(3, min(10, max(abs(d$z))))
  p_prof <- ggplot(d, aes(x = sample, y = partner, fill = pmax(pmin(z, zcap), -zcap))) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig_label(padj)), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-zcap, zcap),
                         name = "z-score") +
    labs(title = sprintf("Neighbors of %s: enrichment vs %s", FOCAL,
                         if (nul == "otherCD8") "other CD8+ T cells" else "random cells"),
         subtitle = "Stars = FDR across partner types within sample", x = NULL, y = "Neighbor type") +
    theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_pdf_png(p_prof, file.path(OUTPUT_DIR, paste0("neighborhood_profile_vs_", nul)),
               width = 7, height = 7)
}

# 5d. LR dot plot for the focus sample (significant pairs)
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

# 5e. CXCL16 -> CXCR6 (partner ligand, focal receptor) across samples and partners
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

# 5f. Spatial map of the focus sample
fm <- meta[meta$sample == FOCUS_SAMPLE, ]
fxy <- xy[rownames(fm), ]
fm$group <- case_when(fm$label == FOCAL ~ FOCAL,
                      fm$label %in% subtype_levels ~ "other CD8+ T",
                      fm$is_tumor ~ "Tumor",
                      fm$is_cxcl16 ~ "Cxcl16+ cell",
                      TRUE ~ "other")
fm$group <- factor(fm$group, levels = c("other", "Tumor", "Cxcl16+ cell", "other CD8+ T", FOCAL))
df_map <- data.frame(x = fxy[, 1], y = fxy[, 2], group = fm$group)[order(fm$group), ]
p_map <- ggplot(df_map, aes(x, y, color = group, size = group)) +
  geom_point(stroke = 0) +
  scale_color_manual(values = c("other" = "grey90", "Tumor" = "#F4A6A0", "Cxcl16+ cell" = "#7FBF7B",
                                "other CD8+ T" = "#6A6A6A", setNames("#7B2CBF", FOCAL))) +
  scale_size_manual(values = c("other" = 0.15, "Tumor" = 0.2, "Cxcl16+ cell" = 0.3,
                               "other CD8+ T" = 0.5, setNames(0.9, FOCAL)), guide = "none") +
  coord_fixed() + theme_void() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(title = sprintf("%s: %s, tumor and Cxcl16+ cells", FOCUS_SAMPLE, FOCAL), color = NULL) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(OUTPUT_DIR, paste0("spatial_map_", FOCUS_SAMPLE, ".png")), p_map,
       width = 8, height = 8, dpi = 300)

# 5g. DKO vs WT: per-sample colocalization metrics with test results
main_targets <- c("Tumor", "Cxcl16_expressing")
glmm_lab <- cond_glmm_df %>%
  mutate(metric = ifelse(grepl("^enrichment", model), "log2_enrich_vs_random", "logOR_vs_otherCD8")) %>%
  select(target, metric, p_glmm = p_lrt)
d_cond <- target_df %>% filter(target %in% main_targets) %>%
  select(sample, condition, target, all_of(names(COND_METRICS))) %>%
  pivot_longer(all_of(names(COND_METRICS)), names_to = "metric", values_to = "value") %>%
  mutate(condition = factor(condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION)))
lab_cond <- cond_sample_df %>% filter(target %in% main_targets) %>%
  left_join(glmm_lab, by = c("target", "metric")) %>%
  mutate(label = sprintf("Welch p=%.2g\nexact p=%.2g%s\nw/o %s p=%.2g", p_welch, p_exact_perm,
                         ifelse(is.na(p_glmm), "", sprintf("\nGLMM p=%.2g", p_glmm)),
                         FOCUS_SAMPLE, p_welch_without_focus))
facet_labs <- labeller(metric = as_labeller(COND_METRICS, label_wrap_gen(28)))
p_cond <- ggplot(d_cond, aes(x = condition, y = value)) +
  geom_boxplot(aes(fill = condition), alpha = 0.4, outlier.shape = NA, width = 0.55) +
  geom_point(aes(color = condition), size = 2.2, position = position_jitter(width = 0.08, seed = SEED)) +
  geom_text(data = d_cond %>% filter(sample == FOCUS_SAMPLE), aes(label = sample),
            nudge_x = 0.28, size = 2.8) +
  geom_text(data = lab_cond, aes(x = 1.5, y = Inf, label = label), vjust = 1.1, size = 2.6,
            lineheight = 0.9) +
  facet_grid(metric ~ target, scales = "free_y", labeller = facet_labs) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.55))) +
  scale_fill_manual(values = CONDITION_COLORS) + scale_color_manual(values = CONDITION_COLORS) +
  labs(title = sprintf("%s colocalization: DKO vs WT (one point per sample)", FOCAL), x = NULL, y = NULL) +
  theme_publication(base_size = 10) + theme(legend.position = "none")
save_pdf_png(p_cond, file.path(OUTPUT_DIR, "condition_DKOvsWT_colocalization"), width = 7.5, height = 9)

# 5h. DKO vs WT: change in the neighborhood profile per partner type
d_cp <- cond_profile_df %>%
  mutate(partner = factor(partner, levels = rev(label_levels)),
         metric = factor(metric, levels = c("log2_enrich_vs_random", "log2_enrich_vs_otherCD8"),
                         labels = c("vs random cells", "vs other CD8+ T")),
         label = sig_label(padj_welch)) %>%
  filter(!is.na(partner))
p_cp <- ggplot(d_cp, aes(x = diff_DKO_minus_WT, y = partner)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = loo_diff_min, xmax = loo_diff_max), height = 0, color = "grey60") +
  geom_point(aes(size = -log10(p_welch), color = diff_DKO_minus_WT)) +
  geom_text(aes(label = label), nudge_y = 0.35, size = 3) +
  facet_wrap(~ metric, scales = "free_x") +
  scale_color_gradient2(low = "#2166AC", mid = "grey70", high = "#B2182B", midpoint = 0,
                        name = "DKO - WT") +
  scale_size_continuous(name = "-log10(Welch p)", range = c(1.5, 5)) +
  labs(title = sprintf("Change in %s neighborhood enrichment, DKO - WT (log2 obs/exp)", FOCAL),
       subtitle = "Grey bar = range over leave-one-sample-out; stars = FDR across partner types",
       x = "DKO - WT", y = "Neighbor type") +
  theme_publication(base_size = 10)
save_pdf_png(p_cp, file.path(OUTPUT_DIR, "condition_DKOvsWT_neighborhood_profile"), width = 10, height = 6.5)

# 5i. Per sample: focal and comparison subtypes vs random cells
#     (log2 enrichment of target contacts; below 0 = excluded / avoids target)
vs_random <- bind_rows(
  target_df %>% transmute(sample, condition, target, subtype = FOCAL, n_cells = n_focal,
                          frac_links_target = focal_frac_links_target, random_frac_links_target,
                          log2_enrich_vs_random, log2_null_lo, log2_null_hi, p_vs_random_two_sided),
  bind_rows(compare_res)) %>%
  filter(target %in% main_targets) %>%
  group_by(target, subtype) %>%
  mutate(padj_vs_random = p.adjust(p_vs_random_two_sided, "fdr")) %>% ungroup() %>%
  mutate(sample = factor(sample, levels = samples),
         subtype = factor(subtype, levels = c(FOCAL, COMPARE_SUBTYPES)),
         pct_links_to_target = 100 * frac_links_target,
         random_pct_links_to_target = 100 * random_frac_links_target,
         fold_vs_random = 2^log2_enrich_vs_random,
         call = case_when(padj_vs_random >= FDR_THRESHOLD ~ "no difference",
                          log2_enrich_vs_random > 0 ~ "closer than random",
                          TRUE ~ "farther than random (excluded)"))
write.csv(vs_random, file.path(OUTPUT_DIR, "subtypes_vs_random_per_sample.csv"), row.names = FALSE)
cat("\n== CD8+ T subtypes vs random cells, per sample ==\n")
print(as.data.frame(vs_random %>% select(subtype, sample, target, n_cells, pct_links_to_target,
                                         random_pct_links_to_target, fold_vs_random,
                                         padj_vs_random, call)), digits = 3)

vs_random_cond <- bind_rows(lapply(split(vs_random, list(vs_random$target, vs_random$subtype), drop = TRUE),
  function(d) data.frame(target = d$target[1], subtype = d$subtype[1],
                         cond_test(d$log2_enrich_vs_random, d$condition, as.character(d$sample)))))
write.csv(vs_random_cond, file.path(OUTPUT_DIR, "subtypes_vs_random_DKOvsWT.csv"), row.names = FALSE)
cat("\n== DKO vs WT on log2 enrichment vs random cells ==\n")
print(vs_random_cond %>% select(target, subtype, mean_WT, mean_DKO, diff_DKO_minus_WT,
                                p_welch, p_exact_perm), digits = 3)

p_vr <- ggplot(vs_random, aes(x = sample, y = log2_enrich_vs_random)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = log2_null_lo, ymax = log2_null_hi), color = "grey70", width = 0.3) +
  geom_point(aes(color = condition), size = 3.2) +
  geom_text(aes(label = sig_label(padj_vs_random),
                vjust = ifelse(log2_enrich_vs_random > 0, -0.8, 1.8)), size = 4.5) +
  geom_text(aes(y = -Inf, label = sprintf("%.1f%%\nvs %.1f%%", pct_links_to_target,
                                          random_pct_links_to_target)),
            vjust = -0.2, size = 2.5, lineheight = 0.85, color = "grey30") +
  facet_grid(target ~ subtype, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.35, 0.15))) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "Are CD8+ T subtypes next to the target more than random cells? (per sample)",
       subtitle = paste0("y = log2(% of focal links reaching target / % expected for random cells); ",
                         "above 0 = closer, below 0 = farther\n",
                         "grey bar = 95% range for random cells; stars = FDR across samples; ",
                         "bottom = % subtype links vs % random"),
       x = NULL, y = "log2 enrichment vs random cells", color = NULL) +
  theme_publication() + theme(plot.subtitle = element_text(size = 8))
save_pdf_png(p_vr, file.path(OUTPUT_DIR, "subtypes_vs_random_per_sample"), width = 13, height = 7)

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))
cat("\nOutputs saved to:", OUTPUT_DIR, "\nDone.\n")
