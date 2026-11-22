# 02_09_hclust optimization.R
# Alternative domain merge, compared with the current 02_03 merge (k = MERGE_K).
#
# 02_03 merges the Banksy clusters with Ward.D2 on the unweighted distances
# between cluster centroids, so a 200-cell cluster pulls on the tree as hard as a
# 20,000-cell one and small, sample-specific clusters can flip the shape of the
# tree. MERGE_METHOD picks the alternative, cut at K_NEW domains:
#   unweighted  02_03 as is
#   weighted    each centroid carries its cluster size: the starting
#               dissimilarity between clusters i and j is the Ward merge cost
#                   d_ij = sqrt(2 * n_i * n_j / (n_i + n_j)) * ||c_i - c_j||
#               and hclust(members = n) applies the size-aware Lance-Williams
#               update, which makes the tree the upper part of a cell-level
#               Ward.D2 tree whose bottom merges are the Banksy clusters
#   minsize     clusters smaller than MIN_CLUSTER_FRAC of all cells are dissolved
#               and their cells take the majority label of their SNN neighbours
#               (the graph the clusters came from); the unweighted 02_03 merge is
#               then run on the remaining clusters
#
# New domains are named after the current domain they overlap most; when a
# current domain is split, the parts are suffixed by tumour-cell fraction
# (provisional names, local to this script; 00_parameters.R is not changed).
# Sections 2-9 mirror 02_03 on the new merge; section 10 compares it with the
# current merge; section 11 (optional) re-runs Harmony with several seeds and
# compares the reproducibility of the two merges.
#
# Read-only: nothing is saved except CSV tables and plots, all in a folder used
# only by this script (no RDS output; DKO_banksy_merged.rds is never overwritten).

library(Seurat)
library(harmony)
library(Matrix)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dplyr)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

# Source shared color palette
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

# --- Parameters ----------------------------------------------------------------
MERGE_METHOD     <- "minsize"     # "unweighted" (= 02_03), "weighted" or "minsize"
K_NEW            <- 5             # domains in the new merge (current merge: MERGE_K)
MIN_CLUSTER_FRAC <- 0.005         # minsize: clusters below this fraction of cells are dissolved
RUN_SEED_CHECK   <- FALSE         # section 11: re-run Harmony + clustering (slow)
HARMONY_SEEDS    <- c(42, 1, 2, 3, 4)   # 42 = the seed used in 02_01

stopifnot(MERGE_METHOD %in% c("unweighted", "weighted", "minsize"))
METHOD_LABELS <- c(unweighted = "Unweighted Ward.D2",
                   weighted   = "Size-weighted Ward.D2",
                   minsize    = sprintf("Min-size %.1f%% + SNN reassignment, unweighted Ward.D2",
                                        100 * MIN_CLUSTER_FRAC))
CURRENT_LABEL <- sprintf("Current: unweighted Ward.D2, k = %d (02_03)", MERGE_K)
NEW_LABEL     <- sprintf("New: %s, k = %d", METHOD_LABELS[MERGE_METHOD], K_NEW)

DOMAIN_COL <- DOMAIN_SOURCE_COL
NEW_TAG    <- sprintf("%s_k%d", MERGE_METHOD, K_NEW)

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "02_banksy_lambda08", "hclust_optimization", NEW_TAG)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Helpers -------------------------------------------------------------------
ari <- function(a, b) mclust::adjustedRandIndex(a, b)

# Ward.D2 tree on cluster centroids; weighted = TRUE gives every centroid its
# cluster size, weighted = FALSE reproduces 02_03.
centroid_tree <- function(emb, clusters, weighted) {
  clusters <- as.character(clusters)
  ids <- sort(unique(clusters))
  n   <- as.numeric(table(clusters)[ids])
  cen <- rowsum(emb, clusters)[ids, , drop = FALSE] / n
  d   <- dist(cen, method = MERGE_DIST_METHOD)
  if (!weighted) return(hclust(d, method = MERGE_LINK_METHOD))
  d <- as.matrix(d) * sqrt(2 * outer(n, n) / outer(n, n, "+"))
  hclust(as.dist(d), method = MERGE_LINK_METHOD, members = setNames(n, ids))
}

# Dissolve clusters with fewer than min_size cells: each of their cells takes the
# label most strongly connected to it in the SNN graph among the kept clusters.
# Repeated so cells whose neighbours were all in small clusters are reached
# through already-reassigned cells; any cell left unconnected goes to the
# nearest kept centroid in `emb`.
reassign_small <- function(clusters, graph, emb, min_size) {
  cl    <- as.character(clusters)
  sizes <- table(cl)
  keep  <- sort(names(sizes)[sizes >= min_size])
  todo  <- !(cl %in% keep)
  while (any(todo)) {
    done  <- which(!todo)
    lab   <- sparseMatrix(i = done, j = match(cl[done], keep), x = 1,
                          dims = c(length(cl), length(keep)))
    votes <- graph[todo, , drop = FALSE] %*% lab
    has   <- rowSums(votes) > 0
    if (!any(has)) break
    idx <- which(todo)[has]
    cl[idx]   <- keep[max.col(as.matrix(votes[has, , drop = FALSE]), ties.method = "first")]
    todo[idx] <- FALSE
  }
  if (any(todo)) {
    cen <- rowsum(emb[!todo, , drop = FALSE], cl[!todo])[keep, , drop = FALSE] /
           as.numeric(table(cl[!todo])[keep])
    d2  <- outer(rowSums(emb[todo, , drop = FALSE]^2), rowSums(cen^2), "+") -
           2 * emb[todo, , drop = FALSE] %*% t(cen)
    cl[todo] <- keep[max.col(-d2, ties.method = "first")]
  }
  cl
}

# Build one merge: the cluster labels the tree is built on, the tree, and the
# cell-level domains at k.
merge_variant <- function(method, emb, clusters, graph, k) {
  labels <- as.character(clusters)
  if (method == "minsize") {
    labels <- reassign_small(labels, graph, emb, MIN_CLUSTER_FRAC * length(labels))
  }
  tree <- centroid_tree(emb, labels, weighted = method == "weighted")
  list(labels = labels, tree = tree, domains = unname(cutree(tree, k = k)[labels]))
}

# --- 1. Load object (read-only) ------------------------------------------------
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

# Recode condition and set factor levels
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$Tier1_celltype <- factor(seu$Tier1_celltype,
                             levels = intersect(TIER1_ORDER, unique(seu$Tier1_celltype)))
seu$condition <- factor(seu$condition,
                        levels = intersect(names(CONDITION_COLORS), unique(seu$condition)))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

harmony_emb    <- Embeddings(seu, "harmony")
domain_labels  <- as.character(seu@meta.data[[DOMAIN_COL]])
domain_current <- as.character(seu@meta.data[[COL_DOMAIN]])
current_names  <- unname(DOMAIN_NAMES[domain_current])
domains        <- sort(unique(domain_labels))
cluster_n      <- table(domain_labels)[domains]
snn_name       <- grep("snn", Graphs(seu), value = TRUE)[1]
snn_graph      <- as(seu[[snn_name]], "dgCMatrix")
cat(sprintf("Original clusters: %d (sizes %d-%d cells); SNN graph: %s\n", length(domains),
            min(cluster_n), max(cluster_n), snn_name))

is_immune <- seu$Tier1_celltype %in% IMMUNE_CELLTYPES
is_tumour <- seu$Tier1_celltype == "Tumor"
is_gland  <- seu$Tier1_celltype == "Basal/Luminal"

# --- 2. Build the current and new merges ---------------------------------------
current <- merge_variant("unweighted", harmony_emb, domain_labels, snn_graph, k = MERGE_K)
new     <- merge_variant(MERGE_METHOD, harmony_emb, domain_labels, snn_graph, k = K_NEW)

# Reference check: the unweighted k = MERGE_K merge must reproduce the 02_03 domains
stopifnot(isTRUE(all.equal(ari(current$domains, domain_current), 1)))

if (MERGE_METHOD == "minsize") {
  small <- cluster_n < MIN_CLUSTER_FRAC * ncol(seu)
  cat(sprintf("minsize: %d clusters (%d cells) below %.1f%% of cells reassigned via SNN; %d clusters kept\n",
              sum(small), sum(cluster_n[small]), 100 * MIN_CLUSTER_FRAC, sum(!small)))
}

# Name each new domain after the current domain it overlaps most; parts of a
# split domain are told apart by tumour-cell fraction.
tab_parent  <- table(new$domains, current_names)
parent      <- setNames(colnames(tab_parent)[max.col(unclass(tab_parent), ties.method = "first")],
                        rownames(tab_parent))
tumour_frac <- tapply(is_tumour, new$domains, mean)
new_names   <- parent
for (p in unique(parent[duplicated(parent)])) {
  ids <- names(parent)[parent == p]
  ids <- ids[order(-tumour_frac[ids])]
  sfx <- if (length(ids) == 2) c(" (tumour-rich)", " (tumour-poor)") else paste0(" ", seq_along(ids))
  new_names[ids] <- paste0(p, sfx)
}
domain_new <- unname(new_names[as.character(new$domains)])

# Display order and colours: parts of a split domain follow their parent, the
# first keeps the parent colour and the others get lighter shades of it.
NEW_ORDER  <- unlist(lapply(DOMAIN_ORDER, function(p) {
  ids <- names(parent)[parent == p]
  unname(new_names[ids[order(-tumour_frac[ids])]])
}))
NEW_COLORS <- unlist(lapply(DOMAIN_ORDER, function(p) {
  nm <- NEW_ORDER[NEW_ORDER %in% new_names[parent == p]]
  if (length(nm) == 0) return(NULL)
  setNames(colorRampPalette(c(DOMAIN_COLORS[[p]], "white"))(length(nm) + 1)[seq_along(nm)], nm)
}))
ALL_COLORS <- c(DOMAIN_COLORS, NEW_COLORS[setdiff(names(NEW_COLORS), names(DOMAIN_COLORS))])

cat(sprintf("\n%s\n", NEW_LABEL))
print(data.frame(cut_id = names(new_names), domain = unname(new_names), parent = unname(parent),
                 n_cells = as.integer(table(new$domains)[names(new_names)]),
                 tumour_frac = round(as.numeric(tumour_frac[names(new_names)]), 3)), row.names = FALSE)

# --- 3. Plot dendrograms -------------------------------------------------------
plot_dendro <- function(tree, title, k) {
  plot(tree, main = title, xlab = "Cluster", ylab = "Ward distance", sub = "", cex = 0.6)
  cut_height <- mean(sort(tree$height, decreasing = TRUE)[c(k - 1, k)])
  abline(h = cut_height, col = "red", lty = 2)
  text(x = length(tree$labels) * 0.8, y = cut_height * 1.05,
       labels = paste("k =", k), col = "red", cex = 0.9)
}
pdf(file.path(OUTPUT_DIR, "domain_dendrograms_current_vs_new.pdf"), width = 14, height = 10)
par(mfrow = c(2, 1))
plot_dendro(current$tree, CURRENT_LABEL, MERGE_K)
plot_dendro(new$tree, NEW_LABEL, K_NEW)
dev.off()
cat("Dendrograms saved.\n")

# --- 4. Cell-type composition heatmap per cluster (new dendrogram order) -------
# Clusters are those the new tree is built on (for minsize: kept clusters after
# reassignment).
hc       <- new$tree
ct_table <- table(new$labels, seu$Tier1_celltype)
ct_prop  <- prop.table(ct_table, margin = 1)
ct_prop_ordered <- ct_prop[hc$labels[hc$order], ]
tree_n   <- table(new$labels)
tree_dom <- tapply(domain_new, new$labels, `[`, 1)[rownames(ct_prop_ordered)]

pdf(file.path(OUTPUT_DIR, "heatmap_domain_celltype_composition.pdf"), width = 12, height = 14)
col_fun <- colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#6BAED6", "#08519C"))
ht <- Heatmap(ct_prop_ordered,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              clustering_method_columns = "ward.D2",
              col = col_fun,
              name = "Proportion",
              row_names_gp = gpar(fontsize = 6),
              column_names_rot = 45,
              column_names_side = "bottom",
              row_names_side = "left",
              border = TRUE,
              left_annotation = rowAnnotation(domain = tree_dom, col = list(domain = NEW_COLORS)),
              right_annotation = rowAnnotation(
                n_cells = anno_barplot(as.numeric(tree_n[rownames(ct_prop_ordered)]))),
              column_title = paste("Cell-type composition per Banksy cluster,", NEW_LABEL))
draw(ht)
dev.off()
cat("Cell-type composition heatmap saved.\n")

# --- 5. Merge map --------------------------------------------------------------
# Domain of each original cluster = majority domain of its cells (for minsize,
# dissolved clusters can be split across domains; frac_majority shows how much).
cl_dom <- table(domain_labels, domain_new)[domains, , drop = FALSE]
merge_map <- data.frame(
  original_domain = domains,
  n_cells         = as.integer(cluster_n),
  dissolved       = MERGE_METHOD == "minsize" & as.integer(cluster_n) < MIN_CLUSTER_FRAC * ncol(seu),
  merged_domain   = colnames(cl_dom)[max.col(cl_dom, ties.method = "first")],
  frac_majority   = apply(cl_dom, 1, max) / rowSums(cl_dom))
merge_map <- merge_map[order(match(merge_map$merged_domain, NEW_ORDER),
                             as.numeric(merge_map$original_domain)), ]
write.csv(merge_map, file.path(OUTPUT_DIR, "domain_merge_mapping.csv"), row.names = FALSE)

cat(sprintf("\nMerged %d clusters into %d domains (* = dissolved, assigned by majority of its cells):\n",
            length(domains), K_NEW))
for (m in NEW_ORDER) {
  sub <- merge_map[merge_map$merged_domain == m, ]
  cat(sprintf("  %s <- original: %s\n", m,
              paste0(sub$original_domain, ifelse(sub$dissolved, "*", ""), collapse = ", ")))
}

# --- 6. Condition composition per merged domain --------------------------------
ct_cond <- as.data.frame(prop.table(table(domain_new, seu$condition), margin = 1))
colnames(ct_cond) <- c("Domain", "Condition", "Proportion")
ct_cond$Domain <- factor(ct_cond$Domain, levels = NEW_ORDER)

p_cond <- ggplot(ct_cond, aes(x = Domain, y = Proportion, fill = Condition)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = "Merged Domain", y = "Proportion", title = "Condition composition per domain") +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "barplot_domain_condition.pdf"), p_cond, width = 7, height = 5)

# --- 7. Cell-type composition per merged domain --------------------------------
ct_merged <- table(domain_new, seu$Tier1_celltype)
ct_merged_prop <- prop.table(ct_merged, margin = 1)[NEW_ORDER, ]

pdf(file.path(OUTPUT_DIR, "heatmap_merged_domain_celltype.pdf"), width = 10, height = 4.5)
col_fun <- colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#A8D8E8", "#4A9FBC"))
ht <- Heatmap(ct_merged_prop,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              clustering_method_columns = "ward.D2",
              col = col_fun,
              name = "Proportion",
              column_names_rot = 45,
              column_names_side = "bottom",
              row_names_side = "left",
              border = TRUE,
              show_column_dend = FALSE,
              cell_fun = function(j, i, x, y, w, h, fill) {
                grid.text(sprintf("%.2f", ct_merged_prop[i, j]), x, y, gp = gpar(fontsize = 9, col = "black"))
              },
              column_title = "Cell-type composition per domain (new merge)")
draw(ht)
dev.off()

write.csv(as.data.frame.matrix(ct_merged_prop),
          file.path(OUTPUT_DIR, "merged_domain_celltype_proportions.csv"))

# --- 8. Cell-type distribution across merged domains (inverse perspective) -----
# Shows: what % of each cell type is in each domain; one heatmap for all cells,
# then one per condition and per sample in the same column order.
inv_heatmap <- function(idx, title, file_stem, col_order = NULL) {
  tab  <- table(seu$Tier1_celltype[idx], factor(domain_new[idx], levels = NEW_ORDER))
  full <- prop.table(tab, margin = 1)
  full[is.nan(full)] <- 0
  full <- unclass(full)
  if (!is.null(col_order)) full <- full[col_order, , drop = FALSE]
  mat <- t(full)

  pdf(file.path(OUTPUT_DIR, paste0(file_stem, ".pdf")), width = 10, height = 4.5)
  ht <- Heatmap(mat,
                cluster_rows = FALSE,
                cluster_columns = is.null(col_order),
                clustering_method_columns = "ward.D2",
                col = colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#A8D8E8", "#4A9FBC")),
                name = "Proportion",
                column_names_rot = 45,
                column_names_side = "bottom",
                row_names_side = "left",
                show_column_dend = FALSE,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  grid.text(sprintf("%.2f", mat[i, j]), x, y, gp = gpar(fontsize = 9, col = "black"))
                },
                border = TRUE,
                column_title = title)
  ht <- draw(ht)
  dev.off()
  write.csv(as.data.frame.matrix(full), file.path(OUTPUT_DIR, paste0(file_stem, ".csv")))
  colnames(mat)[column_order(ht)]
}

ref_col_order <- inv_heatmap(seq_len(ncol(seu)), "Distribution of cell types across domains",
                             "celltype_distribution_across_domains")
for (cond in levels(seu$condition)) {
  inv_heatmap(which(seu$condition == cond),
              paste("Distribution of cell types across domains in", cond),
              paste0("celltype_distribution_across_domains_", cond), ref_col_order)
}
for (samp in levels(seu$sample)) {
  inv_heatmap(which(seu$sample == samp),
              paste("Distribution of cell types across domains in", samp),
              paste0("celltype_distribution_across_domains_", samp), ref_col_order)
}
cat("Cell-type distribution heatmaps saved.\n")

# --- 9. Spatial plots: current vs new ------------------------------------------
all_spatial <- Embeddings(seu, "spatial")
max_x_span <- max(tapply(all_spatial[, 1], seu$sample, function(x) diff(range(x))))
max_y_span <- max(tapply(all_spatial[, 2], seu$sample, function(x) diff(range(x))))
half_span  <- max(max_x_span, max_y_span) / 2 * 1.05

# Plotted straight from the coordinates with scattermore (rasterized points)
# rather than subset() + DimPlot per sample, which is slow on this object.
set.seed(1)
shuffle <- sample(ncol(seu))
spatial_df <- bind_rows(
  data.frame(x = all_spatial[, 1], y = all_spatial[, 2], sample = seu$sample,
             domain = current_names, merge = CURRENT_LABEL)[shuffle, ],
  data.frame(x = all_spatial[, 1], y = all_spatial[, 2], sample = seu$sample,
             domain = domain_new, merge = NEW_LABEL)[shuffle, ]) %>%
  mutate(merge  = factor(merge, levels = c(CURRENT_LABEL, NEW_LABEL)),
         domain = factor(domain, levels = union(DOMAIN_ORDER, NEW_ORDER)))
p_spatial <- ggplot(spatial_df, aes(x, y, color = domain)) +
  scattermore::geom_scattermore(pointsize = 1.5, pixels = c(600, 600)) +
  scale_color_manual(values = ALL_COLORS, name = "Domain") +
  facet_grid(merge ~ sample, labeller = labeller(merge = label_wrap_gen(22))) +
  theme_void() +
  theme(strip.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 3)))
p_spatial <- if (CENTER_COORDS) {
  p_spatial + coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span))
} else {
  p_spatial + coord_fixed()
}
ggsave(file.path(OUTPUT_DIR, "spatial_domains_current_vs_new.png"),
       p_spatial, width = 3 * nlevels(seu$sample) + 1, height = 7, dpi = 200, bg = "white")
rm(spatial_df)
cat("Spatial plots saved.\n")

# --- 10. Comparison with the current merge -------------------------------------
# Collapsing split domains back onto their parent shows how much of the change
# is the extra split and how much is cells moving between current domains.
new_collapsed <- unname(parent[as.character(new$domains)])
summary_df <- data.frame(
  comparison = c("new vs current", "new (split domains collapsed) vs current"),
  ari        = c(ari(domain_new, current_names), ari(new_collapsed, current_names)),
  frac_cells_same_domain = c(NA, mean(new_collapsed == current_names)))
write.csv(summary_df, file.path(OUTPUT_DIR, "comparison_summary.csv"), row.names = FALSE)
cat("\nNew vs current merge:\n")
print(summary_df, digits = 3)

crosstab <- as.data.frame.matrix(table(current = factor(current_names, levels = DOMAIN_ORDER),
                                       new     = factor(domain_new, levels = NEW_ORDER)))
write.csv(crosstab, file.path(OUTPUT_DIR, "domain_crosstab_current_vs_new.csv"))
cat("\nCells: current domain (rows) x new domain (columns):\n")
print(crosstab)

# Per domain: size and immune / tumour / basal-luminal fraction
domain_comp <- bind_rows(
  data.frame(merge = "current", domain = current_names, immune = is_immune, tumour = is_tumour, gland = is_gland),
  data.frame(merge = "new", domain = domain_new, immune = is_immune, tumour = is_tumour, gland = is_gland)) %>%
  group_by(merge, domain) %>%
  summarise(n_cells = n(), immune_frac = mean(immune), tumour_frac = mean(tumour),
            basal_luminal_frac = mean(gland), .groups = "drop") %>%
  arrange(merge, match(domain, union(DOMAIN_ORDER, NEW_ORDER)))
write.csv(domain_comp, file.path(OUTPUT_DIR, "domain_composition_current_vs_new.csv"), row.names = FALSE)
cat("\nComposition per domain:\n")
print(as.data.frame(domain_comp), digits = 3)

# Per sample: cells in each domain, both merges
domain_by_sample <- bind_rows(
  as.data.frame(table(sample = seu$sample, domain = current_names), responseName = "n_cells") %>%
    mutate(merge = "current"),
  as.data.frame(table(sample = seu$sample, domain = domain_new), responseName = "n_cells") %>%
    mutate(merge = "new")) %>%
  group_by(merge, sample) %>% mutate(frac_of_sample = n_cells / sum(n_cells)) %>% ungroup()
write.csv(domain_by_sample, file.path(OUTPUT_DIR, "domain_cells_by_sample_current_vs_new.csv"),
          row.names = FALSE)
cat("\nCells per sample and domain:\n")
print(as.data.frame(tidyr::pivot_wider(
  domain_by_sample %>%
    mutate(col = paste0(merge, ": ", domain)) %>%
    arrange(merge, match(domain, union(DOMAIN_ORDER, NEW_ORDER))) %>%
    select(sample, col, n_cells),
  names_from = col, values_from = n_cells)))

# Per original cluster: does the new merge put it in a different domain?
samp_frac <- prop.table(table(domain_labels, seu$sample), margin = 1)[domains, ]
cluster_cmp <- merge_map %>%
  transmute(cluster = original_domain, n_cells, dissolved,
            domain_new = merged_domain, frac_majority,
            parent_new = unname(parent[match(merged_domain, new_names)])) %>%
  mutate(max_sample_frac = apply(samp_frac[cluster, ], 1, max),
         top_sample      = colnames(samp_frac)[apply(samp_frac[cluster, ], 1, which.max)],
         domain_current  = current_names[match(cluster, domain_labels)],
         changed_parent  = domain_current != parent_new) %>%
  arrange(desc(changed_parent), desc(n_cells))
write.csv(cluster_cmp, file.path(OUTPUT_DIR, "cluster_reassignment_vs_current.csv"), row.names = FALSE)
cat(sprintf("\n%d of %d clusters move to a different current domain (ignoring the split)\n",
            sum(cluster_cmp$changed_parent), nrow(cluster_cmp)))

# --- 11. Seed reproducibility: current vs new merge ----------------------------
# Same settings as 02_01 (Harmony on all PCs with default parameters, SNN on
# 1:BANKSY_NPC_USE, resolution BANKSY_RESOLUTION); only the Harmony seed changes.
# Each clustering is merged both ways (each run uses its own SNN graph for the
# minsize reassignment) and compared (a) pairwise across runs (ARI) and (b) with
# the stored merges domain by domain: every run's domains are matched one-to-one
# to the stored domains (Hungarian) and the cell overlap of each pair is scored
# as Jaccard index.
if (RUN_SEED_CHECK) {
  runs <- list(stored = list(clusters = domain_labels, current = current$domains, new = new$domains))
  for (s in HARMONY_SEEDS) {
    cat(sprintf("Seed check: Harmony + clustering with seed %d...\n", s))
    name <- paste0("seed", s)
    set.seed(s)
    seu <- RunHarmony(seu, group.by.vars = HARMONY_GROUP_VAR, reduction.use = "pca",
                      dims.use = seq_len(BANKSY_NPCS), reduction.save = name, verbose = FALSE)
    set.seed(42)
    seu <- FindNeighbors(seu, reduction = name, dims = 1:BANKSY_NPC_USE,
                         graph.name = paste0(name, c("_nn", "_snn")), verbose = FALSE)
    seu <- FindClusters(seu, graph.name = paste0(name, "_snn"), resolution = BANKSY_RESOLUTION,
                        cluster.name = paste0(name, "_clusters"), verbose = FALSE)
    emb <- Embeddings(seu, name)
    cls <- as.character(seu@meta.data[[paste0(name, "_clusters")]])
    g   <- as(seu[[paste0(name, "_snn")]], "dgCMatrix")
    runs[[name]] <- list(clusters = cls,
                         current  = merge_variant("unweighted", emb, cls, g, k = MERGE_K)$domains,
                         new      = merge_variant(MERGE_METHOD, emb, cls, g, k = K_NEW)$domains)
    seu[[name]] <- NULL
    seu[[paste0(name, "_nn")]] <- NULL
    seu[[paste0(name, "_snn")]] <- NULL
    rm(emb, g); invisible(gc())
  }

  # (a) Pairwise agreement between runs
  pairs <- combn(names(runs), 2, simplify = FALSE)
  seed_cmp <- bind_rows(lapply(pairs, function(p) data.frame(
    run_a = p[1], run_b = p[2],
    n_clusters_a = length(unique(runs[[p[1]]]$clusters)),
    n_clusters_b = length(unique(runs[[p[2]]]$clusters)),
    cluster_ari        = ari(runs[[p[1]]]$clusters, runs[[p[2]]]$clusters),
    domain_ari_current = ari(runs[[p[1]]]$current, runs[[p[2]]]$current),
    domain_ari_new     = ari(runs[[p[1]]]$new, runs[[p[2]]]$new))))
  write.csv(seed_cmp, file.path(OUTPUT_DIR, "seed_reproducibility_current_vs_new.csv"), row.names = FALSE)
  cat("\nPairwise agreement across Harmony seeds:\n")
  print(seed_cmp, digits = 3)
  cat(sprintf("Mean domain ARI over pairs: current %.3f, new %.3f\n",
              mean(seed_cmp$domain_ari_current), mean(seed_cmp$domain_ari_new)))

  # (b) Domain-by-domain agreement with the stored merges
  match_to_ref <- function(x, ref) {
    tab <- table(x, ref)
    asg <- clue::solve_LSAP(unclass(tab), maximum = TRUE)
    unname(setNames(colnames(tab)[as.integer(asg)], rownames(tab))[as.character(x)])
  }
  jaccard <- function(a, b, lv) sapply(lv, function(l) sum(a == l & b == l) / sum(a == l | b == l))
  seed_runs  <- setdiff(names(runs), "stored")
  matched    <- lapply(runs[seed_runs], function(r) list(current = match_to_ref(r$current, current_names),
                                                         new     = match_to_ref(r$new, domain_new)))
  per_domain <- bind_rows(lapply(seed_runs, function(r) bind_rows(
    data.frame(run = r, merge = "current", domain = DOMAIN_ORDER,
               jaccard = jaccard(matched[[r]]$current, current_names, DOMAIN_ORDER)),
    data.frame(run = r, merge = "new", domain = NEW_ORDER,
               jaccard = jaccard(matched[[r]]$new, domain_new, NEW_ORDER)))))
  write.csv(per_domain, file.path(OUTPUT_DIR, "seed_domain_jaccard_vs_stored.csv"), row.names = FALSE)
  cat("\nJaccard overlap of each domain with the stored merge (rows: runs):\n")
  print(as.data.frame(tidyr::pivot_wider(per_domain %>% mutate(col = paste0(merge, ": ", domain)) %>%
                                           select(run, col, jaccard),
                                         names_from = col, values_from = jaccard)), digits = 2)

  # Tumor Core cells per sample in every run (new merge)
  tc_by_run <- sapply(c(list(stored = domain_new), lapply(matched, `[[`, "new")),
                      function(d) tapply(d == "Tumor Core", seu$sample, sum))
  write.csv(tc_by_run, file.path(OUTPUT_DIR, "seed_tumor_core_cells_by_sample.csv"))
  cat("\nTumor Core cells per sample, new merge, each run:\n")
  print(tc_by_run)

  merge_labels <- c(current = sprintf("current\n(k = %d)", MERGE_K),
                    new     = sprintf("%s\n(k = %d)", MERGE_METHOD, K_NEW))
  p_seed <- seed_cmp %>%
    mutate(pair = paste(run_a, "vs", run_b)) %>%
    tidyr::pivot_longer(c(domain_ari_current, domain_ari_new), names_to = "merge", values_to = "ari",
                        names_prefix = "domain_ari_") %>%
    mutate(merge = factor(merge_labels[merge], levels = merge_labels)) %>%
    ggplot(aes(x = merge, y = ari, group = pair)) +
    geom_line(color = "grey60") +
    geom_point(aes(color = merge), size = 2.5) +
    scale_color_manual(values = c("#7570B3", "#1B9E77"), guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "Domain merge", y = "Domain ARI between runs", title = "Pairwise, across Harmony seeds") +
    theme_publication()
  p_dom <- per_domain %>%
    mutate(merge = factor(merge_labels[merge], levels = merge_labels),
           domain = factor(domain, levels = union(DOMAIN_ORDER, NEW_ORDER))) %>%
    ggplot(aes(x = domain, y = jaccard, color = domain)) +
    geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 2) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, color = "black", linewidth = 0.3) +
    scale_color_manual(values = ALL_COLORS, guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    facet_grid(~ merge, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = "Jaccard overlap with stored domain", title = "Per domain, each seed vs stored") +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
  ggsave(file.path(OUTPUT_DIR, "seed_reproducibility_current_vs_new.pdf"),
         p_seed + p_dom + plot_layout(widths = c(1, 2.5)), width = 12, height = 5.5)
}

cat("Outputs saved to:", OUTPUT_DIR, "\n")
cat("Done.\n")
