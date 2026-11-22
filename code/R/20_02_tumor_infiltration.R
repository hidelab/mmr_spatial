# 20_02_tumor_infiltration.R
# Exploratory: CD8+ T cell infiltration into tumor, measured relative to the
# tumor edge so that it does not depend on tumor size.
#
# Tumor region: cells whose 19_01 tumor score (heat-diffusion local tumor
# fraction, in [0,1], comparable across slides) is >= TUMOR_REGION_THRESHOLD.
# Signed distance to the tumor edge, per sample: for a cell outside the region,
# + distance to the nearest region cell; for a cell inside, - distance to the
# nearest non-region cell.
#
# Metrics for CD8+ T subtypes (03_01 clusters) and all QC-pass CD8+ T cells:
#   depth profile        fraction of cells that belong to the group in bins of
#                        signed distance, divided by the group's slide-wide
#                        fraction (1 = slide average)
#   infiltration ratio   log2[ fraction of group among cells within BAND um
#                        inside the edge / same within BAND um outside ];
#                        > 0 = denser inside the tumor border than just outside.
#                        Only the border zones are compared, so large tumor
#                        interiors do not dilute it.
#   deep fraction        share of the group's border-zone cells lying deeper
#                        than BAND / 2 inside the tumor
# Per-sample tests treat cells as independent (Fisher's exact; approximate).
# DKO vs WT: Welch t-test and exact permutation test over condition labels,
# samples as replicates. Thresholds in TUMOR_REGION_SENSITIVITY are re-run for
# the infiltration ratio.
#
# Reads CSV outputs only; no RDS is read or written.

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(FNN)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "20_tumor_infiltration")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette and helpers
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

# =============================================================================
# Parameters
# =============================================================================
TUMOR_CORE_CSV  <- file.path(PROJECT_DIR, "output", "R", "19_tumor_core_domain",
                             "diffusion_filled", "tumor_core_assignment.csv")   # from 19_01
TCELL_ASSIGN_CSV <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters",
                              "tcell_cluster_assignments.csv")                  # from 03_01
TCELL_CLUSTER_COL <- "tcell_clust_res.0.5"
GROUP_CLUSTERS    <- c("2", "1")      # subtypes analysed (plus all CD8+ T)
TUMOR_REGION_THRESHOLD   <- 0.5
TUMOR_REGION_SENSITIVITY <- c(0.2, 0.75)
BAND        <- 100     # um, border band on each side of the tumor edge
DEPTH_MAX   <- 300     # um, range of the depth profile
DEPTH_BIN   <- 25      # um
MIN_BAND_CELLS <- 150  # min cells in each border band to report a sample
MIN_GROUP_CELLS <- 10  # min group cells across both bands to report a ratio

GROUPS <- c(paste0("CD8T_C", GROUP_CLUSTERS), "CD8T_all")

# =============================================================================
# 1. Load inputs
# =============================================================================
stopifnot(file.exists(TUMOR_CORE_CSV), file.exists(TCELL_ASSIGN_CSV))
writeLines(c(sprintf("20_02 tumor infiltration, run %s", format(Sys.time())),
             sprintf("tumor score source: %s (modified %s)", TUMOR_CORE_CSV,
                     format(file.mtime(TUMOR_CORE_CSV))),
             sprintf("TCELL_CLUSTER_COL = %s; GROUPS = %s", TCELL_CLUSTER_COL, paste(GROUPS, collapse = ", ")),
             sprintf("TUMOR_REGION_THRESHOLD = %g; sensitivity = %s", TUMOR_REGION_THRESHOLD,
                     paste(TUMOR_REGION_SENSITIVITY, collapse = ", ")),
             sprintf("BAND = %g; DEPTH_MAX = %g; DEPTH_BIN = %g; MIN_BAND_CELLS = %d; MIN_GROUP_CELLS = %d",
                     BAND, DEPTH_MAX, DEPTH_BIN, MIN_BAND_CELLS, MIN_GROUP_CELLS)),
           file.path(OUTPUT_DIR, "run_parameters.txt"))

cells <- read.csv(TUMOR_CORE_CSV, stringsAsFactors = FALSE)
cells$sample <- factor(cells$sample, levels = intersect(names(SAMPLE_COLORS), unique(cells$sample)))
cells$condition <- factor(cells$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
samples <- levels(cells$sample)

tc <- read.csv(TCELL_ASSIGN_CSV, stringsAsFactors = FALSE, check.names = FALSE)
cells$cd8_subtype <- NA_character_
cells$cd8_subtype[match(tc$cell, cells$cell)] <- paste0("CD8T_C", tc[[TCELL_CLUSTER_COL]])
stopifnot(sum(!is.na(cells$cd8_subtype)) == nrow(tc))
in_group <- function(g) if (g == "CD8T_all") !is.na(cells$cd8_subtype) else cells$cd8_subtype %in% g
for (g in GROUPS) cells[[g]] <- in_group(g)
cat(sprintf("Cells: %d; QC-pass CD8+ T with subtype: %d\n", nrow(cells), nrow(tc)))

# =============================================================================
# 2. Signed distance to the tumor edge
# =============================================================================
signed_distance <- function(threshold) {
  sd <- rep(NA_real_, nrow(cells))
  for (s in samples) {
    idx <- which(cells$sample == s)
    inside <- cells$tumor_score[idx] >= threshold
    if (sum(inside) == 0 || sum(!inside) == 0) next
    xy <- as.matrix(cells[idx, c("x", "y")])
    d_out <- get.knnx(xy[inside, , drop = FALSE], xy[!inside, , drop = FALSE], k = 1)$nn.dist[, 1]
    d_in  <- get.knnx(xy[!inside, , drop = FALSE], xy[inside, , drop = FALSE], k = 1)$nn.dist[, 1]
    sd[idx[!inside]] <- d_out
    sd[idx[inside]]  <- -d_in
  }
  sd
}
cells$signed_dist <- signed_distance(TUMOR_REGION_THRESHOLD)

region_summary <- cells %>% group_by(sample, condition) %>%
  summarise(n_cells = n(), n_tumor_region = sum(signed_dist < 0, na.rm = TRUE),
            n_inside_band = sum(signed_dist < 0 & signed_dist >= -BAND, na.rm = TRUE),
            n_outside_band = sum(signed_dist >= 0 & signed_dist <= BAND, na.rm = TRUE),
            .groups = "drop")
write.csv(region_summary, file.path(OUTPUT_DIR, "tumor_region_summary.csv"), row.names = FALSE)
cat("\nTumor region and border bands per sample:\n"); print(region_summary)

# =============================================================================
# 3. Infiltration ratio per sample (border bands)
# =============================================================================
infiltration_table <- function(sdist, threshold) {
  bind_rows(lapply(samples, function(s) bind_rows(lapply(GROUPS, function(g) {
    i <- cells$sample == s & !is.na(sdist)
    ins <- i & sdist < 0 & sdist >= -BAND
    out <- i & sdist >= 0 & sdist <= BAND
    deep <- i & sdist < -BAND / 2 & sdist >= -BAND
    gi <- sum(cells[[g]][ins]); go <- sum(cells[[g]][out])
    ok <- sum(ins) >= MIN_BAND_CELLS && sum(out) >= MIN_BAND_CELLS && gi + go >= MIN_GROUP_CELLS
    ft <- if (ok) fisher.test(matrix(c(gi, sum(ins) - gi, go, sum(out) - go), 2)) else NULL
    data.frame(threshold = threshold, sample = s, condition = as.character(cells$condition[i][1]),
               group = g, n_inside_band = sum(ins), n_outside_band = sum(out),
               group_inside = gi, group_outside = go,
               frac_inside = gi / max(sum(ins), 1), frac_outside = go / max(sum(out), 1),
               log2_infiltration_ratio = if (ok) log2(((gi + 0.5) / sum(ins)) / ((go + 0.5) / sum(out))) else NA,
               deep_fraction = if (ok) sum(cells[[g]][deep]) / max(gi + go, 1) else NA,
               p_fisher = if (ok) ft$p.value else NA)
  }))))
}

infil <- infiltration_table(cells$signed_dist, TUMOR_REGION_THRESHOLD) %>%
  group_by(group) %>% mutate(padj_fisher = p.adjust(p_fisher, "fdr")) %>% ungroup()
write.csv(infil, file.path(OUTPUT_DIR, "infiltration_ratio_per_sample.csv"), row.names = FALSE)
cat("\nInfiltration ratio (log2 inside band / outside band):\n")
print(as.data.frame(infil %>% select(group, sample, n_inside_band, group_inside, group_outside,
                                     log2_infiltration_ratio, padj_fisher)), digits = 3)

# DKO vs WT, samples as replicates
cond_test <- function(v, cond) {
  ok <- is.finite(v); v <- v[ok]; cond <- cond[ok]
  dko <- cond == TEST_CONDITION
  if (sum(dko) < 2 || sum(!dko) < 2) return(data.frame(n_WT = sum(!dko), n_DKO = sum(dko)))
  diff_fun <- function(g) mean(v[g]) - mean(v[!g])
  d_obs <- diff_fun(dko)
  d_perm <- apply(combn(length(v), sum(!dko)), 2, function(w) {
    g <- rep(TRUE, length(v)); g[w] <- FALSE; diff_fun(g) })
  data.frame(n_WT = sum(!dko), n_DKO = sum(dko), mean_WT = mean(v[!dko]), mean_DKO = mean(v[dko]),
             diff_DKO_minus_WT = d_obs, p_welch = t.test(v[dko], v[!dko])$p.value,
             p_exact_perm = mean(abs(d_perm) >= abs(d_obs) - 1e-12))
}

sens <- bind_rows(lapply(TUMOR_REGION_SENSITIVITY, function(th)
  infiltration_table(signed_distance(th), th)))
infil_all <- bind_rows(infil %>% select(-padj_fisher), sens)
cond_df <- bind_rows(lapply(split(infil_all, list(infil_all$threshold, infil_all$group), drop = TRUE),
  function(d) bind_rows(
    data.frame(threshold = d$threshold[1], group = d$group[1], metric = "log2_infiltration_ratio",
               cond_test(d$log2_infiltration_ratio, d$condition)),
    data.frame(threshold = d$threshold[1], group = d$group[1], metric = "deep_fraction",
               cond_test(d$deep_fraction, d$condition)))))
write.csv(sens, file.path(OUTPUT_DIR, "infiltration_ratio_sensitivity_thresholds.csv"), row.names = FALSE)
write.csv(cond_df, file.path(OUTPUT_DIR, "infiltration_DKOvsWT.csv"), row.names = FALSE)
cat("\nDKO vs WT (samples as replicates):\n")
print(cond_df %>% arrange(metric, group, threshold), digits = 3)

# =============================================================================
# 4. Depth profiles
# =============================================================================
breaks <- seq(-DEPTH_MAX, DEPTH_MAX, by = DEPTH_BIN)
cells$depth_bin <- cut(cells$signed_dist, breaks, right = FALSE)
bin_mid <- setNames(head(breaks, -1) + DEPTH_BIN / 2, levels(cells$depth_bin))

slide_frac <- cells %>% group_by(sample) %>%
  summarise(across(all_of(GROUPS), mean), .groups = "drop") %>%
  pivot_longer(all_of(GROUPS), names_to = "group", values_to = "slide_frac")
profile <- cells %>% filter(!is.na(depth_bin)) %>%
  group_by(sample, condition, depth_bin) %>%
  summarise(n_cells = n(), across(all_of(GROUPS), sum), .groups = "drop") %>%
  pivot_longer(all_of(GROUPS), names_to = "group", values_to = "n_group") %>%
  left_join(slide_frac, by = c("sample", "group")) %>%
  mutate(depth = bin_mid[as.character(depth_bin)],
         frac = n_group / n_cells, enrichment = frac / slide_frac)
write.csv(profile, file.path(OUTPUT_DIR, "depth_profile.csv"), row.names = FALSE)

# =============================================================================
# 5. Plots
# =============================================================================
sig_label <- function(p) case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
                                   p < 0.05 ~ "*", p < FDR_THRESHOLD ~ "^", TRUE ~ "")
group_labs <- c(setNames(paste0("C", GROUP_CLUSTERS), paste0("CD8T_C", GROUP_CLUSTERS)),
                CD8T_all = "all CD8+ T")

# 5a. Depth profiles: one line per sample, mean per condition
prof_plot <- profile %>% filter(n_cells >= 50) %>%
  mutate(group = factor(group, levels = GROUPS, labels = group_labs[GROUPS]),
         sample = factor(sample, levels = samples))
prof_mean <- prof_plot %>% group_by(condition, group, depth) %>%
  summarise(enrichment = mean(enrichment), .groups = "drop")
p_prof <- ggplot(prof_plot, aes(x = depth, y = enrichment, color = condition)) +
  annotate("rect", xmin = -DEPTH_MAX, xmax = 0, ymin = -Inf, ymax = Inf, fill = "#F4C7C3", alpha = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, color = "grey30") +
  geom_line(aes(group = sample), alpha = 0.35, linewidth = 0.5) +
  geom_line(data = prof_mean, linewidth = 1.3) +
  facet_wrap(~ group, nrow = 1, scales = "free_y") +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "CD8+ T density across the tumor edge",
       subtitle = sprintf(paste0("Shaded = inside tumor region (local tumor fraction >= %.2f); ",
                                 "thin = per sample, thick = condition mean; 1 = slide average"),
                          TUMOR_REGION_THRESHOLD),
       x = "Signed distance to tumor edge (um; negative = inside tumor)",
       y = "Density relative to slide average", color = NULL) +
  theme_publication()
save_pdf_png(p_prof, file.path(OUTPUT_DIR, "depth_profile_tumor_edge"), width = 12, height = 4.5)

# 5b. Infiltration ratio per sample by condition
lab_c <- cond_df %>% filter(threshold == TUMOR_REGION_THRESHOLD, metric == "log2_infiltration_ratio") %>%
  mutate(group = factor(group, levels = GROUPS, labels = group_labs[GROUPS]),
         label = sprintf("Welch p=%.2g\nexact p=%.2g", p_welch, p_exact_perm))
d_ir <- infil %>% filter(!is.na(log2_infiltration_ratio)) %>%
  mutate(group = factor(group, levels = GROUPS, labels = group_labs[GROUPS]))
p_ir <- ggplot(d_ir, aes(x = condition, y = log2_infiltration_ratio)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(aes(fill = condition), alpha = 0.4, outlier.shape = NA, width = 0.55) +
  geom_point(aes(color = condition), size = 2.3) +
  geom_text(aes(label = paste0(sample, sig_label(padj_fisher))), nudge_x = 0.3, size = 2.6) +
  geom_text(data = lab_c, aes(x = 1.5, y = Inf, label = label), vjust = 1.2, size = 2.8) +
  facet_wrap(~ group, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.3))) +
  scale_fill_manual(values = CONDITION_COLORS) + scale_color_manual(values = CONDITION_COLORS) +
  labs(title = sprintf("Infiltration ratio: density within %g um inside vs outside the tumor edge", BAND),
       subtitle = "> 0 = denser inside the tumor border; one point per sample; stars after names = per-sample Fisher FDR",
       x = NULL, y = "log2 infiltration ratio") +
  theme_publication() + theme(legend.position = "none")
save_pdf_png(p_ir, file.path(OUTPUT_DIR, "infiltration_ratio_by_condition"), width = 10, height = 4.8)

# 5c. Sensitivity of the infiltration ratio to the tumor-region threshold
d_sens <- infil_all %>% filter(!is.na(log2_infiltration_ratio)) %>%
  mutate(group = factor(group, levels = GROUPS, labels = group_labs[GROUPS]),
         sample = factor(sample, levels = samples))
p_sens <- ggplot(d_sens, aes(x = factor(threshold), y = log2_infiltration_ratio, group = sample, color = condition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(alpha = 0.6) + geom_point(size = 1.8) +
  facet_wrap(~ group, nrow = 1) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "Infiltration ratio vs tumor-region threshold", x = "Local tumor fraction threshold",
       y = "log2 infiltration ratio", color = NULL) +
  theme_publication()
save_pdf_png(p_sens, file.path(OUTPUT_DIR, "infiltration_ratio_threshold_sensitivity"), width = 10, height = 4)

# 5d. Spatial maps: tumor region, border bands and CD8+ T groups
max_span <- max(tapply(cells$x, cells$sample, function(v) diff(range(v))),
                tapply(cells$y, cells$sample, function(v) diff(range(v))))
half_span <- max_span / 2 * 1.05
plots_map <- lapply(samples, function(s) {
  d <- cells[cells$sample == s, ]
  d$zone <- case_when(is.na(d$signed_dist) ~ "outside",
                      d$signed_dist < -BAND ~ "tumor interior",
                      d$signed_dist < 0 ~ "inner band",
                      d$signed_dist <= BAND ~ "outer band",
                      TRUE ~ "outside")
  tcd <- d[!is.na(d$cd8_subtype) & d$cd8_subtype %in% paste0("CD8T_C", GROUP_CLUSTERS), ]
  p <- ggplot() +
    geom_point(data = d, aes(x, y, color = zone), size = 0.03) +
    scale_color_manual(values = c("outside" = "grey92", "outer band" = "#C6DBEF",
                                  "inner band" = "#FCBBA1", "tumor interior" = "#FB6A4A"),
                       name = NULL)
  p <- p + geom_point(data = tcd, aes(x, y, shape = cd8_subtype), size = 0.5, color = "black") +
    scale_shape_manual(values = setNames(c(16, 17)[seq_along(GROUP_CLUSTERS)],
                                         paste0("CD8T_C", GROUP_CLUSTERS)), name = NULL) +
    ggtitle(s) + theme_void() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (CENTER_COORDS) p + coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span))
  else p + coord_fixed()
})
p_maps <- wrap_plots(plots_map, nrow = 1, guides = "collect") &
  guides(color = guide_legend(override.aes = list(size = 3)), shape = guide_legend(override.aes = list(size = 2.5)))
ggsave(file.path(OUTPUT_DIR, "spatial_tumor_region_bands.png"), p_maps,
       width = 4 * length(samples) + 1.5, height = 5, dpi = 250)

cat("\nOutputs saved to:", OUTPUT_DIR, "\nDone.\n")
