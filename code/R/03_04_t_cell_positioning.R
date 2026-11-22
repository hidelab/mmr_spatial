# 03_04_t_cell_positioning.R
# Positioning of CD8+ T cell subtypes (clusters from 03_01) relative to the tumor
# core, and whether it differs between DKO and WT.
#
# Tumor proximity comes from 19_01_tumor_core_domain.R (per-cell
# tumor_core_assignment.csv): the heat-diffusion local tumor fraction
# ("tumor_score", in [0,1], comparable across slides) and the score-defined
# tumor core ("tumor_zone"). Metrics per CD8+ T cell:
#   score      log10(tumor_score + SCORE_EPS); higher = closer to tumor
#   score_pct  percentile of tumor_score among all non-tumor cells of the same
#              sample; 0.5 = positioned like a typical non-tumor cell. Removes
#              between-slide differences in tumor burden.
#   dist_core  distance (spatial units, microns) to the nearest tumor-core cell,
#              0 inside the core; NA in samples with < MIN_CORE_CELLS core cells
#   in_core    cell lies in the tumor core
#
# Statistics are run on sample-level summaries (one value per sample x subtype),
# so the 3 WT + 4 DKO samples, not cells, are the replicates:
#   (1) offset = subtype summary - all-CD8 summary of the same sample.
#       offset != 0 (t-test across samples): subtype sits closer / farther than
#       CD8+ T cells overall. offset DKO vs WT (Welch t-test): positioning changes
#       between conditions, independently of overall CD8 / tumor shifts.
#   (2) absolute summary DKO vs WT per subtype (Welch t-test).
#   (3) global mixed model: summary ~ subtype * condition + (1 | sample).
# P values are FDR-adjusted across subtypes within each metric and test.
#
# Reads CSV outputs only; no RDS is read or written.

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(FNN)
library(lme4)
library(lmerTest)

# Source centralized parameters (defines PROJECT_DIR and all shared constants)
source("./code/R/00_parameters.R")

TCELL_DIR  <- file.path(PROJECT_DIR, "output", "R", "03_t_cell_clusters")
OUTPUT_DIR <- file.path(TCELL_DIR, "positioning")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Source shared color palette and helpers
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))
source(file.path(PROJECT_DIR, "code", "R", "00_utils.R"))

# =============================================================================
# Parameters
# =============================================================================
TCELL_ASSIGN_CSV       <- file.path(TCELL_DIR, "tcell_cluster_assignments.csv")   # from 03_01
TUMOR_CORE_DIR         <- file.path(PROJECT_DIR, "output", "R", "19_tumor_core_domain",
                                    "diffusion_filled")                            # from 19_01
TCELL_CLUSTER_COL      <- "tcell_clust_res.0.4"
TCELL_EXCLUDE_CLUSTERS <- character(0)   # e.g. "C0" to drop a contaminated cluster
TARGET_CELLTYPE        <- "Tumor"
ZONE_CORE              <- "Tumor Core"
SCORE_EPS              <- 1e-3   # floor for log10(tumor_score)
MIN_CORE_CELLS         <- 500    # min core cells in a sample to compute dist_core
SUMMARY_FUN            <- median # per sample x subtype summary of continuous metrics

METRICS <- c(score     = "log10 tumor score (local tumor fraction)",
             score_pct = "Tumor score percentile (within sample)",
             dist_core = "Distance to tumor core",
             in_core   = "Fraction in tumor core")

# =============================================================================
# 1. Load inputs
# =============================================================================
tumor_csv <- file.path(TUMOR_CORE_DIR, "tumor_core_assignment.csv")
stopifnot(file.exists(TCELL_ASSIGN_CSV), file.exists(tumor_csv))
cat("Tumor score source:", tumor_csv, "\n  modified:",
    format(file.mtime(tumor_csv)), "\n")

# Record the 19_01 run that produced the score (19_01 may be re-run later)
params_19 <- file.path(TUMOR_CORE_DIR, "run_parameters.txt")
if (file.exists(params_19)) {
  writeLines(c(sprintf("source: %s (modified %s)", tumor_csv, format(file.mtime(tumor_csv))),
               "", readLines(params_19)),
             file.path(OUTPUT_DIR, "tumor_score_source_19_01.txt"))
}

all_cells <- read.csv(tumor_csv, stringsAsFactors = FALSE)
all_cells$sample <- factor(all_cells$sample,
                           levels = intersect(names(SAMPLE_COLORS), unique(all_cells$sample)))
cat("All cells:", nrow(all_cells), "\n")

tc <- read.csv(TCELL_ASSIGN_CSV, stringsAsFactors = FALSE, check.names = FALSE)
SUBTYPE_ORDER <- paste0("C", sort(unique(tc[[TCELL_CLUSTER_COL]])))
CLUSTER_COLORS <- tcell_subtype_colors(SUBTYPE_ORDER)
tc$subtype <- factor(paste0("C", tc[[TCELL_CLUSTER_COL]]), levels = SUBTYPE_ORDER)
tc <- tc[!tc$subtype %in% TCELL_EXCLUDE_CLUSTERS, ]
tc$subtype <- droplevels(tc$subtype)
subtypes <- levels(tc$subtype)

cd8 <- merge(tc[, c("cell", "subtype")],
             all_cells[, c("cell", "sample", "condition", "x", "y", "tumor_score", "tumor_zone")],
             by = "cell")
cat(sprintf("CD8+ T cells: %d in 03_01 assignments, %d matched to 19_01 scores\n",
            nrow(tc), nrow(cd8)))
stopifnot(nrow(cd8) == nrow(tc))
cd8$condition <- factor(cd8$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))

# =============================================================================
# 2. Per-cell positioning metrics
# =============================================================================
cd8$score   <- log10(cd8$tumor_score + SCORE_EPS)
cd8$in_core <- as.numeric(cd8$tumor_zone == ZONE_CORE)
cd8$score_pct <- NA_real_
cd8$dist_core <- NA_real_

core_n <- table(all_cells$sample[all_cells$tumor_zone == ZONE_CORE])
core_n <- setNames(as.integer(core_n[levels(all_cells$sample)]), levels(all_cells$sample))
core_n[is.na(core_n)] <- 0L
cat("\nTumor core cells per sample:\n"); print(core_n)

for (s in levels(all_cells$sample)) {
  bg <- all_cells[all_cells$sample == s, ]
  idx <- which(cd8$sample == s)
  if (length(idx) == 0) next

  # Percentile among all non-tumor cells of the sample
  ref_ecdf <- ecdf(bg$tumor_score[bg$celltype != TARGET_CELLTYPE])
  cd8$score_pct[idx] <- ref_ecdf(cd8$tumor_score[idx])

  # Distance to nearest tumor-core cell (0 inside the core)
  if (core_n[[s]] >= MIN_CORE_CELLS) {
    core_xy <- as.matrix(bg[bg$tumor_zone == ZONE_CORE, c("x", "y")])
    nn <- get.knnx(core_xy, as.matrix(cd8[idx, c("x", "y")]), k = 1)
    cd8$dist_core[idx] <- ifelse(cd8$in_core[idx] == 1, 0, nn$nn.dist[, 1])
  }
}
dist_samples <- names(core_n)[core_n >= MIN_CORE_CELLS]
cat(sprintf("dist_core computed for samples with >= %d core cells: %s\n",
            MIN_CORE_CELLS, paste(dist_samples, collapse = ", ")))

write.csv(cd8, file.path(OUTPUT_DIR, "CD8T_positioning_per_cell.csv"), row.names = FALSE)

# =============================================================================
# 3. Sample-level summaries and offsets
# =============================================================================
summarise_metric <- function(v, metric) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return(NA_real_)
  if (metric == "in_core") mean(v) else SUMMARY_FUN(v)
}

long_cells <- cd8 %>%
  select(cell, sample, condition, subtype, all_of(names(METRICS))) %>%
  pivot_longer(all_of(names(METRICS)), names_to = "metric", values_to = "value")

samp_sub <- long_cells %>%
  group_by(metric, sample, condition, subtype) %>%
  summarise(n_cells = n(), value = summarise_metric(value, first(metric)), .groups = "drop")
samp_all <- long_cells %>%
  group_by(metric, sample, condition) %>%
  summarise(all_cd8 = summarise_metric(value, first(metric)), .groups = "drop")

samp_sub <- samp_sub %>%
  left_join(samp_all, by = c("metric", "sample", "condition")) %>%
  mutate(value = ifelse(n_cells >= MIN_CELLS_PER_SAMPLE, value, NA_real_),
         offset = value - all_cd8)
write.csv(samp_sub, file.path(OUTPUT_DIR, "CD8T_positioning_per_sample_subtype.csv"), row.names = FALSE)

# =============================================================================
# 4. Statistics
# =============================================================================
safe_t <- function(x, y = NULL) {
  x <- x[!is.na(x)]
  if (is.null(y)) {
    if (length(x) < 3) return(c(est = NA, p = NA))
    tt <- t.test(x)
    return(c(est = unname(tt$estimate), p = tt$p.value))
  }
  y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return(c(est = NA, p = NA))
  tt <- t.test(y, x)   # y - x
  c(est = unname(tt$estimate[1] - tt$estimate[2]), p = tt$p.value)
}

stats_df <- samp_sub %>%
  group_by(metric, subtype) %>%
  group_modify(function(d, key) {
    wt  <- d$condition == REFERENCE_CONDITION
    off_all <- safe_t(d$offset)
    off_cnd <- safe_t(d$offset[wt], d$offset[!wt])
    abs_cnd <- safe_t(d$value[wt], d$value[!wt])
    data.frame(
      n_samples_WT = sum(!is.na(d$value[wt])), n_samples_DKO = sum(!is.na(d$value[!wt])),
      mean_WT = mean(d$value[wt], na.rm = TRUE), mean_DKO = mean(d$value[!wt], na.rm = TRUE),
      offset_mean = off_all[["est"]], offset_p = off_all[["p"]],
      offset_mean_WT = mean(d$offset[wt], na.rm = TRUE),
      offset_mean_DKO = mean(d$offset[!wt], na.rm = TRUE),
      offset_diff_DKOvsWT = off_cnd[["est"]], offset_diff_p = off_cnd[["p"]],
      abs_diff_DKOvsWT = abs_cnd[["est"]], abs_diff_p = abs_cnd[["p"]]
    )
  }) %>%
  group_by(metric) %>%
  mutate(offset_padj      = p.adjust(offset_p, method = "fdr"),
         offset_diff_padj = p.adjust(offset_diff_p, method = "fdr"),
         abs_diff_padj    = p.adjust(abs_diff_p, method = "fdr")) %>%
  ungroup()
write.csv(stats_df, file.path(OUTPUT_DIR, "CD8T_positioning_stats.csv"), row.names = FALSE)

cat("\n== Subtype offsets vs all CD8+ T (positive = closer to tumor for score / pct / in_core,",
    "farther for dist_core) ==\n")
print(as.data.frame(stats_df %>%
  select(metric, subtype, offset_mean, offset_padj, offset_diff_DKOvsWT, offset_diff_padj,
         abs_diff_DKOvsWT, abs_diff_padj)), digits = 3)

# Global mixed model on sample-level summaries
global_df <- bind_rows(lapply(names(METRICS), function(m) {
  d <- samp_sub %>% filter(metric == m, !is.na(value))
  if (length(unique(d$sample)) < 4) return(NULL)
  fit <- lmer(value ~ subtype * condition + (1 | sample), data = d)
  a <- as.data.frame(anova(fit))
  data.frame(metric = m, term = rownames(a), F = a$`F value`,
             df_num = a$NumDF, df_den = a$DenDF, p = a$`Pr(>F)`)
}))
write.csv(global_df, file.path(OUTPUT_DIR, "CD8T_positioning_mixed_model_anova.csv"), row.names = FALSE)
cat("\n== Mixed model: summary ~ subtype * condition + (1 | sample) ==\n")
print(global_df, digits = 3)

# =============================================================================
# 5. Plots
# =============================================================================
sig_label <- function(p) {
  case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
            p < 0.05 ~ "*", p < FDR_THRESHOLD ~ "^", TRUE ~ "")
}
metric_lab <- function(m) {
  if (m == "dist_core") paste0(METRICS[[m]], " (", paste(dist_samples, collapse = ", "), ")")
  else METRICS[[m]]
}

# 5a. Per-cell distributions by subtype and condition, with sample summaries
for (m in setdiff(names(METRICS), "in_core")) {
  dc <- long_cells %>% filter(metric == m, !is.na(value))
  ds <- samp_sub %>% filter(metric == m, !is.na(value))
  lab <- stats_df %>% filter(metric == m) %>%
    mutate(label = sig_label(abs_diff_padj), y = max(dc$value) * 1.02 + 0.02 * diff(range(dc$value)))
  p <- ggplot(dc, aes(x = subtype, y = value, fill = condition)) +
    geom_violin(scale = "width", alpha = 0.5, color = NA, position = position_dodge(0.8)) +
    geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(0.8), alpha = 0.8) +
    geom_point(data = ds, aes(color = condition), position = position_dodge(0.8),
               size = 2, shape = 21, fill = "white", stroke = 1) +
    geom_text(data = lab, aes(x = subtype, y = y, label = label), inherit.aes = FALSE, size = 5) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(title = metric_lab(m),
         subtitle = "Circles = per-sample median; stars = DKO vs WT on sample medians (FDR)",
         x = NULL, y = METRICS[[m]], fill = NULL, color = NULL) +
    theme_publication() + theme(legend.position = "top")
  save_pdf_png(p, file.path(OUTPUT_DIR, paste0("violin_", m, "_by_subtype")), width = 7, height = 5)
}

# 5b. Offsets vs all CD8+ T per sample
p_off <- lapply(names(METRICS), function(m) {
  d <- samp_sub %>% filter(metric == m, !is.na(offset))
  lab <- stats_df %>% filter(metric == m) %>%
    mutate(label = paste0(sig_label(offset_padj),
                          ifelse(sig_label(offset_diff_padj) == "", "",
                                 paste0("\n", sig_label(offset_diff_padj), " cond"))),
           y = max(d$offset) + 0.1 * diff(range(d$offset)))
  ggplot(d, aes(x = subtype, y = offset)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_boxplot(aes(fill = condition), outlier.shape = NA, alpha = 0.6,
                 position = position_dodge(0.8)) +
    geom_point(aes(color = condition), position = position_dodge(0.8), size = 1.8) +
    geom_text(data = lab, aes(x = subtype, y = y, label = label), size = 3.5, lineheight = 0.8) +
    scale_fill_manual(values = CONDITION_COLORS) +
    scale_color_manual(values = CONDITION_COLORS) +
    labs(title = metric_lab(m), x = NULL, y = "Subtype - all CD8+ T", fill = NULL, color = NULL) +
    theme_publication(base_size = 10)
})
save_pdf_png(wrap_plots(p_off, ncol = 2, guides = "collect") +
               plot_annotation(
                 title = "CD8+ T subtype positioning relative to all CD8+ T cells (per sample)",
                 subtitle = "Stars: offset != 0 across samples; 'cond': offset differs DKO vs WT (FDR)"),
             file.path(OUTPUT_DIR, "offset_vs_all_CD8T"), width = 11, height = 8)

# 5c. ECDF of within-sample percentile per subtype
p_ecdf <- ggplot(cd8, aes(x = score_pct, color = condition)) +
  stat_ecdf(aes(group = sample), alpha = 0.35, linewidth = 0.4) +
  stat_ecdf(linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  facet_wrap(~ subtype, nrow = 1) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(title = "Tumor-score percentile of CD8+ T subtypes",
       subtitle = "Thin = per sample, thick = pooled; diagonal = typical non-tumor cell; curves to the right = closer to tumor",
       x = "Tumor score percentile (within sample)", y = "Cumulative fraction", color = NULL) +
  theme_publication()
save_pdf_png(p_ecdf, file.path(OUTPUT_DIR, "ecdf_score_pct_by_subtype"),
             width = 3 * length(subtypes) + 1, height = 3.8)

# 5d. Fraction in tumor core per sample x subtype
d_core <- samp_sub %>% filter(metric == "in_core", !is.na(value))
p_core <- ggplot(d_core, aes(x = sample, y = value, fill = subtype)) +
  geom_col(position = position_dodge(0.85), width = 0.8) +
  scale_fill_manual(values = CLUSTER_COLORS) +
  labs(title = "Fraction of each CD8+ T subtype inside the tumor core",
       subtitle = sprintf("Tumor core cells per sample: %s",
                          paste(names(core_n), core_n, sep = "=", collapse = ", ")),
       x = NULL, y = "Fraction in tumor core", fill = "Subtype") +
  theme_publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf_png(p_core, file.path(OUTPUT_DIR, "fraction_in_core_per_sample"), width = 9, height = 4.5)

# 5e. Spatial maps: tumor core + CD8+ T subtypes
max_span <- max(tapply(all_cells$x, all_cells$sample, function(v) diff(range(v))),
                tapply(all_cells$y, all_cells$sample, function(v) diff(range(v))))
half_span <- max_span / 2 * 1.05
plots_spatial <- lapply(levels(all_cells$sample), function(s) {
  bg <- all_cells[all_cells$sample == s, ]
  p <- ggplot() +
    geom_point(data = bg, aes(x, y), color = "grey90", size = 0.03) +
    geom_point(data = bg[bg$tumor_zone == ZONE_CORE, ], aes(x, y),
               color = "#F4C7C3", size = 0.03) +
    geom_point(data = cd8[cd8$sample == s, ], aes(x, y, color = subtype), size = 0.35) +
    scale_color_manual(values = CLUSTER_COLORS, drop = FALSE) +
    ggtitle(s) + theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  if (CENTER_COORDS) {
    p + coord_fixed(xlim = c(-half_span, half_span), ylim = c(-half_span, half_span))
  } else {
    p + coord_fixed()
  }
})
p_spatial <- wrap_plots(plots_spatial, nrow = 1, guides = "collect") &
  guides(color = guide_legend(override.aes = list(size = 3)))
ggsave(file.path(OUTPUT_DIR, "spatial_CD8T_subtypes_tumor_core.png"), p_spatial,
       width = 4 * length(plots_spatial) + 1, height = 5, dpi = 300)

# =============================================================================
# Done
# =============================================================================
cat("\n========================================\n")
cat("CD8+ T subtype positioning complete!\n")
cat("========================================\n")
cat(sprintf("Clustering: %s; subtypes: %s\n", TCELL_CLUSTER_COL, paste(subtypes, collapse = ", ")))
cat("Output directory:", OUTPUT_DIR, "\n")
