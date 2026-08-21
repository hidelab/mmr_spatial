# 14_03_knn_purity_by_domain.R
# Combine k-NN tumor purity with spatial domain assignments.
# Ask: does tumor purity differ across domains, and specifically between
# Tumor Core vs. Immune Engulfing? Does condition modulate this?

library(Seurat)
library(ggplot2)
library(patchwork)
library(FNN)
library(dplyr)
library(glmmTMB)

# Source centralized parameters
source("./code/R/00_parameters.R")
source(file.path(PROJECT_DIR, "code", "R", "00_color_palette.R"))

OUTPUT_DIR <- file.path(PROJECT_DIR, "output", "R", "14_knn_purity_by_domain")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Parameters
# =============================================================================
K_NEIGHBORS <- 30
TARGET_CELLTYPE <- "Tumor"
DOMAINS_OF_INTEREST <- c("Tumor Core", "Immune Engulfing")

# =============================================================================
# Load data
# =============================================================================
seu <- readRDS(file.path(PROJECT_DIR, "output", "R", "DKO_banksy_merged.rds"))
cat("Loaded:", ncol(seu), "cells\n")

seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
seu$condition <- factor(seu$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
seu$sample <- factor(seu$sample,
                     levels = intersect(names(SAMPLE_COLORS), unique(seu$sample)))

# Domain names
if ("banksy_domain_name" %in% colnames(seu@meta.data)) {
  seu$domain_name <- seu$banksy_domain_name
} else {
  seu$domain_name <- DOMAIN_NAMES[as.character(seu@meta.data[[COL_DOMAIN]])]
}
seu$domain_name <- factor(seu$domain_name, levels = DOMAIN_ORDER)

# =============================================================================
# Compute k-NN tumor purity per sample
# =============================================================================
all_coords <- Embeddings(seu, "spatial")
celltypes <- seu$Tier1_celltype
samples <- seu$sample
is_tumor <- as.integer(celltypes == TARGET_CELLTYPE)

cat(sprintf("Computing %d-NN purity...\n", K_NEIGHBORS))

purity_scores <- rep(NA_real_, ncol(seu))

for (s in levels(samples)) {
  idx <- which(samples == s)
  if (length(idx) < K_NEIGHBORS + 1) next

  coords_s <- all_coords[idx, , drop = FALSE]
  knn_result <- get.knn(coords_s, k = K_NEIGHBORS)
  is_tumor_s <- is_tumor[idx]

  neighbor_is_tumor <- matrix(is_tumor_s[knn_result$nn.index],
                              nrow = nrow(knn_result$nn.index),
                              ncol = K_NEIGHBORS)
  purity_scores[idx] <- rowMeans(neighbor_is_tumor)
  cat(sprintf("  %s: %d cells\n", s, length(idx)))
}

# =============================================================================
# Build tumor-cell data frame with domain info
# =============================================================================
tumor_idx <- which(celltypes == TARGET_CELLTYPE)
df <- data.frame(
  cell = colnames(seu)[tumor_idx],
  sample = samples[tumor_idx],
  condition = seu$condition[tumor_idx],
  domain = seu$domain_name[tumor_idx],
  purity = purity_scores[tumor_idx],
  stringsAsFactors = FALSE
)
df <- df[!is.na(df$purity) & !is.na(df$domain), ]
df$domain <- factor(df$domain, levels = DOMAIN_ORDER)

cat(sprintf("\nTumor cells with purity + domain: %d\n", nrow(df)))
cat("Per domain:\n")
print(table(df$domain, df$condition))

# =============================================================================
# Summary statistics: purity by domain and condition
# =============================================================================
domain_summary <- df %>%
  group_by(domain, condition) %>%
  summarise(
    n_cells = n(),
    mean_purity = mean(purity),
    median_purity = median(purity),
    sd_purity = sd(purity),
    .groups = "drop"
  )

cat("\nPurity by domain x condition:\n")
print(as.data.frame(domain_summary))

write.csv(domain_summary, file.path(OUTPUT_DIR, "purity_by_domain_condition.csv"),
          row.names = FALSE)

# =============================================================================
# Plot 1: Violin of purity by domain, colored by condition
# =============================================================================
p_domain_cond <- ggplot(df, aes(x = domain, y = purity, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width", position = position_dodge(0.8)) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8,
               position = position_dodge(0.8)) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "k-NN Tumor Purity",
       title = sprintf("Tumor Purity by Domain and Condition (k = %d)", K_NEIGHBORS)) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "violin_purity_domain_x_condition.pdf"),
       p_domain_cond, width = 8, height = 5)

# =============================================================================
# Plot 2: Violin of purity by domain (collapsed across conditions)
# =============================================================================
p_domain <- ggplot(df, aes(x = domain, y = purity, fill = domain)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.8) +
  scale_fill_manual(values = DOMAIN_COLORS) +
  labs(x = NULL, y = "k-NN Tumor Purity",
       title = sprintf("Tumor Purity by Spatial Domain (k = %d)", K_NEIGHBORS)) +
  theme_publication() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "violin_purity_by_domain.pdf"),
       p_domain, width = 6, height = 5)

# =============================================================================
# Plot 3: Per-sample dot plot for the two domains of interest
# =============================================================================
df_doi <- df[df$domain %in% DOMAINS_OF_INTEREST, ]
df_doi$domain <- droplevels(df_doi$domain)

sample_doi_summary <- df_doi %>%
  group_by(sample, condition, domain) %>%
  summarise(median_purity = median(purity), n = n(), .groups = "drop")

p_doi_samples <- ggplot(sample_doi_summary,
                        aes(x = domain, y = median_purity, color = condition)) +
  geom_jitter(width = 0.15, size = 3) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               aes(group = condition), position = position_dodge(0.4),
               color = "black", linewidth = 0.4) +
  scale_color_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "Median Tumor Purity (per sample)",
       title = "Tumor Core vs. Immune Engulfing") +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "dotplot_purity_TC_vs_IE.pdf"),
       p_doi_samples, width = 5, height = 5)

# =============================================================================
# Statistical model: Beta GLMM — domain * condition interaction
# =============================================================================
eps <- 1e-4
df_doi$purity_beta <- pmin(pmax(df_doi$purity, eps), 1 - eps)
df_doi$domain <- relevel(df_doi$domain, ref = "Tumor Core")

cat("\nFitting beta GLMM: purity ~ domain * condition + (1 | sample) ...\n")
fit_interaction <- glmmTMB(purity_beta ~ domain * condition + (1 | sample),
                           data = df_doi,
                           family = beta_family(link = "logit"))

cat("\n--- Beta GLMM: domain * condition interaction ---\n")
print(summary(fit_interaction))

# Simpler model: domain + condition (no interaction)
fit_additive <- glmmTMB(purity_beta ~ domain + condition + (1 | sample),
                        data = df_doi,
                        family = beta_family(link = "logit"))

# LRT for interaction
lrt <- anova(fit_additive, fit_interaction)
cat("\nLikelihood ratio test for interaction:\n")
print(lrt)

# Extract key coefficients
coefs <- summary(fit_interaction)$coefficients$cond
cat("\n--- Model coefficients (logit scale) ---\n")
print(coefs)

# Predicted marginal means
newdata <- expand.grid(
  domain = levels(df_doi$domain),
  condition = levels(df_doi$condition)
)
# Predict on logit scale (re = 0 for population-level)
newdata$pred_logit <- predict(fit_interaction, newdata = newdata,
                              type = "link", re.form = ~ 0)
newdata$pred_purity <- plogis(newdata$pred_logit)

cat("\nPredicted marginal means (population-level):\n")
print(newdata)

write.csv(newdata, file.path(OUTPUT_DIR, "predicted_means_domain_condition.csv"),
          row.names = FALSE)

# Save full model output
sink(file.path(OUTPUT_DIR, "beta_glmm_domain_condition.txt"))
cat("Beta GLMM: purity ~ domain * condition + (1 | sample)\n")
cat("Family: beta (logit link)\n")
cat(sprintf("Reference: domain = Tumor Core, condition = %s\n", REFERENCE_CONDITION))
cat(sprintf("N cells = %d, N samples = %d\n\n", nrow(df_doi), length(unique(df_doi$sample))))
print(summary(fit_interaction))
cat("\n\nLikelihood ratio test (interaction vs additive):\n")
print(lrt)
cat("\n\nPredicted marginal means:\n")
print(newdata)
sink()

# =============================================================================
# Plot 4: Predicted means bar plot with error bars and significance lines
# =============================================================================

# Compute SE on response scale via delta method
# SE(plogis(eta)) ≈ dlogis(eta) * SE(eta)
pred_obj <- predict(fit_interaction, newdata = newdata, type = "link",
                    re.form = ~ 0, se.fit = TRUE)
newdata$se_logit <- pred_obj$se.fit
newdata$ci_lo <- plogis(newdata$pred_logit - 1.96 * newdata$se_logit)
newdata$ci_hi <- plogis(newdata$pred_logit + 1.96 * newdata$se_logit)

# P-values for significance annotations
coef_table <- summary(fit_interaction)$coefficients$cond
p_domain <- coef_table["domainImmune Engulfing", "Pr(>|z|)"]
p_condition <- coef_table["conditionDKO", "Pr(>|z|)"]
p_interaction <- coef_table["domainImmune Engulfing:conditionDKO", "Pr(>|z|)"]

format_pval <- function(p) {
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  if (p < 0.1) return(paste0("p=", formatC(p, digits = 2, format = "f")))
  return("ns")
}

# Position significance lines above the upper CI
ci_max_tc <- max(newdata$ci_hi[newdata$domain == "Tumor Core"])
ci_max_ie <- max(newdata$ci_hi[newdata$domain == "Immune Engulfing"])

sig_y_tc <- ci_max_tc + 0.06
sig_y_ie <- ci_max_ie + 0.06
sig_y_domain <- max(sig_y_tc, sig_y_ie) + 0.10

p_pred <- ggplot(newdata, aes(x = domain, y = pred_purity, fill = condition)) +
  geom_col(position = position_dodge(0.7), width = 0.6, alpha = 0.85) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                position = position_dodge(0.7), width = 0.2, linewidth = 0.4) +
  # Significance: WT vs DKO within Tumor Core
  annotate("segment", x = 0.825, xend = 1.175, y = sig_y_tc, yend = sig_y_tc, linewidth = 0.4) +
  annotate("text", x = 1, y = sig_y_tc + 0.02, label = format_pval(p_condition), size = 3.5) +
  # Significance: WT vs DKO within Immune Engulfing
  annotate("segment", x = 1.825, xend = 2.175, y = sig_y_ie, yend = sig_y_ie, linewidth = 0.4) +
  annotate("text", x = 2, y = sig_y_ie + 0.02,
           label = format_pval(p_interaction), size = 3.5) +
  # Significance: domain effect (Tumor Core vs Immune Engulfing, bracket across)
  annotate("segment", x = 1, xend = 2, y = sig_y_domain, yend = sig_y_domain, linewidth = 0.4) +
  annotate("segment", x = 1, xend = 1, y = sig_y_domain - 0.02, yend = sig_y_domain, linewidth = 0.4) +
  annotate("segment", x = 2, xend = 2, y = sig_y_domain - 0.02, yend = sig_y_domain, linewidth = 0.4) +
  annotate("text", x = 1.5, y = sig_y_domain + 0.02, label = format_pval(p_domain), size = 3.5) +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(x = NULL, y = "Predicted Tumor Purity",
       title = "Model-Predicted Purity: Domain x Condition",
       subtitle = sprintf("Interaction LRT p = %.3g", lrt[["Pr(>Chisq)"]][2])) +
  coord_cartesian(ylim = c(0, sig_y_domain + 0.08)) +
  theme_publication()

ggsave(file.path(OUTPUT_DIR, "barplot_predicted_purity.pdf"),
       p_pred, width = 5, height = 5)

# =============================================================================
# Plot 5: Violin by domain, faceted by sample
# =============================================================================
p_by_sample <- ggplot(df_doi, aes(x = domain, y = purity, fill = domain)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.8) +
  facet_wrap(~ sample, nrow = 1) +
  scale_fill_manual(values = DOMAIN_COLORS) +
  labs(x = NULL, y = "k-NN Tumor Purity",
       title = sprintf("Tumor Purity by Domain per Sample (k = %d)", K_NEIGHBORS)) +
  theme_publication() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 9))

ggsave(file.path(OUTPUT_DIR, "violin_purity_domain_by_sample.pdf"),
       p_by_sample, width = 12, height = 5)

cat("\nDone. Outputs saved to:", OUTPUT_DIR, "\n")
