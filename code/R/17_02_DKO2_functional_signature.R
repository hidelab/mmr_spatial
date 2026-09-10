# 17_02_DKO2_functional_signature.R
# Stronger, functional evidence on whether DKO2 is a failed/leaky Mlh1/Msh2 KO.
#
# Rationale: the KO transcripts sit near the detection floor on this 680-gene
# panel (see 17_01), so the transcript-level argument is inherently weak. A
# *successful* MMR knockout, however, causes microsatellite instability ->
# hypermutation -> neoantigens -> a robustly-expressed, MMR-deficiency-specific
# IMMUNOGENIC program that we CAN measure well:
#   - tumor-cell-intrinsic IFN response (Irf1/7, Gbp3/4/5, Ifit1/3, Usp18, Cxcl9/10/11)
#   - antigen-presentation machinery (Nlrc5 = MHC-I transactivator, Ciita = MHC-II, Psmb9)
#   - cGAS-STING sensing of genomic instability
#   - checkpoint induction (Cd274/PD-L1 — IFN-driven, hallmark of MSI tumors)
#   - cytotoxic lymphocyte infiltration (CD8+ T / NK composition)
# If DKO2's tumor is immunologically WT-like (cold, low IFN/antigen presentation)
# while the other DKOs are hot, that is strong DIRECT evidence its KO is inactive.
#
# Also re-tests the KO transcripts with proper DEPTH CONTROL (DKO2 has the highest
# counts/cell, which inflated its raw KO-gene counts in 17_01).

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

TUMOR_CELLTYPE <- "Tumor"

# --- Functional modules (in-panel genes only) --------------------------------
MODULES <- list(
  IFN_response      = c("Irf1","Irf7","Gbp3","Gbp4","Gbp5","Ifit1","Ifit3",
                        "Usp18","Cxcl9","Cxcl10","Cxcl11","Ifngr1"),
  Antigen_presentation = c("Nlrc5","Ciita","Psmb9"),
  cGAS_STING        = c("Cgas","Sting1","Ifih1","Tbk1","Irf3"),
  Checkpoint        = c("Cd274","Pdcd1","Ctla4","Lag3","Tigit"),
  Cytotoxicity      = c("Gzmb","Gzmk","Prf1","Nkg7","Klrk1","Ifng","Gzma")
)
KO_GENES <- c("Mlh1","Msh2","Msh6")

# =============================================================================
# Load + normalize
# =============================================================================
seu <- readRDS(SEURAT_RDS)
seu$condition <- gsub("^KO$", TEST_CONDITION, seu$condition)
sample_levels <- c(WT_SAMPLES, DKO_SAMPLES)
sample_levels <- sample_levels[sample_levels %in% unique(as.character(seu$sample))]
seu$sample <- factor(seu$sample, levels = sample_levels)
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)

cond_of <- function(s) ifelse(s %in% WT_SAMPLES, REFERENCE_CONDITION, TEST_CONDITION)

# =============================================================================
# (A) Functional signatures in TUMOR cells, per sample
# =============================================================================
tum <- subset(seu, subset = Tier1_celltype == TUMOR_CELLTYPE)
cat(sprintf("Tumor cells: %d\n", ncol(tum)))

mods_present <- lapply(MODULES, function(g) intersect(g, rownames(tum)))
# nbin lowered (small 680-gene panel: default 24 bins leaves <30 genes/bin)
tum <- AddModuleScore(tum, features = mods_present, name = "MOD_",
                      ctrl = 20, nbin = 12, seed = 42)
score_cols <- paste0("MOD_", seq_along(mods_present))
names(score_cols) <- names(mods_present)

md <- tum@meta.data
md$sample <- factor(md$sample, levels = sample_levels)

# per-sample mean of each module score in tumor cells
sig_tbl <- do.call(rbind, lapply(names(score_cols), function(m) {
  v <- md[[score_cols[m]]]
  data.frame(module = m,
             sample = levels(md$sample),
             mean_score = tapply(v, md$sample, mean),
             n_tumor = as.integer(table(md$sample)),
             row.names = NULL)
}))
sig_tbl$condition <- factor(cond_of(as.character(sig_tbl$sample)),
                            levels = c(REFERENCE_CONDITION, TEST_CONDITION))
sig_tbl$sample <- factor(sig_tbl$sample, levels = sample_levels)
write.csv(sig_tbl, file.path(OUTPUT_DIR, "tumor_functional_signatures_by_sample.csv"),
          row.names = FALSE)

cat("\n--- Tumor-cell functional signature scores (per-sample mean) ---\n")
print(reshape(sig_tbl[, c("module","sample","mean_score")],
              idvar = "module", timevar = "sample", direction = "wide"),
      digits = 2, row.names = FALSE)

# per-sample mean bars, faceted by module; DKO2 outlined
sig_tbl$is_dko2 <- sig_tbl$sample == "DKO2"
p_sig <- ggplot(sig_tbl, aes(sample, mean_score, fill = condition)) +
  geom_col(aes(color = is_dko2), linewidth = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_color_manual(values = c(`FALSE` = NA, `TRUE` = "black"), guide = "none") +
  labs(title = "Tumor-cell functional signatures (per-sample mean module score)",
       subtitle = "MMR-deficient tumors should be IFN/antigen-presentation HIGH; DKO2 outlined",
       y = "mean AddModuleScore", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "tumor_functional_signatures.pdf"),
       p_sig, width = 3.4 * length(mods_present), height = 4.8)

# per-cell violins for the two most MMR-specific modules
viol_df <- do.call(rbind, lapply(c("IFN_response","Antigen_presentation"), function(m) {
  data.frame(module = m, sample = md$sample, condition = cond_of(as.character(md$sample)),
             score = md[[score_cols[m]]])
}))
viol_df$condition <- factor(viol_df$condition, levels = c(REFERENCE_CONDITION, TEST_CONDITION))
viol_df$sample <- factor(viol_df$sample, levels = sample_levels)
p_viol <- ggplot(viol_df, aes(sample, score, fill = condition)) +
  geom_violin(scale = "width", linewidth = 0.2) +
  stat_summary(fun = mean, geom = "point", size = 1.5, color = "black") +
  facet_wrap(~ module, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = CONDITION_COLORS) +
  labs(title = "Per-tumor-cell signature distribution", y = "module score", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "tumor_signature_violins.pdf"), p_viol, width = 9, height = 7)

# =============================================================================
# (B) Cd274 (PD-L1) — single most diagnostic IFN-driven MSI marker, in tumor
# =============================================================================
if ("Cd274" %in% rownames(tum)) {
  v <- as.numeric(GetAssayData(tum, layer = "data")["Cd274", ])
  pdl1 <- data.frame(sample = levels(md$sample),
                     mean_expr = tapply(v, md$sample, mean),
                     pct_pos   = tapply(v, md$sample, function(x) 100 * mean(x > 0)),
                     row.names = NULL)
  pdl1$condition <- factor(cond_of(as.character(pdl1$sample)),
                           levels = c(REFERENCE_CONDITION, TEST_CONDITION))
  pdl1$sample <- factor(pdl1$sample, levels = sample_levels)
  write.csv(pdl1, file.path(OUTPUT_DIR, "tumor_Cd274_by_sample.csv"), row.names = FALSE)
  cat("\n--- Cd274/PD-L1 in tumor cells ---\n"); print(pdl1, digits = 3, row.names = FALSE)
}

# =============================================================================
# (C) Cytotoxic immune infiltration, per sample (composition of all cells)
# =============================================================================
imm_types <- c("CD8+ T", "NK", "Treg")
comp <- as.data.frame.matrix(table(seu$sample, seu$Tier1_celltype))
comp_frac <- sweep(comp, 1, rowSums(comp), "/") * 100
infil <- data.frame(sample = rownames(comp_frac),
                    CD8_T = comp_frac[["CD8+ T"]],
                    NK    = if ("NK" %in% colnames(comp_frac)) comp_frac[["NK"]] else NA,
                    Treg  = if ("Treg" %in% colnames(comp_frac)) comp_frac[["Treg"]] else NA)
infil$cytotoxic <- rowSums(infil[, c("CD8_T","NK")], na.rm = TRUE)
infil$condition <- factor(cond_of(infil$sample), levels = c(REFERENCE_CONDITION, TEST_CONDITION))
infil$sample <- factor(infil$sample, levels = sample_levels)
write.csv(infil, file.path(OUTPUT_DIR, "immune_infiltration_by_sample.csv"), row.names = FALSE)
cat("\n--- Immune infiltration (% of all cells) ---\n"); print(infil, digits = 3, row.names = FALSE)

p_inf <- ggplot(infil, aes(sample, cytotoxic, fill = condition)) +
  geom_col(aes(color = sample == "DKO2"), linewidth = 0.7) +
  scale_fill_manual(values = CONDITION_COLORS) +
  scale_color_manual(values = c(`FALSE` = NA, `TRUE` = "black"), guide = "none") +
  labs(title = "Cytotoxic lymphocyte infiltration (CD8+ T + NK)",
       subtitle = "MSI/MMR-deficient tumors are cytotoxic-hot; DKO2 outlined",
       y = "% of all cells", x = NULL) +
  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUTPUT_DIR, "immune_infiltration.pdf"), p_inf, width = 6, height = 4.5)

# =============================================================================
# (D) Depth-controlled KO-transcript test in tumor cells
#     detected(gene) ~ sample + log(total_counts)  [logistic]  -> DKO2 vs other-DKO
# =============================================================================
co  <- GetAssayData(seu, assay = "RNA", layer = "counts")
tci <- which(seu$Tier1_celltype == TUMOR_CELLTYPE)
tc_smp   <- droplevels(seu$sample[tci])
tc_depth <- log1p(seu$total_counts[tci])
other_dko <- setdiff(DKO_SAMPLES, "DKO2")

ko_test <- do.call(rbind, lapply(KO_GENES, function(g) {
  det <- as.integer(co[g, tci] > 0)
  grp <- ifelse(as.character(tc_smp) == "DKO2", "DKO2",
         ifelse(as.character(tc_smp) %in% other_dko, "otherDKO", "WT"))
  grp <- factor(grp, levels = c("otherDKO", "DKO2", "WT"))   # otherDKO = reference
  fit <- tryCatch(glm(det ~ grp + tc_depth, family = binomial()), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  s <- summary(fit)$coefficients
  data.frame(gene = g,
             logOR_DKO2_vs_otherDKO = s["grpDKO2", "Estimate"],
             p_DKO2_vs_otherDKO     = s["grpDKO2", "Pr(>|z|)"],
             logOR_WT_vs_otherDKO   = s["grpWT",   "Estimate"],
             p_WT_vs_otherDKO       = s["grpWT",   "Pr(>|z|)"],
             row.names = NULL)
}))
write.csv(ko_test, file.path(OUTPUT_DIR, "KO_transcript_depthcontrolled_test.csv"),
          row.names = FALSE)
cat("\n--- Depth-controlled KO-transcript detection (logistic; ref = other DKOs) ---\n")
cat("    positive logOR = MORE detected than other DKOs (WT-like retention)\n")
print(ko_test, digits = 3, row.names = FALSE)

cat("\nDone. Outputs ->", OUTPUT_DIR, "\n")
