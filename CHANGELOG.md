# Changelog

All notable changes to this project are documented here.
Newest entries appear at the top.

---

## 2026-08-27

### Added
- `code/R/15_05_spatial_gradient_masks.R`: New script (adapted from `08_01_spatial_masks_08lambda.R`) that plots the 15_04 DKO-tumor-core-surge trajectory genes (Nt5e, Bcl2, Icos, Irf1, Gzmb, Prf1, Cd274, Fasl, Il7r, Cx3cl1, Vegfa, Nfatc1, Lat2, Grem1, Rnf213) in tissue space. One multi-panel figure per gene: masked cells (default `MASK_CELLTYPE = "CD8+ T"`) colored by expression, overlaid on the tumor mass drawn in grey (`SHOW_TUMOR_BACKGROUND`) as spatial context; samples ordered WT→DKO. A **shared per-gene color scale** (clipped at `EXPR_CLIP_Q`, default 0.99) across all panels makes WT-vs-DKO directly comparable, and high-expressers are drawn last so they aren't occluded. Confirms visually that these genes pick up specifically in CD8+ T cells inside the DKO tumor core (e.g. Gzmb dense/bright across DKO1/2/4, dim in WT). Outputs `spatial_<gene>.{pdf,png}` to `output/R/15_spatial_gradient_masks/`.
- `code/R/15_04_tumor_trajectory_working.R`: New **significant-trajectory heatmap** (`trajectory_heatmap_significant_<model>.pdf`) — a readable replacement for the 100 per-gene GAM scatter+spline PDFs. Genes are filtered to those passing **BH-FDR < `SIG_FDR`** (default 0.05) on `p_shape` **or** `p_level` (new `fdr_shape`/`fdr_level` columns added to `gam_condition_test_<model>.csv`). The kept genes are shown as a WT | DKO side-by-side ComplexHeatmap of the **GAM-fitted** population-level curves along the tumor-proximity axis (periphery→core), per-gene z-scored across *both* conditions on a shared scale so a *level* shift reads as one panel uniformly brighter and a *shape* change as the two panels diverging. Rows are cascade-ordered by fitted-peak position; a `class` strip labels each gene shape/level/both and capped `-log10 FDR` barplots (clipped at 10 so a single ~0 p-value like Bcl2 doesn't flatten the axis) annotate significance. Added `PLOT_INDIVIDUAL_GAMS` (default `FALSE`, skips the per-gene PDFs), `SIG_FDR`, and `HEATMAP_ROW_ORDER` (`"cluster"` — hierarchical clustering of the fitted curves, default; or `"shape"` — ranked by `p_shape`) parameters. Panels are 6.5 cm wide with a compact per-row height so the figure is wide-and-short. For CD8+ T, 32/100 genes pass (28 shape-only, 3 both — Cd3e/Nkg7/Tnfrsf1b, 1 level-only — Thy1); effector/checkpoint genes (Prf1, Gzmb, Cd274, Fasl, Icos, Irf1) show the DKO-specific surge toward the tumor core that WT lacks.

### Changed
- `code/R/15_04_tumor_trajectory_working.R`: Set `TOP_N_GENES <- 200` and added `limitsize = FALSE` to the two gene-count-scaled ggsave heatmaps (`heatmap_trajectory_combined.pdf`, `normalized_expression_heatmap.pdf`) so they render past ggplot's 50-inch height guard at large gene counts (200 genes → 50 in). Re-ran pseudobulk + ranef on 200 genes: pseudobulk 61 significant-trajectory genes, ranef 102 (FDR<0.1).
- `code/R/15_04_tumor_trajectory_working.R`: Promoted **`ranef`** (single-cell GAM with a per-sample random intercept, `s(sample, bs="re")`) to a full secondary spec run alongside `pseudobulk` on **all** `top_genes` (100). Added `FULL_GENE_MODELS` (default `c("pseudobulk","ranef")`) — models fit on the full gene set; any other GAM model (e.g. `fs`) is still restricted to the top `N_GENES_OTHER_MODELS` pseudobulk hits. Changed `PLOT_INDIVIDUAL_GAMS` from a logical to a **vector of model names** for which to emit the per-gene `gam_<gene>.pdf` spline curves (default `c("ranef")`; `character(0)` disables) — so ranef now gets both the 100 per-gene spline PDFs and its own significant-trajectory heatmap. Dropped `fs` from the default `GAM_MODELS`. For CD8+ T, ranef flags 73 genes at FDR<0.1 (62 shape-only, 9 both, 2 level-only) vs pseudobulk's 48 — ranef is more anticonservative (random *intercept* only absorbs per-sample level, not per-sample trajectory shape), consistent with the pseudoreplication trade-off; pseudobulk remains the primary/defensible test.
- `code/R/15_04_tumor_trajectory_working.R`: Added an optional **QC filter** (`APPLY_QC`, matching the `08_01`/`15_05` convention — `total_counts >= 50` and `n_genes_by_counts >= 10`), computed over all cells but applied only at per-cell-type selection so it never alters the tumor-proximity score field. After review the filter is left **OFF by default** (`APPLY_QC <- FALSE`): the analysis uses all CD8+ T cells (12612). Also raised the significant-trajectory-heatmap cutoff `SIG_FDR` from `0.05` to `0.1` (a gene enters if `fdr_shape < 0.1` OR `fdr_level < 0.1`) to admit borderline trajectory-shape genes for visualization; 48 genes now pass (41 shape-only, 7 both, 0 level-only).
- `code/R/15_05_spatial_gradient_masks.R`: Removed the QC subset (reviewed) — all CD8+ T cells (12612) are now plotted, consistent with 15_04.
- `code/R/15_04_tumor_trajectory_working.R`: Switched the default tumor score from KDE to **graph heat-diffusion** (`SCORE_METHOD <- "diffusion"`, `DIFFUSION_VARIANT <- "heat"`) and added a `DIFFUSION_HEAT_STEPS` parameter (default 20) that is now forwarded to `score_diffusion()` via `compute_tumor_score(..., heat_steps=)`. The heat-diffusion field is more sharply concentrated near 0 than KDE (long right tail), emphasizing cells inside/adjacent to tumor; quantile binning adapts so the pipeline is unaffected. CD8+ T top trajectory-shape hits remain biologically coherent and overlap the KDE result (Cd274, Prf1, Il2rb, Lat2, Cxcl10, Fasl, Pdcd1, Cd3e — robust across scoring methods), and the diffusion axis additionally promotes the chemokine axis (Cxcr6, Cxcl10, Cx3cl1, Cxcr3) including the project-of-interest Cxcl16/Cxcr6.
- `code/R/15_04_tumor_trajectory_working.R`: Replaced the condition test in `fit_trajectory_gam()`. The old approach compared two penalized GAMs via `anova(null_ml, full_ml, test="Chisq")` (ML refit), which returned `NA` for ~15% of genes — diagnosis showed this is the well-known unreliability of deviance-based LRTs on penalized smooths (the penalized "full" model can have *higher* deviance than the null, yielding a negative, untestable chi-square), and it conflated overall-level with trajectory-shape differences. Switched to the mgcv-canonical **ordered-factor difference smooth** (Wood 2017): `expr ~ s(score) + s(score, by = cond_ord) + cond_ord + <rand eff>` with `condition` as an ordered factor (WT = reference). This is a single REML fit (also ~1.4× faster than the prior three-fit LRT) and yields two interpretable, never-`NA` p-values per gene: **`p_shape`** (the `s(score):cond_ord` deviation smooth — does the WT→DKO *trajectory shape* differ, the trajectory question) and **`p_level`** (the `cond_ord` parametric offset — overall up/down shift, already covered by bulk DE). The per-gene CSV (`gam_condition_test_<model>.csv`) now has `gene,p_shape,p_level`, sorted by `p_shape`; GAM curve subtitles show both. Applied to all three GAM specs (`pseudobulk`/`fs`/`ranef`). For CD8+ T the top trajectory-shape hits are Grem1, Cd274, Irf1, Cd3e, Prf1, Il2rb, Cxcl10, Ctla4, Fasl, Pdcd1 — checkpoint/cytotoxicity/IFN programs whose spatial gradient reshapes in DKO.
- `code/R/15_04_tumor_trajectory_working.R`: Lowered the default `KDE_BW_SCALE` from `2.0` to `1.0` — at `2.0` the tumor-density field was over-smoothed, producing unrealistic per-sample score distributions for the DKO samples (a broad plateau plus an artifactual sharp spike near score≈0.45 in DKO3/DKO4). At `1.0` (the raw `bandwidth.nrd` rule-of-thumb) every sample shows the expected right-skewed decay (most cells at the periphery, fewer deep inside tumor). Also set `TOP_N_GENES` to `100` and confirmed `GAM_MODELS <- c("pseudobulk")` for the current CD8+ T analysis.

### Added
- `code/R/00_utils.R`: New shared utilities file holding all tumor-proximity scoring options, refactored out of the 15_xx/16_01 scripts: `score_knn_purity()` (k-NN tumor purity), `score_kde()` (`kde2d` + bilinear interpolation), `score_diffusion()` (graph-Laplacian diffusion with `"implicit"`/`"heat"`/`"ppr"` variants, plus `cg_solve()`/`.diffuse_field()` helpers), and a `compute_tumor_score()` dispatcher that loops per sample, calls the chosen scorer via `method=` + `...`, and min-max normalizes to [0,1]. Lets downstream scripts pick a scoring method with a one-line switch.
- `code/R/15_04_tumor_trajectory_working.R`: New "working" tumor-trajectory script that addresses pseudoreplication in the earlier single-cell GAMs (3 WT + 4 DKO samples treated as thousands of independent cells). Scores cells via `compute_tumor_score()` (default KDE). Provides a switchable menu of **sample-aware GAM specifications** — `"fs"` (factor-smooth random effect `s(score, sample, bs="fs")`), `"ranef"` (random intercept `s(sample, bs="re")`), and `"pseudobulk"` (sample×bin weighted GAM) — fit and saved separately per model, with the condition-effect p-value from an ML-refit likelihood-ratio test. Adds binning-quality diagnostics (score histogram with decile breaks, cells-per-sample×bin tile, usable-bins-per-sample bar, counts CSV) and a normalized-expression-along-gradient "wave" plot (each gene min-max scaled to [0,1]). Gene selection is switchable (`binned_variance` default, `spearman`, and an `svg` hook for a future top-200 spatially-variable-gene set). Outputs to `output/R/15_tumor_trajectory_working`.

---

## 2026-08-26

### Added
- `code/R/16_01_tumor_contours.R`: Added a `DIFFUSION_METHOD` parameter (`"implicit"`, `"heat"`, `"ppr"`) selecting how the binary tumor indicator is diffused over the spatial k-NN graph. `"heat"` (Option 1) applies `HEAT_STEPS` backward-Euler steps for a long total diffusion time (genuine long-range spread); `"ppr"` (Option 2) uses random-walk personalized PageRank `(1-β)(I-βP)^{-1}s_raw` with tunable `PPR_BETA`. Factored the diffusion math into `compute_diffusion()`/`cg_solve()` helpers and replaced the inline Laplacian solve in the per-sample loop. The prior single-step behavior is preserved as `"implicit"`. Motivation: raising `ALPHA_DIFFUSION` alone left scores strongly skewed toward zero (the implicit solve is only a mild low-pass filter).
- `code/R/16_01_tumor_contours.R`: Added a per-sample histogram of the diffused tumor score (`tumor_score_histogram_per_sample.{pdf,png}`) emitted right after score computation for diagnostics, and an isolated contour map of the tumor-score field (`spatial_tumor_score_contour.{pdf,png}`) via `akima::interp` + `geom_contour_filled`.
- `code/R/15_02_tumor_trajectory_KDE.R`: Added a `KDE_BW_SCALE` parameter (default 2.0) that multiplies the `kde2d` bandwidth, letting the tumor-density field be spread wider (>1) or sharpened (<1). Applied to both auto-estimated (`bandwidth.nrd`) and user-supplied bandwidths; the effective scaled bandwidth is now reported in the per-sample log.

---

## 2026-07-27

### Changed
- `code/R/02_03_merge_domains_08lambda.R`: Replaced `pheatmap` with `ComplexHeatmap` for all five heatmaps to enable proper 45° column-name rotation (`column_names_rot = 45`), improving readability. Unified color function using `colorRamp2(c(0, 0.5, 1), c("#FFFFFF", "#A8D8E8", "#4A9FBC"))` for consistency and better contrast. Added `cell_fun` with `grid.text()` to display cell values on numeric heatmaps. Replaced `pheatmap` library import with `ComplexHeatmap`, `circlize`, and `grid`.

### Added
- `code/R/11_01_spatial_dirichlet.R`: Added a third Dirichlet regression analysis testing whether **spatial domain** composition (Banksy merged domains via `COL_DOMAIN`, mapped to `DOMAIN_NAMES`/`DOMAIN_ORDER`) differs between DKO and WT. Mirrors the existing immune and epithelial/fibroblast blocks: builds a sample×domain count table over all QC-passing cells, fits `DirichReg(AL ~ condition)` with an LRT vs the null, and emits results CSVs plus stacked-bar, collapsed-bar, boxplot-with-significance, and dot-plot figures (`dirichlet_DKOvsWT_domains*`, `domain_*`). Updated the header docstring and end-of-run summary accordingly.

### Changed
- `code/R/11_01_spatial_dirichlet.R`: Refactored to load the Seurat object **once** (via a single `readRDS()` call) rather than three independent times. All three analyses (immune, epithelial/fibroblast, spatial domain) now filter copies of `seu_full` (after universal QC), eliminating redundant I/O on the large ~328k-cell object. Reorganized sections and updated comments to reflect the new structure.

### Changed
- `code/R/12_pathway_activity.R`: Replaced the naive mouse→human gene mapping (which simply uppercased mouse symbols) with a proper orthology-based conversion using the `babelgene` package (offline HCOP consensus of 12 databases). This fixes both the forward mapping (mouse symbol → human ENSEMBL for PanomiR input) and the reverse mapping (human ENSEMBL → mouse symbol for heatmaps). Uppercasing silently missed real orthologs whose human symbol differs (e.g. `Trp53`→`TP53`, `H2-K1`→`HLA-A`) and risked spurious collisions with unrelated human genes. Forward mapping now resolves to 1:1 (strongest support) to avoid duplicate ENSEMBL row names; dropped the now-unused `org.Hs.eg.db`/`AnnotationDbi` imports.

---

## 2026-07-24

### Added
- `code/R/05_02_DEG_DESeq_cellxdomains.R`: New script containing all plotting code (volcano, MA, heatmap, summary bar) extracted from `05_DEG_DESeq_cellxdomains.R`. Loops over cell types and reads per-celltype CSV results to produce visualizations independently of the DEG computation.
- `code/R/04_02_DEG_DESeq_merged_domains.R`: New script containing all plotting code (volcano, MA, heatmap, summary bar) extracted from `04_DEG_DESeq_merged_domains.R`. Reads saved CSV results and produces all visualizations independently of the DEG computation.

### Changed
- `code/R/05_DEG_DESeq_cellxdomains.R`: Removed plotting sections (volcano, MA, heatmap, summary bar). Script now only performs DEG calculation and saves result tables. Plotting is handled by `05_02_DEG_DESeq_cellxdomains.R`.
- `code/R/04_DEG_DESeq_merged_domains.R`: Removed plotting sections (volcano, MA, heatmap, summary bar, pseudobulk QC plots). Script now only performs DEG calculation and saves result tables. Plotting is handled by `04_02_DEG_DESeq_merged_domains.R`.

---

## 2026-07-23

### Changed
- `code/R/07_02_plot_subjects_merged.R`: Reworked per-subject domain panel layout for publication figures.
  - The "All Domains" merged panel now appears as the **leftmost** plot in each subject figure (was previously rightmost).
  - Per-domain highlight panels are now ordered by `DOMAIN_ORDER` from `00_parameters.R` (Tumor Core, Immune Engulfing, Stroma, Mammary Glands) instead of numeric domain ID.
  - Per-domain panels keep non-domain cells in gray and highlight domain cells in their original domain color.
  - Refactored the per-sample plotting into a reusable `build_sample_panels()` helper (also robust to `subset()` inside a function via `cells =`).

### Fixed
- `code/R/07_02_plot_subjects_merged.R`: Domain-highlight color bug — the `"Domain"` color came from a *named* element (`DOMAIN_COLORS_NUM[dom]`), which mangled the color-map key to `"Domain.<id>"` so `scale_color_manual` failed to match and domain cells rendered as grey `na.value`. Now `unname()`-ed so domain cells show their true domain color.
- `code/R/07_02_plot_subjects_merged.R`: Point-size mismatch between the All Domains panel and per-domain panels. In raster mode Seurat's on-screen dot size depends on each panel's rendered geometry (not just `pt.size`), so the legend on the All Domains panel plus per-panel coord limits made its dots render smaller. Now every panel is legend-free, shares one `PT_SIZE` constant, and pins identical `coord_fixed(xlim, ylim)` from the sample's spatial range — guaranteeing identical dot size across panels.

### Added
- `code/R/07_02_plot_subjects_merged.R`: New stacked all-samples figure `spatial_subjects_all_stacked.pdf` — one row of panels per sample, column titles on the top row, sample name as the left row label.

---

## 2026-07-21

### Changed
- `code/R/04_DEG_DESeq_merged_domains.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Removed local QC filters, gene-level filters, thresholds, condition labels, `DOMAIN_COL`, `DOMAIN_NAMES`.
  - Fixed `DOMAIN_NAMES` inconsistency (underscores vs spaces) by using centralized version.
  - Volcano top genes uses `VOLCANO_TOP_N`.
- `code/R/05_DEG_DESeq_cellxdomains.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - `CELLTYPES_OF_INTEREST` → `IMMUNE_CELLTYPES_DE`.
  - `COMPARISONS` built dynamically from `DOMAIN_COMPARISONS`.
  - Volcano top genes uses `VOLCANO_TOP_N`.
- `code/R/06_DEG_DESeq_cellxdomain_DKOvsWT.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Same pattern as 05: `IMMUNE_CELLTYPES_DE`, `COL_DOMAIN`, `VOLCANO_TOP_N`.
- `code/R/07_01_plot_domains_merged_08lambda.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - `DOMAIN_COLORS_NUM` derived from centralized `DOMAIN_COLORS` + `DOMAIN_NAMES` mapping.
- `code/R/07_02_plot_subjects_merged.R`: Same refactor as 07_01.
- `code/R/07_03_plot_cell_types.R`: Replaced `PROJECT_DIR` and RDS path with `SEURAT_RDS`.
- `code/R/08_01_spatial_masks_08lambda.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Mask parameters from `MASK_DOMAIN_DEFAULT`, `MASK_CELLTYPE_DEFAULT`.
  - Gene list uses `GENES_SPATIAL_VIS`.
- `code/R/08_02_violin_masks.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Mask parameters from `MASK_DOMAIN_VIOLIN`, `MASK_CELLTYPE_VIOLIN`.
  - `PANEL_NCOL` for multi-panel layout.
- `code/R/08_03_spatial_features.R`: Replaced `PROJECT_DIR` and gene list with `GENES_SPATIAL_VIS`.
- `code/R/08_04_violin_x_cells.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - `MASK_DOMAINS_VIOLIN_X`, `IMMUNE_CELLTYPES_DE`, `PANEL_NCOL`.
  - Script-specific DEG genes appended to `GENES_SPATIAL_VIS`.
- `code/R/09_proportion_analysis.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - `DOMAIN_COLORS_NAMED` → uses centralized `DOMAIN_COLORS` directly.
  - `IMMUNE_CELLTYPES` from parameters; domain factor levels from `DOMAIN_ORDER`.
- `code/R/10_minor_diagnostics.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Removed local `DC_TYPES`, `DOMAIN_NAMES`, `DOMAIN_COLORS_NAMED` definitions.
  - Domain factor levels from `DOMAIN_ORDER`.
- `code/R/11_spatial_dirichlet.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - `IMMUNE_CELLTYPES`, `EPITHELIAL_STROMAL_CELLTYPES`.
  - Significance threshold uses `FDR_THRESHOLD`.
- `code/R/12_pathway_activity.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - MSigDB paths: `MSIGDB_PATHWAYS`, `MSIGDB_NAMES`.
  - `MIN_PATH_SIZE`, `TOP_PATHWAYS_PLOT`, `PATHWAY_NAME_TRUNC`, `FDR_THRESHOLD`.
  - Contrast condition built dynamically from `TEST_CONDITION`/`REFERENCE_CONDITION`.
- `code/R/13_DE_pathways_immuneXdomains.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Same pathway parameters as 12; `IMMUNE_CELLTYPES_DE` for cell types.
  - Domain subset kept in-script (excludes Mammary Glands).

- `code/R/02_03_merge_domains_08lambda.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Removed local `DOMAIN_COLORS`, `DOMAIN_NAMES`, `DOMAIN_ORDER`, `MERGE_K` definitions (now from parameters).
  - Domain column uses `DOMAIN_SOURCE_COL`.
  - Hierarchical clustering uses `MERGE_DIST_METHOD` and `MERGE_LINK_METHOD`.
  - Paths: `BANKSY_RDS` input, derived `BANKSY_MERGED_RDS` output.
  - Condition recode: `TEST_CONDITION`.
  - Fixed domain color mismatch in spatial/UMAP plots: created `DOMAIN_COLORS_NUM` mapping numeric factor levels to colors.
- `code/R/02_02_reclustering_resolutions_08lambda.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Resolution sweep uses `BANKSY_RESOLUTIONS_SWEEP`.
  - Primary resolution reference uses `BANKSY_RESOLUTION` and `BANKSY_NPC_USE`.
  - Paths: `BANKSY_RDS` for input/output.
  - Condition recode: `TEST_CONDITION`.
- `code/R/02_01_reclustering_seurat_banksy_08lambda.R`: Replaced hardcoded parameters with `00_parameters.R` constants.
  - Banksy settings: `BANKSY_LAMBDA`, `BANKSY_K_GEOM`, `BANKSY_NPCS`, `BANKSY_NPC_USE`, `BANKSY_RESOLUTION`.
  - Harmony: `HARMONY_GROUP_VAR`.
  - Paths: `SEURAT_RDS` for input, derived `BANKSY_RDS` for output.
  - Condition recode: `TEST_CONDITION`.
  - Primary cluster column name derived from `BANKSY_NPC_USE` and `BANKSY_RESOLUTION`.
- `code/R/01_data_import.R`: Replaced all hardcoded parameters with references to `00_parameters.R`.
  - `PROJECT_DIR` now sourced from `00_parameters.R` instead of being redefined.
  - Input path uses `H5AD_INPUT` constant.
  - Condition recode uses `TEST_CONDITION` instead of literal `"DKO"`.
  - Sample-to-condition mapping uses `WT_SAMPLES`, `REFERENCE_CONDITION`, `TEST_CONDITION` instead of regex.
  - Seurat save path uses `SEURAT_RDS` constant.