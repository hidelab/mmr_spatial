# Changelog

All notable changes to this project are documented here.
Newest entries appear at the top.

---

## 2026-07-27

### Fixed
- `code/R/02_03_merge_domains_08lambda.R`: Fixed readability of heatmap text by replacing dark navy-blue color palette (`"white", "steelblue", "navy"`) with a lighter palette (`"#FFFFFF", "#E8F4F8", "#A8D8E8", "#4A9FBC"`) in three heatmaps (merged domain cell-type composition, cell-type distribution across domains, and per-condition variants). The original palette made black text on dark-blue cells nearly invisible; the new palette maintains visual gradient while keeping text legible.

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