# ==============================================================================
# 00_parameters.R — Centralized user-defined parameters
# ==============================================================================
# Source this file at the top of every R analysis script:
#   source(file.path(PROJECT_DIR, "code", "R", "00_parameters.R"))
#
# All user-configurable values live here. Individual scripts should NOT
# redefine these unless there is a script-specific override clearly noted.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. PROJECT PATHS
# ------------------------------------------------------------------------------
PROJECT_DIR <- "." # assuming that the scripts are executed from the project root.
# Used in: ALL scripts (01–13)

# Key input files
H5AD_INPUT <- file.path(PROJECT_DIR, "input", "DKO_results", "DKO_Tier1.h5ad")
# Used in: 01_data_import.R (loaded)

SEURAT_RDS <- file.path(PROJECT_DIR, "output", "R", "DKO_Tier1_seurat.rds")
# Generated in: 01_data_import.R (saved)
# Used in: 02_01 (loaded), 03, 04, 05, 06 (loaded)

# Pathway databases
MSIGDB_PATHWAYS <- file.path(PROJECT_DIR, "input", "MSigDBPathGeneTab2024.RDS")
MSIGDB_NAMES    <- file.path(PROJECT_DIR, "input", "MSigDB_display_names_v2024.csv")
# Used in: 12_pathway_activity.R, 13_DE_pathways_immuneXdomains.R

# ------------------------------------------------------------------------------
# 2. EXPERIMENTAL DESIGN
# ------------------------------------------------------------------------------
REFERENCE_CONDITION <- "WT"
TEST_CONDITION      <- "DKO"
# Used in: 03, 04, 06 (DESeq2 reference level), 11, 12, 13

WT_SAMPLES  <- c("WT1", "WT2", "WT3")
DKO_SAMPLES <- c("DKO1", "DKO2", "DKO3", "DKO4")
# Used in: throughout all scripts for filtering/grouping

# Metadata column names
COL_SAMPLE    <- "sample"
COL_CONDITION <- "condition"
COL_CELLTYPE  <- "Tier1_celltype"
COL_DOMAIN    <- "banksy_domain_merged"
# Used in: ALL scripts for subsetting metadata

# ------------------------------------------------------------------------------
# 3. CELL-LEVEL QC FILTERS
# ------------------------------------------------------------------------------
MIN_TOTAL_COUNTS <- 50
# Used in: 03, 04, 05, 06, 08_01, 08_02, 08_04, 09, 11, 13

MIN_GENES_DETECTED <- 10
# Used in: 03, 04, 05, 06, 08_01, 08_02, 08_04, 09, 11, 13

MIN_CELLS_PER_SAMPLE <- 10
# Used in: 03, 04, 05, 06, 08_04, 13

# ------------------------------------------------------------------------------
# 4. GENE-LEVEL FILTERS (for DE analysis)
# ------------------------------------------------------------------------------
MIN_GENE_EXPR_FRAC <- 0.05       # Gene expressed in >= 5% of cells
# Used in: 03, 04, 05, 06

MIN_COUNTS_PER_SAMPLE <- 5       # Minimum counts in a sample
# Used in: 03, 04, 05, 06

MIN_SAMPLES_EXPRESSED <- 3       # Gene detected in >= 3 samples
# Used in: 03, 04, 05, 06

# ------------------------------------------------------------------------------
# 5. STATISTICAL THRESHOLDS
# ------------------------------------------------------------------------------
FDR_THRESHOLD <- 0.1
# Used in: 03, 04, 05, 06 (DE significance), 11 (Dirichlet)

LFC_THRESHOLD <- 0.25
# Used in: 03, 04, 05, 06 (volcano plot labeling threshold)

# Significance annotation (for Dirichlet / pathway plots)
SIG_LEVELS  <- c(0.1, 0.05, 0.01, 0.001)
SIG_SYMBOLS <- c("^", "*", "**", "***")
# Used in: 11_spatial_dirichlet.R

# ------------------------------------------------------------------------------
# 6. BANKSY / SPATIAL CLUSTERING
# ------------------------------------------------------------------------------
BANKSY_LAMBDA   <- 0.8    # Spatial context weight
BANKSY_K_GEOM   <- 30     # Neighborhood size
BANKSY_NPCS     <- 50     # PCA components to compute
BANKSY_NPC_USE  <- 30     # PCs for nearest neighbors / UMAP
BANKSY_RESOLUTION <- 0.8  # Seurat clustering resolution
# Used in: 02_01_reclustering_seurat_banksy_08lambda.R

# Additional resolutions for sweep
BANKSY_RESOLUTIONS_SWEEP <- c(0.05, 0.3, 0.5)
# Used in: 02_02_reclustering_resolutions_08lambda.R

# Harmony batch correction variable
HARMONY_GROUP_VAR <- "sample"
# Used in: 02_01_reclustering_seurat_banksy_08lambda.R

# ------------------------------------------------------------------------------
# 7. DOMAIN MERGING
# ------------------------------------------------------------------------------
DOMAIN_SOURCE_COL <- "banksy_res_0.8"
# Used in: 02_03, 04, 05, 06 (identifies source clustering column)

MERGE_K <- 4       # Number of merged domains
# Used in: 02_03_merge_domains_08lambda.R

# Domain naming and colors
DOMAIN_NAMES <- c(
  "1" = "Mammary Glands",
  "2" = "Tumor Core",
  "3" = "Immune Engulfing",
  "4" = "Stroma"
)
# Generated in: 02_03 (cluster-to-name mapping)
# Used in: 04, 05, 06, 07_01, 07_02, 08_04, 09, 10

DOMAIN_COLORS <- c(
  "Mammary Glands"   = "#F9A825",
  "Tumor Core"       = "#B71C1C",
  "Immune Engulfing" = "#4A148C",
  "Stroma"           = "#4b5c09" #"#BDBDBD"
)
# Used in: 02_03, 07_01, 07_02, 08_04, 09, 10

# Display order (for plots)
DOMAIN_ORDER <- c("Tumor Core", "Immune Engulfing", "Stroma", "Mammary Glands")
# Used in: 04, 05, 06, 07_01, 07_02, 08_04, 09

# Hierarchical clustering settings for merging
MERGE_DIST_METHOD <- "euclidean"
MERGE_LINK_METHOD <- "ward.D2"
# Used in: 02_03_merge_domains_08lambda.R

# ------------------------------------------------------------------------------
# 8. CELL TYPE DEFINITIONS
# ------------------------------------------------------------------------------

# All immune cell types (11)
IMMUNE_CELLTYPES <- c(
  "CD8+ T", "Treg", "NK",
  "Macrophage (Cxcl16+)", "Macrophage (Mrc1+)", "Macrophage (Arg1+)",
  "cDC1", "pDC", "DC (Ccr7+)",
  "Neutrophil", "Mast Cell"
)
# Used in: 09_proportion_analysis.R, 11_spatial_dirichlet.R

# Immune cell types used in domain-specific DE (8)
IMMUNE_CELLTYPES_DE <- c(
  "CD8+ T", "Treg", "NK",
  "Macrophage (Cxcl16+)", "cDC1", "pDC",
  "DC (Ccr7+)", "Neutrophil"
)
# Used in: 05, 06, 08_04, 13

# Epithelial / stromal cell types
EPITHELIAL_STROMAL_CELLTYPES <- c(
  "Tumor", "Basal/Luminal", "Fibroblast",
  "Myofibroblast", "Endothelial Cell", "Adipocyte"
)
# Used in: 11_spatial_dirichlet.R

# DC subtypes
DC_TYPES <- c("cDC1", "pDC", "DC (Ccr7+)")
# Used in: 10_minor_diagnostics.R

# ------------------------------------------------------------------------------
# 9. GENES OF INTEREST
# ------------------------------------------------------------------------------
GENES_CYTOTOXICITY <- c("Klrk1", "Gzmb", "Prf1", "Nkg7")
GENES_CHEMOKINE    <- c("Cxcl16", "Cxcr6")
GENES_CYTOKINE     <- c("Il15", "Osm", "Osmr", "Il24")
GENES_IMMUNE       <- c("Cd8a", "Cd4", "Foxp3", "Cd3d", "Ptprc")
GENES_MACROPHAGE   <- c("Mrc1", "Cd163", "Nos2", "Arg1")
GENES_TUMOR        <- c("Mki67", "Top2a", "Tigit")
GENES_DEGS         <- c("Pdgfra", "Fcgr4", "Lig1")
# Used in: 08_01, 08_02, 08_03, 08_04 (spatial/violin plots)

# Combined gene list for spatial visualization
GENES_SPATIAL_VIS <- c(
  GENES_CYTOTOXICITY, GENES_CHEMOKINE, GENES_CYTOKINE,
  GENES_IMMUNE, GENES_MACROPHAGE, GENES_TUMOR, GENES_DEGS
)
# Used in: 08_01, 08_03 (full spatial feature maps)

# ------------------------------------------------------------------------------
# 10. DOMAIN COMPARISONS (for pairwise DE)
# ------------------------------------------------------------------------------
DOMAIN_COMPARISONS <- list(
  c("Tumor_Core", "Immune_Engulfing"),
  c("Immune_Engulfing", "Stroma"),
  c("Tumor_Core", "Stroma")
)
# Used in: 05_DEG_DESeq_cellxdomains.R, 06_DEG_DESeq_cellxdomain_DKOvsWT.R

# ------------------------------------------------------------------------------
# 11. SPATIAL MASK DEFAULTS
# ------------------------------------------------------------------------------
# Script 08_01: mask tumor core, show all cell types
MASK_DOMAIN_DEFAULT   <- c("2")
MASK_CELLTYPE_DEFAULT <- c()   # empty = all
# Used in: 08_01_spatial_masks_08lambda.R

# Script 08_02: mask stroma, keep tumor cells
MASK_DOMAIN_VIOLIN   <- "4"
MASK_CELLTYPE_VIOLIN <- "Tumor"
# Used in: 08_02_violin_masks.R

# Script 08_04: domains for cell-type x domain violins
MASK_DOMAINS_VIOLIN_X <- c("2", "4")
# Used in: 08_04_violin_x_cells.R

# ------------------------------------------------------------------------------
# 12. PATHWAY ANALYSIS
# ------------------------------------------------------------------------------
MIN_PATH_SIZE <- 10            # Minimum genes per pathway
# Used in: 12_pathway_activity.R, 13_DE_pathways_immuneXdomains.R

TOP_PATHWAYS_PLOT <- 30        # Number of pathways to display
# Used in: 12, 13

PATHWAY_NAME_TRUNC <- 70      # Max characters for pathway names
# Used in: 12, 13

# ------------------------------------------------------------------------------
# 13. VISUALIZATION DEFAULTS
# ------------------------------------------------------------------------------
VOLCANO_TOP_N <- 20            # Number of top genes to label in volcano plots
# Used in: 03, 04, 05, 06

PANEL_NCOL <- 4                # Default columns for multi-panel plots
# Used in: 08_02, 08_04

THEME_BASE_SIZE <- 11          # Base font size for theme_publication()
# Used in: 00_color_palette.R (theme_publication default arg)

# Heatmap settings
HEATMAP_TOP_N         <- 30
HEATMAP_CLUSTER_ROWS  <- TRUE
HEATMAP_CLUSTER_COLS  <- FALSE
HEATMAP_LOG_TRANSFORM <- TRUE   # log2(CPM + 1)
HEATMAP_SCALE_ROWS    <- TRUE
# Used in: 03, 04 (DEG heatmaps)
