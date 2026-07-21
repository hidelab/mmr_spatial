# 00_color_palette.R
# Shared color definitions for the MMR Spatial project.
# Source this file in any R script that produces plots for consistent styling.
#
# Design: colors are grouped by biological lineage so that related cell types
# share the same hue family but differ in saturation/brightness.

# --- Tier1 cell type palette --------------------------------------------------
# Grouped by lineage:
#   Tumor/Epithelial  = reds
#   Lymphoid (cyto)   = blues
#   Lymphoid (reg)    = teal
#   Macrophages       = greens
#   Dendritic cells   = purples
#   Myeloid (other)   = oranges/ambers
#   Stromal           = earth tones

TIER1_COLORS <- c(
  # Tumor / Epithelial
  "Tumor"                 = "#C62828",
  "Basal/Luminal"         = "#EF9A9A",

  # Lymphoid
  "CD8+ T"               = "#1565C0",
  "NK"                   = "#64B5F6",
  "Treg"                 = "#81D4FA",

  # Macrophages
  "Macrophage (Arg1+)" = "#2E7D32",
  "Macrophage (Mrc1+)"   = "#81C784",
  "Macrophage (Cxcl16+)"   = "#A5D6A7",

  # Dendritic cells
  "cDC1"                 = "#6A1B9A",
  "pDC"                  = "#BA68C8",
  "DC (Ccr7+)"           = "#CE93D8",

  # Granulocytes
  "Neutrophil"           = "#E65100",
  "Mast Cell"            = "#FFB74D",

  # Stromal / Structural
  "Fibroblast"           = "#4E342E",
  "Myofibroblast"        = "#A1887F",
  "Endothelial Cell"     = "#78909C",
  "Adipocyte"            = "#827717"
)

# Canonical ordering for legends (matches the grouping above, dark → light within groups)
TIER1_ORDER <- names(TIER1_COLORS)

# --- Condition palette --------------------------------------------------------
CONDITION_COLORS <- c(

  "WT"  = "#2196F3",
  "DKO" = "#F44336"
)

# --- Sample palette -----------------------------------------------------------
SAMPLE_COLORS <- c(
  "WT1"  = "#90CAF9",
  "WT2"  = "#42A5F5",
  "WT3"  = "#1565C0",
  "DKO1" = "#EF9A9A",
  "DKO2" = "#EF5350",
  "DKO3" = "#C62828",
  "DKO4" = "#880E4F"
)

# --- Publication theme --------------------------------------------------------
theme_publication <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0.5),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.title = element_text(size = base_size, face = "bold"),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size - 1),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.background = element_blank()
    )
}
