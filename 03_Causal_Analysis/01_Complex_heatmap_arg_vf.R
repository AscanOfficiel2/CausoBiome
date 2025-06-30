# ==============================================
# heatmap_consistent_biomarkers.R
# ==============================================

# 0) Install missing packages if needed
pkgs <- c("ComplexHeatmap","circlize","dplyr","readr","tibble","grid")
install_if_missing <- function(p) if(!requireNamespace(p,quietly=TRUE)) install.packages(p)
invisible(sapply(pkgs, install_if_missing))

# 1) Load libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)
library(tibble)
library(grid)  # for unit()

# 2) Read input files
expr     <- read.csv("combined_ARG_VFDB_final_DATA.csv",        stringsAsFactors=FALSE, check.names=FALSE)
meta     <- read.csv("Metadata_Aligned_VF_ARG_CountMatrix.csv", stringsAsFactors=FALSE)
filtered <- read.csv("Filtered_biomarker_annotated.csv",        stringsAsFactors=FALSE)

# 3) Rename first column to "Feature"
colnames(expr)[1] <- "Feature"

# 4) Subset to filtered biomarkers
present <- intersect(filtered$Feature, expr$Feature)
if(length(present)==0) stop("No filtered features found!")
expr_sub <- expr[expr$Feature %in% present, ]

# 5) Extract the "Genes" code (second underscore segment)
expr_sub$Genes <- sapply(strsplit(expr_sub$Feature, "_"), `[`, 2)

# 6) Build numeric matrix (rows=Feature, cols=samples)
sample_cols <- setdiff(colnames(expr_sub), c("Feature","Genes"))
mat <- as.matrix(expr_sub[, sample_cols])
rownames(mat) <- expr_sub$Feature

# 7) Z-score normalize each row
mat_z <- t(scale(t(mat)))

# 8) Prepare row annotation DF
annot_df <- data.frame(
  Feature   = rownames(mat_z),
  Genes     = expr_sub$Genes[match(rownames(mat_z), expr_sub$Feature)],
  Mechanism = filtered$Mechanism[match(rownames(mat_z), filtered$Feature)],
  stringsAsFactors = FALSE
)
annot_df$Genes     <- factor(annot_df$Genes,   levels=unique(annot_df$Genes))
annot_df$Mechanism <- factor(annot_df$Mechanism, levels=unique(annot_df$Mechanism))

row_anno <- rowAnnotation(
  Genes     = annot_df$Genes,
  Mechanism = annot_df$Mechanism,
  col = list(
    Genes     = setNames(rainbow(nlevels(annot_df$Genes)),   levels(annot_df$Genes)),
    Mechanism = setNames(rainbow(nlevels(annot_df$Mechanism)), levels(annot_df$Mechanism))
  ),
  width = unit(4, "cm")
)

# 9) Prepare column annotation (sample Group)
meta$Sample_ID <- as.character(meta$Sample_ID)
samples_keep  <- colnames(mat_z)
group_df      <- meta %>%
  filter(Sample_ID %in% samples_keep) %>%
  select(Sample_ID, Group) %>%
  column_to_rownames("Sample_ID")
group_df      <- group_df[samples_keep, , drop=FALSE]
group_df$Group <- factor(group_df$Group, levels=c("Healthy","Adenoma","Cancer"))

top_anno <- HeatmapAnnotation(
  Group = group_df$Group,
  col = list(
    Group = c(
      Healthy = "#ffcc00",
      Adenoma = "#984ea3",
      Cancer  = "#4daf4a"
    )
  ),
  annotation_name_side = "left"
)

# 10) Draw AND SAVE the heatmap as high-res PNG
png(
  filename = "heatmap_consistent_genes.png",
  width    = 12,       # inches
  height   = 10,       # inches
  units    = "in",
  res      = 600       # DPI
)

ht <- Heatmap(
  mat_z,
  name               = "Z-score",
  show_row_names     = TRUE,
  show_column_names  = FALSE,
  cluster_rows       = TRUE,
  cluster_columns    = TRUE,
  row_names_gp       = gpar(fontsize=8),
  row_title          = "Consistent Genes",
  top_annotation     = top_anno,
  left_annotation    = row_anno
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()
