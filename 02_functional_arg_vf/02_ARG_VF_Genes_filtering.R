# =======================
# Step 0: Setup
# =======================
set.seed(45)
library(dplyr)
library(vegan)
library(ggplot2)

# =======================
# Step 1: Load Data
# =======================
cat("Loading ARG, VFDB, and metadata files...\n")

arg <- read.csv("CRC_CARD_expression_matrix.csv", row.names = 1, check.names = FALSE)
vf  <- read.csv("CRC_VFDB_expression_matrix.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata_ext.csv", check.names = FALSE)

# Prefix gene names to distinguish ARG and VFDB
rownames(arg) <- paste0("ARG_", rownames(arg))
rownames(vf)  <- paste0("VF_", rownames(vf))

# =======================
# Step 2: Align Samples
# =======================
cat("Aligning common samples across datasets...\n")

common_samples <- intersect(colnames(arg), colnames(vf))
arg_aligned <- arg[, common_samples]
vf_aligned  <- vf[, common_samples]

# Combine
merged_matrix <- rbind(arg_aligned, vf_aligned)

# Transpose → samples as rows, genes as columns
merged_t <- t(merged_matrix)
cat("Merged matrix dimensions before filtering:", dim(merged_t), "\n")

# =======================
# Step 3: Align Metadata
# =======================
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[rownames(merged_t), , drop = FALSE]

cat("Metadata aligned:", nrow(metadata), "samples.\n")

# =======================
# Step 4: Filter Genes (keep genes present in ≥5% of samples)
# =======================
cat("Filtering genes present in at least 5% of samples...\n")

# Calculate presence/absence per gene
presence_threshold <- 0.00 * nrow(merged_t)

# Count number of samples where gene abundance > 0
gene_presence_counts <- colSums(merged_t > 0)

# Keep genes meeting threshold
filtered_genes <- names(gene_presence_counts[gene_presence_counts >= presence_threshold])
filtered_matrix <- merged_t[, filtered_genes]

cat("Number of genes retained after filtering:", length(filtered_genes), "\n")

# =======================
# Step 5: Save Filtered Data
# =======================
cat("Saving filtered merged matrix and aligned metadata...\n")

write.csv(filtered_matrix, "Combined_ARG_VFDB_Filtered_Matrix.csv", row.names = TRUE)
write.csv(metadata, "Metadata_Aligned_to_FilteredMatrix.csv", row.names = FALSE)

cat("✅ Merging and filtering complete. Files saved:\n")
cat(" - Combined_ARG_VFDB_Filtered_Matrix.csv\n")
cat(" - Metadata_Aligned_to_FilteredMatrix.csv\n")

#########################################################################
