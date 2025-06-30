set.seed(45)
# Load required libraries
library(vegan)
library(ggplot2)
library(car)
library(gridExtra)
library(lme4)
library(MuMIn)
library(lattice)
library(scales)
library(limma)
library(umap)
library(Rtsne)
library(dplyr)
library(shadowtext)

########################################
# STEP 1: Load and Merge Hellinger ARG + VF Matrices
########################################

arg <- read.csv("CRC_CARD_expression_matrix_hellinger.csv", row.names = 1)
vf  <- read.csv("CRC_VFDB_expression_matrix_hellinger.csv", row.names = 1)

# Prefix gene names to distinguish ARG and VF
rownames(arg) <- paste0("ARG_", rownames(arg))
rownames(vf)  <- paste0("VF_", rownames(vf))

# Keep only samples present in both datasets
common_samples <- intersect(colnames(arg), colnames(vf))
arg_aligned <- arg[, common_samples]
vf_aligned  <- vf[, common_samples]

# Merge and transpose
counts <- rbind(arg_aligned, vf_aligned)
count_t <- t(counts)  # Samples as rows, genes as columns

########################################
# STEP 2: Load and Align Metadata
########################################

metadata <- read.csv("CRC_vf_arg_meta_main.csv")
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[rownames(count_t), ]

#  Save the Aligned Metadata
metadata_export <- metadata
metadata_export$Sample_ID <- rownames(metadata_export)

# Save
write.csv(metadata_export, "Metadata_Aligned_to_CountMatrix.csv", row.names = FALSE)
#########################################################################
# ===== Confounding Check (Cramér’s V) =====
group_batch_table <- table(metadata$Group, metadata$Project)
write.csv(as.data.frame.matrix(group_batch_table), "Group_vs_Project_table.csv")
cramer_v <- CramerV(group_batch_table)
cat("Cramér's V between Group and Project:", round(cramer_v, 3), "\n")

if (cramer_v < 0.10) {
  interpretation <- "Negligible association — batch correction is not needed."
  action <- "Skip batch correction."
} else if (cramer_v >= 0.10 & cramer_v < 0.30) {
  interpretation <- "Weak to moderate association — batch correction decision should be based on PERMANOVA."
  action <- "⚠️ Run PERMANOVA to assess if Project explains more variance than Group."
} else if (cramer_v >= 0.30 & cramer_v < 0.50) {
  interpretation <- "Moderate association — batch correction is recommended."
  action <- "Proceed with batch correction."
} else {
  interpretation <- "Strong association — batch correction is strongly recommended."
  action <- "Proceed with batch correction."
}

# Print results
cat("\n Cramér's V:", round(cramer_v, 3), "\n")
cat(" Interpretation:", interpretation, "\n")
cat(" Recommended action:", action, "\n")

########################################
# STEP 3: Initial PCA + PERMANOVA (Uncorrected)
########################################

# Remove zero-variance genes
count_t_clean <- count_t[, apply(count_t, 2, var) > 0]

# PCA
pca_before <- prcomp(count_t_clean, scale. = TRUE)
pca_data_before <- as.data.frame(pca_before$x)
pca_data_before$Group <- metadata$Group


# Individual PERMANOVAs
adonis_project    <- adonis2(count_t_clean ~ Project, data = metadata, method = "euclidean")
adonis_instrument <- adonis2(count_t_clean ~ Instrument, data = metadata, method = "euclidean")
adonis_center_name     <- adonis2(count_t_clean ~ Center_Name, data = metadata, method = "euclidean")
adonis_country    <- adonis2(count_t_clean ~ Country, data = metadata, method = "euclidean")
adonis_continent  <- adonis2(count_t_clean ~ Continent, data = metadata, method = "euclidean")
adonis_bmi     <- adonis2(count_t_clean ~ BMI , data = metadata, method = "euclidean")
adonis_age  <- adonis2(count_t_clean ~ Age, data = metadata, method = "euclidean")
adonis_sex   <- adonis2(count_t_clean ~ Sex, data = metadata, method = "euclidean")

# Extract each result into a named data frame with metadata
extract_adonis <- function(result, factor_name) {
  df <- as.data.frame(result)
  df$Factor <- factor_name
  return(df)
}

# Combine all into one dataframe
full_batch_adonis_df <- rbind(
  extract_adonis(adonis_project, "Project"),
  extract_adonis(adonis_instrument, "Instrument"),
  extract_adonis(adonis_center_name, "Center_Name"),
  extract_adonis(adonis_country, "Country"),
  extract_adonis(adonis_continent, "Continent"),
  extract_adonis(adonis_age, "Age"),
  extract_adonis(adonis_bmi, "BMI"),
  extract_adonis(adonis_sex, "Sex")
)

# Reorder columns for clarity
full_batch_adonis_df <- full_batch_adonis_df[, c("Factor", names(full_batch_adonis_df)[1:5])]

# Save full result to CSV
write.csv(full_batch_adonis_df, "Batch_Factor_VF_ARG_FULL_PERMANOVA.csv", row.names = TRUE)


########################################
# STEP 4: Batch Correction
########################################

# Transpose for limma (genes as rows, samples as columns)
counts_matrix <- t(count_t_clean)

# Design matrix to preserve Group (biological signal)
design <- model.matrix(~ metadata$Group)

covariates <- model.matrix(~ metadata$Instrument + metadata$Center_Name)[, -1]

# Apply batch correction (remove Project_ID effect)
counts_corrected <- removeBatchEffect(
  counts_matrix,
  batch = metadata$Project,
  covariates = covariates,
  design = design
)

# Save the corrected counts
write.csv(counts_corrected, "Combined_ARG_VFDB_batch_corrected.csv", row.names = TRUE)

# Transpose back to samples x genes
counts_corrected <- t(counts_corrected)

#For UMAP PLOT
counts_matrix_corrected <- t(counts_corrected)
############################################
# === PERMANOVA BEFORE Batch Correction ===
permanova_before <- adonis2(count_t_clean ~ Group, data = metadata, method = "euclidean")
cat("PERMANOVA BEFORE batch correction:\n")
print(permanova_before)

# Save to CSV
permanova_before_df <- as.data.frame(permanova_before)
write.csv(permanova_before_df, "permanova_before_batch_correction.csv", row.names = TRUE)
################################################################################
# === PERMANOVA AFTER Batch Correction ===
permanova_after <- adonis2(counts_corrected ~ Group, data = metadata, method = "euclidean")
cat("PERMANOVA AFTER batch correction:\n")
print(permanova_after)

# Save to CSV
permanova_after_df <- as.data.frame(permanova_after)
write.csv(permanova_after_df, "permanova_after_batch_correction.csv", row.names = TRUE)

######################################################
# Step 1: Perform UMAP (Before and After Batch Correction)

# UMAP for data before batch correction
set.seed(42)
umap_before <- umap(t(counts_matrix), n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_data_before <- as.data.frame(umap_before$layout)
umap_data_before$Group <- metadata$Group
colnames(umap_data_before) <- c("UMAP1", "UMAP2", "Group")

# UMAP for data after batch correction
set.seed(42)
umap_after <- umap(t(counts_matrix_corrected), n_neighbors = 15, min_dist = 0.1, metric = "euclidean")
umap_data_after <- as.data.frame(umap_after$layout)
umap_data_after$Group <- metadata$Group
colnames(umap_data_after) <- c("UMAP1", "UMAP2", "Group")

####################################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)

# Add type labels for plotting
umap_data_before$Type <- "Before Batch Correction"
umap_data_after$Type  <- "After Batch Correction"

# Combine both into a single dataframe
umap_combined <- rbind(umap_data_before, umap_data_after)

# Ensure proper column order and factors
colnames(umap_combined)[1:2] <- c("UMAP1", "UMAP2")
umap_combined$Type <- factor(umap_combined$Type, levels = c("Before Batch Correction", "After Batch Correction"))
umap_combined$Group <- factor(umap_combined$Group, levels = c("Healthy", "Adenoma", "Cancer"))

# Define custom colors
group_colors <- c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")

# Plot and save as TIFF
tiff("UMAP_Before_After_Correction.tiff", width = 8, height = 4.5, units = "in", res = 600)
ggplot(umap_combined, aes(x = UMAP1, y = UMAP2, fill = Group)) +
  geom_point(shape = 21, size = 2.8, stroke = 0.3, color = "black", alpha = 0.85) +
  facet_wrap(~ Type) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 16) +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    strip.background = element_rect(fill = "#f0f0f0", color = NA),
    strip.text = element_text(face = "bold", size = 14),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right"
  )
dev.off()


######Corrected the batch corrected file replacing negatives with zeros for following steps

combined_ARG_VFDB <- read.csv("Combined_ARG_VFDB_batch_corrected.csv", row.names = 1)

# Set negative values to zero
combined_ARG_VFDB <- combined_ARG_VFDB
combined_ARG_VFDB[combined_ARG_VFDB < 0] <- 0

# Save the results to a single CSV file
write.csv(combined_ARG_VFDB, file = "combined_ARG_VFDB_final_DATA.csv", row.names = TRUE)
########################