##################### INITIAL SETUP #####################
set.seed(45)

library(vegan)
library(ggplot2)
library(car)
library(gridExtra)
library(lme4)
library(MuMIn)
library(lattice)
library(scales)
library(limma)
library(DescTools)
library(sva)
library(readr)
library(tibble)
library(umap)
library(ggpubr)
library(dplyr)

# ===== Load Data =====
exp_matrix <- read.csv("Species_subset_matrix.csv", row.names = 1)
metadata <- read.csv("crc_meta.csv")  

#######################################
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

# ===== Hellinger Transformation =====
hellinger_matrix <- decostand(exp_matrix, method = "hellinger")
hellinger_matrix_t <- t(hellinger_matrix)

# ===== Align Metadata =====
rownames(metadata) <- metadata$Sample_ID
common_samples <- intersect(rownames(hellinger_matrix_t), rownames(metadata))
hellinger_matrix_t <- hellinger_matrix_t[common_samples, ]
metadata <- metadata[common_samples, ]
write.csv(metadata, "Aligned_metadata_Taxonomy.csv", row.names = TRUE)

##################### PCA Diagnostics #####################
# Remove zero-variance columns
zero_var_cols <- apply(hellinger_matrix_t, 2, function(x) var(x, na.rm = TRUE) == 0)
counts_t_filtered <- hellinger_matrix_t[ , !zero_var_cols]

# PCA
pca <- prcomp(counts_t_filtered, scale. = TRUE)
pca_data <- as.data.frame(pca$x)
pca_data <- cbind(pca_data, metadata)

# Reusable PCA Plot Function
create_pca_plot <- function(data, group_var, title) {
  xyplot(PC2 ~ PC1, data = data, groups = group_var,
         auto.key = FALSE,
         xlab = list("PC1", cex = 0.6, font = 2),
         ylab = list("PC2", cex = 0.6, font = 2),
         scales = list(x = list(rot = 45, cex = 0.6), y = list(cex = 0.6)),
         par.settings = list(axis.line = list(col = "black")),
         axis.components = list(top = list(tck = 0.5), right = list(tck = 0.5),
                                left = list(tck = 0.5), bottom = list(tck = 1)),
         xlab.top = list(title, cex = 0.6, col = "black", border = TRUE, font = 2),
         panel = function(...) {
           panel.grid(h = -1, v = -1, lty = 2, col = alpha("black", 0.3))
           panel.superpose(...)
         })
}

# Save PCA plots
tiff("CRC_species_contaminants_MetaData_PCA_batch_sex.tiff", width = 8, height = 10, units = "in", res = 600)
print(create_pca_plot(pca_data, pca_data$Group, "Group"),       split = c(1, 1, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$Age, "Age"),           split = c(2, 1, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$Sex, "Sex"),           split = c(3, 1, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$BMI, "BMI"),           split = c(1, 2, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$Project, "Project"),   split = c(2, 2, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$Continent, "Continent"), split = c(3, 2, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$Country, "Country"),   split = c(1, 3, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$Instrument, "Instrument"), split = c(2, 3, 3, 3), more = TRUE)
print(create_pca_plot(pca_data, pca_data$Center_Name, "Center_Name"), split = c(3, 3, 3, 3), more = FALSE)
dev.off()

#############################################
variables_to_test <- c("Project", "Center_Name", "Instrument", "Country", "Continent", "Age", "BMI", "Sex")
permanova_results <- lapply(variables_to_test, function(var) {
  formula <- as.formula(paste("hellinger_matrix_t ~", var))
  result <- adonis2(formula, data = metadata, method = "bray")
  data.frame(
    Variable = var,
    R2 = round(result$R2[1], 4),
    P_value = round(result$`Pr(>F)`[1], 4)
  )
})
permanova_df <- do.call(rbind, permanova_results)
permanova_df <- permanova_df[order(-permanova_df$R2), ]
print(permanova_df)
write.csv(permanova_df, "PERMANOVA_Metadata_Association.csv", row.names = FALSE)

# ===== Load Required Libraries =====Batch correction##############
library(vegan)       # For Hellinger transformation
library(sva)         # For ComBat
library(limma)       # Optional alternative for batch correction
library(dplyr)       # Data wrangling
library(readr)       # Faster CSV
library(tibble)      # For rownames_to_column

# ===== Load Data =====
species_matrix <- read_csv("subset_species_Ecology_matrix.csv")
metadata <- read_csv("crc_meta.csv")

# ===== Prepare Species Matrix =====
species_matrix <- as.data.frame(species_matrix)
rownames(species_matrix) <- species_matrix$Taxonomy_name
species_matrix <- species_matrix[ , !(names(species_matrix) %in% "Taxonomy_name")]
species_matrix <- t(species_matrix)  # transpose: rows = samples

# Clean sample IDs
rownames(species_matrix) <- trimws(rownames(species_matrix))
metadata$Sample_ID <- trimws(as.character(metadata$Sample_ID))

# Remove NA and duplicates from metadata
metadata <- metadata[!is.na(metadata$Sample_ID), ]
metadata <- metadata[!duplicated(metadata$Sample_ID), ]
rownames(metadata) <- metadata$Sample_ID

# Align metadata and species matrix
common_samples <- intersect(rownames(species_matrix), metadata$Sample_ID)
species_matrix <- species_matrix[common_samples, ]
metadata <- metadata[common_samples, ]
###########################################################################################
# ===== PERMANOVA Before Correction =====
permanova_before <- adonis2(species_matrix ~ Group, data = metadata, method = "bray")
# Print to console
cat("PERMANOVA BEFORE batch correction:\n")
print(permanova_before)
# Save to CSV
permanova_before_df <- as.data.frame(permanova_before)
write.csv(permanova_before_df, "permanova_before_correction.csv", row.names = TRUE)


# ===== Batch Correction with ComBat (on Hellinger matrix) =====
modcombat <- model.matrix(~ Group, data = metadata)
#batch <- metadata$Project
metadata$Combined_Batch <- factor(paste(metadata$Project, metadata$Center_Name, sep = "_"))
batch <- metadata$Combined_Batch
combat_corrected <- ComBat(dat = t(species_matrix), batch = batch, mod = modcombat, par.prior = TRUE)
combat_corrected <- t(combat_corrected)

# Remove negative values for Bray–Curtis (if needed)
min_val <- min(combat_corrected)
if (min_val < 0) {
  combat_corrected <- combat_corrected + abs(min_val)
}

# ===== PERMANOVA After Correction =====
permanova_after <- adonis2(combat_corrected ~ Group, data = metadata, method = "bray")
print("PERMANOVA AFTER batch correction:")
print(permanova_after)
# Save to CSV
permanova_after_df <- as.data.frame(permanova_after)
write.csv(permanova_after_df, "permanova_after_correction.csv", row.names = TRUE)

# ===== Hellinger Transformation =====
hellinger_matrix <- decostand(combat_corrected, method = "hellinger")

hellinger_matrix_t <- t(hellinger_matrix)


# ===== Save Outputs =====
write.csv(hellinger_matrix_t, "hellinger_transformed_combat.csv")
write.csv(combat_corrected, "combat_corrected_matrix.csv")


# ===== Load Required Libraries =====
library(umap)
library(ggplot2)
library(dplyr)
library(ggpubr)

# ===== Run UMAP BEFORE Correction =====
set.seed(123)
umap_before <- umap(species_matrix)
umap_df_before <- as.data.frame(umap_before$layout)
umap_df_before$Group <- metadata$Group
umap_df_before$Type <- "Before Batch Correction"

# ===== Run UMAP AFTER Correction =====
set.seed(123)
umap_after <- umap(combat_corrected)
umap_df_after <- as.data.frame(umap_after$layout)
umap_df_after$Group <- metadata$Group
umap_df_after$Type <- "After Batch Correction"

# ===== Combine & Rename =====
umap_combined <- rbind(
  rename(umap_df_before, UMAP1 = V1, UMAP2 = V2),
  rename(umap_df_after, UMAP1 = V1, UMAP2 = V2)
)

# ===== Ensure Facet Order =====
umap_combined$Type <- factor(umap_combined$Type,
                             levels = c("Before Batch Correction", "After Batch Correction"))

# ===== Define Colors =====
group_colors <- c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")

# ===== Plot UMAP with Enhanced Aesthetics =====
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
#############################End ########################################################

