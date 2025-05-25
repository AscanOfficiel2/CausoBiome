################
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
#######################################
# Load data files
exp_matrix <- read.csv("General_species_matrix.csv", row.names = 1)
metadata <- read.csv("CRC_Species_contaminants_metadata.csv")

# Data Normalisation using vegan
counts <- decostand(exp_matrix, method =  "hellinger")


# Set Sample_ID as row names in metadata
rownames(metadata) <- metadata$Sample_ID

# Transpose Hellinger matrix so rows = samples
counts_t <- t(counts)

# Align metadata and expression matrix based on common Sample_IDs
common_samples <- intersect(rownames(counts_t), rownames(metadata))

# Subset both matrices to shared samples
counts_t <- counts_t[common_samples, ]
metadata <- metadata[common_samples, ]


write.csv(metadata, "Aligned_metadata.csv", row.names = TRUE)



######## use this when there is error ###################
#zero_var_rows <- apply(counts, 1, var) == 0
#counts <- counts[!zero_var_rows, ]
#counts <- counts[, metadata$Sample_ID, drop = FALSE]
###########################################################################################
#####Inspecting the original Data

# Step 2: Perform PCA (Before)
pca <- prcomp(counts_t , scale. = TRUE)

# Combine PCA results with metadata for plotting
pca_data <- as.data.frame(pca$x)
pca_data$Project <- metadata$Project
pca_data$Country <- metadata$Country
pca_data$Center_Name <- metadata$Center_Name
pca_data$Age <- metadata$Age
pca_data$Group <- metadata$Group
pca_data$Continent <- metadata$Continent  
pca_data$Instrument <- metadata$Instrument
pca_data$BMI <- metadata$BMI
pca_data$Sex <- metadata$Sex


# Create PCA plots
pca_plot_project <- xyplot(PC2 ~ PC1, data = pca_data, groups = Project,
                           auto.key = FALSE, main = "Project", xlab = "PC1", ylab = "PC2",
                           par.settings = list(axis.line = list(col = "black")))

pca_plot_country <- xyplot(PC2 ~ PC1, data = pca_data, groups = Country,
                           auto.key = FALSE, main = "Country", xlab = "PC1", ylab = "PC2",
                           par.settings = list(axis.line = list(col = "black")))

pca_plot_center_name <- xyplot(PC2 ~ PC1, data = pca_data, groups = Center_Name,
                               auto.key = FALSE, main = "Center_Name", xlab = "PC1", ylab = "PC2",
                               par.settings = list(axis.line = list(col = "black")))

pca_plot_continent <- xyplot(PC2 ~ PC1, data = pca_data, groups = Continent,
                             auto.key = FALSE, main = "Continent", xlab = "PC1", ylab = "PC2",
                             par.settings = list(axis.line = list(col = "black")))

pca_plot_instrument <- xyplot(PC2 ~ PC1, data = pca_data, groups = Instrument,
                              auto.key = FALSE, main = "Instrument", xlab = "PC1", ylab = "PC2",
                              par.settings = list(axis.line = list(col = "black")))

pca_plot_age <- xyplot(PC2 ~ PC1, data = pca_data, groups = Age,
                       auto.key = FALSE, main = "Age", xlab = "PC1", ylab = "PC2",
                       par.settings = list(axis.line = list(col = "black")))

pca_plot_group <- xyplot(PC2 ~ PC1, data = pca_data, groups = Group,
                         auto.key = FALSE, main = "Group", xlab = "PC1", ylab = "PC2",
                         par.settings = list(axis.line = list(col = "black")))

pca_plot_BMI <- xyplot(PC2 ~ PC1, data = pca_data, groups = BMI,
                       auto.key = FALSE, main = "BMI", xlab = "PC1", ylab = "PC2",
                       par.settings = list(axis.line = list(col = "black")))

pca_plot_sex <- xyplot(PC2 ~ PC1, data = pca_data, groups = Sex,
                       auto.key = FALSE, main = "Sex", xlab = "PC1", ylab = "PC2",
                       par.settings = list(axis.line = list(col = "black")))


# Helper function to create plots with centered titles in a box, grid lines, tick modifications, slanted x-axis labels, and reduced axis label sizes
create_pca_plot <- function(data, group_var, title) {
  xyplot(PC2 ~ PC1, data = data, groups = group_var,
         auto.key = FALSE,
         xlab = list("PC1", cex = 0.6, font = 2),  # Reduce x-axis title font size
         ylab = list("PC2", cex = 0.6, font = 2),  # Reduce y-axis title font size
         scales = list(
           x = list(rot = 45, cex = 0.6, font = 2),  # Rotate x-axis labels and reduce font size
           y = list(cex = 0.6, font = 2)  # Reduce y-axis labels font size
         ),
         par.settings = list(axis.line = list(col = "black")),
         axis.components = list(
           top = list(tck = 0.5),  # No tick marks on top axis
           right = list(tck = 0.5),  # No tick marks on right axis
           left = list(tck = 0.5),  # Maintain tick marks on left axis
           bottom = list(tck = 1)  # Maintain tick marks on bottom axis
         ),
         xlab.top = list(title, cex = 0.6, col = "black", border = TRUE, font = 2),  # Reduce font size and make not bold
         panel = function(...) {
           panel.grid(h = -1, v = -1, lty = 2, col = alpha("black", 0.3))  # Dashed grid lines with increased alpha
           panel.superpose(...)
         })
}

# Save the plot as a tiff file
tiff("CRC_species_contaminants_MetaData_PCA_batch_sex.tiff", width = 8, height = 10, units = "in", res = 600)


# Create PCA plots
pca_plot_project <- create_pca_plot(pca_data, pca_data$Project, "Project")
pca_plot_country <- create_pca_plot(pca_data, pca_data$Country, "Country")
pca_plot_age <- create_pca_plot(pca_data, pca_data$Age, "Age")
pca_plot_continent <- create_pca_plot(pca_data, pca_data$Continent, "Continent")
pca_plot_instrument <- create_pca_plot(pca_data, pca_data$Instrument, "Instrument")
pca_plot_center_name<- create_pca_plot(pca_data, pca_data$Center_Name, "Center_Name")
pca_plot_BMI <- create_pca_plot(pca_data, pca_data$BMI, "BMI")
pca_plot_Group <- create_pca_plot(pca_data, pca_data$Group, "Group")
pca_plot_sex <- create_pca_plot(pca_data, pca_data$Sex, "Sex")



# Arrange all PCA plots in a 3x3 grid using lattice's print + split layout
print(pca_plot_Group,       split = c(1, 1, 3, 3), more = TRUE)
print(pca_plot_age,         split = c(2, 1, 3, 3), more = TRUE)
print(pca_plot_sex,         split = c(3, 1, 3, 3), more = TRUE)
print(pca_plot_BMI,         split = c(1, 2, 3, 3), more = TRUE)
print(pca_plot_project,     split = c(2, 2, 3, 3), more = TRUE)
print(pca_plot_continent,   split = c(3, 2, 3, 3), more = TRUE)
print(pca_plot_country,     split = c(1, 3, 3, 3), more = TRUE)
print(pca_plot_instrument,  split = c(2, 3, 3, 3), more = TRUE)
print(pca_plot_center_name, split = c(3, 3, 3, 3), more = FALSE)  # Last plot: more = FALSE


#Close
dev.off()

# ===== Load Required Libraries =====Batch correction##############
library(vegan)       # For Hellinger transformation
library(sva)         # For ComBat
library(limma)       # Optional alternative for batch correction
library(dplyr)       # Data wrangling
library(readr)       # Faster CSV
library(tibble)      # For rownames_to_column

# ===== Load Data =====
species_matrix <- read_csv("General_species_matrix.csv")
metadata <- read_csv("CRC_Species_contaminants_metadata.csv")

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



# ===== PERMANOVA Before Correction =====
permanova_before <- adonis2(species_matrix ~ Group, data = metadata, method = "bray")
print("PERMANOVA BEFORE batch correction:")
print(permanova_before)

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


# ===== Hellinger Transformation =====
hellinger_matrix <- decostand(combat_corrected, method = "hellinger")

hellinger_matrix_t <- t(hellinger_matrix)


# ===== Save Outputs =====
write.csv(hellinger_matrix_t, "hellinger_transformed.csv")
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