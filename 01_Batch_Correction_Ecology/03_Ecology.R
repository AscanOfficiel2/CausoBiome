###############################################################
## Continue with subse
########################################################################################################
# ============================
# Alpha Diversity: R Script
# ============================

# ==== Load Required Libraries ====
library(vegan)        # for diversity functions
library(ggplot2)      # for plotting
library(dplyr)        # for data wrangling
library(readr)        # for reading CSV
library(ggpubr)       # for publication-style plots

# ==== Load Expression Matrix ====
expr <- read_csv("Expression_Matrix_Aligned.csv")
expr <- as.data.frame(expr)

# Fix column naming and set sample IDs as rownames
colnames(expr)[1] <- "Sample_ID"
rownames(expr) <- expr$Sample_ID
expr$Sample_ID <- NULL

# Convert all values to numeric (required for diversity calculations)
expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))

# ==== Load Metadata ====
meta <- read_csv("Aligned_metadata.csv")
meta <- meta %>% filter(Sample_ID %in% rownames(expr))
rownames(meta) <- meta$Sample_ID

# Ensure row alignment between metadata and expression matrix
expr <- expr[meta$Sample_ID, ]

# ==== Compute Alpha Diversity Metrics ====
alpha_df <- data.frame(
  Sample_ID = rownames(expr),
  Richness = specnumber(expr),
  Shannon = diversity(expr, index = "shannon"),
  Simpson = diversity(expr, index = "simpson")
)

# ==== Merge with Metadata ====
alpha_df <- left_join(alpha_df, meta, by = "Sample_ID")


# ==== Print Kruskal-Wallis Results ====
print(kruskal.test(Shannon ~ Group, data = alpha_df))
print(kruskal.test(Simpson ~ Group, data = alpha_df))
print(kruskal.test(Richness ~ Group, data = alpha_df))

# ==== Run Kruskal-Wallis Tests ====
shannon_test <- kruskal.test(Shannon ~ Group, data = alpha_df)
simpson_test <- kruskal.test(Simpson ~ Group, data = alpha_df)
richness_test <- kruskal.test(Richness ~ Group, data = alpha_df)
View(alpha_df)
summary(alpha_df)
write.csv(alpha_df,"Alpha_df.csv", row.names = FALSE )

# ==== Create a Summary Table ====
kw_results <- data.frame(
  Metric = c("Shannon", "Simpson", "Richness"),
  Chi_Squared = c(shannon_test$statistic, simpson_test$statistic, richness_test$statistic),
  DF = c(shannon_test$parameter, simpson_test$parameter, richness_test$parameter),
  P_Value = c(shannon_test$p.value, simpson_test$p.value, richness_test$p.value)
)

# ==== Save to CSV ====
write.csv(kw_results, "Kruskal_Wallis_Alpha_Diversity_Results.csv", row.names = FALSE)
##############################################POSTHOC ANALYSIS##############################
# ==== Install and Load Required Packages for Post Hoc ====
if (!require(FSA)) install.packages("FSA")
library(FSA)

# ==== Post Hoc: Dunn's Test ====
# Adjusted for multiple comparisons using Holm method (can also use "bonferroni", "BH", etc.)
dunn_shannon <- dunnTest(Shannon ~ Group, data = alpha_df, method = "holm")
dunn_simpson <- dunnTest(Simpson ~ Group, data = alpha_df, method = "holm")
dunn_richness <- dunnTest(Richness ~ Group, data = alpha_df, method = "holm")

# ==== View and Save Results ====
# Extract tables from the list structure
shannon_posthoc <- dunn_shannon$res
simpson_posthoc <- dunn_simpson$res
richness_posthoc <- dunn_richness$res

# Save to CSV
write.csv(shannon_posthoc, "Dunn_Posthoc_Shannon.csv", row.names = FALSE)
write.csv(simpson_posthoc, "Dunn_Posthoc_Simpson.csv", row.names = FALSE)
write.csv(richness_posthoc, "Dunn_Posthoc_Richness.csv", row.names = FALSE)


####################PLots for Alpha diversity diversity#########################

# === Required Libraries ===
library(ggplot2)
library(dplyr)

# === Set Your Custom Colors for All Plots ===
group_colors <- c("Adenoma" = "#984ea3", "Cancer" = "#4daf4a", "Healthy" = "#ffcc00")


# Function with full axis border
plot_alpha_individual <- function(metric, y_label, filename) {
  p <- ggplot(alpha_df, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.5,
      color = "black",
      size = 1.2,
      alpha = 0.9
    ) +
    geom_jitter(
      width = 0.15,
      size = 1.5,
      alpha = 0.6,
      color = "black"
    ) +
    scale_fill_manual(values = group_colors) +
    theme_classic(base_size = 18) +  # Full axis lines
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 12,color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.text.y = element_text(size = 12, color = "black"),
      panel.grid = element_blank(),
      legend.position = "none",
      panel.border = element_rect(color = "black", size = 1.5, fill = NA),  # Full box
      plot.margin = unit(c(1, 1, 1, 1), "cm")
    ) +
    labs(y = y_label)
  
  # Export as TIFF
  tiff(filename, width = 4, height = 4, units = "in", res = 600)
  print(p)
  dev.off()
}

# === Save Individual Plots ===
plot_alpha_individual("Shannon", "Shannon Diversity", "Alpha_Shannon_BorderBox.tiff")
plot_alpha_individual("Simpson", "Simpson Diversity", "Alpha_Simpson_BorderBox.tiff")
plot_alpha_individual("Richness", "Species Richness", "Alpha_Richness_BorderBox.tiff")

cat("✅ Plots saved with fully enclosed border using axis box.\n")
#####################################################################################################
## ================================
# Beta Diversity: NMDS + PERMANOVA + PERMDISP
# ================================

group_colors <- c("Adenoma" = "#984ea3", "Cancer" = "#4daf4a", "Healthy" = "#ffcc00")

# ==== Compute Bray–Curtis Dissimilarity ====
bray_dist <- vegdist(expr, method = "bray")

# ==== Run NMDS Ordination ====
set.seed(42)
nmds <- metaMDS(expr, k = 3, trymax = 100, try = 10, distance = "bray", no.share = TRUE, parallel = detectCores() - 1) 

NMDS_result <- nmds

NMDS_result$stress  # take the value of the stress result here. if below 0.05 = excellent. if >0.2 then bad if  between 0.05 to 0.1 = good fit of the NMDS model.
# Make a stress plot for the NMDS.The stressplot evaluates how well the ordination represented the complexity in the data.They help to show you how closely the ordination (y-axis) represents the dissimilarities calculated (x-axis). The points around the red stair steps are the communities, and the distance from the line represents the “stress”, or how they are pulled from their original position to be represented in their ordination.
stressplot(nmds) # visualize how well the model fits here.

tiff("Species_NMDS_Stress_plot.tiff", width = 4, height = 4, units = "in", res = 600)
# Customize the plot
par(
  font.axis = 2,   # Bold axis text
  font.lab = 2,    # Bold axis labels
  lwd = 2          # Thickness of the axis lines and border
)

# Re-plot with customizations
stressplot(NMDS_result)

# ==== Save Stress Plot ====
tiff("Species_NMDS_Stress_plot.tiff", width = 4, height = 4, units = "in", res = 600)
par(
  font.axis = 2,
  font.lab = 2,
  lwd = 2
)
stressplot(NMDS_result)
box(lwd = 2)
dev.off()

# Extract NMDS site scores
nmds_points <- as.data.frame(scores(nmds, display = "sites"))
nmds_points$Group <- meta$Group

# Define group color palette
group_colors <- c("Adenoma" = "#984ea3", "Cancer" = "#4daf4a", "Healthy" = "#ffcc00")

# Prepare NMDS points (already aligned)
nmds_points <- as.data.frame(scores(nmds, display = "sites"))
nmds_points$Group <- meta$Group

# Plot
tiff("Species_NMDS_plot_dispersion1.tiff", width = 6, height = 5, units = "in", res = 600)

ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Group, fill = Group)) +
  geom_point(size = 2.5, alpha = 0.85, stroke = 0.4, shape = 21, color = "black") +
  stat_ellipse(geom = "polygon", alpha = 0.05, aes(fill = Group), show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    axis.line = element_line(size = 1.2, color = "black"),
    axis.ticks = element_line(size = 0.8, color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    legend.position = c(0.12, 0.14),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = "white", size = 0.3)
  )

dev.off()

####################################################################
# ==== Run ANOSIM ====
anosim_result <- anosim(bray_dist, grouping = meta$Group, permutations = 999)

# Print ANOSIM results
cat("\n✅ ANOSIM Results:\n")
print(anosim_result)

# Save ANOSIM summary to file
anosim_summary <- data.frame(
  Statistic_R = anosim_result$statistic,
  P_Value = anosim_result$signif,
  Permutations = anosim_result$permutations
)
write.csv(anosim_summary, "ANOSIM_Results.csv", row.names = FALSE)

cat("✅ ANOSIM analysis completed and results saved to 'ANOSIM_Results.csv'\n")
######################PLot ANOSIM ######################################################
tiff("ANOSIM_species.tiff", width = 9, height = 5, units = "in", res = 600)
plot(anosim_result, col = c("grey", "#984ea3", "#4daf4a",  "#ffcc00"),
     ylab = "Dissimilarity Rank Value", xlab = "",
     cex.lab = 1.5, cex.axis = 1.2, # Increase axis label and text size
     lwd = 2) # Make line width thicker
box(lwd = 2.5) # Add a thick border around the plot
dev.off()
#################################################################################################
# ==== PERMDISP: Test homogeneity of group dispersions ====
bd <- betadisper(bray_dist, meta$Group)               # Calculate distances to group centroids
permdisp_result <- permutest(bd, permutations = 999)  # Permutation test

# Print result
cat("\n✅ PERMDISP Results (dispersion assumption for PERMANOVA):\n")
print(permdisp_result)

# Save test summary
write.csv(as.data.frame(permdisp_result$tab), "PERMDISP_Results.csv")
########################## PLot Permdisp #####################################################################
# ==== Extract distances and group labels ====
disp_df <- data.frame(
  Sample_ID = names(bd$distances),
  Distance = bd$distances
)

# Merge with metadata to get group labels
disp_df <- merge(disp_df, meta[, c("Sample_ID", "Group")], by = "Sample_ID")

# Ensure Group is treated as a factor (not numeric or character)
disp_df$Group <- factor(disp_df$Group, levels = c("Healthy", "Adenoma", "Cancer"))

# ==== Create Violin Plot ====
tiff("PERMDISP_ViolinPlot.tiff", width = 6, height = 5, units = "in", res = 600)

ggplot(disp_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "black") +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 3, fill = "white", color = "black") +
  scale_fill_manual(values = c("Adenoma" = "#984ea3", "Cancer" = "#4daf4a", "Healthy" = "#ffcc00")) +
  theme_classic(base_size = 16) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 1.2)
  ) +
  labs(y = "Distance to Group Centroid", title = "")

dev.off()
##############################################################################################
# ==== PERMANOVA: Group differences in Bray–Curtis space ====
permanova_result <- adonis2(bray_dist ~ Group, data = meta, permutations = 999, parallel = detectCores() - 1)
print(permanova_result)
# ==== Save PERMANOVA Output ====
write.csv(as.data.frame(permanova_result), "PERMANOVA_Results.csv")

# If not installed:
#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)

# Run pairwise PERMANOVA on Bray-Curtis distances
pairwise_results <- pairwise.adonis2(
  bray_dist ~ Group,
  data = meta,
  permutations = 999,
  p.adjust.m = "bonferroni"  # or "fdr"
)

# View and save results
print(pairwise_results)

write.csv(pairwise_results, "Pairwise_PERMANOVA_Results.csv", row.names = FALSE)

cat("✅ Pairwise PERMANOVA completed and results saved to 'Pairwise_PERMANOVA_Results.csv'\n")

#############################################
# === Step 1: Extract R² values from pairwise PERMANOVA result ===
# Skip the first element (metadata/parent call)
pairwise_r2 <- sapply(pairwise_results[-1], function(x) {
  x["Model", "R2"]
})

# Extract pair names
pair_names <- names(pairwise_r2)

# Create a square matrix (3x3) and fill upper or lower triangle
group_levels <- unique(meta$Group)
r2_matrix <- matrix(NA, nrow = length(group_levels), ncol = length(group_levels))
rownames(r2_matrix) <- colnames(r2_matrix) <- group_levels

# Populate matrix using pairwise comparisons
for (name in pair_names) {
  parts <- unlist(strsplit(name, "_vs_"))
  if (length(parts) == 2 && all(parts %in% group_levels)) {
    r2_matrix[parts[1], parts[2]] <- pairwise_r2[[name]]
    r2_matrix[parts[2], parts[1]] <- pairwise_r2[[name]]  # symmetric
  }
}

# === Step 2: Convert to data frame for heatmap ===
library(reshape2)
r2_df <- melt(r2_matrix, na.rm = FALSE)

# === Step 3: Plot Heatmap ===
tiff("PERMANOVA_R2_Heatmap_Auto.tiff", width = 6, height = 4, units = "in", res = 600)
ggplot(r2_df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = ifelse(is.na(value), "", round(value, 3))), size = 5) +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  labs(x = "", y = "", fill = expression(R^2), title = "Pairwise PERMANOVA (R²)") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()

cat("✅ Automatic PERMANOVA R² heatmap saved as 'PERMANOVA_R2_Heatmap_Auto.tiff'\n")

###################CCA ANALYSIS ##############################################################
# ==== Load CCA library if not already loaded ====
library(vegan)

# ==== Prepare input matrices ====
# Y = species data (samples × taxa) — already your expr matrix
# X = environmental/clinical variables — from meta

# ==== Ensure Covariates Are Categorical/Numeric ====
meta$Group <- factor(meta$Group)
meta$Sex <- factor(meta$Sex)
meta$Age <- as.numeric(meta$Age)
meta$BMI <- as.numeric(meta$BMI)

# ==== Run Adjusted CCA Model (Group + Age + BMI + Sex) ==== Evaluates the effect of the Group + covariates on the microbiome structure.
cca_model_adj <- cca(expr ~ Group + Age + BMI + Sex, data = meta)

# ==== Global and Term-level Tests ====
anova_global_adj <- anova(cca_model_adj, permutations = 999)
anova_terms_adj <- anova(cca_model_adj, by = "term", permutations = 999)

# ==== Save Adjusted Model Outputs ====
sink("CCA_Adjusted_Summary.txt")
print(summary(cca_model_adj))
sink()

write.csv(as.data.frame(anova_global_adj), "CCA_Adjusted_Global_Significance.csv")
write.csv(as.data.frame(anova_terms_adj), "CCA_Adjusted_Term_Significance.csv")

# ==== Run Partial CCA Model (Group | Age + BMI + Sex) ==== Isolates the effect of the Group on the microbiome structure independent of the covariates.
cca_model_partial <- cca(expr ~ Group + Condition(Age + BMI + Sex), data = meta)

# ==== Global and Term-level Tests for Partial Model ====
anova_global_partial <- anova(cca_model_partial, permutations = 999)
anova_terms_partial <- anova(cca_model_partial, by = "term", permutations = 999)

# ==== Save Partial Model Outputs ====
sink("CCA_Partial_Summary.txt")
print(summary(cca_model_partial))
sink()

write.csv(as.data.frame(anova_global_partial), "CCA_Partial_Global_Significance.csv")
write.csv(as.data.frame(anova_terms_partial), "CCA_Partial_Term_Significance.csv")

cat("✅ CCA analysis (adjusted and partial) completed. Results saved.\n")
##########################################################################################################################

# === Multivariable CCA Axis Labels ===
eig_adj <- summary(cca_model_adj)$concont$importance
cca1_lab <- paste0("CCA1 (", round(eig_adj["Proportion Explained", 1] * 100, 1), "%)")
cca2_lab <- paste0("CCA2 (", round(eig_adj["Proportion Explained", 2] * 100, 1), "%)")

tiff("CCA_Multivariable_Group_Age_BMI_Sex.tiff", width = 6, height = 5, units = "in", res = 600)

plot(cca_model_adj,
     display = c("sites", "bp"),
     xlab = cca1_lab,
     ylab = cca2_lab,
     col = "gray40",
     cex = 0.9)
box(lwd = 2.8)

dev.off()
#########################################################################################################

# === Partial CCA Axis Labels ===
eig_partial <- summary(cca_model_partial)$concont$importance
cca1_lab_p <- paste0("CCA1 (", round(eig_partial["Proportion Explained", 1] * 100, 1), "%)")
cca2_lab_p <- paste0("CCA2 (", round(eig_partial["Proportion Explained", 2] * 100, 1), "%)")

tiff("CCA_Partial_Group_only_Adjusted.tiff", width = 6.2, height = 5, units = "in", res = 600)

plot(cca_model_partial,
     display = c("sites", "bp"),
     xlab = cca1_lab_p,
     ylab = cca2_lab_p,
     col = "gray40",
     cex = 0.9)
box(lwd = 2.8)

dev.off()

