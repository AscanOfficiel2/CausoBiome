# === Load Packages ===
library(vegan)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(effsize)
library(broom)
library(patchwork)
library(ggeffects)

# === Load Data ===
metadata <- read.csv("Metadata_Aligned_to_CountMatrix.csv", row.names = 1)
gene_matrix <- read.csv("combined_ARG_VFDB_final_DATA.csv", row.names = 1)

# Transpose gene matrix to samples × genes
gene_matrix_t <- t(gene_matrix)
common_samples <- intersect(rownames(metadata), rownames(gene_matrix_t))
metadata <- metadata[common_samples, ]
gene_matrix_t <- gene_matrix_t[common_samples, ]

# === Filter Metadata with Complete Covariates ===
metadata <- metadata[!is.na(metadata$Sex) & !is.na(metadata$Age) & !is.na(metadata$BMI), ]

# === Split into ARG and VF Matrices ===
arg_cols <- grep("^ARG_", colnames(gene_matrix_t), value = TRUE)
vf_cols  <- grep("^VF_", colnames(gene_matrix_t), value = TRUE)
arg_matrix <- gene_matrix_t[rownames(metadata), arg_cols]
vf_matrix  <- gene_matrix_t[rownames(metadata), vf_cols]

# === Alpha Diversity Computation ===
arg_alpha <- data.frame(
  Sample   = rownames(arg_matrix),
  Richness = specnumber(arg_matrix),
  Shannon  = diversity(arg_matrix, index = "shannon"),
  Simpson  = diversity(arg_matrix, index = "simpson"),
  Group    = metadata$Group,
  Age      = metadata$Age,
  BMI      = metadata$BMI,
  Sex      = factor(metadata$Sex)
)

vf_alpha <- data.frame(
  Sample   = rownames(vf_matrix),
  Richness = specnumber(vf_matrix),
  Shannon  = diversity(vf_matrix, index = "shannon"),
  Simpson  = diversity(vf_matrix, index = "simpson"),
  Group    = metadata$Group,
  Age      = metadata$Age,
  BMI      = metadata$BMI,
  Sex      = factor(metadata$Sex)
)

# === Statistical Models ===
# ARG
lm_arg_rich <- lm(Richness ~ Group + Age + BMI + Sex, data = arg_alpha)
lm_arg_shan <- lm(Shannon  ~ Group + Age + BMI + Sex, data = arg_alpha)
lm_arg_simp <- lm(Simpson  ~ Group + Age + BMI + Sex, data = arg_alpha)

# VF
lm_vf_rich <- lm(Richness ~ Group + Age + BMI + Sex, data = vf_alpha)
lm_vf_shan <- lm(Shannon  ~ Group + Age + BMI + Sex, data = vf_alpha)
lm_vf_simp <- lm(Simpson  ~ Group + Age + BMI + Sex, data = vf_alpha)

# === Save LM Summaries ===
sink("alpha_diversity_linear_models.txt")
cat("=== Linear Models: ARG Alpha Diversity ===\n\n")
cat("ARG Richness:\n"); print(summary(lm_arg_rich)); cat("\n")
cat("ARG Shannon:\n"); print(summary(lm_arg_shan)); cat("\n")
cat("ARG Simpson:\n"); print(summary(lm_arg_simp)); cat("\n")

cat("=== Linear Models: VF Alpha Diversity ===\n\n")
cat("VF Richness:\n"); print(summary(lm_vf_rich)); cat("\n")
cat("VF Shannon:\n"); print(summary(lm_vf_shan)); cat("\n")
cat("VF Simpson:\n"); print(summary(lm_vf_simp)); cat("\n")
sink()

# === Plotting Functions ===
plot_box_jitter <- function(df, metric, label, colors) {
  ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, size = 1.1, alpha = 0.9) +
    geom_jitter(width = 0.2, shape = 21, size = 2, stroke = 0.4, color = "black", fill = "black") +
    scale_fill_manual(values = colors) +
    labs(title = label, y = metric, x = NULL) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          legend.position = "none")
}

plot_model_pred <- function(model, title_label, x_var = "Group") {
  pred <- ggpredict(model, terms = x_var)
  ggplot(pred, aes(x = x, y = predicted)) +
    geom_line(group = 1, color = "black") +
    geom_point(size = 3, color = "#0072B2") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
    labs(title = title_label, x = x_var, y = "Predicted Value") +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

# === Define Group Colors ===
group_colors <- c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")

# ARG panels
arg_box_rich <- plot_box_jitter(arg_alpha, "Richness", "ARG Richness", group_colors)
arg_box_shan <- plot_box_jitter(arg_alpha, "Shannon", "ARG Shannon", group_colors)
arg_box_simp <- plot_box_jitter(arg_alpha, "Simpson", "ARG Simpson", group_colors)

arg_pred_rich <- plot_model_pred(lm_arg_rich, "Linear Model Richness")
arg_pred_shan <- plot_model_pred(lm_arg_shan, "Linear Model Shannon")
arg_pred_simp <- plot_model_pred(lm_arg_simp, "Linear Model Simpson")

arg_panel <- (arg_box_rich | arg_pred_rich) /
  (arg_box_shan | arg_pred_shan) /
  (arg_box_simp | arg_pred_simp)

ggsave("Alpha_ARG_Box_vs_Model.png", arg_panel, width = 10, height = 10, dpi = 600)

# VF panels
vf_box_rich <- plot_box_jitter(vf_alpha, "Richness", "VF Richness", group_colors)
vf_box_shan <- plot_box_jitter(vf_alpha, "Shannon", "VF Shannon", group_colors)
vf_box_simp <- plot_box_jitter(vf_alpha, "Simpson", "VF Simpson", group_colors)

vf_pred_rich <- plot_model_pred(lm_vf_rich, "Linear Model Richness")
vf_pred_shan <- plot_model_pred(lm_vf_shan, "Linear Model Shannon")
vf_pred_simp <- plot_model_pred(lm_vf_simp, "Linear Model Simpson")

vf_panel <- (vf_box_rich | vf_pred_rich) /
  (vf_box_shan | vf_pred_shan) /
  (vf_box_simp | vf_pred_simp)

ggsave("Alpha_VF_Box_vs_Model.png", vf_panel, width = 10, height = 10, dpi = 600)

#################################################################################################

# === Separate ARGs and VFs ===
arg_cols <- grep("^ARG_", colnames(gene_matrix_t), value = TRUE)
vf_cols  <- grep("^VF_", colnames(gene_matrix_t), value = TRUE)
arg_matrix <- gene_matrix_t[, arg_cols]
vf_matrix  <- gene_matrix_t[, vf_cols]

# === Bray-Curtis dissimilarity ===
bray_arg <- vegdist(arg_matrix, method = "bray")
bray_vf  <- vegdist(vf_matrix, method = "bray")

# === PERMANOVA controlling only for Age + BMI + Sex ===
adonis_arg_adj <- adonis2(bray_arg ~ Group + Age + BMI + Sex, data = metadata, permutations = 999)
adonis_vf_adj  <- adonis2(bray_vf  ~ Group + Age + BMI + Sex, data = metadata, permutations = 999)

write.csv(as.data.frame(adonis_arg_adj), "adonis_ARG_adjusted.csv")
write.csv(as.data.frame(adonis_vf_adj),  "adonis_VF_adjusted.csv")

################################################################################################################

# === Beta-dispersion ===
disp_arg <- betadisper(bray_arg, metadata$Group)
disp_vf  <- betadisper(bray_vf, metadata$Group)
anova_disp_arg <- anova(disp_arg)
anova_disp_vf  <- anova(disp_vf)
write.csv(as.data.frame(anova_disp_arg), "dispersion_ARG_results.csv")
write.csv(as.data.frame(anova_disp_vf),  "dispersion_VF_results.csv")

# === Prepare data for plotting/stats ===
arg_df <- data.frame(Sample = names(disp_arg$distances),
                     Distance = disp_arg$distances,
                     Group = metadata[names(disp_arg$distances), "Group"])
vf_df <- data.frame(Sample = names(disp_vf$distances),
                    Distance = disp_vf$distances,
                    Group = metadata[names(disp_vf$distances), "Group"])
arg_df <- na.omit(arg_df)
vf_df  <- na.omit(vf_df)

# === Define group colors ===
group_colors <- c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")

# === Kruskal-Wallis and Wilcoxon tests ===
kruskal_arg <- kruskal.test(Distance ~ Group, data = arg_df)
kruskal_vf  <- kruskal.test(Distance ~ Group, data = vf_df)
pairwise_arg <- pairwise.wilcox.test(arg_df$Distance, arg_df$Group, p.adjust.method = "BH")
pairwise_vf  <- pairwise.wilcox.test(vf_df$Distance, vf_df$Group, p.adjust.method = "BH")

# === Cliff's Delta Effect Sizes ===
cliff_arg_ha <- cliff.delta(Distance ~ Group, data = filter(arg_df, Group %in% c("Healthy", "Adenoma")))
cliff_arg_hc <- cliff.delta(Distance ~ Group, data = filter(arg_df, Group %in% c("Healthy", "Cancer")))
cliff_arg_ac <- cliff.delta(Distance ~ Group, data = filter(arg_df, Group %in% c("Adenoma", "Cancer")))
cliff_vf_ha  <- cliff.delta(Distance ~ Group, data = filter(vf_df, Group %in% c("Healthy", "Adenoma")))
cliff_vf_hc  <- cliff.delta(Distance ~ Group, data = filter(vf_df, Group %in% c("Healthy", "Cancer")))
cliff_vf_ac  <- cliff.delta(Distance ~ Group, data = filter(vf_df, Group %in% c("Adenoma", "Cancer")))

# === Save all stats to TXT ===
sink("kruskal_pairwise_statistics.txt")
cat("=== ARG Kruskal-Wallis ===\n"); print(kruskal_arg)
cat("\n=== ARG Pairwise Wilcoxon ===\n"); print(pairwise_arg)
cat("\n=== VF Kruskal-Wallis ===\n"); print(kruskal_vf)
cat("\n=== VF Pairwise Wilcoxon ===\n"); print(pairwise_vf)
cat("\n=== Cliff's Delta Effect Sizes (ARGs) ===\n")
print(cliff_arg_ha); print(cliff_arg_hc); print(cliff_arg_ac)
cat("\n=== Cliff's Delta Effect Sizes (VFs) ===\n")
print(cliff_vf_ha); print(cliff_vf_hc); print(cliff_vf_ac)
sink()

# === Helper functions ===
safe_extract <- function(mat, g1, g2) {
  if (g1 %in% rownames(mat) && g2 %in% colnames(mat)) return(mat[g1, g2])
  if (g2 %in% rownames(mat) && g1 %in% colnames(mat)) return(mat[g2, g1])
  return(NA)
}
interpret_cliff <- function(d) {
  abs_d <- abs(d)
  if (abs_d < 0.147) return("negligible")
  else if (abs_d < 0.33) return("small")
  else if (abs_d < 0.474) return("medium")
  else return("large")
}

# === Create summary table ===
summary_table <- data.frame(
  Comparison = c("Healthy vs Adenoma", "Healthy vs Cancer", "Adenoma vs Cancer"),
  ARG_Wilcox_p = c(
    safe_extract(pairwise_arg$p.value, "Healthy", "Adenoma"),
    safe_extract(pairwise_arg$p.value, "Healthy", "Cancer"),
    safe_extract(pairwise_arg$p.value, "Adenoma", "Cancer")
  ),
  ARG_Cliffs_Delta = c(cliff_arg_ha$estimate, cliff_arg_hc$estimate, cliff_arg_ac$estimate),
  ARG_Effect_Size = c(
    interpret_cliff(cliff_arg_ha$estimate),
    interpret_cliff(cliff_arg_hc$estimate),
    interpret_cliff(cliff_arg_ac$estimate)
  ),
  VF_Wilcox_p = c(
    safe_extract(pairwise_vf$p.value, "Healthy", "Adenoma"),
    safe_extract(pairwise_vf$p.value, "Healthy", "Cancer"),
    safe_extract(pairwise_vf$p.value, "Adenoma", "Cancer")
  ),
  VF_Cliffs_Delta = c(cliff_vf_ha$estimate, cliff_vf_hc$estimate, cliff_vf_ac$estimate),
  VF_Effect_Size = c(
    interpret_cliff(cliff_vf_ha$estimate),
    interpret_cliff(cliff_vf_hc$estimate),
    interpret_cliff(cliff_vf_ac$estimate)
  )
)
summary_table <- summary_table %>% mutate(across(where(is.numeric), round, 4))
write.csv(summary_table, "pairwise_summary_table.csv", row.names = FALSE)

# === Boxplots ===
arg_plot <- ggplot(arg_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  labs(title = "ARG Beta-dispersion", y = "Distance to Centroid") +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(size = 1.5),
    axis.line = element_line(size = 1),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

vf_plot <- ggplot(vf_df, aes(x = Group, y = Distance, fill = Group)) +
  geom_boxplot(size = 1.2, outlier.shape = NA) +
  scale_fill_manual(values = group_colors) +
  labs(title = "VF Beta-dispersion", y = "Distance to Centroid") +
  theme_bw(base_size = 16) +
  theme(
    axis.title.x = element_blank(),
    panel.border = element_rect(size = 1.5),
    axis.line = element_line(size = 1),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("ARG_Beta_dispersion.png", arg_plot, width = 5, height = 5, dpi = 600)
ggsave("VF_Beta_dispersion.png", vf_plot, width = 5, height = 5, dpi = 600)

# === NMDS Plots ===
set.seed(123)
nmds_arg <- metaMDS(arg_matrix, distance = "bray", k = 2, trymax = 100)
nmds_arg_df <- as.data.frame(nmds_arg$points)
nmds_arg_df$Group <- metadata[rownames(nmds_arg_df), "Group"]

arg_nmds_plot <- ggplot(nmds_arg_df, aes(MDS1, MDS2, color = Group, fill = Group)) +
  stat_ellipse(type = "norm", size = 1.2, linetype = "solid") +
  geom_point(shape = 21, size = 3.5, stroke = 1.2, color = "black") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(title = "NMDS of ARG Composition", x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )
ggsave("NMDS_ARG_plot.png", arg_nmds_plot, width = 6, height = 5, dpi = 600)

set.seed(123)
nmds_vf <- metaMDS(vf_matrix, distance = "bray", k = 2, trymax = 100)
nmds_vf_df <- as.data.frame(nmds_vf$points)
nmds_vf_df$Group <- metadata[rownames(nmds_vf_df), "Group"]

vf_nmds_plot <- ggplot(nmds_vf_df, aes(MDS1, MDS2, color = Group, fill = Group)) +
  stat_ellipse(type = "norm", size = 1.2, linetype = "solid") +
  geom_point(shape = 21, size = 3.5, stroke = 1.2, color = "black") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(title = "NMDS of VF Composition", x = "NMDS1", y = "NMDS2") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(size = 1.2),
    axis.ticks = element_line(size = 1.2),
    panel.border = element_rect(colour = "black", fill = NA, size = 1.5),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )
ggsave("NMDS_VF_plot.png", vf_nmds_plot, width = 6, height = 5, dpi = 600)

# === Shepard Plot for ARGs ===
png("Stressplot_ARG_Better.png", width = 1800, height = 1600, res = 300)

par(mar = c(5, 5, 4, 2), lwd = 2, cex.main = 2, cex.lab = 1.6, cex.axis = 1.3)
stressplot(nmds_arg, main = "Stress Plot: ARG")

dev.off()

png("Stressplot_VF_Better.png", width = 1800, height = 1600, res = 300)

par(mar = c(5, 5, 4, 2), lwd = 2, cex.main = 2, cex.lab = 1.6, cex.axis = 1.3)
stressplot(nmds_vf, main = "Stress Plot: VF")

dev.off()

png("Stressplot_ARG_VF_SideBySide.png", width = 3600, height = 1600, res = 300)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2), lwd = 2, cex.main = 2, cex.lab = 1.6, cex.axis = 1.3)

stressplot(nmds_arg, main = "Stress Plot: ARG")
stressplot(nmds_vf,  main = "Stress Plot: VF")

dev.off()


######################################################################################
# === Procrustes Analysis ===
proc_fit <- protest(nmds_arg, nmds_vf, permutations = 999)
print(proc_fit)  # Shows r and p-value

# Save stats
proc_stats <- data.frame(
  Procrustes_r = round(proc_fit$t0, 3),
  p_value = proc_fit$signif
)
write.csv(proc_stats, "procrustes_statistics.csv", row.names = FALSE)

#######################################################################################
# === Variance partitioning (capscale) for ARGs ===
arg_rda <- capscale(arg_matrix ~ Group + Age + BMI + Sex, distance = "bray", data = metadata)
arg_rda_results <- anova(arg_rda, by = "margin", permutations = 999)
write.csv(as.data.frame(arg_rda_results), "capscale_ARG_results.csv")

# === Variance partitioning (capscale) for VFs ===
vf_rda <- capscale(vf_matrix ~ Group + Age + BMI + Sex, distance = "bray", data = metadata)
vf_rda_results <- anova(vf_rda, by = "margin", permutations = 999)
write.csv(as.data.frame(vf_rda_results), "capscale_VF_results.csv")

#######################################################################################################
# Convert abundance matrices to presence/absence
arg_pa <- decostand(arg_matrix, method = "pa")
vf_pa  <- decostand(vf_matrix, method = "pa")

# Run null model test using checkerboard score (C-score)
set.seed(123)
null_arg <- oecosimu(arg_pa, nestedchecker, method = "swap", nsimul = 999)
null_vf  <- oecosimu(vf_pa,  nestedchecker, method = "swap", nsimul = 999)

# Print to console
cat("=== Oecosimu: ARGs ===\n")
print(null_arg)

cat("\n=== Oecosimu: VFs ===\n")
print(null_vf)

# Save to files
sink("oecosimu_ARG_results.txt")
print(null_arg)
sink()

sink("oecosimu_VF_results.txt")
print(null_vf)
sink()
