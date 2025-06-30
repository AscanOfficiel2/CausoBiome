######################### SETUP #########################
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(effsize)
library(patchwork)
library(ggeffects)
library(broom)
library(grid)
library(gridGraphics)

# === Load data ===
hellinger <- read.csv("hellinger_transformed_combat.csv", row.names = 1)
metadata  <- read.csv("Aligned_metadata_Taxonomy.csv", row.names = 1)

# Transpose to sample × feature
hellinger_t <- t(hellinger)

# Align samples
common_samples <- intersect(rownames(hellinger_t), rownames(metadata))
hellinger_t <- hellinger_t[common_samples, ]
metadata    <- metadata[common_samples, ]

# Ensure factors
metadata$Group <- factor(metadata$Group, levels = c("Healthy", "Adenoma", "Cancer"))
metadata$Sex   <- factor(metadata$Sex)

#########################################################
# === Alpha Diversity ===
alpha_div <- data.frame(
  Sample   = rownames(hellinger_t),
  Richness = specnumber(hellinger_t),
  Shannon  = diversity(hellinger_t, index = "shannon"),
  Simpson  = diversity(hellinger_t, index = "simpson"),
  Group    = metadata$Group,
  Age      = metadata$Age,
  BMI      = metadata$BMI,
  Sex      = metadata$Sex
)

# Kruskal-Wallis + Wilcoxon
kw_rich <- kruskal.test(Richness ~ Group, data = alpha_div)
kw_shan <- kruskal.test(Shannon  ~ Group, data = alpha_div)
kw_simp <- kruskal.test(Simpson  ~ Group, data = alpha_div)
pair_rich <- pairwise.wilcox.test(alpha_div$Richness, alpha_div$Group, p.adjust.method = "BH")
pair_shan <- pairwise.wilcox.test(alpha_div$Shannon,  alpha_div$Group, p.adjust.method = "BH")
pair_simp <- pairwise.wilcox.test(alpha_div$Simpson,  alpha_div$Group, p.adjust.method = "BH")

# Linear models
lm_rich <- lm(Richness ~ Group + Age + BMI + Sex, data = alpha_div)
lm_shan <- lm(Shannon  ~ Group + Age + BMI + Sex, data = alpha_div)
lm_simp <- lm(Simpson  ~ Group + Age + BMI + Sex, data = alpha_div)

# Save
sink("alpha_diversity_statistics_ecology.txt")
cat("=== Kruskal-Wallis Tests ===\n"); print(kw_rich); print(kw_shan); print(kw_simp)
cat("\n=== Pairwise Wilcoxon Tests ===\n"); print(pair_rich); print(pair_shan); print(pair_simp)
cat("\n=== Adjusted Linear Models ===\n")
cat("\nRichness:\n"); print(summary(lm_rich))
cat("\nShannon:\n"); print(summary(lm_shan))
cat("\nSimpson:\n"); print(summary(lm_simp))
sink()

##############################################################################
plot_alpha_box_nature <- function(df, metric, title_label) {
  ggplot(df, aes(x = Group, y = .data[[metric]])) +  # <-- Removed fill = Group
    geom_boxplot(aes(fill = Group), outlier.shape = NA, alpha = 0.6, color = "black", width = 0.55, size = 1) +
    geom_jitter(shape = 21, size = 2.0, stroke = 0.2, width = 0.15,alpha = 0.7,  color = "black", fill = "black") +  
    scale_fill_manual(values = c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")) +
    labs(title = title_label, x = NULL, y = metric) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      axis.title.y = element_text(face = "bold", size = 14),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "none",
      axis.ticks.length = unit(0.3, "cm"),
      axis.line = element_line(size = 1),
      axis.ticks = element_line(size = 0.8)
    )
}


alpha_shan_plot <- plot_alpha_box_nature(alpha_div, "Shannon", "")
alpha_simp_plot <- plot_alpha_box_nature(alpha_div, "Simpson", "")
alpha_panel <- alpha_shan_plot + alpha_simp_plot + plot_layout(ncol = 2)
ggsave("Alpha_Diversity_Panel.png", alpha_panel, width = 10, height = 4.5, dpi = 600)

#########################################################
# === Beta Diversity ===
bray_dist <- vegdist(hellinger_t, method = "bray")

# Unadjusted PERMANOVA
adonis_unadj <- adonis2(bray_dist ~ Group, data = metadata, permutations = 999)
write.csv(as.data.frame(adonis_unadj), "adonis_unadjusted.csv")

# Adjusted PERMANOVA
adonis_adj <- adonis2(bray_dist ~ Group + Age + BMI + Sex, data = metadata, permutations = 999)
write.csv(as.data.frame(adonis_adj), "adonis_adjusted.csv")

##############################################################################################################################

# ==== Run ANOSIM ====
anosim_result <- anosim(bray_dist, grouping = metadata$Group, permutations = 999)

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

######################PLot ANOSIM ######################################################
tiff("ANOSIM_species.tiff", width = 9, height = 5, units = "in", res = 600)
plot(anosim_result, col = c("grey", "#984ea3", "#4daf4a",  "#ffcc00"),
     ylab = "Dissimilarity Rank Value", xlab = "",
     cex.lab = 1.5, cex.axis = 1.2, # Increase axis label and text size
     lwd = 2) # Make line width thicker
box(lwd = 2.5) # Add a thick border around the plot
dev.off()

###########################################################################################################
# Beta-dispersion
disp <- betadisper(bray_dist, metadata$Group)
anova_disp <- anova(disp)
write.csv(as.data.frame(anova_disp), "dispersion_results.csv")

# Distance to centroid
beta_df <- data.frame(Sample = names(disp$distances),
                      Distance = disp$distances,
                      Group = metadata[names(disp$distances), "Group"])
beta_df <- na.omit(beta_df)

# Kruskal-Wallis and pairwise
kruskal_beta <- kruskal.test(Distance ~ Group, data = beta_df)
pair_beta <- pairwise.wilcox.test(beta_df$Distance, beta_df$Group, p.adjust.method = "BH")

# Cliff's delta
safe_cliff <- function(data, g1, g2) {
  sub_df <- filter(data, Group %in% c(g1, g2))
  sub_df$Group <- droplevels(sub_df$Group)
  if (nlevels(sub_df$Group) != 2) return(NA)
  tryCatch({
    cliff.delta(Distance ~ Group, data = sub_df)
  }, error = function(e) NA)
}
cliff_ha <- safe_cliff(beta_df, "Healthy", "Adenoma")
cliff_hc <- safe_cliff(beta_df, "Healthy", "Cancer")
cliff_ac <- safe_cliff(beta_df, "Adenoma", "Cancer")

sink("beta_diversity_statistics.txt")
print(kruskal_beta)
print(pair_beta)
print(cliff_ha); print(cliff_hc); print(cliff_ac)
sink()

safe_extract <- function(mat, g1, g2) {
  if (g1 %in% rownames(mat) && g2 %in% colnames(mat)) return(mat[g1, g2])
  if (g2 %in% rownames(mat) && g1 %in% colnames(mat)) return(mat[g2, g1])
  return(NA)
}
interpret_cliff <- function(d) {
  abs_d <- abs(d)
  if (is.na(abs_d)) return("NA")
  if (abs_d < 0.147) return("negligible")
  else if (abs_d < 0.33) return("small")
  else if (abs_d < 0.474) return("medium")
  else return("large")
}

summary_table <- data.frame(
  Comparison = c("Healthy vs Adenoma", "Healthy vs Cancer", "Adenoma vs Cancer"),
  Wilcox_p = c(
    safe_extract(pair_beta$p.value, "Healthy", "Adenoma"),
    safe_extract(pair_beta$p.value, "Healthy", "Cancer"),
    safe_extract(pair_beta$p.value, "Adenoma", "Cancer")
  ),
  Cliff_Delta = c(cliff_ha$estimate, cliff_hc$estimate, cliff_ac$estimate),
  Effect_Size = c(
    interpret_cliff(cliff_ha$estimate),
    interpret_cliff(cliff_hc$estimate),
    interpret_cliff(cliff_ac$estimate)
  )
)
write.csv(summary_table, "beta_diversity_summary.csv", row.names = FALSE)

#########################################################
# === NMDS Ordination + Stress ===
set.seed(123)
nmds <- metaMDS(hellinger_t, distance = "bray", k = 2, trymax = 100)
stress_val <- round(nmds$stress, 3)

# Prepare NMDS plot data
nmds_df <- as.data.frame(nmds$points)
nmds_df$Group <- metadata[rownames(nmds_df), "Group"]

# === NMDS Plot (separate)
nmds_plot <- ggplot(nmds_df, aes(MDS1, MDS2, color = Group, fill = Group)) +
  stat_ellipse(type = "norm", size = 1.2, alpha = 0.25) +
  geom_point(shape = 21, size = 3.5, stroke = 1.2, color = "black") +
  scale_color_manual(values = c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")) +
  scale_fill_manual(values = c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")) +
  labs(
    title = paste0("NMDS Ordination (Stress = ", stress_val, ")"),
    x = "NMDS1", y = "NMDS2"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank()
  )

# Save NMDS plot
ggsave("NMDS_Only.png", nmds_plot, width = 6, height = 5, dpi = 600)


# === Shepard (Stress) Plot (base R)
png("Stressplot_Only.png", width = 1800, height = 1600, res = 300)
par(mar = c(5, 5, 4, 2), lwd = 2, cex.main = 2, cex.lab = 1.6, cex.axis = 1.3)
stressplot(nmds, main = "Stress Plot")
dev.off()

