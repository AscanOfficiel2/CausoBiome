################################################################################
##Script Created BY Suleiman Aminu and AbdulAziz Ascandari
### 25th September 2025
####Title: Batch Correction using Limma

# --- Setting seed for reproducibility ---
set.seed(43)


# --- Step 1: Load libraries ---
library(limma)
library(compositions)   # CLR transform
library(vegan)          # PERMANOVA
library(ggplot2)        # plotting
library(gridExtra)      # arrange plots
library(umap)           # UMAP
library(Rtsne)          # t-SNE


############################################################################
# --- Step 2: Load TPM matrix and metadata ---
tpm <- read.csv("CRC_species_matrix_filtered.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("CRC_Cleaned_MetaData.csv")

# Match metadata order to TPM columns
samples <- colnames(tpm)
meta <- meta[match(samples, meta$Sample_ID), ]

# --- Step 3: Add pseudocount and CLR transform ---
tpm_pseudo <- tpm + 1
tpm_clr <- t(apply(tpm_pseudo, 2, function(x) clr(x)))
write.csv(tpm_clr, "TPM_clr.csv", quote = FALSE)
cat("✅ CLR-transformed matrix saved to TPM_clr.csv\n")


# If samples are in rows instead of columns → transpose
if(all(rownames(tpm_clr) %in% meta$Sample_ID | rownames(tpm) %in% meta$Sample_ID)) {
  tpm_clr <- t(tpm_clr)
}


# Match metadata to samples
samples <- colnames(tpm_clr)
meta <- meta[match(samples, meta$Sample_ID), ]

# --- Filter out zero-variance taxa ---
tpm_clr_filtered <- tpm_clr[apply(tpm_clr, 1, var) > 0, ]
cat("✅ Filtered CLR matrix:", nrow(tpm_clr_filtered), "taxa retained\n")



################################################
####Diagnostics  ######################

# --- PCA ---
pca <- prcomp(t(tpm_clr_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Group <- meta$Group


p_before <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "PCA Before Batch Correction")

ggsave("PCA_Before_Correction.png", p_before, width = 6, height = 5, dpi = 600, bg = "white")

# --- PERMANOVA (for multiple metadata factors) ---
counts_t <- t(tpm_clr_filtered)  # samples as rows

adonis_results <- list(
  Group = adonis2(counts_t ~ meta$Group, method = "euclidean"),
  Project = adonis2(counts_t ~ meta$Project, method = "euclidean"),
  Country = adonis2(counts_t ~ meta$Country, method = "euclidean"),
  Center = adonis2(counts_t ~ meta$Center_Name, method = "euclidean")
)

adonis_df <- do.call(rbind, lapply(names(adonis_results), function(x) {
  df <- as.data.frame(adonis_results[[x]])
  df$Factor <- x
  return(df)
}))
write.csv(adonis_df, "PERMANOVA_before.csv", row.names = FALSE)

# --- UMAP ---
umap_before <- umap(counts_t)
umap_df <- as.data.frame(umap_before$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Group <- meta$Group

u_before <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "UMAP Before Batch Correction")

ggsave("UMAP_Before_Correction.png", u_before, width = 6, height = 5, dpi = 600, bg = "white")

# --- t-SNE ---
set.seed(42)

# Remove duplicate rows (samples)
counts_unique <- counts_t[!duplicated(counts_t), ]

tsne_before <- Rtsne(counts_unique, dims = 2, perplexity = 30)
tsne_df <- as.data.frame(tsne_before$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")

# Match back to metadata (remove duplicates there too)
meta_unique <- meta[!duplicated(counts_t), ]
tsne_df$Group <- meta_unique$Group

ts_before <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "t-SNE Before Batch Correction")

ggsave("tSNE_Before_Correction.png", ts_before, width = 6, height = 5, dpi = 600, bg = "white")

######################################################################################################

####### Batch Correction using Limma ################################


# --- Step 4a: Correction with Project + Center_Name only ---
design <- model.matrix(~ meta$Group)   # preserve biological signal
covariates <- model.matrix(~ meta$Project + meta$Center_Name)[, -1]

counts_corrected_v1 <- removeBatchEffect(
  as.matrix(tpm_clr_filtered),
  covariates = covariates,
  design = design
)

write.csv(counts_corrected_v1, "TPM_clr_batch_corrected_v1.csv", row.names = TRUE)
cat("✅ Version 1 correction done (Instrument + Center only).\n")


#########################################################################################################


# --- Function to run diagnostics after correction ---
run_diagnostics <- function(counts_corrected, meta, prefix) {
  # PCA
  pca <- prcomp(t(counts_corrected), scale. = TRUE)
  pca_df <- as.data.frame(pca$x)
  pca_df$Group <- meta$Group
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2) + theme_minimal() + labs(title = paste("PCA After Correction", prefix))
  ggsave(paste0("PCA_After_Correction_", prefix, ".png"), p, width = 6, height = 5, dpi = 600, bg = "white")
  
  # PERMANOVA
  counts_t <- t(counts_corrected)
  adonis_results <- list(
    Group = adonis2(counts_t ~ meta$Group, method = "euclidean"),
    Center = adonis2(counts_t ~ meta$Center_Name, method = "euclidean"),
    Project = adonis2(counts_t ~ meta$Project, method = "euclidean")
  )
  adonis_df <- do.call(rbind, lapply(names(adonis_results), function(x) {
    df <- as.data.frame(adonis_results[[x]])
    df$Factor <- x
    return(df)
  }))
  write.csv(adonis_df, paste0("PERMANOVA_after_", prefix, ".csv"), row.names = FALSE)
  
  # UMAP
  umap_res <- umap(counts_t)
  umap_df <- as.data.frame(umap_res$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Group <- meta$Group
  u <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
    geom_point(size = 2) + theme_minimal() + labs(title = paste("UMAP After Correction", prefix))
  ggsave(paste0("UMAP_After_Correction_", prefix, ".png"), u, width = 6, height = 5, dpi = 600, bg = "white")
  
  # t-SNE
  set.seed(42)
  counts_unique <- counts_t[!duplicated(counts_t), ]
  tsne_res <- Rtsne(counts_unique, dims = 2, perplexity = 30)
  tsne_df <- as.data.frame(tsne_res$Y)
  colnames(tsne_df) <- c("tSNE1", "tSNE2")
  meta_unique <- meta[!duplicated(counts_t), ]
  tsne_df$Group <- meta_unique$Group
  ts <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Group)) +
    geom_point(size = 2) + theme_minimal() + labs(title = paste("t-SNE After Correction", prefix))
  ggsave(paste0("tSNE_After_Correction_", prefix, ".png"), ts, width = 6, height = 5, dpi = 600, bg = "white")
}

# --- Run diagnostics for both versions ---
tpm_corr_v1 <- read.csv("TPM_clr_batch_corrected_v1.csv", row.names = 1, check.names = FALSE)


run_diagnostics(tpm_corr_v1, meta, "v1")  # Instrument + Center
run_diagnostics(tpm_corr_v2, meta, "v2")  # Instrument + Center + Project_ID



library(ggplot2)
library(dplyr)

# --- Load PERMANOVA results ---
before <- read.csv("PERMANOVA_before.csv")
after_v1 <- read.csv("PERMANOVA_after_v1.csv")
after_v2 <- read.csv("PERMANOVA_after_v2.csv")

# Add condition labels
before$Condition <- "Before"
after_v1$Condition <- "After_v1"


# Combine
permanova_all <- bind_rows(before, after_v1)

# --- Keep only actual test rows ---
permanova_all <- permanova_all %>%
  filter(!is.na(F), !is.na(Pr..F.))   # drop residual/total rows

# --- Select key factors ---
key_factors <- c("Group", "Instrument", "Center", "Project")
permanova_all <- permanova_all %>%
  filter(Factor %in% key_factors)

# --- Barplot ---
ggplot(permanova_all, aes(x = Factor, y = R2, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  labs(
    title = "Variance explained (PERMANOVA R²)",
    y = "R² (Variance Explained)",
    x = "Metadata Factor"
  ) +
  scale_fill_manual(values = c("Before" = "#d95f02", "After_v1" = "#1b9e77", "After_v2" = "#7570b3")) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    legend.title = element_blank(),
    legend.position = "top"
  )

ggsave("PERMANOVA_R2_comparison_fixed.png", width = 8, height = 6, dpi = 600, bg = "white")




#################################
##### Figures For Publications


#------  (1) BarPlot  ----------------------------------------------------

library(ggplot2)
library(dplyr)


p <- ggplot(permanova_all, aes(x = Factor, y = R2, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black") +  # add bar border
  geom_text(aes(label = round(R2, 2)), 
            position = position_dodge(width = 0.8), 
            vjust = -0.4, size = 3.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "",
    y = expression(R^2~"(Variance Explained)"),
    x = "Metadata Factor"
  ) +
  scale_fill_manual(values = c("Before" = "#d95f02", 
                               "After_v1" = "#1b9e77", 
                               "After_v2" = "#7570b3")) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 15, hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.major = element_blank(),   # remove internal lines
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # add border
    axis.line = element_blank()
  )

# Save high-res TIFF
tiff("Barplot_variance_explained.tiff", width = 5, height = 5, units = "in", res = 600)
print(p)
dev.off()





#######################################################################

# --- (2) t-SNE Before and After Correction  ---

library(ggplot2)
library(dplyr)
library(shadowtext)

# Define group colors (consistent across plots)
group_colors <- c(
  "Healthy" = "#ffcc00",            
  "Adenoma" = "#984ea3", 
  "Cancer" = "#4daf4a"
)



###################
# Before Correction (tsne_df from BEFORE pipeline)
###################

tsne_plot_before <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Group)) +
  geom_point(size = 2, shape = 16, alpha = 1.6, stroke = 0) +
  scale_color_manual(values = group_colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

ggsave("tSNE_Before_Batch_Correction.tiff",
       plot = tsne_plot_before,
       width = 4, height = 4, dpi = 600, bg = "white")

###################
# After Correction v1 (tsne_df_after from v1 correction pipeline)
# ============================================
# STEP 1. Load libraries
# ============================================
library(ggplot2)
library(dplyr)
library(Rtsne)
library(shadowtext)

set.seed(42)

# ============================================
# STEP 2. Load your corrected matrix
# ============================================
tpm_corr_v1 <- read.csv("TPM_clr_batch_corrected_v1.csv", row.names = 1, check.names = FALSE)
counts_t_after <- t(tpm_corr_v1)
counts_unique_after <- counts_t_after[!duplicated(counts_t_after), ]

# Run t-SNE
tsne_after <- Rtsne(counts_unique_after, dims = 2, perplexity = 30, verbose = TRUE)
tsne_df_after <- as.data.frame(tsne_after$Y)
colnames(tsne_df_after) <- c("tSNE1", "tSNE2")

# Match metadata
meta_unique_after <- meta[!duplicated(counts_t_after), ]
tsne_df_after$Group <- meta_unique_after$Group

# ============================================
# STEP 3. Define a clean, print-friendly palette
# ============================================
group_colors <- c(
  "Healthy" = "#ffcc00",            
  "Adenoma" = "#984ea3", 
  "Cancer" = "#4daf4a"
)



# ============================================
# STEP 4. Compute centroids for labeling
# ============================================
centroids <- tsne_df_after %>%
  group_by(Group) %>%
  summarize(tSNE1 = median(tSNE1), tSNE2 = median(tSNE2))

# ============================================
# STEP 5. Nature-style t-SNE with black outlines
# ============================================
tsne_axes <- ggplot(tsne_df_after, aes(tSNE1, tSNE2)) +
  # Add black-stroked points
  geom_point(aes(fill = Group), size = 3, shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
  # Group labels
  geom_shadowtext(
    data = centroids, aes(label = Group, color = Group),
    bg.color = "white", bg.r = 0.3, fontface = "bold", size = 3.5, show.legend = FALSE
  ) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    title = "t-SNE after batch correction (CLR)"
  ) +
  theme_classic(base_family = "Helvetica") +
  theme(
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.ticks = element_line(color = "black", size = 0.4),
    axis.line = element_line(color = "gray40", size = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
    plot.margin = margin(6, 6, 6, 6)
  )

# ============================================
# STEP 6. Save publication-grade figure
# ============================================
ggsave(
  "tSNE_After_Batch_Correction_BlackOutline.tiff",
  plot = tsne_axes,
  width = 4.5, height = 4.5,
  dpi = 600, compression = "lzw", bg = "white"
)










