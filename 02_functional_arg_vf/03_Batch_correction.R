#########################################################################
## Diagnostics and Figures for Combined ARG + VFDB Dataset (CLR pipeline)
#########################################################################

set.seed(45)
library(vegan)
library(ggplot2)
library(limma)
library(Rtsne)
library(dplyr)
library(compositions)
library(shadowtext)

# =========================
# Step 1. Load matrix and metadata
# =========================
expr <- read.csv("Combined_ARG_VFDB_Filtered_Matrix.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("Metadata_Aligned_to_FilteredMatrix.csv", check.names = FALSE)
rownames(meta) <- meta$Sample_ID
meta <- meta[rownames(expr), , drop = FALSE]

if (!"Group" %in% colnames(meta))
  stop("Metadata must contain a 'Group' column!")

# =========================
# Step 2. CLR transformation
# =========================
cat("Performing CLR transformation...\n")

# Add pseudocount
min_pos <- suppressWarnings(min(expr[expr > 0], na.rm = TRUE))
if (!is.finite(min_pos)) min_pos <- 1e-6
pseudo <- min(min_pos / 2, 1e-6)
expr_pseudo <- expr + pseudo

# Apply CLR column-wise (taxa), keeping samples as rows
tpm_clr <- t(apply(expr_pseudo, 1, clr))
tpm_clr <- as.data.frame(tpm_clr)
rownames(tpm_clr) <- rownames(expr)

# =========================
# Step 3. Filter zero-variance taxa
# =========================
# Filter columns (genes) with variance == 0
tpm_clr_filtered <- tpm_clr[, apply(tpm_clr, 2, function(x) var(as.numeric(x), na.rm = TRUE)) > 0]
cat("✅ Filtered CLR matrix:", ncol(tpm_clr_filtered), "taxa retained\n")

# =========================
# Step 4. PCA before correction
# =========================
pca <- prcomp(tpm_clr_filtered, scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df$Group <- meta$Group
ggsave("PCA_Before_Correction.png",
       ggplot(pca_df, aes(PC1, PC2, color = Group)) +
         geom_point(size = 2) + theme_minimal() +
         labs(title = "PCA Before Batch Correction"),
       width = 6, height = 5, dpi = 600, bg = "white")

# =========================
# Step 5. PERMANOVA before correction
# =========================
adonis_df <- bind_rows(
  as.data.frame(adonis2(tpm_clr_filtered ~ meta$Group, method = "euclidean")) %>% mutate(Factor = "Group"),
  as.data.frame(adonis2(tpm_clr_filtered ~ meta$Project, method = "euclidean")) %>% mutate(Factor = "Project"),
  as.data.frame(adonis2(tpm_clr_filtered ~ meta$Country, method = "euclidean")) %>% mutate(Factor = "Country"),
  as.data.frame(adonis2(tpm_clr_filtered ~ meta$Center_Name, method = "euclidean")) %>% mutate(Factor = "Center_Name")
)
write.csv(adonis_df, "PERMANOVA_before.csv", row.names = FALSE)

## =========================
# Step 6. t-SNE before correction (auto-tuned perplexity)
# =========================
set.seed(42)

n_samples <- nrow(tpm_clr_filtered)
# choose perplexity safely: between 5 and one-third of (n_samples-1)
perp <- max(5, min(10, floor((n_samples - 1) / 3)))
cat("Using perplexity =", perp, "for", n_samples, "samples.\n")

counts_unique <- tpm_clr_filtered[!duplicated(tpm_clr_filtered), ]

tsne_before <- Rtsne(counts_unique, dims = 2, perplexity = perp, check_duplicates = FALSE, verbose = TRUE)
tsne_df <- as.data.frame(tsne_before$Y)
colnames(tsne_df) <- c("tSNE1", "tSNE2")

meta_unique <- meta[!duplicated(tpm_clr_filtered), ]
tsne_df$Group <- meta_unique$Group

ggsave("tSNE_Before_Correction.png",
       ggplot(tsne_df, aes(tSNE1, tSNE2, color = Group)) +
         geom_point(size = 2) + theme_minimal() +
         labs(title = paste0("t-SNE Before Batch Correction (perplexity=", perp, ")")),
       width = 6, height = 5, dpi = 600, bg = "white")


# =========================
# Step 7. Batch correction (Project + Center_Name)
# =========================
cat("Applying batch correction using limma...\n")
counts_matrix <- t(tpm_clr_filtered)  # genes × samples
design <- model.matrix(~ meta$Group)
covariates <- model.matrix(~ meta$Project + meta$Center_Name)[, -1]

counts_corrected <- removeBatchEffect(counts_matrix,
                                      covariates = covariates,
                                      design = design)
counts_corrected <- t(counts_corrected)
write.csv(counts_corrected, "Combined_ARG_VFDB_CLR_batch_corrected.csv", row.names = TRUE)
cat("✅ Batch correction done (Project + Center_Name).\n")

# =========================
# Step 8. PCA & PERMANOVA after correction
# =========================
pca_after <- prcomp(counts_corrected, scale. = TRUE)
pca_df_after <- as.data.frame(pca_after$x)
pca_df_after$Group <- meta$Group
ggsave("PCA_After_Correction_CLR.png",
       ggplot(pca_df_after, aes(PC1, PC2, color = Group)) +
         geom_point(size = 2) + theme_minimal() +
         labs(title = "PCA After Batch Correction"),
       width = 6, height = 5, dpi = 600, bg = "white")

adonis_after_df <- bind_rows(
  as.data.frame(adonis2(counts_corrected ~ meta$Group, method = "euclidean")) %>% mutate(Factor = "Group"),
  as.data.frame(adonis2(counts_corrected ~ meta$Center_Name, method = "euclidean")) %>% mutate(Factor = "Center_Name"),
  as.data.frame(adonis2(counts_corrected ~ meta$Project, method = "euclidean")) %>% mutate(Factor = "Project")
)
write.csv(adonis_after_df, "PERMANOVA_after_CLR.csv", row.names = FALSE)

# =========================
# Step 9. t-SNE after correction (publication-style, auto-tuned perplexity)
# =========================
set.seed(42)

n_samples_after <- nrow(counts_corrected)
perp <- max(5, min(10, floor((n_samples_after - 1) / 3)))
cat("Using perplexity =", perp, "for", n_samples_after, "samples (after correction).\n")

counts_unique_after <- counts_corrected[!duplicated(counts_corrected), ]

tsne_after <- Rtsne(counts_unique_after, dims = 2,
                    perplexity = perp, check_duplicates = FALSE, verbose = TRUE)

# Create data frame
tsne_df_after <- as.data.frame(tsne_after$Y)
colnames(tsne_df_after) <- c("tSNE1", "tSNE2")

meta_unique_after <- meta[!duplicated(counts_corrected), ]
tsne_df_after$Group <- meta_unique_after$Group

# Compute centroids
centroids <- tsne_df_after %>%
  group_by(Group) %>%
  summarize(tSNE1 = median(tSNE1), tSNE2 = median(tSNE2))

# Colors
group_colors <- c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")


# Plot
tsne_axes <- ggplot(tsne_df_after, aes(tSNE1, tSNE2)) +
  # Points
  geom_point(aes(fill = Group), size = 3, shape = 21,
             color = "black", stroke = 0.4, alpha = 0.9) +
  
  # Smaller and offset labels (shift them slightly up and smaller font)
  geom_shadowtext(
    data = centroids,
    aes(label = Group, color = Group),
    bg.color = "white", bg.r = 0.3,
    fontface = "bold", size = 2.4,       # smaller labels
    show.legend = FALSE,
    alpha = 0.8,
    nudge_y = 3                           # vertical offset (move labels above the cluster)
  ) +
  
  # Colors and styling
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    #title = paste0("t-SNE after batch correction (CLR, perplexity=", perp, ")")
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

# Save
ggsave(
  "tSNE_After_Batch_Correction_BlackOutline_CLR.tiff",
  plot = tsne_axes,
  width = 4.5, height = 4.5,
  dpi = 600, compression = "lzw", bg = "white"
)

cat("\n✅ Publication-grade t-SNE (after correction) saved successfully.\n")

