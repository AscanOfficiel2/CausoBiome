# === Load Required Libraries ===
library(ordinalForest)
library(readr)
library(dplyr)
library(ggplot2)
library(pROC)
library(tidyr)
library(scales)

# === Clean Workspace ===
rm(list = ls())
gc()

# === 1. Load & Transpose gene × sample matrix ===
expr_raw <- as.data.frame(read_csv("combined_ARG_VFDB_final_DATA.csv", show_col_types = FALSE))
rownames(expr_raw) <- expr_raw[[1]]
expr_raw <- expr_raw[, -1, drop = FALSE]
expr_mat <- as.matrix(expr_raw)
expr_t   <- t(expr_mat)
rownames(expr_t) <- colnames(expr_mat)
colnames(expr_t) <- rownames(expr_mat)
expr_df <- as.data.frame(expr_t)

# === 2. Load Metadata ===
metadata <- read_csv("Metadata_Aligned_VF_ARG_CountMatrix.csv", show_col_types = FALSE)
metadata$Sample_ID <- trimws(as.character(metadata$Sample_ID))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Sample_ID
metadata <- metadata[, -which(names(metadata) == "Sample_ID"), drop = FALSE]

# === 3. Subset to Overlapping Samples ===
common_samples <- intersect(rownames(expr_df), rownames(metadata))
expr_df  <- expr_df[common_samples, , drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

# === 4. Assemble Final Data Frame ===
df <- expr_df
df$Group <- factor(metadata$Group, levels = c("Healthy", "Adenoma", "Cancer"), ordered = TRUE)
df <- as.data.frame(df)

# === 5. Remove Zero-Variance Predictors ===
gene_matrix <- df[, -ncol(df), drop = FALSE]
nzv <- vapply(gene_matrix, function(x) length(unique(x)) > 1, logical(1))
df <- cbind(gene_matrix[, nzv, drop = FALSE], Group = df$Group)
df <- as.data.frame(df)
df$Group <- factor(df$Group, levels = c("Healthy", "Adenoma", "Cancer"), ordered = TRUE)

# === 6. Fit Ordinal Forest Model ===
set.seed(42)
ordforest_model <- ordfor(
  depvar = "Group",
  data = df,
  nsets = 1000,
  ntreeperdiv = 100,
  ntreefinal = 500,
  perffunction = "probability"
)

# === 7. Print Model Summary ===
cat("\n=== Ordinal Forest Model Summary ===\n")
print(ordforest_model$classification.table)
cat(sprintf("\nOverall OOB Accuracy: %.3f\n\n", ordforest_model$overall.classification.rate))

# === 8. Save Variable Importance to CSV ===
importance_df <- data.frame(
  Gene = names(ordforest_model$varimp),
  Importance = ordforest_model$varimp
) %>% arrange(desc(Importance))

write_csv(importance_df, "ordinal_forest_feature_importance.csv")

# === 9. Plot Top 20 Genes by Importance ===
if (nrow(importance_df) > 0) {
  top_n <- min(20, nrow(importance_df))
  top_plot <- ggplot(head(importance_df, top_n),
                     aes(x = reorder(Gene, Importance), y = Importance)) +
    geom_col(fill = "#2C3E50") +
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Top", top_n, "Genes by Ordinal Forest Importance"),
      x = "Gene",
      y = "Permutation Importance"
    )
  print(top_plot)
} else {
  message("importance_df is empty; no barplot generated.")
}

# === 10. Multiclass ROC Curve ===
# Get predicted probabilities
pred_probs <- predict(ordforest_model, newdata = df, type = "probabilities")$classprobs
colnames(pred_probs) <- levels(df$Group)  # Assign class names manually

true_labels <- df$Group
classes <- levels(true_labels)

# Reshape to long format
roc_data <- cbind(as.data.frame(pred_probs), true = true_labels) %>%
  pivot_longer(cols = -true, names_to = "Class", values_to = "Probability")

# Compute ROC for each class
roc_list <- lapply(classes, function(cls) {
  bin_true <- as.numeric(true_labels == cls)
  roc_obj <- roc(response = bin_true, predictor = pred_probs[, cls])
  list(class = cls, roc = roc_obj)
})

# Plot ROC curves
plot(roc_list[[1]]$roc, col = hue_pal()(3)[1], lwd = 2.5, legacy.axes = TRUE,
     main = "", print.auc = TRUE)

for (i in 2:length(roc_list)) {
  plot(roc_list[[i]]$roc, col = hue_pal()(3)[i], lwd = 2.5, add = TRUE)
}

legend("bottomright",
       legend = paste0("Class: ", sapply(roc_list, `[[`, "class"),
                       " (AUC = ", sapply(roc_list, function(x) round(auc(x$roc), 3)), ")"),
       col = hue_pal()(3), lwd = 3, cex = 0.9, box.lty = 0)



###############################################################################
# === 9. Plot Top 20 Genes by Importance ===
if (nrow(importance_df) > 0) {
  top_n <- min(20, nrow(importance_df))
  
  # Open PNG device
  png("ordinal_forest_top20_features.png", width = 7, height = 5, units = "in", res = 600)
  
  # Create barplot
  top_plot <- ggplot(head(importance_df, top_n),
                     aes(x = reorder(Gene, Importance), y = Importance)) +
    geom_col(fill = "blue", width = 0.7, color = "blue", linewidth = 0.7) +  # thick bar border
    coord_flip() +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(color = "black", size = 2.5),  # thick outer plot border
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5)
    ) +
    labs(
      title = paste("Top", top_n, "Genes by Ordinal Forest Importance"),
      x = "Gene",
      y = "Permutation Importance"
    )
  print(top_plot)
  dev.off()
} else {
  message("importance_df is empty; no barplot generated.")
}

# === 10. Multiclass ROC Curve ===

# Get predicted probabilities
pred_probs <- predict(ordforest_model, newdata = df, type = "probabilities")$classprobs
colnames(pred_probs) <- levels(df$Group)

true_labels <- df$Group
classes <- levels(true_labels)

# Compute ROC for each class
roc_list <- lapply(classes, function(cls) {
  bin_true <- as.numeric(true_labels == cls)
  roc_obj <- roc(response = bin_true, predictor = pred_probs[, cls], quiet = TRUE)
  list(class = cls, roc = roc_obj)
})

# Open PNG device
png("ordinal_forest_multiclass_ROC.png", width = 5, height = 5, units = "in", res = 600)

# Set margins and plot first ROC
par(mar = c(5, 5, 4, 2) + 0.1)
plot(roc_list[[1]]$roc, col = hue_pal()(3)[1], lwd = 2.5, legacy.axes = TRUE,
     main = "", print.auc = TRUE,
     xlab = "1 - Specificity", ylab = "Sensitivity", cex.lab = 1.3, cex.main = 1.4)

# Plot remaining ROCs
for (i in 2:length(roc_list)) {
  plot(roc_list[[i]]$roc, col = hue_pal()(3)[i], lwd = 2.5, add = TRUE)
}

# Diagonal reference line
abline(a = 0, b = 1, lty = 2, col = "gray60", lwd = 1.2)


# Add legend
legend("bottomright",
       legend = paste0("Class: ", sapply(roc_list, `[[`, "class"),
                       " (AUC = ", sapply(roc_list, function(x) round(auc(x$roc), 3)), ")"),
       col = hue_pal()(3), lwd = 3, cex = 0.95, box.lty = 0)

dev.off()

