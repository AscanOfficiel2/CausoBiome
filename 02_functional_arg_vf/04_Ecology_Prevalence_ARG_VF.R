######################################################################
# ARG + VFDB Metagenomic Ecology Analysis
# Total ARG/VF Load + Host Factor Influence (with posthoc & regression)
######################################################################
set.seed(42)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(FSA)
  library(tibble)
  library(tidyr)
})

######################################################################
# Helper Functions
######################################################################

clean_ids <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_\\-]", "", x)
  toupper(x)
}

align_matrix_to_meta <- function(mat, meta, sample_col = "Sample_ID", name = "matrix") {
  meta[[sample_col]] <- clean_ids(meta[[sample_col]])
  rownames(meta) <- meta[[sample_col]]
  rn <- clean_ids(rownames(mat))
  cn <- clean_ids(colnames(mat))
  rownames(mat) <- rn
  colnames(mat) <- cn
  
  shared_r <- intersect(rownames(mat), rownames(meta))
  shared_c <- intersect(colnames(mat), rownames(meta))
  if (length(shared_r) == 0 && length(shared_c) == 0) {
    stop(paste0("❌ No shared samples between ", name, " and metadata."))
  }
  if (length(shared_c) > length(shared_r)) {
    mat <- t(mat)
    rownames(mat) <- clean_ids(rownames(mat))
    shared <- intersect(rownames(mat), rownames(meta))
  } else {
    shared <- shared_r
  }
  
  mat  <- mat[shared, , drop = FALSE]
  meta <- meta[shared, , drop = FALSE]
  stopifnot(all(rownames(mat) == rownames(meta)))
  
  colnames(meta) <- make.names(colnames(meta), unique = TRUE)
  cat("✅", name, "aligned →", length(shared), "samples and", ncol(mat), "features.\n")
  list(mat = mat, meta = meta)
}

detect_feature_sets <- function(feat_names) {
  arg_patterns <- c("^ARG", "ARG", "AMR", "RESIST", "RESISTANCE", "ARO", "CARD", "RESFAM")
  vf_patterns  <- c("^VF", "VFDB", "VIRULENCE", "VIR", "TOX", "TOXIN", "ADHESIN", "VFG")
  pick_cols <- function(pats) {
    idx <- logical(length(feat_names))
    for (p in pats) idx <- idx | grepl(p, feat_names, ignore.case = TRUE)
    which(idx)
  }
  list(
    arg_idx = pick_cols(arg_patterns),
    vf_idx  = pick_cols(vf_patterns)
  )
}

######################################################################
# 1️⃣ Load Data and Calculate Total ARG/VF Load
######################################################################
cat("📂 Loading TPM matrix and metadata...\n")
tpm_raw <- read.csv("Combined_ARG_VFDB_Filtered_Matrix.csv", row.names = 1, check.names = FALSE)
meta_in  <- read.csv("Metadata_Aligned_to_FilteredMatrix.csv", check.names = FALSE)
aln <- align_matrix_to_meta(tpm_raw, meta_in, sample_col = "Sample_ID", name = "TPM")
tpm <- aln$mat
meta <- aln$meta
group_colors <- c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")

feat_sets <- detect_feature_sets(colnames(tpm))
arg_idx <- feat_sets$arg_idx
vf_idx  <- feat_sets$vf_idx
tpm_ARG <- if (length(arg_idx) > 0) tpm[, arg_idx, drop = FALSE] else NULL
tpm_VF  <- if (length(vf_idx)  > 0) tpm[, vf_idx,  drop = FALSE] else NULL
cat("🔎 Feature split (TPM) → ARG:", ncol(tpm_ARG), "| VF:", ncol(tpm_VF), "\n")

# --- Total Load Calculation ---
meta$Total_ARG_Load <- if (!is.null(tpm_ARG)) rowSums(tpm_ARG, na.rm = TRUE) else NA
meta$Total_VF_Load  <- if (!is.null(tpm_VF)) rowSums(tpm_VF,  na.rm = TRUE) else NA
write.csv(meta, "Metadata_with_Total_Loads.csv", row.names = TRUE)
cat("✅ Total ARG & VF loads calculated and saved.\n")

######################################################################
# 2️⃣ Boxplots of Total ARG & VF Load
######################################################################
# Make sure order is Healthy, Adenoma, Cancer
meta$Group <- factor(meta$Group, levels = c("Healthy", "Adenoma", "Cancer"))


plot_total_load <- function(df, load_col, ylab, filename) {
  p <- ggplot(df, aes(x = Group, y = .data[[load_col]], fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, color = "black", linewidth = 1.1) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, color = "black") +
    scale_fill_manual(values = group_colors) +
    labs(x = "", y = ylab, title = "") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1.3),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      legend.position = "none"
    )
  ggsave(filename, p, width = 3.5, height = 3.5, dpi = 600)
  cat(paste0("✅ Saved: ", filename, "\n"))
}

plot_total_load(meta, "Total_ARG_Load", "Total ARG load", "Boxplot_Total_ARG_Load.tiff")
plot_total_load(meta, "Total_VF_Load", "Total VF load", "Boxplot_Total_VF_Load.tiff")

######################################################################
# 3️⃣ Host Factor Influence on Total ARG/VF Load (Kruskal/Spearman + Posthoc)
######################################################################
cat("\n📊 Testing host-factor influence on Total ARG/VF load...\n")

meta$Total_ARG_Load <- as.numeric(meta$Total_ARG_Load)
meta$Total_VF_Load  <- as.numeric(meta$Total_VF_Load)

host_factors_cat  <- c("Country", "Sex")
host_factors_cont <- c("Age", "BMI")

results <- data.frame(
  Factor = character(),
  Type   = character(),
  Test   = character(),
  Effect = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

# ---------- A) Categorical host factors (Kruskal–Wallis + Dunn’s post-hoc) ----------
for (factor in host_factors_cat) {
  if (factor %in% colnames(meta)) {
    cat(paste0("\n▶ ", factor, " effect:\n"))
    if (length(unique(na.omit(meta[[factor]]))) > 1) {
      # ARG load
      kw_arg <- tryCatch(kruskal.test(Total_ARG_Load ~ as.factor(meta[[factor]]), data = meta),
                         error = function(e) NULL)
      if (!is.null(kw_arg)) {
        cat(paste0("  ARG load — Kruskal-Wallis p = ", round(kw_arg$p.value, 4), "\n"))
        results <- rbind(results, data.frame(
          Factor = factor, Type = "Categorical", Test = "Kruskal-Wallis",
          Effect = "ARG load", p_value = kw_arg$p.value
        ))
        if (kw_arg$p.value < 0.05) {
          dunn_arg <- FSA::dunnTest(Total_ARG_Load ~ as.factor(meta[[factor]]),
                                    data = meta, method = "bh")$res
          file <- paste0("PostHoc_Dunn_ARG_", factor, ".csv")
          write.csv(dunn_arg, file, row.names = FALSE)
          cat("    ✅ Saved Dunn post-hoc for ARG → ", file, "\n")
        }
      }
      # VF load
      kw_vf <- tryCatch(kruskal.test(Total_VF_Load ~ as.factor(meta[[factor]]), data = meta),
                        error = function(e) NULL)
      if (!is.null(kw_vf)) {
        cat(paste0("  VF load — Kruskal-Wallis p = ", round(kw_vf$p.value, 4), "\n"))
        results <- rbind(results, data.frame(
          Factor = factor, Type = "Categorical", Test = "Kruskal-Wallis",
          Effect = "VF load", p_value = kw_vf$p.value
        ))
        if (kw_vf$p.value < 0.05) {
          dunn_vf <- FSA::dunnTest(Total_VF_Load ~ as.factor(meta[[factor]]),
                                   data = meta, method = "bh")$res
          file <- paste0("PostHoc_Dunn_VF_", factor, ".csv")
          write.csv(dunn_vf, file, row.names = FALSE)
          cat("    ✅ Saved Dunn post-hoc for VF → ", file, "\n")
        }
      }
    }
  }
}

# ---------- B) Continuous host factors (Spearman correlations + regression plots) ----------
for (factor in host_factors_cont) {
  if (factor %in% colnames(meta)) {
    cat(paste0("\n▶ Continuous factor: ", factor, "\n"))
    
    cor_arg <- suppressWarnings(cor.test(meta[[factor]], meta$Total_ARG_Load, method = "spearman"))
    cor_vf  <- suppressWarnings(cor.test(meta[[factor]], meta$Total_VF_Load,  method = "spearman"))
    
    results <- rbind(results,
                     data.frame(Factor = factor, Type = "Continuous", Test = "Spearman",
                                Effect = paste0("ARG load (rho=", round(cor_arg$estimate, 3), ")"),
                                p_value = cor_arg$p.value),
                     data.frame(Factor = factor, Type = "Continuous", Test = "Spearman",
                                Effect = paste0("VF load (rho=", round(cor_vf$estimate, 3), ")"),
                                p_value = cor_vf$p.value))
    
    cat(paste0("  ARG load — rho = ", round(cor_arg$estimate, 3),
               ", p = ", round(cor_arg$p.value, 4), "\n"))
    cat(paste0("  VF load  — rho = ", round(cor_vf$estimate, 3),
               ", p = ", round(cor_vf$p.value, 4), "\n"))
    
    # ----- Regression plots -----
    # ARG
    p1 <- ggplot(meta, aes(x = .data[[factor]], y = Total_ARG_Load)) +
      geom_point(alpha = 0.7, color = "#377eb8") +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      scale_y_log10() +
      labs(x = factor, y = "Total ARG load (log10)",
           title = paste0("ARG load vs ", factor)) +
      theme_bw(base_size = 12) +
      theme(panel.border = element_rect(color = "black", linewidth = 1.2))
    ggsave(paste0("Regression_ARG_", factor, ".tiff"), p1, width = 4, height = 3.5, dpi = 600)
    
    # VF
    p2 <- ggplot(meta, aes(x = .data[[factor]], y = Total_VF_Load)) +
      geom_point(alpha = 0.7, color = "#e41a1c") +
      geom_smooth(method = "lm", se = TRUE, color = "black") +
      scale_y_log10() +
      labs(x = factor, y = "Total VF load (log10)",
           title = paste0("VF load vs ", factor)) +
      theme_bw(base_size = 12) +
      theme(panel.border = element_rect(color = "black", linewidth = 1.2))
    ggsave(paste0("Regression_VF_", factor, ".tiff"), p2, width = 4, height = 3.5, dpi = 600)
  }
}
#####################################################################################################

######################################################################
# ARG + VFDB Metagenomic Ecology Analysis
# Total ARG/VF Load + Host Factor Influence (with posthoc & regression)
######################################################################
set.seed(42)
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(FSA)
  library(tibble)
  library(tidyr)
})

######################################################################
# Helper Functions
######################################################################

clean_ids <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9_\\-]", "", x)
  toupper(x)
}

align_matrix_to_meta <- function(mat, meta, sample_col = "Sample_ID", name = "matrix") {
  meta[[sample_col]] <- clean_ids(meta[[sample_col]])
  rownames(meta) <- meta[[sample_col]]
  rn <- clean_ids(rownames(mat))
  cn <- clean_ids(colnames(mat))
  rownames(mat) <- rn
  colnames(mat) <- cn
  
  shared_r <- intersect(rownames(mat), rownames(meta))
  shared_c <- intersect(colnames(mat), rownames(meta))
  if (length(shared_r) == 0 && length(shared_c) == 0) {
    stop(paste0("❌ No shared samples between ", name, " and metadata."))
  }
  if (length(shared_c) > length(shared_r)) {
    mat <- t(mat)
    rownames(mat) <- clean_ids(rownames(mat))
    shared <- intersect(rownames(mat), rownames(meta))
  } else {
    shared <- shared_r
  }
  
  mat  <- mat[shared, , drop = FALSE]
  meta <- meta[shared, , drop = FALSE]
  stopifnot(all(rownames(mat) == rownames(meta)))
  
  colnames(meta) <- make.names(colnames(meta), unique = TRUE)
  cat("✅", name, "aligned →", length(shared), "samples and", ncol(mat), "features.\n")
  list(mat = mat, meta = meta)
}

detect_feature_sets <- function(feat_names) {
  arg_patterns <- c("^ARG", "ARG", "AMR", "RESIST", "RESISTANCE", "ARO", "CARD", "RESFAM")
  vf_patterns  <- c("^VF", "VFDB", "VIRULENCE", "VIR", "TOX", "TOXIN", "ADHESIN", "VFG")
  pick_cols <- function(pats) {
    idx <- logical(length(feat_names))
    for (p in pats) idx <- idx | grepl(p, feat_names, ignore.case = TRUE)
    which(idx)
  }
  list(
    arg_idx = pick_cols(arg_patterns),
    vf_idx  = pick_cols(vf_patterns)
  )
}

######################################################################
# 1️⃣ Load Data and Calculate Total ARG/VF Load
######################################################################
cat("📂 Loading TPM matrix and metadata...\n")
tpm_raw <- read.csv("Combined_ARG_VFDB_Filtered_Matrix.csv", row.names = 1, check.names = FALSE)
meta_in  <- read.csv("Metadata_Aligned_to_FilteredMatrix.csv", check.names = FALSE)
aln <- align_matrix_to_meta(tpm_raw, meta_in, sample_col = "Sample_ID", name = "TPM")
tpm <- aln$mat
meta <- aln$meta
group_colors <- c("Healthy" = "#ffcc00", "Adenoma" = "#984ea3", "Cancer" = "#4daf4a")

feat_sets <- detect_feature_sets(colnames(tpm))
arg_idx <- feat_sets$arg_idx
vf_idx  <- feat_sets$vf_idx
tpm_ARG <- if (length(arg_idx) > 0) tpm[, arg_idx, drop = FALSE] else NULL
tpm_VF  <- if (length(vf_idx)  > 0) tpm[, vf_idx,  drop = FALSE] else NULL
cat("🔎 Feature split (TPM) → ARG:", ncol(tpm_ARG), "| VF:", ncol(tpm_VF), "\n")

# --- Total Load Calculation ---
meta$Total_ARG_Load <- if (!is.null(tpm_ARG)) rowSums(tpm_ARG, na.rm = TRUE) else NA
meta$Total_VF_Load  <- if (!is.null(tpm_VF)) rowSums(tpm_VF,  na.rm = TRUE) else NA
write.csv(meta, "Metadata_with_Total_Loads.csv", row.names = TRUE)
cat("✅ Total ARG & VF loads calculated and saved.\n")

######################################################################
# 2️⃣ Boxplots of Total ARG & VF Load
######################################################################
# Make sure order is Healthy, Adenoma, Cancer
meta$Group <- factor(meta$Group, levels = c("Healthy", "Adenoma", "Cancer"))


plot_total_load <- function(df, load_col, ylab, filename) {
  p <- ggplot(df, aes(x = Group, y = .data[[load_col]], fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85, color = "black", linewidth = 1.1) +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.7, color = "black") +
    scale_fill_manual(values = group_colors) +
    labs(x = "", y = ylab, title = "") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", linewidth = 1.3),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      legend.position = "none"
    )
  ggsave(filename, p, width = 3.5, height = 3.5, dpi = 600)
  cat(paste0("✅ Saved: ", filename, "\n"))
}

plot_total_load(meta, "Total_ARG_Load", "Total ARG load", "Boxplot_Total_ARG_Load.tiff")
plot_total_load(meta, "Total_VF_Load", "Total VF load", "Boxplot_Total_VF_Load.tiff")

###########################################################################################
# Load necessary libraries
library(dplyr)
library(FSA)

# Load the metadata containing Total ARG and VF loads
meta <- read.csv("Metadata_with_Total_Loads.csv")

# Ensure "Group" is a factor with the correct order
meta$Group <- factor(meta$Group, levels = c("Healthy", "Adenoma", "Cancer"))

# Kruskal-Wallis Test for ARG load
kw_ARG <- kruskal.test(Total_ARG_Load ~ Group, data = meta)
cat("ARG Load - Kruskal-Wallis test:\n")
print(kw_ARG)

# Dunn's post-hoc test if significant
if (kw_ARG$p.value < 0.05) {
  dunn_ARG <- dunnTest(Total_ARG_Load ~ Group, data = meta, method = "bh")
  cat("\nARG Load - Dunn's post-hoc test (FDR corrected):\n")
  print(dunn_ARG$res)
  write.csv(dunn_ARG$res, "ARG_Load_Dunn_PostHoc.csv", row.names = FALSE)
}

# Kruskal-Wallis Test for VF load
kw_VF <- kruskal.test(Total_VF_Load ~ Group, data = meta)
cat("\nVF Load - Kruskal-Wallis test:\n")
print(kw_VF)

# Dunn's post-hoc test if significant
if (kw_VF$p.value < 0.05) {
  dunn_VF <- dunnTest(Total_VF_Load ~ Group, data = meta, method = "bh")
  cat("\nVF Load - Dunn's post-hoc test (FDR corrected):\n")
  print(dunn_VF$res)
  write.csv(dunn_VF$res, "VF_Load_Dunn_PostHoc.csv", row.names = FALSE)
}

# Save Kruskal-Wallis results
kw_results <- data.frame(
  Load_Type = c("ARG", "VF"),
  KW_statistic = c(kw_ARG$statistic, kw_VF$statistic),
  KW_p_value = c(kw_ARG$p.value, kw_VF$p.value)
)
write.csv(kw_results, "Kruskal_Wallis_ARG_VF_Results.csv", row.names = FALSE)

cat("\n✅ Kruskal–Wallis and Dunn's post-hoc analyses completed and saved.\n")


######################################################################
# 4️⃣ Save All Results
######################################################################
results <- results %>% arrange(Factor, p_value)
write.csv(results, "HostFactor_TotalLoad_Associations.csv", row.names = FALSE)
cat("\n✅ Host factor analysis complete.\n✅ Kruskal-Wallis, Dunn post-hoc, and regression outputs saved.\n")

