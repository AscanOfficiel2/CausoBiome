---
title: "Module 02 — Machine Learning and Feature Selection"
output:
  pdf_document: default
  html_document:
    df_print: paged
date: "2025-06-30"
---

``` {r setup, include=FALSE}
 knitr::opts_chunk$set(echo = TRUE)
```

# Module 02 — Machine Learning and Feature Selection

**Part of the CausoBiome Pipeline**  
Lead Developer: AbdulAziz Ascandari  
Affiliation: Mohammed VI Polytechnic University  
Version: 1.0

## Objective

This module evaluates microbial features (species, ARGs, and VFs) for classifying colorectal cancer (CRC) progression stages. It performs ML benchmarking, hyperparameter tuning, feature importance analysis, synthetic cohort validation, and robust biomarker extraction.

---

## Scripts Overview

### `01_ml_benchmark_ARG-VF.sh`
- **Purpose**: Benchmarks 5 classifiers (RF, SVM, LR, KNN, GB) using microbial data, with and without metadata.
- **Outputs**:
  - `foldwise_scores_microbial_only.csv`
  - `foldwise_scores_microbial_plus_metadata.csv`
  - `full_model_benchmark_comparison.csv`

---

### `02_RF_TUNING_ARG-VF.sh`
- **Purpose**: Hyperparameter tuning of Random Forest using grid search on microbial data only.
- **Outputs**:
  - `rf_microbial_tuning_results.csv`
  - `best_rf_model_microbial.pkl`
  - `rf_confusion_matrix_plot.png`

---

### `03_ARG_VF_data_generation.py`
- **Purpose**: Generates synthetic datasets using PCA-enhanced sampling (balanced and real-ratio).
- **Outputs**:
  - `PCA_synthetic_hellinger_balanced_PCA.csv`
  - `PCA_synthetic_metadata_balanced_PCA.csv`
  - `ks_test_top_100_features.csv`, silhouette scores
  - `pca_real_vs_synthetic.png`, `tsne_real_vs_synthetic.png`

---

### `04_validation_model.sh`
- **Purpose**: Evaluates trained RF model on synthetic datasets.
- **Outputs**:
  - `classification_report_balanced.csv`, `roc_curve_balanced.png`
  - `confusion_matrix_balanced.csv`, `rf_external_validation_results.csv`

---

### `05_Get_feature_importance.sh`
- **Purpose**: Calculates permutation importance and bootstrap stability of features.
- **Outputs**:
  - `robust_biomarkers.csv`
  - `top_robust_biomarkers.png`

---

### `06_retrain_top_20.sh`
- **Purpose**: Retrains RF model on top 20 robust biomarkers and evaluates performance.
- **Outputs**:
  - `rf_reduced_biomarker_performance.csv`

---

### `07_foldwise_Comparison_ARG-VF.R`
- **Purpose**: Visualizes classifier performance across metrics.
- **Output**: `model_comparison_plot_with_metadata.tiff`

---

### `08_ARG_VF_Ordinal_forest.R`
- **Purpose**: Fits ordinal forest model on ARG/VF data; ranks features and plots top 20.
- **Outputs**:
  - `ordinal_forest_feature_importance.csv`
  - `ordinal_forest_top20_features.png`
  - `ordinal_forest_multiclass_ROC.png`

---

### `09_taxon_Machine_learning_feature.sh`
- **Purpose**: Full ML pipeline on species-level features with bootstrap feature stability.
- **Outputs**:
  - `clf_benchmark_metrics.csv`, `rf_gridsearch_results.csv`
  - `robust_biomarkers.csv`, `top_robust_biomarkers.png`

---

### `10_Taxon_Ordinal_forest.R`
- **Purpose**: Ordinal forest on species-level data, with ROC and top 20 feature plots.
- **Outputs**:
  - `ordinal_forest_feature_importance.csv`
  - `ordinal_forest_top20_features.png`
  - `ordinal_forest_multiclass_ROC.png`

---

### `11_taxon_stats_correlations.py`
- **Purpose**: Statistical evaluation of top species using Kruskal-Wallis, Dunn’s test, and correlation networks.
- **Outputs**:
  - `kruskal_results.csv`, `dunn_posthoc_results.csv`
  - `crc_trend_clustered_biomarkers.csv`, `crc_genus_beeswarm.png`
  - `crc_topspecies.png`, `cooccurrence_network_curved_edges.png`

---

## Output Summary Table

| Step | Output                                  | Description |
|------|-----------------------------------------|-------------|
| 01   | `foldwise_scores_*.csv`                 | Cross-validation scores |
| 02   | `rf_confusion_matrix_plot.png`          | Tuning results visual |
| 03   | `ks_test_top_100_features.csv`          | Synthetic vs. real validation |
| 04   | `roc_curve_balanced.png`                | External validation ROC |
| 05   | `top_robust_biomarkers.png`             | Robust features |
| 06   | `rf_reduced_biomarker_performance.csv`  | Retrained model results |
| 07   | `model_comparison_plot_with_metadata.tiff` | Classifier boxplot |
| 08   | `ordinal_forest_top20_features.png`     | ARG/VF ranking |
| 10   | `ordinal_forest_top20_features.png`     | Species ranking |
| 11   | `crc_genus_beeswarm.png`, `cooccurrence_network_curved_edges.png` | Statistics + correlations |

---

*End of Module 02 README*