
#  CRC Microbiome Classification & Feature Importance Pipeline

This module contains scripts for supervised machine learning modeling, hyperparameter tuning, model benchmarking, feature selection, and performance visualization for microbiome-based classification of colorectal cancer (CRC) stages.

---

## Scripts Overview

### 04_RF_Hyperparametertuning.sh

**Purpose:**  
Tunes a Random Forest classifier for microbiome-only data using grid search and 10-fold cross-validation.

**Inputs:**
- `Species_matrix_transposed.csv`
- `Aligned_metadata.csv`

**Outputs:**
- Best model parameters
- `confusion_matrix.png`
- `roc_curve_without_confounders.png`
- `best_RF_MODEL.pkl`

---

### 05_Benchmarking_ML_models.sh

**Purpose:**  
Benchmarks five machine learning classifiers on microbiome data using stratified 10-fold cross-validation.

**Models:**
- Random Forest
- Logistic Regression
- Support Vector Machine (SVM)
- k-Nearest Neighbors (KNN)
- Gradient Boosting

**Scoring Metrics:**
- Accuracy
- F1 Macro
- Cohen’s Kappa
- Matthews Correlation Coefficient (MCC)

**Outputs:**
- `CRC_model_comparison_results.csv`
- `CRC_model_foldwise_scores.csv`

---

### 06_RF_Permutation_Stability_Importance.sh

**Purpose:**  
Performs permutation importance and bootstrap stability analysis on Random Forest features.

**Key Outputs:**
- `permutation_importance_species_only.csv`
- `bootstrap_feature_stability_species_only.csv`
- `youden_index_species_only.csv`
- `top20_species_permutation_barplot.png`
- `consistent_species_feature_summary.csv`

---

### 07_Comparison_Script.R

**Purpose:**  
Visualizes fold-wise model performance using tidy boxplots.

**Input:**  
- `CRC_model_foldwise_scores_with_metadata.csv`

**Output:**  
- `model_comparison_plot_with_metadata.tiff`

---

## Dependencies

### Python
```bash
pip install pandas numpy scikit-learn matplotlib seaborn joblib
```

### R
```r
install.packages(c("ggplot2", "tidyverse"))
```

---

## Workflow Summary

```text
04_RF_Hyperparametertuning.sh
  └── Train and optimize RF model

05_Benchmarking_ML_models.sh
  └── Benchmark RF, SVM, KNN, GB, and LR

06_RF_Permutation_Stability_Importance.sh
  └── Feature selection: permutation + bootstrap analysis

07_Comparison_Script.R
  └── Visualize model performance (boxplots of metrics)
```

---

##  Notes

- All scripts assume a SLURM-based HPC environment.
- Update `base_path` and `file_path` variables as per your directory structure.
- Ensure Python 3.11+ and a configured virtual environment are available for SLURM compute nodes.

---

