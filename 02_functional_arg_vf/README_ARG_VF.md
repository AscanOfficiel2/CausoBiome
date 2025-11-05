# **CausoBiome ARG/VF Analysis and Causal Modeling Pipeline**

This repository contains reproducible scripts for the **antimicrobial resistance (ARG)** and **virulence factor (VF)** functional analysis of colorectal cancer (CRC) microbiomes.  
It integrates **gene-level preprocessing**, **batch correction**, **functional ecology**, **causal inference modeling (PLS–DML)**, and **external validation**.

---

## **Overview**

Scripts are modular and numbered (`01–08`) for sequential execution, covering data preparation to causal modeling and validation.

| Step | Script | Language | Description |
|------|---------|-----------|--------------|
| 01 | `01_arg_vf_preprocessing.py` | Python | Normalize CARD and VFDB gene abundances; create expression matrices |
| 02 | `02_ARG_VF_Genes_filtering.R` | R | Filter, align, and merge ARG + VF genes across cohorts |
| 03 | `03_Batch_correction.R` | R | Apply CLR transformation and batch correction with diagnostic plots |
| 04 | `04_Ecology_Prevalence_ARG_VF.R` | R | Analyze total ARG/VF loads, host-factor correlations, and ecology statistics |
| 05 | `05_Nuisance_models.py` | Python | Train machine-learning models (RF/GB) for feature benchmarking and importance |
| 06 | `06_causal_review.py` | Python | Perform PLS–Double Machine Learning causal analysis with bootstrap inference |
| 07 | `07_risk_score.py` | Python | Develop logistic risk score (DAI_logit) from PLS₁ and PLS₃; test overfitting |
| 08 | `08_external_validation_review.py` | Python | Validate PLS–DML signatures and risk models using external CRC datasets |

---

##  **Environment Setup**

### Python (≥3.9)
```bash
pip install pandas numpy seaborn matplotlib scikit-learn scipy openpyxl econml tqdm statsmodels networkx adjustText
```

###  R (≥4.3)
```r
install.packages(c(
  "limma", "compositions", "vegan", "ggplot2", "Rtsne",
  "FSA", "dplyr", "tibble", "tidyr", "ComplexHeatmap",
  "circlize", "RColorBrewer", "shadowtext"
))
```

---

## **Pipeline Summary**

### **Step 01 — ARG/VF Preprocessing**
**Script:** `01_arg_vf_preprocessing.py`  
Processes **CARD** and **VFDB** alignments, normalizes counts by total genes per sample, and generates relative abundance matrices.

**Outputs:**  
- `CRC_CARD_expression_matrix.csv`  
- `CRC_VFDB_expression_matrix.csv`

---

### **Step 02 — Gene Filtering and Alignment**
**Script:** `02_ARG_VF_Genes_filtering.R`  
Combines ARG and VF matrices, aligns with metadata, filters genes present in ≥5% of samples, and saves unified dataset.

**Outputs:**  
- `Combined_ARG_VFDB_Filtered_Matrix.csv`  
- `Metadata_Aligned_to_FilteredMatrix.csv`

---

### **Step 03 — Batch Correction (CLR)**
**Script:** `03_Batch_correction.R`  
Applies CLR transformation, removes batch effects (Project, Center), and visualizes PCA, PERMANOVA, and t-SNE before/after correction.

**Outputs:**  
- `Combined_ARG_VFDB_CLR_batch_corrected.csv`  
- `PCA_After_Correction_CLR.png`, `tSNE_After_Batch_Correction_BlackOutline_CLR.tiff`

---

### **Step 04 — ARG/VF Ecology and Prevalence**
**Script:** `04_Ecology_Prevalence_ARG_VF.R`  
Computes total ARG/VF load per sample, visualizes distributions (boxplots), and performs host-factor correlation analysis (BMI, Age, Country, Sex).

**Outputs:**  
- `Metadata_with_Total_Loads.csv`  
- `Boxplot_Total_ARG_Load.tiff`, `Boxplot_Total_VF_Load.tiff`  
- `HostFactor_TotalLoad_Associations.csv`

---

### **Step 05 — Machine Learning Benchmarking**
**Script:** `05_Nuisance_models.py`  
Trains five classifiers (RF, GB, LR, SVM, DT) for multi-class CRC stage prediction; performs nested CV, stability runs, and feature ranking.

**Outputs:**  
- `NestedCV_summary_RF_GB.csv`  
- `Consensus_Feature_Ranking_RF_GB.csv`  
- `Consensus_Top30_Barplot_Colorblind.tiff`

---

### **Step 06 — Causal Inference (PLS–DML Framework)**
**Script:** `06_causal_review.py`  
Integrates Partial Least Squares (PLS) regression with Double Machine Learning (DML) to estimate functional causal effects.  
Includes bootstrap inference, E-value sensitivity, gene-level enrichment, and mechanistic mapping (ARG/VF mechanisms).

**Outputs:**  
- `Stable_PLS_Scores.csv`, `Stable_PLS_Loadings.csv`  
- `DML_PLS_Bootstrap_Summary.csv`  
- `PLS_TopGenes_Ranked_with_ComponentStats.csv`  
- `PLS1_PLS3_Functional_Enrichment.csv`, `PLS1_PLS3_Network_Metrics_AllGenes_NoMechanism.csv`

---

### **Step 07 — Functional Risk Score (PLS₁ + PLS₃)**
**Script:** `07_risk_score.py`  
Builds a logistic regression–based risk model (`DAI_logit`) from PLS₁ (pathogenic axis) and PLS₃ (protective axis).  
Evaluates cross-validation AUC, optimism, permutation significance, and calibration.

**Outputs:**  
- `DAIlogit_OOF.csv`, `RiskTable_OOF.csv`  
- `ROC_DAIlogit_OOF.png`, `Calibration_DAIlogit_OOF.png`  
- `Overfitting_Report.txt`

---

### **Step 08 — External Validation**
**Script:** `08_external_validation_review.py`  
Validates gene-based risk model across external CRC cohorts, performs binary classification (CRC vs Control), and tests concordance of PLS₁/₃ gene directionality.

**Outputs:**  
- `External_LogReg_Publication_Metrics.csv`  
- `External_ROC_LogisticRegression_Publication.png`  
- `PLS_DML_RobustGene_Concordance_Scatter_Publication.png`  
- `PLS1_PLS3_Component_PowerSummary.csv`

---

## **Execution Order**

Run scripts sequentially:

```bash
python 01_arg_vf_preprocessing.py
Rscript 02_ARG_VF_Genes_filtering.R
Rscript 03_Batch_correction.R
Rscript 04_Ecology_Prevalence_ARG_VF.R
python 05_Nuisance_models.py
python 06_causal_review.py
python 07_risk_score.py
python 08_external_validation_review.py
```

---

##  **Notes**

- Run each step after verifying input files from the previous stage.  
- The causal inference model (Step 06) is exploratory and hypothesis-generating.  
- Figures are generated in **TIFF/PNG (600 dpi)** for publication quality.  
- Steps 05–08 require ≥64 GB RAM for large metagenomic matrices.

---

## **Authors & Acknowledgements**

Developed by  
**AbdulAziz Ascandari et al. (2025)**  
*Mohammed VI Polytechnic University, Morocco*  
️ *November 2025*  

