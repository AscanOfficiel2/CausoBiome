# CausoBiome: A Microbiome Causality and Biomarker Discovery framework for Colorectal Cancer

## 🧠 Purpose

CausoBiome is a modular, stage-aware, and functionally grounded microbiome analysis pipeline purpose-built to uncover **microbial features (species, ARGs, and VFs)** that are **functionally and causally linked** to human disease progression most notably **colorectal cancer (CRC)**.

Unlike conventional pipelines that stop at descriptive comparisons, CausoBiome is designed to **bridge microbiome discovery with translational insight** by leveraging modern statistical learning, causal inference, and external validation strategies.CausoBiome builds upon high-quality, genome-resolved metagenomic 

preprocessing pipelines from the [genome-resolved-urban-microbiome-biosurveillance](https://github.com/SuleimanAminu/genome-resolved-urban-microbiome-biosurveillance) repository, specifically `Modules 01_Bioinformatics and 02_Quality_Batch_subsetting`

---

## 🎯 Design Philosophy

CausoBiome is built around five central goals:

1. **Stage-awareness**  
   It models microbiome dynamics along **clinically meaningful transitions**, e.g.  
   `Healthy → Adenoma → Cancer`, capturing directional microbial shifts rather than static contrasts.

2. **Causality over correlation**  
   Through **Double Machine Learning (DML)** CausoBiome estimates the **average treatment effect (ATE)** of each feature while controlling for key confounders (e.g., Age, BMI, Sex).

3. **Functional resolution**  
   It analyzes both **taxonomic** (species-level) and **functional** (ARGs, VFs) features to capture mechanisms of microbial influence, including antibiotic resistance, immune modulation, and virulence.

4. **Feature robustness**  
   Using **bootstrap stability**, **permutation importance**, and **ordinal forest modeling**, the pipeline identifies **robust biomarkers** that generalize across datasets.

5. **Generalizability**  
   By incorporating **external cohort validation**, trend consistency, and statistical replication, CausoBiome ensures that its findings are not dataset-specific artifacts.

---


### 🔹 Step 1: Preprocessing and Normalization  
CausoBiome builds upon genome-resolved upstream modules (adapted from [genome-resolved-urban-microbiome-biosurveillance](https://github.com/SuleimanAminu/genome-resolved-urban-microbiome-biosurveillance)) starting from the `01_Bioinformatics` module and then proceed to the `02_Quality_batch_subsetting` module by running specifically the script `normalize_species_counts.py`.Then transition into **CausoBiome** starting from module `01_Quality_Normalization_Batch_ecology from the script: `01_mag_quality_metrics_analysis.py` to deliver:

- MAG binning & QC (completeness, contamination)
- Species abundance normalization (Bracken)
- ARG/VF normalization (CARD, VFDB + Hellinger transform)
- Batch correction using ComBat / limma

### 🔹 Step 2: Ecology and Diversity  
It performs alpha/beta diversity analyses across disease stages using:

- Shannon, Simpson, Richness indices  
- Bray–Curtis dissimilarity and NMDS ordinations  
- PERMANOVA and ANOSIM statistical tests  
- Dispersion tests and Cliff’s delta effect size metrics

### 🔹 Step 3: Machine Learning for Feature Selection  
ML classifiers (Random Forest, SVM, LR, GB, Ordinal Forest) are benchmarked using:

- Foldwise cross-validation  
- Model comparison using Accuracy, Kappa, MCC, and F1 Macro  
- Feature importance via permutation + bootstrap ranking  
- Signature score projection and retraining on top biomarkers

### 🔹 Step 4: Causal Inference via DML  
Using the **econML** framework, CausoBiome applies **LinearDML** to:

- Estimate ATE per microbial gene on CRC stage (Healthy–>Cancer)
- Control for Age, BMI, Sex using machine-learned nuisance functions
- Compute **confidence intervals**, **E-values**, and **required sample sizes**

### 🔹 Step 5: Microbial Interaction Modeling  
Microbial feature interactions are analyzed through:

- Pairwise causal effect modeling (LinearDML on gene-gene products)
- Identification of **synergistic vs. antagonistic** effects
- Network construction with **node centrality metrics**
- Visualization of causal hubs and clusters

### 🔹 Step 6: External Validation  
CausoBiome validates biomarker generalizability via:

- Mann–Whitney tests in both internal and external datasets  
- Directional trend comparison (CRC vs. Control)  
- Concordance barplots and scatterplots  
- Venn diagrams of statistical overlap

---

## 🧬 Who Should Use CausoBiome?

CausoBiome is ideal for:

- Microbiome researchers aiming for **causal inference beyond correlation**  
- Cancer biologists exploring **functional microbial signatures**  
- Clinical bioinformaticians validating **microbial biomarkers across cohorts**  
- Systems biologists modeling **microbial interactions and networks**

## Pipeline Structure

```text
CausoBiome/
├── 01_Quality_Normalization_Batch_ecology/
│   ├── 01_mag_quality_metrics_analysis.py
│   ├── 02_Taxon_counts_normalization.py
│   ├── 03_Taxon_subset_pathogenic_samples.py
│   ├── 04_taxon_batch_correction.R
│   ├── 05_taxon_ecology.R
│   ├── 06_arg_vf_normalization.py
│   ├── 07_arg_vf_batch_correction.R
│   ├── 08_arg_vf_ecology.R
│   └── 01_README_Module_QC_Normalization.md
│
├── 02_Machine_learning_Feature_selection/
│   ├── 01_ml_benchmark_ARG-VF.sh
│   ├── 02_RF_TUNING_ARG-VF.sh
│   ├── 03_ARG_VF_data_generation.py
│   ├── 04_validation_model.sh
│   ├── 05_Get_feature_importance.sh
│   ├── 06_retrain_top_20.sh
│   ├── 07_foldwise_Comparison_ARG-VF.R
│   ├── 08_ARG_VF_Ordinal_forest.R
│   ├── 09_taxon_Machine_learning_feature.sh
│   ├── 10_Taxon_Ordinal_forest.R
│   ├── 11_taxon_stats_correlations.py
│   └── 02_README_Module_Machine_Learning_FULL.md
│
├── 03_Causal_Analysis/
│   ├── 01_Complex_heatmap_arg_vf.R
│   ├── 02_causal_analysis.py
│   └── 03_README_Module_Causal_Analysis.md
│
├── 04_External_cohort_validation/
│   ├── External_cohort_Validation.ipynb
│   └── 04_README_Module_External_Validation.md
│
├── data/
│   ├── combined_ARG_VFDB_final_DATA.csv
│   ├── Metadata_Aligned_VF_ARG_CountMatrix.csv
│   └── ... (external validation matrices)
│
└── README_CausoBiome_Detailed.md

```

---

# 🚀 What Makes CausoBiome Novel?

CausoBiome differs from typical microbiome pipelines by addressing not just *what is different* but *what functionally drives disease progression*. Its novelty lies in several aspects:

## 1. Stage-Aware Ordinal Modeling

- Directly models disease trajectory: **Healthy → Adenoma → Cancer**
- Employs **Ordinal Forests** and **Double Machine Learning (DML)** to capture progression-aware microbial signals

## 2. Causal Inference Core

- Implements **DML via econML** to estimate **Average Treatment Effects (ATEs)** per feature
- Controls for covariates like **Age, Sex, BMI** using flexible machine-learned nuisance models
- Computes:
  - **E-values** for sensitivity to unmeasured confounding  
  - **Bootstrap confidence intervals** for robustness  
  - **Required sample sizes** for validation studies

## 3. Synthetic Cohort Generation

- Produces **PCA-enhanced synthetic metagenomes**
- Enables model testing under both **balanced** and **realistic class distributions**

## 4. Intervention-Ready Outputs

- Identifies **synergistic** and **antagonistic** microbial gene interactions
- Constructs **causal and co-occurrence networks** from inferred ATE effects
- Annotates ARGs and VFs with known mechanisms and **matched microbial species**

---

# ⚙️ Key Functional Highlights

| Component             | Description                                                                 |
|----------------------|-----------------------------------------------------------------------------|
| **Ecological Analysis**  | Rich alpha/beta diversity metrics, NMDS ordination, dispersion tests         |
| **Machine Learning**      | Foldwise benchmarking of RF, SVM, KNN, GB, Logistic Regression             |
| **Feature Robustness**    | Combined permutation importance + bootstrap stability (Top 20 features)    |
| **Ordinal Forest**        | Identifies rank-aware discriminative genes and taxa                       |
| **Signature Score**       | Mean and PCA1 scores for individual-level CRC burden assessment           |
| **Causal Estimation**     | ATEs from LinearDML with Random Forests as base learners                  |
| **Interaction Effects**   | Gene × gene product modeling for combined causal effects                  |
| **External Validation**   | Tests trend consistency in independent cohort (e.g., PRJEB10878)          |
| **Heatmaps & Networks**   | Visualizes consistent biomarkers and their interaction hubs               |

---

# 🔄 Extensibility

CausoBiome is designed with **modularity and disease-agnostic flexibility**:

- Supports any disease with **ordered clinical stages** (e.g., liver fibrosis, IBD, NAFLD)
- Compatible with **metabolomic**, **proteomic**, or **transcriptomic** feature matrices
- Adaptable to **longitudinal** or **time-series** microbiome datasets
---

## Pipeline Origin

CausoBiome is an extension of the upstream **genome-resolved-urban-microbiome-biosurveillance** workflow in:

**GitHub**: [genome-resolved-urban-microbiome-biosurveillance](https://github.com/SuleimanAminu/genome-resolved-urban-microbiome-biosurveillance)

- Users should start from the `01_Bioinformatics` module and then proceed to the `02_Quality_batch_subsetting` module by running specifically scripts `normalize_species_counts.py`
- Then transition into **CausoBiome** starting from module `01_Quality_Normalization_Batch_ecology from the script: `01_mag_quality_metrics_analysis.py`.

---

## ️ Requirements

- **SLURM-based HPC** environment (for `.sh` scripts)
- Python ≥ 3.11 and R ≥ 4.0
- Pip packages: `pandas`, `numpy`, `scikit-learn`, `matplotlib`, `seaborn`, `joblib`, `econml`, `networkx`, `statsmodels`
- R packages: `vegan`, `sva`, `ggplot2`, `umap`, `pairwiseAdonis`, `FSA`, etc.

##  Tools/Databases

-  **VFDB** – Virulence Factor Database  
- **CARD** – Comprehensive Antibiotic Resistance Database  
- **ComBat** – Batch correction via the `sva` R package  
- **EconML** – Causal inference library for treatment effect estimation  
- **scikit-learn** – Model training, permutation importance, ROC/AUC scoring

---

## Output Highlights

- Diversity metrics, ordination plots
- Classification metrics (F1, MCC, AUROC)
- Feature importance/stability plots
- Microbial features causal estimates (DML ATE)
- Robustness plots (E-values, bootstraps)
- Microbial interaction networks (weighted, annotated)


---

## Use Cases

- Functional microbiome biomarker discovery
- Ecological profiling of CRC microbiomes
- Translational microbiome-based risk stratification
- Design of synthetic consortia or microbial interventions

---


## License


MIT License — free to use, adapt, and cite with attribution.

This Framework, otherwise referred to CausoBiome is currently part of a manuscript under peer review.
This repository is shared under the MIT License to promote transparency and reproducibility.

We kindly request that you do not republish or repackage this methodology before journal publication.



## Citation

If you use **CausoBiome**, please cite the following manuscript:

> **Ascandari, A., Aminu, S., Benhida, R., & Daoud, R.** (2025).  
> *A Core Genome-Resolved Microbial Resistome–Virulome Hub Causally Drives Colorectal Cancer Progression*.  


## Submitted Articles Related to the Framework

> **Ascandari, A., Aminu, S., Benhida, R., & Daoud, R.** (2025).  
> *A Core Genome-Resolved Microbial Resistome–Virulome Hub Causally Drives Colorectal Cancer Progression* (under review).


## Contact
For questions, feedback, or collaboration regarding this framework, please reach out:

AbdulAziz Ascandari, PhD Researcher, Department of Chemical and Biochemical Sciences, University Mohammed VI Polytechnic (UM6P), Morocco, abdulaziz.ascandari@um6p.ma

Prof. Rachid Daoud, Group Leader & Supervisor, Department of Chemical and Biochemical Sciences, University Mohammed VI Polytechnic (UM6P), Morocco, rachid.daoud@um6p.ma
