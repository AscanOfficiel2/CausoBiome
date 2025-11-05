# CausoBiome: A Microbiome Causality and Biomarker Discovery framework for Colorectal Cancer

## 🧠 Purpose

CausoBiome is a modular, stage-aware, and functionally grounded microbiome analysis pipeline purpose-built to uncover **microbial features (species, ARGs, and VFs)** that are **functionally and causally linked** to human disease progression most notably **colorectal cancer (CRC)**.

Unlike conventional pipelines that stop at descriptive comparisons, CausoBiome is designed to **bridge microbiome discovery with translational insight** by leveraging modern statistical learning, causal inference, and external validation strategies.CausoBiome builds upon high-quality, genome-resolved metagenomic 

preprocessing pipelines from the [genome-resolved-urban-microbiome-biosurveillance](https://github.com/SuleimanAminu/genome-resolved-urban-microbiome-biosurveillance) repository, specifically `Module 01_Bioinformatics` and run the scripts `01_run.sh`  and then `02_fastq_screen.sh` before continuing to `module 01` of Causobiome.

---

## 🎯 Design Philosophy

CausoBiome is built around five central goals:

1. **Stage-awareness**  
   It models microbiome dynamics along **clinically meaningful transitions**, e.g.  
   `Healthy → Adenoma → Cancer`, capturing directional microbial shifts rather than static contrasts.

2. **Causality over correlation**  
   Through **Partial Least squares (PLS) dimentionality reduction of features and Double Machine Learning (DML)** CausoBiome estimates the **average treatment effect (ATE)** of each PLS component while controlling for key confounders (e.g., Age, BMI, Sex).

3. **Functional resolution**  
   It analyzes both **taxonomic** (species-level) and **functional** (ARGs, VFs) features to capture mechanisms of microbial influence, including antibiotic resistance, immune modulation, and virulence.

4. **Feature robustness**  
   Using **bootstrap stability**, and  **permutation importance**the pipeline identifies **robust biomarkers** that generalize across datasets.

5. **Generalizability**  
   By incorporating **external cohort validation**, trend consistency, and statistical replication, CausoBiome ensures that its findings are not dataset-specific artifacts.

---



## 🧩 Modular Architecture

CausoBiome comprises two analytical layers:

| Module | Description | Focus |
|--------|-------------|--------|
| **01_species_analysis** | Species-level ecological and machine learning analysis | Taxonomic profiling, diversity, ML biomarker discovery |
| **02_functional_arg_vf** | Functional-level (ARG & VF) causal and risk modeling | Functional inference, causal estimation, translational biomarkers |

Each module is fully self-contained and executable independently, with standardized I/O formats and high-resolution publication outputs.

---

## 📦 Module Summaries

### 🔹 **Module 01 — Species-Level Analysis**

Performs **species-level ecological modeling and predictive biomarker discovery**, covering contamination filtering, normalization, ecological diagnostics, and machine learning benchmarking.

**Key Functionalities**
- Contaminant removal and assembly QC (`FastQ Screen`, `QUAST`)
- Cross-cohort taxonomic normalization and merging
- CLR + batch correction with `limma`
- Ecological diversity (Shannon, Simpson, PERMANOVA)
- Indicator species analysis (`IndVal.g`)
- Machine learning benchmarking (RF, GB, LR, SVM, DT)
- Robust feature importance and consensus biomarkers

**Representative Outputs**
- `Species_expression_matrix.csv`  
- `PERMANOVA_R2_comparison_fixed.png`  
- `Alpha_Shannon.tiff`, `NMDS_CLR_Euclidean_AllGroups.tiff`  
- `Consensus_Top30_Barplot.tiff`

**Scripts**
```
01_Fastqscreen_contigs_quality.py
02_Taxon_counts_normalization.py
03_batch_correction_species.R
04_Reference_species_catalog_building.py
05_Overlap_with_reference.R
06_Ecology_prevalence.R
07_Indicator_Specie_Analysis.R
08_ML_training_test.py
09_ML_feature_importance.py
```

---

### 🔹 **Module 02 — Functional ARG/VF Analysis and Causal Modeling**

Performs **functional-level causal inference and translational risk modeling** using antimicrobial resistance genes (ARGs) and virulence factors (VFs).

**Key Functionalities**
- Normalization and filtering of CARD & VFDB gene hits  
- Functional ecology and total load correlation analysis  
- Causal modeling via **PLS–Double Machine Learning (DML)**  
- Identification of progression vs. protective functional axes (PLS₁, PLS₃)  
- Development of a **logistic risk score (DAI_logit)**  
- Cross-cohort validation of functional gene signatures

**Representative Outputs**
- `Combined_ARG_VFDB_CLR_batch_corrected.csv`  
- `Causal_DAG_ARG_VF_to_CRC.pdf`  
- `DML_PLS_Bootstrap_Summary.csv`  
- `ROC_DAIlogit_OOF.png`, `RiskTable_OOF.csv`  
- `External_ROC_LogisticRegression_Publication.png`

**Scripts**
```
01_arg_vf_preprocessing.py
02_ARG_VF_Genes_filtering.R
03_Batch_correction.R
04_Ecology_Prevalence_ARG_VF.R
05_Nuisance_models.py
06_DAG_analysis.R
07_causal_review.py
08_risk_score.py
09_external_validation_review.py
```

---

## ⚙️ Environment Setup

### Python (≥3.9)
```bash
pip install pandas numpy seaborn matplotlib scikit-learn scipy openpyxl econml tqdm statsmodels networkx adjustText upsetplot matplotlib-venn
```

### R (≥4.3)
```r
install.packages(c(
  "limma", "compositions", "vegan", "ggplot2", "gridExtra", "umap", "Rtsne",
  "ComplexHeatmap", "circlize", "igraph", "ggraph", "FSA", "pheatmap",
  "RColorBrewer", "ggvenn", "patchwork", "UpSetR", "dagitty"
))
```

---

## 🚀 Execution Workflow

### 🧩 Species-Level Module
```bash
python 01_Fastqscreen_contigs_quality.py
python 02_Taxon_counts_normalization.py
Rscript 03_batch_correction_species.R
python 04_Reference_species_catalog_building.py
Rscript 05_Overlap_with_reference.R
Rscript 06_Ecology_prevalence.R
Rscript 07_Indicator_Specie_Analysis.R
python 08_ML_training_test.py
python 09_ML_feature_importance.py
```

### 🧬 Functional ARG/VF Module
```bash
python 01_arg_vf_preprocessing.py
Rscript 02_ARG_VF_Genes_filtering.R
Rscript 03_Batch_correction.R
Rscript 04_Ecology_Prevalence_ARG_VF.R
python 05_Nuisance_models.py
Rscript 06_DAG_analysis.R
python 07_causal_review.py
python 08_risk_score.py
python 09_external_validation_review.py
```

---

## 🧾 Outputs Overview

| Category | Example Outputs |
|-----------|----------------|
| **Ecological Metrics** | Shannon/Simpson indices, PERMANOVA tables |
| **Functional Profiles** | ARG/VF load, mechanism abundance tables |
| **Causal Estimates** | ATEs, bootstrap intervals, E-values |
| **ML Benchmarks** | Accuracy, ROC-AUC, feature stability plots |
| **Risk Stratification** | DAI_logit scores, calibration curves |
| **Validation Metrics** | Concordance and directionality across cohorts |

---

## 🧪 Notes

- Each submodule is self-contained and can be executed independently.  
- DAG analysis defines the adjustment set for DML causal estimation.  
- All figures are saved as **TIFF (600 dpi)** ready for publication.  
- Steps involving DML or ML benchmarking require ≥64 GB RAM for large datasets.  
- Causal inference results are exploratory; biological validation is recommended.

---

## Pipeline Origin

CausoBiome is an extension of the upstream **genome-resolved-urban-microbiome-biosurveillance** workflow in:

**GitHub**: [genome-resolved-urban-microbiome-biosurveillance](https://github.com/SuleimanAminu/genome-resolved-urban-microbiome-biosurveillance)

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
- Classification metrics (F1, AUROC etc.)
- Feature importance/stability plots
- Microbial features causal estimates (DML ATE)
- Robustness plots (E-values, bootstraps)
- Gene-gene interaction networks (weighted, annotated)


---

## Use Cases

- Functional microbiome biomarker discovery
- Ecological profiling of CRC microbiomes
- Translational functional microbiome-based risk stratification

---


## License


MIT License — free to use, adapt, and cite with attribution.

This Framework, otherwise referred to CausoBiome is currently part of a manuscript under peer review.
This repository is shared under the MIT License to promote transparency and reproducibility.

We kindly request that you do not republish or repackage this methodology before journal publication.



## Citation

If you use **CausoBiome**, please cite the following manuscript:

> **Ascandari, A., Aminu, S., Benhida, R., & Daoud, R.** (2025).  
> *Genome-resolved metagenomics with causal modelling implicates a resistome and virulome module in colorectal cancer*.  


## Submitted Articles Related to the Framework

> **Ascandari, A., Aminu, S., Benhida, R., & Daoud, R.** (2025).  
> *Genome-resolved metagenomics with causal modelling implicates a resistome and virulome module in colorectal cancer* (under review; npj Biofilms and Microbiomes).


## Contact
For questions, feedback, or collaboration regarding this framework, please reach out:

AbdulAziz Ascandari, PhD Researcher, Department of Chemical and Biochemical Sciences, University Mohammed VI Polytechnic (UM6P), Morocco, abdulaziz.ascandari@um6p.ma

Prof. Rachid Daoud, Group Leader & Supervisor, Department of Chemical and Biochemical Sciences, University Mohammed VI Polytechnic (UM6P), Morocco, rachid.daoud@um6p.ma
