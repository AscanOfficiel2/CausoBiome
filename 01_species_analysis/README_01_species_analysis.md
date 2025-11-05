# 01_species_analysis Module — Species-Level Preprocessing and Machine Learning Workflow

**Author:** AbdulAziz Ascandari, Mohammed VI Polytechnic University  
**Date:** November 2025  

---

## Overview

The **`01_species_analysis`** module forms the first stage of the **CausoBiome** pipeline and prepares species-level features for downstream ecological, functional, and machine-learning analyses.  
It integrates **read contamination assessment**, **taxonomic normalization**, **batch correction**, **reference overlap**, **ecological diversity**, and **predictive modeling** into a unified, reproducible framework.

This document describes each sub-module (01–09), input/output structure, and key visualization outputs.

---

## Module Structure

| Step | Script | Language | Description |
|------|---------|-----------|--------------|
| 01 | `01_Fastqscreen_contigs_quality.py` | Python | Evaluates contaminant reads (pre/post filtering) and contig assembly quality (N50, L50, total length) |
| 02 | `02_Taxon_counts_normalization.py` | Python | Normalizes and merges species-level Bracken counts |
| 03 | `03_batch_correction_species.R` | R | Applies CLR transformation and `limma`-based batch correction |
| 04 | `04_Reference_species_catalog_building.py` | Python | Builds a reference catalog of CRC-associated species |
| 05 | `05_Overlap_with_reference.R` | R | Quantifies overlap with literature markers and visualizes reproducible taxa |
| 06 | `06_Ecology_prevalence.R` | R | Computes alpha/beta diversity, prevalence, and ecological composition |
| 07 | `07_Indicator_Specie_Analysis.R` | R | Identifies statistically significant indicator taxa |
| 08 | `08_ML_training_test.py` | Python | Benchmarks multiple ML models for CRC stage classification |
| 09 | `09_ML_feature_importance.py` | Python | Computes robust feature importance and consensus biomarkers |

---

## Step 01 — FASTQ Screen and Contig Quality Assessment

**Script:** `01_Fastqscreen_contigs_quality.py`  
**Purpose:** Evaluate contaminant removal efficiency and assess assembly quality metrics before proceeding to taxonomic profiling.

### Workflow Summary
1. Reads **FastQ Screen** reports before and after decontamination (`Fastq_screen_before.csv`, `Fastq_screen_after.csv`).
2. Cleans and harmonizes column names.
3. Calculates mapped and unmapped reads across categories (human, mouse, adapters, vectors, *Plasmodium*).
4. Aggregates by cohort and computes percentage change in mapped reads.
5. Generates visual summaries:
   - **Figure 2A:** Bar plot showing contaminant read counts before vs. after filtering.
   - **Supplementary Figure 1B:** Boxplots for QUAST assembly metrics (`N50`, `#Contigs`, `L50`, `Total Length`).

### Key Outputs
| Output File | Description |
|--------------|-------------|
| `CRC_final_corrected_mapped_reads_comparison.csv` | Summary table of mapped/unmapped read improvements |
| `CRC_Contaminants_reads.png` | Bar chart (Before vs After filtering) |
| `CRC_Quast_combined_boxplots.png` | Boxplots of contig assembly quality metrics |

### Scientific Context
This step ensures that host DNA, vector, or adapter contaminants are removed effectively before ecological inference. The quality plots provide a visual confirmation of improved mapping accuracy and contig assembly stability across cohorts.

---

## Step 02 — Taxon Count Normalization

**Script:** `02_Taxon_counts_normalization.py`  
Normalizes Bracken outputs, adjusts for sequencing depth, and merges across cohorts into a unified `Species_expression_matrix.csv`.

---

## Step 03 — Batch Correction

**Script:** `03_batch_correction_species.R`  
Performs CLR transformation, removes batch effects using `limma::removeBatchEffect`, and visualizes pre/post correction through PCA, UMAP, and t-SNE.

---

## Step 04 — Reference Catalog Construction

**Script:** `04_Reference_species_catalog_building.py`  
Builds a unified CRC reference panel integrating *Wirbel et al. (2019)* and *Piccinno et al. (2025)* species markers for benchmarking reproducibility.

---

## Step 05 — Overlap with Literature

**Script:** `05_Overlap_with_reference.R`  
Computes intersection between study-detected species and reference panels, generating reproducibility heatmaps and marker networks.

---

## Step 06 — Microbial Ecology and Prevalence

**Script:** `06_Ecology_prevalence.R`  
Computes diversity indices (Shannon, Simpson, Richness), PERMANOVA tests, and prevalence heatmaps, generating UpSet and bar plots for core and peripheral species.

---

## Step 07 — Indicator Species Analysis

**Script:** `07_Indicator_Specie_Analysis.R`  
Performs indicator value analysis (`IndVal.g`) to identify taxa specifically enriched in Healthy, Adenoma, or Cancer groups.

---

## Step 08 — Machine Learning Benchmarking

**Script:** `08_ML_training_test.py`  
Benchmarks Random Forest, Logistic Regression, SVM, Decision Tree, and Gradient Boosting classifiers using nested cross-validation and holdout evaluation.

---

## Step 09 — Feature Importance and Stability

**Script:** `09_ML_feature_importance.py`  
Computes robust, reproducible feature rankings through repeated model fitting and permutation importance; produces consensus biomarker visualizations.

---

## Environment Setup

### Python (≥3.9)
```bash
pip install pandas numpy seaborn matplotlib scikit-learn scipy openpyxl upsetplot matplotlib-venn
```

### R (≥4.3)
```r
install.packages(c(
  "limma", "compositions", "vegan", "ggplot2", "gridExtra", "umap", "Rtsne",
  "ComplexHeatmap", "circlize", "igraph", "ggraph", "FSA", "pheatmap",
  "RColorBrewer", "ggvenn", "patchwork", "UpSetR"
))
```

---

## Execution Order

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

---

## Citation
If you use this module or figures generated from it, please cite:  
**Ascandari A. et al. (2025).** *A Unified Computational Causal Inference Framework for Reproducible Microbiome-Based Biomarkers to Enhance Precision Therapies.*  
Mohammed VI Polytechnic University, Morocco.
