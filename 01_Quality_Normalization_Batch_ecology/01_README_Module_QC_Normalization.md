---
title: "01_README"
output:
  pdf_document: default
  html_document:
    df_print: paged
date: "2025-06-30"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Module 01 — Quality, Normalization, Batch Correction & Ecological Profiling

**Part of the CausoBiome Pipeline**  
Lead Developer: AbdulAziz Ascandari  
Affiliation: Mohammed VI Polytechnic University  
Version: 1.0

## Objective

This module prepares high-quality, normalized, and batch-corrected microbial taxonomic and functional data for downstream biomarker discovery and causal inference in colorectal cancer (CRC). It integrates read-level quality control, taxonomic normalization, ARG/VF functional normalization, batch harmonization, and ecological diversity assessments.

## Scripts Included in This Module

### 1. `01_mag_quality_metrics_analysis.py`

**Purpose:** Assess the quality of MAGs (completeness, contamination, strain heterogeneity)

**Inputs:**  
- `SELECTED_BINS.csv` — MAG quality metadata across cohorts

**Outputs:**  
- Completeness & contamination KDE plots  
- Strain heterogeneity violin plots  
- Scatterplots and heatmaps of QC correlations  
- Marker lineage distributions  
- Regression analysis of N50 vs assembly features  

### 2. `02_Taxon_counts_normalization.py`

**Purpose:** Normalize species-level Bracken abundance using total reads per sample, then merge across cohorts.

**Inputs:**  
- Bracken outputs (fractional read counts)  
- Sample metadata (`Merged_metadata_all.csv`)

**Outputs:**  
- `Species_expression_matrix.csv` — normalized species counts  
- `Species_metadata.csv` — harmonized sample metadata

**Functionality:**  
- Read scaling, per-sample aggregation, merging of cohort tables, metadata alignment

### 3. `03_Taxon_subset_pathogenic_samples.py`

**Purpose:** Subset species matrix to include only samples positive for ARGs or VFs.

**Inputs:**  
- `General_species_matrix.csv`  
- Sample list of ARG/VF-positive samples (`vf_or_arg_samples`)

**Outputs:**  
- `Species_subset_matrix.csv` — taxon profiles for ARG/VF-carrier samples

**Context:**  
Used prior to batch correction to harmonize inputs for functional ecology.

### 4. `04_taxon_batch_correction.R`

**Purpose:** Perform batch correction on the taxonomic matrix using ComBat.

**Inputs:**  
- `Species_subset_matrix.csv`  
- `crc_meta.csv` — metadata with batch covariates

**Outputs:**  
- `combat_corrected_matrix.csv` — corrected counts  
- `hellinger_transformed_combat.csv` — Hellinger-transformed version  
- PCA, UMAP, and PERMANOVA outputs pre/post correction

**Key Methods:**  
- Cramér’s V for confounding assessment  
- PCA + UMAP visualization  
- PERMANOVA (bray) to validate correction  

### 5. `05_taxon_ecology.R`

**Purpose:** Assess ecological diversity (alpha, beta) on batch-corrected taxonomic profiles.

**Inputs:**  
- `hellinger_transformed_combat.csv`  
- `Aligned_metadata_Taxonomy.csv`

**Outputs:**  
- Alpha diversity statistics and plots (Shannon, Simpson)  
- Beta diversity: NMDS, ANOSIM, betadisper  
- Stressplots and Cliff’s Delta for pairwise contrasts

**Figures:**  
- `Alpha_Diversity_Panel.png`, `NMDS_Only.png`, `ANOSIM_species.tiff`  
- `beta_diversity_summary.csv`, `Stressplot_Only.png`

### 6. `06_arg_vf_normalization.py`

**Purpose:** Normalize CARD and VFDB gene hits per sample and apply Hellinger transformation.

**Inputs:**  
- `Merged_CARD.csv` and `Merged_VFDB.csv`

**Outputs:**  
- `CRC_CARD_expression_matrix_hellinger.csv`  
- `CRC_VFDB_expression_matrix_hellinger.csv`

**Steps:**  
- Presence–absence counting  
- Relative abundance calculation  
- Gene-level aggregation  
- Transformation for ecological analysis

### 7. `07_arg_vf_batch_correction.R`

**Purpose:** Apply batch correction to merged ARG + VF matrices using `limma::removeBatchEffect`.

**Inputs:**  
- Merged Hellinger matrices from step 6  
- Metadata (`CRC_vf_arg_meta_main.csv`)

**Outputs:**  
- `Combined_ARG_VFDB_batch_corrected.csv`  
- Final matrix with negatives zeroed: `combined_ARG_VFDB_final_DATA.csv`  
- UMAP and PERMANOVA diagnostics before/after correction

**Batch variables:** Project, Center, Instrument

### 8. `08_arg_vf_ecology.R`

**Purpose:** Run ecological analyses (alpha/beta diversity, dispersion, NMDS, Procrustes) on batch-corrected ARG/VF data.

**Inputs:**  
- `combined_ARG_VFDB_final_DATA.csv`  
- `Metadata_Aligned_to_CountMatrix.csv`

**Outputs:**  
- Alpha diversity boxplots + LM model predictions  
- `adonis_ARG_adjusted.csv`, `adonis_VF_adjusted.csv`  
- `NMDS_ARG_plot.png`, `NMDS_VF_plot.png`  
- Procrustes correlation: `procrustes_statistics.csv`  
- Null model tests: `oecosimu_ARG_results.txt`, `oecosimu_VF_results.txt`

## Dependencies

### Python (Scripts 1, 2, 3, 6):
- `pandas`, `seaborn`, `matplotlib`, `statsmodels`, `numpy`

### R (Scripts 4, 5, 7, 8):
- `vegan`, `sva`, `limma`, `ggplot2`, `ggpubr`, `patchwork`, `effsize`, `ggeffects`, `umap`, `car`, `broom`, `dplyr`

## Outputs Summary

| Script | Output Type | Example Files |
|--------|-------------|----------------|
| 01     | QC Plots     | `CRC_MAGS_completeness_distribution.png` |
| 02     | Matrix       | `Species_expression_matrix.csv` |
| 03     | Subset       | `Species_subset_matrix.csv` |
| 04     | Corrected Data + PCA/UMAP | `combat_corrected_matrix.csv`, `UMAP_Before_After_Correction.tiff` |
| 05     | Ecology Stats | `adonis_adjusted.csv`, `beta_diversity_summary.csv` |
| 06     | ARG/VF matrices | `CRC_CARD_expression_matrix_hellinger.csv` |
| 07     | Batch-Corrected ARG/VF | `combined_ARG_VFDB_final_DATA.csv` |
| 08     | Ecology Visuals | `NMDS_ARG_plot.png`, `ARG_Beta_dispersion.png` |