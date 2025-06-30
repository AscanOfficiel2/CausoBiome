---
title: "Module 03 — Causal Inference Analysis"
output:
  pdf_document: default
  html_document:
    df_print: paged
date: "2025-06-30"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Module 03 — Causal Inference Analysis

**Part of the CausoBiome Pipeline**  
Lead Developer: AbdulAziz Ascandari  
Affiliation: Mohammed VI Polytechnic University  
Version: 1.0

## Objective

This module identifies potentially causal microbial genes (ARGs/VFs) contributing to colorectal cancer (CRC) progression. It integrates feature annotation, functional distribution analysis, effect estimation using Double Machine Learning (DML), bootstrap robustness, interaction effects, and sensitivity analysis via E-values.

---

## Prerequisite Notes

**Prepare the annotation files** for CARD and VFDB by extracting the following columns:

- **CARD**: `Sample_ID`, `Gene`, `Matched_Species`, `AMR Gene Family`, `Resistance_Mechanism`, `Antibiotics_Class`  
- **VFDB**: `Sample_ID`, `Gene`, `Matched_Species`, `Protein`, `Resistance_Mechanism`, `Class_Mechanism`

These files are required before executing the causal analysis pipeline.

---

## Scripts Included

### `01_Complex_heatmap_arg_vf.R`
- **Purpose**: Generate z-score heatmap of consistent ARG/VF biomarkers with functional and taxonomic annotations.
- **Inputs**:  
  - `combined_ARG_VFDB_final_DATA.csv`  
  - `Metadata_Aligned_VF_ARG_CountMatrix.csv`  
  - `Filtered_biomarker_annotated.csv`
- **Output**: `heatmap_consistent_genes.png`

---

### `02_causal_analysis.py`
- **Purpose**: End-to-end causal analysis and visualization of microbial gene effects on CRC stage.
- **Includes**:
  - Annotation and merging of robust biomarkers
  - Functional pie chart & group differences
  - Post-hoc Pearson residuals for functional shift
  - Gene-wise log2FC + statistical testing (Mann–Whitney + FDR)
  - Signature score analysis (mean & PCA-based)
  - Causal effect estimation (DML with Random Forest)
  - Bootstrap robustness (confidence intervals)
  - E-value sensitivity analysis
  - Interaction network effects and classification (synergistic, antagonistic)
- **Outputs**:
  - `robust_biomarkers_annotated.csv`, `Filtered_biomarker_annotated.csv`
  - `heatmap_log2fc_consistent_biomarkers.png`
  - `DML_Causal_Effects_Per_Species.csv`, `Sensitivity_Analysis_Evalues.csv`
  - `Significant_Microbial_Interactions_Final.png`, `Gene_Interaction_Network_Metrics.csv`

---

## Output Summary

| Analysis Step                     | Key Outputs                                               |
|----------------------------------|------------------------------------------------------------|
| Functional Annotation            | `Filtered_biomarker_annotated.csv`                         |
| Expression Heatmap               | `heatmap_consistent_genes.png`                             |
| Functional Shift Analysis        | `Functional_Group_Differences_Across_Stages.csv`           |
| log2FC + Stats                   | `statistical_comparison_log2fc.csv`, `heatmap_log2fc_...` |
| Signature Score Profiling        | `signature_scores.csv`, boxplots                          |
| Causal Estimation (ATE)          | `DML_Causal_Effects_Per_Species.csv`                      |
| Bootstrap Robustness             | `Bootstrap_Causal_Effects_CI.csv`, `Bootstrap_Robustness_Check.png` |
| E-value Sensitivity              | `Sensitivity_Analysis_Evalues.csv`, evalue plot           |
| Interaction Effects              | `Microbial_Interaction_Effects_Corrected.csv`, `Significant_Microbial_Interactions_Final.png` |
| Interaction Network & Metrics    | `Gene_Interaction_Network_Metrics.csv`                    |

---

*End of Module 03 README*