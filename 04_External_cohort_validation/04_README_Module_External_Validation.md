---
title: "Module 04 — External Cohort Validation"
output:
  pdf_document: default
  html_document:
    df_print: paged
date: "2025-06-30"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Module 04 — External Cohort Validation

**Part of the CausoBiome Pipeline**  
Lead Developer: AbdulAziz Ascandari  
Affiliation: Mohammed VI Polytechnic University  
Version: 1.0

## Objective

This module validates microbial biomarkers in an independent metagenomic dataset (PRJEB10878). It evaluates:
- Consistency of CRC vs. control trends across datasets
- Statistical significance using Mann–Whitney U tests
- Shared directionality and effect size
- Visual agreement of biomarkers using plots and Venn diagrams

---

## Script Overview

### `External_cohort_Validation.ipynb`

**Purpose:**  
Evaluate differential abundance and consistency of microbial biomarkers between internal and external datasets.

**Inputs:**
- `Filtered_biomarker_annotated.csv`  
- Internal: `combined_ARG_VFDB_final_DATA.csv`, `Metadata_Aligned_VF_ARG_CountMatrix.csv`  
- External: `Matrix_PRJEB10878_combined_VF_ARG.csv`, `PRJEB10878_Metadata_Aligned.csv`

**Key Steps:**
1. Recode sample groupings (Case/Control, CRC/Control)
2. Perform Mann–Whitney U tests internally and externally
3. Compute mean differences (CRC - Control) per biomarker
4. Assess trend consistency and plot barplot
5. Venn and barplot of shared significance
6. Annotated scatterplot for consistent biomarkers
7. Save merged results with trends and flags

---

## Output Files

| Output File                            | Description                                        |
|----------------------------------------|----------------------------------------------------|
| `External_MannWhitney_Stats.csv`       | Stats from external CRC vs. control test           |
| `Internal_MannWhitney_Stats.csv`       | Stats from internal case vs. control test          |
| `Trend_Comparison.csv`                 | Mean trend comparison table                        |
| `Barplot_Trend_Comparison.png`         | Barplot of trends for shared biomarkers            |
| `Merged_Significance_Comparison.csv`   | Unified significance table with trend concordance  |
| `Significant_in_Both_Datasets.csv`     | Biomarkers significant in both datasets            |
| `Significance_Comparison_Barplot.png`  | Barplot of overlapping significance labels         |
| `Significance_Comparison_Venn.png`     | Venn diagram of shared significance                |
| `Consistent_Trends.png`                | Annotated scatterplot of consistent biomarker trends |

---

## Visual Highlights

- **Barplot of Trends**: Displays direction of case/control mean differences for both datasets
- **Venn Diagram**: Illustrates overlap of statistically significant biomarkers
- **Annotated Scatterplot**: Shows correlated trends in internal and external datasets with size-encoded magnitude

---

*End of Module 04 README*