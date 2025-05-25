
# CRC Microbiome Batch Correction and Ecology Module

This repository contains scripts to analyze colorectal cancer (CRC) metagenomic data with a focus on batch correction, feature subsetting, ecological diversity analyses, and advanced ordination/statistical methods.

---

## Scripts Overview

### 01_Batch_correction.R
**Purpose:**  
Prepares, normalizes, and batch-corrects species-level abundance matrices across cohorts. Visualizes technical covariates' effects via PCA and UMAP before and after ComBat correction.

**Main Functions:**
- Aligns metadata and feature matrix
- Applies ComBat correction using sequencing project and center
- Performs PCA and UMAP to assess batch effects
- Saves corrected matrices for downstream use

---

### 02_Subset_ARG-VF_Species.py
**Purpose:**  
Filters the global microbiome abundance matrix to retain only species annotated in CARD and VFDB. Transposes and aligns this matrix with metadata.

**Main Functions:**
- Filters species by matching names in pathogen annotation lists
- Aligns metadata and microbial expression matrices
- Saves aligned expression and metadata for downstream analyses

---

### 03_Ecology.R
**Purpose:**  
Performs ecological analyses (alpha & beta diversity), ordination (NMDS, CCA), and statistical testing (PERMANOVA, ANOSIM, PERMDISP) on functionally enriched taxa.

**Main Functions:**
- Computes Shannon, Simpson, and richness indices
- Conducts Kruskal-Wallis and Dunn's post hoc tests
- Calculates Bray–Curtis dissimilarity, NMDS, PERMANOVA, ANOSIM
- Conducts adjusted and partial Canonical Correspondence Analysis (CCA)
- Generates publication-quality figures (boxplots, NMDS, CCA, PERMDISP)

---

##  Software Dependencies

### R packages (install via CRAN):
- vegan
- sva
- ggplot2
- ggpubr
- dplyr
- lattice
- umap
- car
- MuMIn
- FSA
- pairwiseAdonis

### Python packages (install via pip):
- pandas

---

## Usage Instructions

1. Run `01_Batch_correction.R` to clean and correct the raw species matrix.  
2. Run `02_Subset_ARG-VF_Species.py` to subset the functionally annotated taxa.  
3. Run `03_Ecology.R` to perform statistical, ordination, and diversity analyses.  

---

## Execution Tree

```bash
crc_microbiome_pipeline/
├── 01_Batch_correction.R
├── 02_Subset_ARG-VF_Species.py
└── 03_Ecology.R
```

---

## Notes

- Ensure all input files are present in the working directory.
- Modify file paths in scripts if using a different directory structure.
- Figures will be saved as high-resolution `.tiff` files suitable for publication.
