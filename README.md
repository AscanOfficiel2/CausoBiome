# CausoBiome: A Multi-Stage Pipeline for Microbiome-Driven Causal Discovery and Intervention Design

**CausoBiome** is a modular, stage-aware pipeline designed to identify, validate, and prioritize microbial species that are functionally and causally linked to disease progression. It is uniquely structured to model **ordinal clinical transitions**, specifically in Colorectal Cancer from *Healthy* to *Adenoma* to *Cancer*, and to uncover microbial contributors driving disease trajectory.CausoBiome builds upon high-quality, genome-resolved metagenomic preprocessing pipelines from the [`genome-resolved-urban-microbiome-biosurveillance`](https://github.com/SuleimanAminu/genome-resolved-urban-microbiome-biosurveillance) repository, specifically Modules `01_Bioinformatics` and `02_Quality_Batch_subsetting`. By extending these upstream workflows, CausoBiome introduces powerful downstream functionality including **ecological profiling**, **machine learning-based classification**, **causal inference via Double Machine Learning (DML)**, and **in silico microbial intervention modeling**. Purpose-built not merely to describe microbiome differences, CausoBiome infers **stage-specific microbial drivers**, models **inter-species interactions**, and simulates **therapeutic perturbations** bridging the gap between microbial discovery and translational action.
Importantly, **CausoBiome is inherently extensible to other diseases with well-defined ordinal stages**. Its modular architecture, causal modeling core, and intervention simulation engine make it adaptable to any pathology where microbiome dynamics evolve progressively across clinical stages.

---

## Pipeline Structure

```text
CausoBiome/
├── Module 1: Batch Correction & Ecology (from processed feature tables)
│   ├── 01_Batch_correction.R
│   ├── 02_Subset_ARG-VF_Species.py
│   └── 03_Ecology.R

├── Module 2: Classification & Feature Importance
│   ├── 04_RF_Hyperparametertuning.sh
│   ├── 05_Benchmarking_ML_models.sh
│   ├── 06_RF_Permutation_Stability_Importance.sh
│   └── 07_Comparison_Script.R

└── Module 3: Ordinal Causal Inference & Intervention Modeling
    └── 08_Ordinal_Forest_Causal_Intervention_Analysis.py
```

---

## Purpose and Novelty

While prior workflows end at taxonomic profiling or classification, **CausoBiome** extends microbiome analytics into the **causal and translational domain**. It enables:

- High-resolution species-level modeling of disease progression 
- Robust batch correction and functional species selection (VFDB, CARD)
- Benchmarking of supervised classifiers for stage prediction
- Causal estimation of microbial effects via Double Machine Learning (DML)
- Network-based modeling of synergistic and antagonistic microbial interactions
- Simulated interventions for therapeutic prioritization based on risk impact

---

##  Key Functional Highlights

| Module | Functionality |
|--------|---------------|
| Batch Correction & Ecology | ComBat correction, Hellinger transform, NMDS/CCA ordination |
| Classification | RF model tuning, model benchmarking, robust feature selection |
| Causal Inference | ATE estimation, bootstrap & E-value robustness, network synergy modeling |
| Intervention Modeling | ARR, RR, OR simulation, composite prioritization score (CIS) |

---

## Pipeline Origin

CausoBiome is an extension of the upstream **genome-resolved-urban-microbiome-biosurveillance** workflow in:

**GitHub**: [genome-resolved-urban-microbiome-biosurveillance](https://github.com/SuleimanAminu/genome-resolved-urban-microbiome-biosurveillance)

- Users should start from the `01_Bioinformatics` module and then proceed to the `02_Quality_batch_subsetting` module by running specifically scripts `04_mag_quality_metrics_analysis.py` and `05_normalize_species_counts.py`
- Then transition into **CausoBiome** starting from `01_Batch_correction.R`

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
- Per-species causal estimates (DML ATE)
- Robustness plots (E-values, bootstraps)
- Microbial interaction networks (weighted, annotated)
- Intervention outcome simulations and priority tiers

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

A Zenodo DOI has been assigned to ensure formal authorship record.

-We kindly request that you do not republish or repackage this methodology before journal publication.



## Citation

If you use **CausoBiome**, please cite the following manuscript:

> **Ascandari, A., Aminu, S., Benhida, R., & Daoud, R.** (2025).  
> *Causal Inference and Species Interaction Networks Reveal Keystone Microbes in Genome-Resolved CRC Progression*.  
> npj Biofilms and Microbiomes (Submitted).

> Zenodo. [https://zenodo.org/records/15511511]([https://doi.org/10.5281/zenodo.15505402](https://zenodo.org/records/15511511))

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15511511.svg)](https://doi.org/10.5281/zenodo.15511511) 

## Submitted Articles Related to the Framework

> **Ascandari, A., Aminu, S., Benhida, R., & Daoud, R.** (2025).  
> *Causal Inference and Species Interaction Networks Reveal Keystone Microbes in Genome-Resolved CRC Progression*.  
> npj Biofilms and Microbiomes (Submitted).

> **Ascandari, A., Aminu, S., Benhida, R., & Daoud, R.** (2025).  
CausoBiome: A Multi-Stage Pipeline for Microbiome-Driven Causal Discovery and Intervention Design.
Bioinformatics (Submitted).

## Contact
For questions, feedback, or collaboration regarding this framework, please reach out:

AbdulAziz Ascandari, PhD Researcher, Department of Chemical and Biochemical Sciences, University Mohammed VI Polytechnic (UM6P), Morocco, abdulaziz.ascandari@um6p.ma

Prof. Rachid Daoud, Group Leader & Supervisor, Department of Chemical and Biochemical Sciences, University Mohammed VI Polytechnic (UM6P), Morocco, rachid.daoud@um6p.ma
