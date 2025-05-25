
# Ordinal Forest & Causal Inference Pipeline for CRC Microbiome Analysis

This script serves as the core engine of a colorectal cancer (CRC) microbiome analysis framework, integrating feature selection, classification, causal inference, and interventional modeling. It identifies microbial species with stage-specific relevance to CRC, quantifies their causal effects using Double Machine Learning, models synergistic and antagonistic microbial interactions, and simulates abundance-based interventions to estimate risk modification.

---

## Script Overview

### 08_Ordinal_Forest_Causal_Intervention_Analysis.py

**Functions:**

1. **Species-level ranking and validation:**
   - Permutation importance
   - Log-loss-based filtering
   - Random Forest feature stability analysis

2. **Ordinal multiclass classification:**
   - Cross-validation
   - AUROC & F1-macro metrics
   - Negative label shuffling control

3. **Causal inference:**
   - Double Machine Learning (DML) for ATE estimation
   - Confidence intervals and forest plots

4. **Robustness checks:**
   - E-value sensitivity analysis
   - 100× bootstrap resampling of ATEs

5. **Interaction modeling:**
   - Synergistic or antagonistic microbial species pairs
   - Visualization via weighted networks

6. **Intervention simulations:**
   - ARR, RR, OR
   - Composite Intervention Score (CIS)
   - High/Moderate/Low priority tiers

---

## Input Files

| File                               | Description                                             |
|------------------------------------|---------------------------------------------------------|
| `Expression_Matrix_Aligned.csv`   | Sample-wise microbiome abundance matrix (post-Hellinger) |
| `Aligned_metadata.csv`            | Metadata with Group, Age, Sex, BMI                      |
| `Top_Important_Species.csv`       | Top microbial species from prior RF modeling            |
| `consistent_species_feature_summary.csv` | Stable features from permutation + bootstrap     |
| `RF_Log_Loss_Per_Species.csv`     | Species ranked by classification log-loss              |

---

##  Required Libraries

Install using pip:

```bash
pip install pandas numpy seaborn matplotlib scikit-learn econml networkx statsmodels
```

---

## Execution Pipeline

```text
08_Ordinal_Forest_Causal_Intervention_Analysis.py
├── Step 1: Feature Selection & Ranking
│   ├─ Filters species by log-loss (≤ 0.1)
│   └─ Outputs: Ordinal_RF_Full_Evaluation_Metrics.csv, Top_Important_Species.csv
│
├── Step 2: Visualization
│   ├─ PCA and t-SNE
│   └─ Boxplots & heatmaps: CRC_Log2FC_Heatmap_Formatted.png
│
├── Step 3: Causal Inference (DML)
│   └─ DML_Causal_Effects_Per_Species.csv, Forest plots
│
├── Step 4: Sensitivity Analysis
│   └─ E-value computation: Sensitivity_Analysis_Evalues.csv
│
├── Step 5: Bootstrap Validation
│   ├─ 100x resampling
│   └─ Outputs: Bootstrap_Causal_Effects_CI.csv
│
├── Step 6: Microbial Interaction Modeling
│   └─ Synergy/antagonism analysis: Microbial_Interaction_Effects_Corrected.csv
│
└── Step 7: Intervention Simulation
    ├─ Simulated perturbation: +1SD / -1SD
    └─ Outputs: EconML_Intervention_Score_Normalized.csv
```

---

## Key Outputs

| Output File                                      | Description                                             |
|--------------------------------------------------|---------------------------------------------------------|
| `Ordinal_RF_Full_Evaluation_Metrics.csv`        | AUROC, F1, accuracy of ordinal classifier              |
| `DML_Causal_Effects_Per_Species.csv`            | Per-species ATE with confidence intervals              |
| `Sensitivity_Analysis_Evalues.csv`              | E-values for confounding robustness                    |
| `Microbial_Interaction_Effects_Corrected.csv`   | Significant synergistic/antagonistic pairs             |
| `EconML_Intervention_Score_Normalized.csv`      | Ranked list of high-priority species for intervention  |

---

##  Methodological Highlights

- Ordinal modeling reflects clinical disease progression (Healthy → Adenoma → Cancer)
- DML enables robust estimation of causal effects
- E-values quantify resistance to unmeasured confounding
- Bootstrap resampling ensures statistical stability
- Network modeling reveals emergent microbial interaction patterns
- Intervention simulations link model predictions to actionable risk metrics

---

## Reproducibility Tips

- Ensure consistent sample IDs across metadata and expression files
- Pre-clean column headers to remove leading/trailing spaces
- Use a fixed `random_state=42` for deterministic output
- Maintain consistent input sources across scripts (esp. `Top_Important_Species.csv`)

---

