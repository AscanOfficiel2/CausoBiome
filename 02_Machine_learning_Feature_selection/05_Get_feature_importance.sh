#!/bin/bash
#SBATCH --job-name=rf_feature_robustness
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/rf_biomarker_robustness_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/rf_biomarker_robustness_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 - <<EOF

import pandas as pd
import numpy as np
import os
import joblib
from sklearn.inspection import permutation_importance
from sklearn.utils import resample
from collections import Counter
import matplotlib.pyplot as plt

# === CONFIGURATION ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal"
data_file = os.path.join(base_path, "combined_ARG_VFDB_final_DATA.csv")
meta_file = os.path.join(base_path, "Metadata_Aligned_VF_ARG_CountMatrix.csv")
model_file = os.path.join(base_path, "best_rf_model_microbial.pkl")
output_csv = os.path.join(base_path, "robust_biomarkers.csv")
plot_file = os.path.join(base_path, "top_robust_biomarkers.png")
target = "Group"
drop_cols = ["Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex", "Sample_ID"]

# === Load data ===
print(" Loading data and model...")
X_raw = pd.read_csv(data_file, index_col=0).T
X_raw.index.name = "Sample_ID"
X_raw.reset_index(inplace=True)
meta = pd.read_csv(meta_file)
meta["Sample_ID"] = meta["Sample_ID"].astype(str)
df = pd.merge(X_raw, meta, on="Sample_ID", how="inner")
df = df.dropna(subset=[target])

# === Prepare features ===
from sklearn.preprocessing import LabelEncoder
le = LabelEncoder()
df[target] = le.fit_transform(df[target])
X = df.drop(columns=drop_cols + [target], errors="ignore").select_dtypes(include=[np.number])
y = df[target]
feature_names = X.columns

model = joblib.load(model_file)

# === Permutation importance ===
print(" Computing permutation importance...")
perm = permutation_importance(model, X, y, n_repeats=30, random_state=42, n_jobs=-1)
perm_df = pd.DataFrame({
    "Feature": feature_names,
    "Importance_Mean": perm.importances_mean,
    "Importance_STD": perm.importances_std
}).sort_values(by="Importance_Mean", ascending=False)

# === Bootstrap stability ===
print(" Computing bootstrap stability...")
top_features_all = []
for i in range(100):
    X_bs, y_bs = resample(X, y, replace=True, random_state=i)
    model_bs = model.fit(X_bs, y_bs)
    importances = model_bs.feature_importances_
    top_idx = np.argsort(importances)[::-1][:20]
    top_features = [feature_names[i] for i in top_idx]
    top_features_all.extend(top_features)

# Count frequency of top features
stability_counts = Counter(top_features_all)
stability_df = pd.DataFrame.from_dict(stability_counts, orient="index", columns=["Top20_Freq"])
stability_df["Stability (%)"] = 100 * stability_df["Top20_Freq"] / 100
stability_df.index.name = "Feature"

# === Save individual components ===
perm_csv = os.path.join(base_path, "permutation_importance_results.csv")
stab_csv = os.path.join(base_path, "bootstrap_stability_results.csv")

perm_df.to_csv(perm_csv, index=False)
print(" Saved permutation importance to:", perm_csv)

stability_df.to_csv(stab_csv)
print(" Saved bootstrap stability to:", stab_csv)


# === Merge and output ===
merged_df = perm_df.set_index("Feature").join(stability_df, how="outer").fillna(0)
merged_df["Consistency"] = ((merged_df["Importance_Mean"] > 0) & (merged_df["Top20_Freq"] > 0)).astype(int)
merged_df["Avg_Rank"] = merged_df[["Importance_Mean", "Stability (%)"]].rank(ascending=False).mean(axis=1)
merged_df = merged_df.sort_values("Avg_Rank")

merged_df.to_csv(output_csv)
print(" Saved robust biomarker rankings to:", output_csv)

# === Plot top 20 ===
top_plot = merged_df.head(20).sort_values("Importance_Mean", ascending=True)
plt.figure(figsize=(10, 6))

bars = plt.barh(
    top_plot.index,
    top_plot["Importance_Mean"],
    color="skyblue",
    edgecolor="black",        # Add black border around bars
    linewidth=1.8             # Thickness of bar borders
)

plt.xlabel("Mean Permutation Importance", fontsize=12)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)

# Thicken all axes lines
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2.0)

plt.tight_layout()
plt.savefig(plot_file, dpi=600)
print(" Saved top biomarker plot to:", plot_file)


EOF
