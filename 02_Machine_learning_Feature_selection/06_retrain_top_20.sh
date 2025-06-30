#!/bin/bash
#SBATCH --job-name=rf_biomarker_panel
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=10:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/rf_biomarker_panel_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/rf_biomarker_panel_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 - <<EOF

import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import make_scorer, cohen_kappa_score, matthews_corrcoef

# === CONFIGURATION ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal"
data_file = os.path.join(base_path, "combined_ARG_VFDB_final_DATA.csv")
meta_file = os.path.join(base_path, "Metadata_Aligned_VF_ARG_CountMatrix.csv")
biomarker_file = os.path.join(base_path, "robust_biomarkers.csv")
output_metrics = os.path.join(base_path, "rf_reduced_biomarker_performance.csv")
top_n = 20  # adjust as needed
target = "Group"
drop_cols = ["Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex", "Sample_ID"]

# === Load and merge data ===
print("Loading data and top biomarkers...")
gene_df = pd.read_csv(data_file, index_col=0).T
gene_df.index.name = "Sample_ID"
gene_df.reset_index(inplace=True)

meta_df = pd.read_csv(meta_file)
meta_df["Sample_ID"] = meta_df["Sample_ID"].astype(str)

df = pd.merge(gene_df, meta_df, on="Sample_ID", how="inner")
df = df.dropna(subset=[target])

# === Encode target ===
label_encoder = LabelEncoder()
df[target] = label_encoder.fit_transform(df[target])

# === Load top biomarkers ===
biomarker_df = pd.read_csv(biomarker_file)
top_features = biomarker_df.head(top_n)["Feature"].tolist()

# === Prepare features ===
X = df[top_features].copy()
y = df[target]

# === Define custom scorers ===
def kappa_scorer(y_true, y_pred):
    return cohen_kappa_score(y_true, y_pred)

def mcc_scorer(y_true, y_pred):
    return matthews_corrcoef(y_true, y_pred)

scoring = {
    'accuracy': 'accuracy',
    'f1_macro': 'f1_macro',
    'kappa': make_scorer(kappa_scorer),
    'mcc': make_scorer(mcc_scorer)
}

# === Evaluate reduced-feature RF model ===
print(" Evaluating RF model on top {} biomarkers...".format(top_n))
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
model = RandomForestClassifier(n_estimators=300, random_state=42, class_weight="balanced")

scores = cross_validate(model, X, y, cv=cv, scoring=scoring, return_train_score=False, n_jobs=-1)

results = {
    "Top N Biomarkers": top_n,
    "Accuracy": np.mean(scores["test_accuracy"]),
    "F1 Macro": np.mean(scores["test_f1_macro"]),
    "Kappa": np.mean(scores["test_kappa"]),
    "MCC": np.mean(scores["test_mcc"])
}

pd.DataFrame([results]).to_csv(output_metrics, index=False)
print("Results saved to:", output_metrics)

EOF
