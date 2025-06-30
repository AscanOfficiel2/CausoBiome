#!/bin/bash
#SBATCH --job-name=rf_tuning_microbial
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/rf_tuning_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/rf_tuning_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 - <<EOF

import pandas as pd
import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import make_scorer, cohen_kappa_score, matthews_corrcoef, classification_report
import joblib

# === CONFIG ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal"
abundance_file = os.path.join(base_path, "combined_ARG_VFDB_final_DATA.csv")  # genes x samples
metadata_file = os.path.join(base_path, "Metadata_Aligned_VF_ARG_CountMatrix.csv")
output_path = os.path.join(base_path, "rf_microbial_tuning_results.csv")
model_output = os.path.join(base_path, "best_rf_model_microbial.pkl")
target = "Group"
drop_cols = ["Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex", "Sample_ID"]

# === Load and format ===
print("Loading data...")
gene_df = pd.read_csv(abundance_file, index_col=0).T
gene_df.index.name = "Sample_ID"
gene_df.reset_index(inplace=True)

meta_df = pd.read_csv(metadata_file)
meta_df["Sample_ID"] = meta_df["Sample_ID"].astype(str)

df = pd.merge(gene_df, meta_df, on="Sample_ID", how="inner")
df = df.dropna(subset=[target])
print(f"Merged data shape: {df.shape}")

# === Encode target ===
label_encoder = LabelEncoder()
df[target] = label_encoder.fit_transform(df[target])

# === Prepare microbial-only features ===
X = df.drop(columns=[target] + drop_cols, errors="ignore")
X = X.select_dtypes(include=[np.number])
y = df[target]

# === Hyperparameter grid ===
param_grid = {
    'n_estimators': [100, 300, 500],
    'max_depth': [10, 20, 50, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4],
    'max_features': ['sqrt'],
    'class_weight': ['balanced']
}

# === Scorer ===
def kappa_scorer(y_true, y_pred):
    return cohen_kappa_score(y_true, y_pred)

scoring = {
    'accuracy': 'accuracy',
    'f1_macro': 'f1_macro',
    'kappa': make_scorer(kappa_scorer),
    'mcc': make_scorer(matthews_corrcoef)
}

# === Grid Search ===
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
clf = RandomForestClassifier(random_state=42)

print("Running grid search...")
grid = GridSearchCV(
    clf,
    param_grid,
    scoring=scoring,
    refit='kappa',
    cv=cv,
    n_jobs=-1,
    verbose=1
)



grid.fit(X, y)


y_pred = grid.predict(X)

# === Save classification report ===
report_dict = classification_report(y, y_pred, target_names=label_encoder.classes_, output_dict=True)
report_df = pd.DataFrame(report_dict).transpose()
report_path = os.path.join(base_path, "rf_classification_report.csv")
report_df.to_csv(report_path)
print(f"Classification report saved to: {report_path}")

# === Save confusion matrix ===
conf_matrix = confusion_matrix(y, y_pred)
conf_matrix_df = pd.DataFrame(conf_matrix, index=label_encoder.classes_, columns=label_encoder.classes_)
conf_matrix_path = os.path.join(base_path, "rf_confusion_matrix.csv")
conf_matrix_df.to_csv(conf_matrix_path)
print(f"Confusion matrix saved to: {conf_matrix_path}")


# === Save best model and results ===
print(f"Best Parameters Found:\n{grid.best_params_}")
print("\n Classification Report:")
y_pred = grid.predict(X)
print(classification_report(y, y_pred, target_names=label_encoder.classes_))

pd.DataFrame(grid.cv_results_).to_csv(output_path, index=False)
joblib.dump(grid.best_estimator_, model_output)


# === Plot confusion matrix ===
plt.figure(figsize=(6, 5))
sns.heatmap(
    conf_matrix_df,
    annot=True,
    fmt='d',
    cmap='Blues',
    cbar=False,
    linewidths=1.5,
    linecolor='black',
    square=True
)
plt.ylabel("True Label", fontsize=12)
plt.xlabel("Predicted Label", fontsize=12)
plt.title("Confusion Matrix (RF)", fontsize=14)

# Thicken outer spines
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(2)

plt.tight_layout()

# Save plot
conf_matrix_plot = os.path.join(base_path, "rf_confusion_matrix_plot.png")
plt.savefig(conf_matrix_plot, dpi=600)
plt.close()
print(f" Confusion matrix plot saved to: {conf_matrix_plot}")

EOF
