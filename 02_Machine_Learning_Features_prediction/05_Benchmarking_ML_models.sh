#!/bin/bash
#SBATCH --job-name=rf_classifier_env
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/slurm-%j.err

# Load Python
module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

# Run Python block
python3 - <<EOF

import pandas as pd
import numpy as np
import os

from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import (
    make_scorer, cohen_kappa_score, matthews_corrcoef
)
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC, LinearSVC
from sklearn.neighbors import KNeighborsClassifier

# === Set paths ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal"
abundance_file = os.path.join(base_path, "Species_matrix_transposed.csv")
metadata_file = os.path.join(base_path, "Aligned_metadata.csv")



# === Load data ===
print(" Loading data...")
env_data = pd.read_csv(abundance_file)
env_meta_data = pd.read_csv(metadata_file)

# === Combine data ===
def combine_data(data, meta_data, columns_to_drop=[]):
    combined_list = []
    meta_data = meta_data.drop(columns=columns_to_drop)
    for i in range(len(meta_data)):
        row = meta_data.iloc[i:i+1]
        sample_id = row['Sample_ID'].values[0]
        if sample_id not in data.columns:
            continue
        df = data[["Species", sample_id]].rename(columns={sample_id: "relative_ab"})
        repeated_row = pd.concat([row]*len(df), ignore_index=True)
        df_id = pd.concat([df, repeated_row], axis=1)
        df_id = df_id[df_id["relative_ab"] != 0.0]
        combined_list.append(df_id)
    combined_df = pd.concat(combined_list, ignore_index=True)
    combined_df = combined_df.drop(columns=["Sample_ID"])
    return combined_df

# === Drop unneeded metadata ===
drop_cols = ["Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex"]
combined_df = combine_data(env_data, env_meta_data, columns_to_drop=drop_cols)
print(f" Combined shape: {combined_df.shape}")
combined_df.to_csv(os.path.join(base_path, "combined_df.csv"), index=False)

# === Prepare features and target ===
target = "Group"
X = combined_df.drop(columns=[target])
y = combined_df[target]

# Encode feature columns
for col in X.columns:
    if X[col].dtype == 'object':
        X[col] = X[col].astype('category').cat.codes

# Encode target
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)

# === Custom scorers ===
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

# === Classifiers ===
models = {
    'Random Forest': RandomForestClassifier(
        random_state=42,
        max_depth=None,
        max_features='sqrt',
        min_samples_leaf=2,
        min_samples_split=15,
        n_estimators=500,
        class_weight='balanced'
    ),
    'Logistic Regression': LogisticRegression(max_iter=5000),
    'SVM': LinearSVC(max_iter=5000, random_state=42),
    'KNN': KNeighborsClassifier(n_neighbors=5),
    'Gradient Boosting': GradientBoostingClassifier()
}


cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# === Run cross-validation ===
results = {}
all_scores = []

for name, model in models.items():
    print(f"\U0001f50d Evaluating: {name}")
    scores = cross_validate(model, X, y_encoded, cv=cv, scoring=scoring, n_jobs=-1, return_estimator=False)
    
    results[name] = {
        'Accuracy': np.mean(scores['test_accuracy']),
        'F1 Macro': np.mean(scores['test_f1_macro']),
        'Kappa': np.mean(scores['test_kappa']),
        'MCC': np.mean(scores['test_mcc'])
    }
    
    for i in range(len(scores['test_accuracy'])):
        all_scores.append({
            'Model': name,
            'Fold': i+1,
            'Accuracy': scores['test_accuracy'][i],
            'F1 Macro': scores['test_f1_macro'][i],
            'Kappa': scores['test_kappa'][i],
            'MCC': scores['test_mcc'][i]
        })

# === Save outputs ===
results_df = pd.DataFrame(results).T.sort_values(by='Kappa', ascending=False)
folds_df = pd.DataFrame(all_scores)

results_df.to_csv(os.path.join(base_path, "CRC_model_comparison_results.csv"))
folds_df.to_csv(os.path.join(base_path, "CRC_model_foldwise_scores.csv"), index=False)

print("\n Mean Scores Per Model (sorted by Kappa):")
print(results_df)

EOF
