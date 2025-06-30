#!/bin/bash
#SBATCH --job-name=ml_benchmark_wide
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/ml_benchmark_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/ml_benchmark_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 <<EOF
import pandas as pd
import numpy as np
import os
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import make_scorer, cohen_kappa_score, matthews_corrcoef
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.neighbors import KNeighborsClassifier

# === CONFIG ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal"
abundance_file = os.path.join(base_path, "combined_ARG_VFDB_final_DATA.csv")
metadata_file = os.path.join(base_path, "Metadata_Aligned_VF_ARG_CountMatrix.csv")
drop_meta_cols = ["Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex", "Sample_ID"]
target = "Group"

# === Load and prepare data ===
print("Loading data...")
gene_df = pd.read_csv(abundance_file, index_col=0).T
gene_df.index.name = "Sample_ID"
gene_df.reset_index(inplace=True)

meta_df = pd.read_csv(metadata_file)
meta_df["Sample_ID"] = meta_df["Sample_ID"].astype(str)

df = pd.merge(gene_df, meta_df, on="Sample_ID", how="inner")
df = df.dropna(subset=[target])
print(f"Merged data shape: {df.shape}")

label_encoder = LabelEncoder()
df[target] = label_encoder.fit_transform(df[target])

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

def prepare_inputs(include_metadata=True):
    X = df.drop(columns=[target])
    if not include_metadata:
        X = X.drop(columns=[col for col in drop_meta_cols if col in X.columns], errors="ignore")
    else:
        X = X.drop(columns=["Sample_ID"], errors="ignore")
        X = X.select_dtypes(include=[np.number])
    y = df[target]
    return X, y

def evaluate_models(X, y, label):
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
    models = {
        'Random Forest': RandomForestClassifier(n_estimators=300, random_state=42, class_weight='balanced'),
        'Logistic Regression': LogisticRegression(max_iter=5000),
        'SVM': LinearSVC(max_iter=5000, random_state=42),
        'KNN': KNeighborsClassifier(n_neighbors=5),
        'Gradient Boosting': GradientBoostingClassifier()
    }
    results = []
    folds = []
    for name, model in models.items():
        print(f"Evaluating: {name} ({label})")
        scores = cross_validate(model, X, y, cv=cv, scoring=scoring, n_jobs=-1)
        results.append({
            "Model": name,
            "Data": label,
            "Accuracy": np.mean(scores["test_accuracy"]),
            "F1 Macro": np.mean(scores["test_f1_macro"]),
            "Kappa": np.mean(scores["test_kappa"]),
            "MCC": np.mean(scores["test_mcc"])
        })
        for i in range(len(scores["test_accuracy"])):
            folds.append({
                "Model": name,
                "Data": label,
                "Fold": i + 1,
                "Accuracy": scores["test_accuracy"][i],
                "F1 Macro": scores["test_f1_macro"][i],
                "Kappa": scores["test_kappa"][i],
                "MCC": scores["test_mcc"][i]
            })
    return results, folds

# Run evaluations
X_bio, y_bio = prepare_inputs(include_metadata=False)
X_all, y_all = prepare_inputs(include_metadata=True)

results_bio, folds_bio = evaluate_models(X_bio, y_bio, "Microbial Only")
results_all, folds_all = evaluate_models(X_all, y_all, "Microbial + Metadata")

# Save summary + foldwise scores
results_df = pd.DataFrame(results_bio + results_all)
results_df.to_csv(os.path.join(base_path, "full_model_benchmark_comparison.csv"), index=False)

pd.DataFrame(folds_bio).to_csv(os.path.join(base_path, "foldwise_scores_microbial_only.csv"), index=False)
pd.DataFrame(folds_all).to_csv(os.path.join(base_path, "foldwise_scores_microbial_plus_metadata.csv"), index=False)

print("\n Benchmark Results (Summary):")
print(results_df)
print("\n Foldwise scores saved to:")
print("→ foldwise_scores_microbial_only.csv")
print("→ foldwise_scores_microbial_plus_metadata.csv")

EOF
