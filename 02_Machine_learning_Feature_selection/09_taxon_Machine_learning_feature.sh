#!/bin/bash
#SBATCH --job-name=rf_biomarker_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/rf_pipeline_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/rf_pipeline_%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 <<'EOF'
import os
import pandas as pd
import numpy as np
from collections import Counter
import joblib
import matplotlib.pyplot as plt

from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_validate, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, cohen_kappa_score, matthews_corrcoef

from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

# === CONFIGURATION ===
base_path        = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal"
matrix_file      = os.path.join(base_path, "hellinger_subset_combat_afterbatch_matrix.csv")
metadata_file    = os.path.join(base_path, "Aligned_metadata_Taxonomy.csv")

benchmark_csv    = os.path.join(base_path, "clf_benchmark_metrics.csv")
gridsearch_csv   = os.path.join(base_path, "rf_gridsearch_results.csv")
featimp_csv      = os.path.join(base_path, "rf_feature_importances.csv")
stability_csv    = os.path.join(base_path, "bootstrap_stability_results.csv")
robust_csv       = os.path.join(base_path, "robust_biomarkers.csv")
plot_file        = os.path.join(base_path, "top_robust_biomarkers.png")

target_col       = "Group"
drop_meta_cols   = ["Instrument","Project","Center_Name","Continent",
                    "Country","Age","BMI","Sex"]

# === LOAD & MERGE ===
print("📥 Loading data…")
X_hell = pd.read_csv(matrix_file, index_col=0).T
X_hell.index.name = "Sample_ID"

meta  = pd.read_csv(metadata_file)
# if first column is an unnamed index, turn it into Sample_ID then drop
if "Unnamed: 0" in meta.columns:
    meta["Sample_ID"] = meta["Unnamed: 0"].astype(str)
    meta.drop(columns=["Unnamed: 0"], inplace=True)
# ensure there's exactly one Sample_ID column
assert meta.columns.tolist().count("Sample_ID") == 1, "Duplicate Sample_ID columns in metadata"

# merge
df = pd.merge(
    X_hell.reset_index(), 
    meta, 
    on="Sample_ID", 
    how="inner"
).dropna(subset=[target_col])
print(f"✅ Merged data shape: {df.shape}")

# === PREPARE FEATURES & TARGET ===
# Recover relative abundances
print("🔄 Recovering relative abundances by squaring Hellinger data…")
X_rel = df.loc[:, X_hell.columns] ** 2

# encode target
le    = LabelEncoder()
y     = le.fit_transform(df[target_col])

# Final feature matrix
X     = X_rel.select_dtypes(include=[np.number])
feature_names = X.columns.tolist()

# Stratified CV splitter
cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# === BENCHMARKING CLASSIFIERS ===
print("🔍 Benchmarking classifiers…")
scoring = {
    "accuracy":    "accuracy",
    "f1_weighted": "f1_weighted",
    "kappa":       make_scorer(cohen_kappa_score),
    "mcc":         make_scorer(matthews_corrcoef)
}

classifiers = {
    "LogisticRegression": Pipeline([
        ("scale",   StandardScaler()),
        ("clf",     LogisticRegression(class_weight="balanced", max_iter=1_000, random_state=42))
    ]),
    "SVM_simple": Pipeline([
        ("scale",   StandardScaler()),
        ("clf",     SVC(kernel="linear", class_weight="balanced", probability=False, random_state=42))
    ]),
    "KNN": Pipeline([
        ("scale",   StandardScaler()),
        ("clf",     KNeighborsClassifier())
    ]),
    "GradientBoosting": GradientBoostingClassifier(random_state=42),
    "RandomForest":      RandomForestClassifier(class_weight="balanced", random_state=42)
}

rows = []
for name, model in classifiers.items():
    print(f" • {name}")
    cvres = cross_validate(
        model, X, y,
        cv=cv,
        scoring=scoring,
        n_jobs=-1,
        return_train_score=False
    )
    row = {"Classifier": name}
    for metric in scoring:
        row[f"{metric}_mean"] = cvres[f"test_{metric}"].mean()
        row[f"{metric}_std"]  = cvres[f"test_{metric}"].std()
    rows.append(row)

pd.DataFrame(rows).to_csv(benchmark_csv, index=False)
print(f"✅ Benchmark results → {benchmark_csv}")

# === RANDOM FOREST GRID SEARCH ===
print("⚙️  GridSearchCV on RandomForest (F1_weighted)…")
rf = RandomForestClassifier(class_weight="balanced", random_state=42)
param_grid = {
    "n_estimators":       [100, 200, 300],
    "max_depth":          [5, 10, 20],
    "max_features":       ["sqrt", "log2"],
    "min_samples_split":  [2, 5, 10]
}
gs = GridSearchCV(
    rf, param_grid,
    cv=cv,
    scoring="f1_weighted",
    n_jobs=-1,
    verbose=1
)
gs.fit(X, y)

best_rf    = gs.best_estimator_
best_score = gs.best_score_
best_params= gs.best_params_
print(f"   • Best RF params: {best_params}")
print(f"   • Best RF weighted-F1: {best_score:.4f}")

from sklearn.metrics import classification_report, confusion_matrix

# Predict on full data
y_pred = best_rf.predict(X)

# === CLASSIFICATION REPORT ===
print("🧾 Generating classification report…")
report = classification_report(y, y_pred, target_names=le.classes_, output_dict=True)
report_df = pd.DataFrame(report).transpose()
report_file = os.path.join(base_path, "rf_classification_report.csv")
report_df.to_csv(report_file)
print(f"✅ Classification report → {report_file}")

# === CONFUSION MATRIX ===
print("🔄 Computing confusion matrix…")
cm = confusion_matrix(y, y_pred)
cm_df = pd.DataFrame(cm, index=le.classes_, columns=le.classes_)
cm_file = os.path.join(base_path, "rf_confusion_matrix.csv")
cm_df.to_csv(cm_file)
print(f"✅ Confusion matrix → {cm_file}")


# === SAVE THE BEST RF MODEL ===
model_file = os.path.join(base_path, "rf_best_model.joblib")
joblib.dump(best_rf, model_file)
print(f"✅ Saved best RF model to: {model_file}")



# Save detailed grid results
res = pd.DataFrame(gs.cv_results_)[[
    "param_n_estimators","param_max_depth",
    "param_max_features","param_min_samples_split",
    "mean_test_score","std_test_score","rank_test_score"
]]
res.to_csv(gridsearch_csv, index=False)
print(f"✅ Grid search → {gridsearch_csv}")

# === FEATURE IMPORTANCE (IMPURITY‐BASED) ===
print("🔁 Extracting impurity-based feature_importances_…")
imp = best_rf.feature_importances_
imp_df = pd.DataFrame({
    "Feature": feature_names,
    "Importance": imp
}).sort_values("Importance", ascending=False)
imp_df.to_csv(featimp_csv, index=False)
print(f"✅ Feature importances → {featimp_csv}")

# === BOOTSTRAP STABILITY (TOP‐20 FREQ) ===
print("🔁 Computing bootstrap stability (top-20 features)…")
top_feats = []
for i in range(100):
    # resample
    idx = np.random.RandomState(i).choice(len(X), size=len(X), replace=True)
    X_bs, y_bs = X.iloc[idx], y[idx]
    rf_bs = RandomForestClassifier(**best_params, class_weight="balanced", random_state=i)
    rf_bs.fit(X_bs, y_bs)
    imp_bs = rf_bs.feature_importances_
    top20 = np.argsort(imp_bs)[::-1][:20]
    top_feats.extend([feature_names[j] for j in top20])

cnt = Counter(top_feats)
stab_df = (
    pd.DataFrame.from_dict(cnt, orient="index", columns=["Top20_Count"])
      .assign(Stability=lambda d: d["Top20_Count"] / 100 * 100)
      .rename_axis("Feature")
      .reset_index()
)
stab_df.to_csv(stability_csv, index=False)
print(f"✅ Bootstrap stability → {stability_csv}")

# === MERGE & RANK ROBUST BIOMARKERS ===
merged = (
    imp_df.set_index("Feature")
           .join(stab_df.set_index("Feature"), how="outer")
           .fillna(0)
           .assign(
               Consistent=lambda d: ((d["Importance"]>0)&(d["Top20_Count"]>0)).astype(int)
           )
)
# average rank across Importance and Stability
r_imp   = merged["Importance"].rank(ascending=False)
r_stab  = merged["Stability"].rank(ascending=False)
merged["Avg_Rank"] = (r_imp + r_stab) / 2
merged = merged.sort_values("Avg_Rank")
merged.reset_index().to_csv(robust_csv, index=False)
print(f"✅ Robust biomarkers → {robust_csv}")

# === PLOT TOP 20 ===
print("📊 Plotting top 20 robust biomarkers…")
top20 = merged.head(20).sort_values("Importance", ascending=True)
fig, ax = plt.subplots(figsize=(10,6))
bars = ax.barh(top20.index, top20["Importance"], edgecolor="black", linewidth=1.5)
ax.set_xlabel("Impurity-based Feature Importance")
ax.set_ylabel("Feature")
ax.set_title("Top 20 Robust Biomarkers")
for spine in ax.spines.values():
    spine.set_linewidth(2.0)
plt.tight_layout()
plt.savefig(plot_file, dpi=600)
plt.close()
print(f"✅ Plot saved to: {plot_file}")

EOF
