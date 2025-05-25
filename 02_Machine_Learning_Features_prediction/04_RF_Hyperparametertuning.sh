#!/bin/bash
#SBATCH --job-name=rf_classifier_env
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=35:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Colorectal/out/slurm-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Colorectal/err/slurm-%j.err

# Load Python module
module load Python/3.11.5-GCCcore-13.2.0

# (Optional) Activate your own Python environment
source ~/python_env/bin/activate

# Run Python code
python3 - <<EOF

import pandas as pd
import numpy as np
import joblib
import os
import matplotlib
matplotlib.use('Agg')  # For headless plotting (no X11)
import matplotlib.pyplot as plt

from sklearn.preprocessing import LabelEncoder, label_binarize
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, train_test_split, cross_val_score
from sklearn.metrics import (
    classification_report,
    confusion_matrix,
    ConfusionMatrixDisplay,
    roc_curve,
    auc,
    accuracy_score
)


# --- Paths ---
file_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Colorectal"
env_meta_data_file_path = os.path.join(file_path, "Aligned_metadata.csv")
env_data_file_path = os.path.join(file_path, "Species_matrix_transposed.csv")



print("\U0001f504 Loading data...")
env_data = pd.read_csv(env_data_file_path)
env_meta_data = pd.read_csv(env_meta_data_file_path)

def combine_data(data, meta_data, columns_to_drop=[]):
    combine_df = []
    meta_data = meta_data.drop(columns=columns_to_drop)
    for i in range(len(meta_data)):
        row = meta_data.iloc[i:i+1]
        sample_id = row['Sample_ID'].values[0]
        if sample_id not in data.columns:
            continue
        df = data[["Species", sample_id]].rename(columns={sample_id: 'relative_ab'})
        repeated_row = pd.concat([row]*len(df), ignore_index=True)
        df_id = pd.concat([df, repeated_row], axis=1)
        df_id = df_id[df_id['relative_ab'] != 0.0]
        combine_df.append(df_id)
    combine_df = pd.concat(combine_df)
    combine_df = combine_df.drop(columns=['Sample_ID'])
    return combine_df

print("\U0001f6e0\ufe0f Processing data...")
columns_to_drop = ["Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex"]
combined_df = combine_data(env_data, env_meta_data, columns_to_drop=columns_to_drop)

# Data
target = "Group"
X = combined_df.drop(columns=[target])
y = combined_df[target]

# Encode object-type columns
label_encoder = LabelEncoder()
for col in X.columns:
    if col != 'Species' and X[col].dtype == 'object':
        X[col] = label_encoder.fit_transform(X[col])

# One-hot encode 'Species'
if 'Species' in X.columns:
    X = pd.get_dummies(X, columns=['Species'])

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, test_size=0.2, random_state=42
)

print("\U0001f332 Running Random Forest GridSearchCV...")
param_grid = {
    'n_estimators': [100, 200, 300],
    'max_depth': [5, 10, None],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 5],
    'max_features': ['sqrt', 'log2']
}


rf = RandomForestClassifier(random_state=42, class_weight="balanced")
grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=10,
                           n_jobs=-1, scoring='accuracy', verbose=2)
grid_search.fit(X_train, y_train)
best_rf = grid_search.best_estimator_

print("\n\u2705 Best Parameters Found:")
print(grid_search.best_params_)

y_pred = best_rf.predict(X_test)
print("\n\U0001f4c8 Classification Report:")
print(classification_report(y_test, y_pred))

# Confusion Matrix
conf_matrix = pd.DataFrame(confusion_matrix(y_test, y_pred),
                           index=np.unique(y_test),
                           columns=np.unique(y_test))
print("\n\U0001f9ee Confusion Matrix:")
print(conf_matrix)

# Save model
model_path = os.path.join(file_path, "best_RF_MODEL_Without_Confounders.pkl")
joblib.dump(best_rf, model_path)
print(f"\n\U0001f4be Model saved to: {model_path}")

# ---------------------- PLOTTING ----------------------

# Confusion Matrix Plot
disp = ConfusionMatrixDisplay(confusion_matrix=conf_matrix.values,
                              display_labels=conf_matrix.columns)
fig, ax = plt.subplots(figsize=(8, 6))
disp.plot(ax=ax, cmap="Blues", xticks_rotation=45)
#plt.title("Confusion Matrix")
plt.tight_layout()
plt.savefig(os.path.join(file_path, "confusion_matrix.png"), dpi=600)
plt.close()

# ROC Curve (One-vs-Rest)
y_test_bin = label_binarize(y_test, classes=np.unique(y))
y_score = best_rf.predict_proba(X_test)

fpr = dict()
tpr = dict()
roc_auc = dict()
n_classes = y_test_bin.shape[1]

plt.figure(figsize=(7, 6))
for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test_bin[:, i], y_score[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])
    plt.plot(fpr[i], tpr[i], lw=2.5,
             label=f"Class {np.unique(y)[i]} (AUC = {roc_auc[i]:.2f})")

plt.plot([0, 1], [0, 1], "k--", lw=3)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
#plt.title("Multiclass ROC Curve")
plt.legend(loc="lower right")
plt.tight_layout()
plt.savefig(os.path.join(file_path, "roc_curve_without_confounders.png"), dpi=600)
plt.close()


# Overfitting Check: Train vs Test Accuracy
train_preds = best_rf.predict(X_train)
test_preds = best_rf.predict(X_test)
train_acc = accuracy_score(y_train, train_preds)
test_acc = accuracy_score(y_test, test_preds)
print("\n Overfitting Check:")
print(f"Train Accuracy: {train_acc:.4f}")
print(f"Test Accuracy:  {test_acc:.4f}")

# Cross-Validation Check
print("\n Running 10-fold cross-validation on entire data...")
cv_scores = cross_val_score(best_rf, X, y, cv=10, scoring='accuracy')
print(f"Cross-validated Accuracy: {cv_scores.mean():.4f} ± {cv_scores.std():.4f}")



EOF
