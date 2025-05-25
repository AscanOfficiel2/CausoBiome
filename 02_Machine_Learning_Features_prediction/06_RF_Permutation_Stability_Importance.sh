#!/bin/bash
#SBATCH --job-name=rf_species_only_features
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Colorectal/out/feature_importance-%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Colorectal/err/feature_importance-%j.err

module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 <<EOF
import pandas as pd
import numpy as np
import os
import joblib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import (
    classification_report, confusion_matrix, roc_curve,
    auc, accuracy_score
)
from sklearn.preprocessing import LabelEncoder, label_binarize
from sklearn.inspection import permutation_importance
from sklearn.utils import resample
from collections import Counter

file_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/user/Colorectal"
meta_file = os.path.join(file_path, "Aligned_metadata.csv")
data_file = os.path.join(file_path, "Species_matrix_transposed.csv")

# --- Load and preprocess data ---
env_data = pd.read_csv(data_file)
env_meta = pd.read_csv(meta_file)

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

columns_to_drop = ["Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex"]
combined_df = combine_data(env_data, env_meta, columns_to_drop)

# Target and features
target = "Group"
X = combined_df.drop(columns=[target])
y = combined_df[target]

# Encode categorical variables
label_encoder = LabelEncoder()
for col in X.columns:
    if col != 'Species' and X[col].dtype == 'object':
        X[col] = label_encoder.fit_transform(X[col])
if 'Species' in X.columns:
    X = pd.get_dummies(X, columns=['Species'])

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, test_size=0.2, random_state=42
)

# --- Random Forest Training ---
rf = RandomForestClassifier(
    n_estimators=300, max_depth=None, max_features='sqrt',
    min_samples_leaf=1, min_samples_split=10,
    class_weight='balanced', random_state=42
)
rf.fit(X_train, y_train)
y_pred = rf.predict(X_test)
y_score = rf.predict_proba(X_test)
#joblib.dump(rf, os.path.join(file_path, "final_species_only_rf_model.pkl"))

# --- Permutation Importance ---
perm = permutation_importance(rf, X_test, y_test, n_repeats=30, random_state=42, n_jobs=-1)
perm_df = pd.DataFrame({
    'Feature': X_test.columns,
    'Importance_Mean': perm.importances_mean,
    'Importance_STD': perm.importances_std
}).sort_values(by='Importance_Mean', ascending=False)
perm_df.to_csv(os.path.join(file_path, "permutation_importance_species_only.csv"), index=False)

# --- Bootstrap Stability ---
bootstrap_iterations = 100
top_features_all = []
for i in range(bootstrap_iterations):
    X_bs, y_bs = resample(X_train, y_train, replace=True, random_state=i)
    rf_bs = RandomForestClassifier(
        n_estimators=300, max_depth=None, max_features='sqrt',
        min_samples_leaf=1, min_samples_split=10,
        class_weight='balanced', random_state=i
    )
    rf_bs.fit(X_bs, y_bs)
    importances = rf_bs.feature_importances_
    top_indices = np.argsort(importances)[::-1][:20]
    top_features = [X.columns[i] for i in top_indices]
    top_features_all.extend(top_features)

stability_counts = Counter(top_features_all)
stability_df = pd.DataFrame.from_dict(stability_counts, orient='index', columns=['Top20_Freq'])
stability_df['Stability (%)'] = 100 * stability_df['Top20_Freq'] / bootstrap_iterations
stability_df = stability_df.sort_values(by='Stability (%)', ascending=False)
stability_df.to_csv(os.path.join(file_path, "bootstrap_feature_stability_species_only.csv"))

# --- Youden's Index ---
youden_index = {}
y_test_bin = label_binarize(y_test, classes=np.unique(y))
for i, cls in enumerate(np.unique(y)):
    fpr, tpr, thresholds = roc_curve(y_test_bin[:, i], y_score[:, i])
    youden = tpr - fpr
    max_index = np.argmax(youden)
    youden_index[cls] = {
        'Youden': youden[max_index],
        'Best Threshold': thresholds[max_index],
        'Sensitivity': tpr[max_index],
        'Specificity': 1 - fpr[max_index]
    }
pd.DataFrame(youden_index).T.to_csv(os.path.join(file_path, "youden_index_species_only.csv"))

# --- Export files for ROC plotting (e.g., in Colab) ---
pd.DataFrame({'TrueLabel': y_test.values}).to_csv(os.path.join(file_path, "y_test_labels.csv"), index=False)
pd.DataFrame(y_score, columns=rf.classes_).to_csv(os.path.join(file_path, "y_score_probs.csv"), index=False)
pd.Series(rf.classes_, name="Classes").to_csv(os.path.join(file_path, "class_order.csv"), index=False)

# --- Plot Top 20 Permutation Importance Features ---
top20_perm = perm_df.head(20).reset_index()
plt.figure(figsize=(9, 6))
plt.barh(top20_perm['Feature'][::-1], top20_perm['Importance_Mean'][::-1], color="#2b8cbe")
plt.xlabel("Permutation Importance (Mean)", fontsize=12)
plt.title("Top 20 Species by Permutation Importance", fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig(os.path.join(file_path, "top20_species_permutation_barplot.png"), dpi=600)
plt.close()

# --- Merge Permutation + Stability ---
perm_df = perm_df.set_index('Feature')
merged_df = perm_df.join(stability_df, how='outer').fillna(0)
merged_df['Consistency'] = ((merged_df['Importance_Mean'] > 0) & 
                            (merged_df['Top20_Freq'] > 0)).astype(int)
merged_df['Avg_Rank'] = merged_df[['Importance_Mean', 'Stability (%)']].rank(ascending=False).mean(axis=1)
merged_df = merged_df.sort_values('Avg_Rank')

EOF
merged_df.to_csv(os.path.join(file_path, "consistent_species_feature_summary.csv"))

EOF
