#!/bin/bash
#SBATCH --job-name=rf_external_validation
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --time=10:00:00
#SBATCH --partition=compute
#SBATCH --output=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/out/rf_validation_%j.out
#SBATCH --error=/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal/err/rf_validation_%j.err


module load Python/3.11.5-GCCcore-13.2.0
source ~/python_env/bin/activate

python3 <<EOF
import pandas as pd
import numpy as np
import seaborn as sns
import os
import joblib
import matplotlib.pyplot as plt
from sklearn.metrics import (
    classification_report, confusion_matrix, accuracy_score, f1_score,
    cohen_kappa_score, matthews_corrcoef, roc_curve, auc
)
from sklearn.preprocessing import LabelEncoder, label_binarize

# === CONFIGURATION ===
base_path = "/srv/lustre01/project/mmrd-cp3fk69sfrq/morad.mokhtar/Colorectal"
model_file = os.path.join(base_path, "best_rf_model_microbial.pkl")
synthetic_files = {
    "balanced": os.path.join(base_path, "PCA_synthetic/PCA_synthetic_hellinger_balanced_PCA.csv"),
    "real_ratio": os.path.join(base_path, "PCA_synthetic/PCA_synthetic_hellinger_real_ratio.csv")
}
metadata_files = {
    "balanced": os.path.join(base_path, "PCA_synthetic/PCA_synthetic_metadata_balanced_PCA.csv"),
    "real_ratio": os.path.join(base_path, "PCA_synthetic/PCA_synthetic_metadata_real_ratio.csv")
}
output_metrics = os.path.join(base_path, "rf_external_validation_results.csv")
output_figures = os.path.join(base_path, "rf_external_validation_roc")
os.makedirs(output_figures, exist_ok=True)

# === Load trained model ===
print("Loading best RF model...")
model = joblib.load(model_file)
results = []

for tag in synthetic_files:
    print(f" Validating on {tag} synthetic dataset...")
    syn_X = pd.read_csv(synthetic_files[tag])
    syn_meta = pd.read_csv(metadata_files[tag])
    df = pd.merge(syn_X, syn_meta, on="Sample_ID")
    y_true = df["Group"]
    drop_cols = ["Sample_ID", "Group", "Instrument", "Project", "Center_Name", "Continent", "Country", "Age", "BMI", "Sex"]
    X_test = df.drop(columns=drop_cols, errors="ignore").select_dtypes(include=[np.number])

    le = LabelEncoder()
    y_encoded = le.fit_transform(y_true)
    y_binarized = label_binarize(y_encoded, classes=range(len(le.classes_)))
    y_pred = model.predict(X_test)
    y_score = model.predict_proba(X_test)

    # Save classification report
    report_dict = classification_report(y_encoded, y_pred, target_names=le.classes_, output_dict=True)
    pd.DataFrame(report_dict).transpose().to_csv(
        os.path.join(base_path, f"classification_report_{tag}.csv"), index=True
    )

    # Save confusion matrix CSV
    conf_matrix = confusion_matrix(y_encoded, y_pred)
    conf_df = pd.DataFrame(conf_matrix, index=le.classes_, columns=le.classes_)
    conf_csv_path = os.path.join(base_path, f"confusion_matrix_{tag}.csv")
    conf_df.to_csv(conf_csv_path)
    print(f" Confusion matrix saved to: {conf_csv_path}")

    # Plot confusion matrix
    plt.figure(figsize=(6, 5))
    sns.heatmap(
        conf_df,
        annot=True,
        fmt="d",
        cmap="Blues",
        cbar=False,
        linewidths=1.5,
        linecolor="black",
        square=True
    )
    plt.title(f"Confusion Matrix - {tag}", fontsize=14)
    plt.xlabel("Predicted Label", fontsize=12)
    plt.ylabel("True Label", fontsize=12)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    plt.tight_layout()
    conf_fig_path = os.path.join(output_figures, f"confusion_matrix_{tag}.png")
    plt.savefig(conf_fig_path, dpi=600)
    plt.close()
    print(f" Confusion matrix plot saved to: {conf_fig_path}")

    # Print and record results
    print(f" Results on {tag} dataset:")
    print(classification_report(y_encoded, y_pred, target_names=le.classes_))
    print("Confusion Matrix:")
    print(conf_matrix)

    results.append({
        "Dataset": tag,
        "Accuracy": accuracy_score(y_encoded, y_pred),
        "F1 Macro": f1_score(y_encoded, y_pred, average="macro"),
        "Kappa": cohen_kappa_score(y_encoded, y_pred),
        "MCC": matthews_corrcoef(y_encoded, y_pred)
    })

    # === Plot ROC curve ===
    fpr, tpr, roc_auc = {}, {}, {}
    for i in range(len(le.classes_)):
        fpr[i], tpr[i], _ = roc_curve(y_binarized[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    plt.figure(figsize=(6, 5))
    for i, label in enumerate(le.classes_):
        plt.plot(fpr[i], tpr[i], lw=2, label=f"{label} (AUC = {roc_auc[i]:.2f})")
    plt.plot([0, 1], [0, 1], linestyle="--", color="gray")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title(f"Multiclass ROC - {tag}")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(output_figures, f"roc_curve_{tag}.png"), dpi=600)
    plt.close()

# Save overall metrics
results_df = pd.DataFrame(results)
results_df.to_csv(output_metrics, index=False)
print("\n External validation completed. Results saved.")
EOF
