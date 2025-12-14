

!pip install econml

# ============================================================
# Cross-Validated PLS (stable) + Cross-Fitted DML + Bootstrapping
# ============================================================

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from econml.dml import LinearDML
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from statsmodels.stats.multitest import multipletests
from scipy.stats import norm
import warnings

warnings.filterwarnings("ignore")

# ------------------------------------------------------------
# 🌱 Reproducibility
# ------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

# ------------------------------------------------------------
# 1) Load & align data
# ------------------------------------------------------------
tpm = pd.read_csv("Combined_CLR.csv", index_col=0)           # features × samples
meta = pd.read_csv("Metadata_Aligned_to_FilteredMatrix.csv") # has Sample_ID, Group, Age, BMI, Sex

common = tpm.columns.intersection(meta["Sample_ID"])
tpm = tpm[common]
meta = meta[meta["Sample_ID"].isin(common)].reset_index(drop=True)
tpm = tpm[meta["Sample_ID"]]  # enforce order

# Encode outcome (ordinal: Healthy=0, Adenoma=1, Cancer=2)
le = LabelEncoder()
meta["Y"] = le.fit_transform(meta["Group"])

# Confounders
meta["Sex"] = meta["Sex"].astype(str).str.lower().map({"male": 1, "m": 1, "female": 0, "f": 0})
Xconf = (
    meta[["Age", "BMI", "Sex"]]
    .apply(pd.to_numeric, errors="coerce")
    .fillna(meta[["Age", "BMI", "Sex"]].median())
    .values
)
Y = meta["Y"].values

# Feature matrix (samples × features)
X_raw = tpm.T.values
feature_names = tpm.index.tolist()
n_samples, n_features = X_raw.shape
print(f"✅ Data aligned — Samples: {n_samples} | Features: {n_features}")

# ------------------------------------------------------------
# 2) Stable PLS via cross-validated loadings (sign-aligned)
# ------------------------------------------------------------
n_components = 10
n_splits = 5
cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=SEED)

# Standardize features once (global scaler for stability)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_raw)

# Collect x_weights_ across folds, align signs, then average
weights_list = []
ref_weights = None

for fold_idx, (tr, te) in enumerate(cv.split(X_scaled, Y), 1):
    pls = PLSRegression(n_components=n_components)
    pls.fit(X_scaled[tr], Y[tr])
    W = pls.x_weights_.copy()  # shape: (features, n_components)

    if ref_weights is None:
        ref_weights = np.sign(W.sum(axis=0, keepdims=True))  # reference sign per component
        ref_weights[ref_weights == 0] = 1.0

    # Align signs component-wise with reference
    signs = np.sign((W * ref_weights).sum(axis=0, keepdims=True))
    signs[signs == 0] = 1.0
    W_aligned = W * signs
    weights_list.append(W_aligned)

# Average aligned weights
W_avg = np.mean(np.stack(weights_list, axis=2), axis=2)  # (features, n_components)

# Orthonormalize columns (optional but helpful)
# Normalize each component vector to unit norm
W_norm = W_avg / (np.linalg.norm(W_avg, axis=0, keepdims=True) + 1e-12)

# Stable PLS component scores for all samples:
# T = X_scaled @ W_norm
T_pls = X_scaled @ W_norm  # (samples × n_components)

# Save stable loadings and scores
loadings_df = pd.DataFrame(W_norm, index=feature_names, columns=[f"PLS_{i+1}" for i in range(n_components)])
loadings_df.to_csv("Stable_PLS_Loadings.csv")
scores_df = pd.DataFrame(T_pls, index=meta["Sample_ID"], columns=[f"PLS_{i+1}" for i in range(n_components)])
scores_df.to_csv("Stable_PLS_Scores.csv")
print(f"✅ Stable PLS built — {n_components} components (saved: Stable_PLS_Loadings.csv, Stable_PLS_Scores.csv)")

# ------------------------------------------------------------
# 3) Tuned base learners for DML
# ------------------------------------------------------------
gb_params = {
    "n_estimators": 500,
    "learning_rate": 0.1,
    "max_depth": 5,
    "subsample": 0.6,
    "random_state": SEED
}
rf_params = {
    "n_estimators": 300,
    "max_depth": 3,
    "min_samples_split": 2,
    "min_samples_leaf": 1,
    "max_features": "log2",
    "random_state": SEED,
    "n_jobs": -1
}

# ------------------------------------------------------------
# 4) Bootstrapped cross-fitted DML over PLS components
# ------------------------------------------------------------
B = 200                     # number of bootstrap replicates
components = [f"PLS_{i+1}" for i in range(n_components)]
checkpoint_every = 25
os.makedirs("dml_boot_partials", exist_ok=True)

boot_records = []
print(f"🚀 Bootstrapping DML with cross-fitting (B={B}, components={n_components})...")

for b in tqdm(range(1, B + 1)):
    # Bootstrap indices
    idx = np.random.choice(n_samples, size=n_samples, replace=True)
    Xb, Yb, Tb = Xconf[idx], Y[idx], T_pls[idx, :]

    for k in range(n_components):
        Ti = Tb[:, [k]]  # single-component as treatment

        # Cross-fitted DML
        dml = LinearDML(
            model_y=GradientBoostingRegressor(**gb_params),
            model_t=RandomForestRegressor(**rf_params),
            discrete_outcome=False,      # Y is treated as continuous ordinal
            discrete_treatment=False,
            cv=5,                        # cross-fitting
            random_state=SEED + b
        )
        dml.fit(Yb, Ti, X=Xb)
        eff = dml.effect(Xb)
        ate = float(np.mean(eff))

        boot_records.append({
            "Bootstrap": b,
            "Component": f"PLS_{k+1}",
            "ATE": ate
        })

    # checkpoint
    if b % checkpoint_every == 0:
        pd.DataFrame(boot_records).to_csv(f"dml_boot_partials/DML_boot_partial_up_to_{b}.csv", index=False)

boot_df = pd.DataFrame(boot_records)
boot_df.to_csv("DML_Bootstrap_Effects_ByComponent.csv", index=False)

# ------------------------------------------------------------
# 5) Summaries (mean, sd, CI, z, p, FDR)
# ------------------------------------------------------------
summ = (
    boot_df.groupby("Component")["ATE"]
    .agg(["mean", "std"])
    .rename(columns={"mean": "ATE_mean_boot", "std": "ATE_std_boot"})
    .reset_index()
)

# 95% CI (normal approximation from bootstrap SD)
summ["Boot_CI_Lower"] = summ["ATE_mean_boot"] - 1.96 * summ["ATE_std_boot"]
summ["Boot_CI_Upper"] = summ["ATE_mean_boot"] + 1.96 * summ["ATE_std_boot"]

# z, p
summ["z"] = summ["ATE_mean_boot"] / (summ["ATE_std_boot"] + 1e-12)
#summ["pval"] = 2 * (1 - norm.cdf(np.abs(summ["z"])))

# FDR
#summ["qval_FDR_BH"] = multipletests(summ["pval"].values, method="fdr_bh")[1]

# Stability score
summ["Stability_Score"] = np.abs(summ["ATE_mean_boot"]) / (summ["ATE_std_boot"] + 1e-12)

summ = summ.sort_values("ATE_mean_boot", ascending=True)  # negative (protective) on top
summ.to_csv("DML_PLS_Bootstrap_Summary.csv", index=False)

print("\n✅ Bootstrapping complete!")
print("📁 Saved:")
print("  - Stable_PLS_Loadings.csv")
print("  - Stable_PLS_Scores.csv")
print("  - DML_Bootstrap_Effects_ByComponent.csv")
print("  - DML_PLS_Bootstrap_Summary.csv")

print("\nTop summary (by ATE_mean_boot):")
print(summ.head(10))

# ============================================================
# 📊 Extract & Rank Top Genes per PLS Component + Plots
#  - Output table: Component, Rank, Gene, Abs_Loading, ATE_mean_boot, E-value (Point)
#  - Plots: Top genes for most stable PROGRESSION & PROTECTIVE components
# ============================================================

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.cross_decomposition import PLSRegression
from sklearn.preprocessing import StandardScaler

# -----------------------------
# Config
# -----------------------------
N_COMPONENTS = 10
TOP_N_GENES  = 10   # per component for plots/table
SEED = 42

# -----------------------------
# 1) Load data
# -----------------------------
tpm = pd.read_csv("Combined_CLR.csv", index_col=0)  # features x samples
meta = pd.read_csv("Metadata_Aligned_to_FilteredMatrix.csv")

# Align
common = tpm.columns.intersection(meta["Sample_ID"])
tpm = tpm[common]
meta = meta[meta["Sample_ID"].isin(common)].reset_index(drop=True)
tpm = tpm[meta["Sample_ID"]]

# -----------------------------
# 2) Load PLS loadings (or compute)
# -----------------------------
if os.path.exists("Stable_PLS_Loadings.csv"):
    loadings = pd.read_csv("Stable_PLS_Loadings.csv", index_col=0)
    # Ensure column names are PLS_1...PLS_N
    loadings.columns = [f"PLS_{i+1}" for i in range(loadings.shape[1])]
    print(f"✅ Loaded PLS loadings from Stable_PLS_Loadings.csv ({loadings.shape[0]} genes × {loadings.shape[1]} comps)")
else:
    print("ℹ️ Stable_PLS_Loadings.csv not found — refitting PLS for loadings...")
    # Encode labels (Healthy=0, Adenoma=1, Cancer=2)
    grp = meta["Group"].astype('category').cat.codes.values
    # Scale and fit PLS
    X_scaled = StandardScaler().fit_transform(tpm.T.values)  # samples x genes
    pls = PLSRegression(n_components=N_COMPONENTS)
    pls.fit(X_scaled, grp)
    loadings = pd.DataFrame(
        pls.x_loadings_,
        index=tpm.index,
        columns=[f"PLS_{i+1}" for i in range(N_COMPONENTS)]
    )
    loadings.to_csv("Stable_PLS_Loadings.csv")
    print(f"✅ Saved PLS loadings → Stable_PLS_Loadings.csv")

# -----------------------------
# 3) Load component-level stats (ATE + E-value)
# -----------------------------
if os.path.exists("PLS_DML_Evalue_Sensitivity.csv"):
    comp_stats = pd.read_csv("PLS_DML_Evalue_Sensitivity.csv")
    print("✅ Loaded component stats from PLS_DML_Evalue_Sensitivity.csv")
else:
    # Fallback: compute E-values from bootstrap summary if present
    if not os.path.exists("DML_PLS_Bootstrap_Summary.csv"):
        raise FileNotFoundError(
            "Missing PLS_DML_Evalue_Sensitivity.csv and DML_PLS_Bootstrap_Summary.csv. "
            "Please run the DML bootstrap + E-value step first."
        )
    summ = pd.read_csv("DML_PLS_Bootstrap_Summary.csv")
    def compute_evalue(estimate, lower_ci, upper_ci):
        rr = float(np.exp(abs(estimate)))
        rr_ci = float(np.exp(abs(lower_ci if estimate > 0 else upper_ci)))
        def ev(rr_):
            if rr_ <= 1:
                return 1.0
            return rr_ + np.sqrt(rr_ * (rr_ - 1))
        return ev(rr), ev(rr_ci)
    rows = []
    for _, r in summ.iterrows():
        e_point, e_ci = compute_evalue(r["ATE_mean_boot"], r["Boot_CI_Lower"], r["Boot_CI_Upper"])
        rows.append({
            "Component": r["Component"],
            "ATE_mean_boot": r["ATE_mean_boot"],
            "Boot_CI_Lower": r["Boot_CI_Lower"],
            "Boot_CI_Upper": r["Boot_CI_Upper"],
            "E-value (Point)": e_point,
            "E-value (CI)": e_ci,
            "Stability_Score": r.get("Stability_Score", np.nan)
        })
    comp_stats = pd.DataFrame(rows)
    comp_stats.to_csv("PLS_DML_Evalue_Sensitivity.csv", index=False)
    print("✅ Computed E-values → PLS_DML_Evalue_Sensitivity.csv")

# Keep only components that exist in loadings
valid_components = [c for c in comp_stats["Component"] if c in loadings.columns]
comp_stats = comp_stats[comp_stats["Component"].isin(valid_components)].copy()

# -----------------------------
# 4) Rank genes per component by |loading| and merge stats
# -----------------------------
abs_loadings = loadings.abs()

records = []
for comp in valid_components:
    ranked = abs_loadings[comp].sort_values(ascending=False).head(TOP_N_GENES)
    for rank, (gene, abs_w) in enumerate(ranked.items(), start=1):
        records.append({
            "Component": comp,
            "Rank": rank,
            "Gene": gene,
            "Abs_Loading": abs_w
        })

ranked_df = pd.DataFrame(records)

# Merge component-level ATE & E-value into gene-ranked table
ranked_annot = ranked_df.merge(
    comp_stats[["Component", "ATE_mean_boot", "E-value (Point)"]],
    on="Component",
    how="left"
).sort_values(["Component", "Rank"])

ranked_annot.to_csv("PLS_TopGenes_Ranked_with_ComponentStats.csv", index=False)
print(f"✅ Saved → PLS_TopGenes_Ranked_with_ComponentStats.csv ({len(ranked_annot)} rows)")

# -----------------------------
# 5) Pick most stable progression/protective components
#    (stability by Stability_Score, direction by ATE sign)
# -----------------------------
if "Stability_Score" not in comp_stats.columns:
    # compute a proxy if not provided
    if {"ATE_mean_boot", "Boot_CI_Lower", "Boot_CI_Upper"}.issubset(comp_stats.columns):
        comp_stats["ATE_std_boot_proxy"] = (comp_stats["Boot_CI_Upper"] - comp_stats["Boot_CI_Lower"]) / (2*1.96 + 1e-9)
        comp_stats["Stability_Score"] = comp_stats["ATE_mean_boot"].abs() / (comp_stats["ATE_std_boot_proxy"] + 1e-12)
    else:
        comp_stats["Stability_Score"] = comp_stats["ATE_mean_boot"].abs()

# Progression: ATE > 0, highest stability
prog_comp = comp_stats[comp_stats["ATE_mean_boot"] > 0].sort_values("Stability_Score", ascending=False)
# Protective: ATE < 0, highest stability
prot_comp = comp_stats[comp_stats["ATE_mean_boot"] < 0].sort_values("Stability_Score", ascending=False)

selected = []
if len(prog_comp) > 0:
    selected.append(("progression", prog_comp.iloc[0]["Component"]))
if len(prot_comp) > 0:
    selected.append(("protective", prot_comp.iloc[0]["Component"]))

print("🔎 Selected components:")
for kind, comp in selected:
    ate = float(comp_stats.loc[comp_stats["Component"] == comp, "ATE_mean_boot"].iloc[0])
    stab = float(comp_stats.loc[comp_stats["Component"] == comp, "Stability_Score"].iloc[0])
    ev   = float(comp_stats.loc[comp_stats["Component"] == comp, "E-value (Point)"].iloc[0])
    print(f"  • {comp} — {kind} | ATE={ate:.3f}, Stability={stab:.3f}, E-value={ev:.2f}")

# -----------------------------
# 6) Plot top genes for selected components
# -----------------------------
sns.set(style="white", context="talk")
plt.rcParams.update({
    "axes.linewidth": 2,
    "xtick.direction": "out",
    "ytick.direction": "out",
})

for kind, comp in selected:
    sub = ranked_annot[ranked_annot["Component"] == comp].sort_values("Rank").head(TOP_N_GENES).copy()
    direction = "Progression (↑ in disease)" if sub["ATE_mean_boot"].iloc[0] > 0 else "Protective (↓ in disease)"

    plt.figure(figsize=(6.5, 6.0))
    ax = sns.barplot(
        data=sub,
        x="Abs_Loading",
        y="Gene",
        color="#1f77b4" if kind == "protective" else "#d62728",  # blue for protective, red-ish for progression
        edgecolor="black"
    )
    # No numbers on bars, clean look
    ax.set_xlabel("Absolute Loading Strength", fontsize=8, weight="bold")
    ax.set_ylabel("")
    title = (
        f"{comp} — {direction}\n"
        f"ATE={sub['ATE_mean_boot'].iloc[0]:.3f}, "
        f"E-value={sub['E-value (Point)'].iloc[0]:.2f}"
    )
    ax.set_title(title, fontsize=10, weight="bold", pad=10)

    # Make borders thick & black
    for spine in ax.spines.values():
        spine.set_linewidth(2.2)
        spine.set_edgecolor("black")

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, weight="bold")

    sns.despine()
    plt.tight_layout()
    outname = f"TopGenes_{comp}_{'Progression' if kind=='progression' else 'Protective'}.png"
    plt.savefig(outname, dpi=600, bbox_inches="tight")
    plt.show()
    print(f"✅ Saved: {outname}")

print("🎉 Done.")

!pip install scikit-posthocs

# Functional Enrichment Script (Explicit PLS₁ & PLS₃)

import os
import re
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import kruskal

try:
    import scikit_posthocs as sp
    HAS_SCPH = True
except Exception:
    HAS_SCPH = False
    print("⚠️ scikit-posthocs not found. Dunn’s test will be skipped. Install with: pip install scikit-posthocs")

# -----------------------------
# File paths
# -----------------------------
RANKED_GENES_FILE = "PLS_TopGenes_Ranked_with_ComponentStats.csv"
MECH_MAP_FILE = "Combined_ARG_VF_Mechanisms.csv"
TPM_FILE = "Combined_CLR.csv"
META_FILE = "Metadata_Aligned_to_FilteredMatrix.csv"

os.makedirs("figs", exist_ok=True)

# -----------------------------
# Explicitly select components
# -----------------------------
ranked = pd.read_csv(RANKED_GENES_FILE)
selected_components = ["PLS_1", "PLS_3"]

comp_stats = ranked.groupby("Component")["ATE_mean_boot"].mean().reset_index()

print("✅ Selected components:")
for comp in selected_components:
    ate = comp_stats.loc[comp_stats["Component"] == comp, "ATE_mean_boot"].values[0]
    print(f"  • {comp} | ATE={ate:.3f}")

# -----------------------------
# Load mechanism mapping & clean gene names
# -----------------------------
mech_map = pd.read_csv(MECH_MAP_FILE)
required_mech_cols = {"Abbreviated_Gene_Name", "Mechanism", "Source"}
missing_mech_cols = required_mech_cols - set(mech_map.columns)
if missing_mech_cols:
    raise ValueError(f"Mechanism file missing columns: {missing_mech_cols}")

strip_prefix = lambda x: re.sub(r"^(ARG|VF)_", "", str(x)).strip()
mech_map["Gene_clean"] = mech_map["Abbreviated_Gene_Name"].apply(strip_prefix)

# -----------------------------
# Functional enrichment
# -----------------------------
def functional_enrichment(comp_name, out_csv):
    sub = ranked[ranked["Component"] == comp_name].drop_duplicates("Gene")
    sub["Gene_clean"] = sub["Gene"].apply(strip_prefix)
    annot = sub.merge(mech_map[["Gene_clean", "Mechanism", "Source"]], on="Gene_clean", how="left")

    summary = annot.groupby(["Component", "Source", "Mechanism"]).size().reset_index(name="Count")
    summary["Percent"] = summary["Count"] / summary["Count"].sum() * 100
    summary.to_csv(out_csv, index=False)
    print(f"✅ Functional enrichment saved → {out_csv}")
    return sub

genes_pls1 = functional_enrichment("PLS_1", "PLS1_Functional_Enrichment.csv")
genes_pls3 = functional_enrichment("PLS_3", "PLS3_Functional_Enrichment.csv")

# -----------------------------
# Per-sample mechanism abundance
# -----------------------------
tpm = pd.read_csv(TPM_FILE, index_col=0)
meta = pd.read_csv(META_FILE)
common_samples = tpm.columns.intersection(meta["Sample_ID"])
tpm = tpm[common_samples]
meta = meta[meta["Sample_ID"].isin(common_samples)].reset_index(drop=True)

selected_genes = set(genes_pls1["Gene_clean"]).union(genes_pls3["Gene_clean"])

tpm_long = tpm.T.reset_index().melt(id_vars="index", var_name="Gene", value_name="Abundance")
tpm_long.columns = ["Sample_ID", "Gene", "Abundance"]
tpm_long["Gene_clean"] = tpm_long["Gene"].apply(strip_prefix)
tpm_long = tpm_long[tpm_long["Gene_clean"].isin(selected_genes)]

merged = tpm_long.merge(mech_map[["Gene_clean", "Mechanism", "Source"]], on="Gene_clean", how="left")
merged = merged.merge(meta[["Sample_ID", "Group"]], on="Sample_ID", how="left").dropna(subset=["Mechanism"])

sample_mech = merged.groupby(["Sample_ID", "Group", "Mechanism"])["Abundance"].sum().reset_index()
sample_mech.to_csv("Mechanism_Abundance_Per_Sample_PLS13.csv", index=False)
print("✅ Saved mechanism abundances → Mechanism_Abundance_Per_Sample_PLS13.csv")

# -----------------------------
# Kruskal–Wallis tests
# -----------------------------
kw_res = []
for mech, dfm in sample_mech.groupby("Mechanism"):
    groups = [group["Abundance"].values for _, group in dfm.groupby("Group") if len(group) >= 3]
    if len(groups) >= 2:
        stat, p = kruskal(*groups)
        kw_res.append({"Mechanism": mech, "KW_stat": stat, "KW_pval": p})

kw_df = pd.DataFrame(kw_res)
kw_df["qval_FDR_BH"] = multipletests(kw_df["KW_pval"], method="fdr_bh")[1]
kw_df.to_csv("Mechanism_KW_PLS13.csv", index=False)
print("✅ KW results saved → Mechanism_KW_PLS13.csv")

# -----------------------------
# 6) Dunn’s posthoc (if applicable)
# -----------------------------
if HAS_SCPH:
    sig_mechs = kw_df[kw_df["qval_FDR_BH"] < 0.05]["Mechanism"]
    dunn_res = []
    for mech in sig_mechs:
        sub = sample_mech[sample_mech["Mechanism"] == mech]
        ph = sp.posthoc_dunn(sub, val_col="Abundance", group_col="Group", p_adjust="fdr_bh")
        ph["Mechanism"] = mech
        dunn_res.append(ph.reset_index().melt(id_vars=["index", "Mechanism"]))
    if dunn_res:
        dunn_df = pd.concat(dunn_res).rename(columns={"index": "Group1", "variable": "Group2", "value": "p_adj"})
        dunn_df.to_csv("Mechanism_Dunn_PLS13.csv", index=False)
        print("✅ Dunn's test saved → Mechanism_Dunn_PLS13.csv")

# ============================================================
# 🧬 Mechanism Tests per Component (PLS_1 and PLS_3, separately)
#     - Functional enrichment table per component
#     - Per-sample mechanism abundances
#     - Kruskal–Wallis per mechanism (FDR-BH)
#     - Dunn post-hoc on significant mechanisms (optional)
# ============================================================

import os
import re
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from scipy.stats import kruskal

# Optional Dunn's post-hoc
try:
    import scikit_posthocs as sp
    HAS_SCPH = True
except Exception:
    HAS_SCPH = False
    print("⚠️ scikit-posthocs not found. Dunn’s test will be skipped. Install with: pip install scikit-posthocs")

# -----------------------------
# File paths
# -----------------------------
RANKED_GENES_FILE = "PLS_TopGenes_Ranked_with_ComponentStats.csv"
MECH_MAP_FILE     = "Combined_ARG_VF_Mechanisms.csv"
TPM_FILE          = "Combined_CLR.csv"                        # features (genes) × samples (CLR)
META_FILE         = "Metadata_Aligned_to_FilteredMatrix.csv"  # includes Sample_ID, Group

os.makedirs("figs", exist_ok=True)

# -----------------------------
# Load inputs
# -----------------------------
ranked = pd.read_csv(RANKED_GENES_FILE)
mech_map = pd.read_csv(MECH_MAP_FILE)
tpm = pd.read_csv(TPM_FILE, index_col=0)
meta = pd.read_csv(META_FILE)

# Basic checks
req_mech_cols = {"Abbreviated_Gene_Name", "Mechanism", "Source"}
missing_mech = req_mech_cols - set(mech_map.columns)
if missing_mech:
    raise ValueError(f"Mechanism file missing columns: {missing_mech}")

if "Sample_ID" not in meta.columns or "Group" not in meta.columns:
    raise ValueError("Metadata must contain 'Sample_ID' and 'Group' columns.")

# Clean helper
strip_prefix = lambda x: re.sub(r"^(ARG|VF)_", "", str(x)).strip()

# Prepare mechanisms map with cleaned names
mech_map["Gene_clean"] = mech_map["Abbreviated_Gene_Name"].astype(str).apply(strip_prefix)

# Align samples
common_samples = tpm.columns.intersection(meta["Sample_ID"])
tpm = tpm[common_samples]
meta = meta[meta["Sample_ID"].isin(common_samples)].copy()
meta = meta.set_index("Sample_ID").loc[common_samples]  # enforce order

# Components to analyze separately
components = ["PLS_1", "PLS_3"]

def functional_enrichment_for_component(ranked_df, comp_name, mech_map_df, out_csv):
    """Make a mechanism count table for genes in one component."""
    sub = ranked_df[ranked_df["Component"] == comp_name].drop_duplicates("Gene").copy()
    if sub.empty:
        print(f"ℹ️ No genes found for {comp_name}. Skipping enrichment.")
        return pd.DataFrame(), pd.DataFrame()

    sub["Gene_clean"] = sub["Gene"].astype(str).apply(strip_prefix)
    annot = sub.merge(
        mech_map_df[["Gene_clean", "Mechanism", "Source"]],
        on="Gene_clean", how="left"
    )

    # Summary counts
    if annot["Mechanism"].notna().any():
        summary = (
            annot.groupby(["Component", "Source", "Mechanism"], dropna=True)
                 .size().reset_index(name="Count")
        )
        # Percent within component
        summary["Percent"] = summary["Count"] / summary["Count"].sum() * 100.0
        summary.to_csv(out_csv, index=False)
        print(f"✅ {comp_name} functional enrichment → {out_csv}")
    else:
        summary = pd.DataFrame()
        print(f"ℹ️ {comp_name}: no mechanisms mapped.")

    return summary, sub

def subset_tpm_by_clean_gene_names(tpm_df, selected_clean_names):
    """
    Select rows from TPM whose *cleaned* gene names are in selected_clean_names.
    The TPM rows keep their original names (with prefixes) but are filtered by cleaned match.
    """
    mask = [strip_prefix(g) in selected_clean_names for g in tpm_df.index]
    sub = tpm_df.loc[mask]
    return sub

for comp in components:
    print(f"\n🚀 Analyzing {comp}")

    # 1) Enrichment table (counts) for the component
    _, comp_genes = functional_enrichment_for_component(
        ranked, comp, mech_map, out_csv=f"{comp}_Functional_Enrichment.csv"
    )
    if comp_genes.empty:
        print(f"⚠️ {comp}: no genes to analyze; skipping to next component.")
        continue

    # 2) Select TPM rows for this component using CLEANED names
    selected_clean = set(comp_genes["Gene_clean"])
    tpm_sub = subset_tpm_by_clean_gene_names(tpm, selected_clean)

    if tpm_sub.empty:
        print(f"⚠️ {comp}: no matching TPM rows after cleaning. Check naming.")
        continue

    # 3) Melt TPM and merge with mechanisms + metadata
    tpm_long = (
        tpm_sub.T
              .reset_index()
              .melt(id_vars="index", var_name="Gene_raw", value_name="Abundance")
              .rename(columns={"index": "Sample_ID"})
    )
    tpm_long["Gene_clean"] = tpm_long["Gene_raw"].astype(str).apply(strip_prefix)

    merged = (
        tpm_long
        .merge(mech_map[["Gene_clean", "Mechanism", "Source"]], on="Gene_clean", how="left")
        .merge(meta[["Group"]], left_on="Sample_ID", right_index=True, how="left")
        .dropna(subset=["Mechanism", "Group"])
        .reset_index(drop=True)
    )

    if merged.empty:
        print(f"ℹ️ {comp}: merged table is empty after mapping/metadata merge.")
        continue

    # 4) Per-sample mechanism abundance
    sample_mech = (
        merged.groupby(["Sample_ID", "Group", "Mechanism"], as_index=False)["Abundance"].sum()
    )
    out_mech = f"{comp}_Mechanism_Abundance_Per_Sample.csv"
    sample_mech.to_csv(out_mech, index=False)
    print(f"✅ {comp} per-sample mechanism abundances → {out_mech}")

    # 5) Kruskal–Wallis per mechanism (across Group)
    kw_rows = []
    for mech, dfm in sample_mech.groupby("Mechanism"):
        # Ensure at least ~3 obs per group & at least 2 groups
        groups = []
        labels = []
        for g, sub in dfm.groupby("Group"):
            vals = sub["Abundance"].dropna().values
            if len(vals) >= 3:
                groups.append(vals)
                labels.append(g)
        if len(groups) >= 2:
            try:
                stat, p = kruskal(*groups)
                kw_rows.append({"Mechanism": mech, "KW_stat": stat, "KW_pval": p, "Groups": ", ".join(labels)})
            except Exception:
                pass

    kw_df = pd.DataFrame(kw_rows).sort_values("KW_pval") if kw_rows else pd.DataFrame(columns=["Mechanism","KW_stat","KW_pval","Groups"])

    if not kw_df.empty:
        kw_df["qval_FDR_BH"] = multipletests(kw_df["KW_pval"], method="fdr_bh")[1]
    else:
        kw_df["qval_FDR_BH"] = []
    kw_outfile = f"{comp}_Mechanism_KW.csv"
    kw_df.to_csv(kw_outfile, index=False)
    print(f"✅ {comp} Kruskal–Wallis results → {kw_outfile} ({len(kw_df)} mechanisms tested)")

    # 6) Dunn’s post-hoc for mechanisms significant after FDR (if package present)
    if HAS_SCPH and not kw_df.empty:
        sig_mechs = kw_df[kw_df["qval_FDR_BH"] < 0.05]["Mechanism"].unique()
        if len(sig_mechs) > 0:
            dunn_list = []
            for mech in sig_mechs:
                sub = sample_mech[sample_mech["Mechanism"] == mech].dropna(subset=["Abundance"]).copy()
                # Dunn test expects a long data frame with value column and group column
                ph = sp.posthoc_dunn(sub, val_col="Abundance", group_col="Group", p_adjust="fdr_bh")
                ph["Mechanism"] = mech
                dunn_list.append(
                    ph.reset_index().melt(id_vars=["index","Mechanism"], var_name="Group2", value_name="p_adj")
                      .rename(columns={"index":"Group1"})
                )
            dunn_df = pd.concat(dunn_list, ignore_index=True)
            dunn_outfile = f"{comp}_Mechanism_Dunn.csv"
            dunn_df.to_csv(dunn_outfile, index=False)
            print(f"✅ {comp} Dunn’s post-hoc (FDR) → {dunn_outfile}")
        else:
            print(f"ℹ️ {comp}: no mechanisms significant after FDR; Dunn’s skipped.")
    elif not HAS_SCPH:
        print("⚠️ scikit-posthocs not installed; Dunn’s post-hoc skipped.")

# ============================================================
# 🧬 Full PLS₁ + PLS₃ Gene Network (All Genes, All Hubs, No Mechanism Filter)
# ============================================================

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# ------------------------------------------------------------
# 1️⃣ Input files
# ------------------------------------------------------------
RANKED_FILE = "PLS_TopGenes_Ranked_with_ComponentStats.csv"

# ------------------------------------------------------------
# 2️⃣ Load ranked genes (no duplicate filtering)
# ------------------------------------------------------------
ranked = pd.read_csv(RANKED_FILE)
subset = ranked[ranked["Component"].isin(["PLS_1", "PLS_3"])].copy()

# Keep *all* genes belonging to PLS₁ and PLS₃
print(f"✅ Total gene entries from PLS₁ and PLS₃ (including duplicates): {subset.shape[0]}")

# ------------------------------------------------------------
# 3️⃣ Build network (no mechanism filter, simple intra-component linking)
# ------------------------------------------------------------
G = nx.Graph()

# Add all genes as nodes, tagged by component
for _, r in subset.iterrows():
    G.add_node(r["Gene"], component=r["Component"])

# Add simple edges within each component based on ranking proximity
for comp, group in subset.groupby("Component"):
    group = group.sort_values("Rank").reset_index(drop=True)
    genes = group["Gene"].tolist()
    for i in range(len(genes) - 1):
        G.add_edge(genes[i], genes[i + 1], weight=0.5)

print(f"✅ Network built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

# ------------------------------------------------------------
# 4️⃣ Compute network metrics (robust eigenvector centrality)
# ------------------------------------------------------------
deg_cent = nx.degree_centrality(G)

try:
    # robust version (handles disconnected nodes)
    eig_cent = nx.eigenvector_centrality_numpy(G)
except Exception as e:
    print("⚠️ Eigenvector centrality failed; using degree only:", e)
    eig_cent = {n: 0 for n in G.nodes()}

metrics = pd.DataFrame({
    "Gene": list(G.nodes()),
    "Component": [G.nodes[n]["component"] for n in G.nodes()],
    "Degree": [deg_cent[n] for n in G.nodes()],
    "Eigenvector": [eig_cent[n] for n in G.nodes()]
})
metrics.to_csv("PLS1_PLS3_Network_Metrics_AllGenes_NoMechanism.csv", index=False)
print("✅ Saved → PLS1_PLS3_Network_Metrics_AllGenes_NoMechanism.csv")

# Identify top 5 hubs per component
hub_pls1 = metrics[metrics["Component"] == "PLS_1"].nlargest(5, "Eigenvector")
hub_pls3 = metrics[metrics["Component"] == "PLS_3"].nlargest(5, "Eigenvector")

print("\n⭐ Top hubs — PLS₁:")
print(hub_pls1[["Gene", "Eigenvector", "Degree"]])
print("\n⭐ Top hubs — PLS₃:")
print(hub_pls3[["Gene", "Eigenvector", "Degree"]])

hub_genes = set(hub_pls1["Gene"].tolist() + hub_pls3["Gene"].tolist())

# ------------------------------------------------------------
# 5️⃣ Plot network (color by component, gold halo hubs)
# ------------------------------------------------------------
sns.set(style="white")
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "axes.linewidth": 2.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

fig, ax = plt.subplots(figsize=(7, 8))
pos = nx.spring_layout(G, k=1.2, seed=42)

color_map = {"PLS_1": "#E74C3C", "PLS_3": "#2ECC71"}  # red / green
node_colors = [color_map.get(G.nodes[n]["component"], "gray") for n in G.nodes()]
node_sizes = [1000 * deg_cent[n] + 350 for n in G.nodes()]

# Draw edges and nodes
nx.draw_networkx_edges(G, pos, alpha=0.3, edge_color="gray", width=0.8, ax=ax)
nx.draw_networkx_nodes(
    G, pos,
    node_color=node_colors, node_size=node_sizes,
    alpha=0.8, edgecolors="black", linewidths=0.6, ax=ax
)

# Highlight hub genes (gold halo)
hubs_present = [h for h in hub_genes if h in G.nodes()]
nx.draw_networkx_nodes(
    G, pos, nodelist=hubs_present,
    node_color=[color_map[G.nodes[h]['component']] for h in hubs_present],
    node_size=[node_sizes[list(G.nodes()).index(h)] * 1.6 for h in hubs_present],
    edgecolors="gold", linewidths=3.0, ax=ax
)

# Labels (small, bold)
nx.draw_networkx_labels(G, pos, font_size=6.5, font_weight="bold", font_color="black", ax=ax)

plt.title("Full Gene Network — PLS₁ (Red) and PLS₃ (Green)\nTop 5 Hubs per Component Highlighted in Gold",
          fontsize=13, fontweight="bold", pad=10)
plt.axis("off")
plt.tight_layout()
plt.savefig("PLS1_PLS3_AllGenes_NoMechanism_Network.png", dpi=600, bbox_inches="tight", transparent=False)
plt.show()

print("✅ Saved → PLS1_PLS3_AllGenes_NoMechanism_Network.png (all genes + hubs highlighted)")
