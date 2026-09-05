from pathlib import Path
import gzip

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data" / "GSE148158_normalized_counts_MS2.csv.gz"
OUT = ROOT / "results" / "GSE148158"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLE_GROUPS = {
    "BJ_2 [re-analysis]": "BJ_fibroblast", "BJ_2": "BJ_fibroblast",
    "BJ_1 [re-analysis]": "BJ_fibroblast", "BJ_1": "BJ_fibroblast",
    "BJ_3 [re-analysis]": "BJ_fibroblast", "BJ_3": "BJ_fibroblast",
    "BJ_4 [re-analysis]": "BJ_fibroblast", "BJ_4": "BJ_fibroblast",
    "H1_2 [re-analysis]": "hESC", "H1_2": "hESC",
    "H9 [re-analysis]": "hESC", "H9": "hESC",
    "H1 [re-analysis]": "hESC", "H1": "hESC",
    "BJ_GFP48": "GFP_48h", "BJ_GFP48b": "GFP_48h",
    "BJ_GFP72": "GFP_72h", "BJ_GFP72b": "GFP_72h",
    "OSKM48": "OSKM_48h", "OSKM72": "OSKM_72h",
}


def group_for_sample(sample):
    s = str(sample).strip()
    if s in SAMPLE_GROUPS:
        return SAMPLE_GROUPS[s]
    for key, group in SAMPLE_GROUPS.items():
        if s.startswith(key) or key in s:
            return group
    return "unclassified"


def read_expression(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        df = pd.read_csv(f, low_memory=False)
    id_col = df.columns[0]
    x = df.drop(columns=id_col).apply(pd.to_numeric, errors="coerce")
    x = x.loc[:, x.notna().sum() > 0]
    x.index = df[id_col].astype(str).str.strip()
    x = x.replace([np.inf, -np.inf], np.nan).dropna(how="all")
    return x.groupby(level=0).mean()


def save_plot(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=250)
    plt.close(fig)


expr = read_expression(DATA)
expr.to_csv(OUT / "expression.csv")

metadata = pd.DataFrame({"sample": expr.columns})
metadata["group"] = metadata["sample"].map(group_for_sample)
metadata.to_csv(OUT / "01_sample_metadata.csv", index=False)

fig, ax = plt.subplots(figsize=(12, 6))
ax.boxplot([expr[c].dropna() for c in expr.columns], tick_labels=expr.columns, showfliers=False)
ax.set_title("GSE148158 - normalized expression distributions")
ax.set_ylabel("Normalized expression")
ax.tick_params(axis="x", rotation=45)
save_plot(fig, "02_boxplot.png")

sample_summary = pd.DataFrame({"mean": expr.mean(), "median": expr.median(), "sd": expr.std(),
                               "min": expr.min(), "max": expr.max(), "missing": expr.isna().sum()})
sample_summary["group"] = sample_summary.index.map(group_for_sample)
sample_summary.to_csv(OUT / "03_sample_QC.csv")

corr = expr.corr(method="pearson")
corr.to_csv(OUT / "04_sample_correlation.csv")
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(corr, vmin=-1, vmax=1, aspect="auto")
ax.set_xticks(range(len(corr.columns))); ax.set_xticklabels(corr.columns, rotation=45, ha="right")
ax.set_yticks(range(len(corr.index))); ax.set_yticklabels(corr.index)
ax.set_title("GSE148158 - sample correlation")
fig.colorbar(im, ax=ax, label="Pearson r")
save_plot(fig, "05_sample_correlation.png")

log_expr = np.log2(expr.clip(lower=0) + 1)
variance = log_expr.var(axis=1).sort_values(ascending=False)
variance.to_csv(OUT / "06_feature_variance.csv", header=["variance"])
top = variance.head(min(2000, len(variance))).index
X = log_expr.loc[top].T.fillna(log_expr.loc[top].T.mean())
centered = X.to_numpy() - X.to_numpy().mean(axis=0)
U, S, _ = np.linalg.svd(centered, full_matrices=False)
PC = U * S
EV = S**2 / np.sum(S**2)
coords = pd.DataFrame(PC[:, :min(5, PC.shape[1])], index=X.index,
                      columns=[f"PC{i+1}" for i in range(min(5, PC.shape[1]))])
coords["group"] = coords.index.map(group_for_sample)
coords.to_csv(OUT / "07_PCA_coordinates.csv")
if len(EV) >= 2:
    fig, ax = plt.subplots(figsize=(9, 7))
    for group, idx in coords.groupby("group").groups.items():
        positions = [X.index.get_loc(s) for s in idx]
        ax.scatter(PC[positions, 0], PC[positions, 1], s=90, label=group)
        for p in positions:
            ax.annotate(X.index[p], (PC[p, 0], PC[p, 1]), xytext=(6, 6), textcoords="offset points", fontsize=8)
    ax.set_xlabel(f"PC1 ({EV[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({EV[1]*100:.1f}%)")
    ax.set_title("GSE148158 - PCA by biological group")
    ax.legend(fontsize=8)
    save_plot(fig, "08_PCA_by_group.png")

features = variance.head(min(50, len(variance))).index
z = log_expr.loc[features]
z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
fig, ax = plt.subplots(figsize=(10, max(6, len(features)*0.18)))
im = ax.imshow(z, aspect="auto", interpolation="nearest")
ax.set_xticks(range(len(z.columns))); ax.set_xticklabels(z.columns, rotation=45, ha="right")
ax.set_yticks(range(len(z.index))); ax.set_yticklabels(z.index, fontsize=6)
ax.set_title("GSE148158 - top variable features")
fig.colorbar(im, ax=ax, label="row z-score")
save_plot(fig, "09_top_variable_features.png")

group_counts = metadata["group"].value_counts().sort_index()
report = [
    "Dataset: GSE148158",
    "Experiment type: expression profiling by high-throughput sequencing (RNA-seq)",
    f"Features: {expr.shape[0]:,}",
    f"Samples: {expr.shape[1]:,}",
    "Input: GEO normalized counts",
    "Transformation for EDA/PCA: log2(x + 1)",
    "No raw-count RNA-seq normalization was performed.",
    "Sample groups:",
]
report += [f"  {group}: {count}" for group, count in group_counts.items()]
if len(EV) >= 2:
    report += [f"PC1: {EV[0]*100:.2f}%", f"PC2: {EV[1]*100:.2f}%", f"PC1+PC2: {(EV[0]+EV[1])*100:.2f}%"]
(OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")
print("GSE148158 exploration complete")
