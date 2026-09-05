from pathlib import Path
import gzip
import tarfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
ARCHIVE = ROOT / "Data" / "GSE52052_RAW.tar"
OUT = ROOT / "results" / "GSE52052"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLE_GROUPS = {
    "HDF_GFP(+)_day11": "GFP_control",
    "HDF_OSK+control_inh_TRA(+)_day11": "OSK_control_inhibitor",
    "HDF_OSK+let-7_inh_TRA(+)_day11": "OSK_let7_inhibitor",
    "HDF_OSKM_TRA(+)_day11": "OSKM",
    "HDF_OSK+LIN-41_TRA(+)_day11": "OSK_LIN41",
    "H1_hESC": "H1_hESC",
}

extract = ROOT / "GSE52052_extracted"
extract.mkdir(exist_ok=True)
if not any(extract.rglob("*.txt.gz")):
    with tarfile.open(ARCHIVE, "r") as tar:
        tar.extractall(extract)

signals = []
for path in sorted(extract.rglob("*.txt.gz")):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    header_i = next(i for i, line in enumerate(lines) if line.startswith("FEATURES\t"))
    cols = lines[header_i].split("\t")[1:]
    rows = []
    for line in lines[header_i + 1:]:
        if line.startswith("*"):
            break
        if line.startswith("DATA\t"):
            fields = line.split("\t")[1:]
            if len(fields) == len(cols):
                rows.append(fields)
    df = pd.DataFrame(rows, columns=cols)
    ct = pd.to_numeric(df["ControlType"], errors="coerce")
    df = df[ct.fillna(0) == 0].copy()
    df["gProcessedSignal"] = pd.to_numeric(df["gProcessedSignal"], errors="coerce")
    df = df.dropna(subset=["gProcessedSignal"])
    df["ProbeName"] = df["ProbeName"].astype(str).str.strip()
    df = df[(df["ProbeName"] != "") & (df["ProbeName"] != "nan")]
    df.loc[df["gProcessedSignal"] < 0, "gProcessedSignal"] = 0
    sample = path.name[:-7]
    signals.append(df.groupby("ProbeName")["gProcessedSignal"].mean().rename(sample))

expr = pd.concat(signals, axis=1, join="inner")
expr.to_csv(OUT / "expression_raw_processed.csv")

metadata = pd.DataFrame({"sample": expr.columns})
metadata["group"] = metadata["sample"].map(SAMPLE_GROUPS).fillna("unclassified")
metadata.to_csv(OUT / "01_sample_metadata.csv", index=False)

log = np.log2(expr + 1)
log.to_csv(OUT / "expression_log2.csv")

fig, ax = plt.subplots(figsize=(12, 6))
ax.boxplot([log[c].dropna() for c in log], tick_labels=log.columns, showfliers=False)
ax.set_title("GSE52052 - log2 processed signal")
ax.set_ylabel("log2(signal + 1)")
ax.tick_params(axis="x", rotation=45)
fig.tight_layout(); fig.savefig(OUT / "02_boxplot.png", dpi=250); plt.close(fig)

qc = pd.DataFrame({"mean": expr.mean(), "median": expr.median(), "sd": expr.std(), "missing": expr.isna().sum()})
qc["group"] = qc.index.map(SAMPLE_GROUPS).fillna("unclassified")
qc.to_csv(OUT / "03_sample_QC.csv")

corr = log.corr()
corr.to_csv(OUT / "04_sample_correlation.csv")
fig, ax = plt.subplots(figsize=(8, 7))
im = ax.imshow(corr, vmin=-1, vmax=1)
ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns, rotation=45, ha="right")
ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.index)
ax.set_title("GSE52052 - sample correlation")
fig.colorbar(im, ax=ax, label="Pearson r")
fig.tight_layout(); fig.savefig(OUT / "05_sample_correlation.png", dpi=250); plt.close(fig)

# GEO reports GeneSpring percentile-normalized values in the sample tables.
# The exploratory quantile normalization below is therefore an additional EDA transformation,
# not a replacement for the GEO processing pipeline.
x = log.to_numpy(float)
order = np.argsort(x, axis=0)
means = np.sort(x, axis=0).mean(axis=1)
norm = np.empty_like(x)
for j in range(x.shape[1]):
    norm[order[:, j], j] = means
norm = pd.DataFrame(norm, index=log.index, columns=log.columns)
norm.to_csv(OUT / "06_exploratory_quantile_normalized.csv")

var = norm.var(axis=1).sort_values(ascending=False)
var.to_csv(OUT / "07_probe_variance.csv", header=["variance"])
top = var.head(min(2000, len(var))).index
X = norm.loc[top].T
C = X.to_numpy() - X.to_numpy().mean(axis=0)
U, S, _ = np.linalg.svd(C, full_matrices=False)
PC = U * S
EV = S**2 / np.sum(S**2)
n_pc = min(5, PC.shape[1])
coords = pd.DataFrame(PC[:, :n_pc], index=X.index, columns=[f"PC{i+1}" for i in range(n_pc)])
coords["group"] = coords.index.map(SAMPLE_GROUPS).fillna("unclassified")
coords.to_csv(OUT / "08_PCA_coordinates.csv")

if len(EV) >= 2:
    fig, ax = plt.subplots(figsize=(9, 7))
    for group, idx in coords.groupby("group").groups.items():
        positions = [X.index.get_loc(s) for s in idx]
        ax.scatter(PC[positions, 0], PC[positions, 1], s=90, label=group)
        for p in positions:
            ax.annotate(X.index[p], (PC[p, 0], PC[p, 1]), xytext=(6, 6), textcoords="offset points", fontsize=8)
    ax.set_xlabel(f"PC1 ({EV[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({EV[1]*100:.1f}%)")
    ax.set_title("GSE52052 - PCA by experimental group")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "09_PCA_by_group.png", dpi=250); plt.close(fig)

features = var.head(min(50, len(var))).index
z = norm.loc[features]
z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
fig, ax = plt.subplots(figsize=(10, max(6, len(features)*0.18)))
im = ax.imshow(z, aspect="auto")
ax.set_xticks(range(len(z.columns))); ax.set_xticklabels(z.columns, rotation=45, ha="right")
ax.set_yticks(range(len(z.index))); ax.set_yticklabels(z.index, fontsize=6)
ax.set_title("GSE52052 - top variable probes")
fig.colorbar(im, ax=ax, label="row z-score")
fig.tight_layout(); fig.savefig(OUT / "10_top_variable_probes.png", dpi=250); plt.close(fig)

group_counts = metadata["group"].value_counts().sort_index()
report = [
    "Dataset: GSE52052",
    "Experiment type: expression profiling by array",
    "Platform: Agilent-028004 SurePrint G3 Human GE 8x60K (GPL14550)",
    f"Probes: {expr.shape[0]:,}",
    f"Samples: {expr.shape[1]:,}",
    "Input: GSE52052_RAW.tar",
    "Signal: Agilent gProcessedSignal",
    "Control probes: removed",
    "Duplicate ProbeName: averaged",
    "GEO processing: GeneSpring GX percentile normalization",
    "Exploratory transformation: log2(x + 1)",
    "Additional exploratory quantile normalization: applied after log2; not a replacement for GEO processing.",
    "Experimental groups:",
]
report += [f"  {group}: {count}" for group, count in group_counts.items()]
if len(EV) >= 2:
    report += [f"PC1: {EV[0]*100:.2f}%", f"PC2: {EV[1]*100:.2f}%", f"PC1+PC2: {(EV[0]+EV[1])*100:.2f}%"]
(OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")
print("GSE52052 exploration complete")
