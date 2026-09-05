from pathlib import Path
import gzip
import tarfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUT = ROOT / "results" / "GSE28688"
OUT.mkdir(parents=True, exist_ok=True)

SAMPLE_GROUPS = {
    "HFF1-a": "HFF1",
    "HFF1-b": "HFF1",
    "HFF1-24 h post-transduction-a": "24h",
    "HFF1-24 h post-transduction-b": "24h",
    "HFF1-48 h post-transduction-a": "48h",
    "HFF1-48 h post-transduction-b": "48h",
    "HFF1-72 h post-transduction-a": "72h",
    "HFF1-72 h post-transduction-b": "72h",
    "H1": "H1_hESC",
    "H9": "H9_hESC",
    "iPS2 from HFF1-a": "iPS2",
    "iPS2 from HFF1-b": "iPS2",
    "iPS4 from HFF1-a": "iPS4",
    "iPS4 from HFF1-b": "iPS4",
}


def read_table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        sample = f.read(8192)
    sep = "," if sample.count(",") > sample.count("\t") else "\t"
    return pd.read_csv(path, sep=sep, compression="infer", low_memory=False)


def expression(df):
    df = df.dropna(how="all").dropna(axis=1, how="all")
    names = {"gene", "geneid", "gene_id", "symbol", "gene_symbol", "probe", "id", "feature", "feature_id", "id_ref", "probe_id"}
    id_col = next((c for c in df.columns if str(c).strip().lower() in names), df.columns[0])
    x = df.drop(columns=id_col).apply(pd.to_numeric, errors="coerce")
    x = x.loc[:, x.notna().sum() > max(2, int(0.5 * len(df)))]
    if x.shape[1] < 2:
        raise ValueError("Brak co najmniej dwóch kolumn ekspresji.")
    x.index = df[id_col].astype(str).str.strip()
    x = x[(x.index != "") & (x.index != "nan")]
    return x.groupby(level=0).mean().replace([np.inf, -np.inf], np.nan).dropna(how="all")


def extract_raw_archive(path, extract):
    extract.mkdir(exist_ok=True)
    marker = extract / ".archive_extracted"
    if not marker.exists():
        with tarfile.open(path, "r") as tar:
            tar.extractall(extract)
        marker.write_text("ok", encoding="utf-8")
    return [f for f in extract.rglob("*") if f.is_file() and f.name != marker.name]


def inspect_bgx(path, out_dir):
    rows = []
    header = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                lowered = [p.strip().lower() for p in parts]
                if any(x in lowered for x in ("probeid", "probe_id", "array_address_id")):
                    header = parts
                    continue
                if len(parts) > 3:
                    header = parts
                    continue
            if header is not None and len(parts) == len(header):
                rows.append(parts)
            if len(rows) >= 100000:
                break
    if header is None or not rows:
        return {"path": path.name, "parsed": False, "rows": 0, "columns": []}
    df = pd.DataFrame(rows, columns=header)
    df.to_csv(out_dir / "BGX_annotation_sample.csv", index=False)
    return {"path": path.name, "parsed": True, "rows": len(df), "columns": list(df.columns)}


def raw_archive_info(archive_path, extract):
    files = extract_raw_archive(archive_path, extract)
    info_dir = OUT / "RAW_archive"
    info_dir.mkdir(exist_ok=True)
    print("\n=== GSE28688 RAW archive contents ===")
    print(f"Archive members: {len(files)}")
    bgx_result = None
    for i, path in enumerate(files, 1):
        size = path.stat().st_size
        print(f"[{i:02d}] {path.name} | {size:,} bytes")
        if path.name.lower().endswith(".bgx.gz"):
            bgx_result = inspect_bgx(path, info_dir)
            print(f"     BGX parsed: {bgx_result['parsed']}")
            if bgx_result["parsed"]:
                print(f"     annotation rows read: {bgx_result['rows']:,}")
                print(f"     columns: {', '.join(map(str, bgx_result['columns']))}")
    print("=== end RAW archive contents ===")
    report = [
        "Dataset: GSE28688",
        "RAW archive is inspected as supplementary platform/annotation content.",
        f"Archive members: {len(files)}",
    ]
    if bgx_result:
        report += [
            f"BGX file: {bgx_result['path']}",
            f"BGX parsed: {bgx_result['parsed']}",
            f"BGX rows read: {bgx_result['rows']}",
            "BGX columns: " + ", ".join(map(str, bgx_result["columns"])),
        ]
    (info_dir / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")


def qnorm(df):
    a = df.to_numpy(float)
    order = np.argsort(a, axis=0)
    means = np.sort(a, axis=0).mean(axis=1)
    out = np.empty_like(a)
    for j in range(a.shape[1]):
        out[order[:, j], j] = means
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def explore(expr, label):
    out = OUT / label
    out.mkdir(exist_ok=True)
    expr = expr.clip(lower=0)
    expr.to_csv(out / "expression_input.csv")
    log = np.log2(expr + 1)
    log.to_csv(out / "expression_log2.csv")
    norm = qnorm(log)
    norm.to_csv(out / "expression_normalized.csv")

    metadata = pd.DataFrame({"sample": expr.columns})
    metadata["group"] = metadata["sample"].map(SAMPLE_GROUPS).fillna("unclassified")
    metadata.to_csv(out / "01_sample_metadata.csv", index=False)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.boxplot([log[c].dropna() for c in log], tick_labels=log.columns, showfliers=False)
    ax.set_title(f"GSE28688 - {label} - log2 distributions")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(out / "02_boxplot.png", dpi=250); plt.close(fig)

    qc = pd.DataFrame({"mean": expr.mean(), "median": expr.median(), "sd": expr.std(), "missing": expr.isna().sum()})
    qc["group"] = qc.index.map(SAMPLE_GROUPS).fillna("unclassified")
    qc.to_csv(out / "03_sample_QC.csv")

    corr = norm.corr()
    corr.to_csv(out / "04_sample_correlation.csv")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(corr, vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.index)
    ax.set_title(f"GSE28688 - {label} - sample correlation")
    fig.colorbar(im, ax=ax, label="Pearson r")
    fig.tight_layout(); fig.savefig(out / "05_sample_correlation.png", dpi=250); plt.close(fig)

    var = norm.var(axis=1).sort_values(ascending=False)
    var.to_csv(out / "06_feature_variance.csv", header=["variance"])
    top = var.head(min(2000, len(var))).index
    X = norm.loc[top].T.fillna(norm.loc[top].T.mean())
    C = X.to_numpy() - X.to_numpy().mean(axis=0)
    U, S, _ = np.linalg.svd(C, full_matrices=False)
    PC = U * S
    EV = S**2 / np.sum(S**2)
    n_pc = min(5, PC.shape[1])
    coords = pd.DataFrame(PC[:, :n_pc], index=X.index, columns=[f"PC{i + 1}" for i in range(n_pc)])
    coords["group"] = coords.index.map(SAMPLE_GROUPS).fillna("unclassified")
    coords.to_csv(out / "07_PCA_coordinates.csv")

    if len(EV) >= 2:
        fig, ax = plt.subplots(figsize=(9, 7))
        for group, idx in coords.groupby("group").groups.items():
            positions = [X.index.get_loc(s) for s in idx]
            ax.scatter(PC[positions, 0], PC[positions, 1], s=90, label=group)
            for p in positions:
                ax.annotate(X.index[p], (PC[p, 0], PC[p, 1]), xytext=(6, 6), textcoords="offset points", fontsize=8)
        ax.set_xlabel(f"PC1 ({EV[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({EV[1]*100:.1f}%)")
        ax.set_title(f"GSE28688 - {label} - PCA by biological group")
        ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(out / "08_PCA_by_group.png", dpi=250); plt.close(fig)

    features = var.head(min(50, len(var))).index
    z = norm.loc[features]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(10, max(6, len(features) * 0.18)))
    im = ax.imshow(z, aspect="auto")
    ax.set_xticks(range(len(z.columns))); ax.set_xticklabels(z.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(z.index))); ax.set_yticklabels(z.index, fontsize=6)
    ax.set_title(f"GSE28688 - {label} - top variable features")
    fig.colorbar(im, ax=ax, label="row z-score")
    fig.tight_layout(); fig.savefig(out / "09_top_variable_features.png", dpi=250); plt.close(fig)

    group_counts = metadata["group"].value_counts().sort_index()
    report = [
        f"Dataset: GSE28688 ({label})",
        "Experiment type: expression profiling by array",
        "Platform: Illumina HumanRef-8 v3.0 (GPL6883)",
        f"Features: {expr.shape[0]:,}",
        f"Samples: {expr.shape[1]:,}",
        "Input: GSE28688_non-normalized.txt.gz",
        "Transformation: log2(x + 1)",
        "Quantile normalization: exploratory normalization applied after log2 transformation.",
        "Sample groups:",
    ]
    report += [f"  {group}: {count}" for group, count in group_counts.items()]
    if len(EV) >= 2:
        report += [f"PC1: {EV[0]*100:.2f}%", f"PC2: {EV[1]*100:.2f}%", f"PC1+PC2: {(EV[0]+EV[1])*100:.2f}%"]
    (out / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")


p = DATA / "GSE28688_non-normalized.txt.gz"
if p.exists():
    explore(expression(read_table(p)), "non_normalized")
else:
    print("Brak GSE28688_non-normalized.txt.gz")

p = DATA / "GSE28688_RAW.tar"
if p.exists():
    raw_archive_info(p, ROOT / "GSE28688_extracted")
else:
    print("Brak GSE28688_RAW.tar")

print("GSE28688 exploration complete")
