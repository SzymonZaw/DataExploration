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


def raw_single_sample(path):
    df = read_table(path).dropna(how="all").dropna(axis=1, how="all")
    lower = {str(c).strip().lower(): c for c in df.columns}
    id_names = ["id_ref", "probe_id", "probeid", "probe", "id", "ilmn_id", "array_address_id", "feature_id"]
    value_names = ["value", "avg_signal", "signal", "signal_a", "mean_signal", "expr", "expression"]
    id_col = next((lower[n] for n in id_names if n in lower), None)
    value_col = next((lower[n] for n in value_names if n in lower), None)

    if id_col is None:
        id_col = df.columns[0]
    if value_col is None:
        numeric = df.drop(columns=id_col).apply(pd.to_numeric, errors="coerce")
        candidates = [c for c in numeric.columns if numeric[c].notna().sum() >= max(10, int(0.5 * len(df)))]
        if not candidates:
            raise ValueError(f"Nie znaleziono kolumny wartości w {path.name}")
        value_col = candidates[0]

    values = pd.to_numeric(df[value_col], errors="coerce")
    ids = df[id_col].astype(str).str.strip()
    out = pd.DataFrame({"id": ids, "value": values})
    out = out[(out["id"] != "") & (out["id"] != "nan") & out["value"].notna()]
    if len(out) < 10:
        raise ValueError(f"Za mało wartości ekspresji w {path.name}")
    out = out.groupby("id")["value"].mean().to_frame()
    sample_name = path.name
    for suffix in (".txt.gz", ".csv.gz", ".tsv.gz", ".gz", ".txt", ".csv", ".tsv"):
        if sample_name.endswith(suffix):
            sample_name = sample_name[:-len(suffix)]
            break
    return out.rename(columns={"value": sample_name})


def inspect_raw_archive(tar_path):
    """Print a compact inventory of the RAW archive so the next run identifies its real format."""
    print("\n=== GSE28688 RAW archive diagnostic ===")
    with tarfile.open(tar_path, "r") as tar:
        members = [m for m in tar.getmembers() if m.isfile()]
        print(f"Archive members: {len(members)}")
        for i, member in enumerate(members, 1):
            name = member.name
            size = member.size
            magic = ""
            try:
                f = tar.extractfile(member)
                head = f.read(64) if f is not None else b""
                magic = head[:16].hex(" ")
                text = head.decode("utf-8", errors="replace").replace("\r", " ").replace("\n", " ")[:100]
            except Exception as exc:
                text = f"<read error: {exc}>"
            print(f"[{i:02d}] {name} | {size:,} bytes | magic={magic}")
            print(f"     head={text}")
    print("=== end diagnostic ===\n")


def raw_archive_expression(extract):
    files = [f for f in extract.rglob("*") if f.is_file() and (f.suffix in {".txt", ".csv", ".tsv"} or f.name.endswith(".txt.gz") or f.name.endswith(".csv.gz") or f.name.endswith(".tsv.gz"))]
    if not files:
        raise ValueError("Archiwum RAW nie zawiera tabel tekstowych.")

    single_samples = []
    for path in files:
        try:
            one = raw_single_sample(path)
            if one.shape[0] >= 100:
                single_samples.append(one)
        except Exception:
            continue

    if len(single_samples) >= 2:
        expr = pd.concat(single_samples, axis=1, join="outer")
        expr.columns = [str(c) for c in expr.columns]
        return expr

    for path in files:
        try:
            expr = expression(read_table(path))
            if expr.shape[0] >= 10 and expr.shape[1] >= 2:
                return expr
        except Exception:
            continue

    raise ValueError("Nie znaleziono użytecznej macierzy ekspresji w GSE28688_RAW.tar")


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

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.boxplot([log[c].dropna() for c in log], tick_labels=log.columns, showfliers=False)
    ax.set_title(f"GSE28688 - {label} - distributions")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(out / "01_boxplot.png", dpi=250)
    plt.close(fig)

    qc = pd.DataFrame({"mean": expr.mean(), "median": expr.median(), "sd": expr.std(), "missing": expr.isna().sum()})
    qc.to_csv(out / "02_sample_QC.csv")

    corr = norm.corr()
    corr.to_csv(out / "03_sample_correlation.csv")
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(corr, vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(corr)))
    ax.set_yticklabels(corr.index)
    ax.set_title(f"GSE28688 - {label} - correlation")
    fig.colorbar(im, ax=ax, label="Pearson r")
    fig.tight_layout()
    fig.savefig(out / "04_sample_correlation.png", dpi=250)
    plt.close(fig)

    var = norm.var(axis=1).sort_values(ascending=False)
    var.to_csv(out / "05_feature_variance.csv", header=["variance"])
    top = var.head(min(2000, len(var))).index
    X = norm.loc[top].T.fillna(norm.loc[top].T.mean())
    C = X.to_numpy() - X.to_numpy().mean(axis=0)
    U, S, _ = np.linalg.svd(C, full_matrices=False)
    PC = U * S
    EV = S**2 / np.sum(S**2)
    n_pc = min(5, PC.shape[1])
    pd.DataFrame(PC[:, :n_pc], index=X.index, columns=[f"PC{i + 1}" for i in range(n_pc)]).to_csv(out / "06_PCA_coordinates.csv")

    if len(EV) >= 2:
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(PC[:, 0], PC[:, 1], s=90)
        for i, sample in enumerate(X.index):
            ax.annotate(sample, (PC[i, 0], PC[i, 1]), xytext=(6, 6), textcoords="offset points", fontsize=8)
        ax.set_xlabel(f"PC1 ({EV[0] * 100:.1f}%)")
        ax.set_ylabel(f"PC2 ({EV[1] * 100:.1f}%)")
        ax.set_title(f"GSE28688 - {label} - PCA")
        fig.tight_layout()
        fig.savefig(out / "07_PCA.png", dpi=250)
        plt.close(fig)

    features = var.head(min(50, len(var))).index
    z = norm.loc[features]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(10, max(6, len(features) * 0.18)))
    im = ax.imshow(z, aspect="auto")
    ax.set_xticks(range(len(z.columns)))
    ax.set_xticklabels(z.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=6)
    ax.set_title(f"GSE28688 - {label} - variable features")
    fig.colorbar(im, ax=ax, label="row z-score")
    fig.tight_layout()
    fig.savefig(out / "08_top_variable_features.png", dpi=250)
    plt.close(fig)

    report = [
        f"Dataset: GSE28688 ({label})",
        f"Features: {expr.shape[0]:,}",
        f"Samples: {expr.shape[1]:,}",
    ]
    if len(EV) >= 2:
        report += [
            f"PC1: {EV[0] * 100:.2f}%",
            f"PC2: {EV[1] * 100:.2f}%",
            f"PC1+PC2: {(EV[0] + EV[1]) * 100:.2f}%",
        ]
    (out / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")


p = DATA / "GSE28688_non-normalized.txt.gz"
if p.exists():
    explore(expression(read_table(p)), "non_normalized")

p = DATA / "GSE28688_RAW.tar"
if p.exists():
    inspect_raw_archive(p)
    extract = ROOT / "GSE28688_extracted"
    extract.mkdir(exist_ok=True)
    if not any(extract.rglob("*")):
        with tarfile.open(p, "r") as tar:
            tar.extractall(extract)
    try:
        explore(raw_archive_expression(extract), "RAW_archive")
    except ValueError as exc:
        print(f"RAW analysis not completed yet: {exc}")
        print("Powyzej jest diagnostyka archiwum RAW potrzebna do dopasowania parsera.")

print("GSE28688 exploration complete")
