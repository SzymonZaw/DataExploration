from pathlib import Path
import gzip
import tarfile
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)


def quantile_normalize(df):
    x = df.to_numpy(dtype=float)
    order = np.argsort(x, axis=0)
    sorted_x = np.sort(x, axis=0)
    means = sorted_x.mean(axis=1)
    out = np.empty_like(x)
    for j in range(x.shape[1]):
        out[order[:, j], j] = means
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def numeric_expression(df):
    """Find an expression matrix in a tabular file."""
    df = df.dropna(how="all").dropna(axis=1, how="all")
    # Prefer an explicit gene/probe identifier column.
    id_candidates = [c for c in df.columns if str(c).lower() in {
        "gene", "geneid", "gene_id", "symbol", "gene_symbol", "probe",
        "probename", "id", "feature", "featureid", "feature_id"
    }]
    if id_candidates:
        id_col = id_candidates[0]
    else:
        id_col = df.columns[0]

    numeric = df.drop(columns=[id_col]).apply(pd.to_numeric, errors="coerce")
    good = numeric.notna().sum(axis=0) > max(2, int(0.5 * len(numeric)))
    numeric = numeric.loc[:, good]
    if numeric.shape[1] < 2:
        raise ValueError("Nie znaleziono co najmniej dwóch kolumn liczbowych z ekspresją/counts.")

    ids = df[id_col].astype(str).str.strip()
    keep = ids.notna() & (ids != "") & (ids != "nan")
    numeric = numeric.loc[keep]
    ids = ids.loc[keep]
    numeric.index = ids
    numeric = numeric.groupby(level=0).mean()
    numeric = numeric.replace([np.inf, -np.inf], np.nan).dropna(how="all")
    return numeric


def read_table(path):
    """Read common CSV/TXT/GZ expression formats."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        sample = f.read(4096)

    sep = "," if sample.count(",") > sample.count("\t") else "\t"
    return pd.read_csv(path, sep=sep, compression="infer", low_memory=False)


def read_gse52052_tar(path, workdir):
    """Extract and read Agilent Feature Extraction TXT.GZ files."""
    extract_dir = workdir / "GSE52052_extracted"
    extract_dir.mkdir(exist_ok=True)
    files = sorted(extract_dir.rglob("*.txt.gz"))
    if not files:
        with tarfile.open(path, "r") as tar:
            tar.extractall(extract_dir)
        files = sorted(extract_dir.rglob("*.txt.gz"))
    if not files:
        raise ValueError("GSE52052_RAW.tar nie zawiera plików TXT.GZ.")

    signals = []
    for file in files:
        with gzip.open(file, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
        idx = next((i for i, line in enumerate(lines) if line.startswith("FEATURES\t")), None)
        if idx is None:
            continue
        cols = lines[idx].split("\t")[1:]
        rows = []
        for line in lines[idx + 1:]:
            if line.startswith("*"):
                break
            if line.startswith("DATA\t"):
                fields = line.split("\t")[1:]
                if len(fields) == len(cols):
                    rows.append(fields)
        df = pd.DataFrame(rows, columns=cols)
        if "ProbeName" not in df or "gProcessedSignal" not in df:
            continue
        if "ControlType" in df:
            ct = pd.to_numeric(df["ControlType"], errors="coerce")
            df = df[ct.fillna(0) == 0]
        df["gProcessedSignal"] = pd.to_numeric(df["gProcessedSignal"], errors="coerce")
        df = df.dropna(subset=["gProcessedSignal"])
        df["ProbeName"] = df["ProbeName"].astype(str).str.strip()
        df = df[(df["ProbeName"] != "") & (df["ProbeName"] != "nan")]
        df.loc[df["gProcessedSignal"] < 0, "gProcessedSignal"] = 0
        name = file.name[:-7]
        signals.append(df.groupby("ProbeName")["gProcessedSignal"].mean().rename(name))

    if not signals:
        raise ValueError("Nie udało się odczytać danych GSE52052.")
    return pd.concat(signals, axis=1, join="inner")


def read_gse28688_raw(path, workdir):
    """Extract GSE28688 RAW.tar and find the first usable expression table."""
    extract_dir = workdir / "GSE28688_extracted"
    extract_dir.mkdir(exist_ok=True)
    files = sorted(extract_dir.rglob("*"))
    if not files:
        with tarfile.open(path, "r") as tar:
            tar.extractall(extract_dir)
        files = sorted(p for p in extract_dir.rglob("*") if p.is_file())

    tables = [p for p in files if p.suffix in {".txt", ".csv", ".tsv"} or p.name.endswith(".txt.gz")]
    if not tables:
        raise ValueError("Nie znaleziono tabeli ekspresji w GSE28688_RAW.tar.")
    for p in tables:
        try:
            df = read_table(p)
            expr = numeric_expression(df)
            if expr.shape[1] >= 2 and expr.shape[0] >= 10:
                return expr
        except Exception:
            pass
    raise ValueError("Nie znaleziono użytecznej macierzy ekspresji w GSE28688_RAW.tar.")


def plot_box(expr, title, path):
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.boxplot([expr[c].dropna() for c in expr.columns], labels=expr.columns, showfliers=False)
    ax.set_title(title)
    ax.set_ylabel("Expression")
    ax.tick_params(axis="x", rotation=35)
    plt.tight_layout()
    plt.savefig(path, dpi=250)
    plt.close()


def analyze(name, expr, already_normalized=False):
    out = RESULTS / name
    out.mkdir(exist_ok=True)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.replace([np.inf, -np.inf], np.nan).dropna(how="all")
    expr = expr.dropna(axis=1, how="all")

    # Counts/intensities are displayed on log2 scale for comparable EDA.
    if (expr.min().min() < 0):
        raise ValueError("Macierz zawiera wartości ujemne.")
    log_expr = np.log2(expr + 1) if not already_normalized else np.log2(expr.clip(lower=0) + 1)
    log_expr.to_csv(out / "expression_log2.csv")
    expr.to_csv(out / "expression_input.csv")

    plot_box(log_expr, f"{name} - log2 expression", out / "01_boxplot.png")

    norm = log_expr if already_normalized else quantile_normalize(log_expr)
    norm.to_csv(out / "expression_normalized.csv")
    plot_box(norm, f"{name} - normalized log2 expression", out / "02_boxplot_normalized.png")

    corr = norm.corr()
    corr.to_csv(out / "03_sample_correlation.csv")

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(corr.values, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(corr.columns)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(corr.index)))
    ax.set_yticklabels(corr.index)
    ax.set_title(f"{name} - sample correlation")
    plt.colorbar(im, ax=ax, label="Pearson r")
    plt.tight_layout()
    plt.savefig(out / "04_sample_correlation.png", dpi=250)
    plt.close()

    variance = norm.var(axis=1).sort_values(ascending=False)
    variance.to_csv(out / "05_probe_variance.csv", header=["variance"])
    top = variance.head(min(2000, len(variance))).index
    X = norm.loc[top].T.fillna(norm.loc[top].T.mean())
    centered = X.to_numpy() - X.to_numpy().mean(axis=0)
    U, S, _ = np.linalg.svd(centered, full_matrices=False)
    pcs = U * S
    ev = S ** 2 / np.sum(S ** 2)
    pca = pd.DataFrame(pcs[:, :min(5, pcs.shape[1])], index=X.index,
                       columns=[f"PC{i+1}" for i in range(min(5, pcs.shape[1]))])
    pca.to_csv(out / "06_PCA_coordinates.csv")

    if len(ev) >= 2:
        fig, ax = plt.subplots(figsize=(8, 7))
        ax.scatter(pcs[:, 0], pcs[:, 1], s=90)
        for i, sample in enumerate(X.index):
            ax.annotate(sample, (pcs[i, 0], pcs[i, 1]), xytext=(6, 6), textcoords="offset points", fontsize=8)
        ax.set_xlabel(f"PC1 ({ev[0] * 100:.1f}%)")
        ax.set_ylabel(f"PC2 ({ev[1] * 100:.1f}%)")
        ax.set_title(f"{name} - PCA")
        plt.tight_layout()
        plt.savefig(out / "07_PCA.png", dpi=250)
        plt.close()

    top_h = variance.head(min(50, len(variance))).index
    h = norm.loc[top_h]
    z = h.sub(h.mean(axis=1), axis=0).div(h.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(10, max(6, len(top_h) * 0.18)))
    im = ax.imshow(z.values, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(z.columns)))
    ax.set_xticklabels(z.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(z.index)))
    ax.set_yticklabels(z.index, fontsize=6)
    ax.set_title(f"{name} - top variable features")
    plt.colorbar(im, ax=ax, label="row z-score")
    plt.tight_layout()
    plt.savefig(out / "08_top_variable_features.png", dpi=250)
    plt.close()

    report = [
        f"Dataset: {name}",
        f"Features: {expr.shape[0]}",
        f"Samples: {expr.shape[1]}",
        f"Input already normalized: {already_normalized}",
    ]
    if len(ev) >= 2:
        report += [f"PC1: {ev[0] * 100:.2f}%", f"PC2: {ev[1] * 100:.2f}%",
                   f"PC1+PC2: {(ev[0] + ev[1]) * 100:.2f}%"]
    (out / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")
    print(f"{name}: {expr.shape[0]:,} features x {expr.shape[1]} samples")


def main():
    print("DataExploration - EDA wszystkich plików w Data/")
    processed = set()

    # 1. GSE148158: normalized counts matrix
    p = DATA_DIR / "GSE148158_normalized_counts_MS2.csv.gz"
    if p.exists():
        expr = numeric_expression(read_table(p))
        analyze("GSE148158", expr, already_normalized=True)
        processed.add(p.name)

    # 2. GSE28688: non-normalized matrix
    p = DATA_DIR / "GSE28688_non-normalized.txt.gz"
    if p.exists():
        expr = numeric_expression(read_table(p))
        analyze("GSE28688_non_normalized", expr, already_normalized=False)
        processed.add(p.name)

    # 3. GSE28688: raw archive
    p = DATA_DIR / "GSE28688_RAW.tar"
    if p.exists():
        expr = read_gse28688_raw(p, ROOT)
        analyze("GSE28688_RAW", expr, already_normalized=False)
        processed.add(p.name)

    # 4. GSE52052: Agilent raw archive
    p = DATA_DIR / "GSE52052_RAW.tar"
    if p.exists():
        expr = read_gse52052_tar(p, ROOT)
        analyze("GSE52052_RAW", expr, already_normalized=False)
        processed.add(p.name)

    missing = [p.name for p in DATA_DIR.iterdir() if p.is_file() and p.name not in processed]
    print("\nGotowe.")
    print(f"Przetworzone pliki: {len(processed)}")
    if missing:
        print("Nieprzetworzone pliki:", ", ".join(missing))
    print(f"Wyniki: {RESULTS.resolve()}")


if __name__ == "__main__":
    main()
