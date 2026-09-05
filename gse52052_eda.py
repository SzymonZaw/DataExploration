from pathlib import Path
import gzip
import io
import tarfile
import urllib.request
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ================================================================
# GSE52052 - eksploracja danych mikromacierzy Agilent
#
# Źródło danych: Data/GSE52052_RAW.tar
#
# Pipeline:
#   1. rozpakowanie/odczyt plików Agilent TXT.GZ
#   2. QC sygnału i flag Agilent
#   3. macierz ProbeName x sample
#   4. log2(x + 1)
#   5. quantile normalization
#   6. boxploty i rozkłady
#   7. korelacja próbek
#   8. PCA
#   9. hierarchical clustering
#  10. heatmapa najbardziej zmiennych sond
#  11. adnotacja GPL14550, jeśli dostępna
#
# Skrypt NIE wykonuje inferencyjnego differential expression,
# ponieważ GSE52052 ma po jednej próbce na warunek.
# ================================================================

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data"
ARCHIVE = DATA_DIR / "GSE52052_RAW.tar"
EXTRACT_DIR = ROOT / "GSE52052_extracted"
RESULTS = ROOT / "results_GSE52052"
RESULTS.mkdir(exist_ok=True)

N_PCA_PROBES = 2000
N_HEATMAP_PROBES = 50
USE_ONLY_CONTROLTYPE_0 = True


def sample_name(filename):
    name = Path(filename).name
    if name.endswith(".txt.gz"):
        name = name[:-7]
    return name.split("_US")[0] if "_US" in name else name


def extract_archive():
    EXTRACT_DIR.mkdir(exist_ok=True)
    existing = list(EXTRACT_DIR.rglob("*.txt.gz"))
    if existing:
        return sorted(existing)
    if not ARCHIVE.exists():
        raise FileNotFoundError(f"Brak archiwum: {ARCHIVE}")
    print(f"Rozpakowuję: {ARCHIVE}")
    with tarfile.open(ARCHIVE, "r") as tar:
        tar.extractall(EXTRACT_DIR)
    return sorted(EXTRACT_DIR.rglob("*.txt.gz"))


def read_agilent(path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    idx = next((i for i, x in enumerate(lines) if x.startswith("FEATURES\t")), None)
    if idx is None:
        raise ValueError(f"Nie znaleziono sekcji FEATURES: {path.name}")
    cols = lines[idx].split("\t")[1:]
    rows = []
    for line in lines[idx + 1:]:
        if line.startswith("*"):
            break
        if line.startswith("DATA\t"):
            p = line.split("\t")[1:]
            if len(p) == len(cols):
                rows.append(p)
    return pd.DataFrame(rows, columns=cols)


def prepare(df, name):
    if "ProbeName" not in df or "gProcessedSignal" not in df:
        raise ValueError(f"Brak wymaganych kolumn w {name}")
    d = df.copy()
    d["ProbeName"] = d["ProbeName"].astype(str).str.strip()
    d["gProcessedSignal"] = pd.to_numeric(d["gProcessedSignal"], errors="coerce")
    if USE_ONLY_CONTROLTYPE_0 and "ControlType" in d:
        d["ControlType"] = pd.to_numeric(d["ControlType"], errors="coerce")
        d = d[d["ControlType"].fillna(0) == 0]
    d = d[(d["ProbeName"] != "") & (d["ProbeName"] != "nan")]
    d = d.dropna(subset=["gProcessedSignal"])
    d.loc[d["gProcessedSignal"] < 0, "gProcessedSignal"] = 0
    duplicates = int(d["ProbeName"].duplicated().sum())
    signal = d.groupby("ProbeName")["gProcessedSignal"].mean().rename(name)
    return signal, d, duplicates


def quantile_normalize(df):
    x = df.to_numpy(dtype=float)
    order = np.argsort(x, axis=0)
    sorted_x = np.sort(x, axis=0)
    rank_mean = sorted_x.mean(axis=1)
    out = np.zeros_like(x)
    for j in range(x.shape[1]):
        out[order[:, j], j] = rank_mean
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def save_heatmap(data, title, filename, ylabel=""):
    fig, ax = plt.subplots(figsize=(10, max(5, min(16, 0.28 * len(data) + 3))))
    im = ax.imshow(data.values, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(data.columns)))
    ax.set_xticklabels(data.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(data.index)))
    ax.set_yticklabels(data.index, fontsize=7)
    ax.set_title(title)
    if ylabel:
        ax.set_ylabel(ylabel)
    plt.colorbar(im, ax=ax)
    plt.tight_layout()
    plt.savefig(RESULTS / filename, dpi=300)
    plt.close()


def main():
    print("=" * 80)
    print("GSE52052 - EKSPLORACJA DANYCH")
    print("=" * 80)

    files = extract_archive()
    print(f"\nZnaleziono plików: {len(files)}")

    signals = []
    qc_rows = []

    for path in files:
        name = sample_name(path.name)
        print(f"\nCzytam: {path.name}")
        df = read_agilent(path)
        signal, prepared, duplicates = prepare(df, name)
        print(f"  feature'ów: {len(df):,}")
        print(f"  kolumn: {len(df.columns)}")
        print(f"  ProbeName po filtracji: {signal.size:,}")
        print(f"  duplikaty ProbeName: {duplicates:,}")
        signals.append(signal)

        row = {
            "sample": name,
            "features": len(df),
            "probes_after_filter": signal.size,
            "duplicate_probe_rows": duplicates,
            "missing_signal": int(df["gProcessedSignal"].isna().sum()),
            "negative_signal": int((pd.to_numeric(df["gProcessedSignal"], errors="coerce") < 0).sum()),
            "saturated": int((prepared.get("gIsSaturated", pd.Series(dtype=str)).astype(str).isin(["1", "True", "TRUE"]).sum()) if "gIsSaturated" in prepared else 0),
            "signal_median": signal.median(),
            "signal_mean": signal.mean(),
            "signal_sd": signal.std(),
            "signal_max": signal.max(),
        }
        qc_rows.append(row)

    expr = pd.concat(signals, axis=1, join="inner")
    expr.index.name = "ProbeName"
    print(f"\nWspólna macierz: {expr.shape[0]:,} sond x {expr.shape[1]} próbek")
    expr.to_csv(RESULTS / "01_expression_raw_gProcessedSignal.csv")

    qc = pd.DataFrame(qc_rows).set_index("sample")
    qc.to_csv(RESULTS / "02_sample_QC.csv")
    print("\nQC:")
    print(qc.round(3).to_string())

    # log2
    log_expr = np.log2(expr.clip(lower=0) + 1)
    log_expr.to_csv(RESULTS / "03_expression_log2.csv")

    # boxplot before normalization
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.boxplot([log_expr[c].dropna() for c in log_expr.columns], labels=log_expr.columns, showfliers=False)
    ax.set_title("GSE52052 - log2 expression przed normalizacją")
    ax.set_ylabel("log2(gProcessedSignal + 1)")
    ax.tick_params(axis="x", rotation=35)
    plt.tight_layout()
    plt.savefig(RESULTS / "04_boxplot_before_normalization.png", dpi=300)
    plt.close()

    # quantile normalization
    print("\nWykonuję quantile normalization...")
    norm = quantile_normalize(log_expr)
    norm.to_csv(RESULTS / "05_expression_log2_quantile_normalized.csv")

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.boxplot([norm[c].dropna() for c in norm.columns], labels=norm.columns, showfliers=False)
    ax.set_title("GSE52052 - po quantile normalization")
    ax.set_ylabel("normalized log2 expression")
    ax.tick_params(axis="x", rotation=35)
    plt.tight_layout()
    plt.savefig(RESULTS / "06_boxplot_after_normalization.png", dpi=300)
    plt.close()

    # distributions
    fig, ax = plt.subplots(figsize=(12, 7))
    for c in norm.columns:
        ax.hist(norm[c].dropna(), bins=100, density=True, alpha=0.3, label=c)
    ax.set_title("Rozkłady ekspresji po normalizacji")
    ax.set_xlabel("log2 normalized expression")
    ax.set_ylabel("Density")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(RESULTS / "07_expression_distributions.png", dpi=300)
    plt.close()

    # correlation
    corr = norm.corr(method="pearson")
    corr.to_csv(RESULTS / "08_sample_correlation.csv")
    save_heatmap(corr, "Korelacja próbek", "09_sample_correlation_heatmap.png")

    # PCA on top variable probes
    variance = norm.var(axis=1).sort_values(ascending=False)
    variance.rename("variance").to_csv(RESULTS / "10_probe_variance.csv")
    top = variance.head(min(N_PCA_PROBES, len(variance))).index
    X = norm.loc[top].T
    X = X.fillna(X.mean(axis=0))
    centered = X.to_numpy() - X.to_numpy().mean(axis=0)
    U, S, _ = np.linalg.svd(centered, full_matrices=False)
    pcs = U * S
    ev = S**2 / np.sum(S**2)
    pca = pd.DataFrame(pcs[:, :min(5, pcs.shape[1])], index=X.index,
                       columns=[f"PC{i+1}" for i in range(min(5, pcs.shape[1]))])
    pca["sample"] = pca.index
    pca.to_csv(RESULTS / "11_PCA_coordinates.csv", index=False)

    print("\nPCA explained variance:")
    print(f"PC1: {ev[0] * 100:.2f}%")
    if len(ev) > 1:
        print(f"PC2: {ev[1] * 100:.2f}%")

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(pcs[:, 0], pcs[:, 1], s=100)
    for i, name in enumerate(X.index):
        ax.annotate(name, (pcs[i, 0], pcs[i, 1]), xytext=(7, 7), textcoords="offset points", fontsize=9)
    ax.set_xlabel(f"PC1 ({ev[0] * 100:.2f}%)")
    ax.set_ylabel(f"PC2 ({ev[1] * 100:.2f}%)")
    ax.set_title("PCA - GSE52052 po quantile normalization")
    ax.axhline(0, linewidth=0.7)
    ax.axvline(0, linewidth=0.7)
    plt.tight_layout()
    plt.savefig(RESULTS / "12_PCA.png", dpi=300)
    plt.close()

    # hierarchical clustering
    try:
        from scipy.cluster.hierarchy import linkage, dendrogram
        from scipy.spatial.distance import squareform
        distance = 1 - corr
        np.fill_diagonal(distance.values, 0)
        Z = linkage(squareform(distance.values, checks=False), method="average")
        fig, ax = plt.subplots(figsize=(14, 7))
        dendrogram(Z, labels=corr.columns, leaf_rotation=35, ax=ax)
        ax.set_title("Hierarchical clustering próbek")
        ax.set_ylabel("Distance = 1 - Pearson r")
        plt.tight_layout()
        plt.savefig(RESULTS / "13_hierarchical_clustering.png", dpi=300)
        plt.close()
    except ImportError:
        print("Brak scipy - pomijam hierarchical clustering.")

    # top variable heatmap, row z-score
    top_h = variance.head(min(N_HEATMAP_PROBES, len(variance))).index
    h = norm.loc[top_h]
    z = h.sub(h.mean(axis=1), axis=0).div(h.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    save_heatmap(z, f"Top {len(top_h)} najbardziej zmiennych sond", "14_top_variable_probes_heatmap.png")

    # Save simple report
    report = []
    report.append("GSE52052 - EDA REPORT")
    report.append(f"Samples: {expr.shape[1]}")
    report.append(f"Common probes: {expr.shape[0]}")
    report.append(f"PCA probes: {len(top)}")
    report.append(f"PC1: {ev[0] * 100:.2f}%")
    report.append(f"PC2: {ev[1] * 100:.2f}%")
    report.append(f"PC1+PC2: {(ev[0] + ev[1]) * 100:.2f}%")
    report.append("")
    report.append("Sample correlation:")
    report.append(corr.round(4).to_string())
    report.append("")
    report.append("NOTE: GSE52052 has one sample per experimental condition; no inferential DE p-values are calculated.")
    (RESULTS / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")

    print("\n" + "=" * 80)
    print("ANALIZA ZAKOŃCZONA")
    print("=" * 80)
    print(f"\nWyniki: {RESULTS.resolve()}")
    print("\nUtworzone główne pliki:")
    for p in sorted(RESULTS.iterdir()):
        print(f"  {p.name}")


if __name__ == "__main__":
    main()
