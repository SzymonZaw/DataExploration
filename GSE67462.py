from pathlib import Path
import ftplib
import gzip
import shutil
import tarfile
import urllib.request
import ssl

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUT = ROOT / "results" / "GSE67462"
OUT.mkdir(parents=True, exist_ok=True)

ARCHIVE = DATA / "GSE67462_RAW.tar"
MATRIX_FILENAME = "GSE67462_series_matrix.txt.gz"
MATRIX_HTTPS_URL = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE67nnn/GSE67462/matrix/{MATRIX_FILENAME}"
MATRIX_FTP_DIR = "/geo/series/GSE67nnn/GSE67462/matrix"

SAMPLE_GROUPS = {
    "GSM1647454": "day0", "GSM1647455": "day0", "GSM1647456": "day1", "GSM1647457": "day1",
    "GSM1647458": "day3", "GSM1647459": "day3", "GSM1647460": "day5", "GSM1647461": "day5",
    "GSM1647462": "day7", "GSM1647463": "day7", "GSM1647464": "day11", "GSM1647465": "day11",
    "GSM1647466": "day15", "GSM1647467": "day15", "GSM1647468": "day18", "GSM1647469": "day18",
    "GSM1647470": "iPSC", "GSM1647471": "iPSC",
    "2ndMEF_at_day0_rep1": "day0", "2ndMEF_at_day0_rep2": "day0",
    "Reprogramming_cells_at_day1_rep1": "day1", "Reprogramming_cells_at_day1_rep2": "day1",
    "Reprogramming_cells_at_day3_rep1": "day3", "Reprogramming_cells_at_day3_rep2": "day3",
    "Reprogramming_cells_at_day5_rep1": "day5", "Reprogramming_cells_at_day5_rep2": "day5",
    "Reprogramming_cells_at_day7_rep1": "day7", "Reprogramming_cells_at_day7_rep2": "day7",
    "Reprogramming_cells_at_day11_rep1": "day11", "Reprogramming_cells_at_day11_rep2": "day11",
    "Reprogramming_cells_at_day15_rep1": "day15", "Reprogramming_cells_at_day15_rep2": "day15",
    "Reprogrammed_cells_at_day18_rep1": "day18", "Reprogrammed_cells_at_day18_rep2": "day18",
    "iPSC_rep1": "iPSC", "iPSC_rep2": "iPSC",
}


def group_for_sample(sample):
    s = str(sample).strip().strip('"')
    if s in SAMPLE_GROUPS:
        return SAMPLE_GROUPS[s]
    for key, group in SAMPLE_GROUPS.items():
        if key in s or s.startswith(key):
            return group
    return "unclassified"


def inspect_archive(path):
    files = []
    with tarfile.open(path, "r") as tar:
        for member in tar.getmembers():
            if member.isfile():
                files.append((member.name, member.size))
    report = [
        "Dataset: GSE67462",
        "RAW archive contains Affymetrix CEL files.",
        f"Archive members: {len(files)}",
    ]
    for name, size in files:
        report.append(f"{name}\t{size} bytes")
    (OUT / "RAW_archive_inventory.txt").write_text("\n".join(report), encoding="utf-8")
    print("=== GSE67462 RAW archive ===")
    print(f"CEL files: {len(files)}")
    for name, size in files:
        print(f"{name} | {size:,} bytes")
    print("=== end archive ===")


def validate_gzip(path):
    try:
        with gzip.open(path, "rb") as f:
            f.read(1024)
        return True
    except (OSError, EOFError):
        return False


def download_https(target, context=None):
    temp = target.with_suffix(target.suffix + ".part")
    if temp.exists():
        temp.unlink()
    request = urllib.request.Request(MATRIX_HTTPS_URL, headers={"User-Agent": "DataExploration-GSE67462/1.0"})
    with urllib.request.urlopen(request, context=context, timeout=120) as response, open(temp, "wb") as out:
        shutil.copyfileobj(response, out, length=1024 * 1024)
    if not validate_gzip(temp):
        temp.unlink(missing_ok=True)
        raise ValueError("Pobrany plik HTTPS nie jest poprawnym archiwum gzip.")
    temp.replace(target)


def download_ftp(target):
    temp = target.with_suffix(target.suffix + ".part")
    if temp.exists():
        temp.unlink()
    print("Trying NCBI FTP as a fallback...")
    with ftplib.FTP("ftp.ncbi.nlm.nih.gov", timeout=120) as ftp:
        ftp.login()
        ftp.cwd(MATRIX_FTP_DIR)
        with open(temp, "wb") as out:
            ftp.retrbinary(f"RETR {MATRIX_FILENAME}", out.write, blocksize=1024 * 1024)
    if not validate_gzip(temp):
        temp.unlink(missing_ok=True)
        raise ValueError("Pobrany plik FTP nie jest poprawnym archiwum gzip.")
    temp.replace(target)


def download_matrix():
    target = OUT / MATRIX_FILENAME
    if target.exists() and target.stat().st_size > 0 and validate_gzip(target):
        print(f"Using existing valid Series Matrix: {target.name} ({target.stat().st_size:,} bytes)")
        return target
    if target.exists():
        print("Existing Series Matrix is incomplete or invalid; downloading it again...")
        target.unlink()

    print("Downloading GEO Series Matrix for processed expression exploration...")
    errors = []

    try:
        download_https(target)
        print("Download successful via standard HTTPS.")
        return target
    except Exception as e:
        errors.append(f"HTTPS: {e}")
        print(f"Standard HTTPS download failed: {e}")

    try:
        print("Retrying HTTPS with a compatibility SSL context...")
        download_https(target, ssl._create_unverified_context())
        print("Download successful via compatibility SSL context.")
        return target
    except Exception as e:
        errors.append(f"HTTPS compatibility: {e}")
        print(f"Compatibility HTTPS download failed: {e}")

    try:
        download_ftp(target)
        print("Download successful via NCBI FTP.")
        return target
    except Exception as e:
        errors.append(f"FTP: {e}")
        print(f"NCBI FTP download failed: {e}")

    raise RuntimeError("Nie udało się pobrać poprawnego GEO Series Matrix. " + " | ".join(errors))


def read_series_matrix(path):
    rows = []
    columns = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
        else:
            raise ValueError("Nie znaleziono !series_matrix_table_begin w Series Matrix.")

        for line in f:
            if line.startswith("!series_matrix_table_end"):
                break
            line = line.rstrip("\r\n")
            if not line:
                continue
            parts = line.split("\t")
            if columns is None:
                columns = [p.strip('"') for p in parts]
            else:
                rows.append(parts)
        else:
            raise ValueError("Nie znaleziono !series_matrix_table_end w Series Matrix.")

    if columns is None or not rows:
        raise ValueError("Nie udało się odczytać tabeli Series Matrix GSE67462.")

    df = pd.DataFrame(rows, columns=columns)
    df = df.rename(columns={df.columns[0]: "ID_REF"})
    for c in df.columns[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.set_index("ID_REF")
    df = df.dropna(how="all")
    print(f"Series Matrix parsed successfully: {df.shape[0]:,} features x {df.shape[1]:,} samples")
    return df


def explore(expr):
    expr = expr.replace([np.inf, -np.inf], np.nan).dropna(how="all").clip(lower=0)
    expr = expr.groupby(level=0).mean()
    expr.to_csv(OUT / "01_expression_processed.csv")

    metadata = pd.DataFrame({"sample": expr.columns})
    metadata["group"] = metadata["sample"].map(group_for_sample)
    metadata.to_csv(OUT / "02_sample_metadata.csv", index=False)

    log = expr.copy()
    log.to_csv(OUT / "03_expression_for_EDA.csv")

    fig, ax = plt.subplots(figsize=(13, 6))
    ax.boxplot([log[c].dropna() for c in log], tick_labels=log.columns, showfliers=False)
    ax.set_title("GSE67462 - processed expression distributions")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(OUT / "04_boxplot.png", dpi=250); plt.close(fig)

    qc = pd.DataFrame({"mean": expr.mean(), "median": expr.median(), "sd": expr.std(), "missing": expr.isna().sum()})
    qc["group"] = qc.index.map(group_for_sample)
    qc.to_csv(OUT / "05_sample_QC.csv")

    corr = expr.corr()
    corr.to_csv(OUT / "06_sample_correlation.csv")
    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(corr, vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr))); ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(corr))); ax.set_yticklabels(corr.index)
    ax.set_title("GSE67462 - sample correlation")
    fig.colorbar(im, ax=ax, label="Pearson r")
    fig.tight_layout(); fig.savefig(OUT / "07_sample_correlation.png", dpi=250); plt.close(fig)

    var = expr.var(axis=1).sort_values(ascending=False)
    var.to_csv(OUT / "08_feature_variance.csv", header=["variance"])
    top = var.head(min(3000, len(var))).index
    X = expr.loc[top].T.fillna(expr.loc[top].T.mean())
    C = X.to_numpy(float) - X.to_numpy(float).mean(axis=0)
    U, S, _ = np.linalg.svd(C, full_matrices=False)
    PC = U * S
    EV = S**2 / np.sum(S**2)
    n_pc = min(5, PC.shape[1])
    coords = pd.DataFrame(PC[:, :n_pc], index=X.index, columns=[f"PC{i+1}" for i in range(n_pc)])
    coords["group"] = coords.index.map(group_for_sample)
    coords.to_csv(OUT / "09_PCA_coordinates.csv")

    if len(EV) >= 2:
        fig, ax = plt.subplots(figsize=(10, 7))
        for group, idx in coords.groupby("group").groups.items():
            positions = [X.index.get_loc(s) for s in idx]
            ax.scatter(PC[positions, 0], PC[positions, 1], s=90, label=group)
            for p in positions:
                ax.annotate(X.index[p], (PC[p, 0], PC[p, 1]), xytext=(5, 5), textcoords="offset points", fontsize=7)
        ax.set_xlabel(f"PC1 ({EV[0]*100:.1f}%)"); ax.set_ylabel(f"PC2 ({EV[1]*100:.1f}%)")
        ax.set_title("GSE67462 - PCA by reprogramming stage")
        ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig(OUT / "10_PCA_by_group.png", dpi=250); plt.close(fig)

    features = var.head(min(50, len(var))).index
    z = expr.loc[features]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(10, max(6, len(features) * 0.18)))
    im = ax.imshow(z, aspect="auto")
    ax.set_xticks(range(len(z.columns))); ax.set_xticklabels(z.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(z.index))); ax.set_yticklabels(z.index, fontsize=6)
    ax.set_title("GSE67462 - top variable features")
    fig.colorbar(im, ax=ax, label="row z-score")
    fig.tight_layout(); fig.savefig(OUT / "11_top_variable_features.png", dpi=250); plt.close(fig)

    counts = metadata["group"].value_counts().sort_index()
    report = [
        "Dataset: GSE67462",
        "Title: Expression data from OSKM-mediated 2nd reprogramming cells and the corresponding iPS cell line",
        "Organism: Mus musculus",
        "Experiment: expression profiling by array",
        "Platform: Affymetrix Mouse Gene 1.0 ST Array (GPL19972)",
        "RAW archive: Affymetrix CEL files",
        f"Features in processed table: {expr.shape[0]:,}",
        f"Samples: {expr.shape[1]:,}",
        "Processed values were taken from the GEO Series Matrix; the GEO sample table describes them as transformed values.",
        "This script performs exploratory QC, correlation, variance and PCA; it does not re-normalize raw CEL intensities.",
        "Sample groups:",
    ]
    report += [f"  {g}: {n}" for g, n in counts.items()]
    if len(EV) >= 2:
        report += [f"PC1: {EV[0]*100:.2f}%", f"PC2: {EV[1]*100:.2f}%", f"PC1+PC2: {(EV[0]+EV[1])*100:.2f}%"]
    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")
    print("GSE67462 EDA results generated successfully.")


if ARCHIVE.exists():
    inspect_archive(ARCHIVE)
else:
    print("Brak GSE67462_RAW.tar")

try:
    matrix = download_matrix()
    explore(read_series_matrix(matrix))
except Exception as e:
    print(f"Processed expression exploration failed: {e}")
    print("RAW CEL inventory was still generated above.")

print("GSE67462 exploration complete")
