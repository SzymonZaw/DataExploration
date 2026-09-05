"""Exploratory preprocessing for GSE297234 single-cell RNA-seq.

Goal
----
Create a compact sample-level representation that can later be used by
Dynamics.py without pretending that single-cell PCA axes are directly
comparable to bulk/microarray PCA axes.

The script prefers the GEO RAW tar archive containing 10x Genomics H5 files.
It reads the filtered feature-barcode matrices with h5py, creates a
sample-level pseudobulk expression profile, log1p-transforms it, and performs
sample-level PCA. The large Seurat RDS files are intentionally not required.

Expected local input
--------------------
Data/GSE297234_RAW.tar

The GEO series contains 8 human fibroblast samples: aged GM00731 and young
GM23815, each at days 0, 3, 7 and 10 of Sendai OSKM reprogramming.

Important scientific limitation
--------------------------------
This is NOT a cell-level trajectory model. It is a sample-level exploratory
representation of a scRNA-seq dataset. The resulting PCA coordinates are not
assumed to be comparable with PCA coordinates from the other GEO studies.
"""

from pathlib import Path
import re
import tarfile

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUT = ROOT / "results" / "GSE297234"
OUT.mkdir(parents=True, exist_ok=True)

ARCHIVE = DATA / "GSE297234_RAW.tar"

# Known GEO sample design. The mapping is deliberately explicit so that
# biological time is available for the future dynamics stage.
SAMPLE_INFO = {
    "GSM8986586": {"donor": "GM00731", "age_group": "aged", "day": 0},
    "GSM8986587": {"donor": "GM00731", "age_group": "aged", "day": 3},
    "GSM8986588": {"donor": "GM00731", "age_group": "aged", "day": 7},
    "GSM8986589": {"donor": "GM00731", "age_group": "aged", "day": 10},
    "GSM8986590": {"donor": "GM23815", "age_group": "young", "day": 0},
    "GSM8986591": {"donor": "GM23815", "age_group": "young", "day": 3},
    "GSM8986592": {"donor": "GM23815", "age_group": "young", "day": 7},
    "GSM8986593": {"donor": "GM23815", "age_group": "young", "day": 10},
}


def infer_gsm(name):
    match = re.search(r"GSM\d+", name)
    return match.group(0) if match else None


def inspect_archive():
    if not ARCHIVE.exists():
        raise FileNotFoundError(
            f"Nie znaleziono {ARCHIVE}. Umieść GSE297234_RAW.tar w katalogu Data."
        )

    members = []
    with tarfile.open(ARCHIVE, "r") as tar:
        for member in tar.getmembers():
            if member.isfile():
                members.append((member.name, member.size))

    report = [
        "Dataset: GSE297234",
        "RAW archive inventory",
        f"Files: {len(members)}",
        "",
    ]
    report.extend(f"{name}\t{size} bytes" for name, size in members)
    (OUT / "RAW_archive_inventory.txt").write_text("\n".join(report), encoding="utf-8")

    print("=== GSE297234 RAW archive ===")
    for name, size in members:
        print(f"{name} | {size:,} bytes")
    print("=== end archive ===")
    return members


def extract_h5_files(members):
    extracted = []
    with tarfile.open(ARCHIVE, "r") as tar:
        for name, _ in members:
            if not name.lower().endswith(".h5"):
                continue
            gsm = infer_gsm(name)
            if gsm not in SAMPLE_INFO:
                print(f"Skipping H5 without known sample mapping: {name}")
                continue
            target = OUT / Path(name).name
            if not target.exists():
                print(f"Extracting {Path(name).name} ...")
                with tar.extractfile(name) as src, open(target, "wb") as dst:
                    if src is None:
                        raise RuntimeError(f"Could not read archive member: {name}")
                    while True:
                        chunk = src.read(1024 * 1024)
                        if not chunk:
                            break
                        dst.write(chunk)
            extracted.append((gsm, target))
    return extracted


def read_10x_h5(path):
    """Read a standard 10x filtered_feature_bc_matrix H5 file.

    Returns genes, cells and a dense sample-level count vector. We sum counts
    over cells immediately, avoiding construction of a huge dense matrix.
    """
    with h5py.File(path, "r") as h5:
        # Standard 10x layout: matrix/{data,indices,indptr,shape,...}
        matrix = h5["matrix"]
        data = np.asarray(matrix["data"], dtype=np.float64)
        indices = np.asarray(matrix["indices"], dtype=np.int64)
        indptr = np.asarray(matrix["indptr"], dtype=np.int64)
        shape = tuple(np.asarray(matrix["shape"], dtype=np.int64))

        features = matrix["features"]
        if "name" in features:
            genes = [x.decode() if isinstance(x, bytes) else str(x) for x in features["name"][:]]
        elif "gene_names" in features:
            genes = [x.decode() if isinstance(x, bytes) else str(x) for x in features["gene_names"][:]]
        else:
            genes = [str(x) for x in range(shape[0])]

    # 10x stores the sparse matrix in CSC-like form. Sum each row without
    # importing scipy: accumulate each nonzero value to its feature index.
    counts = np.bincount(indices, weights=data, minlength=shape[0])
    return genes, counts, shape


def build_sample_matrix(extracted):
    sample_vectors = []
    reference_genes = None

    for gsm, path in extracted:
        print(f"Reading {gsm}: {path.name}")
        genes, counts, shape = read_10x_h5(path)
        print(f"  matrix: {shape[0]:,} genes x {shape[1]:,} cells")

        if reference_genes is None:
            reference_genes = genes
        elif genes != reference_genes:
            raise ValueError(
                f"Feature order differs between samples; cannot safely aggregate {gsm}."
            )
        sample_vectors.append(pd.Series(counts, index=reference_genes, name=gsm))

    matrix = pd.concat(sample_vectors, axis=1)
    matrix = matrix.groupby(level=0).sum()
    matrix.to_csv(OUT / "01_sample_level_counts.csv")
    return matrix


def create_metadata(columns):
    rows = []
    for sample in columns:
        info = SAMPLE_INFO[sample]
        rows.append({
            "sample": sample,
            "donor": info["donor"],
            "age_group": info["age_group"],
            "day": info["day"],
            "treatment": "Sendai_OSKM",
            "cell_type": "fibroblast",
        })
    metadata = pd.DataFrame(rows).sort_values(["donor", "day"])
    metadata.to_csv(OUT / "02_sample_metadata.csv", index=False)
    return metadata


def run_eda(counts, metadata):
    # Pseudobulk CPM followed by log1p. This is only for exploratory
    # sample-level geometry, not differential-expression testing.
    library_size = counts.sum(axis=0)
    cpm = counts.div(library_size.replace(0, np.nan), axis=1) * 1e6
    log_expr = np.log1p(cpm)
    log_expr.to_csv(OUT / "03_log1p_CPM_sample_expression.csv")

    qc = pd.DataFrame({
        "library_size": library_size,
        "detected_features": (counts > 0).sum(axis=0),
        "mean": counts.mean(axis=0),
        "sd": counts.std(axis=0),
    })
    qc["sample"] = qc.index
    qc["day"] = qc["sample"].map(metadata.set_index("sample")["day"])
    qc["donor"] = qc["sample"].map(metadata.set_index("sample")["donor"])
    qc.to_csv(OUT / "04_sample_QC.csv", index=False)

    corr = log_expr.corr()
    corr.to_csv(OUT / "05_sample_correlation.csv")

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(corr, vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr)))
    ax.set_xticklabels(corr.columns, rotation=45, ha="right")
    ax.set_yticks(range(len(corr)))
    ax.set_yticklabels(corr.index)
    ax.set_title("GSE297234 - sample correlation")
    fig.colorbar(im, ax=ax, label="Pearson r")
    fig.tight_layout()
    fig.savefig(OUT / "06_sample_correlation.png", dpi=250)
    plt.close(fig)

    variance = log_expr.var(axis=1).sort_values(ascending=False)
    variance.to_csv(OUT / "07_feature_variance.csv", header=["variance"])

    top = variance.head(min(3000, len(variance))).index
    X = log_expr.loc[top].T.fillna(0.0)
    X_np = X.to_numpy(float)
    X_np -= X_np.mean(axis=0, keepdims=True)
    U, S, _ = np.linalg.svd(X_np, full_matrices=False)
    pc = U * S
    explained = S**2 / max(np.sum(S**2), 1e-12)
    n_pc = min(5, pc.shape[1])
    coords = pd.DataFrame(
        pc[:, :n_pc],
        index=X.index,
        columns=[f"PC{i + 1}" for i in range(n_pc)],
    )
    coords = coords.reset_index().rename(columns={"index": "sample"})
    coords = coords.merge(metadata, on="sample", how="left")
    coords.to_csv(OUT / "08_PCA_coordinates.csv", index=False)

    if n_pc >= 2:
        fig, ax = plt.subplots(figsize=(9, 7))
        for donor, subset in coords.groupby("donor"):
            ax.plot(subset["PC1"], subset["PC2"], marker="o", label=donor)
            for _, row in subset.iterrows():
                ax.annotate(
                    f"d{int(row['day'])}",
                    (row["PC1"], row["PC2"]),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=8,
                )
        ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
        ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
        ax.set_title("GSE297234 - sample-level PCA")
        ax.legend()
        fig.tight_layout()
        fig.savefig(OUT / "09_PCA_by_donor_and_time.png", dpi=250)
        plt.close(fig)

    report = [
        "Dataset: GSE297234",
        "Title: Prevalent mesenchymal drift in aging and disease is reversed by partial reprogramming [scRNA-seq]",
        "Organism: Homo sapiens",
        "Experiment: 10x Genomics single-cell RNA-seq",
        "Platform: Illumina NovaSeq 6000 / GPL24676",
        "Design: aged GM00731 and young GM23815 fibroblasts, Sendai OSKM, days 0/3/7/10",
        f"Samples processed: {len(metadata)}",
        f"Genes in sample-level matrix: {counts.shape[0]:,}",
        "Representation: cell-level counts were summed within each GEO sample (pseudobulk), then CPM and log1p were used for exploratory sample-level PCA.",
        "This is not a cell-level trajectory and is not a differential-expression pipeline.",
        "PCA coordinates are dataset-local and must not be compared numerically with PCA axes from other studies.",
        "",
        "Biological time points:",
        "  day 0, day 3, day 7, day 10 for each donor",
    ]
    if len(explained) >= 2:
        report.extend([
            f"PC1 explained variance: {explained[0] * 100:.2f}%",
            f"PC2 explained variance: {explained[1] * 100:.2f}%",
            f"PC1+PC2 explained variance: {(explained[0] + explained[1]) * 100:.2f}%",
        ])
    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")


def main():
    members = inspect_archive()
    extracted = extract_h5_files(members)
    if len(extracted) != len(SAMPLE_INFO):
        found = sorted(gsm for gsm, _ in extracted)
        expected = sorted(SAMPLE_INFO)
        raise RuntimeError(
            "Nie znaleziono kompletu 8 macierzy H5.\n"
            f"Expected: {expected}\nFound: {found}\n"
            "Sprawdź RAW_archive_inventory.txt."
        )

    counts = build_sample_matrix(extracted)
    metadata = create_metadata(counts.columns)
    run_eda(counts, metadata)
    print(f"GSE297234 EDA results written to: {OUT}")


if __name__ == "__main__":
    main()
