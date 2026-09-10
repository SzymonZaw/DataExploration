from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

SCATAC_DEFAULT = Path("Data/GSE242421/scATAC.zip")
INTEGRATION_DEFAULT = Path("Data/GSE242421/scATAC_scRNA_integration.zip")
OUT_DEFAULT = Path("results/Dynamics/z4_gse242421_gene_activity")
EXPECTED_SAMPLES = ["D0", "D2", "D4", "D6", "D8", "D10", "D12", "D14", "iPSC"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def member_by_basename(zf: zipfile.ZipFile, basename: str) -> str:
    matches = [n for n in zf.namelist() if Path(n).name == basename]
    if not matches:
        raise FileNotFoundError(f"{basename} not found in archive")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple archive members named {basename}: {matches[:10]}")
    return matches[0]


def read_cells(zf: zipfile.ZipFile, member: str) -> pd.DataFrame:
    with zf.open(member) as fh:
        df = pd.read_csv(fh, sep="\t")
    required = {"barcode", "sample"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"cells.tsv missing columns: {sorted(missing)}")
    return df[["barcode", "sample"]].copy()


def peak_id(chrom: str, start: object, end: object) -> str:
    return f"{chrom}:{int(start)}-{int(end)}"


def read_peaks(zf: zipfile.ZipFile, member: str) -> list[str]:
    rows: list[str] = []
    with zf.open(member) as fh:
        for raw in fh:
            if not raw.strip():
                continue
            parts = raw.decode("utf-8", errors="replace").rstrip("\n\r").split("\t")
            if len(parts) < 3:
                raise ValueError("peaks.bed has fewer than 3 columns")
            rows.append(peak_id(parts[0], parts[1], parts[2]))
    return rows


def load_mtx_from_zip(zf: zipfile.ZipFile, member: str) -> sparse.csc_matrix:
    # scipy.io.mmread accepts a binary file-like object. Convert to CSC because
    # columns correspond to cells and we aggregate columns by sample.
    with zf.open(member) as fh:
        matrix = mmread(fh)
    return sparse.csc_matrix(matrix, dtype=np.float64)


def normalize_header(name: str) -> str:
    return "".join(ch.lower() for ch in name if ch.isalnum())


def choose_column(columns: list[str], aliases: list[str], required: bool = True) -> str | None:
    norm = {normalize_header(c): c for c in columns}
    for alias in aliases:
        if normalize_header(alias) in norm:
            return norm[normalize_header(alias)]
    if required:
        raise ValueError(f"Could not identify required column from {aliases}; columns={columns}")
    return None


def read_peak_gene_links(
    zf: zipfile.ZipFile,
    member: str,
    min_abs_correlation: float,
) -> tuple[pd.DataFrame, dict]:
    with zf.open(member) as fh:
        links = pd.read_csv(fh, sep="\t")

    columns = list(links.columns)
    peak_col = choose_column(columns, ["peak", "peak_id", "peakid", "region", "element"])
    gene_col = choose_column(columns, ["gene", "gene_symbol", "genesymbol", "target_gene", "target"])
    corr_col = choose_column(
        columns,
        ["correlation", "corr", "pearson", "spearman", "r"],
        required=False,
    )

    original_n = len(links)
    if corr_col is not None:
        corr = pd.to_numeric(links[corr_col], errors="coerce")
        valid_corr = corr.notna()
        links = links.loc[valid_corr & (corr.abs() >= min_abs_correlation)].copy()
        links["_weight"] = 1.0
        correlation_filter = True
    else:
        links = links.copy()
        links["_weight"] = 1.0
        correlation_filter = False

    links["_peak"] = links[peak_col].astype(str).str.strip()
    links["_gene"] = links[gene_col].astype(str).str.strip()
    links = links[(links["_peak"] != "") & (links["_gene"] != "")]
    links = links.drop_duplicates(["_peak", "_gene"])

    summary = {
        "link_file": member,
        "columns": columns,
        "peak_column": peak_col,
        "gene_column": gene_col,
        "correlation_column": corr_col,
        "correlation_filter_applied": correlation_filter,
        "min_abs_correlation": min_abs_correlation if corr_col is not None else None,
        "links_original": original_n,
        "links_after_filter": len(links),
        "unique_peaks": int(links["_peak"].nunique()),
        "unique_genes": int(links["_gene"].nunique()),
    }
    return links[["_peak", "_gene", "_weight"]], summary


def write_gzip_tsv(df: pd.DataFrame, path: Path, index_label: str = "gene") -> None:
    df.to_csv(path, sep="\t", compression="gzip", index=True, index_label=index_label)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scatac-zip", type=Path, default=SCATAC_DEFAULT)
    ap.add_argument("--integration-zip", type=Path, default=INTEGRATION_DEFAULT)
    ap.add_argument("--output", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--min-abs-correlation", type=float, default=0.45)
    ap.add_argument("--max-gene-activity-cpm", type=float, default=1e12)
    args = ap.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    summary: dict = {
        "candidate": "GSE242421",
        "assembly": "hg38",
        "source": "https://zenodo.org/records/8313962",
        "scATAC_zip": str(args.scatac_zip),
        "integration_zip": str(args.integration_zip),
        "scATAC_sha256": sha256(args.scatac_zip),
        "integration_sha256": sha256(args.integration_zip),
        "expected_samples": EXPECTED_SAMPLES,
        "min_abs_correlation": args.min_abs_correlation,
        "derivation": "sample-level peak aggregation followed by binary peak-gene incidence multiplication; no frozen module fitting",
    }

    with zipfile.ZipFile(args.scatac_zip) as scatac_zf, zipfile.ZipFile(args.integration_zip) as integration_zf:
        cells_member = member_by_basename(scatac_zf, "cells.tsv")
        peaks_member = member_by_basename(scatac_zf, "peaks.bed")
        matrix_member = member_by_basename(scatac_zf, "cell_x_peak.mtx.gz")
        links_member = member_by_basename(integration_zf, "peak_gene_links_fdr1e-4.tsv")

        cells = read_cells(scatac_zf, cells_member)
        peaks = read_peaks(scatac_zf, peaks_member)
        matrix = load_mtx_from_zip(scatac_zf, matrix_member)
        links, link_summary = read_peak_gene_links(
            integration_zf, links_member, args.min_abs_correlation
        )

    summary["members"] = {
        "cells": cells_member,
        "peaks": peaks_member,
        "matrix": matrix_member,
        "links": links_member,
    }
    summary["matrix_shape"] = list(matrix.shape)
    summary["n_peak_ids"] = len(peaks)
    summary["n_cells"] = len(cells)
    summary["link_summary"] = link_summary

    if matrix.shape[0] != len(peaks):
        raise ValueError(f"Matrix rows {matrix.shape[0]} != peaks {len(peaks)}")
    if matrix.shape[1] != len(cells):
        raise ValueError(f"Matrix columns {matrix.shape[1]} != cells {len(cells)}")

    cells["sample"] = cells["sample"].astype(str)
    unknown_samples = sorted(set(cells["sample"]) - set(EXPECTED_SAMPLES))
    missing_samples = sorted(set(EXPECTED_SAMPLES) - set(cells["sample"]))
    if unknown_samples:
        raise ValueError(f"Unexpected sample labels: {unknown_samples}")
    if missing_samples:
        raise ValueError(f"Missing expected sample labels: {missing_samples}")

    sample_order = EXPECTED_SAMPLES
    sample_to_cols: dict[str, np.ndarray] = {}
    sample_counts: dict[str, int] = {}
    sample_values = cells["sample"].to_numpy()
    for sample in sample_order:
        idx = np.flatnonzero(sample_values == sample)
        sample_to_cols[sample] = idx
        sample_counts[sample] = int(len(idx))

    # Aggregate the peak-by-cell matrix to peak-by-sample mean accessibility.
    peak_sample = np.zeros((len(peaks), len(sample_order)), dtype=np.float64)
    for j, sample in enumerate(sample_order):
        idx = sample_to_cols[sample]
        peak_sample[:, j] = np.asarray(matrix[:, idx].mean(axis=1)).ravel()

    # Normalize peak accessibility to counts per million at the sample level.
    peak_sample_cpm = peak_sample / np.maximum(peak_sample.sum(axis=0, keepdims=True), 1e-12) * 1e6
    peak_df = pd.DataFrame(peak_sample_cpm, index=peaks, columns=sample_order)
    write_gzip_tsv(peak_df, args.output / "02_peak_sample_cpm.tsv.gz", index_label="peak")

    peak_index = {p: i for i, p in enumerate(peaks)}
    kept = links[links["_peak"].isin(peak_index)].copy()
    dropped_peak_links = int(len(links) - len(kept))
    kept["peak_idx"] = kept["_peak"].map(peak_index).astype(int)

    genes = sorted(kept["_gene"].unique())
    gene_index = {g: i for i, g in enumerate(genes)}
    kept["gene_idx"] = kept["_gene"].map(gene_index).astype(int)

    incidence = sparse.coo_matrix(
        (
            np.ones(len(kept), dtype=np.float64),
            (kept["peak_idx"].to_numpy(), kept["gene_idx"].to_numpy()),
        ),
        shape=(len(peaks), len(genes)),
    ).tocsr()

    gene_activity = peak_sample_cpm.T @ incidence
    gene_activity = np.asarray(gene_activity)
    gene_activity_log = np.log1p(gene_activity)
    gene_activity_df = pd.DataFrame(
        gene_activity_log.T,
        index=genes,
        columns=sample_order,
    )
    write_gzip_tsv(gene_activity_df, args.output / "01_gene_activity_log1p_cpm.tsv.gz")

    coverage = gene_activity_df.notna().mean(axis=1)
    summary["sample_cell_counts"] = sample_counts
    summary["peak_gene_links_retained_in_peak_set"] = int(len(kept))
    summary["peak_gene_links_dropped_unmatched_peak"] = dropped_peak_links
    summary["genes_with_activity"] = int(len(genes))
    summary["gene_activity_shape"] = list(gene_activity_df.shape)
    summary["gene_activity_nonzero_fraction"] = float((gene_activity_df.to_numpy() > 0).mean())
    summary["sample_library_peak_cpm"] = {
        sample: float(peak_sample_cpm[:, j].sum())
        for j, sample in enumerate(sample_order)
    }
    summary["decision"] = "GENE_ACTIVITY_DERIVED"

    (args.output / "03_peak_gene_mapping_summary.json").write_text(
        json.dumps(
            {
                **link_summary,
                "peak_ids_in_scATAC": len(peaks),
                "links_retained_in_peak_set": int(len(kept)),
                "links_dropped_unmatched_peak": dropped_peak_links,
                "genes_with_activity": len(genes),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    (args.output / "04_derivation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print("GSE242421 Z4 GENE-ACTIVITY DERIVATION")
    print(f"cells: {len(cells)}")
    print(f"peaks: {len(peaks)}")
    print(f"peak-gene links after filter: {len(links)}")
    print(f"peak-gene links retained in peak set: {len(kept)}")
    print(f"genes with activity: {len(genes)}")
    print("samples: " + ",".join(sample_order))
    print("module fitting/refinement: NONE")
    print("decision: GENE_ACTIVITY_DERIVED")
    print(f"output: {args.output}")


if __name__ == "__main__":
    main()
