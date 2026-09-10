"""Sensitivity audit for promoter versus distal regulatory peak-to-gene mapping.

Diagnostic only. It does not modify frozen Z6 predictive support.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamics import validation
from dynamics.run_z6_gse67462_multimodal_validation_analysis import (
    MODALITIES,
    _discover_files,
    _expression_common,
    _load_tss,
    _parse_peak_file,
    _permutation_test,
    _trajectory_concordance,
)

DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_ROOT = "Data/GSE67520"
DEFAULT_OUTPUT = "results/Dynamics/z6_gse67462_regulatory_element_audit"
MAPPING_RADII = {
    "promoter_2kb": 2000,
    "nearest_tss_25kb": 25000,
    "nearest_tss_50kb": 50000,
    "nearest_tss_100kb": 100000,
}
N_PERMUTATIONS = 1000
SEED = 412


def _nearest_tss_burden(path: Path, tss: pd.DataFrame, radius: int):
    peaks = _parse_peak_file(path)
    by_chrom = {c: g.sort_values("tss") for c, g in tss.groupby("chrom")}
    values = {}
    assigned = 0
    for chrom, start, end, score in peaks:
        if chrom not in by_chrom:
            continue
        midpoint = (start + end) // 2
        g = by_chrom[chrom]
        coords = g["tss"].to_numpy(dtype=float)
        pos = int(np.searchsorted(coords, midpoint))
        candidates = []
        for j in (pos - 1, pos):
            if 0 <= j < len(g):
                row = g.iloc[j]
                distance = abs(float(row.tss) - midpoint)
                if distance <= radius:
                    candidates.append((distance, row))
        if not candidates:
            continue
        _, row = min(candidates, key=lambda x: x[0])
        gene = str(row.gene)
        values[gene] = values.get(gene, 0.0) + max(float(score), 0.0)
        assigned += 1
    return values, len(peaks), assigned


def _build_values(files, tss, radius):
    rows = []
    provenance = []
    for modality, entries in files.items():
        for time, path in entries:
            values, n_peaks, n_assigned = _nearest_tss_burden(path, tss, radius)
            rows.append((modality, time, values))
            provenance.append({"modality": modality, "time": str(time), "file": str(path), "n_peaks": n_peaks, "n_assigned": n_assigned})
    return rows, pd.DataFrame(provenance)


def _classify(result: pd.DataFrame):
    adequate = result["n_genes"] >= 50
    supported = adequate & (result["permutation_p"] < 0.05) & (result["observed_global_spearman"] > result["null_q95"])
    result = result.copy()
    result["adequate_gene_count"] = adequate
    result["supported"] = supported
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--modality-root", default=DEFAULT_ROOT)
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    tss = _load_tss(Path(args.gtf))
    matrix, metadata = validation._load_common_space()
    expr_times, expr, branches = _expression_common(matrix, metadata)
    files = _discover_files(Path(args.modality_root))
    if any(not files[k] for k in MODALITIES):
        missing = [k for k in MODALITIES if not files[k]]
        raise RuntimeError(f"Missing modality files: {missing}")

    all_results = []
    all_provenance = []
    for mapping, radius in MAPPING_RADII.items():
        values, provenance = _build_values(files, tss, radius)
        concordance = _trajectory_concordance(expr, values)
        perm = _permutation_test(expr, values, args.permutations, SEED)
        result = concordance.merge(perm, on="modality", how="left")
        result.insert(0, "mapping", mapping)
        result["radius_bp"] = radius
        all_results.append(_classify(result))
        provenance.insert(0, "mapping", mapping)
        all_provenance.append(provenance)

    result = pd.concat(all_results, ignore_index=True)
    provenance = pd.concat(all_provenance, ignore_index=True)
    result.to_csv(out / "01_mapping_concordance.csv", index=False)
    provenance.to_csv(out / "02_peak_mapping_provenance.csv", index=False)

    supported = result[result["supported"]]
    distal_supported = supported[supported["radius_bp"] > 2000]
    support_by_modality = supported.groupby("modality")["mapping"].nunique().to_dict() if not supported.empty else {}
    distal_modalities = sorted(distal_supported["modality"].unique().tolist())

    if len(support_by_modality) >= 2 and len(distal_supported) >= 2:
        interpretation = "ROBUST_MULTIMODAL_SUPPORT"
    elif len(distal_modalities) >= 2:
        interpretation = "DISTAL_SENSITIVE_SUPPORT"
    elif len(supported) == 1:
        interpretation = "PARTIAL_REGULATORY_SUPPORT"
    elif not result["adequate_gene_count"].any():
        interpretation = "MAPPING_INADEQUATE"
    else:
        interpretation = "NO_REGULATORY_SUPPORT"

    summary = {
        "expression_dataset": "GSE67462",
        "regulatory_dataset": "GSE67520",
        "branches": branches,
        "expression_timepoints": expr_times.tolist(),
        "n_expression_genes": int(expr.shape[1]),
        "n_annotated_tss": int(len(tss)),
        "mapping_radii_bp": MAPPING_RADII,
        "n_supported_rows": int(len(supported)),
        "modalities_with_distal_support": distal_modalities,
        "support_by_modality": support_by_modality,
        "interpretation": interpretation,
        "permutations": int(args.permutations),
        "seed": SEED,
        "frozen_z6_support_modified": False,
    }
    (out / "03_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("GSE67462 regulatory-element audit complete.")
    print(result[["mapping", "modality", "n_genes", "median_gene_spearman", "observed_global_spearman", "null_q95", "permutation_p", "supported"]].to_string(index=False))
    print(f"Interpretation: {interpretation}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
