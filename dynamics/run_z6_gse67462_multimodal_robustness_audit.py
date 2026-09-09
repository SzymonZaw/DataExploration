"""Audit GSE67462/GSE67520 multimodal support after identifier correction.

Diagnostic only. Frozen Z6 predictive thresholds/support are not modified.
The audit reports supported modalities from the mapped analysis and tests
sensitivity to deterministic peak-to-TSS assignment radius using the validated
GPL19972 -> symbol -> TSS expression universe.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamics import validation
from dynamics.run_z6_gse67462_identifier_mapping_audit import (
    _build_mapping_report,
    _norm_symbol,
    _read_soft_platform,
)
from dynamics.run_z6_gse67462_multimodal_validation_analysis import (
    MODALITIES,
    SEED,
    _discover_files,
    _expression_common,
    _load_tss,
    _parse_peak_file,
    _permutation_test,
    _spearman,
)

DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_SOFT = "Data/GPL19972_family.soft.gz"
DEFAULT_ROOT = "Data/GSE67520"
DEFAULT_MAPPED = "results/Dynamics/z6_gse67462_multimodal_validation_mapped/02_mapped_multimodal_concordance.csv"
DEFAULT_OUTPUT = "results/Dynamics/z6_gse67462_multimodal_robustness_audit"
RADII = {"promoter_2kb": 2000, "nearest_tss_25kb": 25000, "nearest_tss_50kb": 50000, "nearest_tss_100kb": 100000}


def _validated_expression(gtf: Path, platform_soft: Path):
    matrix, metadata = validation._load_common_space()
    tss = _load_tss(gtf)
    platform = _read_soft_platform(platform_soft)
    mapping, summary = _build_mapping_report(pd.Index(matrix.index.astype(str)), tss, platform)
    valid = set(mapping.loc[mapping["platform_symbol_tss_match"].fillna(0).astype(int) > 0, "expression_id"].map(_norm_symbol))
    _, expr_raw, branches = _expression_common(matrix, metadata)
    keep = {}
    for col in expr_raw.columns:
        key = _norm_symbol(col)
        if key in valid and key not in keep:
            keep[key] = col
    expr = expr_raw.loc[:, list(keep.values())].copy()
    expr.columns = list(keep.keys())
    return expr, tss, branches, summary


def _nearest_values(files, tss, radius, valid_symbols):
    by_chrom = {c: g.sort_values("tss") for c, g in tss.groupby("chrom")}
    rows = []
    provenance = []
    for modality, entries in files.items():
        for time, path in entries:
            values = {}
            peaks = _parse_peak_file(path)
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
                        dist = abs(float(row.tss) - midpoint)
                        if dist <= radius:
                            candidates.append((dist, row))
                if not candidates:
                    continue
                _, row = min(candidates, key=lambda x: x[0])
                gene = _norm_symbol(row.gene)
                if gene not in valid_symbols:
                    continue
                values[gene] = values.get(gene, 0.0) + max(float(score), 0.0)
                assigned += 1
            rows.append((modality, time, values))
            provenance.append({"modality": modality, "time": str(time), "n_peaks": len(peaks), "n_assigned_valid": assigned})
    return rows, provenance


def _concordance(expr, values):
    rows = []
    for modality in sorted(set(m for m, _, _ in values)):
        entries = sorted([(t, v) for m, t, v in values if m == modality], key=lambda x: x[0])
        valid_entries = [(t, v) for t, v in entries if t in expr.index]
        times = [t for t, _ in valid_entries]
        if len(times) < 5:
            rows.append({"modality": modality, "n_timepoints": len(times), "n_genes": 0, "observed_global_spearman": np.nan})
            continue
        genes = pd.Index(expr.columns.astype(str))
        mod = pd.DataFrame([{g: v.get(g, 0.0) for g in genes} for _, v in valid_entries], index=times)
        common = genes[genes.isin(mod.columns)]
        rhos = [_spearman(expr.loc[times, g].to_numpy(), mod.loc[times, g].to_numpy()) for g in common]
        rhos = np.asarray(rhos, dtype=float)
        rhos = rhos[np.isfinite(rhos)]
        rows.append({
            "modality": modality,
            "n_timepoints": len(times),
            "n_genes": int(len(rhos)),
            "median_gene_spearman": float(np.median(rhos)) if len(rhos) else np.nan,
            "mean_gene_spearman": float(np.mean(rhos)) if len(rhos) else np.nan,
            "observed_global_spearman": _spearman(expr.loc[times, common].to_numpy().ravel(), mod.loc[times, common].to_numpy().ravel()),
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--platform-soft", default=DEFAULT_SOFT)
    parser.add_argument("--modality-root", default=DEFAULT_ROOT)
    parser.add_argument("--mapped-results", default=DEFAULT_MAPPED)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--permutations", type=int, default=1000)
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    mapped = pd.read_csv(args.mapped_results)
    required = {"modality", "n_timepoints", "n_genes", "observed_global_spearman", "null_q95", "permutation_p", "status"}
    missing = required - set(mapped.columns)
    if missing:
        raise RuntimeError(f"Mapped concordance is missing columns: {sorted(missing)}")
    mapped["supported"] = (mapped["status"] == "OK") & (mapped["permutation_p"] < 0.05) & (mapped["observed_global_spearman"] > mapped["null_q95"])
    supported = mapped[mapped["supported"]].copy()
    supported.to_csv(out / "01_supported_modalities.csv", index=False)

    expr, tss, branches, mapping_summary = _validated_expression(Path(args.gtf), Path(args.platform_soft))
    valid_symbols = set(expr.columns)
    files = _discover_files(Path(args.modality_root))
    if any(not files[k] for k in MODALITIES):
        raise RuntimeError(f"Missing modality files: {[k for k in MODALITIES if not files[k]]}")

    robustness = []
    provenance = []
    for mapping, radius in RADII.items():
        values, prov = _nearest_values(files, tss, radius, valid_symbols)
        c = _concordance(expr, values)
        p = _permutation_test(expr, values, args.permutations, SEED)
        r = c.merge(p, on="modality", how="left")
        r.insert(0, "mapping", mapping)
        r["radius_bp"] = radius
        r["supported"] = (r["n_genes"] >= 50) & (r["permutation_p"] < 0.05) & (r["observed_global_spearman"] > r["null_q95"])
        robustness.append(r)
        for x in prov:
            x["mapping"] = mapping
            provenance.append(x)
    robustness = pd.concat(robustness, ignore_index=True)
    pd.DataFrame(provenance).to_csv(out / "02_peak_assignment_provenance.csv", index=False)
    robustness.to_csv(out / "03_assignment_robustness.csv", index=False)

    stable = robustness.groupby("modality")["supported"].sum()
    summary = {
        "dataset_expression": "GSE67462",
        "dataset_regulatory": "GSE67520",
        "validated_expression_genes": int(expr.shape[1]),
        "validated_fraction": float(expr.shape[1] / 11899),
        "branches": branches,
        "mapped_supported_modalities": sorted(supported["modality"].tolist()),
        "n_mapped_supported_modalities": int(len(supported)),
        "support_count_across_assignment_mappings": {str(k): int(v) for k, v in stable.items()},
        "frozen_support_criterion": "permutation_p < 0.05 AND observed_global_spearman > null_q95; minimum 50 genes for mapping robustness audit",
        "interpretation": "ROBUST_MULTIMODAL_SUPPORT" if len(supported) >= 2 and (stable >= 2).sum() >= 2 else "ASSIGNMENT_SENSITIVE_SUPPORT",
        "identifier_provenance": "GSE67462 GPL19972 RefSeq feature -> gene symbol -> mm9.refGene TSS",
        "mapping_summary": mapping_summary,
    }
    (out / "04_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    print("Mapped supported modalities:")
    print(supported[["modality", "n_timepoints", "n_genes", "observed_global_spearman", "null_q95", "permutation_p"]].to_string(index=False))
    print("\nAssignment robustness:")
    print(robustness[["mapping", "modality", "n_genes", "observed_global_spearman", "null_q95", "permutation_p", "supported"]].to_string(index=False))
    print(f"\nInterpretation: {summary['interpretation']}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
