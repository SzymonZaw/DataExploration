"""Run provenance-aware GSE67462/GSE67520 multimodal concordance.

Diagnostic only. Frozen Z6 predictive support is unchanged.

This wrapper fixes the identifier layer used by the earlier multimodal analysis:
GSE67462 common-space expression features are gene symbols, while GPL19972
stores RefSeq-like feature IDs (NM_..._at) and symbols in Description. The
GPL19972 SOFT mapping is therefore used to define the validated expression
symbol universe before comparing it with GSE67520 regulatory trajectories.

The regulatory peak-to-gene assignment remains the existing deterministic
TSS-centered assignment. Gene symbols are normalized case-insensitively before
joining expression and regulatory matrices. This is an identifier/provenance
correction only; no Z6 thresholds or support criteria are changed.
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
    _read_soft_platform,
    _norm_symbol,
)
from dynamics.run_z6_gse67462_multimodal_validation_analysis import (
    DAY_TO_HOURS,
    MODALITIES,
    SEED,
    _build_modality_matrix,
    _discover_files,
    _expression_common,
    _load_tss,
    _permutation_test,
    _spearman,
)

DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_SOFT = "Data/GPL19972_family.soft.gz"


def _normalized_expression(expr: pd.DataFrame, validated_symbols: set[str]) -> pd.DataFrame:
    """Normalize expression columns to upper-case symbols and keep validated ones."""
    keep = {}
    for col in expr.columns:
        key = _norm_symbol(col)
        if key in validated_symbols and key not in keep:
            keep[key] = col
    out = expr.loc[:, list(keep.values())].copy()
    out.columns = list(keep.keys())
    return out


def _normalized_modality_values(values, validated_symbols: set[str]):
    normalized = []
    for modality, time, mapping in values:
        out = {}
        for gene, value in mapping.items():
            key = _norm_symbol(gene)
            if key in validated_symbols:
                out[key] = out.get(key, 0.0) + float(value)
        normalized.append((modality, time, out))
    return normalized


def _trajectory_concordance(expr, modality_values):
    rows = []
    genes = pd.Index(expr.columns.astype(str))
    for modality in sorted(set(m for m, _, _ in modality_values)):
        entries = [(t, v) for m, t, v in modality_values if m == modality and isinstance(t, (int, float, np.integer, np.floating))]
        entries = sorted(entries, key=lambda x: x[0])
        times = [t for t, _ in entries if t in expr.index]
        if len(times) < 5:
            rows.append({"modality": modality, "n_timepoints": len(times), "n_genes": 0,
                         "median_gene_spearman": np.nan, "mean_gene_spearman": np.nan,
                         "observed_global_spearman": np.nan, "status": "INSUFFICIENT_TIMEPOINTS"})
            continue
        modality_df = pd.DataFrame(
            [{g: v.get(g, 0.0) for g in genes} for t, v in entries if t in expr.index],
            index=[t for t, _ in entries if t in expr.index],
        )
        common_genes = genes[genes.isin(modality_df.columns)]
        gene_rhos = np.asarray([
            _spearman(expr.loc[times, g].to_numpy(), modality_df.loc[times, g].to_numpy())
            for g in common_genes
        ])
        gene_rhos = gene_rhos[np.isfinite(gene_rhos)]
        expr_flat = expr.loc[times, common_genes].to_numpy().ravel()
        mod_flat = modality_df.loc[times, common_genes].to_numpy().ravel()
        rows.append({
            "modality": modality,
            "n_timepoints": len(times),
            "n_genes": int(len(gene_rhos)),
            "median_gene_spearman": float(np.nanmedian(gene_rhos)) if len(gene_rhos) else np.nan,
            "mean_gene_spearman": float(np.nanmean(gene_rhos)) if len(gene_rhos) else np.nan,
            "observed_global_spearman": _spearman(expr_flat, mod_flat),
            "status": "OK",
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--modality-root", default="Data/GSE67520")
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--platform-soft", default=DEFAULT_SOFT)
    parser.add_argument("--output", default="results/Dynamics/z6_gse67462_multimodal_validation_mapped")
    parser.add_argument("--permutations", type=int, default=1000)
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    matrix, metadata = validation._load_common_space()
    tss = _load_tss(Path(args.gtf))
    platform = _read_soft_platform(Path(args.platform_soft))
    mapping, mapping_summary = _build_mapping_report(pd.Index(matrix.index.astype(str)), tss, platform)
    validated_symbols = set(mapping.loc[mapping["platform_symbol_tss_match"].fillna(0).astype(int) > 0, "expression_id"].map(_norm_symbol))
    if not validated_symbols:
        raise RuntimeError("No validated expression symbols were mapped through GPL19972 to TSS.")

    expr_times, expr_raw, branches = _expression_common(matrix, metadata)
    expr = _normalized_expression(expr_raw, validated_symbols)
    files = _discover_files(Path(args.modality_root))
    missing = [k for k in MODALITIES if not files[k]]
    if missing:
        raise RuntimeError(f"Missing modality files: {missing}")
    values_raw, provenance = _build_modality_matrix(files, tss)
    values = _normalized_modality_values(values_raw, validated_symbols)

    concordance = _trajectory_concordance(expr, values)
    perm = _permutation_test(expr, values, args.permutations, SEED)
    result = concordance.merge(perm, on="modality", how="left")
    supported = result[
        (result.status == "OK")
        & (result.permutation_p < 0.05)
        & (result.observed_global_spearman > result.null_q95)
    ]
    if len(supported) >= 2:
        interpretation = "MULTIMODAL_DYNAMIC_SUPPORT"
    elif len(supported) == 1 or np.any((result.permutation_p < 0.05).fillna(False)):
        interpretation = "PARTIAL_MULTIMODAL_SUPPORT"
    else:
        interpretation = "EXPRESSION_ONLY"

    provenance.to_csv(out / "01_peak_assignment_provenance.csv", index=False)
    result.to_csv(out / "02_mapped_multimodal_concordance.csv", index=False)
    (out / "03_mapping_summary.json").write_text(
        json.dumps(mapping_summary, indent=2, default=str), encoding="utf-8"
    )
    summary = {
        "dataset_expression": "GSE67462",
        "dataset_multimodal": "GSE67520",
        "branches": branches,
        "expression_timepoints": expr_times.tolist(),
        "n_expression_genes_raw": int(expr_raw.shape[1]),
        "n_expression_genes_validated_via_GPL19972": int(expr.shape[1]),
        "n_annotated_tss": int(len(tss)),
        "n_modalities": len(MODALITIES),
        "n_supported_modalities": int(len(supported)),
        "interpretation": interpretation,
        "identifier_provenance": "GPL19972 SOFT RefSeq feature -> gene symbol -> mm9.refGene TSS; expression/regulatory gene joins normalized case-insensitively",
        "mapping_summary": mapping_summary,
        "frozen_rules": {
            "modify_z6_support": False,
            "model_tuning": False,
            "time_label_permutation": True,
        },
    }
    (out / "04_mapped_multimodal_summary.json").write_text(
        json.dumps(summary, indent=2, default=str), encoding="utf-8"
    )

    print("Provenance-aware GSE67462/GSE67520 multimodal validation complete.")
    print(f"Raw common-space genes: {expr_raw.shape[1]}")
    print(f"Validated genes via GPL19972 -> TSS: {expr.shape[1]}")
    print(f"Validated fraction: {expr.shape[1] / expr_raw.shape[1]:.6f}")
    print(f"Supported modalities: {len(supported)}")
    print(f"Interpretation: {interpretation}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
