"""Test directional regulatory-to-expression effects across GSE67462 replicates.

Diagnostic only. GSE67520 regulatory profiles have one temporal profile per
modality/timepoint, so this tests stability of the expression response across
the two independent GSE67462 expression replicates; it does not claim replicate
stability of the regulatory measurements themselves.
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
    _build_modality_matrix,
    _discover_files,
    _load_tss,
    _spearman,
)
from dynamics.run_z6_within_system_audit import _build_system_data

DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_SOFT = "Data/GPL19972_family.soft.gz"
DEFAULT_ASSIGNMENTS = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_OUT = "results/Dynamics/z6_gse67462_replicate_directional_audit"
ACTIVE = ["h3k27ac", "h3k4me3", "rnapii", "total_oct4"]
NEGATIVE = "h3k27me3"
ALL_TESTED = ACTIVE + [NEGATIVE]
MODULES = [1, 4, 6]
SEED = 67463


def _normalise_values(values, validated: set[str]):
    out = []
    for modality, time, mapping in values:
        if not isinstance(time, (int, float, np.integer, np.floating)):
            continue
        d = {}
        for gene, value in mapping.items():
            key = _norm_symbol(gene)
            if key in validated:
                d[key] = d.get(key, 0.0) + float(value)
        out.append((modality, float(time), d))
    return out


def _matrix_by_modality(values, expression_times, genes):
    matrices = {}
    for modality in ALL_TESTED:
        entries = sorted(
            [(t, d) for m, t, d in values if m == modality and t in expression_times],
            key=lambda x: x[0],
        )
        if entries:
            matrices[modality] = pd.DataFrame(
                [{g: d.get(g, 0.0) for g in genes} for _, d in entries],
                index=[t for t, _ in entries],
            )
    return matrices


def _gene_effects(expr: pd.DataFrame, reg: pd.DataFrame, genes: list[str], module: int, modality: str, replicate: str):
    times = sorted(set(expr.index) & set(reg.index))
    if len(times) < 3:
        return pd.DataFrame()
    expected_sign = -1.0 if modality == NEGATIVE else 1.0
    rows = []
    for gene in genes:
        if gene not in expr.columns or gene not in reg.columns:
            continue
        r = reg.loc[times, gene].to_numpy(float)
        e = expr.loc[times, gene].to_numpy(float)
        concurrent = _spearman(r, e)
        lead = _spearman(r[:-1], e[1:])
        gain = lead - concurrent if np.isfinite(lead) and np.isfinite(concurrent) else np.nan
        dr = np.diff(r)
        de = np.diff(e)
        ok = np.isfinite(dr) & np.isfinite(de)
        directional = float(np.mean(np.sign(dr[ok] * de[ok]) * expected_sign > 0)) if ok.any() else np.nan
        rows.append({
            "module": module,
            "modality": modality,
            "replicate": replicate,
            "gene": gene,
            "n_timepoints": len(times),
            "concurrent_rho": concurrent,
            "lead_rho": lead,
            "lead_gain": gain,
            "directional_agreement": directional,
        })
    return pd.DataFrame(rows)


def _standardise_rank_columns(frame: pd.DataFrame):
    """Return column-wise rank vectors centred and normalised for fast rho."""
    ranked = frame.rank(axis=0, method="average", na_option="keep").to_numpy(float)
    means = np.nanmean(ranked, axis=0)
    centred = ranked - means
    norms = np.sqrt(np.nansum(centred * centred, axis=0))
    good = norms > 0
    scaled = np.full_like(centred, np.nan, dtype=float)
    scaled[:, good] = centred[:, good] / norms[good]
    return scaled, good


def _median_gene_correlations(a: np.ndarray, b: np.ndarray):
    """Column-wise Pearson/Spearman correlations from rank-standardised arrays."""
    valid = np.isfinite(a).all(axis=0) & np.isfinite(b).all(axis=0)
    if not valid.any():
        return np.nan
    corr = np.sum(a[:, valid] * b[:, valid], axis=0)
    return float(np.nanmedian(corr)) if corr.size else np.nan


def _null(expr, reg, genes, observed, permutations, rng):
    """Fast circular-shift null; vectorised across genes and shifts."""
    times = sorted(set(expr.index) & set(reg.index))
    if len(times) < 3 or not np.isfinite(observed):
        return np.nan, np.nan, np.nan

    e = expr.loc[times, genes].astype(float)
    r = reg.loc[times, genes].astype(float)
    e_rank, e_good = _standardise_rank_columns(e)
    r_rank, r_good = _standardise_rank_columns(r)
    good = e_good & r_good
    if not good.any():
        return np.nan, np.nan, np.nan

    e_full = e_rank[:, good]
    r_full = r_rank[:, good]
    null = np.empty(permutations, dtype=float)
    shifts = rng.integers(1, len(times), size=permutations)
    for i, shift in enumerate(shifts):
        order = np.roll(np.arange(len(times)), int(shift))
        rs = r_full[order]
        concurrent = np.sum(rs * e_full, axis=0)
        lead = np.sum(rs[:-1] * e_full[1:], axis=0)
        gains = lead - concurrent
        null[i] = np.nanmedian(gains)

    null = null[np.isfinite(null)]
    if not len(null):
        return np.nan, np.nan, np.nan
    q95 = float(np.quantile(null, .95))
    p = float((1 + np.sum(null >= observed)) / (1 + len(null)))
    return q95, p, float(np.median(null))


def _replicate_assignment(matrix, metadata):
    data, _ = _build_system_data(matrix, metadata, "GSE67462")
    if len(data) != 2:
        raise RuntimeError(f"Expected two GSE67462 expression branches, got {list(data)}")
    branches = []
    for key, (times, values) in sorted(data.items()):
        branches.append((str(key), np.asarray(times, float), values))
    return branches


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", default=DEFAULT_GTF)
    p.add_argument("--platform-soft", default=DEFAULT_SOFT)
    p.add_argument("--assignments", default=DEFAULT_ASSIGNMENTS)
    p.add_argument("--modality-root", default="Data/GSE67520")
    p.add_argument("--output", default=DEFAULT_OUT)
    p.add_argument("--permutations", type=int, default=1000)
    args = p.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    matrix, metadata = validation._load_common_space()
    tss = _load_tss(Path(args.gtf))
    platform = _read_soft_platform(Path(args.platform_soft))
    mapping, mapping_summary = _build_mapping_report(pd.Index(matrix.index.astype(str)), tss, platform)
    validated = set(mapping.loc[mapping["platform_symbol_tss_match"].fillna(0).astype(int) > 0, "expression_id"].map(_norm_symbol))
    if not validated:
        raise RuntimeError("No validated expression symbols mapped through GPL19972 to TSS.")

    branches = _replicate_assignment(matrix, metadata)
    expr_branches = []
    for branch, times, values in branches:
        expr = pd.DataFrame(values, index=times, columns=matrix.index.astype(str))
        keep = {}
        for col in expr.columns:
            key = _norm_symbol(col)
            if key in validated and key not in keep:
                keep[key] = col
        expr = expr.loc[:, list(keep.values())].copy()
        expr.columns = list(keep.keys())
        expr_branches.append((branch, expr))

    files = _discover_files(Path(args.modality_root))
    missing = [m for m in ALL_TESTED if not files[m]]
    if missing:
        raise RuntimeError(f"Missing modality files: {missing}")
    values_raw, _ = _build_modality_matrix(files, tss)
    values = _normalise_values(values_raw, validated)
    all_times = sorted(set().union(*[set(e.index) for _, e in expr_branches]))
    matrices = _matrix_by_modality(values, all_times, sorted(validated))
    assignments = pd.read_csv(args.assignments)
    rng = np.random.default_rng(SEED)

    gene_rows = []
    summary_rows = []
    for module in MODULES:
        module_genes = sorted(set(assignments.loc[assignments.module.astype(int) == module, "gene"].astype(str)) & validated)
        for modality in ALL_TESTED:
            if modality not in matrices:
                continue
            for branch, expr in expr_branches:
                sub = _gene_effects(expr, matrices[modality], module_genes, module, modality, branch)
                if sub.empty:
                    continue
                gene_rows.append(sub)
                observed = float(sub.lead_gain.median()) if sub.lead_gain.notna().any() else np.nan
                q95, pval, null_med = _null(expr, matrices[modality], module_genes, observed, args.permutations, rng)
                summary_rows.append({
                    "module": module,
                    "modality": modality,
                    "replicate": branch,
                    "n_genes": len(sub),
                    "median_concurrent_rho": float(sub.concurrent_rho.median()),
                    "median_lead_rho": float(sub.lead_rho.median()),
                    "median_lead_gain": observed,
                    "median_directional_agreement": float(sub.directional_agreement.median()),
                    "null_q95": q95,
                    "null_median": null_med,
                    "permutation_p": pval,
                    "above_null": bool(np.isfinite(observed) and np.isfinite(q95) and observed > q95),
                })

    gene_df = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    summary = pd.DataFrame(summary_rows)
    replicate_rows = []
    for module in MODULES:
        for modality in ALL_TESTED:
            s = summary[(summary.module == module) & (summary.modality == modality)]
            if len(s) != 2:
                continue
            gains = s.median_lead_gain.to_numpy(float)
            replicate_rows.append({
                "module": module,
                "modality": modality,
                "replicate_1": str(s.iloc[0].replicate),
                "replicate_2": str(s.iloc[1].replicate),
                "gain_rep1": gains[0],
                "gain_rep2": gains[1],
                "same_gain_sign": bool(gains[0] * gains[1] > 0),
                "both_above_null": bool(s.above_null.all()),
                "both_p_lt_0_05": bool((s.permutation_p < .05).all()),
                "replicate_concordant_directional": bool(gains[0] * gains[1] > 0),
            })
    replicate_df = pd.DataFrame(replicate_rows)

    gene_df.to_csv(out / "01_gene_replicate_directional_effects.csv", index=False)
    summary.to_csv(out / "02_replicate_modality_summary.csv", index=False)
    replicate_df.to_csv(out / "03_replicate_concordance.csv", index=False)
    manifest = {
        "modules": MODULES,
        "active_modalities": ACTIVE,
        "negative_control": NEGATIVE,
        "expression_replicates": [b for b, _ in expr_branches],
        "validated_expression_genes": len(validated),
        "permutations": args.permutations,
        "null": "gene-wise circular time shifts excluding shift 0, separately per expression replicate",
        "null_implementation": "vectorised rank-correlation across genes; one random non-zero circular shift per permutation",
        "replicate_scope": "expression response only; GSE67520 regulatory profiles are not independently replicated",
        "diagnostic_only": True,
        "frozen_z6_support_unchanged": True,
        "mapping_provenance": mapping_summary,
    }
    (out / "04_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")

    lines = [
        "# GSE67462 Z6 — replicate-level directional falsification",
        "",
        "This audit asks whether aggregate regulatory-to-expression temporal associations are stable across the two independent GSE67462 expression replicates.",
        "",
        "> Important: GSE67520 supplies one regulatory profile per modality/timepoint. Therefore this is replicate stability of the **expression response**, not replication of the regulatory measurement.",
        "",
        "## Replicate concordance",
        "",
        "| Module | Modality | Gain rep1 | Gain rep2 | Same sign | Both above null | Both p<0.05 |",
        "|---:|---|---:|---:|:---:|:---:|:---:|",
    ]
    for _, r in replicate_df.iterrows():
        lines.append(f"| M{int(r.module)} | {r.modality} | {r.gain_rep1:.3f} | {r.gain_rep2:.3f} | {str(bool(r.same_gain_sign))} | {str(bool(r.both_above_null))} | {str(bool(r.both_p_lt_0_05))} |")
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "Same-sign effects strengthen reproducibility of the expression response, but do not establish causality. Opposite signs argue against a stable replicate-level directional interpretation.",
        "",
        "A complete OCT4 -> chromatin -> expression chain cannot be replicated here because the regulatory profiles do not have independent replicate trajectories in the available GSE67520 series.",
        "",
        "## Frozen boundary",
        "",
        "No Z6 predictive support threshold, validation split, temporal module assignment or multimodal-support decision was changed.",
    ]
    (out / "GSE67462_Z6_replicate_directional_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("GSE67462 Z6 REPLICATE-LEVEL DIRECTIONAL AUDIT")
    print(f"validated genes: {len(validated)}")
    print(f"expression branches: {[b for b, _ in expr_branches]}")
    print(f"modules tested: {MODULES}")
    print(f"outputs: {out}")
    print("frozen Z6 support unchanged: true")


if __name__ == "__main__":
    main()
