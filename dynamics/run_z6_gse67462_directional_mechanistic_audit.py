"""Run a directional temporal-order audit for GSE67462/GSE67520.

Diagnostic only: frozen Z6 predictive support and all prior validation decisions
remain unchanged.
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
    _build_modality_matrix,
    _discover_files,
    _expression_common,
    _load_tss,
    _spearman,
)

DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_SOFT = "Data/GPL19972_family.soft.gz"
DEFAULT_ASSIGNMENTS = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_OUT = "results/Dynamics/z6_gse67462_directional_mechanistic_audit"
ACTIVE = ["h3k27ac", "h3k4me3", "rnapii", "total_oct4"]
NEGATIVE = "h3k27me3"
ALL_TESTED = ACTIVE + [NEGATIVE]
MODULES = [1, 4, 6]
SEED = 67462


def _normalise_expression(expr: pd.DataFrame, validated: set[str]) -> pd.DataFrame:
    keep = {}
    for col in expr.columns:
        key = _norm_symbol(col)
        if key in validated and key not in keep:
            keep[key] = col
    out = expr.loc[:, list(keep.values())].copy()
    out.columns = list(keep.keys())
    return out


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
        entries = sorted([(t, d) for m, t, d in values if m == modality and t in expression_times], key=lambda x: x[0])
        if not entries:
            continue
        times = [t for t, _ in entries]
        matrices[modality] = pd.DataFrame(
            [{g: d.get(g, 0.0) for g in genes} for _, d in entries],
            index=times,
        )
    return matrices


def _gene_effects(expr: pd.DataFrame, reg: pd.DataFrame, genes: list[str], module: int, modality: str) -> pd.DataFrame:
    times = sorted([t for t in reg.index if t in expr.index])
    if len(times) < 3:
        return pd.DataFrame()
    rows = []
    expected_sign = -1.0 if modality == NEGATIVE else 1.0
    for gene in genes:
        if gene not in expr.columns or gene not in reg.columns:
            continue
        r = reg.loc[times, gene].to_numpy(dtype=float)
        e = expr.loc[times, gene].to_numpy(dtype=float)
        concurrent = _spearman(r, e)
        lead = _spearman(r[:-1], e[1:])
        gain = lead - concurrent if np.isfinite(lead) and np.isfinite(concurrent) else np.nan
        dr = np.diff(r)
        de_next = np.diff(e)
        valid = np.isfinite(dr) & np.isfinite(de_next)
        if valid.any():
            signed = np.sign(dr[valid] * de_next[valid]) * expected_sign
            directional = float(np.mean(signed > 0))
        else:
            directional = np.nan
        rows.append({
            "module": module,
            "modality": modality,
            "gene": gene,
            "n_timepoints": len(times),
            "concurrent_rho": concurrent,
            "lead_rho": lead,
            "lead_gain": gain,
            "directional_agreement": directional,
        })
    return pd.DataFrame(rows)


def _circular_null(expr: pd.DataFrame, reg: pd.DataFrame, genes: list[str], observed: float, permutations: int, rng: np.random.Generator) -> tuple[float, float, float]:
    times = sorted([t for t in reg.index if t in expr.index])
    if len(times) < 3:
        return np.nan, np.nan, np.nan
    expr2 = expr.loc[times, genes]
    reg2 = reg.loc[times, genes]
    gains = []
    n = len(times)
    for _ in range(permutations):
        shift = int(rng.integers(1, n))
        shifted = reg2.iloc[np.roll(np.arange(n), shift)]
        gene_gain = []
        for gene in genes:
            a = shifted[gene].to_numpy(dtype=float)
            b = expr2[gene].to_numpy(dtype=float)
            c = _spearman(a, b)
            l = _spearman(a[:-1], b[1:])
            if np.isfinite(c) and np.isfinite(l):
                gene_gain.append(l - c)
        if gene_gain:
            gains.append(float(np.median(gene_gain)))
    if not gains or not np.isfinite(observed):
        return np.nan, np.nan, np.nan
    arr = np.asarray(gains)
    q95 = float(np.quantile(arr, 0.95))
    p = float((1 + np.sum(arr >= observed)) / (1 + len(arr)))
    return q95, p, float(np.median(arr))


def _summarise(permutations: int, rng: np.random.Generator, expr: pd.DataFrame, matrices: dict, assignments: pd.DataFrame):
    gene_rows = []
    mod_rows = []
    perm_rows = []
    for module in MODULES:
        module_genes = sorted(set(assignments.loc[assignments.module.astype(int) == module, "gene"].astype(str)) & set(expr.columns))
        for modality in ALL_TESTED:
            if modality not in matrices:
                continue
            sub = _gene_effects(expr, matrices[modality], module_genes, module, modality)
            if sub.empty:
                continue
            gene_rows.append(sub)
            observed = float(sub.lead_gain.median()) if sub.lead_gain.notna().any() else np.nan
            q95, p, null_median = _circular_null(expr, matrices[modality], module_genes, observed, permutations, rng)
            mod_rows.append({
                "module": module,
                "modality": modality,
                "n_genes": int(len(sub)),
                "median_concurrent_rho": float(sub.concurrent_rho.median()),
                "median_lead_rho": float(sub.lead_rho.median()),
                "median_lead_gain": observed,
                "median_directional_agreement": float(sub.directional_agreement.median()),
                "lead_gain_null_q95": q95,
                "lead_gain_null_median": null_median,
                "permutation_p": p,
                "lead_gain_exceeds_null": bool(np.isfinite(observed) and np.isfinite(q95) and observed > q95),
            })
            perm_rows.append({"module": module, "modality": modality, "observed_median_lead_gain": observed, "null_q95": q95, "null_median": null_median, "permutation_p": p, "permutations": permutations})
    gene_effects = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    modality_summary = pd.DataFrame(mod_rows)
    module_rows = []
    for module in MODULES:
        s = modality_summary[modality_summary.module == module]
        active = s[s.modality.isin(ACTIVE)]
        module_rows.append({
            "module": module,
            "n_active_modalities_tested": int(len(active)),
            "n_active_modalities_positive_lead_gain": int((active.median_lead_gain > 0).sum()),
            "n_active_modalities_lead_gain_above_null": int(active.lead_gain_exceeds_null.sum()),
            "n_active_modalities_permutation_p_lt_0_05": int((active.permutation_p < 0.05).sum()),
            "median_active_lead_gain": float(active.median_lead_gain.median()) if not active.empty else np.nan,
            "median_active_directional_agreement": float(active.median_directional_agreement.median()) if not active.empty else np.nan,
            "h3k27me3_lead_gain": float(s.loc[s.modality == NEGATIVE, "median_lead_gain"].iloc[0]) if not s.loc[s.modality == NEGATIVE].empty else np.nan,
        })
    return gene_effects, modality_summary, pd.DataFrame(module_rows), pd.DataFrame(perm_rows)


def _report(out: Path, modality: pd.DataFrame, modules: pd.DataFrame) -> None:
    lines = [
        "# GSE67462 Z6 — directional mechanistic audit",
        "",
        "This diagnostic tests temporal ordering between regulatory trajectories and subsequent expression. Positive lead gain means the regulatory signal at time t is more concordant with expression at t+1 than with expression at t.",
        "",
        "## Modality summary",
        "",
        "| Module | Modality | Genes | Concurrent rho | Lead rho | Lead gain | Directional agreement | Null q95 | p |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in modality.sort_values(["module", "modality"]).iterrows():
        lines.append(f"| M{int(r.module)} | {r.modality} | {int(r.n_genes)} | {r.median_concurrent_rho:.3f} | {r.median_lead_rho:.3f} | {r.median_lead_gain:.3f} | {r.median_directional_agreement:.3f} | {r.lead_gain_null_q95:.3f} | {r.permutation_p:.3f} |")
    lines += ["", "## Module-level synthesis", "", "| Module | Active modalities tested | Positive lead gain | Above null | p<0.05 | Median active gain | Median directional agreement | H3K27me3 gain |", "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for _, r in modules.iterrows():
        lines.append(f"| M{int(r.module)} | {int(r.n_active_modalities_tested)} | {int(r.n_active_modalities_positive_lead_gain)} | {int(r.n_active_modalities_lead_gain_above_null)} | {int(r.n_active_modalities_permutation_p_lt_0_05)} | {r.median_active_lead_gain:.3f} | {r.median_active_directional_agreement:.3f} | {r.h3k27me3_lead_gain:.3f} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "A positive lead gain supports a temporal-order hypothesis in which the regulatory trajectory is aligned with later expression change. This remains observational evidence and does not establish causality.",
        "",
        "H3K27me3 is retained as a negative-control/reference modality and is not included in the active multimodal synthesis.",
        "",
        "Because GSE67462 contains sparse bulk timepoints, absence of a lead effect cannot be interpreted as evidence that a regulatory mechanism is absent.",
        "",
        "## Frozen boundary",
        "",
        "No Z6 predictive support threshold, validation split, temporal module assignment or multimodal-support decision was changed.",
    ]
    (out / "GSE67462_Z6_directional_mechanistic_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
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
    _, expr_raw, _ = _expression_common(matrix, metadata)
    expr = _normalise_expression(expr_raw, validated)
    files = _discover_files(Path(args.modality_root))
    missing = [m for m in ALL_TESTED if not files[m]]
    if missing:
        raise RuntimeError(f"Missing modality files: {missing}")
    values_raw, _ = _build_modality_matrix(files, tss)
    values = _normalise_values(values_raw, validated)
    matrices = _matrix_by_modality(values, list(expr.index), list(expr.columns))
    assignments = pd.read_csv(args.assignments)
    rng = np.random.default_rng(SEED)
    gene_effects, modality_summary, module_summary, perm_summary = _summarise(permutations=args.permutations, rng=rng, expr=expr, matrices=matrices, assignments=assignments)
    gene_effects.to_csv(out / "01_gene_directional_effects.csv", index=False)
    modality_summary.to_csv(out / "02_modality_summary.csv", index=False)
    module_summary.to_csv(out / "03_module_directional_summary.csv", index=False)
    perm_summary.to_csv(out / "04_permutation_summary.csv", index=False)
    manifest = {
        "modules": MODULES,
        "active_modalities": ACTIVE,
        "negative_control": NEGATIVE,
        "validated_expression_genes": int(expr.shape[1]),
        "permutations": args.permutations,
        "null": "gene-wise circular time shifts excluding shift 0",
        "diagnostic_only": True,
        "frozen_z6_support_unchanged": True,
        "mapping_provenance": mapping_summary,
    }
    (out / "05_manifest.json").write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    _report(out, modality_summary, module_summary)
    print("GSE67462 Z6 DIRECTIONAL MECHANISTIC AUDIT")
    print(f"validated genes: {expr.shape[1]}")
    print(f"modules tested: {MODULES}")
    print(f"outputs: {out}")
    print("frozen Z6 support unchanged: true")


if __name__ == "__main__":
    main()
