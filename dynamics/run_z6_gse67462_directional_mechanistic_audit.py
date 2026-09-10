"""Audit directional temporal ordering in the GSE67462/GSE67520 multimodal signal.

Diagnostic only. It compares contemporaneous versus one-step-ahead
regulatory-to-expression concordance, checks trajectory direction, and
contrasts active modalities with H3K27me3. Frozen Z6 decisions are unchanged.
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
DEFAULT_ROOT = "Data/GSE67520"
DEFAULT_SCORES = "results/Dynamics/z6_gse67462_mechanistic_interpretation/01_gene_multimodal_scores.csv"
DEFAULT_ASSIGNMENTS = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_OUT = "results/Dynamics/z6_gse67462_directional_mechanistic_audit"
ACTIVE = ["h3k27ac", "h3k4me3", "rnapii", "total_oct4"]
NEGATIVE = "h3k27me3"


def _num_time(x):
    try:
        v = float(x)
        return v if np.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def _norm_expr(expr_raw: pd.DataFrame, validated: set[str]) -> pd.DataFrame:
    keep = {}
    for c in expr_raw.columns:
        k = _norm_symbol(c)
        if k in validated and k not in keep:
            keep[k] = c
    out = expr_raw.loc[:, list(keep.values())].copy()
    out.columns = list(keep.keys())
    out.index = pd.to_numeric(out.index, errors="coerce")
    out = out.loc[out.index.notna()].copy()
    out.index = out.index.astype(float)
    return out.sort_index()


def _norm_values(values, validated):
    out = []
    for modality, time, mapping in values:
        t = _num_time(time)
        if t is None:
            continue
        d = {}
        for gene, value in mapping.items():
            k = _norm_symbol(gene)
            if k in validated:
                d[k] = d.get(k, 0.0) + float(value)
        out.append((modality, t, d))
    return out


def _frames(values, genes):
    frames = {}
    for modality in sorted({m for m, _, _ in values}):
        rows = [(t, d) for m, t, d in values if m == modality]
        rows.sort(key=lambda x: x[0])
        if not rows:
            continue
        frames[modality] = pd.DataFrame(
            [{g: d.get(g, 0.0) for g in genes} for _, d in rows],
            index=[t for t, _ in rows],
            columns=genes,
        )
    return frames


def _slope_sign(values):
    if len(values) < 3:
        return 0
    x = np.arange(len(values), dtype=float)
    y = np.asarray(values, dtype=float)
    if np.allclose(y, y[0]):
        return 0
    slope = np.polyfit(x, y, 1)[0]
    return int(np.sign(slope))


def _audit(expr, frames, assignments, scores):
    common_times = sorted(expr.index)
    module_map = assignments.set_index("gene")["module"].to_dict()
    score_map = scores.set_index("gene")
    rows = []
    for modality, frame in frames.items():
        for gene in expr.columns:
            times = [t for t in common_times if t in frame.index]
            if len(times) < 5:
                continue
            x_now = np.asarray([frame.loc[t, gene] for t in times], dtype=float)
            y_now = np.asarray([expr.loc[t, gene] for t in times], dtype=float)
            lead_times = [t for t in times if any(tt > t for tt in common_times)]
            next_times = []
            for t in lead_times:
                later = [tt for tt in common_times if tt > t]
                if later and later[0] in frame.index:
                    next_times.append((t, later[0]))
            if len(next_times) < 4:
                continue
            lx = np.asarray([frame.loc[t, gene] for t, _ in next_times], dtype=float)
            ly = np.asarray([expr.loc[tt, gene] for _, tt in next_times], dtype=float)
            concurrent = _spearman(x_now, y_now)
            lagged = _spearman(lx, ly)
            gain = lagged - concurrent if np.isfinite(lagged) and np.isfinite(concurrent) else np.nan
            mod_slope = _slope_sign(x_now)
            expr_next_slope = _slope_sign(ly)
            direction = int(mod_slope != 0 and expr_next_slope != 0 and mod_slope == expr_next_slope)
            core = False
            if gene in score_map.index:
                core = str(score_map.loc[gene].get("mechanistic_class", "")) == "MULTIMODAL_CORE"
            rows.append({
                "gene": gene,
                "module": int(module_map.get(gene, -1)),
                "modality": modality,
                "contemporaneous_rho": concurrent,
                "lagged_rho": lagged,
                "lead_gain": gain,
                "directional_agreement": direction,
                "modality_slope_sign": mod_slope,
                "subsequent_expression_slope_sign": expr_next_slope,
                "is_multimodal_core": core,
                "n_timepoints": len(times),
                "n_lag_pairs": len(next_times),
            })
    return pd.DataFrame(rows)


def _summary(audit):
    if audit.empty:
        return pd.DataFrame()
    a = audit.copy()
    a["lead_candidate"] = (a.lead_gain >= 0.10) & (a.lagged_rho >= 0.30)
    rows = []
    for (module, modality), g in a.groupby(["module", "modality"]):
        rows.append({
            "module": int(module),
            "modality": modality,
            "n_genes": len(g),
            "n_lead_candidates": int(g.lead_candidate.sum()),
            "lead_fraction": float(g.lead_candidate.mean()),
            "median_contemporaneous_rho": float(g.contemporaneous_rho.median()),
            "median_lagged_rho": float(g.lagged_rho.median()),
            "median_lead_gain": float(g.lead_gain.median()),
            "directional_agreement_fraction": float(g.directional_agreement.mean()),
            "core_fraction": float(g.is_multimodal_core.mean()),
            "candidate_core_fraction": float(g.loc[g.lead_candidate, "is_multimodal_core"].mean()) if g.lead_candidate.any() else np.nan,
        })
    return pd.DataFrame(rows).sort_values(["module", "lead_fraction"], ascending=[True, False])


def _global_summary(audit):
    a = audit.copy()
    a["lead_candidate"] = (a.lead_gain >= 0.10) & (a.lagged_rho >= 0.30)
    rows = []
    for modality, g in a.groupby("modality"):
        rows.append({
            "modality": modality,
            "n_genes": len(g),
            "lead_candidates": int(g.lead_candidate.sum()),
            "lead_fraction": float(g.lead_candidate.mean()),
            "median_lead_gain": float(g.lead_gain.median()),
            "median_lagged_rho": float(g.lagged_rho.median()),
            "directional_agreement_fraction": float(g.directional_agreement.mean()),
            "core_fraction": float(g.is_multimodal_core.mean()),
        })
    return pd.DataFrame(rows).sort_values("lead_fraction", ascending=False)


def _report(out, global_s, module_s, audit, raw_n, validated_n):
    a = audit.copy()
    a["lead_candidate"] = (a.lead_gain >= 0.10) & (a.lagged_rho >= 0.30)
    core = a[a.is_multimodal_core]
    active = a[a.modality.isin(ACTIVE)]
    neg = a[a.modality == NEGATIVE]
    active_lead = float(((active.lead_gain >= 0.10) & (active.lagged_rho >= 0.30)).mean()) if not active.empty else np.nan
    neg_lead = float(((neg.lead_gain >= 0.10) & (neg.lagged_rho >= 0.30)).mean()) if not neg.empty else np.nan
    lines = [
        "# GSE67462 Z6 — directional mechanistic ordering audit",
        "",
        f"Validated expression universe: **{validated_n:,}/{raw_n:,} genes**.",
        "",
        "The audit tests whether regulatory/chromatin-associated measurements at time t have stronger association with expression at the subsequent sampled time than with contemporaneous expression. It is diagnostic and does not modify frozen Z6 decisions.",
        "",
        "## Global modality summary",
        "",
        "| Modality | Genes | Lead candidates | Lead fraction | Median lead gain | Median lagged rho | Directional agreement |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in global_s.iterrows():
        lines.append(f"| {r.modality} | {int(r.n_genes)} | {int(r.lead_candidates)} | {r.lead_fraction:.3f} | {r.median_lead_gain:.3f} | {r.median_lagged_rho:.3f} | {r.directional_agreement_fraction:.3f} |")
    lines += [
        "",
        f"Active modalities combined lead-candidate fraction: **{active_lead:.3f}**.",
        f"H3K27me3 negative-control lead-candidate fraction: **{neg_lead:.3f}**.",
        "",
        "## M1 / M4 / M6 module summary",
        "",
        "| Module | Modality | Lead fraction | Median lead gain | Directional agreement | Core fraction | Candidate core fraction |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for _, r in module_s[module_s.module.isin([1,4,6])].iterrows():
        lines.append(f"| M{int(r.module)} | {r.modality} | {r.lead_fraction:.3f} | {r.median_lead_gain:.3f} | {r.directional_agreement_fraction:.3f} | {r.core_fraction:.3f} | {r.candidate_core_fraction:.3f} |")
    lines += [
        "",
        "## Interpretation",
        "",
        "A positive active-modality lag signal supports a temporal-ordering hypothesis: regulatory/chromatin state at an earlier sampled time is more concordant with later expression than with expression at the same time.",
        "",
        "Directional agreement is an additional trajectory-level check, not a causal test. The experiment has eight principal reprogramming timepoints and unequal intervals, so this audit cannot establish a continuous causal delay.",
        "",
        "The H3K27me3 comparison is a falsification control, but failure to separate it from active modalities does not prove that the active signal is non-biological; it would indicate that temporal ordering is not specific enough for mechanistic interpretation.",
        "",
        "## Scientific boundary",
        "",
        "The GSE67462 experiment was designed as a bulk time series with two replicates per stage and jointly profiled gene expression, Oct4 binding and histone modifications. The original study reports stage-specific regulatory signatures and enhancer activation preceding later pluripotency-network activation. This audit tests a narrower directional consequence of that model; it does not reproduce or prove the original mechanistic claims.",
        "",
        "## Frozen boundary",
        "",
        "No Z6 predictive support threshold, multimodal support decision, temporal module assignment or validation split was changed.",
    ]
    (out / "GSE67462_Z6_directional_mechanistic_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", default=DEFAULT_GTF)
    p.add_argument("--platform-soft", default=DEFAULT_SOFT)
    p.add_argument("--modality-root", default=DEFAULT_ROOT)
    p.add_argument("--scores", default=DEFAULT_SCORES)
    p.add_argument("--assignments", default=DEFAULT_ASSIGNMENTS)
    p.add_argument("--output", default=DEFAULT_OUT)
    args = p.parse_args()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    matrix, metadata = validation._load_common_space()
    tss = _load_tss(Path(args.gtf))
    platform = _read_soft_platform(Path(args.platform_soft))
    mapping, mapping_summary = _build_mapping_report(pd.Index(matrix.index.astype(str)), tss, platform)
    validated = set(mapping.loc[mapping["platform_symbol_tss_match"].fillna(0).astype(int) > 0, "expression_id"].map(_norm_symbol))
    expr_times, expr_raw, branches = _expression_common(matrix, metadata)
    expr = _norm_expr(expr_raw, validated)

    files = _discover_files(Path(args.modality_root))
    missing = [k for k in MODALITIES if not files[k]]
    if missing:
        raise RuntimeError(f"Missing modality files: {missing}")
    values_raw, provenance = _build_modality_matrix(files, tss)
    values = _norm_values(values_raw, validated)
    frames = _frames(values, expr.columns)
    scores = pd.read_csv(args.scores)
    scores["gene"] = scores["gene"].astype(str).map(_norm_symbol)
    assignments = pd.read_csv(args.assignments)
    assignments["gene"] = assignments["gene"].astype(str).map(_norm_symbol)

    audit = _audit(expr, frames, assignments, scores)
    global_s = _global_summary(audit)
    module_s = _summary(audit)

    audit.to_csv(out / "01_gene_directional_metrics.csv", index=False)
    global_s.to_csv(out / "02_modality_summary.csv", index=False)
    module_s.to_csv(out / "03_module_directional_summary.csv", index=False)
    (out / "04_mapping_summary.json").write_text(json.dumps(mapping_summary, indent=2, default=str), encoding="utf-8")
    manifest = {
        "validated_genes": int(expr.shape[1]),
        "raw_common_space_genes": int(expr_raw.shape[1]),
        "active_modalities": ACTIVE,
        "negative_control": NEGATIVE,
        "lead_gain_threshold": 0.10,
        "lagged_rho_threshold": 0.30,
        "diagnostic_only": True,
        "frozen_z6_support_unchanged": True,
    }
    (out / "05_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    _report(out, global_s, module_s, audit, int(expr_raw.shape[1]), int(expr.shape[1]))
    print("GSE67462 Z6 DIRECTIONAL MECHANISTIC ORDERING AUDIT")
    print(f"validated genes: {expr.shape[1]}/{expr_raw.shape[1]} ({expr.shape[1] / expr_raw.shape[1]:.2%})")
    print(f"outputs: {out}")
    print("frozen Z6 support unchanged: true")


if __name__ == "__main__":
    main()
