"""Run a diagnostic mechanistic interpretation audit for GSE67462/GSE67520.

The audit ranks genes by multimodal temporal coherence, tests a hypothesis-generating
adjacent-time regulatory lead signal, and optionally performs versioned GMT enrichment.
It never changes the frozen Z6 predictive or multimodal support decisions.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from dynamics import validation
from dynamics.run_z6_gse67462_identifier_mapping_audit import (
    _build_mapping_report,
    _norm_symbol,
    _read_soft_platform,
)
from dynamics.run_z6_gse67462_multimodal_validation_analysis import (
    DAY_TO_HOURS,
    MODALITIES,
    SEED,
    _build_modality_matrix,
    _discover_files,
    _expression_common,
    _load_tss,
    _spearman,
)

DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_SOFT = "Data/GPL19972_family.soft.gz"
DEFAULT_GENE_SETS = None
DEFAULT_OUT = "results/Dynamics/z6_gse67462_mechanistic_interpretation"
ACTIVE_MODALITIES = ["h3k27ac", "h3k4me3", "rnapii", "total_oct4"]
NEGATIVE_MODALITY = "h3k27me3"


def _numeric_time(value):
    try:
        x = float(value)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _normalized_expression(expr_raw: pd.DataFrame, validated: set[str]) -> pd.DataFrame:
    keep: dict[str, str] = {}
    for col in expr_raw.columns:
        key = _norm_symbol(col)
        if key in validated and key not in keep:
            keep[key] = col
    out = expr_raw.loc[:, list(keep.values())].copy()
    out.columns = list(keep.keys())
    out.index = pd.to_numeric(out.index, errors="coerce")
    out = out.loc[out.index.notna()].copy()
    out.index = out.index.astype(float)
    return out.sort_index()


def _normalize_modality_values(values, validated: set[str]):
    result = []
    for modality, time, mapping in values:
        t = _numeric_time(time)
        if t is None:
            continue
        out: dict[str, float] = {}
        for gene, value in mapping.items():
            key = _norm_symbol(gene)
            if key in validated:
                out[key] = out.get(key, 0.0) + float(value)
        result.append((modality, t, out))
    return result


def _modality_frames(values, genes: pd.Index) -> dict[str, pd.DataFrame]:
    frames: dict[str, pd.DataFrame] = {}
    for modality in sorted(set(m for m, _, _ in values)):
        entries = sorted([(t, v) for m, t, v in values if m == modality], key=lambda x: x[0])
        if not entries:
            continue
        rows = [{g: v.get(g, 0.0) for g in genes} for _, v in entries]
        frames[modality] = pd.DataFrame(rows, index=[t for t, _ in entries], columns=genes)
    return frames


def _gene_level_scores(expr: pd.DataFrame, frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    genes = expr.columns
    for gene in genes:
        r: dict[str, object] = {"gene": gene}
        active = []
        all_rhos = []
        for modality, frame in frames.items():
            times = sorted(set(expr.index).intersection(frame.index))
            if len(times) < 5:
                rho = np.nan
            else:
                rho = _spearman(expr.loc[times, gene].to_numpy(), frame.loc[times, gene].to_numpy())
            r[f"rho_{modality}"] = rho
            if np.isfinite(rho):
                all_rhos.append(rho)
            if modality in ACTIVE_MODALITIES and np.isfinite(rho) and rho >= 0.30:
                active.append(modality)
        positive = [r.get(f"rho_{m}") for m in ACTIVE_MODALITIES if np.isfinite(r.get(f"rho_{m}", np.nan)) and r[f"rho_{m}"] > 0]
        positive = [float(x) for x in positive]
        r["n_active_modalities_ge_0_30"] = len(active)
        r["active_modalities_ge_0_30"] = ";".join(sorted(active))
        r["median_positive_active_rho"] = float(np.median(positive)) if positive else np.nan
        r["mean_active_rho"] = float(np.mean([r[f"rho_{m}"] for m in ACTIVE_MODALITIES if np.isfinite(r.get(f"rho_{m}", np.nan))])) if any(np.isfinite(r.get(f"rho_{m}", np.nan)) for m in ACTIVE_MODALITIES) else np.nan
        if len(active) >= 3 and r["median_positive_active_rho"] >= 0.30:
            cls = "MULTIMODAL_CORE"
        elif len(active) >= 2:
            cls = "MULTIMODAL_SUPPORT"
        elif len(active) == 1:
            cls = "SINGLE_MODAL"
        else:
            cls = "UNRESOLVED"
        r["mechanistic_class"] = cls
        rows.append(r)
    return pd.DataFrame(rows)


def _lead_candidates(expr: pd.DataFrame, frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    common_times = sorted(expr.index)
    for modality in ACTIVE_MODALITIES:
        if modality not in frames:
            continue
        frame = frames[modality]
        for gene in expr.columns:
            concurrent_x, concurrent_y, lead_x, lead_y = [], [], [], []
            for i, t in enumerate(common_times):
                if t not in frame.index:
                    continue
                concurrent_x.append(frame.loc[t, gene])
                concurrent_y.append(expr.loc[t, gene])
                if i + 1 < len(common_times) and common_times[i + 1] in frame.index:
                    lead_x.append(frame.loc[t, gene])
                    lead_y.append(expr.loc[common_times[i + 1], gene])
            concurrent = _spearman(np.asarray(concurrent_x), np.asarray(concurrent_y)) if len(concurrent_x) >= 5 else np.nan
            lead = _spearman(np.asarray(lead_x), np.asarray(lead_y)) if len(lead_x) >= 4 else np.nan
            delta = lead - concurrent if np.isfinite(lead) and np.isfinite(concurrent) else np.nan
            if np.isfinite(delta) and delta >= 0.10 and np.isfinite(lead) and lead >= 0.30:
                rows.append({"gene": gene, "modality": modality, "concurrent_rho": concurrent, "lead_rho": lead, "lead_minus_concurrent": delta, "n_pairs": len(lead_x)})
    if not rows:
        return pd.DataFrame(columns=["gene", "modality", "concurrent_rho", "lead_rho", "lead_minus_concurrent", "n_pairs"])
    return pd.DataFrame(rows).sort_values(["lead_minus_concurrent", "lead_rho"], ascending=False)


def _read_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = {_norm_symbol(g) for g in parts[2:] if _norm_symbol(g)}
            if genes:
                sets[name] = genes
    return sets


def _enrichment(core: set[str], universe: set[str], gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    M = len(universe)
    N = len(core)
    rows = []
    if not core or not universe:
        return pd.DataFrame(columns=["term", "overlap", "core_size", "term_size", "universe_size", "p_value", "fdr_bh"])
    for term, genes in gene_sets.items():
        gs = genes & universe
        k = len(core & gs)
        if k == 0:
            continue
        p = float(hypergeom.sf(k - 1, M, len(gs), N))
        rows.append({"term": term, "overlap": k, "core_size": N, "term_size": len(gs), "universe_size": M, "p_value": p})
    if not rows:
        return pd.DataFrame(columns=["term", "overlap", "core_size", "term_size", "universe_size", "p_value", "fdr_bh"])
    out = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    m = len(out)
    out["fdr_bh"] = (out["p_value"] * m / (np.arange(m) + 1)).clip(upper=1.0)
    out["fdr_bh"] = np.minimum.accumulate(out["fdr_bh"][::-1])[::-1]
    return out


def _write_report(out: Path, scores: pd.DataFrame, leads: pd.DataFrame, enrichment: pd.DataFrame, nuisance: pd.DataFrame, raw_n: int, validated_n: int) -> None:
    core = scores[scores.mechanistic_class == "MULTIMODAL_CORE"].copy()
    lead_genes = leads["gene"].nunique() if not leads.empty else 0
    lines = [
        "# Z6 GSE67462 — mechanistic interpretation audit",
        "",
        "## Finding",
        "",
        f"The validated GSE67462 expression space contains **{validated_n:,} of {raw_n:,}** common-space genes. "
        f"The multimodal core ranking identifies **{len(core):,} genes** with positive temporal agreement (rho ≥ 0.30) across at least three of the four supported active modalities: H3K27ac, H3K4me3, RNAPII and OCT4.",
        "",
        "This is evidence of multimodal molecular coherence, not causal proof. H3K27me3 is excluded from the active core score and retained as a contrasting/negative-control modality.",
        "",
        "## Core candidate genes",
        "",
        "| Rank | Gene | Active modalities ≥0.30 | Median positive rho | Class |",
        "|---:|---|---|---:|---|",
    ]
    top = core.sort_values(["median_positive_active_rho", "n_active_modalities_ge_0_30"], ascending=False).head(50)
    for i, (_, row) in enumerate(top.iterrows(), 1):
        lines.append(f"| {i} | {row.gene} | {row.active_modalities_ge_0_30} | {row.median_positive_active_rho:.3f} | {row.mechanistic_class} |")
    lines += [
        "",
        "## Regulatory-lead diagnostic",
        "",
        f"The adjacent-time diagnostic identified **{lead_genes} genes** where an active regulatory modality had a subsequent-expression correlation at least 0.10 higher than its contemporaneous correlation and lead rho ≥ 0.30. With only eight principal timepoints and unequal intervals, this is hypothesis-generating rather than evidence of causal ordering.",
        "",
        "## Functional enrichment",
        "",
    ]
    if enrichment.empty:
        lines.append("No GMT gene-set file was supplied or no enriched terms were detected. No pathway claim is made from this audit.")
    else:
        lines.append("| Term | Overlap | p-value | FDR BH |")
        lines.append("|---|---:|---:|---:|")
        for _, row in enrichment.head(30).iterrows():
            lines.append(f"| {row.term} | {int(row.overlap)} | {row.p_value:.3g} | {row.fdr_bh:.3g} |")
    lines += [
        "",
        "## Nuisance/context audit",
        "",
    ]
    if nuisance.empty:
        lines.append("No nuisance/context GMT was supplied. Biological specificity is therefore unresolved.")
    else:
        lines.append("Nuisance/context enrichment is reported diagnostically; enrichment or depletion does not by itself establish or refute biological specificity.")
        lines.append("")
        lines.append("| Term | Overlap | p-value | FDR BH |")
        lines.append("|---|---:|---:|---:|")
        for _, row in nuisance.head(30).iterrows():
            lines.append(f"| {row.term} | {int(row.overlap)} | {row.p_value:.3g} | {row.fdr_bh:.3g} |")
    lines += [
        "",
        "## Interpretation boundary",
        "",
        "The appropriate scientific statement is: **the GSE67462 transition contains a multimodally coherent molecular core that is suitable for mechanistic hypothesis generation.** The analysis does not establish that OCT4, chromatin remodeling or any individual pathway causes the transition, and it does not establish transferability to other datasets.",
        "",
        "## Reproducibility",
        "",
        "All candidate rankings are generated from machine-readable local inputs. Optional enrichment uses an explicitly supplied GMT file; no external pathway database is silently downloaded.",
    ]
    (out / "GSE67462_Z6_mechanistic_interpretation_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--platform-soft", default=DEFAULT_SOFT)
    parser.add_argument("--modality-root", default="Data/GSE67520")
    parser.add_argument("--gene-sets", default=DEFAULT_GENE_SETS, help="Optional versioned GMT file for pathway enrichment")
    parser.add_argument("--nuisance-gene-sets", default=None, help="Optional versioned GMT file for nuisance/context enrichment")
    parser.add_argument("--output", default=DEFAULT_OUT)
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    matrix, metadata = validation._load_common_space()
    tss = _load_tss(Path(args.gtf))
    platform = _read_soft_platform(Path(args.platform_soft))
    mapping, mapping_summary = _build_mapping_report(pd.Index(matrix.index.astype(str)), tss, platform)
    validated = set(mapping.loc[mapping["platform_symbol_tss_match"].fillna(0).astype(int) > 0, "expression_id"].map(_norm_symbol))
    expr_times, expr_raw, branches = _expression_common(matrix, metadata)
    expr = _normalized_expression(expr_raw, validated)

    files = _discover_files(Path(args.modality_root))
    missing = [k for k in MODALITIES if not files[k]]
    if missing:
        raise RuntimeError(f"Missing modality files: {missing}")
    values_raw, provenance = _build_modality_matrix(files, tss)
    values = _normalize_modality_values(values_raw, validated)
    frames = _modality_frames(values, expr.columns)
    scores = _gene_level_scores(expr, frames)
    leads = _lead_candidates(expr, frames)

    core_genes = set(scores.loc[scores.mechanistic_class == "MULTIMODAL_CORE", "gene"])
    universe = set(expr.columns)
    enrichment = pd.DataFrame()
    nuisance = pd.DataFrame()
    if args.gene_sets:
        enrichment = _enrichment(core_genes, universe, _read_gmt(Path(args.gene_sets)))
    if args.nuisance_gene_sets:
        nuisance = _enrichment(core_genes, universe, _read_gmt(Path(args.nuisance_gene_sets)))

    scores.sort_values(["mechanistic_class", "median_positive_active_rho"], ascending=[True, False]).to_csv(out / "01_gene_multimodal_scores.csv", index=False)
    leads.to_csv(out / "02_regulatory_lead_candidates.csv", index=False)
    enrichment.to_csv(out / "03_pathway_enrichment.csv", index=False)
    nuisance.to_csv(out / "04_nuisance_enrichment.csv", index=False)
    provenance.to_csv(out / "05_peak_assignment_provenance.csv", index=False)

    summary = {
        "dataset_expression": "GSE67462",
        "dataset_multimodal": "GSE67520",
        "raw_common_space_genes": int(expr_raw.shape[1]),
        "validated_genes": int(expr.shape[1]),
        "validated_fraction": float(expr.shape[1] / expr_raw.shape[1]) if expr_raw.shape[1] else 0.0,
        "multimodal_core_genes": int(len(core_genes)),
        "regulatory_lead_candidates": int(leads["gene"].nunique()) if not leads.empty else 0,
        "active_modalities": ACTIVE_MODALITIES,
        "negative_control": NEGATIVE_MODALITY,
        "gene_set_source": str(args.gene_sets) if args.gene_sets else None,
        "nuisance_gene_set_source": str(args.nuisance_gene_sets) if args.nuisance_gene_sets else None,
        "mapping_summary": mapping_summary,
        "frozen_z6_support_unchanged": True,
        "interpretation": "MULTIMODAL_MECHANISTIC_HYPOTHESIS_GENERATION",
    }
    (out / "06_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    _write_report(out, scores, leads, enrichment, nuisance, int(expr_raw.shape[1]), int(expr.shape[1]))

    print("GSE67462 MECHANISTIC INTERPRETATION AUDIT")
    print(f"validated genes: {expr.shape[1]}/{expr_raw.shape[1]} ({expr.shape[1] / expr_raw.shape[1]:.2%})")
    print(f"multimodal core genes: {len(core_genes)}")
    print(f"regulatory lead candidates: {leads['gene'].nunique() if not leads.empty else 0}")
    print("negative control: H3K27me3 (excluded from active core score)")
    print("interpretation: MULTIMODAL_MECHANISTIC_HYPOTHESIS_GENERATION")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
