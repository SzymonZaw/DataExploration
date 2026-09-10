"""Run a bounded mechanistic falsification audit for GSE67462 Z6 modules.

This audit is diagnostic only. It compares competing interpretations of M1, M4 and
M6 using the already validated multimodal gene scores, temporal module assignments,
and specificity categories. It does not change frozen Z6 decisions.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_SCORES = "results/Dynamics/z6_gse67462_mechanistic_interpretation/01_gene_multimodal_scores.csv"
DEFAULT_ASSIGNMENTS = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_MODULES = "results/Dynamics/z6_gse67462_temporal_modules/02_temporal_modules.csv"
DEFAULT_SPECIFICITY = "results/Dynamics/z6_gse67462_specificity_audit/03_category_gene_coverage.csv"
DEFAULT_LEADING = "results/Dynamics/z6_gse67462_specificity_audit/02_leading_genes.csv"
DEFAULT_OUT = "results/Dynamics/z6_gse67462_mechanistic_falsification"

ACTIVE = ["h3k27ac", "h3k4me3", "rnapii", "total_oct4"]
COMPETING = [
    "ECM_MESENCHYMAL",
    "EPITHELIALIZATION",
    "PLURIPOTENCY_STEM",
    "PROLIFERATION_CELL_CYCLE",
    "STRESS_INFLAMMATION",
    "METABOLISM",
]


def _read_optional(path: Path, columns: list[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=columns)
    return pd.read_csv(path)


def _gene_score_table(scores: pd.DataFrame) -> pd.DataFrame:
    out = scores.copy()
    for m in ACTIVE:
        c = f"rho_{m}"
        if c not in out.columns:
            out[c] = np.nan
    out["active_count"] = sum((out[f"rho_{m}"] >= 0.30).fillna(False) for m in ACTIVE)
    out["active_median_rho"] = out[[f"rho_{m}" for m in ACTIVE]].where(lambda x: x >= 0).median(axis=1)
    return out


def _module_summary(assignments: pd.DataFrame, scores: pd.DataFrame, coverage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    s = scores.set_index("gene")
    for module, g in assignments.groupby("module"):
        genes = set(g["gene"].astype(str))
        sub = s.reindex(sorted(genes)).dropna(how="all")
        row = {"module": int(module), "n_genes": len(genes)}
        for m in ACTIVE:
            v = pd.to_numeric(sub[f"rho_{m}"], errors="coerce")
            row[f"median_rho_{m}"] = float(v.median()) if v.notna().any() else np.nan
            row[f"fraction_rho_{m}_ge_0_30"] = float((v >= 0.30).mean()) if v.notna().any() else np.nan
        row["fraction_3plus_active_modalities"] = float((sub["active_count"] >= 3).mean()) if len(sub) else np.nan
        for cat in COMPETING:
            c = coverage[(coverage.module == int(module)) & (coverage.category == cat)]
            row[f"coverage_{cat}"] = float(c.iloc[0].leading_fraction_of_module) if not c.empty else 0.0
            row[f"coverage_rho_{cat}"] = float(c.iloc[0].median_multimodal_rho) if not c.empty else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values("module")


def _m6_score(summary: pd.DataFrame) -> pd.DataFrame:
    r = summary[summary.module == 6].copy()
    if r.empty:
        return pd.DataFrame()
    row = r.iloc[0]
    # These are deliberately transparent diagnostic scores, not probabilities.
    # Category coverage is weighted by multimodal coherence; OCT4 is reported separately.
    scores = []
    for cat in COMPETING:
        cov = float(row.get(f"coverage_{cat}", 0.0))
        rho = row.get(f"coverage_rho_{cat}", np.nan)
        coherence = float(rho) if np.isfinite(rho) else 0.0
        score = cov * max(0.0, coherence)
        scores.append({"module": 6, "hypothesis": cat, "coverage_fraction": cov, "coherence_rho": coherence, "diagnostic_score": score})
    out = pd.DataFrame(scores).sort_values("diagnostic_score", ascending=False)
    oct4 = float(row.get("median_rho_total_oct4", np.nan))
    out["module_median_oct4_rho"] = oct4
    return out


def _leading_contrast(leading: pd.DataFrame) -> pd.DataFrame:
    if leading.empty:
        return pd.DataFrame(columns=["module", "category", "genes", "mean_rank_score", "median_rho", "fraction_3plus"])
    rows = []
    for (module, category), g in leading.groupby(["module", "category"]):
        rows.append({
            "module": int(module),
            "category": category,
            "genes": int(g.gene.nunique()),
            "mean_rank_score": float(g.rank_score.mean()),
            "median_rho": float(g.median_positive_active_rho.median()),
            "fraction_3plus": float((g.n_active_modalities_ge_0_30 >= 3).mean()),
        })
    return pd.DataFrame(rows).sort_values(["module", "mean_rank_score"], ascending=[True, False])


def _decisions(summary: pd.DataFrame, m6: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for module in [1, 4, 6]:
        r = summary[summary.module == module]
        if r.empty:
            continue
        row = r.iloc[0]
        vals = {cat: float(row.get(f"coverage_{cat}", 0.0)) for cat in COMPETING}
        ranked = sorted(vals.items(), key=lambda x: x[1], reverse=True)
        top_cat, top_cov = ranked[0]
        second_cov = ranked[1][1] if len(ranked) > 1 else 0.0
        core = float(row.get("fraction_3plus_active_modalities", 0.0))
        if module == 6:
            # Do not call pluripotency supported unless it is competitive with epithelialization.
            epi = vals.get("EPITHELIALIZATION", 0.0)
            plur = vals.get("PLURIPOTENCY_STEM", 0.0)
            if epi > 0 and plur > 0 and epi >= plur * 1.5:
                decision = "MIXED_EVIDENCE_EPITHELIALIZATION_DOMINANT"
            elif plur > 0 and plur >= epi * 0.67 and core >= 0.50:
                decision = "UNRESOLVED_PLURIPOTENCY_VS_EPITHELIALIZATION"
            else:
                decision = "MIXED_EVIDENCE"
        elif core >= 0.50 and top_cov >= 0.15 and top_cov >= second_cov * 1.15:
            decision = "SUPPORTED_AS_HYPOTHESIS"
        else:
            decision = "MIXED_EVIDENCE"
        rows.append({"module": module, "top_category": top_cat, "top_coverage": top_cov, "second_coverage": second_cov, "multimodal_core_fraction": core, "decision": decision})
    return pd.DataFrame(rows)


def _write_report(out: Path, summary: pd.DataFrame, contrast: pd.DataFrame, m6: pd.DataFrame, decisions: pd.DataFrame) -> None:
    lines = [
        "# GSE67462 Z6 — mechanistic falsification audit",
        "",
        "This audit compares competing biological interpretations of M1, M4 and M6 using the frozen multimodal core and the previously generated specificity audit. It is diagnostic only.",
        "",
        "## Decision summary",
        "",
        "| Module | Dominant diagnostic category | Coverage | 2nd category | Multimodal core fraction | Decision |",
        "|---:|---|---:|---|---:|---|",
    ]
    for _, r in decisions.iterrows():
        lines.append(f"| M{int(r.module)} | {r.top_category} | {r.top_coverage:.3f} | {r.second_coverage:.3f} | {r.multimodal_core_fraction:.3f} | {r.decision} |")
    lines += ["", "## M6 competing hypotheses", ""]
    if m6.empty:
        lines.append("M6 scorecard unavailable.")
    else:
        lines += ["| Hypothesis | Coverage | Multimodal coherence | Diagnostic score |", "|---|---:|---:|---:|"]
        for _, r in m6.iterrows():
            lines.append(f"| {r.hypothesis} | {r.coverage_fraction:.3f} | {r.coherence_rho:.3f} | {r.diagnostic_score:.3f} |")
        oct4 = float(m6.module_median_oct4_rho.iloc[0])
        lines += ["", f"M6 median OCT4-module correlation: **{oct4:.3f}**.", ""]
        lines.append("OCT4 correlation is treated as an associated transition signal, not proof of pluripotency.")
    lines += ["", "## Leading-gene contrast", "", "| Module | Category | Genes | Mean rank score | Median rho | Fraction 3+ modalities |", "|---:|---|---:|---:|---:|---:|"]
    for _, r in contrast.head(60).iterrows():
        lines.append(f"| M{int(r.module)} | {r.category} | {int(r.genes)} | {r.mean_rank_score:.3f} | {r.median_rho:.3f} | {r.fraction_3plus:.3f} |")
    lines += [
        "",
        "## Falsification interpretation",
        "",
        "### M1",
        "",
        "M1 is considered supported only as an ECM/mesenchymal-remodeling hypothesis. The audit does not distinguish EMT from MET direction because the module contains genes associated with both epithelial and mesenchymal processes.",
        "",
        "### M4",
        "",
        "M4 is treated as a structural-remodeling hypothesis rather than a generic proliferation program. RHO/RAC/RAB and integrin/ECM enrichment is interpreted as convergent context, not causal evidence.",
        "",
        "### M6",
        "",
        "M6 is deliberately not labelled a pluripotency module. Its late trajectory and OCT4 association are compatible with reprogramming, but epithelialization/cornification and cholesterol metabolism provide competing explanations. The appropriate result is therefore a bounded competing-hypothesis statement rather than a single mechanistic label.",
        "",
        "## Non-claims",
        "",
        "This audit does not establish causality, perturbation effects, pathway necessity, or transferability across datasets. Those require independent perturbation or intervention data and/or prospective validation.",
        "",
        "## Frozen boundary",
        "",
        "No Z6 predictive support threshold, multimodal support decision, temporal module assignment or validation split was changed.",
    ]
    (out / "GSE67462_Z6_mechanistic_falsification_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scores", default=DEFAULT_SCORES)
    p.add_argument("--assignments", default=DEFAULT_ASSIGNMENTS)
    p.add_argument("--modules", default=DEFAULT_MODULES)
    p.add_argument("--specificity", default=DEFAULT_SPECIFICITY)
    p.add_argument("--leading", default=DEFAULT_LEADING)
    p.add_argument("--output", default=DEFAULT_OUT)
    args = p.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)

    scores = _gene_score_table(pd.read_csv(args.scores))
    assignments = pd.read_csv(args.assignments)
    coverage = _read_optional(Path(args.specificity), ["module", "category", "leading_genes", "module_size", "leading_fraction_of_module", "median_multimodal_rho", "fraction_with_3plus_modalities"])
    leading = _read_optional(Path(args.leading), ["module", "category", "gene", "significant_term_count", "median_positive_active_rho", "n_active_modalities_ge_0_30", "rank_score"])

    summary = _module_summary(assignments, scores, coverage)
    contrast = _leading_contrast(leading)
    m6 = _m6_score(summary)
    decisions = _decisions(summary, m6)

    summary.to_csv(out / "01_module_mechanistic_summary.csv", index=False)
    contrast.to_csv(out / "02_leading_gene_contrast.csv", index=False)
    m6.to_csv(out / "03_m6_competing_hypotheses.csv", index=False)
    decisions.to_csv(out / "04_mechanistic_decisions.csv", index=False)
    manifest = {
        "modules_tested": [1, 4, 6],
        "active_modalities": ACTIVE,
        "negative_control_excluded_from_active_score": "h3k27me3",
        "diagnostic_only": True,
        "frozen_z6_support_unchanged": True,
        "interpretation": "MECHANISTIC_FALSIFICATION_AUDIT",
    }
    (out / "05_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    _write_report(out, summary, contrast, m6, decisions)

    print("GSE67462 Z6 MECHANISTIC FALSIFICATION AUDIT")
    print(f"modules tested: {manifest['modules_tested']}")
    print(f"outputs: {out}")
    print("frozen Z6 support unchanged: true")


if __name__ == "__main__":
    main()
