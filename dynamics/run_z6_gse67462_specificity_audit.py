"""Audit biological specificity of GSE67462 temporal modules.

Diagnostic only. Combines temporal-module assignments, multimodal gene scores,
and versioned ORA results to distinguish ECM/mesenchymal, epithelialization,
pluripotency, proliferation, stress/inflammation and metabolism signals.
No frozen Z6 support decision is modified.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_ENRICHMENT = "results/Dynamics/z6_gse67462_temporal_enrichment/01_enrichment_all.csv"
DEFAULT_SCORES = "results/Dynamics/z6_gse67462_mechanistic_interpretation/01_gene_multimodal_scores.csv"
DEFAULT_ASSIGNMENTS = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_GO = "Data/GeneSets/GO_BP_mouse.gmt"
DEFAULT_REACTOME = "Data/GeneSets/Reactome_mouse.gmt"
DEFAULT_HALLMARK = "Data/GeneSets/Hallmark_mouse.gmt"
DEFAULT_OUT = "results/Dynamics/z6_gse67462_specificity_audit"

CATEGORIES = {
    "ECM_MESENCHYMAL": [
        r"extracellular.?matrix", r"cell.?adhesion", r"wound", r"tissue.?remodel", r"mesenchym",
        r"integrin", r"collagen", r"fibroblast", r"RHO", r"RAC", r"RAB", r"growth.?factor",
        r"receptor tyrosine kinase", r"EMT",
    ],
    "EPITHELIALIZATION": [
        r"epithelial", r"keratin", r"cornified", r"apical.?surface", r"cell.?junction",
        r"epiderm", r"epithelial.?mesenchymal",
    ],
    "PLURIPOTENCY_STEM": [
        r"stem.?cell", r"pluripot", r"embryonic", r"reprogram", r"self.?renew", r"cell.?fate",
        r"developmental potential",
    ],
    "PROLIFERATION_CELL_CYCLE": [
        r"cell.?cycle", r"mitotic", r"DNA replication", r"pre.?replicative", r"proliferat",
        r"chromosome segregation", r"S phase", r"G2/M",
    ],
    "STRESS_INFLAMMATION": [
        r"stress", r"inflamm", r"TNFA", r"NF.?kB", r"hypoxia", r"apoptosis", r"interferon",
        r"response to cytokine", r"UV_RESPONSE",
    ],
    "METABOLISM": [
        r"metabolism", r"glycolysis", r"cholesterol", r"lipid", r"fatty acid", r"steroid",
        r"oxidative phosphorylation", r"respiratory", r"mitochond",
    ],
}


def norm(x: object) -> str:
    return re.sub(r"[^A-Z0-9]+", "_", str(x).upper()).strip("_")


def read_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            p = line.rstrip("\r\n").split("\t")
            if len(p) >= 3:
                genes = {str(x).strip().upper() for x in p[2:] if str(x).strip()}
                if genes:
                    sets[p[0].strip()] = genes
    return sets


def assign_categories(term: str) -> list[str]:
    n = norm(term)
    return [cat for cat, patterns in CATEGORIES.items() if any(re.search(p, n, flags=re.I) for p in patterns)]


def bh(p: pd.Series) -> pd.Series:
    a = p.to_numpy(float)
    order = np.argsort(a)
    q = a[order] * len(a) / np.arange(1, len(a) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return pd.Series(out, index=p.index)


def term_categories(enrichment: pd.DataFrame) -> pd.DataFrame:
    e = enrichment.copy()
    e["category"] = e["term"].map(assign_categories)
    e = e.explode("category").dropna(subset=["category"])
    return e


def build_category_summary(e: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (module, category), g in e.groupby(["module", "category"]):
        sig = g[g["fdr_bh_global"] < 0.05].copy()
        sig_module = g[g["fdr_bh_module"] < 0.05].copy()
        rows.append({
            "module": int(module),
            "category": category,
            "all_matching_terms": len(g),
            "global_fdr_significant_terms": len(sig),
            "module_fdr_significant_terms": len(sig_module),
            "best_global_fdr": float(g["fdr_bh_global"].min()),
            "best_p": float(g["p_value"].min()),
            "max_overlap": int(g["overlap"].max()),
            "total_overlap_counts": int(g["overlap"].sum()),
        })
    return pd.DataFrame(rows).sort_values(["module", "best_global_fdr", "category"])


def load_term_genes(paths: list[Path]) -> dict[str, set[str]]:
    merged: dict[str, set[str]] = {}
    for p in paths:
        if not p.exists():
            raise FileNotFoundError(p)
        merged.update(read_gmt(p))
    return merged


def leading_genes(enrichment: pd.DataFrame, scores: pd.DataFrame, assignments: pd.DataFrame, gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    # Use only terms with global FDR < 0.05 and module FDR < 0.05; then rank genes by
    # the number of independent significant terms covering them plus multimodal coherence.
    sig = enrichment[(enrichment.fdr_bh_global < 0.05) & (enrichment.fdr_bh_module < 0.05)].copy()
    score_cols = [c for c in ["median_positive_active_rho", "mean_active_rho", "n_active_modalities_ge_0_30"] if c in scores.columns]
    score = scores.set_index("gene")
    module_genes = {int(m): set(g["gene"].astype(str).str.upper()) for m, g in assignments.groupby("module")}
    rows = []
    for module in sorted(module_genes):
        msig = sig[sig.module == module]
        if msig.empty:
            continue
        for category in CATEGORIES:
            terms = msig[msig.term.map(lambda x: category in assign_categories(x))]
            if terms.empty:
                continue
            genes: dict[str, int] = {}
            for term in terms.term:
                for gene in gene_sets.get(term, set()):
                    if gene in module_genes[module]:
                        genes[gene] = genes.get(gene, 0) + 1
            for gene, n_terms in genes.items():
                if gene not in score.index:
                    continue
                row = score.loc[gene]
                rows.append({
                    "module": module,
                    "category": category,
                    "gene": gene,
                    "significant_term_count": n_terms,
                    "median_positive_active_rho": float(row.get("median_positive_active_rho", np.nan)),
                    "n_active_modalities_ge_0_30": int(row.get("n_active_modalities_ge_0_30", 0)),
                })
    out = pd.DataFrame(rows)
    if out.empty:
        return pd.DataFrame(columns=["module", "category", "gene", "significant_term_count", "median_positive_active_rho", "n_active_modalities_ge_0_30", "rank_score"])
    out["rank_score"] = (
        out["significant_term_count"]
        * out["median_positive_active_rho"].fillna(0)
        * (1 + 0.25 * out["n_active_modalities_ge_0_30"])
    )
    return out.sort_values(["module", "category", "rank_score"], ascending=[True, True, False]).reset_index(drop=True)


def category_gene_coverage(leading: pd.DataFrame, assignments: pd.DataFrame, scores: pd.DataFrame) -> pd.DataFrame:
    if leading.empty:
        return pd.DataFrame()
    module_sizes = assignments.groupby("module")["gene"].nunique().to_dict()
    rows = []
    for (m, cat), g in leading.groupby(["module", "category"]):
        genes = set(g.gene)
        sg = scores[scores.gene.isin(genes)]
        rows.append({
            "module": int(m),
            "category": cat,
            "leading_genes": len(genes),
            "module_size": int(module_sizes.get(m, 0)),
            "leading_fraction_of_module": len(genes) / max(1, module_sizes.get(m, 0)),
            "median_multimodal_rho": float(sg.median_positive_active_rho.median()) if not sg.empty else np.nan,
            "fraction_with_3plus_modalities": float((sg.n_active_modalities_ge_0_30 >= 3).mean()) if not sg.empty else np.nan,
        })
    return pd.DataFrame(rows).sort_values(["module", "leading_fraction_of_module"], ascending=[True, False])


def write_report(out: Path, summary: pd.DataFrame, coverage: pd.DataFrame, leading: pd.DataFrame, raw_enrichment: pd.DataFrame) -> None:
    lines = [
        "# GSE67462 Z6 — temporal module specificity audit",
        "",
        "This diagnostic audit evaluates whether temporal-module enrichment is better described as ECM/mesenchymal remodeling, epithelialization, pluripotency/stem-cell programs, proliferation/cell-cycle, stress/inflammation, or metabolism.",
        "",
        "The categories are keyword-based hypothesis screens over explicitly versioned GO-BP/Reactome/Hallmark results. They are not causal annotations and do not alter frozen Z6 decisions.",
        "",
        "## Category summary",
        "",
        "| Module | Category | Significant global terms | Best global FDR | Max overlap |",
        "|---:|---|---:|---:|---:|",
    ]
    if summary.empty:
        lines.append("| — | No category-matching terms | — | — | — |")
    else:
        for _, r in summary.iterrows():
            lines.append(f"| M{int(r.module)} | {r.category} | {int(r.global_fdr_significant_terms)} | {r.best_global_fdr:.3g} | {int(r.max_overlap)} |")
    lines += ["", "## Leading genes", "", "Leading genes are ranked by the number of significant terms covering the gene and its multimodal coherence score.", "", "| Module | Category | Gene | Significant terms | Modalities ≥0.30 | Median positive rho |", "|---:|---|---|---:|---:|---:|"]
    if leading.empty:
        lines.append("| — | — | No leading genes | — | — | — |")
    else:
        for _, r in leading.groupby(["module", "category"], sort=True).head(10).iterrows():
            lines.append(f"| M{int(r.module)} | {r.category} | {r.gene} | {int(r.significant_term_count)} | {int(r.n_active_modalities_ge_0_30)} | {r.median_positive_active_rho:.3f} |")
    lines += ["", "## Category coverage", "", "| Module | Category | Leading genes | Fraction of module | Median multimodal rho | Fraction 3+ modalities |", "|---:|---|---:|---:|---:|---:|"]
    if coverage.empty:
        lines.append("| — | — | — | — | — | — |")
    else:
        for _, r in coverage.iterrows():
            lines.append(f"| M{int(r.module)} | {r.category} | {int(r.leading_genes)} | {r.leading_fraction_of_module:.3f} | {r.median_multimodal_rho:.3f} | {r.fraction_with_3plus_modalities:.3f} |")
    lines += [
        "",
        "## Scientific interpretation boundary",
        "",
        "M1 should be treated primarily as an ECM/mesenchymal remodeling and receptor/growth-factor signaling module; its very strong EMT Hallmark enrichment is supportive of similarity to that state program, not proof of EMT/MET direction.",
        "",
        "M4 should be treated as a late remodeling/trafficking/cytoskeletal module, given RHO/RAC/RAB, integrin and ECM signals.",
        "",
        "M6 has evidence compatible with late epithelialization and OCT4-associated transition, but keratinization/cornified-envelope and metabolic terms make a specific pluripotency claim premature.",
        "",
        "M3 has weaker and less coherent pathway evidence and should remain exploratory.",
        "",
        "Generic proliferation, stress/inflammation and metabolism are reported as specificity qualifiers rather than filtered away. The absence of a category does not prove absence of the biological program.",
        "",
        "## Frozen-analysis boundary",
        "",
        "This audit does not change the Z6 predictive support criterion, multimodal support decision, module assignments, or temporal clustering.",
    ]
    (out / "GSE67462_Z6_specificity_audit_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--enrichment", default=DEFAULT_ENRICHMENT)
    p.add_argument("--scores", default=DEFAULT_SCORES)
    p.add_argument("--assignments", default=DEFAULT_ASSIGNMENTS)
    p.add_argument("--go", default=DEFAULT_GO)
    p.add_argument("--reactome", default=DEFAULT_REACTOME)
    p.add_argument("--hallmark", default=DEFAULT_HALLMARK)
    p.add_argument("--output", default=DEFAULT_OUT)
    args = p.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    enrichment = pd.read_csv(args.enrichment)
    scores = pd.read_csv(args.scores)
    assignments = pd.read_csv(args.assignments)
    assignments["gene"] = assignments["gene"].astype(str).str.upper()
    scores["gene"] = scores["gene"].astype(str).str.upper()

    required = {"collection", "module", "term", "overlap", "p_value", "fdr_bh_global", "fdr_bh_module"}
    missing = required - set(enrichment.columns)
    if missing:
        raise RuntimeError(f"Enrichment table missing columns: {sorted(missing)}")

    e = term_categories(enrichment)
    summary = build_category_summary(e)
    gene_sets = load_term_genes([Path(args.go), Path(args.reactome), Path(args.hallmark)])
    leading = leading_genes(enrichment, scores, assignments, gene_sets)
    coverage = category_gene_coverage(leading, assignments, scores)

    e.to_csv(out / "01_category_term_mapping.csv", index=False)
    summary.to_csv(out / "02_category_summary.csv", index=False)
    leading.to_csv(out / "03_leading_genes.csv", index=False)
    coverage.to_csv(out / "04_category_gene_coverage.csv", index=False)

    manifest = {
        "enrichment": str(Path(args.enrichment)),
        "scores": str(Path(args.scores)),
        "assignments": str(Path(args.assignments)),
        "gene_sets": [str(Path(x)) for x in [args.go, args.reactome, args.hallmark]],
        "categories": CATEGORIES,
        "global_fdr_threshold": 0.05,
        "module_fdr_threshold": 0.05,
        "frozen_z6_support_unchanged": True,
    }
    (out / "05_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    write_report(out, summary, coverage, leading, enrichment)

    print("GSE67462 Z6 SPECIFICITY AUDIT")
    print(f"category-matched significant terms: {len(e[(e.fdr_bh_global < 0.05) & (e.fdr_bh_module < 0.05)])}")
    print(f"leading gene-category pairs: {len(leading)}")
    print(f"outputs: {out}")


if __name__ == "__main__":
    main()
