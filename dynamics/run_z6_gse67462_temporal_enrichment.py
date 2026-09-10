"""Run reproducible enrichment for GSE67462 temporal modules.

Diagnostic only. Uses explicit GMT files and the provenance-validated expression
universe as background. Never changes frozen Z6 support decisions.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from dynamics import validation
from dynamics.run_z6_gse67462_identifier_mapping_audit import _build_mapping_report, _norm_symbol, _read_soft_platform
from dynamics.run_z6_gse67462_multimodal_validation_analysis import _load_tss

DEFAULT_ASSIGNMENTS = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_SOFT = "Data/GPL19972_family.soft.gz"
DEFAULT_OUT = "results/Dynamics/z6_gse67462_temporal_enrichment"
DEFAULT_MIN_GENES = 20
DEFAULT_MIN_OVERLAP = 5


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = line.rstrip("\n\r").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = {_norm_symbol(x) for x in parts[2:] if _norm_symbol(x)}
            if name and genes:
                sets[name] = genes
    return sets


def _bh(values: pd.Series) -> pd.Series:
    p = values.to_numpy(dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(p), dtype=float)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return pd.Series(out, index=values.index)


def _enrich(assignments: pd.DataFrame, universe: set[str], gene_sets: dict[str, set[str]], collection: str, min_overlap: int) -> pd.DataFrame:
    rows = []
    for module, block in assignments.groupby("module"):
        module_genes = set(block["gene"].map(_norm_symbol)) & universe
        n = len(module_genes)
        if n == 0:
            continue
        for term, raw_genes in gene_sets.items():
            gs = raw_genes & universe
            if not gs:
                continue
            overlap = len(module_genes & gs)
            if overlap < min_overlap:
                continue
            p = float(hypergeom.sf(overlap - 1, len(universe), len(gs), n))
            rows.append({
                "collection": collection,
                "module": int(module),
                "module_genes": n,
                "term": term,
                "term_genes_in_background": len(gs),
                "overlap": overlap,
                "universe_genes": len(universe),
                "p_value": p,
            })
    columns = ["collection", "module", "module_genes", "term", "term_genes_in_background", "overlap", "universe_genes", "p_value", "fdr_bh_global", "fdr_bh_module"]
    if not rows:
        return pd.DataFrame(columns=columns)
    out = pd.DataFrame(rows).sort_values(["collection", "p_value", "module"]).reset_index(drop=True)
    out["fdr_bh_global"] = out.groupby("collection")["p_value"].transform(_bh)
    out["fdr_bh_module"] = out.groupby(["collection", "module"])["p_value"].transform(_bh)
    return out[columns]


def _write_report(out: Path, enrichment: pd.DataFrame, module_sizes: dict[int, int], manifest: dict) -> None:
    primary_modules = sorted(x for x, n in module_sizes.items() if n >= manifest["min_module_genes"])
    lines = [
        "# GSE67462 Z6 temporal-module enrichment",
        "",
        f"Primary modules (>= {manifest['min_module_genes']} genes): {', '.join('M'+str(x) for x in primary_modules)}.",
        "",
        f"The background is the provenance-validated GSE67462 expression universe (**{manifest['validated_universe_size']:,} genes**). Enrichment is over-representation testing with a one-sided hypergeometric test and Benjamini-Hochberg FDR.",
        "",
    ]
    if enrichment.empty:
        lines.append("No enrichment results were generated. Supply at least one explicit GMT collection.")
    else:
        for collection in sorted(enrichment.collection.unique()):
            lines += [f"## {collection}", "", "| Module | Term | Overlap | p-value | FDR global | FDR module |", "|---:|---|---:|---:|---:|---:|"]
            sub = enrichment[(enrichment.collection == collection) & (enrichment.fdr_bh_global < 0.05)].sort_values(["module", "fdr_bh_global", "p_value"]).head(50)
            if sub.empty:
                lines.append("| — | No term passed FDR < 0.05 | — | — | — | — |")
            else:
                for _, r in sub.iterrows():
                    lines.append(f"| M{int(r.module)} | {r.term} | {int(r.overlap)} | {r.p_value:.3g} | {r.fdr_bh_global:.3g} | {r.fdr_bh_module:.3g} |")
            lines.append("")
    lines += [
        "## Interpretation boundary",
        "",
        "Enriched terms are candidate functional descriptions of temporal modules, not causal pathways. Module labels remain trajectory-derived and should not be replaced by pathway names unless the enrichment evidence supports that interpretation.",
        "",
        "Nuisance/context gene sets, when supplied, are reported separately and are used to qualify biological specificity rather than as an automatic exclusion criterion.",
        "",
        "## Provenance manifest",
        "",
        "```json",
        json.dumps(manifest, indent=2),
        "```",
    ]
    (out / "GSE67462_Z6_temporal_enrichment_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--assignments", default=DEFAULT_ASSIGNMENTS)
    p.add_argument("--gtf", default=DEFAULT_GTF)
    p.add_argument("--platform-soft", default=DEFAULT_SOFT)
    p.add_argument("--go", default=None, help="Versioned GO Biological Process GMT")
    p.add_argument("--reactome", default=None, help="Versioned Reactome GMT")
    p.add_argument("--kegg", default=None, help="Versioned KEGG GMT")
    p.add_argument("--nuisance", default=None, help="Versioned nuisance/context GMT")
    p.add_argument("--min-module-genes", type=int, default=DEFAULT_MIN_GENES)
    p.add_argument("--min-overlap", type=int, default=DEFAULT_MIN_OVERLAP)
    p.add_argument("--output", default=DEFAULT_OUT)
    args = p.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    assignments = pd.read_csv(args.assignments)
    if "gene" not in assignments.columns or "module" not in assignments.columns:
        raise RuntimeError("Assignments file must contain gene and module columns")
    assignments["gene"] = assignments["gene"].map(_norm_symbol)

    matrix, _metadata = validation._load_common_space()
    tss = _load_tss(Path(args.gtf))
    platform = _read_soft_platform(Path(args.platform_soft))
    mapping, mapping_summary = _build_mapping_report(pd.Index(matrix.index.astype(str)), tss, platform)
    validated = set(mapping.loc[mapping["platform_symbol_tss_match"].fillna(0).astype(int) > 0, "expression_id"].map(_norm_symbol))
    if len(validated) < 1000:
        raise RuntimeError(f"Validated background unexpectedly small: {len(validated)} genes")

    module_sizes = {int(k): int(v) for k, v in assignments.groupby("module").size().items()}
    primary = [m for m, n in module_sizes.items() if n >= args.min_module_genes]
    assignments = assignments[assignments.module.isin(primary)].copy()

    sources = {"GO_BP": args.go, "REACTOME": args.reactome, "KEGG": args.kegg, "NUISANCE": args.nuisance}
    results = []
    manifest = {
        "validated_universe_size": len(validated),
        "validated_expression_platform_symbol_to_tss_matches": mapping_summary.get("expression_platform_symbol_to_tss_matches"),
        "gtf": {"path": str(Path(args.gtf)), "sha256": _sha256(Path(args.gtf))},
        "platform_soft": {"path": str(Path(args.platform_soft)), "sha256": _sha256(Path(args.platform_soft))},
        "module_sizes": module_sizes,
        "primary_modules": sorted(primary),
        "min_module_genes": args.min_module_genes,
        "min_overlap": args.min_overlap,
        "fdr_method": "Benjamini-Hochberg",
        "test": "one-sided hypergeometric over-representation",
        "collections": {},
        "frozen_z6_support_unchanged": True,
    }
    for collection, raw_path in sources.items():
        if not raw_path:
            continue
        path = Path(raw_path)
        if not path.exists():
            raise FileNotFoundError(path)
        sets = _read_gmt(path)
        manifest["collections"][collection] = {"path": str(path), "sha256": _sha256(path), "terms": len(sets)}
        results.append(_enrich(assignments, validated, sets, collection, args.min_overlap))

    enrichment = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not enrichment.empty:
        enrichment.to_csv(out / "01_enrichment_all.csv", index=False)
        enrichment.sort_values(["collection", "fdr_bh_global", "p_value"]).to_csv(out / "02_enrichment_ranked.csv", index=False)
    else:
        pd.DataFrame(columns=["collection", "module", "term", "overlap", "p_value", "fdr_bh_global", "fdr_bh_module"]).to_csv(out / "01_enrichment_all.csv", index=False)
        pd.DataFrame().to_csv(out / "02_enrichment_ranked.csv", index=False)
    (out / "03_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    _write_report(out, enrichment, module_sizes, manifest)

    print("GSE67462 TEMPORAL MODULE ENRICHMENT AUDIT")
    print(f"validated background genes: {len(validated)}")
    print(f"primary modules: {sorted(primary)}")
    print(f"collections analyzed: {[x for x in sources if sources[x]]}")
    print(f"significant global-FDR terms: {int((enrichment.fdr_bh_global < 0.05).sum()) if not enrichment.empty else 0}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
