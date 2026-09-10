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

from dynamics.run_z6_gse67462_identifier_mapping_audit import _norm_symbol

DEFAULT_ASSIGNMENTS = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_MAPPING = "results/Dynamics/z6_gse67462_identifier_mapping_audit/02_mapping_report.json"
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
    if not rows:
        return pd.DataFrame(columns=["collection", "module", "module_genes", "term", "term_genes_in_background", "overlap", "universe_genes", "p_value", "fdr_bh_global", "fdr_bh_module"])
    out = pd.DataFrame(rows).sort_values(["collection", "p_value", "module"]).reset_index(drop=True)
    out["fdr_bh_global"] = out.groupby("collection")["p_value"].transform(_bh)
    out["fdr_bh_module"] = out.groupby(["collection", "module"])["p_value"].transform(_bh)
    return out


def _write_report(out: Path, enrichment: pd.DataFrame, assignments: pd.DataFrame, manifest: dict) -> None:
    primary_modules = sorted(int(x) for x in assignments.loc[assignments["n_genes"] >= 20, "module"].unique())
    lines = [
        "# GSE67462 Z6 temporal-module enrichment",
        "",
        f"Primary modules (>=20 genes): {', '.join('M'+str(x) for x in primary_modules)}.",
        "",
        "The background is the provenance-validated GSE67462 expression universe. Enrichment is over-representation testing with a one-sided hypergeometric test and Benjamini-Hochberg FDR.",
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
    p.add_argument("--validated-universe", default=None, help="Optional text file with one validated gene symbol per line; otherwise derive from mapping report")
    p.add_argument("--mapping-report", default=DEFAULT_MAPPING)
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

    if args.validated_universe:
        universe = {_norm_symbol(x.strip()) for x in Path(args.validated_universe).read_text(encoding="utf-8").splitlines() if _norm_symbol(x.strip())}
        universe_source = args.validated_universe
    else:
        report = json.loads(Path(args.mapping_report).read_text(encoding="utf-8"))
        if isinstance(report, dict) and "validated_expression_symbols" in report:
            universe = {_norm_symbol(x) for x in report["validated_expression_symbols"] if _norm_symbol(x)}
        else:
            raise RuntimeError("Mapping report does not contain validated_expression_symbols; supply --validated-universe")
        universe_source = args.mapping_report

    universe &= set(assignments.gene)
    if len(universe) < 1000:
        raise RuntimeError(f"Validated background unexpectedly small: {len(universe)} genes")

    primary = assignments.groupby("module").size()
    assignments = assignments[assignments.module.isin(primary[primary >= args.min_module_genes].index)].copy()

    sources = {"GO_BP": args.go, "REACTOME": args.reactome, "KEGG": args.kegg, "NUISANCE": args.nuisance}
    results = []
    manifest = {
        "validated_universe_size": len(universe),
        "universe_source": universe_source,
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
        results.append(_enrich(assignments, universe, sets, collection, args.min_overlap))

    enrichment = pd.concat(results, ignore_index=True) if results else pd.DataFrame()
    if not enrichment.empty:
        enrichment.to_csv(out / "01_enrichment_all.csv", index=False)
        top = enrichment.sort_values(["collection", "fdr_bh_global", "p_value"])
        top.to_csv(out / "02_enrichment_ranked.csv", index=False)
    else:
        pd.DataFrame(columns=["collection", "module", "term", "overlap", "p_value", "fdr_bh_global", "fdr_bh_module"]).to_csv(out / "01_enrichment_all.csv", index=False)
        pd.DataFrame().to_csv(out / "02_enrichment_ranked.csv", index=False)
    (out / "03_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    _write_report(out, enrichment, assignments, manifest)

    print("GSE67462 TEMPORAL MODULE ENRICHMENT AUDIT")
    print(f"validated background genes: {len(universe)}")
    print(f"primary modules: {sorted(primary[primary >= args.min_module_genes].index.astype(int).tolist())}")
    print(f"collections analyzed: {[x for x in sources if sources[x]]}")
    print(f"significant global-FDR terms: {int((enrichment.fdr_bh_global < 0.05).sum()) if not enrichment.empty else 0}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
