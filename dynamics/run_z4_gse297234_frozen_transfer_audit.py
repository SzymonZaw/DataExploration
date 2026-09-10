"""Frozen Z4 transfer audit for GSE297234.

This audit transfers the frozen GSE67462 temporal-module gene sets into the
independent human GSE297234 scRNA-seq experiment after human->mouse mapping.
It deliberately treats time/module concordance as a transfer diagnostic only:
GSE297234 does not currently provide a clearly independent biological endpoint
that can satisfy the full Z4-A/B/C acceptance rule.

The Seurat object is read by the companion R script so that the Python audit
never depends on a Python Seurat reader.
"""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DEFAULT_RDS = "Data/GSE297234_HFIB_COMBINED_SEVOSKM.rds"
DEFAULT_MAPPING = "results/Dynamics/z4_frozen_orthology_mapping/01_frozen_human_mouse_mapping.csv"
DEFAULT_MODULES = "results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv"
DEFAULT_METADATA = "results/Dynamics/z4_external_metadata_audit/GSE297234/metadata.csv"
DEFAULT_OUT = "results/Dynamics/z4_gse297234_frozen_transfer_audit"
R_SCRIPT = Path(__file__).with_name("run_z4_gse297234_frozen_transfer_audit.R")


def _spearman(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3 or np.std(a[mask]) == 0 or np.std(b[mask]) == 0:
        return np.nan
    return float(spearmanr(a[mask], b[mask]).statistic)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rds", default=DEFAULT_RDS)
    ap.add_argument("--mapping", default=DEFAULT_MAPPING)
    ap.add_argument("--modules", default=DEFAULT_MODULES)
    ap.add_argument("--metadata", default=DEFAULT_METADATA)
    ap.add_argument("--output", default=DEFAULT_OUT)
    args = ap.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)

    pseudobulk = out / "01_gse297234_pseudobulk_human.csv"
    meta_out = out / "02_gse297234_sample_metadata.csv"
    subprocess.run(["Rscript", str(R_SCRIPT), "--rds", args.rds,
                    "--output-expression", str(pseudobulk),
                    "--output-metadata", str(meta_out)], check=True)

    expr = pd.read_csv(pseudobulk, index_col=0)
    meta = pd.read_csv(meta_out)
    mapping = pd.read_csv(args.mapping)
    modules = pd.read_csv(args.modules)

    mapping = mapping.dropna(subset=["human_gene", "mouse_gene"]).drop_duplicates()
    mapping = mapping.drop_duplicates("human_gene")
    human_to_mouse = dict(zip(mapping["human_gene"].astype(str), mapping["mouse_gene"].astype(str)))

    # The expression matrix is human; transfer each frozen mouse module back to
    # its human orthologues. No module discovery is performed on validation data.
    module_results = []
    transferred_genes = []
    for module, g in modules.groupby("module"):
        mouse_genes = set(g["gene"].astype(str))
        human_genes = [h for h, m in human_to_mouse.items() if m in mouse_genes and h in expr.index]
        transferred_genes.extend(human_genes)
        if not human_genes:
            continue
        values = expr.loc[sorted(set(human_genes))]
        score = values.mean(axis=0)
        for sample, value in score.items():
            module_results.append({"module": int(module), "sample": sample,
                                   "score": float(value), "n_genes": len(set(human_genes))})

    scores = pd.DataFrame(module_results)
    scores.to_csv(out / "03_frozen_module_scores.csv", index=False)

    # Merge with explicit GEO metadata where sample IDs can be resolved.
    if Path(args.metadata).exists():
        geo = pd.read_csv(args.metadata)
        geo_cols = [c for c in ["sample_id", "title", "source", "platform", "characteristics", "inferred_time_days"] if c in geo.columns]
        geo = geo[geo_cols].copy()
        scores = scores.merge(geo, left_on="sample", right_on="sample_id", how="left")
    else:
        scores["inferred_time_days"] = np.nan

    scores.to_csv(out / "04_frozen_module_scores_with_metadata.csv", index=False)

    summary_rows = []
    for module, g in scores.groupby("module"):
        pivot = g.pivot_table(index="inferred_time_days", values="score", aggfunc="mean").dropna()
        rho = _spearman(pivot.index.to_numpy(), pivot["score"].to_numpy()) if len(pivot) >= 3 else np.nan
        summary_rows.append({"module": int(module), "n_genes": int(g["n_genes"].max()),
                             "n_samples": int(g["sample"].nunique()),
                             "n_timepoints": int(len(pivot)),
                             "time_spearman": rho})
    summary = pd.DataFrame(summary_rows).sort_values("module")
    summary.to_csv(out / "05_transfer_summary.csv", index=False)

    result = {
        "candidate": "GSE297234",
        "frozen_mapping": str(Path(args.mapping)),
        "frozen_modules": str(Path(args.modules)),
        "validated_transfer_genes": int(len(set(transferred_genes))),
        "modules_evaluated": int(summary.shape[0]),
        "endpoint_independence": "NOT_ESTABLISHED",
        "z4_A_external_target_association": "NOT_EVALUATED",
        "z4_B_target_vs_nontarget_specificity": "NOT_EVALUATED",
        "z4_C_context_robustness": "NOT_EVALUATED",
        "z4_D_replication": "DIAGNOSTIC_ONLY",
        "interpretation": "EXTERNAL_TRANSFER_DIAGNOSTIC_ONLY",
        "decision": "Z4_UNRESOLVED",
        "reason": "GSE297234 provides an independent OSKM time series, but the currently available metadata/Seurat objects do not establish an independent biological endpoint suitable for the frozen Z4 acceptance rule.",
        "frozen_z4_rule_unchanged": True,
    }
    (out / "06_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("GSE297234 Z4 FROZEN TRANSFER AUDIT")
    print(f"validated transfer genes: {result['validated_transfer_genes']}")
    print(f"modules evaluated: {result['modules_evaluated']}")
    print("endpoint independence: NOT ESTABLISHED")
    print("decision: Z4_UNRESOLVED")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
