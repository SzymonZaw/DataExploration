"""Biological characterization of frozen Z6 module 6 in Z4.

This is an interpretation audit, not a validation test. It keeps the mouse module
frozen, maps only through the frozen 1:1 orthology table, and uses g:Profiler for
GO Biological Process enrichment. Target-side genes are ranked by Spearman
correlation with D0-D14 gene-activity trajectories in GSE242421; no target-side
module fitting or gene selection is performed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import ssl
import urllib.request
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DEFAULT_MODULES = Path("results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv")
DEFAULT_MAPPING = Path("results/Dynamics/z4_frozen_orthology_mapping/01_frozen_human_mouse_mapping.csv")
DEFAULT_ACTIVITY = Path("results/Dynamics/z4_gse242421_gene_activity/01_gene_activity_log1p_cpm.tsv.gz")
DEFAULT_OUT = Path("results/Dynamics/z4_module6_biological_enrichment")
TIME_ORDER = ["D0", "D2", "D4", "D6", "D8", "D10", "D12", "D14"]


def norm(x: object) -> str:
    return str(x).strip().upper()


def gprofiler(genes: list[str], organism: str, background: list[str] | None, sources: list[str], user_threshold: float = 0.05, insecure_ssl: bool = False) -> pd.DataFrame:
    payload = {"organism": organism, "query": genes, "sources": sources, "user_threshold": user_threshold, "no_evidences": False, "combined": False}
    if background:
        payload["background"] = background
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", data=data,
        headers={"Content-Type": "application/json", "User-Agent": "DataExploration-Z4-Module6/1.1"}, method="POST")
    context = ssl._create_unverified_context() if insecure_ssl else None
    with urllib.request.urlopen(req, timeout=60, context=context) as r:
        obj = json.loads(r.read().decode("utf-8"))
    rows = []
    for r in (obj.get("result") or []):
        rows.append({"source": r.get("source"), "native": r.get("native"), "name": r.get("name"),
                     "p_value": r.get("p_value"), "significant": r.get("significant"),
                     "intersection_size": r.get("intersection_size"), "term_size": r.get("term_size"),
                     "effective_domain_size": r.get("effective_domain_size"),
                     "intersection": ";".join(r.get("intersections") or [])})
    return pd.DataFrame(rows)


def temporal_core(activity: pd.DataFrame, genes: list[str], direction: float, quantile: float = 0.75) -> pd.DataFrame:
    rows = []
    for g in genes:
        if g not in activity.index:
            continue
        y = activity.loc[g, TIME_ORDER].to_numpy(float)
        x = np.arange(len(TIME_ORDER), dtype=float)
        if not np.all(np.isfinite(y)) or np.std(y) == 0:
            continue
        r = spearmanr(x, y).statistic
        rows.append({"gene": g, "spearman_time": float(r), "signed_spearman": float(direction * r)})
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    cutoff = float(df.signed_spearman.quantile(quantile))
    df["core_temporal"] = df.signed_spearman >= cutoff
    df["core_cutoff_quantile"] = quantile
    df["core_cutoff_value"] = cutoff
    return df.sort_values("signed_spearman", ascending=False)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--modules", type=Path, default=DEFAULT_MODULES)
    ap.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    ap.add_argument("--activity", type=Path, default=DEFAULT_ACTIVITY)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--module", type=int, default=6)
    ap.add_argument("--core-quantile", type=float, default=0.75)
    ap.add_argument("--skip-gprofiler", action="store_true")
    ap.add_argument("--insecure-ssl", action="store_true",
                    help="Retry g:Profiler with certificate verification disabled; use only when the local CA bundle cannot validate the service certificate.")
    a = ap.parse_args()
    out = a.output; out.mkdir(parents=True, exist_ok=True)
    for p in [a.modules, a.mapping, a.activity]:
        if not p.exists(): raise FileNotFoundError(f"Required input missing: {p}")

    modules = pd.read_csv(a.modules)
    mapping = pd.read_csv(a.mapping)
    activity = pd.read_csv(a.activity, sep="\t", compression="gzip", index_col=0)
    activity.index = activity.index.map(norm)
    activity = activity[~activity.index.duplicated(keep="first")]

    if not {"module", "gene"}.issubset(modules.columns): raise ValueError("Module table must contain module and gene")
    if not {"human_gene", "mouse_gene"}.issubset(mapping.columns): raise ValueError("Frozen mapping must contain human_gene and mouse_gene")
    if any(c not in activity.columns for c in TIME_ORDER): raise ValueError("Target gene-activity matrix lacks required D0-D14 samples")

    modules["gene_norm"] = modules.gene.map(norm)
    m6_mouse = sorted(modules.loc[pd.to_numeric(modules.module, errors="coerce") == a.module, "gene_norm"].unique())
    mapping["mouse_norm"] = mapping.mouse_gene.map(norm)
    mapping["human_norm"] = mapping.human_gene.map(norm)
    frozen = mapping.dropna(subset=["mouse_norm", "human_norm"]).drop_duplicates("mouse_norm")
    mouse_to_human = dict(zip(frozen.mouse_norm, frozen.human_norm))
    m6_human = sorted({mouse_to_human[g] for g in m6_mouse if g in mouse_to_human and mouse_to_human[g] in activity.index})
    core = temporal_core(activity, m6_human, 1.0, a.core_quantile)
    core_genes = core.loc[core.core_temporal, "gene"].tolist() if not core.empty else []
    mouse_background = sorted(frozen.mouse_norm.unique())
    human_background = sorted(set(frozen.human_norm) & set(activity.index))

    summary = {"candidate": "GSE242421", "module": a.module, "frozen_mouse_genes": len(m6_mouse),
               "validated_human_genes": len(m6_human), "target_samples": TIME_ORDER,
               "target_module_fitting": False, "target_gene_selection": False,
               "core_definition": f"top {(1-a.core_quantile)*100:.0f}% by signed D0-D14 Spearman within frozen validated module genes",
               "interpretation_only": True, "z4_decision_unchanged": "Z4_ORTHO_THRESHOLD_SENSITIVE",
               "gprofiler": not a.skip_gprofiler, "insecure_ssl_requested": a.insecure_ssl}

    core.to_csv(out / "01_module6_target_gene_ranking.csv", index=False)
    pd.DataFrame({"mouse_gene": m6_mouse, "human_gene": [mouse_to_human.get(g, "") for g in m6_mouse]}).to_csv(out / "02_module6_frozen_orthology.csv", index=False)
    pd.DataFrame({"gene": core_genes}).to_csv(out / "03_module6_core_temporal_genes.csv", index=False)

    if not a.skip_gprofiler:
        results = {}
        for label, genes, organism, bg in [("mouse_frozen_module", m6_mouse, "mmusculus", mouse_background),
                                           ("human_validated_module", m6_human, "hsapiens", human_background),
                                           ("human_target_core", core_genes, "hsapiens", human_background)]:
            if len(genes) >= 3:
                try:
                    df = gprofiler(genes, organism, bg, ["GO:BP", "REAC"], insecure_ssl=a.insecure_ssl)
                    df.to_csv(out / f"04_{label}_enrichment.csv", index=False)
                    results[label] = {"status": "OK", "n_genes": len(genes), "n_terms": len(df)}
                except Exception as e:
                    results[label] = {"status": "ERROR", "n_genes": len(genes), "error": str(e)}
            else:
                results[label] = {"status": "SKIPPED", "n_genes": len(genes)}
        summary["enrichment_results"] = results
    (out / "05_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("GSE242421 Z4 MODULE 6 BIOLOGICAL ENRICHMENT AUDIT")
    print(f"frozen mouse genes: {len(m6_mouse)}")
    print(f"validated human genes: {len(m6_human)}")
    print(f"target core temporal genes: {len(core_genes)}")
    print("decision unchanged: Z4_ORTHO_THRESHOLD_SENSITIVE")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
