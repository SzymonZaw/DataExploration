"""Cluster the GSE67462 multimodal mechanistic core into temporal modules.

Diagnostic only. Does not modify frozen Z6 predictive or multimodal support.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

from dynamics import validation
from dynamics.run_z6_gse67462_mechanistic_interpretation import _normalized_expression, _numeric_time
from dynamics.run_z6_gse67462_identifier_mapping_audit import _build_mapping_report, _norm_symbol, _read_soft_platform
from dynamics.run_z6_gse67462_multimodal_validation_analysis import _build_modality_matrix, _discover_files, _expression_common, _load_tss, _spearman

DEFAULT_OUT = "results/Dynamics/z6_gse67462_temporal_modules"
DEFAULT_MECH = "results/Dynamics/z6_gse67462_mechanistic_interpretation/01_gene_multimodal_scores.csv"
DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
DEFAULT_SOFT = "Data/GPL19972_family.soft.gz"
ACTIVE = ["h3k27ac", "h3k4me3", "rnapii", "total_oct4"]
MODALITY_LABELS = {"h3k27ac":"H3K27ac", "h3k4me3":"H3K4me3", "rnapii":"RNAPII", "total_oct4":"OCT4", "h3k27me3":"H3K27me3"}


def _frames(values, genes):
    out = {}
    for modality in sorted(set(m for m, _, _ in values)):
        entries = sorted([(float(t), v) for m, t, v in values if m == modality and _numeric_time(t) is not None], key=lambda x:x[0])
        if entries:
            out[modality] = pd.DataFrame([{g:v.get(g,0.0) for g in genes} for _,v in entries], index=[t for t,_ in entries], columns=genes)
    return out


def _replicate_mean(expr_raw, metadata):
    times = []
    for x in metadata.get("time_hours", []):
        t = _numeric_time(x)
        times.append(t)
    expr = expr_raw.copy()
    expr.index = times
    expr = expr.loc[[x is not None for x in expr.index]]
    grouped = expr.groupby(expr.index).mean()
    return grouped.sort_index()


def _zscore_rows(df):
    a = df.to_numpy(dtype=float)
    mean = np.nanmean(a, axis=1, keepdims=True)
    sd = np.nanstd(a, axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    return pd.DataFrame((a-mean)/sd, index=df.index, columns=df.columns)


def _trajectory_label(c):
    x=np.asarray(c,float)
    if len(x)<4:return "mixed"
    early=np.mean(x[:3]); late=np.mean(x[-3:]); peak=int(np.argmax(x)); trough=int(np.argmin(x))
    if late-early > 0.8 and peak >= len(x)-2:return "late_rising"
    if early-late > 0.8 and trough <= 2:return "early_declining"
    if peak <= 2 and late < peak-0.4:return "early_transient"
    if peak >= len(x)-3 and early < late-0.4:return "late_rising"
    if abs(late-early)<0.4 and np.max(x)-np.min(x)>1.0:return "transient"
    if np.mean(np.abs(np.diff(x)))<0.25:return "sustained"
    return "mixed"


def _module_table(labels, z, expr, frames):
    rows=[]
    for module in sorted(np.unique(labels)):
        genes=z.index[labels==module]
        centroid=z.loc[genes].mean(axis=0).to_numpy()
        within=[]
        for g in genes:
            r=_spearman(z.loc[g].to_numpy(), centroid)
            if np.isfinite(r):within.append(r)
        modality_scores={}
        for m in frames:
            f=frames[m].reindex(index=expr.index).interpolate(limit_direction="both")
            vals=[]
            for g in genes:
                if g in f.columns:
                    r=_spearman(centroid,f[g].to_numpy())
                    if np.isfinite(r): vals.append(r)
            modality_scores[m]=float(np.median(vals)) if vals else np.nan
        rows.append({"module":int(module),"n_genes":len(genes),"trajectory_label":_trajectory_label(centroid),"within_module_median_rho":float(np.median(within)) if within else np.nan,**{f"module_rho_{m}":v for m,v in modality_scores.items()},**{f"centroid_t{i}":float(v) for i,v in enumerate(centroid)}})
    return pd.DataFrame(rows).sort_values("module")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--mechanistic-scores",default=DEFAULT_MECH)
    p.add_argument("--gtf",default=DEFAULT_GTF)
    p.add_argument("--platform-soft",default=DEFAULT_SOFT)
    p.add_argument("--modality-root",default="Data/GSE67520")
    p.add_argument("--output",default=DEFAULT_OUT)
    p.add_argument("--k-min",type=int,default=4)
    p.add_argument("--k-max",type=int,default=8)
    args=p.parse_args(); out=Path(args.output); out.mkdir(parents=True,exist_ok=True)

    matrix, metadata=validation._load_common_space()
    tss=_load_tss(Path(args.gtf)); platform=_read_soft_platform(Path(args.platform_soft))
    mapping,_=_build_mapping_report(pd.Index(matrix.index.astype(str)),tss,platform)
    validated=set(mapping.loc[mapping["platform_symbol_tss_match"].fillna(0).astype(int)>0,"expression_id"].map(_norm_symbol))
    _, expr_raw, _=_expression_common(matrix,metadata)
    expr=_normalized_expression(expr_raw,validated)
    expr=_replicate_mean(expr, {"time_hours":expr.index.tolist()})
    scores=pd.read_csv(args.mechanistic_scores)
    core=set(scores.loc[scores["mechanistic_class"]=="MULTIMODAL_CORE","gene"].astype(str))
    core=core & set(expr.columns.astype(str)); expr=expr.loc[:,sorted(core)]
    z=_zscore_rows(expr.T)
    corr=z.T.corr(method="spearman").fillna(0.0); dist=1.0-corr.to_numpy(); np.fill_diagonal(dist,0.0)
    link=linkage(squareform(np.clip(dist,0,2),checks=False),method="average")
    records=[]; assignments={}
    for k in range(args.k_min,args.k_max+1):
        lab=fcluster(link,k,criterion="maxclust"); assignments[k]=lab
        for k2 in range(args.k_min,args.k_max+1):
            if k2>k:
                lab2=fcluster(link,k2,criterion="maxclust"); records.append({"k_a":k,"k_b":k2,"adjusted_rand_index":float(adjusted_rand_score(lab,lab2))})
    chosen=int(round((args.k_min+args.k_max)/2)); labels=assignments[chosen]
    files=_discover_files(Path(args.modality_root)); values_raw,_=_build_modality_matrix(files,tss); values=[]
    for m,t,v in values_raw:
        tt=_numeric_time(t)
        if tt is None:continue
        values.append((m,tt,{_norm_symbol(g):float(x) for g,x in v.items() if _norm_symbol(g) in core}))
    frames=_frames(values,expr.columns)
    modules=_module_table(labels,z,expr,frames)
    gene_module=pd.DataFrame({"gene":z.index,"module":labels}).sort_values(["module","gene"])
    centroids=z.copy(); centroids["module"]=labels; centroids=centroids.groupby("module").mean()
    pd.DataFrame(records).to_csv(out/"01_cluster_stability.csv",index=False)
    modules.to_csv(out/"02_temporal_modules.csv",index=False)
    gene_module.to_csv(out/"03_gene_module_assignments.csv",index=False)
    centroids.to_csv(out/"04_module_centroids.csv")
    summary={"raw_common_space_genes":int(expr_raw.shape[1]),"validated_genes":int(len(validated)),"multimodal_core_genes":int(len(core)),"timepoints":expr.index.tolist(),"tested_k":list(range(args.k_min,args.k_max+1)),"chosen_k":chosen,"module_count":int(len(modules)),"median_pairwise_ari":float(np.median([r["adjusted_rand_index"] for r in records])) if records else None,"interpretation":"TEMPORAL_MODULES_FOR_MECHANISTIC_INTERPRETATION","frozen_z6_support_unchanged":True}
    (out/"05_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    lines=["# GSE67462 Z6 temporal modules","",f"The multimodal core contains **{len(core):,} genes** after provenance validation. Hierarchical clustering was evaluated for k={args.k_min}..{args.k_max}; the descriptive report uses k={chosen}.","","| Module | Genes | Shape | Within-module rho | H3K27ac | H3K4me3 | RNAPII | OCT4 |","|---:|---:|---|---:|---:|---:|---:|---:|"]
    for _,r in modules.iterrows():lines.append(f"| {int(r.module)} | {int(r.n_genes)} | {r.trajectory_label} | {r.within_module_median_rho:.3f} | {r.get('module_rho_h3k27ac',np.nan):.3f} | {r.get('module_rho_h3k4me3',np.nan):.3f} | {r.get('module_rho_rnapii',np.nan):.3f} | {r.get('module_rho_total_oct4',np.nan):.3f} |")
    lines += ["","## Interpretation","","These modules describe recurring temporal expression shapes within the multimodal core. They are not causal pathways. Functional enrichment should be applied only with an explicitly versioned gene-set collection.",""]
    (out/"GSE67462_Z6_temporal_modules_report.md").write_text("\n".join(lines),encoding="utf-8")
    print("GSE67462 TEMPORAL MODULE AUDIT")
    print(f"multimodal core genes: {len(core)}")
    print(f"chosen k: {chosen}")
    print(f"modules: {len(modules)}")
    print(f"median pairwise ARI: {summary['median_pairwise_ari']}")
    print(f"Outputs: {out}")

if __name__=="__main__":main()
