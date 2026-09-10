"""Biological characterization of frozen Z6 module 6 in Z4.

Interpretation-only audit. The mouse module and 1:1 orthology are frozen.
"""
from __future__ import annotations
import argparse, json, ssl, urllib.error, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DEFAULT_MODULES = Path("results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv")
DEFAULT_MAPPING = Path("results/Dynamics/z4_frozen_orthology_mapping/01_frozen_human_mouse_mapping.csv")
DEFAULT_ACTIVITY = Path("results/Dynamics/z4_gse242421_gene_activity/01_gene_activity_log1p_cpm.tsv.gz")
DEFAULT_OUT = Path("results/Dynamics/z4_module6_biological_enrichment")
TIME_ORDER = ["D0","D2","D4","D6","D8","D10","D12","D14"]

def norm(x: object) -> str: return str(x).strip().upper()

def flatten_values(x):
    if x is None: return []
    if isinstance(x, (list, tuple)):
        out=[]
        for v in x: out.extend(flatten_values(v))
        return out
    if isinstance(x, str): return [v for v in x.split(";") if v]
    return [str(x)]

def _query_mapping_from_meta(obj, query_name, fallback_query):
    """Recover original submitted IDs aligned to g:Profiler intersections.

    The API's intersections are aligned to the mapped Ensembl IDs in
    meta.genes_metadata.query[query_name].ensgs, not necessarily to the raw
    input list: unmapped/duplicate identifiers can change query_size.
    The corresponding `mapping` field maps those Ensembl IDs back to input IDs.
    """
    try:
        qmeta=obj["meta"]["genes_metadata"]["query"][query_name]
        mapping=qmeta.get("mapping")
        ensgs=qmeta.get("ensgs")
        if isinstance(mapping, list) and len(mapping)==len(ensgs):
            return [str(x) for x in mapping]
        # Some API versions expose mapping as a dict/list-of-lists. Flatten only
        # enough to recover a one-to-one aligned identifier list.
        if isinstance(mapping, dict):
            vals=[]
            for x in ensgs:
                v=mapping.get(x, x)
                vals.append(str(v[0] if isinstance(v,list) and v else v))
            return vals
    except Exception:
        pass
    return [str(x) for x in fallback_query]

def gprofiler(genes, organism, background, sources, insecure_ssl=False):
    query_genes=[str(g) for g in genes]
    payload={"organism":organism,"query":query_genes,"sources":sources,"user_threshold":0.05,"no_evidences":False,"combined":False}
    if background:
        payload["background"]=[str(g) for g in background]; payload["domain_scope"]="custom"
    data=json.dumps(payload).encode("utf-8")
    req=urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",data=data,headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"DataExploration-Z4-Module6/1.7"},method="POST")
    ctx=ssl._create_unverified_context() if insecure_ssl else None
    try:
        with urllib.request.urlopen(req,timeout=60,context=ctx) as r: obj=json.loads(r.read().decode("utf-8"))
    except urllib.error.HTTPError as e:
        body=e.read().decode("utf-8",errors="replace"); raise RuntimeError(f"g:Profiler HTTP {e.code}: {body}") from e
    rows=[]
    for r in (obj.get("result") or []):
        intersections=r.get("intersections") or []
        query_name=r.get("query") or "query_1"
        mapped_query=_query_mapping_from_meta(obj,query_name,query_genes)
        if len(intersections)==len(mapped_query):
            genes_hit=[mapped_query[i] for i,item in enumerate(intersections) if item not in (None, [], "")]
        else:
            genes_hit=[]
        rows.append({"source":r.get("source"),"native":r.get("native"),"name":r.get("name"),"p_value":r.get("p_value"),"significant":r.get("significant"),"intersection_size":r.get("intersection_size"),"term_size":r.get("term_size"),"effective_domain_size":r.get("effective_domain_size"),"intersection":";".join(genes_hit)})
    return pd.DataFrame(rows)

def temporal_core(activity, genes, quantile=0.75):
    rows=[]
    for g in genes:
        if g not in activity.index: continue
        y=activity.loc[g,TIME_ORDER].to_numpy(float)
        if not np.all(np.isfinite(y)) or np.std(y)==0: continue
        r=spearmanr(np.arange(len(TIME_ORDER),dtype=float),y).statistic
        rows.append({"gene":g,"spearman_time":float(r),"signed_spearman":float(r)})
    df=pd.DataFrame(rows)
    if df.empty: return df
    cutoff=float(df.signed_spearman.quantile(quantile)); df["core_temporal"]=df.signed_spearman>=cutoff; df["core_cutoff_quantile"]=quantile; df["core_cutoff_value"]=cutoff
    return df.sort_values("signed_spearman",ascending=False)

def read_enrichment(path):
    if not path.exists() or path.stat().st_size==0: return pd.DataFrame()
    try: return pd.read_csv(path)
    except pd.errors.EmptyDataError: return pd.DataFrame()

def build_term_gene_audit(enrichment_files,m6_mouse,m6_human,core_genes,mouse_to_human,out):
    m6_mouse_set=set(m6_mouse); m6_human_set=set(m6_human); core_set=set(core_genes); rows=[]
    for label,path in enrichment_files.items():
        df=read_enrichment(path)
        for _,r in df.iterrows():
            hits=set(norm(x) for x in flatten_values(r.get("intersection")) if str(x).strip())
            if label=="mouse_frozen_module":
                module_hits=hits&m6_mouse_set
                human_hits={mouse_to_human[g] for g in module_hits if g in mouse_to_human and mouse_to_human[g] in m6_human_set}
                core_hits={mouse_to_human[g] for g in module_hits if g in mouse_to_human and mouse_to_human[g] in core_set}
            else:
                module_hits=hits&m6_human_set; human_hits=module_hits; core_hits=module_hits&core_set
            rows.append({"analysis":label,"source":r.get("source"),"native":r.get("native"),"name":r.get("name"),"p_value":r.get("p_value"),"significant":r.get("significant"),"reported_intersection_size":r.get("intersection_size"),"intersection_genes":";".join(sorted(hits)),"module6_gene_hits":";".join(sorted(module_hits)),"module6_gene_hit_count":len(module_hits),"validated_human_hits":";".join(sorted(human_hits)),"validated_human_hit_count":len(human_hits),"target_core_hits":";".join(sorted(core_hits)),"target_core_hit_count":len(core_hits)})
    columns=["analysis","source","native","name","p_value","significant","reported_intersection_size","intersection_genes","module6_gene_hits","module6_gene_hit_count","validated_human_hits","validated_human_hit_count","target_core_hits","target_core_hit_count"]
    audit=pd.DataFrame(rows,columns=columns); audit.to_csv(out/"06_module6_term_gene_audit.csv",index=False); return audit

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--modules",type=Path,default=DEFAULT_MODULES); ap.add_argument("--mapping",type=Path,default=DEFAULT_MAPPING); ap.add_argument("--activity",type=Path,default=DEFAULT_ACTIVITY); ap.add_argument("--output",type=Path,default=DEFAULT_OUT); ap.add_argument("--module",type=int,default=6); ap.add_argument("--core-quantile",type=float,default=0.75); ap.add_argument("--skip-gprofiler",action="store_true"); ap.add_argument("--insecure-ssl",action="store_true"); a=ap.parse_args(); out=a.output; out.mkdir(parents=True,exist_ok=True)
    for p in [a.modules,a.mapping,a.activity]:
        if not p.exists(): raise FileNotFoundError(f"Required input missing: {p}")
    modules=pd.read_csv(a.modules); mapping=pd.read_csv(a.mapping); activity=pd.read_csv(a.activity,sep="\t",compression="gzip",index_col=0); activity.index=activity.index.map(norm); activity=activity[~activity.index.duplicated()]; modules["gene_norm"]=modules.gene.map(norm)
    m6_mouse=sorted(modules.loc[pd.to_numeric(modules.module,errors="coerce")==a.module,"gene_norm"].unique()); mapping["mouse_norm"]=mapping.mouse_gene.map(norm); mapping["human_norm"]=mapping.human_gene.map(norm); frozen=mapping.dropna(subset=["mouse_norm","human_norm"]).drop_duplicates("mouse_norm"); mouse_to_human=dict(zip(frozen.mouse_norm,frozen.human_norm)); m6_human=sorted({mouse_to_human[g] for g in m6_mouse if g in mouse_to_human and mouse_to_human[g] in activity.index}); core=temporal_core(activity,m6_human,a.core_quantile); core_genes=core.loc[core.core_temporal,"gene"].tolist() if not core.empty else []; mouse_bg=sorted(frozen.mouse_norm.unique()); human_bg=sorted(set(frozen.human_norm)&set(activity.index))
    summary={"candidate":"GSE242421","module":a.module,"frozen_mouse_genes":len(m6_mouse),"validated_human_genes":len(m6_human),"target_core_temporal_genes":len(core_genes),"target_samples":TIME_ORDER,"target_module_fitting":False,"target_gene_selection":False,"core_definition":f"top {(1-a.core_quantile)*100:.0f}% by signed D0-D14 Spearman within frozen validated module genes","interpretation_only":True,"z4_decision_unchanged":"Z4_ORTHO_THRESHOLD_SENSITIVE","gprofiler":not a.skip_gprofiler,"insecure_ssl_requested":a.insecure_ssl}
    core.to_csv(out/"01_module6_target_gene_ranking.csv",index=False); pd.DataFrame({"mouse_gene":m6_mouse,"human_gene":[mouse_to_human.get(g,"") for g in m6_mouse]}).to_csv(out/"02_module6_frozen_orthology.csv",index=False); pd.DataFrame({"gene":core_genes}).to_csv(out/"03_module6_core_temporal_genes.csv",index=False)
    enrichment_files={"mouse_frozen_module":out/"04_mouse_frozen_module_enrichment.csv","human_validated_module":out/"04_human_validated_module_enrichment.csv","human_target_core":out/"04_human_target_core_enrichment.csv"}
    if not a.skip_gprofiler:
        results={}
        for label,genes,org,bg in [("mouse_frozen_module",m6_mouse,"mmusculus",mouse_bg),("human_validated_module",m6_human,"hsapiens",human_bg),("human_target_core",core_genes,"hsapiens",human_bg)]:
            try:
                df=gprofiler(genes,org,bg,["GO:BP","REAC"],a.insecure_ssl) if len(genes)>=3 else pd.DataFrame(columns=["source","native","name","p_value","significant","intersection_size","term_size","effective_domain_size","intersection"]); df.to_csv(enrichment_files[label],index=False); results[label]={"status":"OK","n_genes":len(genes),"n_terms":len(df)}
            except Exception as e: results[label]={"status":"ERROR","n_genes":len(genes),"error":str(e)}
        summary["enrichment_results"]=results
    audit=build_term_gene_audit(enrichment_files,m6_mouse,m6_human,core_genes,mouse_to_human,out); summary["term_gene_audit"]={"status":"OK","n_rows":len(audit),"terms_with_module6_hits":int((audit.module6_gene_hit_count>0).sum()) if not audit.empty else 0,"terms_with_target_core_hits":int((audit.target_core_hit_count>0).sum()) if not audit.empty else 0}; (out/"05_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print("GSE242421 Z4 MODULE 6 BIOLOGICAL ENRICHMENT AUDIT"); print(f"frozen mouse genes: {len(m6_mouse)}"); print(f"validated human genes: {len(m6_human)}"); print(f"target core temporal genes: {len(core_genes)}"); print("decision unchanged: Z4_ORTHO_THRESHOLD_SENSITIVE"); print(f"output: {out}")
if __name__=="__main__": main()
