"""Stage 2.9.12: sensitivity analysis of training-gene selection.

This diagnostic asks whether the Stage 2.9.10 enrichment failure for some
held-out datasets is caused by the particular 2,500-gene stability cutoff.
It does not relax thresholds in the validation stage and does not fit an
ODE/state-space model.

For each valid Stage 2.9.10 fold, enrichment is evaluated for several nested
training-gene sets: top 1,000 / 2,500 / 5,000 stable genes and all genes that
pass the Stage 2.9.1 stability-selection rule.  The cached 2,500-gene result
is reused where possible; other sizes are queried once and cached.
"""
from pathlib import Path
import json, os, ssl, time, urllib.error, urllib.request
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
COMMON=ROOT/"results"/"Dynamics"/"stage2_6"
IN=ROOT/"results"/"Dynamics"/"stage2_9_1"
OUT=ROOT/"results"/"Dynamics"/"stage2_9_12"
CACHE=OUT/"cache"
OUT.mkdir(parents=True,exist_ok=True); CACHE.mkdir(parents=True,exist_ok=True)

def log(x): print(f"Stage 2.9.12: {x}",flush=True)

def sslctx(verify=True):
    if not verify:return ssl._create_unverified_context()
    try:
        import certifi; return ssl.create_default_context(cafile=certifi.where())
    except ImportError:return ssl.create_default_context()

def request(payload,label):
    path=CACHE/f"{label}.json"
    if path.exists():
        try:return json.loads(path.read_text(encoding="utf-8"))
        except Exception:pass
    body=json.dumps(payload).encode(); insecure=os.environ.get("STAGE298_INSECURE_SSL","")=="1"
    for attempt in range(3):
        try:
            req=urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",data=body,headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"DataExploration-stage2.9.12/1.0"},method="POST")
            with urllib.request.urlopen(req,timeout=90,context=sslctx(True)) as r:data=json.loads(r.read().decode())
            path.write_text(json.dumps(data),encoding="utf-8"); return data
        except urllib.error.URLError as e:
            if insecure and attempt==0 and isinstance(e.reason,ssl.SSLCertVerificationError):
                log(f"{label}: verified TLS failed; retrying with verification disabled")
                try:
                    req=urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",data=body,headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"DataExploration-stage2.9.12/1.0"},method="POST")
                    with urllib.request.urlopen(req,timeout=90,context=sslctx(False)) as r:data=json.loads(r.read().decode())
                    path.write_text(json.dumps(data),encoding="utf-8"); return data
                except Exception:pass
            if attempt<2:time.sleep(2**attempt)
        except Exception:
            if attempt<2:time.sleep(2**attempt)
    log(f"{label}: enrichment request failed"); return None

def load_space():
    from dynamics.validation import _load_common_space
    return _load_common_space()

def load_fold_genes():
    stab=pd.read_csv(IN/"03_gene_stability.csv"); status=pd.read_csv(IN/"06_fold_status.csv")
    genes=pd.read_csv(COMMON/"06_common_human_gene_matrix.csv",index_col=0).index.astype(str).tolist(); out={}
    for _,r in status.iterrows():
        ds=str(r.held_out_dataset)
        if str(r.get("reason",""))!="ok":continue
        s=stab[stab.held_out_dataset.astype(str)==ds].copy()
        s["stability"]=pd.to_numeric(s.median_pairwise_correlation,errors="coerce")
        s=s[s.stability.notna()]
        q=s[s.n_datasets>=3].copy()
        if len(q)<20:q=s[s.n_datasets>=2].copy()
        if len(q)<20:continue
        threshold=float(q.stability.quantile(.75)); q=q[q.stability>=threshold]
        q=q.sort_values(["stability","n_datasets"],ascending=False)
        out[ds]=[genes[i] for i in q.gene_index.astype(int) if 0<=i<len(genes)]
    return out

def enrich(genes,background,label):
    payload={"organism":"hsapiens","query":genes,"sources":["GO:BP","REAC","KEGG"],"domain_scope":"custom","background":background,"user_threshold":0.05,"significance_threshold_method":"fdr","no_evidences":False,"output":"json"}
    log(f"{label}: g:Profiler on {len(genes):,} training genes")
    return request(payload,label)

def term_rows(data,held,size):
    rows=[]
    if not data:return rows
    for r in data.get("result",[]):
        p=pd.to_numeric(r.get("p_value"),errors="coerce")
        if not np.isfinite(p):continue
        src=str(r.get("source",r.get("source_id","")))
        rows.append({"held_out_dataset":held,"gene_set_size":size,"source":src,"term_id":r.get("native"),"term_name":r.get("name"),"p_value":float(p),"intersection_size":int(r.get("intersection_size",0) or 0)})
    return rows

def run():
    log("starting gene-selection sensitivity diagnostic; no ODE/state-space model")
    matrix,meta=load_space(); background=matrix.index.astype(str).tolist(); folds=load_fold_genes()
    log(f"loaded {len(background):,} background genes and {len(folds)} training folds")
    sizes=[1000,2500,5000]
    all_terms=[]; diagnostics=[]
    for held,all_genes in sorted(folds.items()):
        log(f"FOLD {held}: {len(all_genes):,} genes available after Stage 2.9.1 selection")
        for size in sizes:
            query=all_genes[:min(size,len(all_genes))]
            label=f"gprofiler_{held}_{len(query)}"
            data=enrich(query,background,label)
            terms=term_rows(data,held,len(query)); all_terms.extend(terms)
            fdr=[x for x in terms if x["p_value"]<.05]
            by={s:sum(1 for x in fdr if x["source"]==s) for s in ("GO:BP","REAC","KEGG")}
            diagnostics.append({"held_out_dataset":held,"gene_set_size":len(query),"available_selected_genes":len(all_genes),"cache_or_result":bool(data),"n_result_terms":len(terms),"n_fdr05_terms":len(fdr),"n_go_bp_fdr05":by["GO:BP"],"n_reactome_fdr05":by["REAC"],"n_kegg_fdr05":by["KEGG"]})
        # All selected genes is deliberately tested separately: it is a diagnostic, not a new validation cutoff.
        query=all_genes; label=f"gprofiler_{held}_all"
        data=enrich(query,background,label); terms=term_rows(data,held,len(query)); all_terms.extend(terms)
        fdr=[x for x in terms if x["p_value"]<.05]; by={s:sum(1 for x in fdr if x["source"]==s) for s in ("GO:BP","REAC","KEGG")}
        diagnostics.append({"held_out_dataset":held,"gene_set_size":len(query),"available_selected_genes":len(all_genes),"cache_or_result":bool(data),"n_result_terms":len(terms),"n_fdr05_terms":len(fdr),"n_go_bp_fdr05":by["GO:BP"],"n_reactome_fdr05":by["REAC"],"n_kegg_fdr05":by["KEGG"]})
    diag=pd.DataFrame(diagnostics); terms=pd.DataFrame(all_terms)
    diag.to_csv(OUT/"01_gene_set_sensitivity.csv",index=False); terms.to_csv(OUT/"02_enriched_terms_by_gene_set.csv",index=False)
    if len(terms):
        sig=terms[terms.p_value<.05].copy(); cross=(sig.groupby(["source","term_id","term_name"]).held_out_dataset.nunique().reset_index(name="n_folds").sort_values(["n_folds","term_name"],ascending=[False,True])); cross.to_csv(OUT/"03_cross_fold_enrichment_recurrence.csv",index=False)
    else: cross=pd.DataFrame(columns=["source","term_id","term_name","n_folds"]);cross.to_csv(OUT/"03_cross_fold_enrichment_recurrence.csv",index=False)
    rows=[]
    for held,g in diag.groupby("held_out_dataset"):
        g=g.sort_values("gene_set_size"); base=int(g.iloc[0].n_fdr05_terms) if len(g) else 0; largest=int(g.iloc[-1].n_fdr05_terms) if len(g) else 0
        rows.append({"held_out_dataset":held,"n_gene_set_sizes":len(g),"min_fdr05_terms":base,"max_fdr05_terms":largest,"first_nonzero_size":next((int(r.gene_set_size) for _,r in g.iterrows() if int(r.n_fdr05_terms)>0),np.nan),"has_any_enrichment":bool((g.n_fdr05_terms>0).any())})
    summary=pd.DataFrame(rows);summary.to_csv(OUT/"04_fold_interpretation.csv",index=False)
    overall=pd.DataFrame([{"n_folds":len(folds),"n_folds_with_any_enrichment":int(summary.has_any_enrichment.sum()) if len(summary) else 0,"n_folds_with_enrichment_at_2500":int((diag[diag.gene_set_size==2500].groupby("held_out_dataset").n_fdr05_terms.max()>0).sum()) if len(diag) else 0,"n_cross_fold_terms":int(len(cross)),"n_terms_recurrent_2plus_folds":int((cross.n_folds>=2).sum()) if len(cross) else 0}]);overall.to_csv(OUT/"05_stage2912_summary.csv",index=False)
    log("complete")
    print("\nStage 2.9.12 gene-set sensitivity:",flush=True);print(diag.to_string(index=False),flush=True)
    print("\nStage 2.9.12 fold interpretation:",flush=True);print(summary.to_string(index=False),flush=True)
    print("\nStage 2.9.12 overall:",flush=True);print(overall.to_string(index=False),flush=True)
    return overall

if __name__=="__main__":run()
