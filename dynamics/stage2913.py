"""Stage 2.9.13: bootstrap stability of biological enrichment.

Tests whether enrichment observed for selected training genes is robust to
subsampling rather than to one arbitrary top-N cutoff.  This is diagnostic
only: gene/program discovery is not used for ODE fitting and no held-out
sample is used to define programs.
"""
from pathlib import Path
import json, os, ssl, time, urllib.error, urllib.request
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
COMMON=ROOT/"results"/"Dynamics"/"stage2_6"
IN=ROOT/"results"/"Dynamics"/"stage2_9_1"
OUT=ROOT/"results"/"Dynamics"/"stage2_9_13"
CACHE=OUT/"cache"
OUT.mkdir(parents=True,exist_ok=True); CACHE.mkdir(parents=True,exist_ok=True)

def log(x): print(f"Stage 2.9.13: {x}",flush=True)
def _ctx(verify=True):
    if not verify:return ssl._create_unverified_context()
    try:
        import certifi; return ssl.create_default_context(cafile=certifi.where())
    except ImportError:return ssl.create_default_context()

def _request(genes,bg,label):
    path=CACHE/f"{label}.json"
    if path.exists():
        try:return json.loads(path.read_text(encoding="utf-8"))
        except Exception:pass
    payload={"organism":"hsapiens","query":genes,"sources":["GO:BP","REAC","KEGG"],"domain_scope":"custom","background":bg,"user_threshold":0.05,"significance_threshold_method":"fdr","no_evidences":False,"output":"json"}
    body=json.dumps(payload).encode(); insecure=os.environ.get("STAGE298_INSECURE_SSL","")=="1"
    for attempt in range(3):
        try:
            req=urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",data=body,headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"DataExploration-stage2.9.13/1.0"},method="POST")
            with urllib.request.urlopen(req,timeout=90,context=_ctx(True)) as r:data=json.loads(r.read().decode())
            path.write_text(json.dumps(data),encoding="utf-8");return data
        except urllib.error.URLError as e:
            if insecure and attempt==0 and isinstance(e.reason,ssl.SSLCertVerificationError):
                log(f"{label}: verified TLS failed; retrying with verification disabled")
                try:
                    with urllib.request.urlopen(req,timeout=90,context=_ctx(False)) as r:data=json.loads(r.read().decode())
                    path.write_text(json.dumps(data),encoding="utf-8");return data
                except Exception:pass
            if attempt<2:time.sleep(2**attempt)
        except Exception:
            if attempt<2:time.sleep(2**attempt)
    return None

def _stable():
    stab=pd.read_csv(IN/"03_gene_stability.csv"); status=pd.read_csv(IN/"06_fold_status.csv")
    matrix=pd.read_csv(COMMON/"06_common_human_gene_matrix.csv",index_col=0); genes=matrix.index.astype(str).tolist();out={}
    for _,r in status.iterrows():
        ds=str(r.held_out_dataset)
        if str(r.get("reason",""))!="ok":continue
        s=stab[stab.held_out_dataset.astype(str)==ds].copy();s["stability"]=pd.to_numeric(s.median_pairwise_correlation,errors="coerce");s=s[s.stability.notna()]
        q=s[s.n_datasets>=3].copy()
        if len(q)<20:q=s[s.n_datasets>=2].copy()
        if len(q)<20:continue
        q=q[q.stability>=float(q.stability.quantile(.75))].sort_values(["stability","n_datasets"],ascending=False)
        out[ds]=[genes[i] for i in q.gene_index.astype(int) if 0<=i<len(genes)]
    return out,genes

def _terms(data):
    rows=[]
    if not data:return rows
    for r in data.get("result",[]):
        p=pd.to_numeric(r.get("p_value"),errors="coerce")
        if np.isfinite(p):rows.append((str(r.get("source","")),str(r.get("native","")),str(r.get("name","")),float(p),int(r.get("intersection_size",0) or 0)))
    return rows

def run(n_boot=30,subsample=1000,seed=2913):
    log("starting enrichment bootstrap stability diagnostic; no ODE/state-space model")
    matrix=pd.read_csv(COMMON/"06_common_human_gene_matrix.csv",index_col=0);bg=matrix.index.astype(str).tolist();folds,_=_stable();log(f"loaded {len(bg):,} background genes and {len(folds)} training folds")
    rng=np.random.default_rng(seed); summary=[]; term_rows=[]; overlap_rows=[]
    for held,genes in sorted(folds.items()):
        if len(genes)<subsample: subsample_use=len(genes)
        else:subsample_use=subsample
        log(f"FOLD {held}: {len(genes):,} stable genes; {n_boot} bootstrap samples of {subsample_use:,}")
        sets=[]
        for b in range(n_boot):
            q=list(rng.choice(np.array(genes,dtype=object),size=subsample_use,replace=False))
            data=_request(q,bg,f"{held}_boot_{b:03d}_{subsample_use}"); ts=_terms(data); sig={(a,b) for a,b,_,p,_ in ts if p<.05};sets.append(sig)
            for src,tid,name,p,inter in ts:
                if p<.05:term_rows.append({"held_out_dataset":held,"bootstrap":b,"source":src,"term_id":tid,"term_name":name,"p_value":p,"intersection_size":inter})
            if (b+1)%5==0:log(f"{held}: bootstrap {b+1}/{n_boot}")
        counts={}
        for s in sets:
            for key in s:counts[key]=counts.get(key,0)+1
        freq=np.array(list(counts.values()),dtype=float)/n_boot if counts else np.array([])
        stable_terms={k:v for k,v in counts.items() if v/n_boot>=.5}
        pair=[]
        for i in range(len(sets)):
            for j in range(i+1,len(sets)):
                u=sets[i]|sets[j]; pair.append(len(sets[i]&sets[j])/len(u) if u else np.nan)
        summary.append({"held_out_dataset":held,"n_stable_genes":len(genes),"bootstrap_subsample_size":subsample_use,"n_bootstrap":n_boot,"mean_fdr05_terms_per_bootstrap":float(np.mean([len(s) for s in sets])) if sets else np.nan,"median_fdr05_terms_per_bootstrap":float(np.median([len(s) for s in sets])) if sets else np.nan,"n_unique_fdr05_terms":len(counts),"n_terms_recurrence_ge50pct":len(stable_terms),"max_term_recurrence":float(freq.max()) if len(freq) else 0.0,"mean_pairwise_jaccard":float(np.nanmean(pair)) if pair else np.nan})
        for (src,tid),n in sorted(counts.items(),key=lambda x:x[1],reverse=True)[:100]:overlap_rows.append({"held_out_dataset":held,"source":src,"term_id":tid,"bootstrap_count":n,"bootstrap_frequency":n/n_boot})
    s=pd.DataFrame(summary);tr=pd.DataFrame(term_rows);ov=pd.DataFrame(overlap_rows)
    s.to_csv(OUT/"01_bootstrap_summary.csv",index=False);tr.to_csv(OUT/"02_bootstrap_enriched_terms.csv",index=False);ov.to_csv(OUT/"03_term_bootstrap_frequency.csv",index=False)
    cross=(ov[ov.bootstrap_frequency>=.5].groupby(["source","term_id"]).held_out_dataset.nunique().reset_index(name="n_folds") if len(ov) else pd.DataFrame(columns=["source","term_id","n_folds"]))
    cross.to_csv(OUT/"04_cross_fold_stable_terms.csv",index=False)
    overall=pd.DataFrame([{"n_folds":len(folds),"n_folds_with_bootstrap_stable_terms":int((s.n_terms_recurrence_ge50pct>0).sum()) if len(s) else 0,"mean_pairwise_jaccard":float(s.mean_pairwise_jaccard.mean()) if len(s) else np.nan,"n_cross_fold_stable_terms":len(cross),"n_cross_fold_stable_terms_2plus":int((cross.n_folds>=2).sum()) if len(cross) else 0}]);overall.to_csv(OUT/"05_stage2913_summary.csv",index=False)
    log("complete");print("\nStage 2.9.13 bootstrap summary:");print(s.to_string(index=False));print("\nStage 2.9.13 overall:");print(overall.to_string(index=False));return overall

if __name__=="__main__":run()
