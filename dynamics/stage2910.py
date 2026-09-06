"""Stage 2.9.10: leakage-free biological program validation.

For each held-out dataset, stable genes are taken only from the corresponding
Stage 2.9.1 training-fold discovery. Biological enrichment is then repeated on
that training gene set, and program gene sets are projected to the held-out
dataset without using held-out data for program discovery or scaling.

This stage does not fit an ODE/state-space model.
"""
from pathlib import Path
import json, os, ssl, time, urllib.error, urllib.request
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; COMMON=ROOT/"results"/"Dynamics"/"stage2_6"; IN=ROOT/"results"/"Dynamics"/"stage2_9_1"; OUT=ROOT/"results"/"Dynamics"/"stage2_9_10"; CACHE=OUT/"cache"
OUT.mkdir(parents=True,exist_ok=True); CACHE.mkdir(parents=True,exist_ok=True)

def log(x): print(f"Stage 2.9.10: {x}",flush=True)
def corr(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(np.corrcoef(a[ok],b[ok])[0,1])
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
            req=urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",data=body,headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"DataExploration-stage2.9.10/1.2"},method="POST")
            with urllib.request.urlopen(req,timeout=90,context=sslctx(True)) as r:data=json.loads(r.read().decode())
            path.write_text(json.dumps(data),encoding="utf-8"); return data
        except urllib.error.URLError as e:
            if insecure and attempt==0 and isinstance(e.reason,ssl.SSLCertVerificationError):
                log(f"{label}: verified TLS failed; retrying with verification disabled")
                try:
                    req=urllib.request.Request("https://biit.cs.ut.ee/gprofiler/api/gost/profile/",data=body,headers={"Content-Type":"application/json","Accept":"application/json","User-Agent":"DataExploration-stage2.9.10/1.2"},method="POST")
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
        s=stab[stab.held_out_dataset.astype(str)==ds].copy(); s=s[pd.to_numeric(s.median_pairwise_correlation,errors="coerce").notna()]
        q=s[s.n_datasets>=3].copy()
        if len(q)<20:q=s[s.n_datasets>=2].copy()
        if len(q)<20:continue
        threshold=float(q.median_pairwise_correlation.quantile(.75)); q=q[q.median_pairwise_correlation>=threshold]
        q=q.sort_values(["median_pairwise_correlation","n_datasets"],ascending=False).head(2500)
        out[ds]=[genes[i] for i in q.gene_index.astype(int) if 0<=i<len(genes)]
    return out

def enrich(genes,background,label):
    payload={"organism":"hsapiens","query":genes,"sources":["GO:BP","REAC","KEGG"],"domain_scope":"custom","background":background,"user_threshold":0.05,"significance_threshold_method":"fdr","no_evidences":False,"output":"json"}
    log(f"fold {label}: g:Profiler on {len(genes):,} training genes"); return request(payload,f"gprofiler_{label}")

def programs(data,query,background,max_programs=8):
    if not data:return []
    bg={str(x).upper() for x in background}; rows=[]
    for r in data.get("result",[]):
        p=float(r.get("p_value",np.nan)); n=int(r.get("intersection_size",0) or 0); name=str(r.get("name","") or "")
        if not np.isfinite(p) or p>=.05 or n<5:continue
        low=name.lower()
        if any(k in low for k in ("hiv","viral messenger","metabolism of rna","gene expression","protein metabolic process","rna polymerase ii transcription","processing of capped intron-containing pre-mrna")):continue
        mask=r.get("intersections",[]); gs=[str(query[i]).upper() for i,x in enumerate(mask) if x and i<len(query)]; gs=sorted(set(gs)&bg)
        if len(gs)>=5:rows.append((p,n,name,r.get("native"),gs))
    rows.sort(key=lambda x:(x[0],-x[1])); sel=[]
    for row in rows:
        s=set(row[4])
        if any(len(s&set(q[4]))/len(s|set(q[4]))>=.75 for q in sel):continue
        sel.append(row)
        if len(sel)>=max_programs:break
    return sel

def activity(matrix,ps,fit_cols):
    idx={str(g).upper():g for g in matrix.index}; fit=[str(c) for c in fit_cols if str(c) in matrix.columns]; out={}
    for pi,p in enumerate(ps,1):
        gs=[idx[g] for g in p[4] if g in idx]
        if len(gs)<5:continue
        train=matrix.loc[gs,fit].to_numpy(float); vals=matrix.loc[gs].to_numpy(float)
        mu=np.nanmedian(train,axis=1); mad=1.4826*np.nanmedian(np.abs(train-mu[:,None]),axis=1); mad[~np.isfinite(mad)|(mad<1e-8)]=1.0
        z=(np.where(np.isfinite(vals),vals,mu[:,None])-mu[:,None])/mad[:,None]
        out[f"P{pi:02d}"]=np.nanmedian(z,axis=0)
    return pd.DataFrame(out,index=matrix.columns)

def run():
    log("starting leakage-free biological program validation; no ODE/state-space model")
    matrix,meta=load_space(); background=matrix.index.astype(str).tolist(); fold_genes=load_fold_genes(); datasets=sorted(fold_genes)
    log(f"valid training folds available: {len(datasets)}"); all_results=[]; defs=[]; trajectories_out=[]
    for i,held in enumerate(datasets,1):
        log(f"FOLD {i}/{len(datasets)}: held out {held}"); query=fold_genes[held]
        data=enrich(query,background,held); ps=programs(data,query,background); log(f"fold {held}: selected {len(ps)} training-derived biological programs")
        for j,p in enumerate(ps,1):defs.append({"held_out_dataset":held,"program_id":f"P{j:02d}","term_id":p[3],"term_name":p[2],"p_value":p[0],"n_genes":len(p[4]),"genes":",".join(p[4])})
        if len(ps)<2:continue
        train_ds=[d for d in meta.dataset.astype(str).unique() if d!=held]; fit_cols=meta[meta.dataset.astype(str).isin(train_ds)&meta.time_hours.notna()]["matrix_column"].astype(str).tolist(); act=activity(matrix,ps,fit_cols)
        trajectories={}
        for ds,g in meta.groupby("dataset",sort=True):
            g=g[g.time_hours.notna()&g.matrix_column.notna()]
            if g.time_hours.nunique()<3:continue
            rec=[]
            for t,gt in g.groupby("time_hours",sort=True):
                cols=[str(c) for c in gt.matrix_column if str(c) in act.index]
                if cols:rec.append((float(t),act.loc[cols].mean(axis=0).to_numpy(float)))
            if len(rec)>=3:trajectories[ds]=(np.array([x[0] for x in rec]),np.vstack([x[1] for x in rec]))
        if held not in trajectories:continue
        tt,tv=trajectories[held]; train={d:v for d,v in trajectories.items() if d!=held}
        for j,t in enumerate(tt):
            preds=[]
            for _,(x,y) in train.items():
                if x.min()<=t<=x.max():preds.append(np.array([np.interp(t,x,y[:,k]) for k in range(y.shape[1])]))
            if not preds:continue
            pred=np.mean(preds,axis=0); true=tv[j]; all_results.append({"held_out_dataset":held,"time_hours":t,"n_programs":len(ps),"n_training_datasets":len(preds),"rmse":float(np.sqrt(np.mean((true-pred)**2))),"mae":float(np.mean(np.abs(true-pred))),"profile_correlation":corr(true,pred)})
        for t,y in zip(tt,tv):
            for k,val in enumerate(y,1):trajectories_out.append({"held_out_dataset":held,"time_hours":t,"program_id":f"P{k:02d}","activity":val})
    pd.DataFrame(defs).to_csv(OUT/"01_fold_program_definitions.csv",index=False); res=pd.DataFrame(all_results);res.to_csv(OUT/"02_lodo_by_timepoint.csv",index=False)
    if len(res):summary=res.groupby("held_out_dataset").agg(n_timepoints=("time_hours","count"),mean_rmse=("rmse","mean"),mean_mae=("mae","mean"),mean_profile_correlation=("profile_correlation","mean")).reset_index()
    else:summary=pd.DataFrame(columns=["held_out_dataset","n_timepoints","mean_rmse","mean_mae","mean_profile_correlation"])
    summary.to_csv(OUT/"03_lodo_summary.csv",index=False);pd.DataFrame(trajectories_out).to_csv(OUT/"04_heldout_program_trajectories.csv",index=False)
    overall=pd.DataFrame([{ "n_folds":len(datasets),"n_valid_folds":int(summary.shape[0]),"mean_rmse":summary.mean_rmse.mean() if len(summary) else np.nan,"mean_profile_correlation":summary.mean_profile_correlation.mean() if len(summary) else np.nan}]);overall.to_csv(OUT/"05_overall_summary.csv",index=False)
    log("complete; enrichment, program definitions and activity scaling were training-fold only")
    print("\nStage 2.9.10 leakage-free LODO:",flush=True);print(summary.to_string(index=False),flush=True);print("\nStage 2.9.10 overall:",flush=True);print(overall.to_string(index=False),flush=True);return summary

if __name__=="__main__":run()
