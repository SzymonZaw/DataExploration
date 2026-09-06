"""Stage 2.9.7: biologically anchored consensus state validation."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]; COMMON=ROOT/"results"/"Dynamics"/"stage2_6"; IN=ROOT/"results"/"Dynamics"/"stage2_9_6"; OUT=ROOT/"results"/"Dynamics"/"stage2_9_7"; OUT.mkdir(parents=True,exist_ok=True)

def _corr(a,b,method="spearman"):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method=method))

def _load():
    m=pd.read_csv(COMMON/"06_common_human_gene_matrix.csv",index_col=0).apply(pd.to_numeric,errors="coerce"); meta=pd.read_csv(COMMON/"07_common_gene_sample_metadata.csv"); rec=pd.read_csv(IN/"02_gene_recurrence.csv"); return m,meta,rec

def _time(meta,ds,c):
    q=meta[(meta["dataset"].astype(str)==str(ds))&(meta["matrix_column"].astype(str)==str(c))]
    if q.empty:return np.nan
    if "time_hours" in q.columns and pd.notna(q.iloc[0]["time_hours"]):return float(q.iloc[0]["time_hours"])
    try:
        from .validation import _time_hours_for_validation,_strip_dataset_prefix
        sample=q.iloc[0]["sample"]; local=int(q.index[0])
        if str(ds)=="GSE28688":local=int(meta.index[meta["dataset"].astype(str)==str(ds)].tolist().index(q.index[0]))
        v=_time_hours_for_validation(str(ds),_strip_dataset_prefix(sample),local if str(ds)=="GSE28688" else None);return float(v) if pd.notna(v) else np.nan
    except Exception:return np.nan

def _trajectory(ds,genes):
    g=META[(META["dataset"].astype(str)==str(ds))&META["matrix_column"].notna()].copy();g["t"]=[_time(META,ds,c) for c in g["matrix_column"]];g=g[np.isfinite(g["t"].to_numpy(float))]
    if len(g)<3 or g["t"].nunique()<3:return None
    genes=[x for x in genes if x in MATRIX.index];
    if not genes:return None
    f=MATRIX.loc[genes,g["matrix_column"].astype(str).tolist()].T;f.index=g["t"].to_numpy(float);f=f.groupby(level=0).mean().sort_index();
    if len(f)<3:return None
    t=f.index.to_numpy(float);return (t-t[0])/(t[-1]-t[0]),f.to_numpy(float)

def _stable_genes(train,recurrence=2,min_temporal_corr=.5,max_genes=1000):
    genes=set(REC.loc[REC["n_discovery_folds"]>=recurrence,"gene"].astype(str));vals=[]
    for g in sorted(genes):
        cors=[]
        for ds in train:
            tr=_trajectory(ds,[g])
            if tr is not None:cors.append(abs(_corr(tr[0],tr[1][:,0],"spearman")))
        if cors and np.mean(cors)>=min_temporal_corr:vals.append((g,float(np.mean(cors)),len(cors)))
    vals=sorted(vals,key=lambda x:(x[1],x[2]),reverse=True)[:max_genes];return [x[0] for x in vals],pd.DataFrame(vals,columns=["gene","mean_abs_temporal_spearman","n_training_datasets"])

def _modules(genes,n_modules=8):
    profiles=[];used=[];grid=np.linspace(0,1,9)
    for g in genes:
        arr=[]
        for ds in TRAIN:
            tr=_trajectory(ds,[g])
            if tr is not None:arr.append(np.interp(grid,tr[0],tr[1][:,0]))
        if len(arr)>=2:
            v=np.nanmean(np.asarray(arr),axis=0);v=(v-v.mean())/(v.std()+1e-12);profiles.append(v);used.append(g)
    if len(used)<n_modules:return None
    score=np.mean(np.asarray(profiles)*np.linspace(-1,1,9),axis=1);order=np.argsort(score);chunks=np.array_split(order,min(n_modules,len(order)));return {f"module_{i+1:02d}":[used[j] for j in ch] for i,ch in enumerate(chunks) if len(ch)}

def _activity(ds,modules):
    genes=sorted(set(g for gs in modules.values() for g in gs));tr=_trajectory(ds,genes)
    if tr is None:return None
    t,y=tr;idx={g:i for i,g in enumerate(genes)};out={"t":t}
    for name,gs in modules.items():
        ids=[idx[g] for g in gs];v=y[:,ids];mu=np.nanmedian(v,axis=0);mad=1.4826*np.nanmedian(np.abs(v-mu),axis=0);mad[~np.isfinite(mad)|(mad<1e-8)]=1;z=(np.where(np.isfinite(v),v,mu)-mu)/mad;out[name]=np.nanmedian(z,axis=1)
    return pd.DataFrame(out).set_index("t")

def _fit_axis(train,modules):
    frames=[]
    for ds in train:
        a=_activity(ds,modules)
        if a is not None and len(a)>=3:
            grid=np.linspace(0,1,9);frames.append(np.column_stack([np.interp(grid,a.index.to_numpy(float),a[c].to_numpy(float)) for c in a.columns]))
    if len(frames)<2:return None
    X=np.concatenate(frames);sc=StandardScaler().fit(X);pc=PCA(n_components=1,random_state=42).fit(sc.transform(X));z=pc.transform(sc.transform(X))[:,0];target=np.tile(np.linspace(0,1,9),len(frames));r=_corr(z,target,"pearson");return sc,pc,(1 if np.isfinite(r) and r>=0 else -1),z,r

def run(recurrence=2,min_temporal_corr=.5,n_modules=8,max_genes=1000):
    global MATRIX,META,REC,TRAIN
    MATRIX,META,REC=_load();folds=sorted(META["dataset"].astype(str).unique());rows=[];membership=[];stable_rows=[]
    for ds in folds:
        TRAIN=[x for x in folds if x!=ds];genes,diag=_stable_genes(TRAIN,recurrence,min_temporal_corr,max_genes);diag.insert(0,"held_out_dataset",ds);stable_rows.append(diag);modules=_modules(genes,n_modules);test=_trajectory(ds,genes)
        if modules is None or test is None:rows.append({"held_out_dataset":ds,"status":"insufficient_consensus_modules","n_stable_genes":len(genes),"n_modules":0});continue
        for name,gs in modules.items():
            for g in gs:membership.append({"held_out_dataset":ds,"module":name,"gene":g})
        axis=_fit_axis(TRAIN,modules)
        if axis is None:rows.append({"held_out_dataset":ds,"status":"insufficient_training_state","n_stable_genes":len(genes),"n_modules":len(modules)});continue
        sc,pc,sign,ref,train_corr=axis;act=_activity(ds,modules);X=act.to_numpy(float);med=np.nanmedian(X,axis=0);X=np.where(np.isfinite(X),X,med);X=np.where(np.isfinite(X),X,0);z=sign*pc.transform(sc.transform(X))[:,0];lo,hi=np.nanpercentile(ref,[2.5,97.5]);zn=np.clip((z-lo)/(hi-lo),0,1) if hi>lo else np.full(len(z),.5);t=act.index.to_numpy(float);tn=(t-t.min())/(t.max()-t.min())
        rows.append({"held_out_dataset":ds,"status":"ok","n_stable_genes":len(genes),"n_modules":len(modules),"train_axis_time_pearson":train_corr,"progress_spearman":_corr(tn,zn),"progress_pearson":_corr(tn,zn,"pearson"),"progress_rmse":float(np.sqrt(np.mean((tn-zn)**2))),"progress_monotonicity":float(np.mean(np.diff(zn)>=-.02))})
    result=pd.DataFrame(rows);result.to_csv(OUT/"01_fold_results.csv",index=False);pd.concat(stable_rows,ignore_index=True).to_csv(OUT/"02_stable_gene_diagnostics.csv",index=False);pd.DataFrame(membership).to_csv(OUT/"03_module_membership.csv",index=False);ok=result[result.status=="ok"]
    summary=pd.DataFrame([{"n_folds":len(result),"n_valid_folds":len(ok),"mean_progress_spearman":ok.progress_spearman.mean() if len(ok) else np.nan,"mean_progress_pearson":ok.progress_pearson.mean() if len(ok) else np.nan,"mean_progress_rmse":ok.progress_rmse.mean() if len(ok) else np.nan,"mean_monotonicity":ok.progress_monotonicity.mean() if len(ok) else np.nan,"mean_train_axis_time_pearson":ok.train_axis_time_pearson.mean() if len(ok) else np.nan,"interpretation":"recurrence + temporal-coherence consensus modules; no ODE fitted"}]);summary.to_csv(OUT/"04_overall_summary.csv",index=False);print("Stage 2.9.7 complete.");print(result.to_string(index=False));print(summary.to_string(index=False));return summary

if __name__=="__main__":run()
