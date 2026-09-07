"""Stage 2.9.30: time-warped shared trajectory validation."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results"/"Dynamics"/"stage2_9_30"; OUT.mkdir(parents=True,exist_ok=True)
TARGET=["GSE67462","GSE28688","GSE297234"]; GRID=np.linspace(0,1,40); MIN_GENES=100; MAX_GENES=3000; N_PERM=100; N_PC=10

def log(x): print(f"Stage 2.9.30: {x}",flush=True)
def finite(X):
 X=np.asarray(X,float).copy()
 for j in range(X.shape[1]):
  ok=np.isfinite(X[:,j]); X[~ok,j]=np.median(X[ok,j]) if ok.any() else 0.0
 return X
def normalize_time(t):
 t=np.asarray(t,float); lo,hi=float(t.min()),float(t.max()); return np.zeros_like(t) if hi<=lo else (t-lo)/(hi-lo)
def resample(t,X):
 tn=normalize_time(t); X=finite(X); Y=np.empty((len(GRID),X.shape[1]))
 for j in range(X.shape[1]): Y[:,j]=np.interp(GRID,tn,X[:,j])
 return Y
def normalize_with_stats(Y,baseline=None,amp=None):
 Y=finite(Y); baseline=Y[0].copy() if baseline is None else baseline; Z=Y-baseline[None,:]
 if amp is None: amp=np.sqrt(np.mean(Z*Z,axis=0)); amp=np.where(amp>1e-8,amp,1.0)
 return Z/amp,baseline,amp
def gene_indices(genes,all_genes):
 pos={str(g):i for i,g in enumerate(all_genes)}; return [pos[str(g)] for g in genes if str(g) in pos]
def select_genes(train,all_genes):
 A=np.stack([normalize_with_stats(resample(t,X))[0] for t,X,_ in train.values()]); mean=A.mean(0); shared=np.var(mean,0); hetero=np.mean((A-mean[None])**2,axis=(0,1)); score=shared/(shared+hetero+1e-12); order=np.argsort(score)[::-1]; chosen=order[score[order]>=0.5]
 if len(chosen)<MIN_GENES: chosen=order[:min(MIN_GENES,len(order))]
 return np.asarray(all_genes)[chosen[:MAX_GENES]]
def load_trajectories():
    # Reuse Stage 2.7's validated sample/time/matrix mapping.
    from dynamics.validation import _load_common_space
    m,meta=_load_common_space()
    meta=meta[meta["dataset"].astype(str).isin(TARGET)].copy()
    out={}
    for ds in TARGET:
        g=meta[(meta["dataset"].astype(str)==ds)&meta["matrix_column"].notna()&meta["time_hours"].notna()].copy()
        g["time_hours"]=pd.to_numeric(g["time_hours"],errors="coerce")
        g=g[np.isfinite(g["time_hours"].to_numpy(float))]
        if g["time_hours"].nunique()<3:
            log(f"{ds}: skipped; timed samples={len(g)}, unique times={g['time_hours'].nunique()}")
            continue
        expr=m.loc[:,g["matrix_column"].astype(str).tolist()].T.copy()
        expr.index=g["time_hours"].to_numpy(float)
        expr=expr.groupby(level=0,sort=True).mean()
        log(f"{ds}: {len(g)} timed samples -> {len(expr)} unique timepoints")
        out[ds]=(expr.index.to_numpy(float),expr.to_numpy(float),list(m.index))
    return out
def fit_artifact(train,all_genes):
 genes=select_genes(train,all_genes); idx=gene_indices(genes,all_genes); states=[normalize_with_stats(resample(t,X[:,idx]))[0] for t,X,_ in train.values()]; Z=np.vstack(states); pca=PCA(n_components=min(N_PC,Z.shape[0],Z.shape[1]),random_state=2930).fit(Z); template=np.mean(np.stack([pca.transform(Y) for Y in states]),0); return genes,pca,template
def dtw_path(A,B):
 n,m=len(A),len(B); dp=np.full((n+1,m+1),np.inf); dp[0,0]=0; prev=np.full((n+1,m+1),-1,dtype=np.int8)
 for i in range(1,n+1):
  for j in range(1,m+1):
   d=float(np.sqrt(np.mean((A[i-1]-B[j-1])**2))); c=(dp[i-1,j-1],dp[i-1,j],dp[i,j-1]); k=int(np.argmin(c)); dp[i,j]=d+c[k]; prev[i,j]=k
 i,j=n,m; path=[]
 while i and j:
  path.append((i-1,j-1)); k=prev[i,j]; i,j=(i-1,j-1) if k==0 else (i-1,j) if k==1 else (i,j-1)
 return path[::-1],float(dp[n,m]/max(1,len(path)))
def infer_progress(prefix_pc,T):
 path,cost=dtw_path(prefix_pc,T); mp={i:[] for i in range(len(prefix_pc))}
 for i,j in path: mp[i].append(j)
 q=[];p=[]
 for i,js in mp.items():
  if js:q.append(GRID[i]);p.append(float(np.mean(GRID[js])))
 q,p=np.asarray(q),np.asarray(p)
 if len(q)<3:return np.nan,cost,np.nan
 slope,intercept=np.polyfit(q,p,1); return float(np.clip(slope+intercept,0,1)),cost,float(slope)
def interp_state(T,p): return np.array([np.interp(float(np.clip(p,0,1)),GRID,T[:,j]) for j in range(T.shape[1])])
def metrics(a,b): return float(np.sqrt(np.mean((a-b)**2))),float(pd.Series(a).corr(pd.Series(b)))
def prepare_fold(heldout,train_names,traj):
 train={d:traj[d] for d in train_names}; t,X,all_genes=traj[heldout]; genes,pca,T=fit_artifact(train,all_genes); idx=gene_indices(genes,all_genes); prefix_raw=resample(t[:-1],X[:-1,idx]); prefix_state,base,amp=normalize_with_stats(prefix_raw); full_state,_,_=normalize_with_stats(resample(t,X[:,idx]),base,amp); prefix_pc=pca.transform(prefix_state); true_pc=pca.transform(full_state)[-1]; persistence=prefix_pc[-1]; return {"genes":genes,"pca":pca,"template":T,"prefix_pc":prefix_pc,"true_pc":true_pc,"persistence":persistence}
def evaluate(a,seed=None):
 prefix=a["prefix_pc"]
 if seed is not None: prefix=prefix[np.random.default_rng(seed).permutation(len(prefix))]
 prog,cost,slope=infer_progress(prefix,a["template"])
 if not np.isfinite(prog):return None
 pred=interp_state(a["template"],prog); true=a["true_pc"]; rw,cw=metrics(pred,true); rp,cp=metrics(a["persistence"],true); ru,cu=metrics(a["template"][-1],true)
 return {"warp_progress_at_final_time":prog,"warp_slope":slope,"dtw_prefix_cost":cost,"rmse_time_warp":rw,"rmse_persistence":rp,"rmse_unwarped_shared":ru,"pearson_time_warp":cw,"pearson_persistence":cp,"pearson_unwarped_shared":cu,"improvement_vs_persistence":rp-rw,"improvement_vs_unwarped_shared":ru-rw}
def run():
 traj=load_trajectories(); names=[d for d in TARGET if d in traj]; log(f"trajectory datasets: {', '.join(names)}")
 if len(names)<3:raise RuntimeError("Stage 2.9.30 requires at least three trajectory datasets.")
 arts={}; rows=[]
 for heldout in names:
  train=[d for d in names if d!=heldout]; log(f"LODO held out {heldout}; training on {', '.join(train)}"); a=prepare_fold(heldout,train,traj); arts[heldout]=a; r=evaluate(a)
  if r:r.update({"heldout_dataset":heldout,"n_selected_genes":len(a["genes"])}); rows.append(r)
 R=pd.DataFrame(rows); R.to_csv(OUT/"01_lodo_time_warp_validation.csv",index=False)
 if R.empty:raise RuntimeError("No valid Stage 2.9.30 LODO folds.")
 null=[]
 for p in range(N_PERM):
  vals=[]
  for di,ds in enumerate(names):
   r=evaluate(arts[ds],700000+p*len(names)+di)
   if r:vals.append(r["improvement_vs_persistence"])
  null.append({"permutation":p+1,"mean_improvement_vs_persistence":float(np.mean(vals)) if vals else np.nan})
 P=pd.DataFrame(null); P.to_csv(OUT/"02_time_permutation_null.csv",index=False); obs=float(R["improvement_vs_persistence"].mean()); nv=P["mean_improvement_vs_persistence"].dropna().to_numpy(float); pv=float((1+np.sum(nv>=obs))/(len(nv)+1)) if len(nv) else np.nan
 S=pd.DataFrame([{"n_trajectory_datasets":len(names),"n_valid_lodo_folds":len(R),"mean_rmse_time_warp":float(R["rmse_time_warp"].mean()),"mean_rmse_persistence":float(R["rmse_persistence"].mean()),"mean_rmse_unwarped_shared":float(R["rmse_unwarped_shared"].mean()),"mean_improvement_vs_persistence":obs,"mean_improvement_vs_unwarped_shared":float(R["improvement_vs_unwarped_shared"].mean()),"mean_pearson_time_warp":float(R["pearson_time_warp"].mean()),"time_warp_permutation_p":pv,"time_warp_predictive_support":bool(obs>0 and np.isfinite(pv) and pv<0.05),"stage3_readiness":False,"interpretation":"Training-only shared trajectory with held-out-prefix DTW time warping; final-point LODO prediction; representation/predictive diagnostic only; no ODE/state-space model"}]); S.to_csv(OUT/"03_stage2930_summary.csv",index=False); log("overall:"); print(S.to_string(index=False)); return S
if __name__=="__main__":run()
