"""Stage 2.9.16: residual biological-state validation.

Question: do fixed biological programs contain information beyond explicit time?
For each program, remove its smooth dependence on normalized time using a
training-only linear fit, then test whether the residual program state carries
cross-dataset signal. This is a diagnostic gate before Stage 3, not an ODE fit.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"results"/"Dynamics"/"stage2_9_14"
OUT=ROOT/"results"/"Dynamics"/"stage2_9_16"
OUT.mkdir(parents=True,exist_ok=True)

def log(x): print(f"Stage 2.9.16: {x}",flush=True)
def corr(a,b,method="pearson"):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method="spearman" if method=="spearman" else "pearson"))
def _load():
    p=SRC/"03_program_trajectories.csv"
    if not p.exists(): raise RuntimeError("Stage 2.9.14 trajectories missing. Run: python validate_pipeline.py --stage2914")
    x=pd.read_csv(p); x["dataset"]=x.dataset.astype(str); x["program_id"]=x.program_id.astype(str)
    return x

def _residuals(df):
    rows=[]
    for (ds,pid),g in df.groupby(["dataset","program_id"],sort=True):
        g=g.sort_values("time_hours").copy(); t=g.time_hours.to_numpy(float); y=g.activity.to_numpy(float)
        if len(g)<3 or np.ptp(t)<=0: continue
        tn=(t-t.min())/np.ptp(t)
        X=np.column_stack([np.ones(len(tn)),tn]); beta=np.linalg.lstsq(X,y,rcond=None)[0]; fit=X@beta
        for ti,yi,fi,ri,tn_i in zip(t,y,fit,y-fit,tn):
            rows.append({"dataset":ds,"program_id":pid,"time_hours":float(ti),"normalized_time":float(tn_i),"activity":float(yi),"time_fitted_activity":float(fi),"residual_activity":float(ri),"time_slope":float(beta[1])})
    return pd.DataFrame(rows)

def _lodo_residual(df):
    rows=[]
    for held,gtest in df.groupby("dataset",sort=True):
        train=df[df.dataset!=held]
        common=sorted(set(gtest.program_id)&set(train.program_id))
        for pid in common:
            tr=train[train.program_id==pid].dropna(subset=["normalized_time","residual_activity"])
            te=gtest[gtest.program_id==pid].dropna(subset=["normalized_time","residual_activity"])
            if len(tr)<3 or len(te)<3: continue
            # Cross-dataset residual calibration: fit residual-vs-time on training
            # data and project the held-out dataset. A non-zero slope is itself
            # diagnostic; no biological dynamics model is inferred here.
            X=np.column_stack([np.ones(len(tr)),tr.normalized_time]); b=np.linalg.lstsq(X,tr.residual_activity,rcond=None)[0]
            pred=np.column_stack([np.ones(len(te)),te.normalized_time])@b
            true=te.residual_activity.to_numpy(float)
            rows.append({"held_out_dataset":held,"program_id":pid,"n_test_timepoints":len(te),"rmse":float(np.sqrt(np.mean((true-pred)**2))),"mae":float(np.mean(np.abs(true-pred))),"profile_pearson":corr(true,pred),"profile_spearman":corr(true,pred,"spearman"),"train_residual_time_slope":float(b[1])})
    return pd.DataFrame(rows)

def _cross_dataset_residual_signal(df):
    rows=[]
    for pid,g in df.groupby("program_id",sort=True):
        dsvals=[]
        for ds,h in g.groupby("dataset"):
            r=h.residual_activity.to_numpy(float); t=h.normalized_time.to_numpy(float)
            if len(r)>=3: dsvals.append((ds,r,t))
        if len(dsvals)<2: continue
        # Compare trajectories after removing each dataset's own linear time trend.
        # Correlation is computed on overlapping normalized-time coordinates by
        # interpolation, requiring at least three points per dataset.
        grid=np.linspace(0,1,9); curves=[]
        for ds,r,t in dsvals:
            order=np.argsort(t); t=t[order];r=r[order]
            if np.unique(t).size<3: continue
            curves.append(np.interp(grid,t,r))
        if len(curves)>=2:
            vals=[]
            for i in range(len(curves)):
                for j in range(i+1,len(curves)): vals.append(corr(curves[i],curves[j],"pearson"))
            rows.append({"program_id":pid,"n_datasets":len(curves),"mean_pairwise_residual_pearson":float(np.nanmean(vals)) if vals else np.nan})
    return pd.DataFrame(rows)

def _permutation(df,n=500,seed=2916):
    rng=np.random.default_rng(seed); observed=[]; null=[]
    for _,g in df.groupby(["dataset","program_id"]):
        observed.append(abs(corr(g.residual_activity,g.normalized_time,"spearman")))
    observed=float(np.nanmean([x for x in observed if np.isfinite(x)])) if any(np.isfinite(observed)) else np.nan
    for _ in range(n):
        vals=[]
        for _,g in df.groupby(["dataset","program_id"]):
            p=rng.permutation(g.residual_activity.to_numpy(float)); vals.append(abs(corr(p,g.normalized_time,"spearman")))
        null.append(float(np.nanmean([x for x in vals if np.isfinite(x)])) if any(np.isfinite(vals)) else np.nan)
    null=np.asarray(null); p=float((1+np.sum(null>=observed))/(1+np.isfinite(null).sum())) if np.isfinite(observed) else np.nan
    return pd.DataFrame([{"observed_mean_abs_residual_time_spearman":observed,"permutation_n":n,"empirical_p":p,"null_mean":float(np.nanmean(null))}]),pd.DataFrame({"permutation":np.arange(n),"mean_abs_residual_time_spearman":null})

def run(permutations=500):
    log("starting residual-state analysis; explicit time effect removed within dataset")
    df=_load(); r=_residuals(df); r.to_csv(OUT/"01_residual_program_state.csv",index=False); log(f"computed residuals for {r.dataset.nunique()} datasets and {r.program_id.nunique()} programs")
    lodo=_lodo_residual(r); lodo.to_csv(OUT/"02_residual_lodo.csv",index=False)
    if len(lodo): summ=lodo.groupby("held_out_dataset").agg(n_programs=("program_id","count"),mean_rmse=("rmse","mean"),mean_mae=("mae","mean"),mean_profile_pearson=("profile_pearson","mean"),mean_profile_spearman=("profile_spearman","mean")).reset_index()
    else: summ=pd.DataFrame(columns=["held_out_dataset","n_programs","mean_rmse","mean_mae","mean_profile_pearson","mean_profile_spearman"])
    summ.to_csv(OUT/"03_residual_lodo_summary.csv",index=False)
    sig=_cross_dataset_residual_signal(r); sig.to_csv(OUT/"04_cross_dataset_residual_signal.csv",index=False)
    perm,null=_permutation(r,permutations);perm.to_csv(OUT/"05_residual_time_permutation_summary.csv",index=False);null.to_csv(OUT/"06_residual_time_permutation_null.csv",index=False)
    overall=pd.DataFrame([{"n_datasets":int(r.dataset.nunique()),"n_programs":int(r.program_id.nunique()),"n_lodo_datasets":int(summ.held_out_dataset.nunique()) if len(summ) else 0,"mean_residual_lodo_spearman":float(summ.mean_profile_spearman.mean()) if len(summ) else np.nan,"mean_cross_dataset_residual_pearson":float(sig.mean_pairwise_residual_pearson.mean()) if len(sig) else np.nan,"residual_time_permutation_p":float(perm.iloc[0].empirical_p) if len(perm) else np.nan,"residual_state_supported":bool(len(lodo)>=3 and np.isfinite(summ.mean_profile_spearman.mean()) and summ.mean_profile_spearman.mean()>0.2 and np.isfinite(perm.iloc[0].empirical_p) and perm.iloc[0].empirical_p<0.05)}])
    overall.to_csv(OUT/"07_STAGE2_9_16_SUMMARY.csv",index=False)
    log("complete")
    print("\nStage 2.9.16 residual LODO:");print(summ.to_string(index=False));print("\nStage 2.9.16 cross-dataset residual signal:");print(sig.to_string(index=False));print("\nStage 2.9.16 gate:");print(overall.to_string(index=False));return overall

if __name__=="__main__":run()
