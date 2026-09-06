"""Stage 2.9.15: pre-Stage-3 readiness gate for biological program state.

Adds the missing decision tests to Stage 2.9.14:
- nearest-time baseline for LODO prediction;
- program-state improvement over that baseline;
- bootstrap stability of marker-panel activities;
- time-label permutation test;
- explicit STAGE3_READINESS gate.

This stage does not fit an ODE/state-space model. A positive gate requires
predictive information beyond nearest time and stable program activities.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"stage2_9_15"
OUT.mkdir(parents=True,exist_ok=True)

from dynamics.stage2914 import PROGRAMS, load_space, activity


def log(x): print(f"Stage 2.9.15: {x}",flush=True)


def corr(a,b,method="spearman"):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method=method))


def trajectories(act,meta):
    out={}
    for ds,g in meta.groupby("dataset",sort=True):
        g=g[g.matrix_column.notna() & g.time_hours.notna()].copy()
        if g.time_hours.nunique()<3: continue
        g["matrix_column"]=g.matrix_column.astype(str)
        a=act[act.matrix_column.astype(str).isin(set(g.matrix_column))].copy()
        a["matrix_column"]=a.matrix_column.astype(str)
        a=a.merge(g[["matrix_column","time_hours"]],on="matrix_column",how="inner")
        rec=a.groupby(["time_hours","program_id"],as_index=False).activity.mean()
        piv=rec.pivot(index="time_hours",columns="program_id",values="activity").sort_index()
        if piv.shape[0]>=3 and piv.notna().all().all():
            out[ds]=(piv.index.to_numpy(float),piv.to_numpy(float),list(piv.columns))
    return out


def _program_lodo(traj):
    rows=[]
    for held,(tt,tv,pids) in traj.items():
        train={d:x for d,x in traj.items() if d!=held and x[2]==pids}
        if not train: continue
        for j,t in enumerate(tt):
            preds=[]; nearest=[]
            for _,(xt,xv,_) in train.items():
                if xt.min()<=t<=xt.max():
                    preds.append(np.array([np.interp(t,xt,xv[:,k]) for k in range(xv.shape[1])]))
                # nearest-time baseline: state at the closest observed training time.
                idx=int(np.argmin(np.abs(xt-t)))
                nearest.append(xv[idx])
            if not preds: continue
            true=tv[j]; pred=np.mean(preds,axis=0); base=np.mean(nearest,axis=0)
            rmse=float(np.sqrt(np.mean((true-pred)**2))); base_rmse=float(np.sqrt(np.mean((true-base)**2)))
            rows.append({"held_out_dataset":held,"time_hours":float(t),"n_training_datasets":len(preds),"program_rmse":rmse,"nearest_time_rmse":base_rmse,"rmse_improvement":base_rmse-rmse,"program_profile_spearman":corr(true,pred),"nearest_profile_spearman":corr(true,base)})
    return pd.DataFrame(rows)


def _bootstrap_markers(matrix,meta,n=200,seed=2915):
    rng=np.random.default_rng(seed); base=activity(matrix); rows=[]
    lookup={str(g).upper():g for g in matrix.index}
    for b in range(n):
        ranks=matrix.rank(axis=0,method="average",pct=True)
        frames=[]
        for pid,spec in PROGRAMS.items():
            pos=[lookup[g] for g in spec["positive"] if g in lookup]
            neg=[lookup[g] for g in spec["negative"] if g in lookup]
            if len(pos)<3: continue
            k=max(3,int(round(0.8*len(pos))))
            sel=list(rng.choice(pos,size=min(k,len(pos)),replace=False))
            score=ranks.loc[sel].mean(axis=0)
            if neg:
                nk=max(2,int(round(0.8*len(neg))))
                nsel=list(rng.choice(neg,size=min(nk,len(neg)),replace=False))
                score=score-ranks.loc[nsel].mean(axis=0)
            frames.append(pd.DataFrame({"program_id":pid,"activity":score.values},index=matrix.columns))
        boot=pd.concat(frames).reset_index(names="matrix_column")
        for ds,(t,v,pids) in trajectories(boot,meta).items():
            nt=(t-t.min())/(t.max()-t.min())
            vals=[abs(corr(v[:,k],nt)) for k in range(v.shape[1])]
            rows.append({"bootstrap":b,"dataset":ds,"mean_abs_time_spearman":float(np.nanmean(vals))})
        if (b+1)%50==0: log(f"marker bootstrap {b+1}/{n}")
    return pd.DataFrame(rows)


def _permutation(traj,n=500,seed=2915):
    rng=np.random.default_rng(seed); obs=[]; null=[]
    for _,(t,v,_) in traj.items():
        nt=(t-t.min())/(t.max()-t.min())
        obs.extend([abs(corr(v[:,k],nt)) for k in range(v.shape[1])])
    observed=float(np.nanmean(obs))
    for b in range(n):
        vals=[]
        for _,(t,v,_) in traj.items():
            nt=(t-t.min())/(t.max()-t.min())
            p=rng.permutation(nt)
            vals.extend([abs(corr(v[:,k],p)) for k in range(v.shape[1])])
        null.append(float(np.nanmean(vals)))
        if (b+1)%100==0: log(f"time permutation {b+1}/{n}")
    null=np.asarray(null); p=float((1+np.sum(null>=observed))/(1+n))
    return observed,p,null


def run(bootstrap_n=200,permutations=500):
    log("starting pre-Stage-3 readiness gate; no ODE/state-space")
    matrix,meta=load_space(); act=activity(matrix); traj=trajectories(act,meta)
    log(f"usable trajectories: {', '.join(traj) if traj else 'none'}")
    lodo=_program_lodo(traj); lodo.to_csv(OUT/"01_lodo_vs_nearest_time.csv",index=False)
    if len(lodo):
        s=lodo.groupby("held_out_dataset").agg(n_timepoints=("time_hours","count"),program_rmse=("program_rmse","mean"),nearest_time_rmse=("nearest_time_rmse","mean"),rmse_improvement=("rmse_improvement","mean"),program_profile_spearman=("program_profile_spearman","mean"),nearest_profile_spearman=("nearest_profile_spearman","mean")).reset_index()
    else: s=pd.DataFrame()
    s.to_csv(OUT/"02_lodo_summary.csv",index=False)
    boot=_bootstrap_markers(matrix,meta,bootstrap_n);boot.to_csv(OUT/"03_marker_bootstrap.csv",index=False)
    if len(boot):
        bs=boot.groupby("dataset").agg(n_bootstrap=("bootstrap","nunique"),mean_abs_time_spearman=("mean_abs_time_spearman","mean"),p05=("mean_abs_time_spearman",lambda x:x.quantile(.05)),p95=("mean_abs_time_spearman",lambda x:x.quantile(.95))).reset_index()
    else: bs=pd.DataFrame()
    bs.to_csv(OUT/"04_marker_bootstrap_summary.csv",index=False)
    observed,p,null=_permutation(traj,permutations);pd.DataFrame([{"observed_mean_abs_spearman":observed,"permutation_n":permutations,"empirical_p":p,"null_mean":float(np.mean(null))}]).to_csv(OUT/"05_time_permutation_summary.csv",index=False)
    pd.DataFrame({"permutation":np.arange(len(null)),"mean_abs_spearman":null}).to_csv(OUT/"06_time_permutation_null.csv",index=False)
    n_lodo=len(s); positive_improvement=int((s.rmse_improvement>0).sum()) if len(s) else 0
    stable_boot=int((bs.p05>0.25).sum()) if len(bs) else 0
    overall=pd.DataFrame([{
        "n_programs":len(PROGRAMS),"n_trajectory_datasets":len(traj),"n_lodo_datasets":n_lodo,
        "n_lodo_better_than_nearest_time":positive_improvement,
        "mean_rmse_improvement":float(s.rmse_improvement.mean()) if len(s) else np.nan,
        "mean_lodo_profile_spearman":float(s.program_profile_spearman.mean()) if len(s) else np.nan,
        "n_marker_bootstrap_datasets_with_p05_gt_0_25":stable_boot,
        "mean_bootstrap_abs_time_spearman":float(bs.mean_abs_time_spearman.mean()) if len(bs) else np.nan,
        "time_permutation_p":p,
        "STAGE3_READINESS":bool(n_lodo>=2 and positive_improvement>=max(1,n_lodo//2+1) and stable_boot>=2 and p<0.05)
    }])
    overall.to_csv(OUT/"07_STAGE3_READINESS.csv",index=False)
    log("complete")
    print("\nStage 2.9.15 LODO vs nearest-time:");print(s.to_string(index=False))
    print("\nStage 2.9.15 marker-bootstrap summary:");print(bs.to_string(index=False))
    print("\nStage 2.9.15 STAGE 3 GATE:");print(overall.to_string(index=False))
    return overall

if __name__=="__main__": run()
