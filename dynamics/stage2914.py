"""Stage 2.9.14: fixed, biologically anchored program-state validation.

Programs are defined independently of the current dataset. Activities are
computed from within-sample gene ranks and validated across all datasets whose
metadata contain at least three resolved time points. No ODE/state-space model
is fit.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"stage2_9_14"
OUT.mkdir(parents=True,exist_ok=True)

PROGRAMS={
"P01_PLURIPOTENCY":{"description":"core pluripotency/reprogramming program","positive":["POU5F1","SOX2","NANOG","LIN28A","LIN28B","DPPA4","UTF1","ESRRB","KLF4","MYC"],"negative":[]},
"P02_PROLIFERATION":{"description":"cell-cycle and proliferative program","positive":["MKI67","PCNA","TOP2A","CCNB1","CCNB2","CCNE1","CDK1","CDC20","UBE2C","TYMS","MCM2","MCM3","MCM5","MCM6","MCM7"],"negative":[]},
"P03_EMT_MESENCHYMAL":{"description":"mesenchymal/EMT versus epithelial state","positive":["VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","FN1","ITGA5","COL1A1","COL1A2","CDH2"],"negative":["CDH1","EPCAM","KRT8","KRT18","KRT19"]},
"P04_STRESS_RESPONSE":{"description":"integrated stress, heat-shock and ER-stress program","positive":["DDIT3","ATF4","HSPA1A","HSPA1B","HMOX1","XBP1","JUN","FOS","DUSP1","PPP1R15A","DNAJB1"],"negative":[]},
"P05_GLYCOLYTIC_METABOLISM":{"description":"glycolytic/metabolic remodeling program","positive":["SLC2A1","HK2","PFKP","ALDOA","GAPDH","ENO1","PKM","LDHA","PGK1","TPI1","PDK1"],"negative":[]},
"P06_FGFR_PI3K_MAPK":{"description":"FGFR receptor and downstream PI3K/MAPK signalling","positive":["FGFR1","FGFR2","FGFR3","FGFR4","FRS2","PLCG1","PIK3CA","PIK3CB","AKT1","AKT2","MAPK1","MAPK3","RAF1","SOS1"],"negative":[]},
"P07_CHROMATIN_EPIGENETIC":{"description":"chromatin/remodelling and epigenetic regulation","positive":["KMT2A","KMT2B","EZH2","DNMT1","DNMT3A","DNMT3B","TET1","TET2","HDAC1","HDAC2","SMARCA4","ARID1A","CHD4","SUZ12"],"negative":[]},
"P08_ECM_ADHESION":{"description":"extracellular matrix and cell-adhesion remodeling","positive":["FN1","ITGA5","ITGB1","COL1A1","COL1A2","COL3A1","SPARC","VCAN","THBS1","LAMC1","LAMA4"],"negative":[]},
}

def log(x): print(f"Stage 2.9.14: {x}",flush=True)
def corr(a,b,method="pearson"):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method="spearman" if method=="spearman" else "pearson"))
def load_space():
    from dynamics.validation import _load_common_space
    return _load_common_space()
def _rank_matrix(matrix): return matrix.rank(axis=0,method="average",pct=True)
def activity(matrix):
    ranks=_rank_matrix(matrix);lookup={str(g).upper():g for g in matrix.index};frames=[]
    for pid,spec in PROGRAMS.items():
        pos=[lookup[g] for g in spec["positive"] if g in lookup];neg=[lookup[g] for g in spec["negative"] if g in lookup]
        if len(pos)<3: continue
        score=ranks.loc[pos].mean(axis=0)
        if neg: score=score-ranks.loc[neg].mean(axis=0)
        frames.append(pd.DataFrame({"program_id":pid,"activity":score.values},index=matrix.columns))
    return pd.concat(frames).reset_index(names="matrix_column")

def _trajectory(activity_df,meta,ds):
    g=meta[meta.dataset.astype(str)==ds].copy()
    # IMPORTANT: _load_common_space has already resolved time_hours using GEO
    # metadata, GSM_TIME and dataset-specific row-order rules.  Do not replace
    # those values by parsing sample names; many GSM labels contain no time.
    if "time_hours" not in g.columns: return None
    g["time_resolved"]=pd.to_numeric(g["time_hours"],errors="coerce")
    g=g[g.time_resolved.notna()&g.matrix_column.notna()]
    g=g[g.matrix_column.astype(str).isin(set(activity_df.matrix_column.astype(str)))]
    if g.time_resolved.nunique()<3:return None
    a=activity_df[activity_df.matrix_column.astype(str).isin(g.matrix_column.astype(str))].copy()
    a["matrix_column"]=a.matrix_column.astype(str);g["matrix_column"]=g.matrix_column.astype(str)
    a=a.merge(g[["matrix_column","time_resolved"]],on="matrix_column",how="inner")
    rec=a.groupby(["time_resolved","program_id"],as_index=False).activity.mean()
    piv=rec.pivot(index="time_resolved",columns="program_id",values="activity").sort_index()
    pids=list(piv.columns)
    return piv.index.to_numpy(float),piv.to_numpy(float),pids

def _lodo(trajectories):
    rows=[]
    for held,(tt,tv,pids) in trajectories.items():
        train={d:v for d,v in trajectories.items() if d!=held and v[2]==pids}
        if not train or len(tt)<3:continue
        for j,t in enumerate(tt):
            preds=[]
            for _,(xt,xv,_) in train.items():
                if xt.min()<=t<=xt.max():preds.append(np.array([np.interp(t,xt,xv[:,k]) for k in range(xv.shape[1])]))
            if not preds:continue
            pred=np.mean(preds,axis=0);true=tv[j]
            rows.append({"held_out_dataset":held,"time_hours":float(t),"n_training_datasets":len(preds),"rmse":float(np.sqrt(np.nanmean((true-pred)**2))),"mae":float(np.nanmean(np.abs(true-pred))),"profile_pearson":corr(true,pred),"profile_spearman":corr(true,pred,"spearman")})
    return pd.DataFrame(rows)

def _time_association(trajectories):
    rows=[]
    for ds,(t,v,pids) in trajectories.items():
        nt=(t-t.min())/(t.max()-t.min()) if t.max()>t.min() else np.full(len(t),np.nan)
        for k,pid in enumerate(pids):rows.append({"dataset":ds,"program_id":pid,"n_timepoints":len(t),"spearman_time":corr(v[:,k],nt,"spearman"),"pearson_time":corr(v[:,k],nt)})
    return pd.DataFrame(rows)

def _permutation_time_null(trajectories,n=500,seed=2914):
    rng=np.random.default_rng(seed);observed=[];null=[]
    for _,(t,v,_) in trajectories.items():
        nt=(t-t.min())/(t.max()-t.min()) if t.max()>t.min() else np.full(len(t),np.nan)
        observed.extend([abs(corr(v[:,k],nt,"spearman")) for k in range(v.shape[1]) if np.isfinite(corr(v[:,k],nt,"spearman"))])
    obs=float(np.nanmean(observed)) if observed else np.nan
    for b in range(n):
        vals=[]
        for _,(t,v,_) in trajectories.items():
            nt=(t-t.min())/(t.max()-t.min()) if t.max()>t.min() else np.full(len(t),np.nan);perm=rng.permutation(nt)
            vals.extend([abs(corr(v[:,k],perm,"spearman")) for k in range(v.shape[1]) if np.isfinite(corr(v[:,k],perm,"spearman"))])
        null.append(float(np.nanmean(vals)) if vals else np.nan)
    null=np.asarray(null,float);p=float((1+np.sum(null>=obs))/(1+np.isfinite(null).sum())) if np.isfinite(obs) else np.nan
    return pd.DataFrame({"observed_mean_abs_spearman":[obs],"permutation_n":[n],"empirical_p":[p],"null_mean":[np.nanmean(null)]}),pd.DataFrame({"permutation":np.arange(n),"mean_abs_spearman":null})

def run(permutations=500):
    log("starting fixed biological program-state validation; no ODE/state-space")
    matrix,meta=load_space();log(f"loaded {len(matrix):,} common genes x {matrix.shape[1]} samples")
    lookup={str(g).upper() for g in matrix.index};defs=[]
    for pid,spec in PROGRAMS.items():
        pos=[g for g in spec["positive"] if g in lookup];neg=[g for g in spec["negative"] if g in lookup]
        defs.append({"program_id":pid,"description":spec["description"],"n_positive_available":len(pos),"n_positive_defined":len(spec["positive"]),"n_negative_available":len(neg),"n_negative_defined":len(spec["negative"]),"positive_genes":";".join(pos),"negative_genes":";".join(neg),"usable":len(pos)>=3})
    pd.DataFrame(defs).to_csv(OUT/"01_program_definitions.csv",index=False)
    act=activity(matrix);act.to_csv(OUT/"02_program_activity_by_sample.csv",index=False)
    trajectories={}
    for ds in sorted(meta.dataset.astype(str).unique()):
        tr=_trajectory(act,meta,ds)
        if tr is not None:trajectories[ds]=tr
    log(f"constructed trajectories for {len(trajectories)} datasets: {', '.join(sorted(trajectories)) if trajectories else 'none'}")
    traj_rows=[]
    for ds,(t,v,pids) in trajectories.items():
        for i,tt in enumerate(t):
            for k,pid in enumerate(pids):traj_rows.append({"dataset":ds,"time_hours":float(tt),"program_id":pid,"activity":float(v[i,k])})
    pd.DataFrame(traj_rows).to_csv(OUT/"03_program_trajectories.csv",index=False)
    assoc=_time_association(trajectories);assoc.to_csv(OUT/"04_program_time_association.csv",index=False)
    lodo=_lodo(trajectories);lodo.to_csv(OUT/"05_lodo_by_timepoint.csv",index=False)
    if len(lodo):summary=lodo.groupby("held_out_dataset").agg(n_timepoints=("time_hours","count"),mean_rmse=("rmse","mean"),mean_mae=("mae","mean"),mean_profile_pearson=("profile_pearson","mean"),mean_profile_spearman=("profile_spearman","mean")).reset_index()
    else:summary=pd.DataFrame(columns=["held_out_dataset","n_timepoints","mean_rmse","mean_mae","mean_profile_pearson","mean_profile_spearman"])
    summary.to_csv(OUT/"06_lodo_summary.csv",index=False)
    main=pd.DataFrame([{"dataset":"GSE67462","program_id":pid,"time_hours":t.tolist(),"activity":v[:,k].tolist()} for t,v,pids in [trajectories["GSE67462"]] for k,pid in enumerate(pids)]) if "GSE67462" in trajectories else pd.DataFrame()
    main.to_csv(OUT/"07_GSE67462_program_trajectories.csv",index=False)
    perm,null=_permutation_time_null(trajectories,permutations);perm.to_csv(OUT/"08_time_permutation_summary.csv",index=False);null.to_csv(OUT/"09_time_permutation_null.csv",index=False)
    overall=pd.DataFrame([{"n_programs_defined":len(PROGRAMS),"n_programs_usable":int(sum(x["usable"] for x in defs)),"n_datasets_with_trajectory":len(trajectories),"n_lodo_datasets":len(summary),"mean_lodo_rmse":float(summary.mean_rmse.mean()) if len(summary) else np.nan,"mean_lodo_profile_spearman":float(summary.mean_profile_spearman.mean()) if len(summary) else np.nan,"GSE67462_present":bool("GSE67462" in trajectories),"time_permutation_empirical_p":float(perm.iloc[0].empirical_p) if len(perm) else np.nan}])
    overall.to_csv(OUT/"10_overall_summary.csv",index=False)
    log("complete")
    print("\nStage 2.9.14 program time association:");print(assoc.to_string(index=False));print("\nStage 2.9.14 LODO:");print(summary.to_string(index=False));print("\nStage 2.9.14 overall:");print(overall.to_string(index=False));return overall

if __name__=="__main__":run()
