"""Stage 2.9.5: bootstrap stability of the dataset-invariant latent axis.

Uses full program-state trajectories persisted by Stage 2.9.2. For each
leave-one-dataset-out fold, PCA is refit after bootstrap resampling of program
dimensions. The held-out trajectory is projected with each axis. Reports axis
agreement, transfer stability, orientation, and individual-program sensitivity.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
IN=ROOT/"results"/"Dynamics"/"stage2_9_2"
OUT=ROOT/"results"/"Dynamics"/"stage2_9_5"
OUT.mkdir(parents=True,exist_ok=True)

def corr(a,b,method="spearman"):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method=method))

def fit_axis(X):
    X=np.asarray(X,float);med=np.nanmedian(X,axis=0);X=np.where(np.isfinite(X),X,med);X=np.where(np.isfinite(X),X,0.0)
    if X.shape[0]<3 or X.shape[1]<2:return None,None
    sc=StandardScaler().fit(X);pc=PCA(n_components=1,random_state=42).fit(sc.transform(X));return sc,pc

def project(X,sc,pc):
    X=np.asarray(X,float);med=np.nanmedian(X,axis=0);X=np.where(np.isfinite(X),X,med);X=np.where(np.isfinite(X),X,0.0);return pc.transform(sc.transform(X))[:,0]

def normalize(z,ref):
    lo,hi=np.nanpercentile(ref,[2.5,97.5])
    if hi<=lo:lo,hi=np.nanmin(ref),np.nanmax(ref)
    if not np.isfinite(lo) or not np.isfinite(hi) or hi<=lo:return np.full(len(z),.5)
    return np.clip((z-lo)/(hi-lo),0,1)

def run(bootstrap_replicates=500,seed=42):
    state_path=IN/"04_program_state_trajectories.csv";traj_path=IN/"02_latent_progress_trajectories.csv"
    if not state_path.exists() or not traj_path.exists():raise RuntimeError("Run Stage 2.9.2 first: python validate_pipeline.py --stage292")
    state=pd.read_csv(state_path);reported=pd.read_csv(traj_path);rng=np.random.default_rng(seed)
    programs=sorted([c for c in state.columns if c.startswith("program_")])
    if len(programs)<2:raise RuntimeError("Stage 2.9.2 did not persist enough program-state columns.")
    boots=[];folds=[];deletions=[]
    for ds in sorted(state.held_out_dataset.astype(str).unique()):
        fold=state[state.held_out_dataset.astype(str)==ds];train=fold[fold.dataset.astype(str)!=ds].sort_values(["dataset","time_hours"]);test=fold[fold.dataset.astype(str)==ds].sort_values("time_hours")
        rep=reported[reported.held_out_dataset.astype(str)==ds].sort_values("time_hours");rz=rep.latent_progress.to_numpy(float);rt=rep.normalized_time.to_numpy(float)
        sc,pc=fit_axis(train[programs].to_numpy(float))
        if pc is None or len(test)<3:continue
        base_tr=project(train[programs],sc,pc);base_te=project(test[programs],sc,pc);base=normalize(base_te,base_tr);base_loading=pc.components_[0]
        base_r=corr(rz,base,"pearson")
        for b in range(bootstrap_replicates):
            idx=rng.integers(0,len(programs),len(programs));scb,pcb=fit_axis(train[programs].to_numpy(float)[:,idx])
            if pcb is None:continue
            load=pcb.components_[0];sim=float(np.dot(load,base_loading)/(np.linalg.norm(load)*np.linalg.norm(base_loading))) if np.linalg.norm(load)>0 else np.nan
            zbtr=project(train[programs].to_numpy(float)[:,idx],scb,pcb);zbte=project(test[programs].to_numpy(float)[:,idx],scb,pcb)
            if np.isfinite(sim) and sim<0:sim=-sim;zbtr=-zbtr;zbte=-zbte
            zn=normalize(zbte,zbtr)
            boots.append({"held_out_dataset":ds,"bootstrap":b,"axis_cosine_similarity":sim,"reported_progress_spearman":corr(rz,zn),"reported_progress_pearson":corr(rz,zn,"pearson"),"time_spearman":corr(rt,zn)})
        for col in programs:
            keep=[c for c in programs if c!=col];scd,pcd=fit_axis(train[keep].to_numpy(float))
            if pcd is None:continue
            zdtr=project(train[keep],scd,pcd);zdte=project(test[keep],scd,pcd);zn=normalize(zdte,zdtr)
            deletions.append({"held_out_dataset":ds,"excluded_program":col,"reported_progress_spearman":corr(rz,zn),"reported_progress_pearson":corr(rz,zn,"pearson")})
        bg=pd.DataFrame([x for x in boots if x["held_out_dataset"]==ds]);folds.append({"held_out_dataset":ds,"n_training_rows":len(train),"n_test_timepoints":len(test),"n_programs":len(programs),"recomputed_pca1_vs_reported_progress_pearson":base_r,"bootstrap_axis_cosine_mean":bg.axis_cosine_similarity.mean(),"bootstrap_axis_cosine_p05":bg.axis_cosine_similarity.quantile(.05),"bootstrap_axis_cosine_p95":bg.axis_cosine_similarity.quantile(.95),"bootstrap_progress_spearman_mean":bg.reported_progress_spearman.mean(),"bootstrap_progress_spearman_p05":bg.reported_progress_spearman.quantile(.05),"bootstrap_progress_spearman_p95":bg.reported_progress_spearman.quantile(.95),"bootstrap_progress_ci_excludes_zero":bool(bg.reported_progress_spearman.quantile(.05)>0)})
        print(f"Stage 2.9.5: fold {ds} ({bootstrap_replicates} program bootstraps)",flush=True)
    boot=pd.DataFrame(boots);fold=pd.DataFrame(folds);delete=pd.DataFrame(deletions)
    boot.to_csv(OUT/"02_program_bootstrap.csv",index=False);fold.to_csv(OUT/"03_fold_axis_stability.csv",index=False);delete.to_csv(OUT/"04_program_deletion_sensitivity.csv",index=False)
    summary=pd.DataFrame([{"n_datasets":len(fold),"mean_axis_cosine":fold.bootstrap_axis_cosine_mean.mean() if not fold.empty else np.nan,"min_axis_cosine_p05":fold.bootstrap_axis_cosine_p05.min() if not fold.empty else np.nan,"n_folds_progress_ci_excludes_zero":int(fold.bootstrap_progress_ci_excludes_zero.sum()) if not fold.empty else 0,"mean_bootstrap_progress_spearman":fold.bootstrap_progress_spearman_mean.mean() if not fold.empty else np.nan,"mean_recomputed_pca1_vs_reported_progress_pearson":fold.recomputed_pca1_vs_reported_progress_pearson.mean() if not fold.empty else np.nan,"interpretation":"program-dimension bootstrap; gene/program discovery itself is not fully refit in every bootstrap"}])
    summary.to_csv(OUT/"01_overall_summary.csv",index=False)
    print("Stage 2.9.5: fold results",flush=True);print(fold.to_string(index=False),flush=True);print("\nStage 2.9.5 summary",flush=True);print(summary.to_string(index=False),flush=True);print("Stage 2.9.5 complete. This is a robustness gate; do not fit ODE/state-space yet if axis stability remains weak.",flush=True);return summary

if __name__=="__main__":run()
