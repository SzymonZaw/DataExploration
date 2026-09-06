"""Stage 2.9.2: dataset-invariant latent reprogramming progress validation."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from .stage291 import _load_common_space, _trajectories, _profile_cube, _stability, _select, _fit_programs, _state, _state_traj

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_2"
OUT.mkdir(parents=True, exist_ok=True)


def _program_grid(traj, grid=17):
    if traj is None or len(traj.index) < 2: return None
    t=traj.index.to_numpy(float); y=traj.to_numpy(float)
    if t[-1] <= t[0]: return None
    u=np.linspace(0,1,grid); tn=(t-t[0])/(t[-1]-t[0]); out=np.full((grid,y.shape[1]),np.nan)
    for j in range(y.shape[1]):
        ok=np.isfinite(y[:,j])
        if ok.sum()>=2: out[:,j]=np.interp(u,tn[ok],y[ok,j])
    return out


def _fit_progress(train_trajs, seed=42, grid=17):
    cubes=[]; used=[]
    for ds in sorted(train_trajs):
        g=_program_grid(train_trajs[ds],grid)
        if g is not None and np.isfinite(g).mean()>=0.7: cubes.append(g); used.append(ds)
    if len(cubes)<2:return None
    X=np.concatenate(cubes,axis=0); med=np.nanmedian(X,axis=0); X=np.where(np.isfinite(X),X,med)
    scaler=StandardScaler().fit(X); Xs=scaler.transform(X); pca=PCA(n_components=min(3,Xs.shape[1]),random_state=seed).fit(Xs)
    scores=pca.transform(Xs)[:,0]; target=np.tile(np.linspace(0,1,grid),len(used)); corr=np.corrcoef(scores,target)[0,1] if np.std(scores)>1e-12 else np.nan
    sign=1.0 if np.isfinite(corr) and corr>=0 else -1.0
    return {"scaler":scaler,"pca":pca,"sign":sign,"training_datasets":used,"train_pc1_time_correlation":float(corr)}


def _fit_progress_with_scores(train_trajs, seed=42, grid=17):
    model=_fit_progress(train_trajs,seed,grid)
    if model is None:return None
    X=np.concatenate([_program_grid(train_trajs[ds],grid) for ds in model["training_datasets"]],axis=0); med=np.nanmedian(X,axis=0); X=np.where(np.isfinite(X),X,med)
    model["training_scores"]=model["sign"]*model["pca"].transform(model["scaler"].transform(X))[:,0]
    return model


def _project_progress(traj,model,grid=17):
    g=_program_grid(traj,grid)
    if g is None:return None
    med=np.nanmedian(g,axis=0); g=np.where(np.isfinite(g),g,med); z=model["sign"]*model["pca"].transform(model["scaler"].transform(g))[:,0]
    ref=model["training_scores"]; lo,hi=np.nanpercentile(ref,[2.5,97.5])
    if hi<=lo: lo,hi=np.nanmin(ref),np.nanmax(ref)
    if hi<=lo:return np.full(len(z),0.5)
    return np.clip((z-lo)/(hi-lo),0,1)


def _fold(matrix,meta,held_out,n_programs=10,max_genes=2500,seed=42,grid=17):
    datasets=sorted(meta.dataset.astype(str).unique()); train=[d for d in datasets if d!=held_out]; tr=_trajectories(matrix,meta,train); te=_trajectories(matrix,meta,[held_out])
    if held_out not in te or len(tr)<2:return None,{"reason":"insufficient_trajectories","n_train_trajectories":len(tr)}
    cube,used=_profile_cube(tr); stab=_stability(cube); ids,selected,min_ds=_select(stab,max_genes)
    if len(ids)<max(2,n_programs):return None,{"reason":"too_few_stable_genes","n_train_trajectories":len(used),"n_stable_genes":len(ids),"min_datasets":min_ds}
    programs=_fit_programs(cube,ids,n_programs,seed); state=_state(matrix,meta,programs,train); trajs=_state_traj(state,meta); train_trajs={d:trajs[d] for d in train if d in trajs}; test=trajs.get(held_out); model=_fit_progress_with_scores(train_trajs,seed,grid)
    if model is None or test is None:return None,{"reason":"no_progress_model","n_train_trajectories":len(used),"n_stable_genes":len(ids)}
    pred=_project_progress(test,model,grid); actual=(test.index.to_numpy(float)-test.index.min())/(test.index.max()-test.index.min())
    if len(pred)<3:return None,{"reason":"too_few_test_timepoints","n_test_timepoints":len(pred)}
    spearman=float(pd.Series(actual).corr(pd.Series(pred),method="spearman")); pear=float(np.corrcoef(actual,pred)[0,1]) if np.std(pred)>1e-12 else np.nan; rmse=float(np.sqrt(np.mean((actual-pred)**2))); mono=float(np.mean(np.diff(pred)>=-0.02))
    row={"held_out_dataset":held_out,"n_timepoints":len(pred),"n_training_trajectories":len(used),"n_stable_genes":len(ids),"n_programs":len(programs),"min_datasets_for_selection":min_ds,"mean_gene_stability":float(selected.median_pairwise_correlation.head(len(ids)).mean()),"progress_spearman":spearman,"progress_pearson":pear,"progress_rmse":rmse,"progress_monotonicity":mono,"train_pc1_time_correlation":model["train_pc1_time_correlation"]}
    path=pd.DataFrame({"held_out_dataset":held_out,"time_hours":test.index.to_numpy(float),"normalized_time":actual,"latent_progress":pred})
    return {"summary":pd.DataFrame([row]),"trajectory":path},{"reason":"ok","n_train_trajectories":len(used),"n_stable_genes":len(ids),"min_datasets":min_ds}


def _permute_meta(meta,seed):
    rng=np.random.default_rng(seed); out=meta.copy(); out["time_hours"]=pd.to_numeric(out["time_hours"],errors="coerce")
    for ds,g in out.groupby("dataset"):
        idx=g.index[g.time_hours.notna()&g.matrix_column.notna()]; vals=out.loc[idx,"time_hours"].to_numpy(float)
        if len(vals)>=3: out.loc[idx,"time_hours"]=rng.permutation(vals)
    return out


def stage2_9_2(n_programs=10,max_genes=2500,permutations=20,seed=42,grid=17):
    print("Stage 2.9.2: loading common gene space...",flush=True); matrix,meta=_load_common_space(); datasets=sorted(meta.dataset.astype(str).unique()); observed=[]; trajectories=[]; statuses=[]
    for i,ds in enumerate(datasets,1):
        print(f"Stage 2.9.2: fold {i}/{len(datasets)}, held out={ds}...",flush=True); result,status=_fold(matrix,meta,ds,n_programs,max_genes,seed,grid); status.update({"held_out_dataset":ds,"mode":"observed"}); statuses.append(status)
        if result is None: print(f"  skipped: {status}",flush=True)
        else: observed.append(result["summary"]); trajectories.append(result["trajectory"])
    summary=pd.concat(observed,ignore_index=True) if observed else pd.DataFrame(); summary.to_csv(OUT/"01_latent_progress_summary.csv",index=False); (pd.concat(trajectories,ignore_index=True) if trajectories else pd.DataFrame()).to_csv(OUT/"02_latent_progress_trajectories.csv",index=False); pd.DataFrame(statuses).to_csv(OUT/"05_fold_status.csv",index=False)
    null=[]
    for p in range(permutations):
        pm=_permute_meta(meta,seed+1000+p); vals=[]
        for ds in datasets:
            r,_=_fold(matrix,pm,ds,n_programs,max_genes,seed+1000+p,grid)
            if r is not None: vals.append(r["summary"])
        if vals:
            q=pd.concat(vals,ignore_index=True); null.append({"permutation":p,"mean_progress_spearman":float(q.progress_spearman.mean()),"mean_progress_pearson":float(q.progress_pearson.mean()),"mean_progress_rmse":float(q.progress_rmse.mean()),"mean_monotonicity":float(q.progress_monotonicity.mean()),"n_cases":len(q)})
        print(f"Stage 2.9.2: permutation {p+1}/{permutations}",flush=True)
    null=pd.DataFrame(null); null.to_csv(OUT/"03_permutation_null.csv",index=False)
    if summary.empty: print("Stage 2.9.2: no valid folds; inspect 05_fold_status.csv.")
    else:
        rho=float(summary.progress_spearman.mean()); pear=float(summary.progress_pearson.mean()); rmse=float(summary.progress_rmse.mean()); mono=float(summary.progress_monotonicity.mean()); print(summary.to_string(index=False)); print(f"Observed mean latent-progress Spearman: {rho:.3f}; Pearson: {pear:.3f}; RMSE: {rmse:.3f}; monotonicity: {mono:.3f}")
        if not null.empty:
            p_rho=(1+int((null.mean_progress_spearman>=rho).sum()))/(len(null)+1); p_pear=(1+int((null.mean_progress_pearson>=pear).sum()))/(len(null)+1); print(f"Permutation empirical p: Spearman={p_rho:.3f}, Pearson={p_pear:.3f}")
        print("Interpret latent progress as a validation target, not yet a mechanistic state. Proceed to ODE/state-space only after cross-dataset progress is reproducible and null-resistant.")
    return summary


if __name__=="__main__":stage2_9_2()
