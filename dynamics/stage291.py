"""Stage 2.9.1: robust leakage-free invariant temporal-program validation."""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_1"
OUT.mkdir(parents=True, exist_ok=True)


def _load_common_space():
    from .validation import _load_common_space as load
    return load()


def _metrics(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 2: return {"rmse": np.nan, "mae": np.nan, "correlation": np.nan}
    a, b = a[ok], b[ok]
    corr = np.corrcoef(a, b)[0, 1] if np.std(a) > 1e-12 and np.std(b) > 1e-12 else np.nan
    return {"rmse": float(np.sqrt(np.mean((a-b)**2))), "mae": float(np.mean(np.abs(a-b))), "correlation": float(corr)}


def _times(meta, override=None):
    out = {}
    for _, r in meta.iterrows():
        if pd.isna(r.matrix_column): continue
        c = str(r.matrix_column); ds = str(r.dataset)
        t = override.get((ds, c), r.time_hours) if override is not None else r.time_hours
        if pd.notna(t): out[(ds, c)] = float(t)
    return out


def _trajectories(matrix, meta, datasets, override=None):
    tm = _times(meta, override); out = {}
    for ds in datasets:
        g = meta[(meta.dataset.astype(str) == str(ds)) & meta.matrix_column.notna()].copy()
        g["_time"] = [tm.get((str(ds), str(c)), np.nan) for c in g.matrix_column]
        g = g[np.isfinite(g["_time"].to_numpy(float))]
        if len(g) < 3 or g["_time"].nunique() < 3: continue
        cols = g.matrix_column.astype(str).tolist()
        frame = matrix[cols].T.copy(); frame.index = g["_time"].to_numpy(float)
        frame = frame.groupby(level=0).mean().sort_index()
        if len(frame) >= 3:
            t = frame.index.to_numpy(float); span = t[-1]-t[0]
            if span > 0: out[str(ds)] = ((t-t[0])/span, frame.to_numpy(float))
    return out


def _profile_cube(trajs, grid=9):
    if not trajs: return None, []
    n_genes = next(iter(trajs.values()))[1].shape[1]; gx = np.linspace(0,1,grid); datasets=list(trajs)
    cube=np.full((n_genes,len(datasets),grid),np.nan)
    for di,ds in enumerate(datasets):
        t,y=trajs[ds]
        for j in range(n_genes):
            v=y[:,j]; ok=np.isfinite(v)
            if ok.sum()>=3: cube[j,di]=np.interp(gx,t[ok],v[ok])
    return cube,datasets


def _stability(cube):
    rows=[]
    for j in range(cube.shape[0]):
        p=cube[j]; valid=[i for i in range(p.shape[0]) if np.isfinite(p[i]).sum()>=3]; cs=[]
        for ai in range(len(valid)):
            for bi in range(ai+1,len(valid)):
                x,y=p[valid[ai]],p[valid[bi]]
                if np.std(x)>1e-10 and np.std(y)>1e-10: cs.append(np.corrcoef(x,y)[0,1])
        rows.append({"gene_index":j,"n_datasets":len(valid),"n_pairs":len(cs),"median_pairwise_correlation":float(np.median(cs)) if cs else np.nan,"mean_pairwise_correlation":float(np.mean(cs)) if cs else np.nan})
    return pd.DataFrame(rows)


def _select(stability,max_genes=2500):
    s=stability[np.isfinite(stability.median_pairwise_correlation)].copy()
    if s.empty:return np.array([],dtype=int),s,0
    # Prefer three independent trajectories; fall back to two only when the
    # training fold cannot supply enough genes for a stable program model.
    for min_ds in (3,2):
        q=s[s.n_datasets>=min_ds].copy()
        if len(q)<20: continue
        threshold=float(q.median_pairwise_correlation.quantile(.75)); q=q[q.median_pairwise_correlation>=threshold]
        q=q.sort_values(["median_pairwise_correlation","n_datasets"],ascending=False).head(max_genes)
        if len(q)>=2:return q.gene_index.astype(int).to_numpy(),q,min_ds
    return np.array([],dtype=int),s.sort_values("median_pairwise_correlation",ascending=False),0


def _fit_programs(cube,ids,n_programs=10,seed=42):
    if len(ids)<max(2,n_programs):raise RuntimeError(f"Too few stable genes ({len(ids)}) for {n_programs} programs.")
    X=cube[ids].reshape(len(ids),-1); med=np.nanmedian(X,axis=0); X=np.where(np.isfinite(X),X,med); X=StandardScaler().fit_transform(X)
    labels=KMeans(n_clusters=min(n_programs,len(ids)),random_state=seed,n_init=20).fit_predict(X)
    return {i:ids[labels==i].tolist() for i in np.unique(labels)}


def _state(matrix,meta,programs,train):
    fit=meta[meta.dataset.astype(str).isin([str(x) for x in train])&meta.matrix_column.notna()&meta.time_hours.notna()]
    state=pd.DataFrame(index=[str(c) for c in matrix.columns]); fit_cols=fit.matrix_column.astype(str).tolist()
    for pi,ids in programs.items():
        genes=[matrix.index[int(i)] for i in ids if 0<=int(i)<len(matrix.index)]
        if not genes or not fit_cols:continue
        vals=matrix.loc[genes,fit_cols].to_numpy(float); mu=np.nanmedian(vals,axis=1); mad=1.4826*np.nanmedian(np.abs(vals-mu[:,None]),axis=1); mad[~np.isfinite(mad)|(mad<1e-8)]=1
        v=matrix.loc[genes,state.index].to_numpy(float); z=(np.where(np.isfinite(v),v,mu[:,None])-mu[:,None])/mad[:,None]; state[f"program_{pi:02d}"]=np.nanmedian(z,axis=0)
    return state


def _state_traj(state,meta,override=None):
    tm=_times(meta,override); out={}
    for ds,g in meta.groupby("dataset"):
        g=g[g.matrix_column.notna()].copy(); g["_time"]=[tm.get((str(ds),str(c)),np.nan) for c in g.matrix_column]; g=g[np.isfinite(g["_time"].to_numpy(float))]
        if len(g)<2 or g["_time"].nunique()<2:continue
        f=state.loc[g.matrix_column.astype(str)].copy();f.index=g["_time"].astype(float).to_numpy();f=f.groupby(level=0).mean().sort_index()
        if len(f)>=2:out[str(ds)]=f
    return out


def _interp(f,t):
    if f is None or t<f.index.min() or t>f.index.max():return None
    x=f.index.to_numpy(float);y=f.to_numpy(float);return np.asarray([np.interp(t,x,y[:,j]) for j in range(y.shape[1])])


def _fold(matrix,meta,held_out,n_programs,max_genes,seed,override=None):
    all_ds=sorted(meta.dataset.astype(str).unique());train=[d for d in all_ds if d!=held_out];tr=_trajectories(matrix,meta,train,override);te=_trajectories(matrix,meta,[held_out],override)
    if held_out not in te or len(tr)<2:return None,None,None,{"reason":"insufficient_trajectories","n_train_trajectories":len(tr)}
    cube,used=_profile_cube(tr);stab=_stability(cube);ids,selected,min_ds=_select(stab,max_genes)
    if len(ids)<max(2,n_programs):return None,stab,None,{"reason":"too_few_stable_genes","n_train_trajectories":len(used),"n_stable_genes":len(ids),"min_datasets":min_ds}
    programs=_fit_programs(cube,ids,n_programs,seed);state=_state(matrix,meta,programs,train);traj=_state_traj(state,meta,override);test=traj.get(str(held_out))
    if test is None:return None,stab,None,{"reason":"no_test_trajectory","n_train_trajectories":len(used),"n_stable_genes":len(ids),"min_datasets":min_ds}
    rows=[]
    for t in test.index:
        preds=[_interp(traj.get(d),float(t)) for d in train if _interp(traj.get(d),float(t)) is not None]
        if not preds:continue
        pred=np.mean(preds,axis=0);actual=test.loc[t].to_numpy(float);near=[f.iloc[int(np.argmin(np.abs(f.index.to_numpy(float)-t)))].to_numpy(float) for d in train if (f:=traj.get(d)) is not None];base=np.mean(near,axis=0) if near else np.full_like(pred,np.nan);m=_metrics(actual,pred);b=_metrics(actual,base)
        rows.append({"held_out_dataset":held_out,"time_hours":float(t),"n_programs":len(programs),"n_genes_selected":len(ids),"n_training_trajectories":len(used),"min_datasets_for_selection":min_ds,"mean_gene_stability":float(selected.median_pairwise_correlation.head(len(ids)).mean()),"program_rmse":m["rmse"],"program_mae":m["mae"],"program_correlation":m["correlation"],"nearest_program_rmse":b["rmse"],"nearest_program_correlation":b["correlation"],"rmse_improvement_vs_nearest":b["rmse"]-m["rmse"]})
    mem=pd.DataFrame([{"held_out_dataset":held_out,"program":int(p),"gene_index":int(i),"gene":str(matrix.index[int(i)])} for p,g in programs.items() for i in g])
    return pd.DataFrame(rows),stab,mem,{"reason":"ok","n_train_trajectories":len(used),"n_stable_genes":len(ids),"min_datasets":min_ds}


def _summary(d):
    if d.empty:return pd.DataFrame()
    return d.groupby("held_out_dataset").agg(n_cases=("time_hours","size"),mean_rmse=("program_rmse","mean"),mean_mae=("program_mae","mean"),mean_correlation=("program_correlation","mean"),mean_nearest_rmse=("nearest_program_rmse","mean"),mean_rmse_improvement=("rmse_improvement_vs_nearest","mean"),mean_gene_stability=("mean_gene_stability","mean")).reset_index()


def _perm_times(meta,seed):
    rng=np.random.default_rng(seed);out={}
    for ds,g in meta.groupby("dataset"):
        q=g[g.time_hours.notna()&g.matrix_column.notna()];vals=q.time_hours.astype(float).to_numpy()
        if len(vals)>=3:
            for c,t in zip(q.matrix_column.astype(str),rng.permutation(vals)):out[(str(ds),str(c))]=float(t)
    return out


def stage2_9_1(n_programs=10,max_genes=2500,stability_quantile=.75,permutations=20,seed=42):
    print("Stage 2.9.1: loading common gene space...",flush=True);matrix,meta=_load_common_space();datasets=sorted(meta.dataset.astype(str).unique());diags=[];members=[];stabs=[];statuses=[]
    for i,ds in enumerate(datasets,1):
        print(f"Stage 2.9.1: fold {i}/{len(datasets)}, held out={ds}...",flush=True);d,s,m,status=_fold(matrix,meta,ds,n_programs,max_genes,seed)
        status.update({"held_out_dataset":ds,"mode":"observed"});statuses.append(status)
        if status["reason"]!="ok":print(f"  skipped: {status}",flush=True)
        if d is not None and not d.empty:diags.append(d);members.append(m)
        if s is not None and not s.empty:s=s.copy();s.insert(0,"held_out_dataset",ds);stabs.append(s)
    diag=pd.concat(diags,ignore_index=True) if diags else pd.DataFrame();mem=pd.concat(members,ignore_index=True) if members else pd.DataFrame();stab=pd.concat(stabs,ignore_index=True) if stabs else pd.DataFrame();diag.to_csv(OUT/"01_fold_summary.csv",index=False);mem.to_csv(OUT/"02_program_membership.csv",index=False);stab.to_csv(OUT/"03_gene_stability.csv",index=False);pd.DataFrame(statuses).to_csv(OUT/"06_fold_status.csv",index=False);summary=_summary(diag);summary.to_csv(OUT/"04_program_state_summary.csv",index=False)
    null=[]
    for p in range(permutations):
        ps=seed+1000+p;pt=_perm_times(meta,ps);pdg=[]
        for ds in datasets:
            d,_,_,_=_fold(matrix,meta,ds,n_programs,max_genes,ps,pt)
            if d is not None and not d.empty:pdg.append(d)
        if pdg:
            q=pd.concat(pdg,ignore_index=True);null.append({"permutation":p,"mean_rmse_improvement":float(q.rmse_improvement_vs_nearest.mean()),"mean_correlation":float(q.program_correlation.mean()),"n_cases":len(q)})
        print(f"Stage 2.9.1: permutation {p+1}/{permutations}",flush=True)
    null=pd.DataFrame(null);null.to_csv(OUT/"05_permutation_null.csv",index=False);print("Stage 2.9.1 complete.")
    if summary.empty:print("No valid observed folds were produced. Inspect stage2_9_1/06_fold_status.csv for the exact reason per dataset.")
    else:
        print(summary.to_string(index=False));imp=float(summary.mean_rmse_improvement.mean());corr=float(summary.mean_correlation.mean());print(f"Observed mean RMSE improvement: {imp:.4f}; correlation: {corr:.3f}")
        if not null.empty:
            p_imp=(1+int((null.mean_rmse_improvement>=imp).sum()))/(len(null)+1);p_corr=(1+int((null.mean_correlation>=corr).sum()))/(len(null)+1);print(f"Permutation empirical p: improvement={p_imp:.3f}, correlation={p_corr:.3f}")
        print("Programs are data-derived invariant temporal programs, not curated pathways. Fit ODE/state-space only if transfer beats baseline and permutation null.")
    return summary


if __name__=="__main__":stage2_9_1()
