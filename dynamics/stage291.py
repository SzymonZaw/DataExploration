"""Stage 2.9.1: robust leakage-free invariant temporal programs."""
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
    corr = np.corrcoef(a, b)[0, 1] if np.std(a) > 0 and np.std(b) > 0 else np.nan
    return {"rmse": float(np.sqrt(np.mean((a-b)**2))), "mae": float(np.mean(np.abs(a-b))), "correlation": float(corr)}


def _trajectories(matrix, meta, datasets, time_override=None):
    out = {}
    for ds in datasets:
        g = meta[(meta.dataset == ds) & meta.matrix_column.notna() & meta.time_hours.notna()].copy()
        if time_override is None:
            times = g.time_hours.astype(float).to_numpy()
        else:
            times = np.asarray([time_override.get((ds, str(c)), np.nan) for c in g.matrix_column.astype(str)], float)
            ok = np.isfinite(times); g = g.loc[ok].copy(); times = times[ok]
        if len(g) < 3 or len(np.unique(times)) < 3: continue
        frame = matrix[g.matrix_column.astype(str)].T.copy(); frame.index = times
        frame = frame.groupby(level=0).mean().sort_index()
        if len(frame) >= 3:
            t = frame.index.to_numpy(float); span = t.max()-t.min()
            if span > 0: out[ds] = ((t-t.min())/span, frame.to_numpy(float))
    return out


def _profile_cube(trajs, grid=9):
    if not trajs: return None, []
    n_genes = next(iter(trajs.values()))[1].shape[1]
    gx = np.linspace(0, 1, grid); datasets = list(trajs)
    cube = np.full((n_genes, len(datasets), grid), np.nan)
    for di, ds in enumerate(datasets):
        t, y = trajs[ds]
        for j in range(n_genes):
            v = y[:, j]; ok = np.isfinite(v)
            if ok.sum() >= 3: cube[j, di] = np.interp(gx, t[ok], v[ok])
    return cube, datasets


def _stability(cube):
    rows=[]
    for j in range(cube.shape[0]):
        p=cube[j]; valid=[i for i in range(p.shape[0]) if np.isfinite(p[i]).sum()>=3]; cs=[]
        for ai in range(len(valid)):
            for bi in range(ai+1,len(valid)):
                x,y=p[valid[ai]],p[valid[bi]]
                if np.std(x)>1e-10 and np.std(y)>1e-10: cs.append(np.corrcoef(x,y)[0,1])
        rows.append({"gene_index":j,"n_datasets":len(valid),"n_pairs":len(cs),
                     "median_pairwise_correlation":float(np.median(cs)) if cs else np.nan,
                     "mean_pairwise_correlation":float(np.mean(cs)) if cs else np.nan})
    return pd.DataFrame(rows)


def _fit_programs(cube, stability, n_programs=10, max_genes=2500, quantile=.75, seed=42):
    s=stability[(stability.n_datasets>=3) & np.isfinite(stability.median_pairwise_correlation)].copy()
    if s.empty: raise RuntimeError("No genes pass temporal stability QC.")
    threshold=float(s.median_pairwise_correlation.quantile(quantile))
    s=s[s.median_pairwise_correlation>=threshold].sort_values("median_pairwise_correlation",ascending=False).head(max_genes)
    ids=s.gene_index.astype(int).to_numpy()
    if len(ids)<n_programs: raise RuntimeError("Too few stable genes for requested programs.")
    X=cube[ids].reshape(len(ids),-1); med=np.nanmedian(X,axis=0); X=np.where(np.isfinite(X),X,med)
    X=StandardScaler().fit_transform(X)
    labels=KMeans(n_clusters=min(n_programs,len(ids)),random_state=seed,n_init=20).fit_predict(X)
    return {i:ids[labels==i].tolist() for i in np.unique(labels)}, s


def _state(matrix, meta, programs, train):
    fit=meta[meta.dataset.isin(train)&meta.matrix_column.notna()&meta.time_hours.notna()]
    state=pd.DataFrame(index=[str(c) for c in matrix.columns])
    for pi,ids in programs.items():
        genes=[matrix.index[int(i)] for i in ids]
        vals=matrix.loc[genes,fit.matrix_column.astype(str).tolist()].to_numpy(float)
        mu=np.nanmedian(vals,axis=1); mad=1.4826*np.nanmedian(np.abs(vals-mu[:,None]),axis=1); mad[~np.isfinite(mad)|(mad<1e-8)]=1
        v=matrix.loc[genes,state.index].to_numpy(float); z=(np.where(np.isfinite(v),v,mu[:,None])-mu[:,None])/mad[:,None]
        state[f"program_{pi:02d}"]=np.nanmedian(z,axis=0)
    return state


def _state_traj(state,meta):
    out={}
    for ds,g in meta.groupby("dataset"):
        g=g[g.matrix_column.notna()&g.time_hours.notna()].copy()
        if len(g)<2 or len(g.time_hours.unique())<2: continue
        f=state.loc[g.matrix_column.astype(str)].copy(); f.index=g.time_hours.astype(float).to_numpy(); f=f.groupby(level=0).mean().sort_index()
        if len(f)>=2: out[ds]=f
    return out


def _interp(f,t):
    if f is None or t<f.index.min() or t>f.index.max(): return None
    x=f.index.to_numpy(float); y=f.to_numpy(float)
    return np.asarray([np.interp(t,x,y[:,j]) for j in range(y.shape[1])])


def _fold(matrix,meta,held_out,n_programs,max_genes,quantile,seed,time_override=None):
    ds_all=sorted(meta.dataset.astype(str).unique()); train=[d for d in ds_all if d!=held_out]
    tr=_trajectories(matrix,meta,train,time_override); te=_trajectories(matrix,meta,[held_out],time_override)
    if held_out not in te or len(tr)<2: return None,None,None
    cube,used=_profile_cube(tr); stab=_stability(cube)
    try: programs,selected=_fit_programs(cube,stab,n_programs,max_genes,quantile,seed)
    except RuntimeError: return None,stab,None
    state=_state(matrix,meta,programs,train); traj=_state_traj(state,meta); test=traj.get(held_out)
    if test is None: return None,stab,None
    rows=[]
    for t in test.index:
        preds=[_interp(traj.get(d),float(t)) for d in train if _interp(traj.get(d),float(t)) is not None]
        if not preds: continue
        pred=np.mean(preds,axis=0); actual=test.loc[t].to_numpy(float)
        near=[_interp(traj.get(d),float(t)) for d in train if traj.get(d) is not None]
        base=np.mean(near,axis=0) if near else np.full_like(pred,np.nan)
        m=_metrics(actual,pred); b=_metrics(actual,base)
        rows.append({"held_out_dataset":held_out,"time_hours":float(t),"n_programs":len(programs),"n_genes_selected":len(selected),
                     "mean_gene_stability":float(selected.median_pairwise_correlation.mean()),"program_rmse":m["rmse"],"program_mae":m["mae"],
                     "program_correlation":m["correlation"],"nearest_program_rmse":b["rmse"],"nearest_program_correlation":b["correlation"],
                     "rmse_improvement_vs_nearest":b["rmse"]-m["rmse"]})
    mem=[{"held_out_dataset":held_out,"program":int(p),"gene_index":int(i),"gene":str(matrix.index[int(i)])} for p,ids in programs.items() for i in ids]
    return pd.DataFrame(rows),stab,pd.DataFrame(mem)


def _summary(d):
    if d.empty:return pd.DataFrame()
    return d.groupby("held_out_dataset").agg(n_cases=("time_hours","size"),mean_rmse=("program_rmse","mean"),mean_mae=("program_mae","mean"),mean_correlation=("program_correlation","mean"),mean_nearest_rmse=("nearest_program_rmse","mean"),mean_rmse_improvement=("rmse_improvement_vs_nearest","mean"),mean_gene_stability=("mean_gene_stability","mean")).reset_index()


def _perm_times(meta,seed):
    rng=np.random.default_rng(seed); out={}
    for ds,g in meta.groupby("dataset"):
        q=g[g.time_hours.notna()&g.matrix_column.notna()]; vals=q.time_hours.astype(float).to_numpy()
        if len(vals)>=3:
            for c,t in zip(q.matrix_column.astype(str),rng.permutation(vals)):out[(ds,c)]=float(t)
    return out


def stage2_9_1(n_programs=10,max_genes=2500,stability_quantile=.75,permutations=20,seed=42):
    print("Stage 2.9.1: loading common gene space...",flush=True); matrix,meta=_load_common_space(); datasets=sorted(meta.dataset.astype(str).unique())
    diags=[]; members=[]; stabs=[]
    for i,ds in enumerate(datasets,1):
        print(f"Stage 2.9.1: fold {i}/{len(datasets)}, held out={ds}...",flush=True)
        d,s,m=_fold(matrix,meta,ds,n_programs,max_genes,stability_quantile,seed)
        if d is not None and not d.empty:diags.append(d);members.append(m)
        if s is not None and not s.empty:s=s.copy();s.insert(0,"held_out_dataset",ds);stabs.append(s)
    diag=pd.concat(diags,ignore_index=True) if diags else pd.DataFrame(); mem=pd.concat(members,ignore_index=True) if members else pd.DataFrame(); stab=pd.concat(stabs,ignore_index=True) if stabs else pd.DataFrame()
    diag.to_csv(OUT/"01_fold_summary.csv",index=False);mem.to_csv(OUT/"02_program_membership.csv",index=False);stab.to_csv(OUT/"03_gene_stability.csv",index=False)
    summary=_summary(diag);summary.to_csv(OUT/"04_program_state_summary.csv",index=False)
    null=[]
    for p in range(permutations):
        ps=seed+1000+p; pt=_perm_times(meta,ps); pdg=[]
        for ds in datasets:
            d,_,_=_fold(matrix,meta,ds,n_programs,max_genes,stability_quantile,ps,pt)
            if d is not None and not d.empty:pdg.append(d)
        if pdg:
            q=pd.concat(pdg,ignore_index=True);null.append({"permutation":p,"mean_rmse_improvement":float(q.rmse_improvement_vs_nearest.mean()),"mean_correlation":float(q.program_correlation.mean()),"n_cases":len(q)})
        print(f"Stage 2.9.1: permutation {p+1}/{permutations}",flush=True)
    null=pd.DataFrame(null);null.to_csv(OUT/"05_permutation_null.csv",index=False)
    print("Stage 2.9.1 complete.")
    if not summary.empty:
        print(summary.to_string(index=False)); imp=float(summary.mean_rmse_improvement.mean());corr=float(summary.mean_correlation.mean())
        print(f"Observed mean RMSE improvement: {imp:.4f}; correlation: {corr:.3f}")
        if not null.empty:
            p_imp=(1+int((null.mean_rmse_improvement>=imp).sum()))/(len(null)+1);p_corr=(1+int((null.mean_correlation>=corr).sum()))/(len(null)+1)
            print(f"Permutation empirical p: improvement={p_imp:.3f}, correlation={p_corr:.3f}")
        print("Programs are data-derived invariant temporal programs, not curated pathways. Fit ODE/state-space only if transfer beats baseline and permutation null.")
    return summary


if __name__=="__main__":stage2_9_1()
