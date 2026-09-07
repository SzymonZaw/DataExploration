"""Stage 2.9.30: time-warped shared trajectory validation.

Tests whether heterogeneous trajectory datasets can be explained by one shared
progression curve traversed at different speeds. The held-out dataset is never
used for gene selection or PCA fitting. Its final state is predicted from a
DTW alignment learned only from the held-out prefix. No ODE/state-space model
is fitted.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_30"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
GRID = np.linspace(0.0, 1.0, 40)
MIN_GENES = 100
MAX_GENES = 3000
N_PERM = 100
N_PC = 10


def log(x): print(f"Stage 2.9.30: {x}", flush=True)


def finite(X):
    X = np.asarray(X, float).copy()
    for j in range(X.shape[1]):
        ok = np.isfinite(X[:, j])
        X[~ok, j] = np.median(X[ok, j]) if ok.any() else 0.0
    return X


def normalize_time(t):
    t = np.asarray(t, float); lo, hi = float(t.min()), float(t.max())
    return np.zeros_like(t) if hi <= lo else (t - lo) / (hi - lo)


def resample(t, X):
    tn = normalize_time(t); X = finite(X)
    Y = np.empty((len(GRID), X.shape[1]))
    for j in range(X.shape[1]): Y[:, j] = np.interp(GRID, tn, X[:, j])
    return Y


def normalize_with_stats(Y, baseline=None, amp=None):
    Y = finite(Y)
    if baseline is None: baseline = Y[0].copy()
    Z = Y - baseline[None, :]
    if amp is None:
        amp = np.sqrt(np.mean(Z * Z, axis=0))
        amp = np.where(amp > 1e-8, amp, 1.0)
    return Z / amp, baseline, amp


def gene_indices(genes, all_genes):
    pos = {str(g): i for i, g in enumerate(all_genes)}
    return [pos[str(g)] for g in genes if str(g) in pos]


def select_genes(train, all_genes):
    blocks = []
    for t, X, _ in train.values():
        Y, _, _ = normalize_with_stats(resample(t, X))
        blocks.append(Y)
    A = np.stack(blocks); mean = A.mean(axis=0)
    shared = np.var(mean, axis=0)
    hetero = np.mean((A - mean[None]) ** 2, axis=(0, 1))
    score = shared / (shared + hetero + 1e-12)
    order = np.argsort(score)[::-1]
    chosen = order[score[order] >= 0.5]
    if len(chosen) < MIN_GENES: chosen = order[:min(MIN_GENES, len(order))]
    return np.asarray(all_genes)[chosen[:MAX_GENES]]


def load_trajectories():
    from dynamics.validation import _load_common_space, _time_hours_for_validation, _strip_dataset_prefix
    m, meta = _load_common_space()
    meta = meta[meta["dataset"].astype(str).isin(TARGET)].copy()
    meta["matrix_column"] = meta["matrix_column"].astype(str)
    times = []
    for i, (_, r) in enumerate(meta.iterrows()):
        ds = str(r["dataset"]); sample = _strip_dataset_prefix(str(r["sample"]))
        times.append(_time_hours_for_validation(ds, sample, i if ds == "GSE28688" else None))
    meta["time_hours"] = pd.to_numeric(times, errors="coerce")
    out = {}
    for ds in TARGET:
        g = meta[(meta["dataset"].astype(str) == ds) & np.isfinite(meta["time_hours"])].copy()
        if g["time_hours"].nunique() < 3: continue
        X = m.loc[:, list(g["matrix_column"])].T.copy()
        X["time_hours"] = g["time_hours"].to_numpy(float)
        X = X.groupby("time_hours", sort=True).mean()
        out[ds] = (X.index.to_numpy(float), X.drop(columns="time_hours").to_numpy(float), list(m.index))
    return out


def fit_artifact(train, all_genes):
    genes = select_genes(train, all_genes); idx = gene_indices(genes, all_genes)
    states = []
    for t, X, _ in train.values():
        Y, _, _ = normalize_with_stats(resample(t, X[:, idx]))
        states.append(Y)
    Z = np.vstack(states)
    ncomp = min(N_PC, Z.shape[0], Z.shape[1])
    pca = PCA(n_components=ncomp, random_state=2930).fit(Z)
    template = np.mean(np.stack([pca.transform(Y) for Y in states]), axis=0)
    return genes, idx, pca, template


def dtw_path(A, B):
    n, m = len(A), len(B)
    dp = np.full((n + 1, m + 1), np.inf); dp[0, 0] = 0.0
    prev = np.full((n + 1, m + 1), -1, dtype=np.int8)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = float(np.sqrt(np.mean((A[i-1] - B[j-1]) ** 2)))
            c = (dp[i-1,j-1], dp[i-1,j], dp[i,j-1]); k = int(np.argmin(c))
            dp[i,j] = d + c[k]; prev[i,j] = k
    i, j = n, m; path = []
    while i and j:
        path.append((i-1,j-1)); k = prev[i,j]
        if k == 0: i -= 1; j -= 1
        elif k == 1: i -= 1
        else: j -= 1
    path.reverse()
    return path, float(dp[n,m] / max(1,len(path)))


def infer_progress(prefix_pc, template):
    path, cost = dtw_path(prefix_pc, template)
    mapping = {i: [] for i in range(len(prefix_pc))}
    for i,j in path: mapping[i].append(j)
    q=[]; p=[]
    for i, js in mapping.items():
        if js: q.append(GRID[i]); p.append(float(np.mean(GRID[js])))
    q,p=np.asarray(q),np.asarray(p)
    if len(q)<3: return np.nan,cost,np.nan
    slope, intercept = np.polyfit(q,p,1)
    return float(np.clip(slope+intercept,0,1)), cost, float(slope)


def interp_state(template, progress):
    return np.array([np.interp(float(np.clip(progress,0,1)), GRID, template[:,j]) for j in range(template.shape[1])])


def metrics(pred, true):
    return float(np.sqrt(np.mean((pred-true)**2))), float(pd.Series(pred).corr(pd.Series(true)))


def prepare_fold(heldout, train_names, traj):
    train={d:traj[d] for d in train_names}; t,X,all_genes=traj[heldout]
    genes,idx,pca,template=fit_artifact(train,all_genes)
    prefix_t,prefix_X=t[:-1],X[:-1]
    prefix_raw=resample(prefix_t,prefix_X[:,idx]); prefix_state,base,amp=normalize_with_stats(prefix_raw)
    full_raw=resample(t,X[:,idx]); full_state,_,_=normalize_with_stats(full_raw,base,amp)
    return {"genes":genes,"idx":idx,"pca":pca,"template":template,"prefix_pc":pca.transform(prefix_state),"true_pc":pca.transform(full_state)[-1],"persistence":pca.transform(prefix_state)[-1],"unwarped":template[-1]}


def evaluate(artifact, perm_seed=None):
    prefix=artifact["prefix_pc"].copy()
    if perm_seed is not None: prefix=prefix[np.random.default_rng(perm_seed).permutation(len(prefix))]
    progress,cost,slope=infer_progress(prefix,artifact["template"])
    if not np.isfinite(progress): return None
    warped=interp_state(artifact["template"],progress); true=artifact["true_pc"]
    rw,cw=metrics(warped,true); rp,cp=metrics(artifact["persistence"],true); ru,cu=metrics(artifact["unwarped"],true)
    return {"warp_progress_at_final_time":progress,"warp_slope":slope,"dtw_prefix_cost":cost,"rmse_time_warp":rw,"rmse_persistence":rp,"rmse_unwarped_shared":ru,"pearson_time_warp":cw,"pearson_persistence":cp,"pearson_unwarped_shared":cu,"improvement_vs_persistence":rp-rw,"improvement_vs_unwarped_shared":ru-rw}


def run():
    traj=load_trajectories(); names=[d for d in TARGET if d in traj]
    log(f"trajectory datasets: {', '.join(names)}")
    if len(names)<3: raise RuntimeError("Stage 2.9.30 requires at least three trajectory datasets.")
    artifacts={}; rows=[]
    for heldout in names:
        train_names=[d for d in names if d!=heldout]
        log(f"LODO held out {heldout}; training on {', '.join(train_names)}")
        art=prepare_fold(heldout,train_names,traj); artifacts[heldout]=art
        r=evaluate(art)
        if r: r.update({"heldout_dataset":heldout,"n_selected_genes":len(art["genes"]) }); rows.append(r)
    R=pd.DataFrame(rows); R.to_csv(OUT/"01_lodo_time_warp_validation.csv",index=False)
    if R.empty: raise RuntimeError("No valid Stage 2.9.30 LODO folds.")

    null_rows=[]
    for p in range(N_PERM):
        vals=[]
        for di,heldout in enumerate(names):
            r=evaluate(artifacts[heldout],700000+p*len(names)+di)
            if r: vals.append(r["improvement_vs_persistence"])
        null_rows.append({"permutation":p+1,"mean_improvement_vs_persistence":float(np.mean(vals)) if vals else np.nan})
    P=pd.DataFrame(null_rows); P.to_csv(OUT/"02_time_permutation_null.csv",index=False)
    obs=float(R["improvement_vs_persistence"].mean()); null=P["mean_improvement_vs_persistence"].dropna().to_numpy(float)
    p=float((1+np.sum(null>=obs))/(len(null)+1)) if len(null) else np.nan
    summary=pd.DataFrame([{"n_trajectory_datasets":len(names),"n_valid_lodo_folds":len(R),"mean_rmse_time_warp":float(R["rmse_time_warp"].mean()),"mean_rmse_persistence":float(R["rmse_persistence"].mean()),"mean_rmse_unwarped_shared":float(R["rmse_unwarped_shared"].mean()),"mean_improvement_vs_persistence":obs,"mean_improvement_vs_unwarped_shared":float(R["improvement_vs_unwarped_shared"].mean()),"mean_pearson_time_warp":float(R["pearson_time_warp"].mean()),"time_warp_permutation_p":p,"time_warp_predictive_support":bool(obs>0 and np.isfinite(p) and p<0.05),"stage3_readiness":False,"interpretation":"Training-only shared trajectory with held-out-prefix DTW time warping; final-point LODO prediction; representation/predictive diagnostic only; no ODE/state-space model"}])
    summary.to_csv(OUT/"03_stage2930_summary.csv",index=False); log("overall:"); print(summary.to_string(index=False)); return summary


if __name__=="__main__": run()
