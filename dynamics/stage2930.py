"""Stage 2.9.30: time-warped shared trajectory validation.

Tests whether heterogeneous trajectory datasets can be explained by one shared
progression curve traversed at different speeds. The held-out dataset is never
used for gene selection, PCA fitting, or the shared template. A monotone DTW
alignment is inferred from the held-out prefix only and used to predict its
final state. This is a predictive diagnostic, not an ODE/state-space fit.
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
N_PERM = 1000
N_PC = 10


def log(x):
    print(f"Stage 2.9.30: {x}", flush=True)


def finite(X):
    X = np.asarray(X, float).copy()
    for j in range(X.shape[1]):
        ok = np.isfinite(X[:, j])
        X[~ok, j] = np.median(X[ok, j]) if ok.any() else 0.0
    return X


def nt(t):
    t = np.asarray(t, float)
    lo, hi = float(t.min()), float(t.max())
    return np.zeros_like(t) if hi <= lo else (t - lo) / (hi - lo)


def resample(t, X):
    tn = nt(t); X = finite(X)
    Y = np.empty((len(GRID), X.shape[1]))
    for j in range(X.shape[1]):
        Y[:, j] = np.interp(GRID, tn, X[:, j])
    return Y


def normalize(Y):
    Y = finite(Y)
    Y = Y - Y[0:1]
    amp = np.sqrt(np.mean(Y * Y, axis=0))
    return Y / np.where(amp > 1e-8, amp, 1.0)


def select_genes(train, all_genes):
    blocks = [normalize(resample(t, X)) for t, X, _ in train.values()]
    A = np.stack(blocks); mean = A.mean(axis=0)
    shared = np.var(mean, axis=0)
    hetero = np.mean((A - mean[None]) ** 2, axis=(0, 1))
    score = shared / (shared + hetero + 1e-12)
    order = np.argsort(score)[::-1]
    chosen = order[score[order] >= 0.5]
    if len(chosen) < MIN_GENES:
        chosen = order[:min(MIN_GENES, len(order))]
    return np.asarray(all_genes)[chosen[:MAX_GENES]]


def gene_indices(genes, all_genes):
    pos = {str(g): i for i, g in enumerate(all_genes)}
    return [pos[str(g)] for g in genes if str(g) in pos]


def state_from_prefix(t, X, idx):
    Y = resample(t, X[:, idx])
    return normalize(Y)


def fit_fold(train, heldout):
    t, X, all_genes = heldout
    genes = select_genes(train, all_genes)
    idx = gene_indices(genes, all_genes)
    train_states = [state_from_prefix(t0, X0, idx) for t0, X0, _ in train.values()]
    Ztrain = np.vstack(train_states)
    ncomp = min(N_PC, Ztrain.shape[0], Ztrain.shape[1])
    pca = PCA(n_components=ncomp, random_state=2930).fit(Ztrain)
    train_pc = [pca.transform(s) for s in train_states]
    template = np.mean(np.stack(train_pc), axis=0)
    return idx, genes, pca, template


def dtw(A, B):
    n, m = len(A), len(B)
    dp = np.full((n + 1, m + 1), np.inf); dp[0, 0] = 0.0
    prev = np.full((n + 1, m + 1), -1, dtype=np.int8)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = float(np.sqrt(np.mean((A[i-1] - B[j-1]) ** 2)))
            c = (dp[i-1, j-1], dp[i-1, j], dp[i, j-1])
            k = int(np.argmin(c)); dp[i, j] = d + c[k]; prev[i, j] = k
    i, j = n, m; path = []
    while i and j:
        path.append((i-1, j-1)); k = prev[i, j]
        if k == 0: i -= 1; j -= 1
        elif k == 1: i -= 1
        else: j -= 1
    path.reverse()
    return path, float(dp[n, m] / max(1, len(path)))


def warp_progress(prefix_pc, template):
    path, cost = dtw(prefix_pc, template)
    q, p = [], []
    mapping = {i: [] for i in range(len(prefix_pc))}
    for i, j in path: mapping[i].append(j)
    for i, js in mapping.items():
        if js: q.append(GRID[i]); p.append(float(np.mean(GRID[js])))
    q, p = np.asarray(q), np.asarray(p)
    if len(q) < 3: return np.nan, cost, np.nan
    slope, intercept = np.polyfit(q, p, 1)
    return float(np.clip(slope + intercept, 0.0, 1.0)), cost, float(slope)


def predict_from_template(template, progress):
    return np.array([np.interp(progress, GRID, template[:, j]) for j in range(template.shape[1])])


def rmse(a, b): return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))


def pearson(a, b): return float(pd.Series(a).corr(pd.Series(b)))


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


def fold_result(heldout, train_names, traj, permute=False, seed=0):
    train = {d: traj[d] for d in train_names}
    idx, genes, pca, template = fit_fold(train, traj[heldout])
    t, X, _ = traj[heldout]
    if len(t) < 4: return None
    prefix_t, prefix_X = t[:-1], X[:-1]
    # Prefix normalization is fitted without the held-out final state.
    prefix_state = state_from_prefix(prefix_t, prefix_X, idx)
    full_state = state_from_prefix(t, X, idx)
    if permute:
        rng = np.random.default_rng(seed)
        prefix_state = prefix_state[rng.permutation(len(prefix_state))]
        full_state = full_state.copy(); full_state[-1] = full_state[-1]
    prefix_pc = pca.transform(prefix_state)
    true_pc = pca.transform(full_state)[-1]
    progress, cost, slope = warp_progress(prefix_pc, template)
    if not np.isfinite(progress): return None
    warped = predict_from_template(template, progress)
    persistence = prefix_pc[-1]
    unwarped = template[-1]
    return {
        "heldout_dataset": heldout, "n_selected_genes": len(genes),
        "warp_progress_at_final_time": progress, "warp_slope": slope,
        "dtw_prefix_cost": cost, "rmse_time_warp": rmse(warped, true_pc),
        "rmse_persistence": rmse(persistence, true_pc),
        "rmse_unwarped_shared": rmse(unwarped, true_pc),
        "pearson_time_warp": pearson(warped, true_pc),
        "pearson_persistence": pearson(persistence, true_pc),
        "pearson_unwarped_shared": pearson(unwarped, true_pc),
        "improvement_vs_persistence": rmse(persistence, true_pc) - rmse(warped, true_pc),
        "improvement_vs_unwarped_shared": rmse(unwarped, true_pc) - rmse(warped, true_pc),
    }


def run():
    traj = load_trajectories(); names = [d for d in TARGET if d in traj]
    log(f"trajectory datasets: {', '.join(names)}")
    if len(names) < 3: raise RuntimeError("Stage 2.9.30 requires at least three trajectory datasets.")
    rows = []
    for heldout in names:
        train_names = [d for d in names if d != heldout]
        log(f"LODO held out {heldout}; training on {', '.join(train_names)}")
        r = fold_result(heldout, train_names, traj)
        if r: rows.append(r)
    R = pd.DataFrame(rows); R.to_csv(OUT / "01_lodo_time_warp_validation.csv", index=False)
    if R.empty: raise RuntimeError("No valid Stage 2.9.30 LODO folds.")

    # Permutation null: reuse each fold's training representation, but destroy
    # temporal ordering in the held-out prefix. This is intentionally a null for
    # the predictive gain of ordered alignment, not a second feature-selection pass.
    null_rows = []
    for p in range(N_PERM):
        vals = []
        for heldout in names:
            train_names = [d for d in names if d != heldout]
            r = fold_result(heldout, train_names, traj, permute=True, seed=700000 + p * len(names) + names.index(heldout))
            if r: vals.append(r["improvement_vs_persistence"])
        null_rows.append({"permutation": p + 1, "mean_improvement_vs_persistence": float(np.mean(vals)) if vals else np.nan})
    P = pd.DataFrame(null_rows); P.to_csv(OUT / "02_time_permutation_null.csv", index=False)
    obs = float(R["improvement_vs_persistence"].mean())
    null = P["mean_improvement_vs_persistence"].dropna().to_numpy(float)
    p = float((1 + np.sum(null >= obs)) / (len(null) + 1)) if len(null) else np.nan
    summary = pd.DataFrame([{
        "n_trajectory_datasets": len(names), "n_valid_lodo_folds": len(R),
        "mean_rmse_time_warp": float(R["rmse_time_warp"].mean()),
        "mean_rmse_persistence": float(R["rmse_persistence"].mean()),
        "mean_rmse_unwarped_shared": float(R["rmse_unwarped_shared"].mean()),
        "mean_improvement_vs_persistence": obs,
        "mean_improvement_vs_unwarped_shared": float(R["improvement_vs_unwarped_shared"].mean()),
        "mean_pearson_time_warp": float(R["pearson_time_warp"].mean()),
        "time_warp_permutation_p": p,
        "time_warp_predictive_support": bool(obs > 0 and np.isfinite(p) and p < 0.05),
        "stage3_readiness": False,
        "interpretation": "Training-only shared trajectory with held-out-prefix DTW time warping; final-point LODO prediction; representation/predictive diagnostic only; no ODE/state-space model",
    }])
    summary.to_csv(OUT / "03_stage2930_summary.csv", index=False)
    log("overall:"); print(summary.to_string(index=False)); return summary


if __name__ == "__main__": run()
