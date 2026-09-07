"""Stage 2.9.30: time-warped shared trajectory validation.

Tests whether trajectory datasets can be explained by one shared progression
curve traversed at dataset-specific speeds. The held-out dataset is never
used to build the shared template or select genes. A monotone DTW alignment
is inferred from the held-out prefix only, then used to predict its final
observed state. This is a predictive validation, not an ODE/state-space fit.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_30"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
GRID = np.linspace(0.0, 1.0, 40)
MIN_GENES = 100
MAX_GENES = 3000
N_PERM = 1000


def log(x):
    print(f"Stage 2.9.30: {x}", flush=True)


def finite_matrix(X):
    X = np.asarray(X, float).copy()
    for j in range(X.shape[1]):
        ok = np.isfinite(X[:, j])
        fill = float(np.median(X[ok, j])) if ok.any() else 0.0
        X[~ok, j] = fill
    return X


def normalize_time(t):
    t = np.asarray(t, float)
    lo, hi = float(np.min(t)), float(np.max(t))
    if hi <= lo:
        return np.zeros_like(t)
    return (t - lo) / (hi - lo)


def resample(t, X):
    t = normalize_time(t)
    X = finite_matrix(X)
    Y = np.empty((len(GRID), X.shape[1]), float)
    for j in range(X.shape[1]):
        Y[:, j] = np.interp(GRID, t, X[:, j])
    return Y


def gene_scores(train_traj):
    blocks = [resample(t, X) for t, X in train_traj.values()]
    A = np.stack(blocks)
    mean = A.mean(axis=0)
    shared = np.var(mean, axis=0)
    hetero = np.mean((A - mean[None, :, :]) ** 2, axis=(0, 1))
    score = shared / (shared + hetero + 1e-12)
    return score


def select_genes(train_traj, gene_names):
    score = gene_scores(train_traj)
    order = np.argsort(score)[::-1]
    chosen = order[score[order] >= 0.5]
    if len(chosen) < MIN_GENES:
        chosen = order[:min(MIN_GENES, len(order))]
    chosen = chosen[:min(MAX_GENES, len(chosen))]
    return np.asarray(gene_names)[chosen]


def trajectory_matrix(t, X, genes, all_genes):
    pos = {str(g): i for i, g in enumerate(all_genes)}
    idx = [pos[str(g)] for g in genes if str(g) in pos]
    return normalize_state(resample(t, X[:, idx]))


def normalize_state(Y):
    Y = finite_matrix(Y)
    Y = Y - Y[0:1, :]
    amp = np.sqrt(np.mean(Y * Y, axis=0))
    amp = np.where(amp > 1e-8, amp, 1.0)
    return Y / amp


def dtw_path(A, B):
    """Classic monotone DTW path between time x feature matrices."""
    A = finite_matrix(A)
    B = finite_matrix(B)
    n, m = len(A), len(B)
    dp = np.full((n + 1, m + 1), np.inf)
    dp[0, 0] = 0.0
    prev = np.full((n + 1, m + 1), -1, dtype=np.int8)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = float(np.sqrt(np.mean((A[i - 1] - B[j - 1]) ** 2)))
            choices = (dp[i - 1, j - 1], dp[i - 1, j], dp[i, j - 1])
            k = int(np.argmin(choices))
            dp[i, j] = d + choices[k]
            prev[i, j] = k
    i, j = n, m
    path = []
    while i > 0 and j > 0:
        path.append((i - 1, j - 1))
        k = prev[i, j]
        if k == 0:
            i -= 1; j -= 1
        elif k == 1:
            i -= 1
        else:
            j -= 1
    path.reverse()
    return path, float(dp[n, m] / max(1, len(path)))


def template(train_traj, genes, all_genes):
    mats = [trajectory_matrix(t, X, genes, all_genes) for t, X in train_traj.values()]
    return np.mean(np.stack(mats), axis=0)


def nearest_grid_prediction(template_matrix, progress):
    p = float(np.clip(progress, 0.0, 1.0))
    return np.array([np.interp(p, GRID, template_matrix[:, j]) for j in range(template_matrix.shape[1])])


def infer_warp_progress(prefix_t, prefix_X, template_matrix, genes, all_genes):
    prefix = trajectory_matrix(prefix_t, prefix_X, genes, all_genes)
    # prefix is already resampled to GRID based on the observed prefix only.
    # Align its complete observed prefix to the training template.
    path, cost = dtw_path(prefix, template_matrix)
    mapping = {i: [] for i in range(len(prefix))}
    for i, j in path:
        mapping[i].append(j)
    q = []
    p = []
    for i, js in mapping.items():
        if js:
            q.append(GRID[i])
            p.append(float(np.mean(GRID[js])))
    q = np.asarray(q, float)
    p = np.asarray(p, float)
    ok = np.isfinite(q) & np.isfinite(p)
    q, p = q[ok], p[ok]
    if len(q) < 3:
        return np.nan, cost, np.nan
    # Fit a monotone-free linear clock on the aligned prefix, then extrapolate
    # only from prefix observations to the known final normalized time 1.0.
    slope, intercept = np.polyfit(q, p, 1)
    pred_progress = float(np.clip(slope * 1.0 + intercept, 0.0, 1.0))
    return pred_progress, cost, float(slope)


def load_trajectories():
    from dynamics.validation import _load_common_space, _time_hours_for_validation, _strip_dataset_prefix
    m, meta = _load_common_space()
    meta = meta[meta["dataset"].astype(str).isin(TARGET)].copy()
    meta["matrix_column"] = meta["matrix_column"].astype(str)
    meta["time_hours"] = [
        _time_hours_for_validation(str(r["dataset"]), _strip_dataset_prefix(str(r["sample"])), None)
        for _, r in meta.iterrows()
    ]
    meta["time_hours"] = pd.to_numeric(meta["time_hours"], errors="coerce")
    out = {}
    for ds in TARGET:
        g = meta[(meta["dataset"].astype(str) == ds) & np.isfinite(meta["time_hours"])].copy()
        if g["time_hours"].nunique() < 3:
            continue
        cols = list(g["matrix_column"])
        X = m.loc[:, cols].T.copy()
        X["time_hours"] = g["time_hours"].to_numpy(float)
        X = X.groupby("time_hours", sort=True).mean()
        out[ds] = (X.index.to_numpy(float), X.drop(columns=["time_hours"]).to_numpy(float), list(m.index))
    return out


def evaluate_fold(heldout, train_names, traj):
    train = {ds: traj[ds] for ds in train_names}
    t, X, genes_all = traj[heldout]
    if len(t) < 4:
        return None
    # Strictly training-only feature selection and template construction.
    genes = select_genes(train, genes_all)
    T = template(train, genes, genes_all)

    prefix_t, prefix_X = t[:-1], X[:-1]
    true_state = trajectory_matrix(np.array([t[-2], t[-1]]), X[-2:], genes, genes_all)[-1]
    pred_progress, warp_cost, slope = infer_warp_progress(prefix_t, prefix_X, T, genes, genes_all)
    if not np.isfinite(pred_progress):
        return None

    warped_pred = nearest_grid_prediction(T, pred_progress)
    unwarped_pred = nearest_grid_prediction(T, 1.0)
    persistence = trajectory_matrix(np.array([t[-2], t[-1]]), X[-2:], genes, genes_all)[0]
    # Nearest-time baseline: map the held-out final normalized time directly to
    # the training template, without learning a dataset-specific warp.
    nearest = nearest_grid_prediction(T, 1.0)

    def rmse(a, b):
        return float(np.sqrt(np.mean((np.asarray(a) - np.asarray(b)) ** 2)))
    def pearson(a, b):
        return float(pd.Series(a).corr(pd.Series(b)))

    return {
        "heldout_dataset": heldout,
        "n_train_datasets": len(train_names),
        "n_train_timepoints": int(sum(len(traj[d][0]) for d in train_names)),
        "n_heldout_timepoints": len(t),
        "n_selected_genes": len(genes),
        "warp_progress_at_final_time": pred_progress,
        "warp_slope": slope,
        "dtw_prefix_cost": warp_cost,
        "rmse_time_warp": rmse(warped_pred, true_state),
        "rmse_persistence": rmse(persistence, true_state),
        "rmse_unwarped_shared": rmse(unwarped_pred, true_state),
        "pearson_time_warp": pearson(warped_pred, true_state),
        "pearson_persistence": pearson(persistence, true_state),
        "pearson_unwarped_shared": pearson(unwarped_pred, true_state),
        "improvement_vs_persistence": rmse(persistence, true_state) - rmse(warped_pred, true_state),
        "improvement_vs_unwarped_shared": rmse(unwarped_pred, true_state) - rmse(warped_pred, true_state),
    }


def permutation_null(traj, observed_rows):
    rng = np.random.default_rng(293000)
    out = []
    names = list(traj)
    for p in range(N_PERM):
        fold_improvements = []
        for heldout in names:
            train_names = [d for d in names if d != heldout]
            # Destroy temporal shape in every dataset while retaining the same
            # time grid and all sample counts. This keeps the prediction task
            # unchanged while removing ordered biological dynamics.
            perm_traj = {}
            for di, ds in enumerate(names):
                t, X, genes = traj[ds]
                perm = np.random.default_rng(700000 + p * len(names) + di).permutation(len(t))
                perm_traj[ds] = (t, X[perm], genes)
            r = evaluate_fold(heldout, train_names, perm_traj)
            if r is not None:
                fold_improvements.append(r["improvement_vs_persistence"])
        out.append({"permutation": p + 1, "mean_improvement_vs_persistence": float(np.mean(fold_improvements)) if fold_improvements else np.nan})
    return pd.DataFrame(out)


def run():
    traj = load_trajectories()
    names = [d for d in TARGET if d in traj]
    log(f"trajectory datasets: {', '.join(names)}")
    if len(names) < 3:
        raise RuntimeError("Stage 2.9.30 requires at least three trajectory datasets.")

    rows = []
    for heldout in names:
        train_names = [d for d in names if d != heldout]
        log(f"LODO held out {heldout}; training on {', '.join(train_names)}")
        r = evaluate_fold(heldout, train_names, traj)
        if r is not None:
            rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "01_lodo_time_warp_validation.csv", index=False)
    if R.empty:
        raise RuntimeError("No valid Stage 2.9.30 LODO folds.")

    P = permutation_null(traj, rows)
    P.to_csv(OUT / "02_time_permutation_null.csv", index=False)
    obs = float(R["improvement_vs_persistence"].mean())
    null = P["mean_improvement_vs_persistence"].dropna().to_numpy(float)
    p = float((1 + np.sum(null >= obs)) / (len(null) + 1)) if len(null) else np.nan
    summary = pd.DataFrame([{
        "n_trajectory_datasets": len(names),
        "n_valid_lodo_folds": len(R),
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
    log("overall:")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    run()
