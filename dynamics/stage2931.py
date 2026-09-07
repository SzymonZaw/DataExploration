"""Stage 2.9.31: shared vs dataset-specific dynamics.

Tests whether a common temporal component plus a dataset-specific residual
predicts a held-out trajectory better than persistence or shared-only.
This is a predictive representation diagnostic, not an ODE/state-space model.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_31"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
GRID = np.linspace(0.0, 1.0, 40)
MIN_GENES = 100
MAX_GENES = 3000
N_PC = 10
N_PERM = 1000


def log(x):
    print(f"Stage 2.9.31: {x}", flush=True)


def finite(X):
    X = np.asarray(X, float).copy()
    for j in range(X.shape[1]):
        ok = np.isfinite(X[:, j])
        X[~ok, j] = np.median(X[ok, j]) if ok.any() else 0.0
    return X


def normalize_time(t):
    t = np.asarray(t, float)
    lo, hi = float(np.min(t)), float(np.max(t))
    return np.zeros_like(t) if hi <= lo else (t - lo) / (hi - lo)


def resample(t, X):
    tn = normalize_time(t)
    X = finite(X)
    Y = np.empty((len(GRID), X.shape[1]))
    for j in range(X.shape[1]):
        Y[:, j] = np.interp(GRID, tn, X[:, j])
    return Y


def normalize_with_stats(Y, baseline=None, amp=None):
    Y = finite(Y)
    baseline = Y[0].copy() if baseline is None else np.asarray(baseline, float)
    Z = Y - baseline[None, :]
    if amp is None:
        amp = np.sqrt(np.mean(Z * Z, axis=0))
        amp = np.where(amp > 1e-8, amp, 1.0)
    return Z / amp, baseline, amp


def gene_indices(genes, all_genes):
    pos = {str(g): i for i, g in enumerate(all_genes)}
    return [pos[str(g)] for g in genes if str(g) in pos]


def select_genes(train, all_genes):
    states = []
    for t, X, _ in train.values():
        states.append(normalize_with_stats(resample(t, X))[0])
    A = np.stack(states)
    mean = A.mean(axis=0)
    shared = np.var(mean, axis=0)
    hetero = np.mean((A - mean[None, :, :]) ** 2, axis=(0, 1))
    score = shared / (shared + hetero + 1e-12)
    order = np.argsort(score)[::-1]
    chosen = order[score[order] >= 0.5]
    if len(chosen) < MIN_GENES:
        chosen = order[:min(MIN_GENES, len(order))]
    return np.asarray(all_genes)[chosen[:MAX_GENES]]


def load_trajectories():
    # Reuse Stage 2.7's validated sample/time/matrix mapping.
    from dynamics.validation import _load_common_space
    matrix, metadata = _load_common_space()
    metadata = metadata[metadata["dataset"].astype(str).isin(TARGET)].copy()
    out = {}
    for ds in TARGET:
        g = metadata[
            (metadata["dataset"].astype(str) == ds)
            & metadata["matrix_column"].notna()
            & metadata["time_hours"].notna()
        ].copy()
        g["time_hours"] = pd.to_numeric(g["time_hours"], errors="coerce")
        g = g[np.isfinite(g["time_hours"].to_numpy(float))]
        if g["time_hours"].nunique() < 3:
            log(f"{ds}: skipped; timed samples={len(g)}, unique times={g['time_hours'].nunique()}")
            continue
        cols = g["matrix_column"].astype(str).tolist()
        expr = matrix.loc[:, cols].T.copy()
        expr.index = g["time_hours"].to_numpy(float)
        expr = expr.groupby(level=0, sort=True).mean()
        log(f"{ds}: {len(g)} timed samples -> {len(expr)} unique timepoints")
        out[ds] = (expr.index.to_numpy(float), expr.to_numpy(float), list(matrix.index))
    return out


def fit_shared_model(train, all_genes):
    genes = select_genes(train, all_genes)
    idx = gene_indices(genes, all_genes)
    states = [normalize_with_stats(resample(t, X[:, idx]))[0] for t, X, _ in train.values()]
    Z = np.vstack(states)
    n_pc = min(N_PC, Z.shape[0], Z.shape[1])
    pca = PCA(n_components=n_pc, random_state=2931).fit(Z)
    pcs = [pca.transform(Y) for Y in states]
    shared_pc = np.mean(np.stack(pcs), axis=0)
    return genes, idx, pca, shared_pc


def fit_residual_extrapolation(residual):
    """Fit residual on the normalized training grid and extrapolate to grid=1."""
    residual = np.asarray(residual, float)
    if residual.shape[0] < 3:
        return np.full(residual.shape[1], np.nan), np.nan
    design = np.column_stack([np.ones(len(GRID)), GRID])
    coef = np.linalg.lstsq(design, residual, rcond=None)[0]
    pred = coef[0] + coef[1]
    return pred, float(np.mean(np.abs(coef[1])))


def metrics(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    rmse = float(np.sqrt(np.mean((a - b) ** 2)))
    pa = pd.Series(a).corr(pd.Series(b))
    return rmse, float(pa) if pd.notna(pa) else np.nan


def prepare_fold(heldout, train_names, traj):
    train = {d: traj[d] for d in train_names}
    t, X, all_genes = traj[heldout]
    if len(t) < 4:
        return None
    genes, idx, pca, shared = fit_shared_model(train, all_genes)

    # Held-out normalization uses only the observed prefix. Both prefix and
    # shared trajectory live on the same normalized 40-point grid.
    prefix_raw = resample(t[:-1], X[:-1, idx])
    prefix_state, baseline, amp = normalize_with_stats(prefix_raw)
    full_raw = resample(t, X[:, idx])
    full_state, _, _ = normalize_with_stats(full_raw, baseline, amp)

    prefix_pc = pca.transform(prefix_state)
    true_pc = pca.transform(full_state)[-1]

    # prefix_pc is already sampled on GRID, so subtraction is dimensionally
    # and temporally aligned with the training-derived shared trajectory.
    residual_prefix = prefix_pc - shared
    residual_final, slope = fit_residual_extrapolation(residual_prefix)
    if not np.all(np.isfinite(residual_final)):
        return None

    pred_full = shared[-1] + residual_final
    pred_shared = shared[-1]
    persistence = prefix_pc[-1]

    # Nearest-time baseline: shared state at the last observed normalized time.
    last_progress = float(normalize_time(t[:-1])[-1])
    pred_nearest = np.array([
        np.interp(last_progress, GRID, shared[:, j])
        for j in range(shared.shape[1])
    ])

    return {
        "heldout_dataset": heldout,
        "n_selected_genes": len(genes),
        "n_pc": shared.shape[1],
        "prefix_pc": prefix_pc,
        "true_pc": true_pc,
        "shared_final": shared[-1],
        "shared_template": shared,
        "pred_full": pred_full,
        "pred_shared": pred_shared,
        "persistence": persistence,
        "pred_nearest": pred_nearest,
        "residual_slope": slope,
    }


def evaluate(artifact):
    true = artifact["true_pc"]
    rw, cw = metrics(artifact["pred_full"], true)
    rs, cs = metrics(artifact["pred_shared"], true)
    rp, cp = metrics(artifact["persistence"], true)
    rn, cn = metrics(artifact["pred_nearest"], true)
    return {
        "heldout_dataset": artifact["heldout_dataset"],
        "n_selected_genes": artifact["n_selected_genes"],
        "rmse_shared_plus_dataset": rw,
        "rmse_shared_only": rs,
        "rmse_persistence": rp,
        "rmse_nearest_shared": rn,
        "improvement_vs_persistence": rp - rw,
        "improvement_vs_shared_only": rs - rw,
        "improvement_vs_nearest": rn - rw,
        "pearson_shared_plus_dataset": cw,
        "pearson_shared_only": cs,
        "pearson_persistence": cp,
        "pearson_nearest_shared": cn,
        "residual_slope": artifact["residual_slope"],
    }


def permutation_null(artifacts, names):
    rows = []
    for p in range(N_PERM):
        vals = []
        for di, ds in enumerate(names):
            if ds not in artifacts:
                continue
            a = artifacts[ds]
            rng = np.random.default_rng(931000 + p * len(names) + di)
            perm = a["prefix_pc"][rng.permutation(len(a["prefix_pc"]))]
            template = a["shared_template"]
            residual = perm - template
            pred_resid, _ = fit_residual_extrapolation(residual)
            if not np.all(np.isfinite(pred_resid)):
                continue
            pred = a["shared_final"] + pred_resid
            true = a["true_pc"]
            rmse = float(np.sqrt(np.mean((pred - true) ** 2)))
            persistence = float(np.sqrt(np.mean((a["persistence"] - true) ** 2)))
            vals.append(persistence - rmse)
        rows.append({
            "permutation": p + 1,
            "mean_improvement_vs_persistence": float(np.mean(vals)) if vals else np.nan,
        })
    return pd.DataFrame(rows)


def run():
    traj = load_trajectories()
    names = [d for d in TARGET if d in traj]
    log(f"trajectory datasets: {', '.join(names)}")
    if len(names) < 3:
        raise RuntimeError("Stage 2.9.31 requires at least three trajectory datasets.")

    artifacts, rows = {}, []
    for heldout in names:
        train = [d for d in names if d != heldout]
        log(f"LODO held out {heldout}; training on {', '.join(train)}")
        a = prepare_fold(heldout, train, traj)
        if a is not None:
            artifacts[heldout] = a
            rows.append(evaluate(a))

    R = pd.DataFrame(rows)
    R.to_csv(OUT / "01_lodo_shared_dataset_dynamics.csv", index=False)
    if R.empty:
        raise RuntimeError("No valid Stage 2.9.31 LODO folds.")

    P = permutation_null(artifacts, names)
    P.to_csv(OUT / "02_time_permutation_null.csv", index=False)
    obs = float(R["improvement_vs_persistence"].mean())
    nv = P["mean_improvement_vs_persistence"].dropna().to_numpy(float)
    p_value = float((1 + np.sum(nv >= obs)) / (len(nv) + 1)) if len(nv) else np.nan

    summary = pd.DataFrame([{
        "n_trajectory_datasets": len(names),
        "n_valid_lodo_folds": len(R),
        "mean_selected_genes": float(R["n_selected_genes"].mean()),
        "mean_rmse_shared_plus_dataset": float(R["rmse_shared_plus_dataset"].mean()),
        "mean_rmse_shared_only": float(R["rmse_shared_only"].mean()),
        "mean_rmse_persistence": float(R["rmse_persistence"].mean()),
        "mean_rmse_nearest_shared": float(R["rmse_nearest_shared"].mean()),
        "mean_improvement_vs_persistence": obs,
        "mean_improvement_vs_shared_only": float(R["improvement_vs_shared_only"].mean()),
        "mean_improvement_vs_nearest": float(R["improvement_vs_nearest"].mean()),
        "mean_pearson_shared_plus_dataset": float(R["pearson_shared_plus_dataset"].mean()),
        "time_permutation_p": p_value,
        "shared_plus_dataset_predictive_support": bool(obs > 0 and np.isfinite(p_value) and p_value < 0.05),
        "stage3_readiness": False,
        "interpretation": "Training-only shared temporal component plus held-out dataset-specific residual extrapolation; LODO final-point prediction; permutation null; no ODE/state-space model",
    }])
    summary.to_csv(OUT / "03_stage2931_summary.csv", index=False)
    log("overall:")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    run()
