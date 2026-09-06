"""Stage 2.9: leakage-free invariant module/state construction.

This stage deliberately does not fit an ODE. It asks whether a lower-dimensional
state built from reproducible gene programs transfers between experiments better
than the 11,899-gene representation.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9"
OUT.mkdir(parents=True, exist_ok=True)


def _metrics(a, b):
    a = np.asarray(a, float).ravel(); b = np.asarray(b, float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 2:
        return {"rmse": np.nan, "mae": np.nan, "correlation": np.nan}
    a, b = a[ok], b[ok]
    corr = np.corrcoef(a, b)[0, 1] if np.std(a) > 0 and np.std(b) > 0 else np.nan
    return {"rmse": float(np.sqrt(np.mean((a - b) ** 2))),
            "mae": float(np.mean(np.abs(a - b))),
            "correlation": float(corr)}


def _load_common_space():
    from .validation import _load_common_space as load
    return load()


def _trajectories(matrix, meta):
    out = {}
    for ds, g in meta.groupby("dataset"):
        g = g[g.matrix_column.notna() & g.time_hours.notna()].copy()
        if len(g) < 2:
            continue
        frame = pd.DataFrame(matrix[g.matrix_column].T.to_numpy(),
                            index=g.time_hours.astype(float).to_numpy(),
                            columns=matrix.index)
        frame = frame.groupby(level=0).mean().sort_index()
        if len(frame) >= 2:
            out[ds] = frame
    return out


def _interp(frame, target):
    if target < frame.index.min() or target > frame.index.max():
        return None
    x = frame.index.to_numpy(float); y = frame.to_numpy(float)
    return np.asarray([np.interp(target, x, y[:, j]) for j in range(y.shape[1])])


def _training_sample_matrix(matrix, meta, datasets):
    rows = []
    labels = []
    for ds in datasets:
        g = meta[(meta.dataset == ds) & meta.matrix_column.notna() & meta.time_hours.notna()]
        for _, r in g.iterrows():
            rows.append(matrix[str(r.matrix_column)].to_numpy(float))
            labels.append(ds)
    return np.asarray(rows), np.asarray(labels)


def _fit_modules(matrix, meta, train_datasets, n_modules=16, max_genes=4000, seed=42):
    """Fit empirical gene programs using training datasets only.

    Each gene is represented by its within-dataset standardized temporal/sample
    profile. K-means clusters genes with similar profiles. This is an empirical
    module representation; names are intentionally numeric rather than implying
    a known pathway or cell identity.
    """
    x, labels = _training_sample_matrix(matrix, meta, train_datasets)
    if len(x) < 4:
        raise RuntimeError("Too few training samples to construct invariant modules.")

    z = np.zeros_like(x, dtype=float)
    for ds in train_datasets:
        idx = np.where(labels == ds)[0]
        if len(idx) == 0:
            continue
        block = x[idx]
        mu = np.nanmean(block, axis=0)
        sd = np.nanstd(block, axis=0)
        sd[~np.isfinite(sd) | (sd < 1e-8)] = 1.0
        z[idx] = (np.nan_to_num(block, nan=mu) - mu) / sd

    gene_var = np.nanvar(z, axis=0)
    keep = np.argsort(np.nan_to_num(gene_var, nan=-np.inf))[::-1][:min(max_genes, z.shape[1])]
    gene_profiles = np.nan_to_num(z[:, keep].T, nan=0.0)
    # Remove sample-specific scale before clustering genes.
    gene_profiles = StandardScaler().fit_transform(gene_profiles)
    k = max(2, min(n_modules, len(keep)))
    km = KMeans(n_clusters=k, random_state=seed, n_init=10)
    module_ids = km.fit_predict(gene_profiles)

    modules = {i: [str(matrix.index[j]) for j, m in zip(keep, module_ids) if m == i]
               for i in range(k)}
    membership = []
    for i, genes in modules.items():
        for gene in genes:
            membership.append({"module": i, "gene": gene})
    return modules, pd.DataFrame(membership)


def _module_state(matrix, meta, modules, fit_datasets):
    """Score samples with training-only gene scaling."""
    cols = [str(c) for c in matrix.columns]
    state = pd.DataFrame(index=cols)
    all_fit = meta[meta.dataset.isin(fit_datasets) & meta.matrix_column.notna() & meta.time_hours.notna()]
    for mi, genes in modules.items():
        genes = [g for g in genes if g in matrix.index]
        if not genes:
            continue
        fit_cols = all_fit.matrix_column.astype(str).tolist()
        vals = matrix.loc[genes, fit_cols].to_numpy(float)
        mu = np.nanmean(vals, axis=1); sd = np.nanstd(vals, axis=1)
        sd[~np.isfinite(sd) | (sd < 1e-8)] = 1.0
        sample_vals = matrix.loc[genes, cols].to_numpy(float)
        sample_z = (np.nan_to_num(sample_vals, nan=mu[:, None]) - mu[:, None]) / sd[:, None]
        state[f"module_{mi:02d}"] = np.nanmean(sample_z, axis=0)
    return state


def _module_trajectories(state, meta):
    out = {}
    for ds, g in meta.groupby("dataset"):
        g = g[g.matrix_column.notna() & g.time_hours.notna()].copy()
        if len(g) < 2:
            continue
        cols = g.matrix_column.astype(str).tolist()
        frame = state.loc[cols].copy()
        frame.index = g.time_hours.astype(float).to_numpy()
        frame = frame.groupby(level=0).mean().sort_index()
        if len(frame) >= 2:
            out[ds] = frame
    return out


def _validate_fold(matrix, meta, held_out, n_modules=16, max_genes=4000, seed=42):
    train = sorted(d for d in meta.dataset.unique() if d != held_out)
    train_traj = _trajectories(matrix, meta)
    train = [d for d in train if d in train_traj]
    test = train_traj.get(held_out)
    if test is None or len(train) < 2:
        return None, None

    modules, membership = _fit_modules(matrix, meta, train, n_modules, max_genes, seed)
    state = _module_state(matrix, meta, modules, train)
    trajectories = _module_trajectories(state, meta)
    rows = []
    for t in test.index:
        preds = [_interp(trajectories[d], float(t)) for d in train if d in trajectories]
        preds = [p for p in preds if p is not None]
        if not preds:
            continue
        actual_cols = meta[(meta.dataset == held_out) & (meta.time_hours == t)].matrix_column.astype(str).tolist()
        actual = state.loc[actual_cols].mean(axis=0).to_numpy(float)
        pred = np.mean(preds, axis=0)
        nearest = []
        for d in train:
            f = trajectories.get(d)
            if f is None: continue
            nearest.append(f.iloc[int(np.argmin(np.abs(f.index.to_numpy(float) - t)))].to_numpy(float))
        base = np.mean(nearest, axis=0) if nearest else np.full_like(pred, np.nan)
        m = _metrics(actual, pred); b = _metrics(actual, base)
        rows.append({"held_out_dataset": held_out, "time_hours": float(t),
                     "n_training_datasets": len(preds), "n_modules": len(modules),
                     "n_genes_used": int(sum(map(len, modules.values()))),
                     "cross_module_rmse": m["rmse"], "cross_module_mae": m["mae"],
                     "cross_module_correlation": m["correlation"],
                     "nearest_module_rmse": b["rmse"], "nearest_module_mae": b["mae"],
                     "nearest_module_correlation": b["correlation"],
                     "rmse_improvement_vs_nearest": b["rmse"] - m["rmse"]})
    return pd.DataFrame(rows), membership


def stage2_9(n_modules=16, max_genes=4000, seed=42):
    print("Stage 2.9: loading common gene space...")
    matrix, meta = _load_common_space()
    datasets = sorted(meta.dataset.unique())
    rows = []
    memberships = []
    for i, ds in enumerate(datasets, 1):
        print(f"Stage 2.9: fold {i}/{len(datasets)}, held out={ds}...", flush=True)
        result, membership = _validate_fold(matrix, meta, ds, n_modules, max_genes, seed)
        if result is not None and not result.empty:
            rows.append(result)
            membership = membership.copy(); membership.insert(0, "held_out_dataset", ds)
            memberships.append(membership)
    diagnostics = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    membership_df = pd.concat(memberships, ignore_index=True) if memberships else pd.DataFrame()
    diagnostics.to_csv(OUT / "01_leave_one_dataset_out_module_state.csv", index=False)
    membership_df.to_csv(OUT / "02_module_membership_by_fold.csv", index=False)

    if diagnostics.empty:
        summary = pd.DataFrame()
    else:
        summary = diagnostics.groupby("held_out_dataset").agg(
            n_cases=("time_hours", "size"),
            mean_rmse=("cross_module_rmse", "mean"),
            mean_mae=("cross_module_mae", "mean"),
            mean_correlation=("cross_module_correlation", "mean"),
            mean_nearest_rmse=("nearest_module_rmse", "mean"),
            mean_rmse_improvement=("rmse_improvement_vs_nearest", "mean"),
        ).reset_index()
    summary.to_csv(OUT / "03_module_state_summary.csv", index=False)

    print("Stage 2.9 complete.")
    if not summary.empty:
        print(summary.to_string(index=False))
        print("\nStage 2.9 interpretation:")
        improvement = float(summary.mean_rmse_improvement.mean())
        corr = float(summary.mean_correlation.mean())
        if improvement > 0:
            print(f"  Module state beats nearest-time baseline on average (RMSE improvement={improvement:.4f}).")
            print(f"  Mean cross-dataset module-state correlation={corr:.3f}.")
            print("  This supports proceeding to a latent-state dynamics model, with uncertainty and controls kept explicit.")
        else:
            print(f"  Module state does not yet beat nearest-time baseline (mean RMSE improvement={improvement:.4f}).")
            print(f"  Mean cross-dataset module-state correlation={corr:.3f}.")
            print("  Do not fit an ODE yet; refine invariant biological programs first.")
    return summary


if __name__ == "__main__":
    stage2_9()
