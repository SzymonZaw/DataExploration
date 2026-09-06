"""Leakage-free Stage 2.4 leave-one-timepoint-out validation.

PCA and feature scaling are fitted only on training timepoints in every fold.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_4"
OUT.mkdir(parents=True, exist_ok=True)


def _fit_transform(train_samples, all_samples, n_components=3, max_features=5000):
    x_train = train_samples.T.apply(pd.to_numeric, errors="coerce")
    x_all = all_samples.T.apply(pd.to_numeric, errors="coerce")
    common = x_train.columns[x_train.notna().sum() >= 2]
    x_train, x_all = x_train[common], x_all[common]
    variances = x_train.var(axis=0, ddof=0).sort_values(ascending=False)
    keep = variances.head(min(max_features, len(variances))).index
    x_train = x_train[keep].fillna(x_train[keep].mean())
    x_all = x_all[keep].fillna(x_train.mean())
    scaler = StandardScaler().fit(x_train)
    z_train = scaler.transform(x_train)
    z_all = scaler.transform(x_all)
    n = min(n_components, z_train.shape[0] - 1, z_train.shape[1])
    if n < 1:
        raise ValueError("Not enough training samples/features for PCA")
    pca = PCA(n_components=n, random_state=42).fit(z_train)
    return pd.DataFrame(pca.transform(z_all), index=x_all.index, columns=[f"PC{i+1}" for i in range(n)])


def leave_one_timepoint_out(expression, sample_times, branch=None):
    """Run leakage-free PCA trajectory interpolation for one dataset."""
    times = pd.Series(sample_times, dtype=float).reindex(expression.columns)
    keep = times.notna()
    if branch is not None:
        keep &= pd.Series(branch, index=expression.columns).astype(bool)
    expression, times = expression.loc[:, keep], times.loc[expression.columns]
    unique_times = sorted(times.unique())
    rows = []
    for held in unique_times:
        train_cols, held_cols = times.index[times != held], times.index[times == held]
        if len(train_cols) < 3 or len(unique_times) < 3:
            continue
        coords = _fit_transform(expression[train_cols], expression)
        train_coords = coords.loc[train_cols]
        held_coords = coords.loc[held_cols].mean(axis=0).to_numpy()
        train_mean = train_coords.groupby(times.loc[train_cols]).mean().sort_index()
        pred = np.asarray([np.interp(held, train_mean.index.to_numpy(float), train_mean[c].to_numpy(float)) for c in train_mean.columns])
        rows.append({"held_out_time_hours":float(held), "n_train_samples":len(train_cols), "n_train_timepoints":len(unique_times)-1, "oos_pca_error":float(np.linalg.norm(held_coords-pred))})
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "02_leakage_free_leave_one_timepoint_out.csv", index=False)
    return result


__all__ = ["leave_one_timepoint_out"]
