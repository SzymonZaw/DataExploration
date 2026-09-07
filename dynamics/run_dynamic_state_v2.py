"""Validated DynamicStateModel v2 experiment.

This is the first serious predictive test of the learned latent dynamics.
It performs leave-one-dataset-out evaluation with:
- training-only feature selection, imputation and scaling,
- explicit delta-t conditioning,
- prefix-to-future forecasting on the held-out dataset,
- persistence, nearest-time and linear baselines,
- multiple random seeds and a time-permutation null.

No biological or causal claim is made by this module.
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch import nn

from .dynamic_state_model import DynamicStateModel
from .validation import _load_common_space

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "dynamic_state_model_v2"
DATASETS = ["GSE67462", "GSE28688", "GSE297234"]
DEFAULT_SEEDS = (211, 212, 213, 214, 215)


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def load_data() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    matrix, meta = _load_common_space()
    data: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for ds in DATASETS:
        g = meta[(meta["dataset"] == ds) & meta["time_hours"].notna() & meta["matrix_column"].notna()].copy()
        if g["time_hours"].nunique() < 3:
            continue
        g = g.sort_values("time_hours")
        X = matrix.loc[:, g["matrix_column"].astype(str).tolist()].T
        X.index = g["time_hours"].to_numpy(float)
        X = X.groupby(level=0, sort=True).mean()
        data[ds] = (X.index.to_numpy(float), X.to_numpy(float))
    if len(data) < 2:
        raise RuntimeError("At least two longitudinal datasets are required.")
    return data


def training_statistics(train, max_genes: int):
    arrays = [X for _, X in train.values()]
    pooled_raw = np.vstack(arrays)
    med = np.nanmedian(np.where(np.isfinite(pooled_raw), pooled_raw, np.nan), axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    pooled = np.vstack([np.where(np.isfinite(X), X, med) for X in arrays])
    var = np.var(pooled, axis=0)
    keep = np.argsort(var)[::-1][: min(max_genes, pooled.shape[1])]
    mean = pooled[:, keep].mean(axis=0)
    std = pooled[:, keep].std(axis=0)
    std = np.where(std > 1e-8, std, 1.0)
    return keep, mean, std


def transform(X, keep, mean, std):
    X = X[:, keep]
    X = np.where(np.isfinite(X), X, mean)
    return (X - mean) / std


def transition_pairs(train, keep, mean, std):
    xs, ys, dts = [], [], []
    for times, X in train.values():
        Z = transform(X, keep, mean, std)
        scale = max(float(times[-1] - times[0]), 1.0)
        tn = (times - times[0]) / scale
        for i in range(len(times) - 1):
            xs.append(Z[i]); ys.append(Z[i + 1]); dts.append(float(tn[i + 1] - tn[i]))
    return (torch.tensor(np.asarray(xs), dtype=torch.float32),
            torch.tensor(np.asarray(ys), dtype=torch.float32),
            torch.tensor(np.asarray(dts)[:, None], dtype=torch.float32))


def build_dt_model(input_dim, state_dim, hidden_dim):
    return DynamicStateModel(input_dim=input_dim, state_dim=state_dim,
                             hidden_dim=hidden_dim, context_dim=1, dropout=0.05)


def train_model(train_data, keep, mean, std, state_dim, hidden_dim, epochs, lr, seed):
    set_seed(seed)
    x, y, dt = transition_pairs(train_data, keep, mean, std)
    model = build_dt_model(x.shape[1], state_dim, hidden_dim)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    best = float("inf")
    patience = max(20, epochs // 8)
    stale = 0
    for _ in range(epochs):
        model.train(); optimizer.zero_grad(set_to_none=True)
        z = model.encode(x)
        z_next = model.encode(y).detach()
        reconstruction = model.decode(z)
        predicted = model.transition(z, context=dt)
        predicted_obs = model.decode(predicted)
        loss = (torch.mean((reconstruction - x) ** 2)
                + 0.25 * torch.mean((predicted - z_next) ** 2)
                + torch.mean((predicted_obs - y) ** 2))
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        value = float(loss.detach())
        if value < best - 1e-6:
            best = value; stale = 0
        else:
            stale += 1
            if stale >= patience:
                break
    return model, best


def predict_step(model, current, dt):
    x = torch.tensor(current[None, :], dtype=torch.float32)
    d = torch.tensor([[float(dt)]], dtype=torch.float32)
    model.eval()
    with torch.no_grad():
        return model.predict_observation(x, context=d).numpy()[0]


def linear_prediction(prefix_t, prefix_X, target_t):
    if len(prefix_t) < 2:
        return prefix_X[-1].copy()
    t0, t1 = float(prefix_t[-2]), float(prefix_t[-1])
    if t1 == t0:
        return prefix_X[-1].copy()
    alpha = (float(target_t) - t1) / (t1 - t0)
    return prefix_X[-1] + alpha * (prefix_X[-1] - prefix_X[-2])


def nearest_prediction(train_data, target_t):
    candidates = []
    for times, X in train_data.values():
        idx = int(np.argmin(np.abs(times - target_t)))
        candidates.append(X[idx])
    return np.mean(candidates, axis=0)


def metrics(y, p):
    e = np.asarray(y) - np.asarray(p)
    return float(np.sqrt(np.mean(e ** 2))), float(np.mean(np.abs(e)))


def evaluate_prefix(model, train_data, heldout, keep, mean, std, prefix_fraction):
    times, Xraw = heldout
    X = transform(Xraw, keep, mean, std)
    n_prefix = max(2, int(np.ceil(len(times) * prefix_fraction)))
    n_prefix = min(n_prefix, len(times) - 1)
    prefix_t, prefix_X = times[:n_prefix], X[:n_prefix]
    preds_model, preds_persist, preds_nearest, preds_linear, true = [], [], [], [], []
    scale = max(float(times[-1] - times[0]), 1.0)
    for j in range(n_prefix, len(times)):
        dt = float((times[j] - times[j - 1]) / scale)
        preds_model.append(predict_step(model, prefix_X[j - 1], dt))
        preds_persist.append(prefix_X[j - 1])
        preds_nearest.append(nearest_prediction(train_data, times[j]))
        preds_linear.append(linear_prediction(prefix_t, prefix_X, times[j]))
        true.append(X[j])
        prefix_t = np.append(prefix_t, times[j])
        prefix_X = np.vstack([prefix_X, X[j]])
    if not true:
        return None
    out = {"n_future_points": len(true)}
    for name, pred in (("model", preds_model), ("persistence", preds_persist),
                       ("nearest", preds_nearest), ("linear", preds_linear)):
        rmse, mae = metrics(np.asarray(true), np.asarray(pred))
        out[f"rmse_{name}"] = rmse; out[f"mae_{name}"] = mae
    out["improvement_vs_persistence"] = out["rmse_persistence"] - out["rmse_model"]
    out["improvement_vs_nearest"] = out["rmse_nearest"] - out["rmse_model"]
    out["improvement_vs_linear"] = out["rmse_linear"] - out["rmse_model"]
    return out


def permutation_null(data, max_genes, state_dim, hidden_dim, epochs, lr, seed, prefix_fraction, n_perm):
    rng = np.random.default_rng(seed)
    observed = []
    for _ in range(n_perm):
        permuted = {}
        for ds, (t, X) in data.items():
            order = rng.permutation(len(t)); permuted[ds] = (t, X[order])
        vals = []
        for held_name in sorted(permuted):
            train = {k: v for k, v in permuted.items() if k != held_name}
            if len(train) < 2: continue
            keep, mean, std = training_statistics(train, max_genes)
            model, _ = train_model(train, keep, mean, std, state_dim, hidden_dim, max(20, epochs // 3), lr, seed + 1000)
            r = evaluate_prefix(model, train, permuted[held_name], keep, mean, std, prefix_fraction)
            if r: vals.append(r["improvement_vs_persistence"])
        if vals: observed.append(float(np.mean(vals)))
    return np.asarray(observed)


def run(max_genes=2000, state_dim=8, hidden_dim=128, epochs=250, lr=1e-3,
        seeds=DEFAULT_SEEDS, prefix_fraction=0.6, n_perm=1000):
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_data(); rows = []
    for seed in seeds:
        for held_name in sorted(data):
            train = {k: v for k, v in data.items() if k != held_name}
            if len(train) < 2: continue
            keep, mean, std = training_statistics(train, max_genes)
            model, loss = train_model(train, keep, mean, std, state_dim, hidden_dim, epochs, lr, seed)
            result = evaluate_prefix(model, train, data[held_name], keep, mean, std, prefix_fraction)
            if result is None: continue
            rows.append({"seed": seed, "heldout_dataset": held_name,
                         "n_training_datasets": len(train), "n_selected_genes": len(keep),
                         "train_loss": loss, **result})
            print(f"DynamicStateModel v2: seed={seed} heldout={held_name} "
                  f"RMSE={result['rmse_model']:.4f} persistence={result['rmse_persistence']:.4f}", flush=True)
    result_df = pd.DataFrame(rows)
    result_df.to_csv(OUT / "01_lodo_prefix_metrics.csv", index=False)
    if result_df.empty:
        raise RuntimeError("No valid LODO prefix-forecast folds.")
    summary = pd.DataFrame([{
        "n_valid_folds": len(result_df),
        "n_seeds": result_df["seed"].nunique(),
        "mean_rmse_model": result_df.rmse_model.mean(),
        "mean_rmse_persistence": result_df.rmse_persistence.mean(),
        "mean_rmse_nearest": result_df.rmse_nearest.mean(),
        "mean_rmse_linear": result_df.rmse_linear.mean(),
        "mean_improvement_vs_persistence": result_df.improvement_vs_persistence.mean(),
        "mean_improvement_vs_nearest": result_df.improvement_vs_nearest.mean(),
        "mean_improvement_vs_linear": result_df.improvement_vs_linear.mean(),
        "median_improvement_vs_persistence": result_df.improvement_vs_persistence.median(),
        "q05_improvement_vs_persistence": result_df.improvement_vs_persistence.quantile(.05),
        "q95_improvement_vs_persistence": result_df.improvement_vs_persistence.quantile(.95),
    }])
    null = permutation_null(data, max_genes, state_dim, hidden_dim, epochs, lr, seeds[0], prefix_fraction, n_perm)
    obs = float(summary.loc[0, "mean_improvement_vs_persistence"])
    p = float((1 + np.sum(null >= obs)) / (1 + len(null))) if len(null) else np.nan
    summary["permutation_p_improvement_vs_persistence"] = p
    summary["dynamic_state_predictive_support"] = bool(
        summary.loc[0, "mean_improvement_vs_persistence"] > 0 and
        summary.loc[0, "mean_improvement_vs_nearest"] > 0 and
        summary.loc[0, "mean_improvement_vs_linear"] > 0 and
        p < 0.05
    )
    summary.to_csv(OUT / "02_summary.csv", index=False)
    pd.DataFrame({"null_improvement_vs_persistence": null}).to_csv(OUT / "03_permutation_null.csv", index=False)
    print(summary.to_string(index=False), flush=True)
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--max-genes", type=int, default=2000)
    p.add_argument("--state-dim", type=int, default=8)
    p.add_argument("--hidden-dim", type=int, default=128)
    p.add_argument("--epochs", type=int, default=250)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--prefix-fraction", type=float, default=.6)
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    a = p.parse_args()
    run(a.max_genes, a.state_dim, a.hidden_dim, a.epochs, a.lr, tuple(a.seeds), a.prefix_fraction, a.permutations)
