"""First end-to-end DynamicStateModel experiment on the validated Stage 2.6 space."""
from __future__ import annotations
import argparse
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from .dynamic_state_model import DynamicStateModel

ROOT = Path(__file__).resolve().parents[1]
COMMON_MATRIX = ROOT / "results" / "Dynamics" / "stage2_6" / "06_common_human_gene_matrix.csv"
COMMON_METADATA = ROOT / "results" / "Dynamics" / "stage2_6" / "07_common_gene_sample_metadata.csv"
OUT = ROOT / "results" / "Dynamics" / "dynamic_state_model"
DATASETS = ["GSE67462", "GSE28688", "GSE297234"]


def set_seed(seed: int) -> None:
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)


def load_data() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Load Stage 2.6 data through the validated Stage 2.7 metadata resolver.

    Stage 2.6 does not store time_hours in its compact metadata. Stage 2.7
    reconstructs times from the authoritative dataset metadata and resolves
    legacy matrix-column mappings, so reuse that logic here rather than
    maintaining a second, potentially inconsistent time resolver.
    """
    if not COMMON_MATRIX.exists() or not COMMON_METADATA.exists():
        raise FileNotFoundError("Common-space files missing; run python validate_pipeline.py --refresh first.")

    from .validation import _load_common_space
    matrix, meta = _load_common_space()
    required = {"dataset", "sample", "matrix_column", "time_hours"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"Resolved metadata missing columns: {sorted(missing)}")

    data = {}
    for ds in DATASETS:
        g = meta[meta["dataset"].astype(str).eq(ds)].copy()
        g["time_hours"] = pd.to_numeric(g["time_hours"], errors="coerce")
        g = g[g["time_hours"].notna() & g["matrix_column"].notna()].copy()
        if g["time_hours"].nunique() < 3:
            continue
        cols = g["matrix_column"].astype(str).tolist()
        X = matrix.loc[:, cols].T.copy()
        X.index = g["time_hours"].to_numpy(float)
        X = X.groupby(level=0, sort=True).mean()
        data[ds] = (X.index.to_numpy(float), X.to_numpy(float))

    if len(data) < 2:
        raise RuntimeError("At least two longitudinal datasets are required.")
    return data


def training_statistics(train, max_genes):
    """Fit feature selection and imputation strictly on training datasets."""
    arrays = [X for _, X in train.values()]
    finite_values = np.concatenate([np.where(np.isfinite(X), X, np.nan) for X in arrays], axis=0)
    med = np.nanmedian(finite_values, axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    pooled = np.vstack([np.where(np.isfinite(X), X, med) for X in arrays])
    var = np.var(pooled, axis=0)
    keep = np.argsort(var)[::-1][:min(max_genes, len(var))]
    mean = np.mean(pooled[:, keep], axis=0)
    std = np.std(pooled[:, keep], axis=0)
    std = np.where(std > 1e-8, std, 1.0)
    return keep, mean, std


def transform(X, keep, mean, std):
    X = X[:, keep]
    X = np.where(np.isfinite(X), X, mean)
    return (X - mean) / std


def pairs(data, keep, mean, std):
    xs, ys = [], []
    for t, X in data.values():
        Z = transform(X, keep, mean, std)
        for i in range(len(t) - 1):
            xs.append(Z[i]); ys.append(Z[i + 1])
    return torch.tensor(np.asarray(xs), dtype=torch.float32), torch.tensor(np.asarray(ys), dtype=torch.float32)


def train(train_data, keep, mean, std, state_dim, hidden_dim, epochs, lr, seed):
    set_seed(seed)
    x, y = pairs(train_data, keep, mean, std)
    model = DynamicStateModel(input_dim=x.shape[1], state_dim=state_dim, hidden_dim=hidden_dim, dropout=0.05)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    best = float("inf")
    for _ in range(epochs):
        model.train(); optimizer.zero_grad(set_to_none=True)
        z = model.encode(x)
        z_next = model.encode(y).detach()
        reconstruction = model.decode(z)
        predicted = model.transition(z)
        predicted_obs = model.decode(predicted)
        loss = (torch.mean((reconstruction - x) ** 2) +
                0.25 * torch.mean((predicted - z_next) ** 2) +
                torch.mean((predicted_obs - y) ** 2))
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        best = min(best, float(loss.detach()))
    return model, best


def evaluate(model, heldout, keep, mean, std):
    _, X = heldout
    Z = transform(X, keep, mean, std)
    x = torch.tensor(Z[:-1], dtype=torch.float32)
    y = torch.tensor(Z[1:], dtype=torch.float32)
    model.eval()
    with torch.no_grad():
        pred = model.predict_observation(x)
        recon = model.decode(model.encode(x))
    pred_np, recon_np, y_np = pred.numpy(), recon.numpy(), y.numpy()
    rmse_model = float(np.sqrt(np.mean((pred_np - y_np) ** 2)))
    rmse_persistence = float(np.sqrt(np.mean((x.numpy() - y_np) ** 2)))
    rmse_reconstruction = float(np.sqrt(np.mean((recon_np - y_np) ** 2)))
    return {"n_test_transitions": len(y_np), "rmse_model": rmse_model,
            "rmse_persistence": rmse_persistence,
            "rmse_reconstruction_to_next": rmse_reconstruction,
            "improvement_vs_persistence": rmse_persistence - rmse_model,
            "mean_abs_error_model": float(np.mean(np.abs(pred_np - y_np)))}


def run(max_genes=2000, state_dim=8, hidden_dim=128, epochs=300, lr=1e-3, seed=211):
    OUT.mkdir(parents=True, exist_ok=True)
    data = load_data(); rows = []
    for heldout_name in sorted(data):
        train_data = {k: v for k, v in data.items() if k != heldout_name}
        if len(train_data) < 2:
            print(f"DynamicStateModel: skipping {heldout_name}; only {len(train_data)} training datasets", flush=True); continue
        keep, mean, std = training_statistics(train_data, max_genes)
        model, train_loss = train(train_data, keep, mean, std, state_dim, hidden_dim, epochs, lr, seed)
        metrics = evaluate(model, data[heldout_name], keep, mean, std)
        rows.append({"heldout_dataset": heldout_name, "n_training_datasets": len(train_data),
                     "n_selected_genes": len(keep), "state_dim": state_dim,
                     "train_loss": train_loss, **metrics})
        print(f"DynamicStateModel LODO: {heldout_name} -> RMSE {metrics['rmse_model']:.4f}, persistence {metrics['rmse_persistence']:.4f}", flush=True)
    result = pd.DataFrame(rows); result.to_csv(OUT / "01_lodo_metrics.csv", index=False)
    if not result.empty:
        summary = pd.DataFrame([{"n_valid_lodo_folds": len(result),
                                 "mean_rmse_model": result.rmse_model.mean(),
                                 "mean_rmse_persistence": result.rmse_persistence.mean(),
                                 "mean_improvement_vs_persistence": result.improvement_vs_persistence.mean(),
                                 "dynamic_state_predictive_support": bool(result.improvement_vs_persistence.mean() > 0)}])
        summary.to_csv(OUT / "02_summary.csv", index=False); print(summary.to_string(index=False), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-genes", type=int, default=2000)
    parser.add_argument("--state-dim", type=int, default=8)
    parser.add_argument("--hidden-dim", type=int, default=128)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--seed", type=int, default=211)
    args = parser.parse_args(); run(args.max_genes, args.state_dim, args.hidden_dim, args.epochs, args.lr, args.seed)
