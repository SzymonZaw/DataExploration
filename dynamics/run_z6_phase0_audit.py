"""Audit Z6 Phase 0 without changing the frozen scientific decision rule.

This wrapper patches the legacy PCA evaluator so predictions are mapped back
through both PCA and StandardScaler into the same training-standardized
observation space used by the truth and baselines. It then runs the existing
LODO benchmark and writes fold-level and dataset-level diagnostics.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

from dynamics import model_benchmark as mb
from dynamics import run_model_benchmark as runner

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "z6_predictive_transition_audit"


def evaluate_pca_corrected(fitted, train_data, heldout_data, keep, mean, std, cfg, order=None):
    scaler, pca, transition = fitted
    times, X, n_prefix = mb._prepare_holdout(heldout_data, keep, mean, std, cfg, order)
    current_z = pca.transform(scaler.transform(X[n_prefix - 1][None, :]))[0]
    pred, true = [], []
    for j in range(n_prefix, len(times)):
        current_z = transition.predict(current_z[None, :])[0]
        # X is in the training-statistics standardized space. PCA was fitted
        # after a second StandardScaler, so both transforms must be inverted.
        x_scaled = pca.inverse_transform(current_z[None, :])
        x_benchmark = scaler.inverse_transform(x_scaled)[0]
        pred.append(x_benchmark)
        true.append(X[j])
    persistence, nearest, linear = mb._static_baselines(
        train_data, times, X, n_prefix, keep, mean, std
    )
    return mb._metrics(true, pred, persistence, nearest, linear)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    mb.evaluate_pca = evaluate_pca_corrected
    mb.permutation_null.__globals__["evaluate_pca"] = evaluate_pca_corrected
    data = runner.load_data()
    cfg = mb.BenchmarkConfig(
        max_genes=2000,
        state_dim=8,
        hidden_dim=128,
        epochs=250,
        lr=1e-3,
        prefix_fraction=.6,
        permutation_n=1000,
    )
    seeds = runner.SEEDS
    df, summary = mb.benchmark(data, cfg, seeds)

    df.to_csv(OUT / "01_fold_metrics.csv", index=False)
    summary.to_csv(OUT / "02_model_summary.csv", index=False)

    dataset_summary = (
        df.groupby(["heldout_dataset", "model"], as_index=False)
        .agg(
            n_runs=("rmse_model", "size"),
            mean_future_points=("n_future_points", "mean"),
            min_future_points=("n_future_points", "min"),
            mean_rmse=("rmse_model", "mean"),
            mean_improvement_vs_persistence=("improvement_vs_persistence", "mean"),
            q05_improvement_vs_persistence=("improvement_vs_persistence", lambda x: x.quantile(.05)),
            mean_improvement_vs_nearest=("improvement_vs_nearest", "mean"),
            q05_improvement_vs_nearest=("improvement_vs_nearest", lambda x: x.quantile(.05)),
            mean_improvement_vs_linear=("improvement_vs_linear", "mean"),
            q05_improvement_vs_linear=("improvement_vs_linear", lambda x: x.quantile(.05)),
        )
    )
    dataset_summary.to_csv(OUT / "03_dataset_model_summary.csv", index=False)

    fold = df[
        [
            "seed", "heldout_dataset", "model", "n_future_points",
            "rmse_model", "rmse_persistence", "rmse_nearest", "rmse_linear",
            "improvement_vs_persistence", "improvement_vs_nearest",
            "improvement_vs_linear", "train_loss",
        ]
    ].sort_values(["heldout_dataset", "model", "seed"])
    fold.to_csv(OUT / "04_fold_diagnostics.csv", index=False)

    protocol = pd.Series({
        "audit": "Z6 Phase 0 implementation and fold-level audit",
        "pca_fix": "pca.inverse_transform followed by scaler.inverse_transform",
        "split": "leave-one-dataset-out",
        "datasets": ",".join(sorted(data)),
        "seeds": ",".join(map(str, seeds)),
        "prefix_fraction": cfg.prefix_fraction,
        "permutation_n": cfg.permutation_n,
        "scientific_support_rule": "unchanged: positive mean and q05 improvement versus all baselines plus permutation p < 0.05",
    })
    protocol.to_json(OUT / "00_audit_protocol.json", indent=2)

    print("\nZ6 Phase 0 audit complete.")
    print("Corrected PCA evaluation: pca.inverse_transform -> scaler.inverse_transform")
    print("\nModel summary:")
    print(summary.to_string(index=False))
    print("\nDataset × model summary:")
    print(dataset_summary.to_string(index=False))
    print(f"\nAudit outputs: {OUT}")


if __name__ == "__main__":
    main()
