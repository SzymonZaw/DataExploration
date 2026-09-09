"""Run Z6 Phase 0.1 within-system branch/donor holdout audit.

The frozen model, preprocessing and predictive-support rule are retained.
Only the validation unit changes from leave-one-dataset-out to leave-one-
branch-out within each biological system.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamics import model_benchmark as mb
from dynamics import validation
from dynamics.run_z6_phase0_audit import evaluate_pca_corrected

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "z6_within_system_audit"
SEEDS = (411, 412, 413, 414, 415)
SYSTEMS = ("GSE28688", "GSE67462", "GSE297234")


def _assign_gse28688_replicates(g: pd.DataFrame) -> pd.Series:
    """Recover the two experimental branches from the GEO sample order.

    GSE28688 records timed HFF1 samples as paired -a/-b measurements at each
    timepoint: GSM710513/514 (0 h), 515/516 (24 h), 517/518 (48 h), and
    519/520 (72 h). The common-space matrix preserves this GEO sample order,
    so even/odd matrix-column indices recover the same a/b pairing without
    relying on the legacy replicate field, which is unresolved for these GSM
    labels in the current validation metadata.
    """
    timed = g[g.time_hours.notna() & g.matrix_column.notna()].copy()
    parsed = pd.to_numeric(
        timed.matrix_column.astype(str).str.extract(r"__(\d+)$")[0],
        errors="coerce",
    )
    if parsed.isna().any():
        raise RuntimeError("GSE28688: could not parse common-space matrix column indices")

    # GEO order is paired a,b at every timed point; validate that invariant
    # before using it as the branch definition.
    tmp = timed.assign(_idx=parsed.astype(int))
    counts = tmp.groupby("time_hours")['_idx'].apply(
        lambda x: sorted(int(v) % 2 for v in x.tolist())
    )
    bad = counts[counts.map(lambda x: x != [0, 1])]
    if not bad.empty:
        raise RuntimeError(
            "GSE28688: expected exactly one odd and one even timed sample per "
            f"timepoint, got {bad.to_dict()}"
        )
    return pd.Series(
        np.where(parsed.to_numpy(dtype=int) % 2 == 0, "a", "b"),
        index=timed.index,
    ).reindex(g.index)


def _build_system_data(matrix: pd.DataFrame, metadata: pd.DataFrame, dataset: str):
    g = metadata[(metadata.dataset == dataset) & metadata.time_hours.notna() & metadata.matrix_column.notna()].copy()
    if dataset == "GSE297234":
        g["branch"] = g.condition.astype(str)
        branch_type = "donor"
    elif dataset == "GSE28688":
        # The legacy replicate field is unresolved for GSM710513-520.
        # Recover the documented GEO -a/-b branches from the preserved sample
        # order in the common-space matrix.
        g["branch"] = _assign_gse28688_replicates(g).astype(str)
        branch_type = "replicate"
    else:
        g["branch"] = g.replicate.astype(str)
        branch_type = "replicate"

    data = {}
    inventory = []
    for branch, b in sorted(g.groupby("branch")):
        if str(branch).lower() == "unknown":
            continue
        b = b.sort_values("time_hours")
        frame = pd.DataFrame(
            matrix[b.matrix_column.tolist()].T.to_numpy(dtype=float),
            index=b.time_hours.to_numpy(dtype=float),
            columns=matrix.index,
        ).groupby(level=0, sort=True).mean()
        if len(frame) < 4:
            continue
        key = f"{dataset}__{branch_type}_{branch}"
        data[key] = (frame.index.to_numpy(dtype=float), frame.to_numpy(dtype=float))
        inventory.append(
            {
                "dataset": dataset,
                "branch_type": branch_type,
                "branch": str(branch),
                "pseudo_dataset": key,
                "n_samples": len(b),
                "n_timepoints": len(frame),
                "time_values": ",".join(map(str, frame.index.tolist())),
            }
        )

    if len(data) < 2:
        raise RuntimeError(f"{dataset}: fewer than two eligible independent branches")
    return data, pd.DataFrame(inventory)


def _run_system(data, cfg, seeds):
    # Keep the corrected PCA implementation used by the Phase 0 audit.
    mb.evaluate_pca = evaluate_pca_corrected
    mb.permutation_null.__globals__["evaluate_pca"] = evaluate_pca_corrected
    return mb.benchmark(data, cfg, seeds)


def main(permutation_n: int = 1000):
    OUT.mkdir(parents=True, exist_ok=True)
    matrix, metadata = validation._load_common_space()

    cfg = mb.BenchmarkConfig(
        max_genes=2000,
        state_dim=8,
        hidden_dim=128,
        epochs=250,
        lr=1e-3,
        prefix_fraction=0.6,
        permutation_n=permutation_n,
    )

    all_rows = []
    all_summaries = []
    inventory_frames = []

    for dataset in SYSTEMS:
        data, inventory = _build_system_data(matrix, metadata, dataset)
        inventory_frames.append(inventory)
        df, summary = _run_system(data, cfg, SEEDS)
        df["system_dataset"] = dataset
        summary["system_dataset"] = dataset
        all_rows.append(df)
        all_summaries.append(summary)
        df.to_csv(OUT / f"fold_metrics_{dataset}.csv", index=False)
        summary.to_csv(OUT / f"model_summary_{dataset}.csv", index=False)

    folds = pd.concat(all_rows, ignore_index=True)
    summaries = pd.concat(all_summaries, ignore_index=True)
    inventory = pd.concat(inventory_frames, ignore_index=True)

    folds.to_csv(OUT / "01_fold_metrics_all_systems.csv", index=False)
    summaries.to_csv(OUT / "02_model_summary_by_system.csv", index=False)
    inventory.to_csv(OUT / "00_branch_inventory.csv", index=False)

    transfer = (
        summaries.groupby("model", as_index=False)
        .agg(
            n_systems=("system_dataset", "nunique"),
            systems_with_positive_mean=("mean_improvement_vs_persistence", lambda x: int((x > 0).sum())),
            systems_with_positive_q05=("q05_improvement_vs_persistence", lambda x: int((x > 0).sum())),
            systems_with_permutation_support=("permutation_p", lambda x: int((x < 0.05).sum())),
        )
    )
    support = summaries[
        [
            "system_dataset", "model", "predictive_support",
            "mean_improvement_vs_persistence", "q05_improvement_vs_persistence",
            "mean_improvement_vs_nearest", "q05_improvement_vs_nearest",
            "mean_improvement_vs_linear", "q05_improvement_vs_linear", "permutation_p",
        ]
    ].copy()
    support.to_csv(OUT / "03_support_decisions_by_system.csv", index=False)
    transfer.to_csv(OUT / "04_transferability_summary.csv", index=False)

    protocol = {
        "protocol": "Z6 Phase 0.1 within-system transferability audit",
        "systems": list(SYSTEMS),
        "validation": "leave-one-branch-out within each biological system",
        "GSE28688_branch": "GEO paired -a/-b timed samples, recovered from preserved sample order",
        "GSE67462_branch": "replicate 1/2",
        "GSE297234_branch": "donor aged/young; not technical replication",
        "seeds": list(SEEDS),
        "config": cfg.__dict__,
        "pca_fix": "pca.inverse_transform followed by scaler.inverse_transform",
        "forecast": "60% observed prefix followed by genuine free rollout",
        "baselines": ["persistence", "nearest-time", "linear"],
        "scientific_support_rule": "unchanged: positive mean and q05 improvement versus all baselines plus permutation p < 0.05",
        "interpretation": {
            "within_positive_lodo_negative": "supports H-Z6b: predictive information may exist but fail cross-dataset transfer",
            "within_negative_lodo_negative": "supports H-Z6a: current representation/dynamics lack robust predictive information even within systems",
        },
    }
    (OUT / "05_protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    print("\nZ6 Phase 0.1 within-system audit complete.")
    print("\nBranch inventory:")
    print(inventory.to_string(index=False))
    print("\nModel summary by system:")
    print(summaries.to_string(index=False))
    print("\nTransferability summary:")
    print(transfer.to_string(index=False))
    print(f"\nOutputs: {OUT}")


if __name__ == "__main__":
    main()
