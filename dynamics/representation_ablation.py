"""Phase 1: ablate the biological state representation.

The forecasting benchmark is kept fixed while the input representation changes:
raw common genes, PROGENy pathway activity and DoRothEA TF activity.
Prior-knowledge networks are fixed external resources; no representation
parameters are fit on the held-out dataset.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .validation import _load_common_space
from .model_benchmark import BenchmarkConfig, benchmark

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "phase1_representation_ablation"
DATASETS = ("GSE67462", "GSE28688", "GSE297234")
SEEDS = (511, 512, 513, 514, 515)


def _trajectory_data(matrix, metadata):
    data = {}
    for ds in DATASETS:
        g = metadata[
            (metadata.dataset == ds)
            & metadata.time_hours.notna()
            & metadata.matrix_column.notna()
        ].copy().sort_values("time_hours")
        if g.time_hours.nunique() < 3:
            continue
        X = matrix.loc[:, g.matrix_column.astype(str).tolist()].T.copy()
        X.index = g.time_hours.to_numpy(float)
        X = X.groupby(level=0, sort=True).mean()
        data[ds] = (X.index.to_numpy(float), X.to_numpy(float), list(X.columns))
    if len(data) != len(DATASETS):
        raise RuntimeError(f"Expected {len(DATASETS)} trajectory datasets, got {len(data)}")
    return data


def _score_network(data, net):
    import decoupler as dc

    scored = {}
    for ds, (times, X, _) in data.items():
        samples_by_gene = pd.DataFrame(X, columns=_.copy())
        samples_by_gene.index = [f"{ds}__{i}" for i in range(len(samples_by_gene))]
        acts, _ = dc.mt.ulm(data=samples_by_gene, net=net)
        acts = acts.apply(pd.to_numeric, errors="coerce")
        acts = acts.replace([np.inf, -np.inf], np.nan)
        scored[ds] = (times, acts.to_numpy(float), list(acts.columns))
    return scored


def _save_representation(data, name):
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    for ds, (times, X, features) in data.items():
        frame = pd.DataFrame(X, index=times, columns=features)
        frame.index.name = "time_hours"
        frame.to_csv(out / f"{ds}.csv")


def run(
    max_features=2000,
    state_dim=8,
    hidden_dim=128,
    epochs=250,
    lr=1e-3,
    prefix_fraction=0.6,
    permutation_n=1000,
    seeds=SEEDS,
):
    OUT.mkdir(parents=True, exist_ok=True)
    matrix, metadata = _load_common_space()
    gene_data_raw = _trajectory_data(matrix, metadata)
    gene_data = {ds: (t, X) for ds, (t, X, _) in gene_data_raw.items()}

    try:
        import decoupler as dc
    except ImportError as exc:
        raise RuntimeError(
            "Phase 1 requires decoupler. Install the AI stack with: "
            "pip install -r requirements-ai.txt"
        ) from exc

    progeny = dc.op.progeny(organism="human", top=100, license="academic")
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"], license="academic")

    pathway_data = {ds: (t, X, f) for ds, (t, X, f) in _score_network(gene_data_raw, progeny).items()}
    tf_data = {ds: (t, X, f) for ds, (t, X, f) in _score_network(gene_data_raw, dorothea).items()}
    pathway_data = {ds: (t, X) for ds, (t, X, _) in pathway_data.items()}
    tf_data = {ds: (t, X) for ds, (t, X, _) in tf_data.items()}

    _save_representation(gene_data, "genes")
    _save_representation(pathway_data, "progeny")
    _save_representation(tf_data, "dorothea")

    cfg = BenchmarkConfig(
        max_genes=max_features,
        state_dim=state_dim,
        hidden_dim=hidden_dim,
        epochs=epochs,
        lr=lr,
        prefix_fraction=prefix_fraction,
        permutation_n=permutation_n,
    )

    summaries = []
    metrics = []
    representations = {"genes": gene_data, "progeny": pathway_data, "dorothea": tf_data}
    for name, rep_data in representations.items():
        df, summary = benchmark(rep_data, cfg, seeds)
        df.insert(0, "representation", name)
        summary.insert(0, "representation", name)
        metrics.append(df)
        summaries.append(summary)
        df.to_csv(OUT / f"{name}_model_metrics.csv", index=False)
        summary.to_csv(OUT / f"{name}_model_summary.csv", index=False)

    all_metrics = pd.concat(metrics, ignore_index=True)
    all_summary = pd.concat(summaries, ignore_index=True)
    all_metrics.to_csv(OUT / "01_all_model_metrics.csv", index=False)
    all_summary.to_csv(OUT / "02_all_model_summary.csv", index=False)

    comparison = all_summary[
        [
            "representation",
            "model",
            "n_runs",
            "mean_rmse",
            "mean_improvement_vs_persistence",
            "q05_improvement_vs_persistence",
            "mean_improvement_vs_nearest",
            "q05_improvement_vs_nearest",
            "mean_improvement_vs_linear",
            "q05_improvement_vs_linear",
            "permutation_p",
            "predictive_support",
        ]
    ].copy()
    comparison.to_csv(OUT / "03_representation_comparison.csv", index=False)

    protocol = {
        "phase": "Phase 1 — biological state representation ablation",
        "representations": ["genes", "progeny", "dorothea"],
        "datasets": list(DATASETS),
        "split": "leave-one-dataset-out",
        "seeds": list(seeds),
        "forecast": "same 60% prefix and genuine free rollout as Phase 0",
        "baselines": ["persistence", "nearest-time", "linear"],
        "permutation": "same held-out-row permutation with fixed time vector",
        "progeny": "human PROGENy, top=100 target genes/pathway, academic license",
        "dorothea": "human DoRothEA confidence levels A/B/C, academic license",
        "representation_parameters_fit_on_holdout": False,
        "scientific_question": "Does biological prior knowledge yield a representation with stronger temporal generalization than the raw gene-level common space?",
    }
    (OUT / "00_protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")
    print("\nPhase 1 representation ablation:")
    print(comparison.to_string(index=False), flush=True)
    return comparison


if __name__ == "__main__":
    run()
