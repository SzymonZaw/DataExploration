"""Phase 1: ablate the biological state representation.

The forecasting benchmark is kept fixed while the input representation changes:
raw common genes, PROGENy pathway activity and DoRothEA TF activity.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from .validation import _load_common_space
from .model_benchmark import BenchmarkConfig, benchmark

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "phase1_representation_ablation"
DATASETS = ("GSE67462", "GSE28688", "GSE297234")
SEEDS = (511, 512, 513, 514, 515)


def _finite_frame(frame):
    """Coerce to finite floats, replacing non-finite values with zero."""
    frame = frame.apply(pd.to_numeric, errors="coerce").astype(float)
    arr = frame.to_numpy(copy=True)
    bad = ~np.isfinite(arr)
    report = {"nan_count": int(np.isnan(arr).sum()), "inf_count": int(np.isinf(arr).sum()), "nonfinite_count": int(bad.sum())}
    if bad.any():
        arr[bad] = 0.0
        frame = pd.DataFrame(arr, index=frame.index, columns=frame.columns)
    return frame, report


def _get_gene_data(matrix, metadata):
    data = {}
    for ds in DATASETS:
        g = metadata[(metadata.dataset == ds) & metadata.time_hours.notna() & metadata.matrix_column.notna()].copy().sort_values("time_hours")
        if g.time_hours.nunique() < 3:
            continue
        columns = g.matrix_column.astype(str).tolist()
        X = matrix.loc[:, columns].T.copy()
        X.index = g.time_hours.to_numpy(float)
        X = X.groupby(level=0, sort=True).mean()
        X, _ = _finite_frame(X)
        data[ds] = (X.index.to_numpy(float), X.to_numpy(float), list(X.columns))
    return data


def _score_network(data, net):
    import decoupler as dc
    scored = {}
    audits = []
    for ds, (times, X, genes) in data.items():
        samples_by_gene = pd.DataFrame(X, columns=genes)
        samples_by_gene.index = [f"{ds}__{i}" for i in range(len(samples_by_gene))]
        samples_by_gene, report = _finite_frame(samples_by_gene)
        report.update({"dataset": ds, "stage": "input", "n_samples": len(samples_by_gene), "n_genes": len(samples_by_gene.columns)})
        audits.append(report)
        acts, _ = dc.mt.ulm(data=samples_by_gene, net=net)
        acts, score_report = _finite_frame(acts)
        score_report.update({"dataset": ds, "stage": "activity", "n_activity_features": len(acts.columns)})
        audits.append(score_report)
        scored[ds] = (times, acts)

    # decoupler can return a slightly different set of activities for different
    # datasets when a network has insufficient measured targets. The frozen
    # benchmark requires one common feature space across all datasets. Use the
    # intersection rather than padding missing activities with arbitrary zeros.
    common_features = None
    for _, acts in scored.values():
        features = set(acts.columns)
        common_features = features if common_features is None else common_features & features
    common_features = sorted(common_features or [])
    if len(common_features) < 3:
        raise ValueError(
            f"Representation has only {len(common_features)} features shared across "
            f"{len(scored)} datasets; cannot run the common-space benchmark."
        )

    aligned = {}
    for ds, (times, acts) in scored.items():
        acts = acts.loc[:, common_features]
        aligned[ds] = (times, acts.to_numpy(float))
    for report in audits:
        if report.get("stage") == "activity":
            report["n_common_activity_features"] = len(common_features)
    audits.append({"stage": "alignment", "n_common_activity_features": len(common_features), "datasets": len(aligned)})

    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audits).to_csv(OUT / "01_representation_scoring_audit.csv", index=False)
    return aligned


def _get_prior_knowledge():
    import decoupler as dc
    progeny = dc.op.progeny(organism="human", top=100)
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
    return progeny, dorothea


def _tag(result, representation):
    out = result.copy()
    out["representation"] = representation
    return out


def run():
    loaded = _load_common_space()
    matrix, metadata = loaded[:2]
    gene_data_raw = _get_gene_data(matrix, metadata)
    progeny, dorothea = _get_prior_knowledge()
    pathway_data = _score_network(gene_data_raw, progeny)
    tf_data = _score_network(gene_data_raw, dorothea)
    # model_benchmark expects exactly (time, X). Feature names are only needed
    # for the decoupler step above.
    gene_data = {ds: (t, X) for ds, (t, X, _) in gene_data_raw.items()}

    cfg = BenchmarkConfig(max_genes=2000, state_dim=8, hidden_dim=128, epochs=250, lr=1e-3, prefix_fraction=0.6, history_len=2, history_dim=16, dropout=0.05, permutation_n=1000, time_scale_hours=168.0)
    representations = (("genes", gene_data), ("PROGENy", pathway_data), ("DoRothEA", tf_data))
    detail_frames, summary_frames = [], []
    for representation, rep_data in representations:
        detail, summary = benchmark(rep_data, cfg=cfg, seeds=SEEDS)
        detail_frames.append(_tag(detail, representation))
        summary_frames.append(_tag(summary, representation))
    detail = pd.concat(detail_frames, ignore_index=True)
    summary = pd.concat(summary_frames, ignore_index=True)
    OUT.mkdir(parents=True, exist_ok=True)
    detail.to_csv(OUT / "02_phase1_results.csv", index=False)
    summary.to_csv(OUT / "03_phase1_summary.csv", index=False)
    return summary
