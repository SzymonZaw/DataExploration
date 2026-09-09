"""Audit reproducibility versus branch/context dependence of Z6 trajectories.

This is a diagnostic following Z6 Phase 0.1. It does not alter the frozen
predictive benchmark. It asks whether the two independent branches within
each system exhibit a reproducible molecular transition.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamics.run_z6_within_system_audit import _build_system_data
from dynamics import validation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "z6_context_specificity_audit"
SYSTEMS = ("GSE28688", "GSE67462", "GSE297234")
SEED = 412
N_PERMUTATIONS = 1000


def _rank_corr(x, y):
    a = pd.Series(np.asarray(x, dtype=float)).rank(method="average").to_numpy()
    b = pd.Series(np.asarray(y, dtype=float)).rank(method="average").to_numpy()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _pearson(x, y):
    a = np.asarray(x, dtype=float)
    b = np.asarray(y, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _align(a, b):
    ta, xa = a
    tb, xb = b
    common = sorted(set(np.asarray(ta).tolist()) & set(np.asarray(tb).tolist()))
    if len(common) < 3:
        raise RuntimeError("Fewer than three shared timepoints between branches")
    ia = {float(t): i for i, t in enumerate(ta)}
    ib = {float(t): i for i, t in enumerate(tb)}
    xa2 = np.vstack([xa[ia[float(t)]] for t in common])
    xb2 = np.vstack([xb[ib[float(t)]] for t in common])
    return np.asarray(common, dtype=float), xa2, xb2


def _branch_metrics(dataset, branch_a, branch_b):
    times, xa, xb = _align(branch_a, branch_b)

    gene_temporal = np.asarray([_rank_corr(xa[:, j], xb[:, j]) for j in range(xa.shape[1])])
    temporal_abs = np.abs(gene_temporal[np.isfinite(gene_temporal)])

    delta_a = xa[-1] - xa[0]
    delta_b = xb[-1] - xb[0]
    finite_delta = np.isfinite(delta_a) & np.isfinite(delta_b)
    directional = np.sign(delta_a[finite_delta]) == np.sign(delta_b[finite_delta])
    nonzero = (delta_a[finite_delta] != 0) & (delta_b[finite_delta] != 0)
    directional_nonzero = directional[nonzero]

    endpoint_corr = _pearson(delta_a, delta_b)
    timepoint_rows = []
    for i, t in enumerate(times):
        timepoint_rows.append(
            {
                "dataset": dataset,
                "time_hours": float(t),
                "spearman_across_genes": _rank_corr(xa[i], xb[i]),
            }
        )

    return {
        "dataset": dataset,
        "n_genes": int(xa.shape[1]),
        "n_shared_timepoints": int(len(times)),
        "time_values": ",".join(map(str, times.tolist())),
        "median_abs_gene_temporal_spearman": float(np.median(temporal_abs)) if len(temporal_abs) else np.nan,
        "median_gene_temporal_spearman": float(np.median(gene_temporal[np.isfinite(gene_temporal)])) if np.isfinite(gene_temporal).any() else np.nan,
        "endpoint_directional_fraction": float(np.mean(directional)) if len(directional) else np.nan,
        "endpoint_directional_fraction_nonzero": float(np.mean(directional_nonzero)) if len(directional_nonzero) else np.nan,
        "endpoint_effect_pearson": endpoint_corr,
        "median_timepoint_spearman": float(np.nanmedian([r["spearman_across_genes"] for r in timepoint_rows])),
        "min_timepoint_spearman": float(np.nanmin([r["spearman_across_genes"] for r in timepoint_rows])),
        "timepoint_rows": timepoint_rows,
        "delta_a": delta_a,
        "delta_b": delta_b,
    }


def _permutation_null(delta_a, delta_b, n=N_PERMUTATIONS, seed=SEED):
    rng = np.random.default_rng(seed)
    observed = _pearson(delta_a, delta_b)
    vals = []
    for i in range(n):
        vals.append(_pearson(delta_a, rng.permutation(delta_b)))
    vals = np.asarray(vals, dtype=float)
    finite = vals[np.isfinite(vals)]
    p = (1 + np.sum(finite >= observed)) / (1 + len(finite)) if np.isfinite(observed) else np.nan
    return observed, vals, float(p)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    matrix, metadata = validation._load_common_space()

    temporal_rows = []
    direction_rows = []
    endpoint_rows = []
    timepoint_rows = []
    null_rows = []
    interpretation_rows = []

    rng_seed = SEED
    for dataset in SYSTEMS:
        data, inventory = _build_system_data(matrix, metadata, dataset)
        keys = sorted(data)
        if len(keys) != 2:
            raise RuntimeError(f"{dataset}: expected exactly two branches, got {keys}")

        metrics = _branch_metrics(dataset, data[keys[0]], data[keys[1]])
        for row in metrics["timepoint_rows"]:
            timepoint_rows.append(row)
        temporal_rows.append({k: v for k, v in metrics.items() if k not in {"timepoint_rows", "delta_a", "delta_b"}})
        direction_rows.append({
            "dataset": dataset,
            "branch_a": keys[0],
            "branch_b": keys[1],
            "endpoint_directional_fraction": metrics["endpoint_directional_fraction"],
            "endpoint_directional_fraction_nonzero": metrics["endpoint_directional_fraction_nonzero"],
        })
        endpoint_rows.append({
            "dataset": dataset,
            "branch_a": keys[0],
            "branch_b": keys[1],
            "endpoint_effect_pearson": metrics["endpoint_effect_pearson"],
        })

        observed, null, p = _permutation_null(metrics["delta_a"], metrics["delta_b"], seed=rng_seed)
        rng_seed += 1
        for i, value in enumerate(null):
            null_rows.append({"dataset": dataset, "permutation": i, "endpoint_effect_pearson": float(value)})

        strong = (
            metrics["median_abs_gene_temporal_spearman"] >= 0.60
            and metrics["endpoint_directional_fraction_nonzero"] >= 0.60
            and metrics["endpoint_effect_pearson"] >= 0.60
            and metrics["median_timepoint_spearman"] >= 0.60
        )
        weak = (
            metrics["median_abs_gene_temporal_spearman"] < 0.40
            or metrics["endpoint_directional_fraction_nonzero"] < 0.50
            or metrics["endpoint_effect_pearson"] < 0.40
        )
        if strong:
            label = "REPRODUCIBLE_TRANSITION"
        elif weak:
            label = "BRANCH_SPECIFIC_STRUCTURE"
        else:
            label = "MIXED_CONTEXT_DEPENDENCE"
        interpretation_rows.append({
            "dataset": dataset,
            "context_label": label,
            "predictive_support_phase0_1": dataset == "GSE67462",
            "positive_models_phase0_1": "autoencoder;PCA" if dataset == "GSE67462" else "",
            "endpoint_effect_permutation_p": p,
            "note": "Diagnostic only; does not alter frozen Z6 support decisions.",
        })

    pd.DataFrame(temporal_rows).to_csv(OUT / "01_branch_temporal_concordance.csv", index=False)
    pd.DataFrame(direction_rows).to_csv(OUT / "02_branch_directional_concordance.csv", index=False)
    pd.DataFrame(endpoint_rows).to_csv(OUT / "03_branch_endpoint_correlation.csv", index=False)
    pd.DataFrame(timepoint_rows).to_csv(OUT / "04_branch_timepoint_correlation.csv", index=False)
    pd.DataFrame(null_rows).to_csv(OUT / "05_permutation_endpoint_null.csv", index=False)
    pd.DataFrame(interpretation_rows).to_csv(OUT / "06_system_interpretation.csv", index=False)

    protocol = {
        "protocol": "Z6 context-specificity audit",
        "systems": list(SYSTEMS),
        "branch_definition": "identical to Z6 Phase 0.1",
        "predictive_support_changed": False,
        "thresholds_are_descriptive_only": True,
        "descriptive_thresholds": {
            "reproducible": "median abs gene temporal Spearman >= 0.60; nonzero endpoint direction >= 0.60; endpoint Pearson >= 0.60; median timepoint Spearman >= 0.60",
            "branch_specific": "median abs gene temporal Spearman < 0.40 OR nonzero endpoint direction < 0.50 OR endpoint Pearson < 0.40",
        },
        "permutation": {"n": N_PERMUTATIONS, "seed": SEED, "endpoint_effect": "shuffle gene identities in branch B"},
        "interpretation": "context reproducibility diagnostic, not biological specificity or causality",
    }
    (OUT / "07_protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")

    print("\nZ6 context-specificity audit complete.")
    print(pd.DataFrame(interpretation_rows).to_string(index=False))
    print(f"\nOutputs: {OUT}")


if __name__ == "__main__":
    main()
