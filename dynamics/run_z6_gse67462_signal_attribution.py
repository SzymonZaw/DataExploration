"""Attribute the reproducible GSE67462 transition signal.

Diagnostic only: does not alter frozen Z6 predictive-support decisions.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from dynamics.run_z6_within_system_audit import _build_system_data
from dynamics import validation

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "z6_gse67462_signal_attribution"
DATASET = "GSE67462"
SEED = 412
N_PERMUTATIONS = 1000


def _zscore(v):
    v = np.asarray(v, dtype=float)
    mu = np.nanmean(v)
    sd = np.nanstd(v)
    return (v - mu) / sd if np.isfinite(sd) and sd > 1e-12 else np.zeros_like(v)


def _spearman(x, y):
    a = pd.Series(np.asarray(x, dtype=float)).rank(method="average").to_numpy()
    b = pd.Series(np.asarray(y, dtype=float)).rank(method="average").to_numpy()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _align(a, b):
    ta, xa = a
    tb, xb = b
    common = sorted(set(map(float, ta)) & set(map(float, tb)))
    ia = {float(t): i for i, t in enumerate(ta)}
    ib = {float(t): i for i, t in enumerate(tb)}
    xa2 = np.vstack([xa[ia[t]] for t in common])
    xb2 = np.vstack([xb[ib[t]] for t in common])
    return np.asarray(common), xa2, xb2


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    matrix, metadata = validation._load_common_space()
    data, inventory = _build_system_data(matrix, metadata, DATASET)
    keys = sorted(data)
    if len(keys) != 2:
        raise RuntimeError(f"Expected two GSE67462 branches, got {keys}")

    times, xa, xb = _align(data[keys[0]], data[keys[1]])
    # Center each gene on day 0: attribution concerns changes, not absolute branch level.
    da = xa - xa[0]
    db = xb - xb[0]
    common = (da + db) / 2.0
    branch = (da - db) / 2.0

    rows = []
    for g in range(common.shape[1]):
        c = common[:, g]
        b = branch[:, g]
        va = da[:, g]
        vb = db[:, g]
        common_energy = float(np.nanmean(c ** 2))
        branch_energy = float(np.nanmean(b ** 2))
        total_energy = float(np.nanmean((np.vstack([va, vb])) ** 2))
        common_fraction = common_energy / total_energy if total_energy > 0 else np.nan
        branch_fraction = branch_energy / total_energy if total_energy > 0 else np.nan
        rho_dynamic = _spearman(va, vb)
        rho_common_endpoint = _spearman(np.asarray([c[-1]]), np.asarray([b[-1]]))
        direction = np.sign(va[-1]) == np.sign(vb[-1]) and va[-1] != 0 and vb[-1] != 0
        rows.append({
            "gene_index": g,
            "common_energy": common_energy,
            "branch_energy": branch_energy,
            "total_change_energy": total_energy,
            "common_fraction": common_fraction,
            "branch_fraction": branch_fraction,
            "branch_dynamic_spearman": rho_dynamic,
            "endpoint_direction_concordant": bool(direction),
            "endpoint_abs_common": float(abs(c[-1])),
            "endpoint_abs_branch": float(abs(b[-1])),
        })

    genes = pd.DataFrame(rows)
    genes["attribution_ratio"] = genes["common_fraction"] / (genes["common_fraction"] + genes["branch_fraction"] + 1e-12)
    genes["common_dynamic_score"] = genes["attribution_ratio"] * genes["branch_dynamic_spearman"].clip(lower=0).fillna(0)
    genes["branch_context_score"] = genes["branch_fraction"] * (1 - genes["branch_dynamic_spearman"].clip(lower=0).fillna(0))
    genes["endpoint_only_score"] = genes["endpoint_abs_common"] / (genes["common_energy"] + 1e-12)

    # Descriptive categories; thresholds are frozen for this diagnostic and are not model-selection criteria.
    genes["component"] = np.select(
        [
            (genes["common_fraction"] >= 0.60) & (genes["branch_dynamic_spearman"] >= 0.60),
            (genes["endpoint_abs_common"] > 0) & (genes["branch_dynamic_spearman"] < 0.40),
            (genes["branch_fraction"] >= 0.60),
        ],
        ["COMMON_DYNAMIC_SIGNAL", "COMMON_ENDPOINT_SIGNAL", "BRANCH_CONTEXT_SIGNAL"],
        default="UNRESOLVED",
    )

    # Endpoint permutation: preserve each branch's temporal profile while destroying gene identity correspondence.
    rng = np.random.default_rng(SEED)
    observed = float(np.corrcoef(
        (da[-1] - da[0])[np.isfinite(da[-1]) & np.isfinite(db[-1])],
        (db[-1] - db[0])[np.isfinite(da[-1]) & np.isfinite(db[-1])],
    )[0, 1])
    null = []
    for _ in range(N_PERMUTATIONS):
        idx = rng.permutation(db.shape[1])
        x = da[-1]
        y = db[-1, idx]
        ok = np.isfinite(x) & np.isfinite(y)
        null.append(float(np.corrcoef(x[ok], y[ok])[0, 1]))
    null = np.asarray(null)
    p = float((1 + np.sum(null >= observed)) / (1 + len(null)))

    summary = {
        "dataset": DATASET,
        "branches": keys,
        "time_values": times.tolist(),
        "n_genes": int(genes.shape[0]),
        "median_common_fraction": float(genes["common_fraction"].median()),
        "median_branch_fraction": float(genes["branch_fraction"].median()),
        "fraction_common_dynamic": float(np.mean(genes["component"] == "COMMON_DYNAMIC_SIGNAL")),
        "fraction_common_endpoint": float(np.mean(genes["component"] == "COMMON_ENDPOINT_SIGNAL")),
        "fraction_branch_context": float(np.mean(genes["component"] == "BRANCH_CONTEXT_SIGNAL")),
        "fraction_unresolved": float(np.mean(genes["component"] == "UNRESOLVED")),
        "endpoint_observed_pearson": observed,
        "endpoint_permutation_p": p,
        "interpretation": "Associative attribution only; technical nuisance is not proven by this expression-only decomposition.",
    }

    genes.to_csv(OUT / "01_gene_level_attribution.csv", index=False)
    pd.DataFrame([summary]).to_csv(OUT / "02_system_summary.csv", index=False)
    pd.DataFrame({"permutation": np.arange(N_PERMUTATIONS), "endpoint_pearson": null}).to_csv(OUT / "03_endpoint_permutation_null.csv", index=False)
    pd.DataFrame({"time_hours": times, "median_common": np.nanmedian(common, axis=1), "median_branch_abs": np.nanmedian(np.abs(branch), axis=1)}).to_csv(OUT / "04_time_resolved_components.csv", index=False)
    (OUT / "05_protocol_snapshot.json").write_text(json.dumps({
        "dataset": DATASET,
        "branches": keys,
        "baseline": "day-0 centered expression trajectories",
        "diagnostic_thresholds": {"common_dynamic": "common_fraction >= 0.60 and branch_dynamic_spearman >= 0.60", "branch_context": "branch_fraction >= 0.60"},
        "permutation_n": N_PERMUTATIONS,
        "seed": SEED,
        "predictive_support_changed": False,
    }, indent=2), encoding="utf-8")

    print("\nGSE67462 signal attribution audit complete.")
    print(pd.DataFrame([summary]).to_string(index=False))
    print(f"\nOutputs: {OUT}")


if __name__ == "__main__":
    main()
