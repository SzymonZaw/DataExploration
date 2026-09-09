from __future__ import annotations

"""Calibrate H3 structural eligibility metrics on the Yamanaka source data.

This is a calibration/audit tool, not an H3 similarity test. It computes the
same structural quantities used for prospective controls so we can determine
whether those quantities are commensurate with the target process. It never
changes thresholds automatically.

Example:
  python dynamics/h3_yamanaka_structural_calibration.py \
      --expression Data/GSE297234_expression.tsv.gz \
      --sample-days "GM00731_D0:0,GM00731_D3:3,GM00731_D7:7,GM00731_D10:10,HFIB_D0:0,HFIB_D3:3,HFIB_D7:7,HFIB_D10:10"

The expression matrix must be gene x sample, with the first column containing
gene identifiers. Multiple biological trajectories may be supplied; values are
averaged within each native timepoint before calculating gene-level temporal
metrics, matching the pointwise aggregation logic used elsewhere in the H3
protocol.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
DEFAULT_CANDIDATES = [ROOT / "Data" / "GSE297234_expression.tsv.gz", ROOT / "Data" / "GSE297234_expression.csv.gz"]

THRESHOLDS = {
    "min_progeny_overlap": 0.90,
    "min_dorothea_overlap": 0.90,
    "min_pc1_abs_spearman": 0.80,
    "min_pc1_adjacent_monotonic_fraction": 0.80,
    "min_directional_gene_fraction_abs_spearman_ge_0_80": 0.60,
    "min_median_abs_gene_spearman": 0.60,
    "min_median_abs_day_final_day_zero_log2fc": 0.50,
    "min_time_span_days": 7.0,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--expression", type=Path, default=None)
    p.add_argument("--sample-days", required=True, help="Comma-separated SAMPLE:DAY pairs; multiple samples may share a day")
    p.add_argument("--out", type=Path, default=OUT / "H3_YAMANAKA_STRUCTURAL_CALIBRATION.json")
    return p.parse_args()


def locate_expression(path: Path | None) -> Path:
    if path is not None:
        if not path.exists():
            raise FileNotFoundError(path)
        return path
    for candidate in DEFAULT_CANDIDATES:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("No GSE297234 expression matrix found. Supply --expression explicitly.")


def read_matrix(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    compression = "gzip" if name.endswith(".gz") else None
    d = pd.read_csv(path, sep="," if name.endswith(".csv") or name.endswith(".csv.gz") else "\t", compression=compression, low_memory=False)
    if d.shape[1] < 5:
        raise ValueError(f"Expression matrix has too few columns: {d.shape}")
    x = d.set_index(d.columns[0]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    x.index = x.index.astype(str).str.strip().str.upper()
    return x.groupby(level=0).sum()


def parse_sample_days(spec: str, columns: pd.Index) -> tuple[list[str], np.ndarray]:
    pairs = []
    for item in spec.split(","):
        item = item.strip()
        if not item or ":" not in item:
            raise ValueError(f"Invalid --sample-days item: {item!r}; expected SAMPLE:DAY")
        sample, day = item.rsplit(":", 1)
        sample = sample.strip()
        if sample not in columns:
            raise ValueError(f"Sample {sample!r} not found in expression matrix")
        pairs.append((sample, float(day)))
    if len(pairs) < 4:
        raise ValueError("At least four ordered sample/time assignments are required")
    days = np.asarray([p[1] for p in pairs], dtype=float)
    if len(np.unique(days)) < 4:
        raise ValueError("At least four distinct timepoints are required")
    if not np.all(np.sort(np.unique(days))[1:] > np.sort(np.unique(days))[:-1]):
        raise ValueError("Timepoints must be strictly ordered")
    return [p[0] for p in pairs], days


def aggregate_by_day(x: pd.DataFrame, samples: list[str], days: np.ndarray) -> tuple[pd.DataFrame, np.ndarray]:
    tmp = x[samples].T.copy()
    tmp["__day__"] = days
    y = tmp.groupby("__day__").mean().T
    ordered_days = y.columns.to_numpy(float)
    return y, ordered_days


def log_cpm(x: pd.DataFrame) -> pd.DataFrame:
    lib = x.sum(axis=0).replace(0, np.nan)
    return np.log2(x.div(lib, axis=1) * 1e6 + 1.0)


def temporal_metrics(x: pd.DataFrame, samples: list[str], days: np.ndarray) -> dict:
    day_mean_counts, ordered_days = aggregate_by_day(x, samples, days)
    lcpm = log_cpm(day_mean_counts)
    rhos = []
    for _, row in lcpm.iterrows():
        y = row.to_numpy(float)
        if np.std(y) > 0:
            rhos.append(float(spearmanr(ordered_days, y).statistic))
    rhos = np.asarray(rhos, dtype=float)
    z = lcpm.to_numpy(float)
    z -= z.mean(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    pc1 = vt[0]
    pc1_rho = float(spearmanr(ordered_days, pc1).statistic)
    diffs = np.diff(pc1)
    mono = float(max(np.mean(diffs >= 0), np.mean(diffs <= 0)))
    transition = lcpm.iloc[:, -1] - lcpm.iloc[:, 0]
    return {
        "n_genes": int(len(rhos)),
        "n_timepoints": int(len(ordered_days)),
        "timepoints": ordered_days.tolist(),
        "pc1_abs_spearman": abs(pc1_rho),
        "pc1_adjacent_monotonic_fraction": mono,
        "directional_gene_fraction_abs_spearman_ge_0_80": float(np.mean(np.abs(rhos) >= 0.80)),
        "median_abs_gene_spearman": float(np.median(np.abs(rhos))),
        "median_abs_day_final_day_zero_log2fc": float(np.median(np.abs(transition.to_numpy(float)))),
        "time_span_days": float(ordered_days[-1] - ordered_days[0]),
    }


def representation_overlap(x: pd.DataFrame) -> dict:
    import decoupler as dc

    out = {}
    for family, net in (("PROGENy", dc.op.progeny(organism="human", top=100)), ("DoRothEA", dc.op.dorothea(organism="human", levels=["A", "B", "C"]))):
        targets = pd.Index(net["target"].astype(str).str.upper().unique())
        common = len(x.index.intersection(targets))
        out[family] = {"network_targets": int(len(targets)), "common_genes": int(common), "overlap_fraction": float(common / max(len(targets), 1))}
    return out


def compare_to_thresholds(metrics: dict, overlap: dict) -> dict:
    checks = {
        "PROGENy_overlap": overlap["PROGENy"]["overlap_fraction"] >= THRESHOLDS["min_progeny_overlap"],
        "DoRothEA_overlap": overlap["DoRothEA"]["overlap_fraction"] >= THRESHOLDS["min_dorothea_overlap"],
        "PC1_abs_spearman": metrics["pc1_abs_spearman"] >= THRESHOLDS["min_pc1_abs_spearman"],
        "PC1_monotonicity": metrics["pc1_adjacent_monotonic_fraction"] >= THRESHOLDS["min_pc1_adjacent_monotonic_fraction"],
        "directional_gene_fraction": metrics["directional_gene_fraction_abs_spearman_ge_0_80"] >= THRESHOLDS["min_directional_gene_fraction_abs_spearman_ge_0_80"],
        "median_abs_gene_spearman": metrics["median_abs_gene_spearman"] >= THRESHOLDS["min_median_abs_gene_spearman"],
        "median_abs_endpoint_log2fc": metrics["median_abs_day_final_day_zero_log2fc"] >= THRESHOLDS["min_median_abs_day_final_day_zero_log2fc"],
        "time_span": metrics["time_span_days"] >= THRESHOLDS["min_time_span_days"],
    }
    return {"checks": checks, "all_pass": bool(all(checks.values()))}


def main() -> None:
    args = parse_args()
    path = locate_expression(args.expression)
    x = read_matrix(path)
    samples, days = parse_sample_days(args.sample_days, x.columns)
    metrics = temporal_metrics(x, samples, days)
    overlap = representation_overlap(x)
    comparison = compare_to_thresholds(metrics, overlap)
    result = {
        "status": "ok",
        "purpose": "Yamanaka structural calibration only; no H3 control similarity is run",
        "source_expression": str(path),
        "samples": samples,
        "sample_days": {s: float(d) for s, d in zip(samples, days)},
        "metrics": metrics,
        "representation_overlap": overlap,
        "locked_control_thresholds": THRESHOLDS,
        "threshold_comparison": comparison,
        "interpretation_guardrail": "A threshold failure does not prove that the H3 gate is wrong. It flags a potential target-control commensurability problem that must be reviewed before changing prospective-control thresholds. Thresholds are never changed automatically.",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
