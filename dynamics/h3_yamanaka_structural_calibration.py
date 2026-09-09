from __future__ import annotations

"""Calibrate H3 structural eligibility metrics on the Yamanaka source data.

This is a calibration/audit tool, not an H3 similarity test. It deliberately
computes the same structural quantities used for prospective controls, so that
we can determine whether those quantities are commensurate with the biological
process being used as the target. It never changes thresholds automatically.

Usage examples:
  python dynamics/h3_yamanaka_structural_calibration.py
  python dynamics/h3_yamanaka_structural_calibration.py --expression Data/GSE297234_expression.tsv.gz --days 0,3,7,10

The expression matrix must be gene x sample, with the first column containing
gene identifiers. Sample columns may be explicitly supplied with --samples.
Otherwise the script attempts to infer four-point day labels from sample names.
For a strict calibration, use the same gene-level input representation that was
used to construct the prospective H3 controls.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
DEFAULT_CANDIDATES = [
    ROOT / "Data" / "GSE297234_expression.tsv.gz",
    ROOT / "Data" / "GSE297234_expression.csv.gz",
    ROOT / "Data" / "GSE297234_raw_counts_matrix.csv.gz",
]

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
    p.add_argument("--days", default="0,3,7,10")
    p.add_argument("--samples", default=None, help="Comma-separated sample columns in chronological order")
    p.add_argument("--gene-column", default=None)
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
    raise FileNotFoundError(
        "No Yamanaka expression matrix found. Supply --expression pointing to the "
        "gene x sample matrix used for the Yamanaka representation."
    )


def read_matrix(path: Path) -> pd.DataFrame:
    suffix = path.name.lower()
    compression = "gzip" if suffix.endswith(".gz") else None
    if suffix.endswith(".csv") or suffix.endswith(".csv.gz"):
        d = pd.read_csv(path, compression=compression, low_memory=False)
    else:
        d = pd.read_csv(path, sep="\t", compression=compression, low_memory=False)
    if d.shape[1] < 5:
        raise ValueError(f"Expression matrix has too few columns: {d.shape}")
    x = d.set_index(d.columns[0]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    x.index = x.index.astype(str).str.strip().str.upper()
    return x.groupby(level=0).sum()


def infer_samples(columns: list[str], days: np.ndarray) -> list[str]:
    # Accept common day encodings: D0, day0, _0, etc. Only retain columns with
    # an unambiguous numeric day that is one of the requested points.
    out = []
    for day in days:
        hits = []
        for c in columns:
            s = str(c)
            m = re.search(r"(?:DAY|D)[_ -]?(\d+(?:\.\d+)?)\b", s, re.I)
            if m is None:
                continue
            if abs(float(m.group(1)) - float(day)) < 1e-9:
                hits.append(c)
        if len(hits) != 1:
            raise ValueError(
                f"Could not infer exactly one sample for day {day}. Found {hits}. "
                "Supply --samples explicitly."
            )
        out.append(hits[0])
    return out


def log_cpm(x: pd.DataFrame) -> pd.DataFrame:
    lib = x.sum(axis=0).replace(0, np.nan)
    return np.log2(x.div(lib, axis=1) * 1e6 + 1.0)


def temporal_metrics(x: pd.DataFrame, samples: list[str], days: np.ndarray) -> dict:
    lcpm = log_cpm(x[samples])
    rhos = []
    for _, row in lcpm.iterrows():
        y = row.to_numpy(float)
        if np.std(y) > 0:
            rhos.append(float(spearmanr(days, y).statistic))
    rhos = np.asarray(rhos, dtype=float)
    z = lcpm.to_numpy(float)
    z -= z.mean(axis=1, keepdims=True)
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    pc1 = vt[0]
    pc1_rho = float(spearmanr(days, pc1).statistic)
    diffs = np.diff(pc1)
    mono = float(max(np.mean(diffs >= 0), np.mean(diffs <= 0)))
    transition = lcpm.iloc[:, -1] - lcpm.iloc[:, 0]
    return {
        "n_genes": int(len(rhos)),
        "pc1_abs_spearman": abs(pc1_rho),
        "pc1_adjacent_monotonic_fraction": mono,
        "directional_gene_fraction_abs_spearman_ge_0_80": float(np.mean(np.abs(rhos) >= 0.80)),
        "median_abs_gene_spearman": float(np.median(np.abs(rhos))),
        "median_abs_day_final_day_zero_log2fc": float(np.median(np.abs(transition.to_numpy(float)))),
        "time_span_days": float(days[-1] - days[0]),
    }


def representation_overlap(x: pd.DataFrame) -> dict:
    import decoupler as dc

    out = {}
    for family, net in (
        ("PROGENy", dc.op.progeny(organism="human", top=100)),
        ("DoRothEA", dc.op.dorothea(organism="human", levels=["A", "B", "C"])),
    ):
        targets = pd.Index(net["target"].astype(str).str.upper().unique())
        common = len(x.index.intersection(targets))
        out[family] = {
            "network_targets": int(len(targets)),
            "common_genes": int(common),
            "overlap_fraction": float(common / max(len(targets), 1)),
        }
    return out


def compare_to_thresholds(metrics: dict, overlap: dict) -> dict:
    checks = {
        "progeny_overlap": overlap["PROGENy"] if False else overlap["PROGENy"]
    }
    # The odd-looking key guard above is avoided below; retain explicit names so
    # a typo in a future network label cannot silently pass the audit.
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
    days = np.asarray([float(v) for v in args.days.split(",")], dtype=float)
    if len(days) < 4 or not np.all(np.diff(days) > 0):
        raise ValueError("--days must contain at least four strictly increasing timepoints")
    path = locate_expression(args.expression)
    x = read_matrix(path)
    samples = [s.strip() for s in args.samples.split(",")] if args.samples else infer_samples(list(x.columns), days)
    if len(samples) != len(days) or any(s not in x.columns for s in samples):
        raise ValueError(f"Invalid --samples: {samples}")
    metrics = temporal_metrics(x, samples, days)
    overlap = representation_overlap(x)
    comparison = compare_to_thresholds(metrics, overlap)
    result = {
        "status": "ok",
        "purpose": "Yamanaka structural calibration only; no H3 control similarity is run",
        "source_expression": str(path),
        "samples": samples,
        "days": days.tolist(),
        "metrics": metrics,
        "representation_overlap": overlap,
        "locked_control_thresholds": THRESHOLDS,
        "threshold_comparison": comparison,
        "interpretation_guardrail": (
            "A threshold failure does not prove that the H3 gate is wrong. It flags a "
            "potential target-control commensurability problem that must be reviewed before "
            "changing any prospective-control threshold. Thresholds are never changed automatically."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
