from __future__ import annotations

"""Audit whether the frozen H3 structural gates are commensurate with Yamanaka.

This is a target-calibration audit, not an H3 similarity test and not a
threshold-optimization procedure. It evaluates the two Yamanaka donor
trajectories separately and the canonical trajectory defined as the pointwise
arithmetic mean of donor-wise normalized trajectories. It also reports the
alternative pooled-count construction used by the earlier calibration script,
so that the aggregation choice is explicit rather than implicit.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"

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


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--expression", type=Path, required=True)
    p.add_argument("--sample-days", required=True)
    p.add_argument("--donors", required=True, help="DONOR:SAMPLE1,SAMPLE2,...;DONOR2:SAMPLE...")
    p.add_argument("--out", type=Path, default=OUT / "H3_TARGET_COMMENSURABILITY_AUDIT_v1.json")
    return p.parse_args()


def read_matrix(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    compression = "gzip" if name.endswith(".gz") else None
    sep = "," if name.endswith(".csv") or name.endswith(".csv.gz") else "\t"
    d = pd.read_csv(path, sep=sep, compression=compression, low_memory=False)
    x = d.set_index(d.columns[0]).apply(pd.to_numeric, errors="coerce").fillna(0.0)
    x.index = x.index.astype(str).str.strip().str.upper()
    x = x.groupby(level=0).sum()
    if x.shape[1] == 0:
        raise ValueError(
            "Expression matrix has no sample columns after using the first column as gene ID. "
            f"Input: {path}"
        )
    return x


def parse_days(spec: str, columns: pd.Index) -> dict[str, float]:
    out: dict[str, float] = {}
    for item in spec.split(","):
        sample, day = item.strip().rsplit(":", 1)
        if sample not in columns:
            raise ValueError(f"Sample not found: {sample}")
        out[sample] = float(day)
    return out


def parse_donors(spec: str, days: dict[str, float]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for group in spec.split(";"):
        donor, raw_samples = group.strip().split(":", 1)
        samples = [s.strip() for s in raw_samples.split(",") if s.strip()]
        if len(samples) < 4:
            raise ValueError(f"Donor {donor} has fewer than 4 samples")
        if len(set(days[s] for s in samples)) != len(samples):
            raise ValueError(f"Donor {donor} must have one sample per timepoint")
        out[donor] = samples
    if len(out) < 2:
        raise ValueError("At least two independent donor trajectories are required")
    return out


def library_sizes(x: pd.DataFrame) -> pd.Series:
    lib = x.sum(axis=0, dtype=float)
    if not np.isfinite(lib.to_numpy()).all():
        bad = lib[~np.isfinite(lib)]
        raise ValueError(f"Non-finite library sizes detected: {bad.to_dict()}")
    if (lib <= 0).any():
        bad = lib[lib <= 0]
        raise ValueError(
            "Non-positive library size detected. This means the selected sample columns "
            "contain no positive counts and cannot be normalized to CPM: "
            f"{bad.to_dict()}"
        )
    return lib


def log_cpm(x: pd.DataFrame) -> pd.DataFrame:
    lib = library_sizes(x)
    out = np.log2(x.div(lib, axis=1) * 1e6 + 1.0)
    if not np.isfinite(out.to_numpy(dtype=float, copy=False)).all():
        bad_cols = out.columns[~np.isfinite(out.to_numpy(dtype=float, copy=False)).all(axis=0)].tolist()
        raise ValueError(f"Non-finite log2(CPM+1) values detected in columns: {bad_cols}")
    return out


def sample_log_cpm(x: pd.DataFrame) -> pd.DataFrame:
    return log_cpm(x)


def pointwise_donor_trajectory(x: pd.DataFrame, samples: list[str], days: dict[str, float]) -> pd.DataFrame:
    ordered = sorted(samples, key=lambda s: days[s])
    return sample_log_cpm(x[ordered])


def _pc1_from_time_covariance(matrix: np.ndarray) -> np.ndarray:
    """Stable PC1 for a gene x time matrix with only a few timepoints.

    This is computationally equivalent to SVD for the temporal score direction,
    but works on the small time-by-time Gram matrix and avoids LAPACK SVD
    convergence failures on ill-conditioned expression matrices.
    """
    z = matrix - matrix.mean(axis=1, keepdims=True)
    gram = z.T @ z
    gram = (gram + gram.T) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    pc1 = eigenvectors[:, int(np.argmax(eigenvalues))]
    if not np.isfinite(pc1).all():
        raise ValueError("PC1 contains non-finite values")
    return pc1


def metrics(lcpm: pd.DataFrame, days: list[float], label: str = "trajectory") -> dict:
    t = np.asarray(days, dtype=float)
    matrix = lcpm.to_numpy(dtype=float, copy=True)
    finite_rows = np.isfinite(matrix).all(axis=1)
    n_total = int(matrix.shape[0])
    n_finite = int(finite_rows.sum())
    if n_finite == 0:
        finite_by_column = np.isfinite(matrix).sum(axis=0).tolist()
        raise ValueError(
            f"No finite gene trajectories available for {label}. "
            f"shape={matrix.shape}; finite_rows={n_finite}/{n_total}; "
            f"finite_values_by_timepoint={finite_by_column}. "
            "This indicates a normalization/input problem upstream of PC1, not an SVD/PC1 problem."
        )
    matrix = matrix[finite_rows]

    rhos = []
    for y in matrix:
        if np.std(y) > 0:
            rho = spearmanr(t, y).statistic
            if np.isfinite(rho):
                rhos.append(float(rho))
    rhos = np.asarray(rhos, dtype=float)
    if len(rhos) == 0:
        raise ValueError(f"No variable finite gene trajectories available for {label}")

    z = matrix - matrix.mean(axis=1, keepdims=True)
    if not np.isfinite(z).all():
        raise ValueError(f"Non-finite values remain after centering for {label}")
    pc1 = _pc1_from_time_covariance(matrix)
    pc1_rho = spearmanr(t, pc1).statistic
    if not np.isfinite(pc1_rho):
        raise ValueError(f"PC1/time Spearman correlation is non-finite for {label}")
    diffs = np.diff(pc1)
    mono = float(max(np.mean(diffs >= 0), np.mean(diffs <= 0)))
    endpoint = matrix[:, -1] - matrix[:, 0]
    abs_rho = np.abs(rhos)
    return {
        "n_genes_with_variable_trajectory": int(len(rhos)),
        "n_genes_used_for_pc1": int(matrix.shape[0]),
        "n_timepoints": int(len(t)),
        "timepoints_days": t.tolist(),
        "time_span_days": float(t[-1] - t[0]),
        "pc1_abs_spearman": float(abs(pc1_rho)),
        "pc1_adjacent_monotonic_fraction": mono,
        "directional_gene_fraction_abs_spearman_ge_0_80": float(np.mean(abs_rho >= 0.80)),
        "median_abs_gene_spearman": float(np.median(abs_rho)),
        "median_abs_endpoint_log2fc": float(np.median(np.abs(endpoint))),
        "gene_abs_spearman_quantiles": {str(q): float(np.quantile(abs_rho, q)) for q in [0.1, 0.25, 0.5, 0.75, 0.9]},
        "endpoint_abs_log2fc_quantiles": {str(q): float(np.quantile(np.abs(endpoint), q)) for q in [0.1, 0.25, 0.5, 0.75, 0.9]},
    }


def compare(metrics_dict: dict) -> dict:
    checks = {
        "PC1_abs_spearman": metrics_dict["pc1_abs_spearman"] >= THRESHOLDS["min_pc1_abs_spearman"],
        "PC1_monotonicity": metrics_dict["pc1_adjacent_monotonic_fraction"] >= THRESHOLDS["min_pc1_adjacent_monotonic_fraction"],
        "directional_gene_fraction": metrics_dict["directional_gene_fraction_abs_spearman_ge_0_80"] >= THRESHOLDS["min_directional_gene_fraction_abs_spearman_ge_0_80"],
        "median_abs_gene_spearman": metrics_dict["median_abs_gene_spearman"] >= THRESHOLDS["min_median_abs_gene_spearman"],
        "median_abs_endpoint_log2fc": metrics_dict["median_abs_endpoint_log2fc"] >= THRESHOLDS["min_median_abs_day_final_day_zero_log2fc"],
        "time_span": metrics_dict["time_span_days"] >= THRESHOLDS["min_time_span_days"],
    }
    return {"checks": checks, "all_temporal_gates_pass": bool(all(checks.values()))}


def representation_overlap(x: pd.DataFrame) -> dict:
    import decoupler as dc
    out = {}
    for family, net in (("PROGENy", dc.op.progeny(organism="human", top=100)), ("DoRothEA", dc.op.dorothea(organism="human", levels=["A", "B", "C"]))):
        targets = pd.Index(net["target"].astype(str).str.upper().unique())
        common = len(x.index.intersection(targets))
        out[family] = {"network_targets": int(len(targets)), "common_genes": int(common), "overlap_fraction": float(common / max(len(targets), 1))}
    return out


def main() -> None:
    a = args()
    x = read_matrix(a.expression)
    days = parse_days(a.sample_days, x.columns)
    donors = parse_donors(a.donors, days)

    selected_samples = [s for samples in donors.values() for s in samples]
    selected_lib = library_sizes(x[selected_samples])
    print(json.dumps({
        "input_shape": [int(x.shape[0]), int(x.shape[1])],
        "selected_samples": selected_samples,
        "selected_library_sizes": {k: float(v) for k, v in selected_lib.items()},
        "selected_nonzero_gene_counts": {k: int((x[k] > 0).sum()) for k in selected_samples},
    }, indent=2))

    donor_traj = {}
    donor_metrics = {}
    for donor, samples in donors.items():
        traj = pointwise_donor_trajectory(x, samples, days)
        ordered_days = [days[s] for s in sorted(samples, key=lambda s: days[s])]
        donor_traj[donor] = traj
        donor_metrics[donor] = metrics(traj, ordered_days, label=f"donor {donor}")

    common_index = donor_traj[next(iter(donor_traj))].index
    if any(not common_index.equals(v.index) for v in donor_traj.values()):
        raise RuntimeError("Donor gene indices are not identical after deterministic mapping")
    canonical = sum(donor_traj.values()) / len(donor_traj)
    canonical_days = donor_metrics[next(iter(donor_metrics))]["timepoints_days"]
    canonical_metrics = metrics(canonical, canonical_days, label="canonical Yamanaka")

    pooled_parts = []
    for day in canonical_days:
        samples = [s for donor_samples in donors.values() for s in donor_samples if days[s] == day]
        pooled_parts.append(x[samples].mean(axis=1))
    pooled_counts = pd.concat(pooled_parts, axis=1)
    pooled_counts.columns = canonical_days
    pooled_lcpm = log_cpm(pooled_counts)
    pooled_metrics = metrics(pooled_lcpm, canonical_days, label="pooled-count alternative")

    overlap = representation_overlap(x)
    result = {
        "status": "ok",
        "purpose": "H3 target commensurability audit v1; no control similarity or threshold optimization",
        "source_expression": str(a.expression),
        "sample_days": days,
        "donor_assignments": donors,
        "normalization_definition": "Each sample is normalized to log2(CPM+1); donor trajectories are then averaged pointwise across donors to form the canonical target.",
        "historical_representation_note": "The historical exploratory PCA used log1p(CPM). This audit uses log2(CPM+1) to match the prospective H3 metric space; the rank-based temporal metrics are invariant to the log base.",
        "donor_metrics": donor_metrics,
        "canonical_yamanaka_metrics": canonical_metrics,
        "canonical_threshold_comparison": compare(canonical_metrics),
        "pooled_count_alternative_metrics": pooled_metrics,
        "canonical_vs_pooled_absolute_metric_difference": {
            key: abs(canonical_metrics[key] - pooled_metrics[key])
            for key in ["pc1_abs_spearman", "pc1_adjacent_monotonic_fraction", "directional_gene_fraction_abs_spearman_ge_0_80", "median_abs_gene_spearman", "median_abs_endpoint_log2fc"]
        },
        "representation_overlap": overlap,
        "locked_control_thresholds": THRESHOLDS,
        "interpretation_guardrails": [
            "This audit does not change any threshold.",
            "A Yamanaka gate failure is a commensurability warning, not evidence that the Yamanaka biology is absent.",
            "The canonical target is donor-wise normalized pointwise averaging, matching the H3 protocol; pooled-count averaging is reported only as a sensitivity diagnostic.",
            "No prospective control should be selected, rejected, or similarity-tested from these metrics alone without the protocol review decision."
        ],
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
