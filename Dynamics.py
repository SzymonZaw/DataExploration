"""Dynamics v2.1: time-aware exploratory dynamics of OSKM reprogramming.

This script reuses the existing EDA PCA outputs. It uses explicit biological time,
preserves replicate observations, and avoids treating PCA axes from different
studies as directly comparable. Derivatives are calculated on the mean state at
each unique time point, so biological replicates do not create duplicate time
coordinates or artificial infinite/NaN velocities.

Scientific status
-----------------
Exploratory dynamical analysis only. It is not a causal model and does not claim
critical transitions, bifurcations, universality, quantum mechanisms, or
relativistic mechanisms. The symbolic-regression table is a scaffold for a later
model-discovery stage; no biological law is fitted here.
"""

from pathlib import Path
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "Dynamics"
OUT.mkdir(parents=True, exist_ok=True)

# Existing PCA outputs from the EDA scripts.
# GSE297234 uses 08_PCA_coordinates.csv: 07 is its feature-variance table.
DATASETS = {
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "07_PCA_coordinates.csv",
    "GSE148158": RESULTS / "GSE148158" / "07_PCA_coordinates.csv",
    "GSE52052": RESULTS / "GSE52052" / "07_PCA_coordinates.csv",
    "GSE67462": RESULTS / "GSE67462" / "07_PCA_coordinates.csv",
    "GSE297234": RESULTS / "GSE297234" / "08_PCA_coordinates.csv",
}

# Biological time is represented in hours. iPSC/hESC are kept as endpoint/control
# observations but are not assigned an artificial elapsed time.
TIME_PATTERNS = {
    "GSE28688": [
        (r"24\s*h", 24.0),
        (r"48\s*h", 48.0),
        (r"72\s*h", 72.0),
    ],
    "GSE148158": [
        (r"48\s*h|\b48\b", 48.0),
        (r"72\s*h|\b72\b", 72.0),
    ],
    "GSE52052": [(r"day\s*11", 264.0)],
    "GSE67462": [
        (r"day\s*0\b", 0.0),
        (r"day\s*1\b", 24.0),
        (r"day\s*3\b", 72.0),
        (r"day\s*5\b", 120.0),
        (r"day\s*7\b", 168.0),
        (r"day\s*11\b", 264.0),
        (r"day\s*15\b", 360.0),
        (r"day\s*18\b", 432.0),
    ],
    "GSE297234": [
        (r"d0\b|day\s*0\b", 0.0),
        (r"d3\b|day\s*3\b", 72.0),
        (r"d7\b|day\s*7\b", 168.0),
        (r"d10\b|day\s*10\b", 240.0),
    ],
}


def load_pca(path):
    """Load PCA coordinates and return available numeric PC1-PC3 columns."""
    if not path.exists():
        return None
    df = pd.read_csv(path, index_col=0)
    pcs = [c for c in ("PC1", "PC2", "PC3") if c in df.columns]
    if not pcs:
        return None
    coords = df[pcs].apply(pd.to_numeric, errors="coerce")
    coords.index = coords.index.astype(str)
    return coords


def infer_time_hours(dataset, sample):
    """Infer numeric experimental time from a sample identifier."""
    text = str(sample).lower().replace("_", " ").replace("-", " ")
    for pattern, hours in TIME_PATTERNS.get(dataset, []):
        if re.search(pattern, text):
            return hours
    return np.nan


def infer_stage(sample, time_hours):
    """Return a readable stage label without inventing time for endpoints."""
    text = str(sample).lower().replace("_", " ").replace("-", " ")
    if re.search(r"ipsc|\bips\b", text):
        return "iPSC"
    if re.search(r"hesc", text):
        return "hESC"
    if pd.notna(time_hours):
        days = time_hours / 24.0
        if days.is_integer():
            return f"day{int(days)}"
        return f"{time_hours:g}h"
    return "unknown"


def parse_replicate(sample):
    """Extract an explicit replicate label when present in the sample name."""
    text = str(sample)
    match = re.search(r"(?:-|_|\s)([ab])$", text, flags=re.IGNORECASE)
    if match:
        return match.group(1).lower()
    match = re.search(r"(?:rep|replicate)[_\s-]*(\d+)", text, flags=re.IGNORECASE)
    if match:
        return match.group(1)
    return "unknown"


def orient_pc(series):
    """Choose a deterministic sign for presentation; sign has no biology."""
    x = pd.Series(series, dtype=float).copy()
    valid = x.dropna()
    if valid.empty:
        return x
    idx = valid.abs().idxmax()
    return -x if x.loc[idx] < 0 else x


def zscore(series):
    """Within-study z-score to remove arbitrary PCA scale."""
    x = pd.Series(series, dtype=float)
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=x.index)
    return (x - x.mean()) / sd


def finite_difference(values, times):
    """Derivative on strictly increasing, unique time points."""
    x = np.asarray(values, dtype=float)
    t = np.asarray(times, dtype=float)
    result = np.full(len(x), np.nan, dtype=float)
    valid = np.isfinite(x) & np.isfinite(t)
    if valid.sum() < 2:
        return result

    xv = x[valid]
    tv = t[valid]
    order = np.argsort(tv)
    xv = xv[order]
    tv = tv[order]

    # This function expects unique times. Duplicates are collapsed before call.
    if np.any(np.diff(tv) <= 0):
        raise ValueError("finite_difference requires strictly increasing unique times")

    if len(xv) == 2:
        velocity = (xv[1] - xv[0]) / (tv[1] - tv[0])
        values_out = np.full(2, velocity)
    else:
        values_out = np.gradient(xv, tv)

    valid_indices = np.where(valid)[0][order]
    result[valid_indices] = values_out
    return result


def rolling_indicators(values, window=5):
    """Exploratory rolling variance and lag-1 autocorrelation."""
    s = pd.Series(values, dtype=float)
    if len(s) < 3:
        return pd.Series(np.nan, index=s.index), pd.Series(np.nan, index=s.index)
    effective_window = min(window, len(s))
    min_periods = 3
    variance = s.rolling(effective_window, min_periods=min_periods).var()
    autocorr = s.rolling(effective_window, min_periods=min_periods).apply(
        lambda x: x.autocorr(lag=1) if x.std() > 0 else np.nan,
        raw=False,
    )
    return variance, autocorr


def candidate_library(x):
    """Interpretable basis for future symbolic regression / GP."""
    x = np.asarray(x, dtype=float)
    return pd.DataFrame({
        "1": np.ones_like(x),
        "x": x,
        "x2": x ** 2,
        "x3": x ** 3,
        "sin_x": np.sin(x),
        "cos_x": np.cos(x),
        "log_abs_x": np.log1p(np.abs(x)),
    })


def process_dataset(dataset, path):
    coords = load_pca(path)
    if coords is None:
        return None, None

    out = coords.copy()
    out.insert(0, "sample", out.index.astype(str))
    out["dataset"] = dataset
    out["time_hours"] = [infer_time_hours(dataset, s) for s in out["sample"]]
    out["stage"] = [infer_stage(s, t) for s, t in zip(out["sample"], out["time_hours"])]
    out["replicate"] = [parse_replicate(s) for s in out["sample"]]

    pcs = [c for c in ("PC1", "PC2", "PC3") if c in coords.columns]
    for pc in pcs:
        out[f"{pc}_oriented"] = orient_pc(out[pc])
        out[f"{pc}_z"] = zscore(out[f"{pc}_oriented"])

    timed = out[out["time_hours"].notna()].copy()
    timed = timed.sort_values("time_hours")

    # Critical fix: replicate samples at the same biological time must not be
    # passed directly to np.gradient. We first average the standardized state
    # across biological replicates at each unique time point.
    if len(timed):
        time_means = timed.groupby("time_hours", as_index=False)[
            [f"{pc}_z" for pc in pcs]
        ].mean()
        time_means = time_means.sort_values("time_hours").reset_index(drop=True)

        for pc in pcs:
            derivative = finite_difference(
                time_means[f"{pc}_z"].to_numpy(),
                time_means["time_hours"].to_numpy(),
            )
            time_means[f"d{pc}_dt"] = derivative

        velocity_cols = [f"d{pc}_dt" for pc in pcs]
        time_means["state_speed"] = np.sqrt(
            np.sum(time_means[velocity_cols].to_numpy() ** 2, axis=1)
        )

        variance, autocorr = rolling_indicators(
            time_means["PC1_z"].to_numpy(), window=5
        )
        time_means["rolling_variance_PC1"] = variance.to_numpy()
        time_means["rolling_autocorrelation_PC1"] = autocorr.to_numpy()

        derivative_cols = ["time_hours"] + [f"d{pc}_dt" for pc in pcs]
        derivative_cols += [
            "state_speed",
            "rolling_variance_PC1",
            "rolling_autocorrelation_PC1",
        ]
        out = out.merge(time_means[derivative_cols], on="time_hours", how="left")
    else:
        for pc in pcs:
            out[f"d{pc}_dt"] = np.nan
        out["state_speed"] = np.nan
        out["rolling_variance_PC1"] = np.nan
        out["rolling_autocorrelation_PC1"] = np.nan

    return out, pcs


def build_trajectory_summary(states):
    """Summarize trajectories using study-internal, time-aware quantities."""
    rows = []
    for dataset, group in states.groupby("dataset"):
        timed = group[group["time_hours"].notna()].copy()
        unique = (
            timed.groupby("time_hours")[[c for c in ["PC1_z", "PC2_z", "PC3_z"] if c in timed]]
            .mean()
            .sort_index()
        )
        row = {
            "dataset": dataset,
            "n_samples": len(group),
            "n_timed_samples": len(timed),
            "n_unique_times": len(unique),
            "time_min_hours": unique.index.min() if len(unique) else np.nan,
            "time_max_hours": unique.index.max() if len(unique) else np.nan,
            "n_replicate_labelled": int((group["replicate"] != "unknown").sum()),
        }
        if len(unique) >= 2:
            row["net_PC1_z_change"] = unique["PC1_z"].iloc[-1] - unique["PC1_z"].iloc[0]
            row["net_PC2_z_change"] = (
                unique["PC2_z"].iloc[-1] - unique["PC2_z"].iloc[0]
                if "PC2_z" in unique else np.nan
            )
            row["net_PC3_z_change"] = (
                unique["PC3_z"].iloc[-1] - unique["PC3_z"].iloc[0]
                if "PC3_z" in unique else np.nan
            )
        else:
            row["net_PC1_z_change"] = np.nan
            row["net_PC2_z_change"] = np.nan
            row["net_PC3_z_change"] = np.nan

        speeds = group["state_speed"].dropna()
        row["mean_state_speed"] = speeds.mean() if len(speeds) else np.nan
        row["max_state_speed"] = speeds.max() if len(speeds) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def build_symbolic_training_table(states):
    """Return derivative observations prepared for later equation discovery."""
    cols = [
        "dataset", "sample", "stage", "replicate", "time_hours",
        "PC1_z", "dPC1_dt", "PC2_z", "dPC2_dt", "PC3_z", "dPC3_dt",
        "state_speed",
    ]
    table = states[[c for c in cols if c in states.columns]].copy()
    table = table[table["time_hours"].notna()].copy()
    return table.sort_values(["dataset", "time_hours", "sample"])


def main():
    availability = []
    processed = []

    for dataset, path in DATASETS.items():
        state, pcs = process_dataset(dataset, path)
        found = state is not None
        availability.append({
            "dataset": dataset,
            "PCA_file_found": found,
            "path": str(path),
            "PCs_available": ",".join(pcs) if pcs else "",
        })
        if found:
            processed.append(state)

    availability_df = pd.DataFrame(availability)
    availability_df.to_csv(OUT / "01_dataset_availability.csv", index=False)

    if not processed:
        print("No existing PCA trajectories found. Run the EDA scripts first.")
        return

    states = pd.concat(processed, ignore_index=True)
    states.to_csv(OUT / "02_time_aware_state_dynamics.csv", index=False)

    summary = build_trajectory_summary(states)
    summary.to_csv(OUT / "03_trajectory_summary.csv", index=False)

    training = build_symbolic_training_table(states)
    training.to_csv(OUT / "04_symbolic_training_table.csv", index=False)

    # Candidate library uses unique time-point means to avoid counting the same
    # biological time multiple times merely because it has replicates.
    if len(training):
        source = (
            training.groupby(["dataset", "time_hours"], as_index=False)
            .agg({"PC1_z": "mean", "dPC1_dt": "mean"})
            .dropna(subset=["PC1_z", "dPC1_dt"])
        )
        if len(source):
            library = candidate_library(source["PC1_z"].to_numpy())
            library.insert(0, "dataset", source["dataset"].to_numpy())
            library.insert(1, "time_hours", source["time_hours"].to_numpy())
            library["target_dPC1_dt"] = source["dPC1_dt"].to_numpy()
            library.to_csv(OUT / "05_symbolic_candidate_library.csv", index=False)

    timed_counts = (
        states[states["time_hours"].notna()]
        .groupby("dataset")["time_hours"]
        .nunique()
        .to_dict()
    )

    report = [
        "Dynamics v2.1 — time-aware exploratory dynamics of OSKM reprogramming",
        "",
        "Implementation status",
        "--------------------",
        "All derivatives are calculated against explicit biological time in hours.",
        "Replicates at the same time point are first averaged for derivative estimation.",
        "The individual replicate rows remain in the state table for uncertainty analysis.",
        "",
        "Datasets",
        "--------",
        f"PCA datasets available: {int(availability_df['PCA_file_found'].sum())}/{len(availability_df)}.",
    ]
    for dataset in DATASETS:
        report.append(f"{dataset}: {timed_counts.get(dataset, 0)} unique numeric time point(s).")

    report += [
        "",
        "Scientific interpretation",
        "--------------------------",
        "The primary object is a time-dependent low-dimensional state trajectory within",
        "each experiment. PC1-PC3 are standardized within study; they are not assumed",
        "to form a universal biological coordinate system across platforms or species.",
        "",
        "Important methodological safeguards",
        "------------------------------------",
        "- Duplicate biological times are aggregated before numerical differentiation.",
        "- No artificial time is assigned to iPSC or hESC endpoint/control samples.",
        "- Derivatives are descriptive and depend on sampling density and preprocessing.",
        "- Rolling variance and autocorrelation are exploratory stability indicators only.",
        "- No causal, universal, critical-transition, quantum, or relativistic claim is made.",
        "",
        "Research target",
        "---------------",
        "Test whether independent OSKM/reprogramming experiments contain a robust",
        "low-dimensional temporal structure that can generalise across experiments.",
        "A strong future test is to learn a dynamical model on one complete study and",
        "validate it on an entirely held-out study.",
        "",
        "Next layer",
        "-----------",
        "1. Define a genuinely shared biological state using harmonised genes/markers",
        "   rather than assuming PCA coordinates are shared.",
        "2. Quantify replicate uncertainty with bootstrap or hierarchical models.",
        "3. Compare multiple state definitions and preprocessing choices.",
        "4. Fit symbolic-regression / Genetic-Programming equations and simple baselines.",
        "5. Hold out an entire study, preferably GSE297234, for external validation.",
        "6. Later connect expression dynamics to GSE67520 regulatory dynamics.",
    ]

    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")

    print(f"Dynamics v2.1 results written to: {OUT}")
    print(f"Datasets with PCA: {int(availability_df['PCA_file_found'].sum())}/{len(availability_df)}")
    print(f"Time-aware observations: {len(training)}")
    print("Replicate-safe derivative calculation: enabled")


if __name__ == "__main__":
    main()
