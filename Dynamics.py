"""Dynamics v2: time-aware exploratory dynamics of OSKM reprogramming.

This script reuses the existing EDA PCA outputs, but unlike the original seed it
uses explicit biological time and preserves replicate structure. It deliberately
does NOT assume that PCA axes are comparable between studies. Instead, each study
is analysed in its own PCA coordinate system and then transformed into a
within-study standardized trajectory representation for cross-study shape
comparison.

Scientific status
-----------------
This is an exploratory dynamical-analysis layer, not a causal model. It is meant
to test whether independently measured reprogramming trajectories show robust,
low-dimensional temporal structure that can later support equation discovery.
Critical-transition metrics are hypothesis-generating only. No claim of
universality, causality, quantum mechanism, or relativistic mechanism is made.
"""

from pathlib import Path
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "Dynamics"
OUT.mkdir(parents=True, exist_ok=True)

# PCA outputs produced by the existing EDA scripts.
# GSE297234 uses 08_PCA_coordinates.csv because 07 is its feature-variance table.
DATASETS = {
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "07_PCA_coordinates.csv",
    "GSE148158": RESULTS / "GSE148158" / "07_PCA_coordinates.csv",
    "GSE52052": RESULTS / "GSE52052" / "07_PCA_coordinates.csv",
    "GSE67462": RESULTS / "GSE67462" / "07_PCA_coordinates.csv",
    "GSE297234": RESULTS / "GSE297234" / "08_PCA_coordinates.csv",
}

# Explicit biological-time metadata where the experimental design is known.
# Times are in hours, with iPSC represented as the final qualitative endpoint
# after the measured time course. It is excluded from finite-difference velocity
# unless a numeric time is available.
TIME_PATTERNS = {
    "GSE28688": [
        (r"24\s*h", 24.0),
        (r"48\s*h", 48.0),
        (r"72\s*h", 72.0),
    ],
    "GSE148158": [
        (r"48", 48.0),
        (r"72", 72.0),
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


# ---------------------------------------------------------------------------
# General helpers
# ---------------------------------------------------------------------------


def load_pca(path):
    """Load PCA coordinates and retain numeric PC columns."""
    if not path.exists():
        return None
    df = pd.read_csv(path, index_col=0)
    cols = [c for c in df.columns if str(c).startswith("PC")]
    if not cols:
        return None
    coords = df[cols].apply(pd.to_numeric, errors="coerce")
    coords.index = coords.index.astype(str)
    return coords


def infer_time_hours(dataset, sample):
    """Infer experimental time from the sample identifier/metadata text."""
    text = str(sample).lower().replace("_", " ").replace("-", " ")
    for pattern, hours in TIME_PATTERNS.get(dataset, []):
        if re.search(pattern, text):
            return hours
    # iPSC is an endpoint rather than a regular elapsed-time measurement.
    if re.search(r"\bips\b|ipsc|hesc|h1 hesc|h9 hesc", text):
        return np.nan
    return np.nan


def infer_stage(dataset, sample, time_hours):
    """Create a readable stage label from sample text/time."""
    text = str(sample).lower()
    if re.search(r"ipsc|\bips\b", text):
        return "iPSC"
    if re.search(r"hesc|h1 hesc|h9 hesc", text):
        return "hESC"
    if pd.notna(time_hours):
        if time_hours % 24 == 0:
            return f"day{int(time_hours // 24)}"
        return f"{int(time_hours)}h"
    return "unknown"


def parse_replicate(dataset, sample):
    """Extract a simple replicate identifier where the naming scheme allows it."""
    text = str(sample)
    # Common replicate suffixes such as -a/-b or _1/_2.
    m = re.search(r"(?:-|_|\s)([ab])$", text, flags=re.IGNORECASE)
    if m:
        return m.group(1).lower()
    m = re.search(r"(?:rep|replicate)[_\s-]*(\d+)", text, flags=re.IGNORECASE)
    if m:
        return m.group(1)
    # GSM IDs themselves are unique but are not biological replicate labels.
    return "unknown"


def standardize_series(x):
    """Within-study z-score; preserves shape but removes arbitrary PCA scale."""
    x = pd.Series(x, dtype=float)
    sd = x.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return x * np.nan
    return (x - x.mean()) / sd


def orient_pc(series):
    """Deterministic orientation only; PC sign has no biological meaning."""
    x = pd.Series(series, dtype=float).copy()
    if x.notna().sum() == 0:
        return x
    i = x.abs().idxmax()
    return -x if x.loc[i] < 0 else x


def finite_difference(values, times):
    """Time-aware local derivative. Returns NaN where time is unavailable."""
    x = np.asarray(values, dtype=float)
    t = np.asarray(times, dtype=float)
    out = np.full(len(x), np.nan, dtype=float)
    valid = np.isfinite(x) & np.isfinite(t)
    if valid.sum() < 2:
        return out
    xv = x[valid]
    tv = t[valid]
    order = np.argsort(tv)
    xv = xv[order]
    tv = tv[order]
    if len(xv) == 2:
        v = np.diff(xv) / np.diff(tv)
        out[np.where(valid)[0][order]] = [v[0], v[0]]
        return out
    grad = np.gradient(xv, tv)
    out[np.where(valid)[0][order]] = grad
    return out


def rolling_indicators(values, window=5):
    """Exploratory rolling variance and lag-1 autocorrelation."""
    s = pd.Series(values, dtype=float)
    min_periods = max(3, min(window, len(s)))
    variance = s.rolling(window, min_periods=min_periods).var()
    autocorr = s.rolling(window, min_periods=min_periods).apply(
        lambda x: pd.Series(x).autocorr(lag=1) if pd.Series(x).std() else np.nan,
        raw=False,
    )
    return variance, autocorr


def candidate_library(x):
    """Small interpretable basis for a later symbolic-regression/GP stage."""
    x = np.asarray(x, dtype=float)
    return pd.DataFrame({
        "1": np.ones_like(x),
        "x": x,
        "x2": x**2,
        "x3": x**3,
        "sin_x": np.sin(x),
        "cos_x": np.cos(x),
        "log_abs_x": np.log1p(np.abs(x)),
    })


# ---------------------------------------------------------------------------
# Dataset processing
# ---------------------------------------------------------------------------


def process_dataset(dataset, path):
    coords = load_pca(path)
    if coords is None:
        return None, None

    # Keep the first three PCs as a compact state representation. Their absolute
    # coordinates remain study-specific; standardized versions are used only for
    # trajectory-shape comparisons.
    pcs = [c for c in ["PC1", "PC2", "PC3"] if c in coords.columns]
    out = coords[pcs].copy()
    out.insert(0, "sample", out.index.astype(str))
    out["dataset"] = dataset
    out["time_hours"] = [infer_time_hours(dataset, s) for s in out["sample"]]
    out["stage"] = [infer_stage(dataset, s, t) for s, t in zip(out["sample"], out["time_hours"])]
    out["replicate"] = [parse_replicate(dataset, s) for s in out["sample"]]

    for pc in pcs:
        out[f"{pc}_oriented"] = orient_pc(out[pc])
        out[f"{pc}_z"] = standardize_series(out[f"{pc}_oriented"])

    # Velocity is calculated only from real numeric time points and independently
    # for each dataset. Samples without explicit biological time (e.g. hESC/iPSC)
    # remain in the state table but are not assigned a fake time.
    order = out["time_hours"].notna()
    if order.any():
        timed = out.loc[order].sort_values("time_hours").copy()
        for pc in pcs:
            timed[f"d{pc}_dt"] = finite_difference(
                timed[f"{pc}_z"].to_numpy(), timed["time_hours"].to_numpy()
            )
        out = out.merge(
            timed[["sample"] + [f"d{pc}_dt" for pc in pcs]],
            on="sample",
            how="left",
        )
    else:
        for pc in pcs:
            out[f"d{pc}_dt"] = np.nan

    # A scalar speed in standardized PCA space.
    velocity_cols = [f"d{pc}_dt" for pc in pcs]
    out["state_speed"] = np.sqrt(np.nansum(out[velocity_cols].to_numpy() ** 2, axis=1))
    out.loc[out[velocity_cols].isna().all(axis=1), "state_speed"] = np.nan

    # Stability indicators are computed only for the ordered numeric-time
    # trajectory. They are deliberately labelled exploratory.
    out["rolling_variance_PC1"] = np.nan
    out["rolling_autocorrelation_PC1"] = np.nan
    if order.sum() >= 3:
        timed = out.loc[order].sort_values("time_hours")
        variance, autocorr = rolling_indicators(timed["PC1_z"].to_numpy(), window=5)
        out.loc[timed.index, "rolling_variance_PC1"] = variance.to_numpy()
        out.loc[timed.index, "rolling_autocorrelation_PC1"] = autocorr.to_numpy()

    return out, pcs


# ---------------------------------------------------------------------------
# Cross-study trajectory descriptors
# ---------------------------------------------------------------------------


def build_trajectory_summary(all_states):
    """Summarize each dataset using only quantities that are safe to compare."""
    rows = []
    for dataset, group in all_states.groupby("dataset"):
        timed = group[group["time_hours"].notna()].sort_values("time_hours")
        row = {
            "dataset": dataset,
            "n_samples": len(group),
            "n_timed_samples": len(timed),
            "n_unique_times": timed["time_hours"].nunique(),
            "time_min_hours": timed["time_hours"].min() if len(timed) else np.nan,
            "time_max_hours": timed["time_hours"].max() if len(timed) else np.nan,
            "n_replicate_labels": int((group["replicate"] != "unknown").sum()),
            "mean_state_speed": timed["state_speed"].mean() if len(timed) else np.nan,
            "max_state_speed": timed["state_speed"].max() if len(timed) else np.nan,
        }
        if len(timed) >= 2:
            row["net_PC1_z_change"] = timed["PC1_z"].iloc[-1] - timed["PC1_z"].iloc[0]
            row["net_PC2_z_change"] = (
                timed["PC2_z"].iloc[-1] - timed["PC2_z"].iloc[0]
                if "PC2_z" in timed else np.nan
            )
        else:
            row["net_PC1_z_change"] = np.nan
            row["net_PC2_z_change"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def build_symbolic_training_table(all_states):
    """Prepare time/state/derivative observations for future equation discovery."""
    cols = [
        "dataset", "sample", "stage", "replicate", "time_hours",
        "PC1_z", "dPC1_dt", "PC2_z", "dPC2_dt", "PC3_z", "dPC3_dt",
        "state_speed",
    ]
    available = [c for c in cols if c in all_states.columns]
    table = all_states[available].copy()
    table = table[table["time_hours"].notna()].copy()
    return table.sort_values(["dataset", "time_hours", "sample"])


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


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

    pd.DataFrame(availability).to_csv(OUT / "01_dataset_availability.csv", index=False)

    if not processed:
        print("No existing PCA trajectories found. Run the EDA scripts first.")
        return

    states = pd.concat(processed, ignore_index=True)
    states.to_csv(OUT / "02_time_aware_state_dynamics.csv", index=False)

    summary = build_trajectory_summary(states)
    summary.to_csv(OUT / "03_trajectory_summary.csv", index=False)

    training = build_symbolic_training_table(states)
    training.to_csv(OUT / "04_symbolic_training_table.csv", index=False)

    # Candidate library is generated from the standardized PC1 state. This is a
    # scaffold for equation discovery, not a fitted biological law.
    first_timed = training.dropna(subset=["PC1_z", "dPC1_dt"])
    if len(first_timed):
        library = candidate_library(first_timed["PC1_z"].to_numpy())
        library.insert(0, "dataset", first_timed["dataset"].to_numpy())
        library.insert(1, "sample", first_timed["sample"].to_numpy())
        library.insert(2, "time_hours", first_timed["time_hours"].to_numpy())
        library["target_dPC1_dt"] = first_timed["dPC1_dt"].to_numpy()
        library.to_csv(OUT / "05_symbolic_candidate_library.csv", index=False)

    report = [
        "Dynamics v2 — time-aware exploratory dynamics of OSKM reprogramming",
        "",
        "What changed from the original prototype",
        "------------------------------------------",
        "1. Biological time is parsed explicitly from sample identifiers where known.",
        "2. Finite-difference velocities use hours rather than arbitrary sample order.",
        "3. Replicate labels are retained instead of collapsing observations.",
        "4. PC1-PC3 are standardized within each study to remove arbitrary PCA scale.",
        "5. Raw PCA coordinates are NOT treated as comparable between studies.",
        "6. A symbolic-regression training table is prepared but no equation is fitted.",
        "",
        "Interpretation",
        "--------------",
        "The primary object is a time-dependent state trajectory within each study.",
        "Cross-study tables compare trajectory descriptors and standardized shape,",
        "not absolute PCA coordinates or raw expression values.",
        "",
        "Scientific target",
        "----------------",
        "Test whether independent OSKM/reprogramming experiments contain a robust",
        "low-dimensional temporal structure that can generalise across experiments.",
        "A stronger future test is to learn a dynamical equation on one complete",
        "study and validate its predictions on an entirely held-out study.",
        "",
        "Important limitations",
        "---------------------",
        "- PCA coordinates remain study-specific; z-scoring does not create a",
        "  biologically shared latent space.",
        "- Time points without explicit numeric elapsed time (e.g. iPSC endpoints)",
        "  are retained but excluded from numerical derivatives.",
        "- The current velocity is descriptive and sensitive to sampling density and",
        "  preprocessing choices.",
        "- Rolling variance/autocorrelation are exploratory indicators, not proof of",
        "  a critical transition or bifurcation.",
        "- No causal interpretation is made.",
        "- No claim is made for a quantum or relativistic biological mechanism.",
        "",
        "Next research layer",
        "-------------------",
        "1. Build a genuinely shared biological latent space from harmonised gene-level",
        "   features/markers rather than assuming PCA axes are shared.",
        "2. Model replicates explicitly and estimate uncertainty by bootstrap.",
        "3. Compare alternative state definitions and preprocessing pipelines.",
        "4. Fit candidate dynamical equations with symbolic regression / Genetic",
        "   Programming and compare against simpler baselines.",
        "5. Hold out GSE297234 or another complete study for external validation.",
        "6. Connect expression dynamics with GSE67520 regulatory dynamics later.",
    ]
    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")

    print(f"Dynamics v2 results written to: {OUT}")
    print(f"Datasets with PCA: {sum(a['PCA_file_found'] for a in availability)}/{len(availability)}")
    print(f"Time-aware observations: {len(training)}")


if __name__ == "__main__":
    main()
