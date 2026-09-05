"""Dynamics v2.2: time-aware exploratory dynamics of OSKM reprogramming.

Uses existing PCA outputs, explicit biological time, replicate-safe derivatives,
and study-internal standardized state coordinates. It does not assume PCA axes
are comparable between studies and does not fit a biological law.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "Dynamics"
OUT.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "07_PCA_coordinates.csv",
    "GSE148158": RESULTS / "GSE148158" / "07_PCA_coordinates.csv",
    "GSE52052": RESULTS / "GSE52052" / "08_PCA_coordinates.csv",
    "GSE67462": RESULTS / "GSE67462" / "09_PCA_coordinates.csv",
    "GSE297234": RESULTS / "GSE297234" / "08_PCA_coordinates.csv",
}

# GEO sample IDs whose biological time is known exactly.
GSM_TIME = {
    # GSE148158
    "GSM4455240": 48.0, "GSM4455241": 48.0,
    "GSM4455242": 72.0, "GSM4455243": 72.0,
    "GSM4455244": 48.0, "GSM4455245": 72.0,
    # GSE28688
    "GSM710515": 24.0, "GSM710516": 24.0,
    "GSM710517": 48.0, "GSM710518": 48.0,
    "GSM710519": 72.0, "GSM710520": 72.0,
    # GSE52052: all samples are day 11
    "GSM1258008": 264.0, "GSM1258009": 264.0, "GSM1258010": 264.0,
    "GSM1258011": 264.0, "GSM1258012": 264.0, "GSM1258013": 264.0,
    # GSE67462: two replicates per stage
    "GSM1647454": 0.0, "GSM1647455": 0.0,
    "GSM1647456": 24.0, "GSM1647457": 24.0,
    "GSM1647458": 72.0, "GSM1647459": 72.0,
    "GSM1647460": 120.0, "GSM1647461": 120.0,
    "GSM1647462": 168.0, "GSM1647463": 168.0,
    "GSM1647464": 264.0, "GSM1647465": 264.0,
    "GSM1647466": 360.0, "GSM1647467": 360.0,
    "GSM1647468": 432.0, "GSM1647469": 432.0,
    # GSE297234
    "GSM8986586": 0.0, "GSM8986587": 72.0,
    "GSM8986588": 168.0, "GSM8986589": 240.0,
    "GSM8986590": 0.0, "GSM8986591": 72.0,
    "GSM8986592": 168.0, "GSM8986593": 240.0,
}

TIME_PATTERNS = {
    "GSE28688": [(r"24\s*h", 24.0), (r"48\s*h", 48.0), (r"72\s*h", 72.0)],
    "GSE148158": [(r"48", 48.0), (r"72", 72.0)],
    "GSE52052": [(r"day\s*11", 264.0)],
    "GSE67462": [
        (r"day\s*0\b", 0.0), (r"day\s*1\b", 24.0), (r"day\s*3\b", 72.0),
        (r"day\s*5\b", 120.0), (r"day\s*7\b", 168.0),
        (r"day\s*11\b", 264.0), (r"day\s*15\b", 360.0),
        (r"day\s*18\b", 432.0),
    ],
    "GSE297234": [
        (r"d0\b|day\s*0\b", 0.0), (r"d3\b|day\s*3\b", 72.0),
        (r"d7\b|day\s*7\b", 168.0), (r"d10\b|day\s*10\b", 240.0),
    ],
}


def load_pca(path):
    if not path.exists():
        return None
    df = pd.read_csv(path, index_col=0)
    pcs = [c for c in ("PC1", "PC2", "PC3") if c in df.columns]
    if not pcs:
        return None
    out = df[pcs].apply(pd.to_numeric, errors="coerce")
    out.index = out.index.astype(str)
    return out


def time_hours(dataset, sample):
    s = str(sample).strip().strip('"')
    if s in GSM_TIME:
        return GSM_TIME[s]
    text = s.lower().replace("_", " ").replace("-", " ")
    for pattern, value in TIME_PATTERNS.get(dataset, []):
        if re.search(pattern, text):
            return value
    return np.nan


def stage(sample, t):
    text = str(sample).lower()
    if re.search(r"ipsc|\bips\b", text):
        return "iPSC"
    if re.search(r"hesc", text):
        return "hESC"
    if pd.notna(t):
        return f"day{int(t / 24)}" if t % 24 == 0 else f"{int(t)}h"
    return "unknown"


def replicate(sample):
    s = str(sample)
    m = re.search(r"(?:-|_|\s)([ab])$", s, re.I)
    if m:
        return m.group(1).lower()
    m = re.search(r"(?:rep|replicate)[_\s-]*(\d+)", s, re.I)
    if m:
        return m.group(1)
    # GSE67462 has paired GSMs: even accession = first replicate, odd = second.
    m = re.fullmatch(r"GSM(\d+)", s)
    if m and 1647454 <= int(m.group(1)) <= 1647469:
        return "1" if int(m.group(1)) % 2 == 0 else "2"
    return "unknown"


def zscore(x):
    x = pd.Series(x, dtype=float)
    sd = x.std(ddof=0)
    return (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else pd.Series(np.nan, index=x.index)


def orient(x):
    x = pd.Series(x, dtype=float).copy()
    valid = x.dropna()
    if valid.empty:
        return x
    i = valid.abs().idxmax()
    return -x if x.loc[i] < 0 else x


def derivative(x, t):
    """Derivative for strictly increasing unique biological times."""
    x = np.asarray(x, float)
    t = np.asarray(t, float)
    if len(x) < 2:
        return np.full(len(x), np.nan)
    if np.any(np.diff(t) <= 0):
        raise ValueError("Duplicate/non-increasing time reached derivative().")
    if len(x) == 2:
        v = (x[1] - x[0]) / (t[1] - t[0])
        return np.array([v, v])
    return np.gradient(x, t)


def process(dataset, path):
    coords = load_pca(path)
    if coords is None:
        return None, []
    pcs = [c for c in ("PC1", "PC2", "PC3") if c in coords.columns]
    out = coords[pcs].copy()
    out.insert(0, "sample", out.index.astype(str))
    out["dataset"] = dataset
    out["time_hours"] = [time_hours(dataset, s) for s in out["sample"]]
    out["stage"] = [stage(s, t) for s, t in zip(out["sample"], out["time_hours"])]
    out["replicate"] = [replicate(s) for s in out["sample"]]

    for pc in pcs:
        out[f"{pc}_z"] = zscore(orient(out[pc]))

    timed = out[out["time_hours"].notna()].copy()
    if timed.empty:
        for pc in pcs:
            out[f"d{pc}_dt"] = np.nan
        out["state_speed"] = np.nan
        out["rolling_variance_PC1"] = np.nan
        out["rolling_autocorrelation_PC1"] = np.nan
        return out, pcs

    # Replicate-safe: first average the state at identical biological times.
    zcols = [f"{pc}_z" for pc in pcs]
    mean_state = timed.groupby("time_hours", as_index=False)[zcols].mean().sort_values("time_hours")
    if len(mean_state) >= 2:
        for pc in pcs:
            mean_state[f"d{pc}_dt"] = derivative(
                mean_state[f"{pc}_z"].to_numpy(),
                mean_state["time_hours"].to_numpy(),
            )
        vcols = [f"d{pc}_dt" for pc in pcs]
        mean_state["state_speed"] = np.sqrt(np.sum(mean_state[vcols].to_numpy() ** 2, axis=1))
    else:
        for pc in pcs:
            mean_state[f"d{pc}_dt"] = np.nan
        mean_state["state_speed"] = np.nan

    mean_state["rolling_variance_PC1"] = np.nan
    mean_state["rolling_autocorrelation_PC1"] = np.nan
    if len(mean_state) >= 3:
        s = mean_state["PC1_z"]
        w = min(5, len(s))
        mean_state["rolling_variance_PC1"] = s.rolling(w, min_periods=3).var().to_numpy()
        mean_state["rolling_autocorrelation_PC1"] = s.rolling(w, min_periods=3).apply(
            lambda q: q.autocorr(lag=1) if q.std() > 0 else np.nan, raw=False
        ).to_numpy()

    merge_cols = ["time_hours"] + [f"d{pc}_dt" for pc in pcs] + [
        "state_speed", "rolling_variance_PC1", "rolling_autocorrelation_PC1"
    ]
    return out.merge(mean_state[merge_cols], on="time_hours", how="left"), pcs


def summary(states):
    rows = []
    for dataset, g in states.groupby("dataset"):
        t = g[g["time_hours"].notna()]
        rows.append({
            "dataset": dataset,
            "n_samples": len(g),
            "n_timed_samples": len(t),
            "n_unique_times": t["time_hours"].nunique(),
            "time_min_hours": t["time_hours"].min() if len(t) else np.nan,
            "time_max_hours": t["time_hours"].max() if len(t) else np.nan,
            "replicate_labelled": int((g["replicate"] != "unknown").sum()),
            "mean_state_speed": t["state_speed"].mean() if len(t) else np.nan,
            "max_state_speed": t["state_speed"].max() if len(t) else np.nan,
        })
    return pd.DataFrame(rows)


def main():
    availability, parts = [], []
    for dataset, path in DATASETS.items():
        state, pcs = process(dataset, path)
        found = state is not None
        availability.append({
            "dataset": dataset,
            "PCA_file_found": found,
            "path": str(path),
            "PCs_available": ",".join(pcs),
            "timed_samples": int(state["time_hours"].notna().sum()) if found else 0,
            "unique_times": int(state["time_hours"].dropna().nunique()) if found else 0,
        })
        if found:
            parts.append(state)

    pd.DataFrame(availability).to_csv(OUT / "01_dataset_availability.csv", index=False)
    if not parts:
        print("No PCA trajectories found.")
        return

    states = pd.concat(parts, ignore_index=True)
    states.to_csv(OUT / "02_time_aware_state_dynamics.csv", index=False)
    summary(states).to_csv(OUT / "03_trajectory_summary.csv", index=False)

    cols = ["dataset", "sample", "stage", "replicate", "time_hours",
            "PC1_z", "dPC1_dt", "PC2_z", "dPC2_dt", "PC3_z", "dPC3_dt", "state_speed"]
    training = states[[c for c in cols if c in states.columns]]
    training = training[training["time_hours"].notna()].sort_values(["dataset", "time_hours", "sample"])
    training.to_csv(OUT / "04_symbolic_training_table.csv", index=False)

    # Candidate library uses one row per unique time point, not one row per replicate.
    source = training.groupby(["dataset", "time_hours"], as_index=False).agg(
        PC1_z=("PC1_z", "mean"), dPC1_dt=("dPC1_dt", "mean")
    ).dropna()
    if len(source):
        x = source["PC1_z"].to_numpy(float)
        lib = pd.DataFrame({
            "dataset": source["dataset"], "time_hours": source["time_hours"],
            "constant": 1.0, "x": x, "x2": x**2, "x3": x**3,
            "sin_x": np.sin(x), "cos_x": np.cos(x),
            "log_abs_x": np.log1p(np.abs(x)),
            "target_dPC1_dt": source["dPC1_dt"],
        })
        lib.to_csv(OUT / "05_symbolic_candidate_library.csv", index=False)

    report = [
        "Dynamics v2.2 — time-aware exploratory dynamics of OSKM reprogramming",
        "",
        "Biological time is recovered from GEO sample IDs and descriptive names.",
        "Replicates at identical time points are averaged before derivatives are estimated.",
        "PCA coordinates are standardized within each study; PCA axes are not assumed comparable between studies.",
        "The symbolic-regression table is a preparation layer, not a fitted biological equation.",
        "",
        "Research target: test whether independent reprogramming experiments share a robust",
        "low-dimensional temporal dynamical structure and whether a model learned on one",
        "complete experiment can predict an entirely held-out experiment.",
        "",
        "Limitations: derivatives depend on sampling density and preprocessing; endpoints such as",
        "iPSC/hESC without explicit elapsed time are excluded from numerical derivatives; current",
        "PC coordinates do not constitute a shared biological latent space; stability indicators",
        "are hypothesis-generating only; no causal, quantum, or relativistic mechanism is claimed.",
    ]
    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")

    print(f"Dynamics v2.2 results written to: {OUT}")
    print(f"Datasets with PCA: {len(parts)}/{len(DATASETS)}")
    print(f"Time-aware observations: {len(training)}")
    print("Replicate-safe derivative calculation: enabled")
    print("Dataset timing summary:")
    for row in availability:
        print(f"  {row['dataset']}: PCA={row['PCA_file_found']}, timed={row['timed_samples']}, unique_times={row['unique_times']}")


if __name__ == "__main__":
    main()
