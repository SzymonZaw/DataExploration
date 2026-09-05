"""Dynamics v3.0: staged exploratory dynamics pipeline for OSKM reprogramming.

Stages implemented in this single script:
1. Data integration and metadata harmonisation.
2. Study-normalized latent state representation.
3. Time-aware trajectory reconstruction.
4. Derivative-based dynamics.
5. Critical-transition / stability indicators.

Important scientific boundary:
- PCA coordinates are standardized within each study.
- Stage 2 is NOT yet a cross-study biological latent space.
- No causal law, quantum mechanism, or relativistic mechanism is claimed.
- Symbolic regression is not fitted here; its training table is prepared for a later stage.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "Dynamics"
OUT.mkdir(parents=True, exist_ok=True)

STAGE_DIRS = {i: OUT / f"stage{i}" for i in range(1, 6)}
for directory in STAGE_DIRS.values():
    directory.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "07_PCA_coordinates.csv",
    "GSE148158": RESULTS / "GSE148158" / "07_PCA_coordinates.csv",
    "GSE52052": RESULTS / "GSE52052" / "08_PCA_coordinates.csv",
    "GSE67462": RESULTS / "GSE67462" / "09_PCA_coordinates.csv",
    "GSE297234": RESULTS / "GSE297234" / "08_PCA_coordinates.csv",
}

GSM_TIME = {
    "GSM4455240": 48.0, "GSM4455241": 48.0,
    "GSM4455242": 72.0, "GSM4455243": 72.0,
    "GSM4455244": 48.0, "GSM4455245": 72.0,
    "GSM710515": 24.0, "GSM710516": 24.0,
    "GSM710517": 48.0, "GSM710518": 48.0,
    "GSM710519": 72.0, "GSM710520": 72.0,
    "GSM1258008": 264.0, "GSM1258009": 264.0, "GSM1258010": 264.0,
    "GSM1258011": 264.0, "GSM1258012": 264.0, "GSM1258013": 264.0,
    "GSM1647454": 0.0, "GSM1647455": 0.0,
    "GSM1647456": 24.0, "GSM1647457": 24.0,
    "GSM1647458": 72.0, "GSM1647459": 72.0,
    "GSM1647460": 120.0, "GSM1647461": 120.0,
    "GSM1647462": 168.0, "GSM1647463": 168.0,
    "GSM1647464": 264.0, "GSM1647465": 264.0,
    "GSM1647466": 360.0, "GSM1647467": 360.0,
    "GSM1647468": 432.0, "GSM1647469": 432.0,
    "GSM8986586": 0.0, "GSM8986587": 72.0,
    "GSM8986588": 168.0, "GSM8986589": 240.0,
    "GSM8986590": 0.0, "GSM8986591": 72.0,
    "GSM8986592": 168.0, "GSM8986593": 240.0,
}

GSE28688_ROW_TIME = [
    0.0, 0.0, 24.0, 24.0, 48.0, 48.0, 72.0, 72.0,
    np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,
]
GSE28688_ROW_SAMPLE = [
    "GSM710513", "GSM710514", "GSM710515", "GSM710516",
    "GSM710517", "GSM710518", "GSM710519", "GSM710520",
    "GSM710521", "GSM710522", "GSM710523", "GSM710524",
    "GSM710525", "GSM710526",
]

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


def stage_label(sample, t):
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


def stage1_data_integration():
    """Stage 1: integrate PCA outputs and harmonise sample metadata."""
    availability, parts = [], []
    for dataset, path in DATASETS.items():
        coords = load_pca(path)
        if coords is None:
            availability.append({
                "dataset": dataset, "PCA_file_found": False, "n_samples": 0,
                "n_timed_samples": 0, "n_unique_times": 0, "role": "unavailable",
                "timing_source": "none", "path": str(path),
            })
            continue

        pcs = [c for c in ("PC1", "PC2", "PC3") if c in coords.columns]
        out = coords[pcs].copy()
        out.insert(0, "sample", out.index.astype(str))
        timing_source = "GSM_or_text"

        if dataset == "GSE28688" and len(out) == len(GSE28688_ROW_SAMPLE):
            out["sample"] = GSE28688_ROW_SAMPLE
            timing_source = "GSE28688_GEO_row_order"

        out["dataset"] = dataset
        out["time_hours"] = [time_hours(dataset, s) for s in out["sample"]]
        if dataset == "GSE28688" and timing_source == "GSE28688_GEO_row_order":
            out["time_hours"] = GSE28688_ROW_TIME
        out["stage"] = [stage_label(s, t) for s, t in zip(out["sample"], out["time_hours"])]
        out["replicate"] = [replicate(s) for s in out["sample"]]
        out["timing_source"] = timing_source
        for pc in pcs:
            out[f"{pc}_z"] = zscore(orient(out[pc]))

        timed = out[out["time_hours"].notna()]
        role = "trajectory" if timed["time_hours"].nunique() >= 2 else "context_only"
        availability.append({
            "dataset": dataset,
            "PCA_file_found": True,
            "n_samples": len(out),
            "n_timed_samples": len(timed),
            "n_unique_times": timed["time_hours"].nunique(),
            "role": role,
            "timing_source": timing_source,
            "path": str(path),
        })
        parts.append(out)

    availability_df = pd.DataFrame(availability)
    states = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    availability_df.to_csv(STAGE_DIRS[1] / "01_dataset_availability.csv", index=False)
    states.to_csv(STAGE_DIRS[1] / "02_master_sample_metadata.csv", index=False)
    return states, availability_df


def stage2_latent_state(states):
    """Stage 2: create a study-normalized low-dimensional state representation."""
    if states.empty:
        return states.copy()

    out = states.copy()
    zcols = [c for c in ("PC1_z", "PC2_z", "PC3_z") if c in out.columns]
    out["latent_1"] = out["PC1_z"] if "PC1_z" in out.columns else np.nan
    out["latent_2"] = out["PC2_z"] if "PC2_z" in out.columns else np.nan
    out["latent_3"] = out["PC3_z"] if "PC3_z" in out.columns else np.nan
    out["latent_space_type"] = "study_normalized_PCA"

    out[["dataset", "sample", "time_hours", "stage", "replicate",
         "latent_1", "latent_2", "latent_3", "latent_space_type"]].to_csv(
        STAGE_DIRS[2] / "01_latent_state_coordinates.csv", index=False
    )
    return out


def stage3_trajectory_reconstruction(states):
    """Stage 3: reconstruct time-ordered trajectories and replicate means."""
    if states.empty:
        return states.copy(), pd.DataFrame()

    timed = states[states["time_hours"].notna()].copy()
    trajectory_rows = []
    for dataset, group in timed.groupby("dataset"):
        role = "trajectory" if group["time_hours"].nunique() >= 2 else "context_only"
        if role != "trajectory":
            continue
        zcols = [c for c in ("latent_1", "latent_2", "latent_3") if c in group.columns]
        mean_state = group.groupby("time_hours", as_index=False)[zcols].mean().sort_values("time_hours")
        mean_state.insert(0, "dataset", dataset)
        mean_state["n_replicates"] = group.groupby("time_hours").size().reindex(
            mean_state["time_hours"]
        ).to_numpy()
        trajectory_rows.append(mean_state)

    trajectories = pd.concat(trajectory_rows, ignore_index=True) if trajectory_rows else pd.DataFrame()
    trajectories.to_csv(STAGE_DIRS[3] / "01_reconstructed_trajectories.csv", index=False)

    summary_rows = []
    for dataset, g in trajectories.groupby("dataset") if not trajectories.empty else []:
        summary_rows.append({
            "dataset": dataset,
            "n_timepoints": g["time_hours"].nunique(),
            "time_min_hours": g["time_hours"].min(),
            "time_max_hours": g["time_hours"].max(),
            "total_replicate_observations": int(g["n_replicates"].sum()),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(STAGE_DIRS[3] / "02_trajectory_summary.csv", index=False)
    return states, trajectories


def stage4_dynamics(trajectories):
    """Stage 4: calculate temporal derivatives, speed and acceleration."""
    if trajectories.empty:
        return trajectories.copy()

    out = trajectories.copy()
    for dataset, idx in out.groupby("dataset").groups.items():
        sub = out.loc[idx].sort_values("time_hours")
        t = sub["time_hours"].to_numpy(float)
        for axis in ("latent_1", "latent_2", "latent_3"):
            if axis not in sub.columns:
                continue
            x = sub[axis].to_numpy(float)
            if len(x) >= 2:
                v = derivative(x, t)
                a = derivative(v, t) if len(x) >= 3 else np.full(len(x), np.nan)
            else:
                v = np.full(len(x), np.nan)
                a = np.full(len(x), np.nan)
            out.loc[sub.index, f"d{axis}_dt"] = v
            out.loc[sub.index, f"d2{axis}_dt2"] = a

    vcols = [c for c in ("dlatent_1_dt", "dlatent_2_dt", "dlatent_3_dt") if c in out.columns]
    acols = [c for c in ("d2latent_1_dt2", "d2latent_2_dt2", "d2latent_3_dt2") if c in out.columns]
    out["state_speed"] = np.sqrt(np.nansum(out[vcols].to_numpy(float) ** 2, axis=1)) if vcols else np.nan
    out["state_acceleration"] = np.sqrt(np.nansum(out[acols].to_numpy(float) ** 2, axis=1)) if acols else np.nan
    out.to_csv(STAGE_DIRS[4] / "01_dynamics.csv", index=False)

    training = out[["dataset", "time_hours", "latent_1", "dlatent_1_dt",
                    "latent_2", "dlatent_2_dt", "latent_3", "dlatent_3_dt",
                    "state_speed", "state_acceleration"]].copy()
    training.to_csv(STAGE_DIRS[4] / "02_symbolic_training_table.csv", index=False)
    return out


def stage5_critical_transitions(dynamics):
    """Stage 5: calculate exploratory stability and critical-transition indicators."""
    if dynamics.empty:
        return dynamics.copy()

    rows = []
    for dataset, group in dynamics.groupby("dataset"):
        g = group.sort_values("time_hours").copy()
        x = g["latent_1"].astype(float)
        n = len(g)
        window = min(5, n)
        g["rolling_variance"] = x.rolling(window, min_periods=3).var()
        g["rolling_autocorrelation"] = x.rolling(window, min_periods=3).apply(
            lambda q: q.autocorr(lag=1) if q.std() > 0 else np.nan, raw=False
        )
        g["trajectory_curvature_proxy"] = np.nan
        if {"dlatent_1_dt", "dlatent_2_dt"}.issubset(g.columns):
            v1 = g["dlatent_1_dt"].to_numpy(float)
            v2 = g["dlatent_2_dt"].to_numpy(float)
            angle = np.arctan2(v2, v1)
            g["trajectory_curvature_proxy"] = np.abs(np.gradient(angle)) if n >= 2 else np.nan
        rows.append(g)

    out = pd.concat(rows, ignore_index=True)
    out.to_csv(STAGE_DIRS[5] / "01_critical_transition_indicators.csv", index=False)

    critical_rows = []
    for dataset, g in out.groupby("dataset"):
        score_columns = ["rolling_variance", "rolling_autocorrelation", "trajectory_curvature_proxy"]
        candidates = []
        for col in score_columns:
            if col not in g.columns or g[col].dropna().empty:
                continue
            s = g[col]
            idx = s.idxmax()
            if pd.notna(s.loc[idx]):
                candidates.append({
                    "dataset": dataset,
                    "indicator": col,
                    "candidate_time_hours": g.loc[idx, "time_hours"],
                    "indicator_value": s.loc[idx],
                })
        critical_rows.extend(candidates)

    candidates_df = pd.DataFrame(critical_rows)
    candidates_df.to_csv(STAGE_DIRS[5] / "02_candidate_critical_points.csv", index=False)
    return out


def write_report(availability):
    rows = availability.to_dict("records") if not availability.empty else []
    report = [
        "Dynamics v3.0 — staged OSKM reprogramming dynamics pipeline", "",
        "STAGE 1: Data integration and metadata harmonisation.",
        "STAGE 2: Study-normalized low-dimensional latent state representation.",
        "STAGE 3: Time-aware trajectory reconstruction with replicate averaging.",
        "STAGE 4: Derivatives, state speed and acceleration.",
        "STAGE 5: Exploratory critical-transition and stability indicators.", "",
        "Dataset timing summary:",
    ]
    for row in rows:
        report.append(
            f"  {row['dataset']}: role={row['role']}, samples={row['n_samples']}, "
            f"timed={row['n_timed_samples']}, unique_times={row['n_unique_times']}, "
            f"timing={row['timing_source']}"
        )
    report += [
        "",
        "Scientific boundary:",
        "Stage 2 is study-normalized PCA, not a validated cross-study biological latent space.",
        "Independent datasets must not be treated as directly comparable PC coordinate systems.",
        "Stage 5 indicators are hypothesis-generating and do not establish bifurcation or causality.",
        "The symbolic training table is preparation for later symbolic regression, not a fitted law.",
        "No quantum or relativistic biological mechanism is assumed.",
    ]
    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")


def main():
    print("=== STAGE 1: DATA INTEGRATION ===")
    states, availability = stage1_data_integration()

    print("=== STAGE 2: LATENT STATE REPRESENTATION ===")
    states = stage2_latent_state(states)

    print("=== STAGE 3: TRAJECTORY RECONSTRUCTION ===")
    states, trajectories = stage3_trajectory_reconstruction(states)

    print("=== STAGE 4: DYNAMICS ===")
    dynamics = stage4_dynamics(trajectories)

    print("=== STAGE 5: CRITICAL TRANSITIONS ===")
    stage5_critical_transitions(dynamics)

    write_report(availability)

    print(f"Dynamics v3.0 results written to: {OUT}")
    print(f"Datasets with PCA: {int(availability['PCA_file_found'].sum())}/{len(DATASETS)}")
    print(f"Time-aware observations: {int(states['time_hours'].notna().sum()) if not states.empty else 0}")
    print(f"Trajectory datasets: {int((availability['role'] == 'trajectory').sum())}")
    print("Replicate-safe derivative calculation: enabled")
    print("Stages 1-5: completed")
    print("Dataset timing summary:")
    for row in availability.to_dict("records"):
        print(
            f"  {row['dataset']}: PCA={row['PCA_file_found']}, role={row['role']}, "
            f"timed={row['n_timed_samples']}, unique_times={row['n_unique_times']}, "
            f"timing={row['timing_source']}"
        )


if __name__ == "__main__":
    main()
