"""Stage 2.7 validation of the biological common gene space.

The primary test is leave-one-dataset-out prediction in the Stage 2.6 common
human gene space. Replicate holdout and a time-permutation null are reported
as secondary controls.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
STAGE26 = ROOT / "results" / "Dynamics" / "stage2_6"
OUT = ROOT / "results" / "Dynamics" / "stage2_7"
OUT.mkdir(parents=True, exist_ok=True)


def _metrics(y_true, y_pred):
    a = np.asarray(y_true, dtype=float).ravel()
    b = np.asarray(y_pred, dtype=float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if not ok.any():
        return {"rmse": np.nan, "mae": np.nan, "correlation": np.nan}
    a, b = a[ok], b[ok]
    corr = np.corrcoef(a, b)[0, 1] if len(a) > 1 and np.std(a) > 0 and np.std(b) > 0 else np.nan
    return {"rmse": float(np.sqrt(np.mean((a - b) ** 2))), "mae": float(np.mean(np.abs(a - b))), "correlation": float(corr)}


def _interpolate(points, times, target):
    times = np.asarray(times, dtype=float)
    points = np.asarray(points, dtype=float)
    if len(times) < 2 or target < times.min() or target > times.max():
        return None
    return np.asarray([np.interp(target, times, points[:, j]) for j in range(points.shape[1])])


def _normalise_name(value):
    s = str(value).strip().strip('"').strip().replace("\\", "/")
    return re.sub(r"\s+", "", s).lower()


def _candidate_column_names(dataset, sample):
    ds, sm = str(dataset).strip().strip('"'), str(sample).strip().strip('"')
    return [sm, f"{ds}__{sm}", f"{ds}_{sm}", f"{ds}/{sm}"]


def _build_matrix_column_map(matrix, metadata):
    actual = list(matrix.columns)
    normalised, ambiguous = {}, set()
    for col in actual:
        key = _normalise_name(col)
        if key in normalised and normalised[key] != col:
            ambiguous.add(key)
        else:
            normalised[key] = col
    resolved, statuses = [], []
    for _, row in metadata.iterrows():
        found, method = None, "unmatched"
        for candidate in _candidate_column_names(row["dataset"], row["sample"]):
            key = _normalise_name(candidate)
            if key in normalised and key not in ambiguous:
                found, method = normalised[key], "exact_or_alias"
                break
        resolved.append(found)
        statuses.append(method)
    metadata = metadata.copy()
    metadata["matrix_column"] = resolved
    metadata["matrix_match_status"] = statuses
    return metadata


def _print_mapping_diagnostics(matrix, metadata):
    rows = []
    for ds, g in metadata.groupby("dataset", sort=True):
        timed = g[g["time_hours"].notna()]
        matched = g[g["matrix_column"].notna()]
        timed_matched = timed[timed["matrix_column"].notna()]
        rows.append({
            "dataset": ds,
            "matrix_columns": int(timed_matched["matrix_column"].nunique()),
            "metadata_samples": len(g),
            "matched_samples": len(matched),
            "timed_samples": len(timed),
            "timed_matched": len(timed_matched),
            "unique_times": int(timed_matched["time_hours"].nunique()),
            "replicates": ",".join(sorted(map(str, timed_matched["replicate"].dropna().unique()))),
        })
    diag = pd.DataFrame(rows)
    diag.to_csv(OUT / "00_mapping_diagnostics.csv", index=False)
    print("\nStage 2.7 sample-to-matrix mapping:")
    print(diag.to_string(index=False))
    return diag


def _strip_dataset_prefix(sample):
    """Recover the original sample ID from Stage 2.6's dataset__sample label."""
    s = str(sample).strip().strip('"')
    return s.split("__", 1)[1] if "__" in s else s


def _time_hours_for_validation(dataset, sample, gse28688_row_index=None):
    """Apply Dynamics.py time parsing to the original sample ID.

    Stage 2.6 writes columns as ``dataset__sample``. The original Dynamics
    time parser expects the underlying GSM/text label, so validation must strip
    that technical prefix before parsing. GSE28688 additionally has four
    untimed GEO rows whose time points are defined by the documented row order.
    """
    from Dynamics import time_hours
    raw = _strip_dataset_prefix(sample)
    if dataset == "GSE28688" and gse28688_row_index is not None:
        from Dynamics import GSE28688_ROW_TIME
        if gse28688_row_index < len(GSE28688_ROW_TIME):
            row_time = GSE28688_ROW_TIME[gse28688_row_index]
            if pd.notna(row_time):
                return float(row_time)
    return time_hours(dataset, raw)


def _load_common_space():
    matrix_path = STAGE26 / "06_common_human_gene_matrix.csv"
    meta_path = STAGE26 / "07_common_gene_sample_metadata.csv"
    if not matrix_path.exists() or not meta_path.exists():
        raise FileNotFoundError("Stage 2.6 outputs are missing; run Dynamics.py first.")

    matrix = pd.read_csv(matrix_path, index_col=0).apply(pd.to_numeric, errors="coerce")
    metadata = pd.read_csv(meta_path)
    metadata["dataset"] = metadata["dataset"].astype(str)
    metadata["sample"] = metadata["sample"].astype(str)

    from Dynamics import condition, replicate

    times = []
    conditions = []
    replicates = []
    for ds, sample in zip(metadata["dataset"], metadata["sample"]):
        raw = _strip_dataset_prefix(sample)
        idx = None
        if ds == "GSE28688":
            m = re.fullmatch(r"GSM(\d+)", raw)
            if m:
                idx = int(m.group(1)) - 710513
        times.append(_time_hours_for_validation(ds, sample, idx))
        conditions.append(condition(ds, raw))
        replicates.append(replicate(raw))
    metadata["time_hours"] = times
    metadata["condition"] = conditions
    metadata["replicate"] = replicates
    metadata = _build_matrix_column_map(matrix, metadata)
    _print_mapping_diagnostics(matrix, metadata)

    if metadata["matrix_column"].notna().sum() == 0:
        raise RuntimeError("Stage 2.7 could not match any metadata sample to the Stage 2.6 matrix columns. See results/Dynamics/stage2_7/00_mapping_diagnostics.csv.")
    return matrix, metadata[metadata["matrix_column"].notna()].copy()


def _trajectories(matrix, metadata, time_override=None):
    out = {}
    for ds, g in metadata.groupby("dataset"):
        g = g[g["matrix_column"].notna() & g["time_hours"].notna()].copy()
        if len(g) < 2:
            continue
        times = g["time_hours"].astype(float).to_numpy()
        if time_override:
            times = np.asarray([time_override.get((ds, c), t) for c, t in zip(g["matrix_column"], times)], dtype=float)
        frame = pd.DataFrame(matrix[g["matrix_column"]].T.to_numpy(), index=times, columns=matrix.index).groupby(level=0).mean().sort_index()
        if len(frame) >= 2:
            out[ds] = (frame.index.to_numpy(float), frame.to_numpy(float))
    return out


def leave_one_dataset_out(matrix, metadata):
    """Predict each held-out dataset from all other datasets at overlapping real times."""
    traj = _trajectories(matrix, metadata)
    rows = []
    for test_ds, (test_times, _) in sorted(traj.items()):
        for target in test_times:
            preds = []
            for train_ds, (times, values) in traj.items():
                if train_ds == test_ds:
                    continue
                p = _interpolate(values, times, target)
                if p is not None:
                    preds.append(p)
            if not preds:
                continue
            cols = metadata[(metadata["dataset"] == test_ds) & (metadata["time_hours"] == target) & metadata["matrix_column"].notna()]["matrix_column"].tolist()
            if not cols:
                continue
            truth = matrix[cols].mean(axis=1).to_numpy(float)
            rows.append({"validation": "leave_one_dataset_out", "test_dataset": test_ds, "time_hours": float(target), "n_training_datasets": len(preds), **_metrics(truth, np.mean(preds, axis=0))})
    return pd.DataFrame(rows)


def leave_one_replicate_out(matrix, metadata):
    """Predict a held-out replicate from other replicates in the same dataset."""
    rows = []
    for ds, g in metadata.groupby("dataset"):
        reps = sorted(r for r in g["replicate"].dropna().unique() if str(r).lower() != "unknown")
        if len(reps) < 2:
            continue
        for held in reps:
            train = g[(g["replicate"] != held) & g["matrix_column"].notna()].copy()
            test = g[(g["replicate"] == held) & g["matrix_column"].notna()].copy()
            train = train[train["time_hours"].notna()]
            for target in sorted(test["time_hours"].dropna().unique()):
                frame = pd.DataFrame(matrix[train["matrix_column"]].T.to_numpy(), index=train["time_hours"].to_numpy(), columns=matrix.index).groupby(level=0).mean().sort_index()
                pred = _interpolate(frame.to_numpy(), frame.index.to_numpy(float), float(target))
                if pred is None:
                    continue
                cols = test[test["time_hours"] == target]["matrix_column"].tolist()
                truth = matrix[cols].mean(axis=1).to_numpy(float)
                rows.append({"validation": "leave_one_replicate_out", "dataset": ds, "held_out_replicate": str(held), "time_hours": float(target), **_metrics(truth, pred)})
    return pd.DataFrame(rows)


def permutation_null(matrix, metadata, n_permutations=25, seed=42):
    """Shuffle training time labels within datasets and repeat holdout prediction."""
    rng = np.random.default_rng(seed)
    rows = []
    for permutation in range(n_permutations):
        override = {}
        for ds, g in metadata.groupby("dataset"):
            g = g[g["time_hours"].notna() & g["matrix_column"].notna()]
            vals = g["time_hours"].to_numpy(float)
            shuffled = rng.permutation(vals)
            for col, value in zip(g["matrix_column"], shuffled):
                override[(ds, col)] = float(value)
        traj = _trajectories(matrix, metadata, override)
        for test_ds, _ in sorted(traj.items()):
            test = metadata[(metadata["dataset"] == test_ds) & metadata["time_hours"].notna() & metadata["matrix_column"].notna()]
            for target in sorted(test["time_hours"].unique()):
                preds = []
                for train_ds, (times, values) in traj.items():
                    if train_ds == test_ds:
                        continue
                    p = _interpolate(values, times, float(target))
                    if p is not None:
                        preds.append(p)
                if not preds:
                    continue
                cols = test[test["time_hours"] == target]["matrix_column"].tolist()
                truth = matrix[cols].mean(axis=1).to_numpy(float)
                rows.append({"validation": "time_permutation_null", "permutation": permutation, "test_dataset": test_ds, "time_hours": float(target), **_metrics(truth, np.mean(preds, axis=0))})
    return pd.DataFrame(rows)


def stage2_7(n_permutations=25, seed=42):
    matrix, metadata = _load_common_space()
    dataset_df = leave_one_dataset_out(matrix, metadata)
    replicate_df = leave_one_replicate_out(matrix, metadata)
    null_df = permutation_null(matrix, metadata, n_permutations=n_permutations, seed=seed)
    dataset_df.to_csv(OUT / "01_leave_one_dataset_out.csv", index=False)
    replicate_df.to_csv(OUT / "02_leave_one_replicate_out.csv", index=False)
    null_df.to_csv(OUT / "03_time_permutation_null.csv", index=False)
    frames = [("leave_one_dataset_out", dataset_df), ("leave_one_replicate_out", replicate_df), ("time_permutation_null", null_df)]
    summary = []
    for name, frame in frames:
        summary.append({"validation": name, "n_cases": len(frame), "mean_rmse": frame.rmse.mean() if not frame.empty else np.nan, "median_rmse": frame.rmse.median() if not frame.empty else np.nan, "mean_mae": frame.mae.mean() if not frame.empty else np.nan, "mean_correlation": frame.correlation.mean() if not frame.empty else np.nan})
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(OUT / "04_validation_summary.csv", index=False)
    pd.DataFrame([{"common_genes": int(matrix.shape[0]), "samples": int(matrix.shape[1]), "datasets": sorted(metadata.dataset.unique().tolist()), "dataset_holdout_cases": len(dataset_df), "replicate_holdout_cases": len(replicate_df), "permutation_cases": len(null_df), "permutations": n_permutations}]).to_json(OUT / "05_stage27_report.json", orient="records", indent=2)
    return summary_df


__all__ = ["stage2_7", "leave_one_dataset_out", "leave_one_replicate_out", "permutation_null"]
