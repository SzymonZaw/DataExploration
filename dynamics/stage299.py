"""Stage 2.9.9: biologically anchored program-state construction.

This stage converts the enrichment evidence from Stage 2.9.8 into a small,
interpretable set of biological programs and computes their sample-level
activities. It also runs an *exploratory* leave-one-dataset-out projection
using fixed program definitions from Stage 2.9.8.

Important: because the Stage 2.9.8 enrichment was computed from the global
recurrence set, this LODO is not a fully leakage-free validation of program
discovery. It is a projection/sensitivity diagnostic. A leakage-free version
must repeat enrichment independently inside every training fold (planned for
Stage 2.9.10).
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
COMMON = ROOT / "results" / "Dynamics" / "stage2_6"
IN = ROOT / "results" / "Dynamics" / "stage2_9_8"
OUT = ROOT / "results" / "Dynamics" / "stage2_9_9"
OUT.mkdir(parents=True, exist_ok=True)


def _log(msg):
    print(f"Stage 2.9.9: {msg}", flush=True)


def _normalise_gene(x):
    return str(x).strip().upper()


def _safe_corr(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _time_column(metadata):
    if "time_hours" in metadata.columns:
        return pd.to_numeric(metadata["time_hours"], errors="coerce")
    from dynamics.validation import _time_hours_for_validation, _strip_dataset_prefix
    vals = []
    counters = {}
    for _, row in metadata.iterrows():
        ds = str(row["dataset"])
        idx = counters.get(ds, 0)
        counters[ds] = idx + 1
        vals.append(_time_hours_for_validation(
            ds, _strip_dataset_prefix(row["sample"]), idx if ds == "GSE28688" else None
        ))
    return pd.Series(vals, index=metadata.index, dtype=float)


def _load_common_space():
    matrix_path = COMMON / "06_common_human_gene_matrix.csv"
    metadata_path = COMMON / "07_common_gene_sample_metadata.csv"
    if not matrix_path.exists() or not metadata_path.exists():
        raise FileNotFoundError("Stage 2.6 common-space files are missing.")
    matrix = pd.read_csv(matrix_path, index_col=0).apply(pd.to_numeric, errors="coerce")
    metadata = pd.read_csv(metadata_path)
    metadata["dataset"] = metadata["dataset"].astype(str)
    metadata["sample"] = metadata["sample"].astype(str)
    metadata["time_hours"] = _time_column(metadata)
    from dynamics.validation import _build_matrix_column_map, _recover_ordered_sample_labels
    metadata = _recover_ordered_sample_labels(metadata)
    metadata["time_hours"] = _time_column(metadata)
    metadata = _build_matrix_column_map(matrix, metadata)
    metadata = metadata[metadata["matrix_column"].notna()].copy()
    return matrix, metadata


def _load_enrichment():
    path = IN / "07_program_summary.csv"
    if not path.exists():
        raise FileNotFoundError("Stage 2.9.8 program summary is missing; run --stage298 first.")
    x = pd.read_csv(path)
    x["term_name"] = x["term_name"].fillna("").astype(str)
    x["p_value"] = pd.to_numeric(x["p_value"], errors="coerce")
    x["intersection_size"] = pd.to_numeric(x["intersection_size"], errors="coerce")
    x["intersection_genes"] = x["intersection_genes"].fillna("").astype(str)
    return x


def _is_generic(term):
    s = term.lower()
    bad = (
        "hiv", "viral messenger", "metabolism of rna", "gene expression",
        "protein metabolic process", "rna polymerase ii transcription",
        "processing of capped intron-containing pre-mrna",
    )
    return any(k in s for k in bad)


def _program_label(term):
    s = term.lower()
    if "fgfr" in s or "pi-3k" in s or "pi3k" in s:
        return "FGFR_PI3K_SIGNALING"
    if any(k in s for k in ("pou5f1", "oct4", "sox2", "nanog")):
        return "OCT4_SOX2_NANOG"
    if "alternative splicing" in s:
        return "FGFR_SPLICING"
    if "nucleus" in s or "localization" in s:
        return "NUCLEAR_LOCALIZATION"
    if "rna polymerase" in s or "transcription" in s:
        return "TRANSCRIPTION"
    return "BIOLOGICAL_PROGRAM"


def _gene_set(text):
    return {_normalise_gene(x) for x in str(text).split(",") if str(x).strip()}


def _jaccard(a, b):
    u = a | b
    return len(a & b) / len(u) if u else 0.0


def _select_programs(enrichment, background):
    x = enrichment.copy()
    x = x[(x["p_value"] < 0.05) & (x["intersection_size"] >= 5)].copy()
    x = x[~x["term_name"].map(_is_generic)].copy()
    x["genes"] = x["intersection_genes"].map(_gene_set)
    bg = {_normalise_gene(g) for g in background}
    x["genes"] = x["genes"].map(lambda s: s & bg)
    x = x[x["genes"].map(len) >= 5].copy()
    x["program_family"] = x["term_name"].map(_program_label)
    x = x.sort_values(["p_value", "intersection_size"], ascending=[True, False])

    selected = []
    family_seen = set()
    for _, row in x.iterrows():
        genes = row["genes"]
        family = row["program_family"]
        # Prefer the most significant term per biological family, while also
        # allowing a second term when it adds substantially different genes.
        redundant = any(_jaccard(genes, q["genes"]) >= 0.75 for q in selected)
        if family in family_seen and redundant:
            continue
        if redundant:
            continue
        selected.append({
            "program_id": f"P{len(selected)+1:02d}",
            "program_family": family,
            "term_id": row["term_id"],
            "term_name": row["term_name"],
            "p_value": row["p_value"],
            "n_genes": len(genes),
            "genes": genes,
        })
        family_seen.add(family)
        if len(selected) >= 8:
            break

    # If the strict non-redundancy filter leaves too few programs, add the
    # strongest remaining terms at a lower overlap threshold.
    if len(selected) < 4:
        for _, row in x.iterrows():
            genes = row["genes"]
            if any(_jaccard(genes, q["genes"]) >= 0.90 for q in selected):
                continue
            selected.append({
                "program_id": f"P{len(selected)+1:02d}",
                "program_family": row["program_family"],
                "term_id": row["term_id"],
                "term_name": row["term_name"],
                "p_value": row["p_value"],
                "n_genes": len(genes),
                "genes": genes,
            })
            if len(selected) >= 4:
                break

    return selected


def _activity_matrix(matrix, programs):
    genes = {_normalise_gene(g): g for g in matrix.index}
    rows = []
    for p in programs:
        matched = [genes[g] for g in p["genes"] if g in genes]
        if not matched:
            continue
        values = matrix.loc[matched].to_numpy(dtype=float)
        # Gene-wise standardisation prevents large-expression genes from
        # dominating a program. NaNs are ignored per sample.
        mean = np.nanmean(values, axis=1, keepdims=True)
        sd = np.nanstd(values, axis=1, keepdims=True)
        z = (values - mean) / np.where(sd > 1e-12, sd, np.nan)
        activity = np.nanmean(z, axis=0)
        rows.append(pd.Series(activity, index=matrix.columns, name=p["program_id"]))
        p["n_matrix_genes"] = len(matched)
    if not rows:
        return pd.DataFrame(index=matrix.columns)
    return pd.DataFrame(rows)


def _orient_programs(activity, metadata):
    meta = metadata.set_index("matrix_column")
    times = meta.loc[activity.columns, "time_hours"].to_numpy(float)
    for program in activity.index:
        a = activity.loc[program].to_numpy(float)
        corr = _safe_corr(a, times)
        if np.isfinite(corr) and corr < 0:
            activity.loc[program] = -a
    return activity


def _trajectory_table(activity, metadata):
    m = metadata.set_index("matrix_column")
    records = []
    for ds, g in metadata.groupby("dataset", sort=True):
        g = g[g["time_hours"].notna()].copy()
        for t, gt in g.groupby("time_hours", sort=True):
            cols = [c for c in gt["matrix_column"] if c in activity.columns]
            if not cols:
                continue
            for program in activity.index:
                vals = pd.to_numeric(activity.loc[program, cols], errors="coerce")
                records.append({
                    "dataset": ds,
                    "time_hours": float(t),
                    "program_id": program,
                    "activity": float(vals.mean()) if vals.notna().any() else np.nan,
                    "n_samples": int(vals.notna().sum()),
                })
    return pd.DataFrame(records)


def _interpolate(times, values, target):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    ok = np.isfinite(times) & np.isfinite(values)
    if ok.sum() < 2 or target < times[ok].min() or target > times[ok].max():
        return np.nan
    order = np.argsort(times[ok])
    return float(np.interp(target, times[ok][order], values[ok][order]))


def _exploratory_lodo(activity, metadata):
    rows = []
    trajectories = {}
    for ds, g in metadata.groupby("dataset", sort=True):
        g = g[g["time_hours"].notna()]
        if g["time_hours"].nunique() < 2:
            continue
        records = []
        for t, gt in g.groupby("time_hours", sort=True):
            cols = [c for c in gt["matrix_column"] if c in activity.columns]
            if cols:
                records.append((float(t), activity[cols].mean(axis=1).to_numpy(float)))
        if len(records) >= 2:
            trajectories[ds] = (np.array([r[0] for r in records]), np.vstack([r[1] for r in records]))

    for held, (test_times, test_values) in sorted(trajectories.items()):
        train = {k: v for k, v in trajectories.items() if k != held}
        for j, target in enumerate(test_times):
            preds = []
            for times, values in train.values():
                if target < times.min() or target > times.max():
                    continue
                pred = np.array([_interpolate(times, values[:, p], target) for p in range(values.shape[1])])
                if np.isfinite(pred).all():
                    preds.append(pred)
            if not preds:
                continue
            pred = np.mean(preds, axis=0)
            true = test_values[j]
            rows.append({
                "held_out_dataset": held,
                "time_hours": float(target),
                "n_training_datasets": len(preds),
                "rmse": float(np.sqrt(np.mean((true - pred) ** 2))),
                "mae": float(np.mean(np.abs(true - pred))),
                "program_profile_correlation": _safe_corr(true, pred),
            })
    return pd.DataFrame(rows)


def run():
    _log("starting biological program construction")
    matrix, metadata = _load_common_space()
    enrichment = _load_enrichment()
    background = matrix.index.astype(str).tolist()
    _log(f"loaded {matrix.shape[0]:,} genes x {matrix.shape[1]} samples")

    programs = _select_programs(enrichment, background)
    _log(f"selected {len(programs)} non-redundant biologically anchored programs")
    if not programs:
        raise RuntimeError("No biological programs could be constructed from Stage 2.9.8 enrichment.")

    activity = _activity_matrix(matrix, programs)
    activity = _orient_programs(activity, metadata)
    _log(f"computed activity for {len(activity)} programs")

    program_rows = []
    for p in programs:
        program_rows.append({
            "program_id": p["program_id"],
            "program_family": p["program_family"],
            "term_id": p["term_id"],
            "term_name": p["term_name"],
            "p_value": p["p_value"],
            "n_enrichment_genes": p["n_genes"],
            "n_matrix_genes": p.get("n_matrix_genes", 0),
            "genes": ",".join(sorted(p["genes"])),
        })
    pd.DataFrame(program_rows).to_csv(OUT / "01_program_definitions.csv", index=False)

    activity_out = activity.copy()
    activity_out.insert(0, "program_id", activity_out.index)
    activity_out.to_csv(OUT / "02_program_activity_by_sample.csv", index=False)

    trajectory = _trajectory_table(activity, metadata)
    trajectory.to_csv(OUT / "03_program_trajectories.csv", index=False)

    lodo = _exploratory_lodo(activity, metadata)
    lodo.to_csv(OUT / "04_exploratory_lodo.csv", index=False)

    if len(lodo):
        summary = lodo.groupby("held_out_dataset").agg(
            n_timepoints=("time_hours", "count"),
            mean_rmse=("rmse", "mean"),
            mean_mae=("mae", "mean"),
            mean_profile_correlation=("program_profile_correlation", "mean"),
        ).reset_index()
    else:
        summary = pd.DataFrame(columns=["held_out_dataset", "n_timepoints", "mean_rmse", "mean_mae", "mean_profile_correlation"])
    summary["validation_type"] = "exploratory_fixed_program_lodo"
    summary.to_csv(OUT / "05_lodo_summary.csv", index=False)

    _log("complete; no ODE/state-space model fitted")
    print("\nStage 2.9.9 programs:", flush=True)
    print(pd.DataFrame(program_rows)[[
        "program_id", "program_family", "term_name", "p_value", "n_matrix_genes"
    ]].to_string(index=False), flush=True)
    print("\nStage 2.9.9 exploratory LODO:", flush=True)
    print(summary.to_string(index=False), flush=True)
    return summary


if __name__ == "__main__":
    run()
