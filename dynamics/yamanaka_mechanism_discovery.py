"""Stage 2.11: perturbation-conditioned, biologically interpretable dynamics.

This is a deliberately conservative research benchmark. It uses the existing
trajectory POC outputs as the starting point and asks whether pathway/TF
activity plus explicit perturbation information improves cross-dataset
prediction of the next observed state.

No causal or lineage claim is made here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np
import pandas as pd


DATASETS = ("GM00731", "HFIB_COMBINED")


def _read_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError(f"Empty matrix: {path}")
    return df


def _find_column(df: pd.DataFrame, tokens: Iterable[str]) -> str | None:
    for col in df.columns:
        name = str(col).lower()
        if any(token in name for token in tokens):
            return col
    return None


def _standardize_activities(df: pd.DataFrame) -> pd.DataFrame:
    """Column-wise z-score while retaining constant columns as zero."""
    values = df.astype(float)
    mean = values.mean(axis=0)
    std = values.std(axis=0, ddof=0).replace(0, 1.0)
    return (values - mean) / std


def _load_inputs(root: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load Stage 2.11 trajectory/representation artifacts if present.

    The function accepts a small number of conventional filenames so the
    benchmark can be run against locally generated Phase 1 outputs without
    coupling the scientific logic to one artifact layout.
    """
    candidates = {
        "activities": [
            root / "results" / "representation_activities.csv",
            root / "results" / "phase1_activities.csv",
            root / "results" / "yamanaka_activities.csv",
        ],
        "metadata": [
            root / "results" / "representation_metadata.csv",
            root / "results" / "phase1_metadata.csv",
            root / "results" / "yamanaka_metadata.csv",
        ],
    }
    found = {}
    for kind, paths in candidates.items():
        for path in paths:
            if path.exists():
                found[kind] = path
                break
    if set(found) != {"activities", "metadata"}:
        raise FileNotFoundError(
            "Stage 2.11 input artifacts not found. Expected an activities CSV "
            "and a metadata CSV under results/. Run the representation/trajectory "
            "pipeline first."
        )
    return _read_matrix(found["activities"]), _read_matrix(found["metadata"])


def _align(activities: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    common = activities.index.intersection(metadata.index)
    if len(common) < 4:
        raise ValueError(f"Only {len(common)} samples/groups can be aligned")
    out = metadata.loc[common].copy()
    acts = activities.loc[common].copy()
    out = out.join(acts.add_prefix("ACT_"))
    return out


def _parse_day(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.lower()
    result = pd.Series(np.nan, index=series.index, dtype=float)
    for day in (0, 3, 7, 10):
        mask = text.str.contains(rf"(?:^|[^0-9]){day}(?:$|[^0-9])", regex=True)
        result.loc[mask] = float(day)
    return result


def _prepare(df: pd.DataFrame) -> pd.DataFrame:
    day_col = _find_column(df, ("day", "time", "hour"))
    if day_col is None:
        raise ValueError("No day/time metadata column found")
    days = _parse_day(df[day_col])
    df = df.copy()
    df["day"] = days
    df = df[df["day"].notna()].copy()
    dataset_col = _find_column(df, ("dataset", "study"))
    if dataset_col is None:
        df["dataset"] = "UNKNOWN"
    else:
        df["dataset"] = df[dataset_col].astype(str)
    perturb_col = _find_column(df, ("perturb", "condition", "treatment", "factor", "group"))
    df["perturbation"] = df[perturb_col].astype(str) if perturb_col else "UNKNOWN"
    return df


def _fit_ridge(X: np.ndarray, Y: np.ndarray, alpha: float = 1.0) -> np.ndarray:
    X1 = np.column_stack([np.ones(len(X)), X])
    reg = np.eye(X1.shape[1])
    reg[0, 0] = 0.0
    return np.linalg.solve(X1.T @ X1 + alpha * reg, X1.T @ Y)


def _predict(beta: np.ndarray, X: np.ndarray) -> np.ndarray:
    return np.column_stack([np.ones(len(X)), X]) @ beta


def _rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def _build_transitions(df: pd.DataFrame, activity_cols: list[str]) -> list[dict]:
    rows = []
    for dataset, sub in df.groupby("dataset"):
        sub = sub.sort_values("day")
        for _, row in sub.iterrows():
            future = sub[sub["day"] > row["day"]]
            if future.empty:
                continue
            nxt = future.iloc[0]
            rows.append({"dataset": dataset, "day": row["day"], "next_day": nxt["day"], "row": row, "next": nxt})
    return rows


def run(root: Path) -> Dict:
    activities, metadata = _load_inputs(root)
    df = _align(activities, metadata)
    df = _prepare(df)

    activity_cols = [c for c in df.columns if c.startswith("ACT_")]
    if not activity_cols:
        raise ValueError("No activity features found")
    df[activity_cols] = _standardize_activities(df[activity_cols])

    transitions = _build_transitions(df, activity_cols)
    if len(transitions) < 2:
        raise ValueError("Not enough temporal transitions for a benchmark")

    # Deterministic leave-one-dataset-out benchmark.
    results = []
    for held_out in sorted({t["dataset"] for t in transitions}):
        train = [t for t in transitions if t["dataset"] != held_out]
        test = [t for t in transitions if t["dataset"] == held_out]
        if not train or not test:
            continue

        def matrix(ts, include_perturbation: bool):
            X, Y = [], []
            perturb_levels = sorted({str(t["row"]["perturbation"]) for t in train})
            for t in ts:
                x = t["row"][activity_cols].to_numpy(dtype=float)
                if include_perturbation:
                    p = np.array([float(str(t["row"]["perturbation"]) == level) for level in perturb_levels])
                    x = np.concatenate([x, p])
                X.append(x)
                Y.append(t["next"][activity_cols].to_numpy(dtype=float))
            return np.asarray(X), np.asarray(Y)

        Xs, Ys = matrix(train, False)
        Xp, Yp = matrix(train, True)
        Xt, Yt = matrix(test, False)
        # For held-out perturbation categories, one-hot columns naturally become zero.
        Xtt, _ = matrix(test, True)

        beta_state = _fit_ridge(Xs, Ys)
        beta_pert = _fit_ridge(Xp, Yp)
        pred_state = _predict(beta_state, Xt)
        pred_pert = _predict(beta_pert, Xtt)

        persistence = np.asarray([t["row"][activity_cols].to_numpy(float) for t in test])
        nearest = persistence.copy()
        # nearest-time baseline is identical to persistence here for strictly
        # ordered sample-level transitions; retained explicitly for auditability.

        results.append({
            "held_out_dataset": held_out,
            "n_test_transitions": len(test),
            "state_only_rmse": _rmse(pred_state, Yt),
            "state_plus_perturbation_rmse": _rmse(pred_pert, Yt),
            "persistence_rmse": _rmse(persistence, Yt),
            "nearest_rmse": _rmse(nearest, Yt),
        })

    return {
        "status": "ok",
        "n_features": len(activity_cols),
        "n_aligned_groups": int(len(df)),
        "n_transitions": len(transitions),
        "folds": results,
        "guardrail": "Predictive improvement is evidence of useful state/perturbation information, not causal mechanism discovery.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    print(json.dumps(run(args.root), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
