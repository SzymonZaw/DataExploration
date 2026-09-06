"""Stage 2.9.18: conservative lead/lag validation with null controls.

Tests whether apparent non-zero program lags exceed a within-dataset
permutation null. This is a diagnostic, not causal inference and not an ODE.
Because trajectories contain only 4-8 time points, results are deliberately
conservative: only small fixed lags on a 9-point interpolation grid are used,
and a pair must show the same non-zero direction in at least two datasets.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "results" / "Dynamics" / "stage2_9_14"
OUT = ROOT / "results" / "Dynamics" / "stage2_9_18"
OUT.mkdir(parents=True, exist_ok=True)

def log(x): print(f"Stage 2.9.18: {x}", flush=True)

def _load():
    p = SRC / "03_program_trajectories.csv"
    if not p.exists(): raise RuntimeError("Run Stage 2.9.14 first")
    x = pd.read_csv(p)
    x["dataset"] = x["dataset"].astype(str)
    x["program_id"] = x["program_id"].astype(str)
    return x

def _curve(g, grid):
    g = g.sort_values("time_hours")
    t = g["time_hours"].to_numpy(float)
    y = g["activity"].to_numpy(float)
    if len(np.unique(t)) < 4 or np.ptp(t) <= 0: return None
    tn = (t - t.min()) / np.ptp(t)
    return np.interp(grid, tn, y)

def _residual(y, degree=2):
    x = np.linspace(0.0, 1.0, len(y))
    degree = min(degree, len(y) - 1)
    X = np.column_stack([x ** k for k in range(degree + 1)])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    return y - X @ beta

def _corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    if len(a) < 4 or np.std(a) < 1e-12 or np.std(b) < 1e-12: return np.nan
    return float(np.corrcoef(a, b)[0, 1])

def _lag_corr(a, b, lag):
    if lag < 0: x, y = a[:lag], b[-lag:]
    elif lag > 0: x, y = a[lag:], b[:-lag]
    else: x, y = a, b
    return _corr(x, y)

def _best(a, b, lags):
    vals = [(lag, _lag_corr(a, b, lag)) for lag in lags]
    vals = [(l, r) for l, r in vals if np.isfinite(r)]
    if not vals: return np.nan, np.nan
    return max(vals, key=lambda z: abs(z[1]))

def _pairs(df):
    grid = np.linspace(0, 1, 9)
    curves = {}
    for (ds, pid), g in df.groupby(["dataset", "program_id"]):
        y = _curve(g, grid)
        if y is not None: curves[(ds, pid)] = _residual(y)
    pids = sorted(df.program_id.unique()); rows = []
    # Fixed non-zero lags only; +/- 1 and +/- 2 grid steps are less extreme
    # than the +/-3 scan used in Stage 2.9.17.
    lags = [-2, -1, 0, 1, 2]
    for ds in sorted(df.dataset.unique()):
        for i, a in enumerate(pids):
            for b in pids[i + 1:]:
                if (ds, a) not in curves or (ds, b) not in curves: continue
                ya, yb = curves[(ds, a)], curves[(ds, b)]
                best_lag, best_r = _best(ya, yb, lags)
                zero = _lag_corr(ya, yb, 0)
                nonzero = [(l, _lag_corr(ya, yb, l)) for l in lags if l != 0]
                nonzero = [(l, r) for l, r in nonzero if np.isfinite(r)]
                nz_lag, nz_r = (max(nonzero, key=lambda z: abs(z[1])) if nonzero else (np.nan, np.nan))
                rows.append({"dataset": ds, "program_a": a, "program_b": b,
                             "best_lag": best_lag, "best_abs_pearson": abs(best_r) if np.isfinite(best_r) else np.nan,
                             "zero_lag_pearson": zero,
                             "best_nonzero_lag": nz_lag,
                             "best_nonzero_abs_pearson": abs(nz_r) if np.isfinite(nz_r) else np.nan,
                             "nonzero_minus_zero_abs": (abs(nz_r) - abs(zero)) if np.isfinite(nz_r) and np.isfinite(zero) else np.nan})
    return pd.DataFrame(rows)

def _null(df, observed, n_perm=1000, seed=2918):
    rng = np.random.default_rng(seed); grid = np.linspace(0, 1, 9); lags = [-2, -1, 0, 1, 2]
    curves = {}
    for (ds, pid), g in df.groupby(["dataset", "program_id"]):
        y = _curve(g, grid)
        if y is not None: curves[(ds, pid)] = _residual(y)
    rows = []
    for idx, r in observed.iterrows():
        a = curves.get((r.dataset, r.program_a)); b = curves.get((r.dataset, r.program_b))
        if a is None or b is None: continue
        vals = []
        for _ in range(n_perm):
            bp = rng.permutation(b)
            nz = [abs(_lag_corr(a, bp, l)) for l in lags if l != 0]
            nz = [v for v in nz if np.isfinite(v)]
            z = abs(_lag_corr(a, bp, 0))
            if nz and np.isfinite(z): vals.append(max(nz) - z)
        obs_delta = r.nonzero_minus_zero_abs
        p = (1 + np.sum(np.asarray(vals) >= obs_delta)) / (1 + len(vals)) if vals and np.isfinite(obs_delta) else np.nan
        rows.append({"dataset": r.dataset, "program_a": r.program_a, "program_b": r.program_b,
                     "observed_nonzero_minus_zero_abs": obs_delta,
                     "permutation_n": len(vals), "empirical_p": p,
                     "null_mean_delta": float(np.mean(vals)) if vals else np.nan,
                     "null_p95_delta": float(np.quantile(vals, .95)) if vals else np.nan})
    return pd.DataFrame(rows)

def run(permutations=1000):
    log("starting conservative lag/null analysis")
    df = _load(); obs = _pairs(df)
    obs.to_csv(OUT / "01_within_dataset_lag_effects.csv", index=False)
    null = _null(df, obs, permutations); null.to_csv(OUT / "02_lag_permutation_null.csv", index=False)
    if len(null):
        key = ["dataset", "program_a", "program_b"]
        merged = obs.merge(null, on=key, how="left")
        merged.to_csv(OUT / "03_lag_effects_with_pvalues.csv", index=False)
        pair = merged.groupby(["program_a", "program_b"], as_index=False).agg(
            n_datasets=("dataset", "nunique"),
            mean_nonzero_minus_zero_abs=("nonzero_minus_zero_abs", "mean"),
            n_significant=("empirical_p", lambda x: int(np.sum(np.asarray(x) < .05))),
            n_nonzero=("best_nonzero_lag", lambda x: int(np.sum(pd.to_numeric(x, errors="coerce").notna()))),
            median_best_nonzero_lag=("best_nonzero_lag", "median"))
        pair["directional_consistency"] = pair.apply(lambda r: np.nan, axis=1)
        # Direction is evaluated only among dataset-level effects that are
        # individually positive and significant; with 3 datasets this avoids
        # treating a single lucky lag as reproducible.
        sig = merged[(merged["empirical_p"] < .05) & (merged["nonzero_minus_zero_abs"] > 0)]
        if len(sig):
            d = sig.groupby(["program_a", "program_b"])["best_nonzero_lag"].apply(lambda x: float(max(np.sum(np.asarray(x)>0), np.sum(np.asarray(x)<0)) / len(x))).rename("directional_consistency")
            pair = pair.drop(columns=["directional_consistency"]).merge(d, on=["program_a", "program_b"], how="left")
        pair.to_csv(OUT / "04_cross_dataset_lag_summary.csv", index=False)
    else:
        merged = pd.DataFrame(); pair = pd.DataFrame()
    supported = bool(len(pair) and np.any((pair["n_significant"] >= 2) & (pair["directional_consistency"] >= .67)))
    overall = pd.DataFrame([{
        "n_datasets": int(df.dataset.nunique()), "n_program_pairs": int(len(pair)),
        "n_pairs_with_ge2_significant_datasets": int(np.sum(pair.n_significant >= 2)) if len(pair) else 0,
        "n_pairs_reproducible_direction": int(np.sum((pair.n_significant >= 2) & (pair.directional_consistency >= .67))) if len(pair) else 0,
        "mean_nonzero_vs_zero_abs_effect": float(merged.nonzero_minus_zero_abs.mean()) if len(merged) else np.nan,
        "lag_signal_supported": supported,
        "stage3_readiness": False
    }])
    overall.to_csv(OUT / "05_STAGE2_9_18_SUMMARY.csv", index=False)
    log("complete")
    print("\nStage 2.9.18 cross-dataset lag summary:")
    print(pair.sort_values(["n_significant", "mean_nonzero_minus_zero_abs"], ascending=False).head(20).to_string(index=False) if len(pair) else "none")
    print("\nStage 2.9.18 gate:"); print(overall.to_string(index=False))
    return overall

if __name__ == "__main__": run()
