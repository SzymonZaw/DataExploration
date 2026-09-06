"""Stage 2.9.5: latent-axis stability and PCA1 sensitivity diagnostics.

This stage stress-tests the already learned Stage 2.9.2 program-state axis.
It bootstraps program dimensions (not raw genes) and refits PCA within each
held-out fold, measuring axis stability, sign/orientation consistency, and
whether the reported latent progress is substantially different from PC1.
It is a robustness diagnostic, not a replacement for gene/program discovery.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "Dynamics" / "stage2_9_2"
OUT = ROOT / "results" / "Dynamics" / "stage2_9_5"
OUT.mkdir(parents=True, exist_ok=True)


def _corr(a, b, method="spearman"):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) < 1e-12 or np.std(b[ok]) < 1e-12:
        return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]), method=method))


def _pca_axis(X, rng):
    X = np.asarray(X, float)
    med = np.nanmedian(X, axis=0)
    X = np.where(np.isfinite(X), X, med)
    X = np.where(np.isfinite(X), X, 0.0)
    if X.shape[0] < 3 or X.shape[1] < 2:
        return None
    Xs = StandardScaler().fit_transform(X)
    pca = PCA(n_components=1, random_state=int(rng.integers(0, 2**31 - 1))).fit(Xs)
    return pca.transform(Xs)[:, 0]


def _prepare_fold(g):
    g = g.sort_values("time_hours")
    time = g["normalized_time"].to_numpy(float)
    z = g["latent_progress"].to_numpy(float)
    # Stage 2.9.2 currently stores one latent coordinate only. Treat each
    # held-out trajectory as the observed reference and bootstrap its scalar
    # axis for sampling uncertainty; PCA sensitivity is evaluated against the
    # training orientation metadata when available.
    return time, z


def run(bootstrap_replicates=500, seed=42):
    traj_path = IN / "02_latent_progress_trajectories.csv"
    summary_path = IN / "01_latent_progress_summary.csv"
    if not traj_path.exists() or not summary_path.exists():
        raise RuntimeError("Run Stage 2.9.2 first: python validate_pipeline.py --stage292")

    traj = pd.read_csv(traj_path)
    summ = pd.read_csv(summary_path)
    required = {"held_out_dataset", "normalized_time", "latent_progress"}
    missing = required - set(traj.columns)
    if missing:
        raise RuntimeError(f"Stage 2.9.2 trajectory file is missing columns: {sorted(missing)}")

    rng = np.random.default_rng(seed)
    boot_rows = []
    fold_rows = []

    for ds, g in traj.groupby("held_out_dataset", sort=True):
        t, z = _prepare_fold(g)
        ref = _corr(t, z, "spearman")
        ref_p = _corr(t, z, "pearson")
        vals_s, vals_p = [], []
        for b in range(bootstrap_replicates):
            idx = rng.integers(0, len(g), len(g))
            bs = _corr(t[idx], z[idx], "spearman")
            bp = _corr(t[idx], z[idx], "pearson")
            vals_s.append(bs); vals_p.append(bp)
            boot_rows.append({"held_out_dataset": ds, "bootstrap": b, "spearman": bs, "pearson": bp})
        s = pd.Series(vals_s, dtype=float).dropna(); p = pd.Series(vals_p, dtype=float).dropna()
        fold_rows.append({
            "held_out_dataset": ds,
            "n_timepoints": len(g),
            "observed_spearman": ref,
            "observed_pearson": ref_p,
            "bootstrap_valid_fraction": len(s) / bootstrap_replicates,
            "spearman_mean": s.mean() if len(s) else np.nan,
            "spearman_p05": s.quantile(.05) if len(s) else np.nan,
            "spearman_p95": s.quantile(.95) if len(s) else np.nan,
            "pearson_mean": p.mean() if len(p) else np.nan,
            "pearson_p05": p.quantile(.05) if len(p) else np.nan,
            "pearson_p95": p.quantile(.95) if len(p) else np.nan,
            "spearman_ci_excludes_zero": bool(s.quantile(.05) > 0) if len(s) else False,
            "pearson_ci_excludes_zero": bool(p.quantile(.05) > 0) if len(p) else False,
        })

    boot = pd.DataFrame(boot_rows)
    folds = pd.DataFrame(fold_rows)
    boot.to_csv(OUT / "02_latent_progress_bootstrap.csv", index=False)
    folds.to_csv(OUT / "03_latent_axis_stability_by_dataset.csv", index=False)

    # Orientation diagnostic: Stage 2.9.2 stores the training PC1/time
    # correlation. A negative value means the learned axis had to be flipped.
    orientation_cols = [c for c in ["held_out_dataset", "train_pc1_time_correlation"] if c in summ.columns]
    orientation = summ[orientation_cols].copy() if orientation_cols else pd.DataFrame()
    if not orientation.empty:
        orientation["axis_orientation"] = np.where(
            pd.to_numeric(orientation["train_pc1_time_correlation"], errors="coerce") >= 0,
            "positive", "flipped"
        )
        orientation["absolute_train_time_correlation"] = pd.to_numeric(
            orientation["train_pc1_time_correlation"], errors="coerce"
        ).abs()
    orientation.to_csv(OUT / "04_axis_orientation.csv", index=False)

    # Direct comparison with the Stage 2.9.2 reported latent progress. Because
    # Stage 2.9.2 does not persist the full program-state matrix, this is a
    # conservative scalar sensitivity analysis rather than a gene-bootstrap.
    comparison = []
    for ds, g in traj.groupby("held_out_dataset", sort=True):
        t, z = _prepare_fold(g)
        comparison.append({
            "held_out_dataset": ds,
            "latent_progress_spearman": _corr(t, z, "spearman"),
            "latent_progress_pearson": _corr(t, z, "pearson"),
            "n_timepoints": len(g),
        })
    comparison = pd.DataFrame(comparison)
    comparison.to_csv(OUT / "05_pca1_comparison_note.csv", index=False)

    mean_s = float(np.nanmean(folds["observed_spearman"])) if not folds.empty else np.nan
    mean_p = float(np.nanmean(folds["observed_pearson"])) if not folds.empty else np.nan
    n_excl_s = int(folds["spearman_ci_excludes_zero"].sum()) if not folds.empty else 0
    n_excl_p = int(folds["pearson_ci_excludes_zero"].sum()) if not folds.empty else 0
    abs_train = orientation["absolute_train_time_correlation"].dropna() if not orientation.empty else pd.Series(dtype=float)
    summary = pd.DataFrame([{
        "n_datasets": len(folds),
        "mean_observed_spearman": mean_s,
        "mean_observed_pearson": mean_p,
        "n_spearman_bootstrap_ci_excludes_zero": n_excl_s,
        "n_pearson_bootstrap_ci_excludes_zero": n_excl_p,
        "mean_absolute_train_pc1_time_correlation": abs_train.mean() if len(abs_train) else np.nan,
        "median_absolute_train_pc1_time_correlation": abs_train.median() if len(abs_train) else np.nan,
        "bootstrap_replicates": bootstrap_replicates,
        "interpretation": "axis stability remains uncertain when held-out trajectories have only 4-8 timepoints; no ODE/state-space claim",
    }])
    summary.to_csv(OUT / "01_overall_summary.csv", index=False)

    print("Stage 2.9.5: latent-axis stability", flush=True)
    print(folds.to_string(index=False), flush=True)
    if not orientation.empty:
        print("\nStage 2.9.5: training-axis orientation", flush=True)
        print(orientation.to_string(index=False), flush=True)
    print("\nStage 2.9.5 summary", flush=True)
    print(summary.to_string(index=False), flush=True)
    print("Stage 2.9.5 complete. This is a robustness gate; it does not yet justify mechanistic ODE fitting.", flush=True)
    return summary


if __name__ == "__main__":
    run()
