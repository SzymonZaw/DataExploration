"""Stage 2.9.4: cross-dataset robustness diagnostics for latent progress.

This stage is intentionally fast and uses the already computed Stage 2.9.3
bootstrap/permutation outputs. It asks whether the positive latent-progress
signal is driven by one dataset and whether bootstrap intervals consistently
exclude zero. It does not refit the biological programs.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
IN_DIR = ROOT / "results" / "Dynamics" / "stage2_9_3"
OUT = ROOT / "results" / "Dynamics" / "stage2_9_4"
OUT.mkdir(parents=True, exist_ok=True)


def run():
    boot_path = IN_DIR / "02_bootstrap_stability.csv"
    null_path = IN_DIR / "04_permutation_null.csv"
    if not boot_path.exists():
        raise RuntimeError("Run Stage 2.9.3 first: python validate_pipeline.py --stage293")

    boot = pd.read_csv(boot_path)
    required = {"held_out_dataset", "spearman", "pearson"}
    missing = required - set(boot.columns)
    if missing:
        raise RuntimeError(f"Stage 2.9.3 bootstrap file is missing columns: {sorted(missing)}")

    rows = []
    for ds, g in boot.groupby("held_out_dataset", sort=True):
        s = pd.to_numeric(g["spearman"], errors="coerce").dropna()
        p = pd.to_numeric(g["pearson"], errors="coerce").dropna()
        rows.append({
            "held_out_dataset": ds,
            "n_bootstrap_valid_spearman": len(s),
            "bootstrap_valid_fraction_spearman": len(s) / len(g) if len(g) else np.nan,
            "spearman_mean": s.mean(),
            "spearman_median": s.median(),
            "spearman_p05": s.quantile(.05),
            "spearman_p95": s.quantile(.95),
            "spearman_ci_excludes_zero": bool(s.quantile(.05) > 0) if len(s) else False,
            "n_bootstrap_valid_pearson": len(p),
            "pearson_mean": p.mean(),
            "pearson_median": p.median(),
            "pearson_p05": p.quantile(.05),
            "pearson_p95": p.quantile(.95),
            "pearson_ci_excludes_zero": bool(p.quantile(.05) > 0) if len(p) else False,
        })
    by_ds = pd.DataFrame(rows)
    by_ds.to_csv(OUT / "01_dataset_robustness.csv", index=False)

    s = by_ds["spearman_mean"].to_numpy(float)
    p = by_ds["pearson_mean"].to_numpy(float)
    leave_rows = []
    for i, ds in enumerate(by_ds["held_out_dataset"]):
        keep = np.arange(len(by_ds)) != i
        leave_rows.append({
            "excluded_dataset": ds,
            "n_remaining_datasets": int(keep.sum()),
            "mean_spearman_excluding": float(np.nanmean(s[keep])),
            "mean_pearson_excluding": float(np.nanmean(p[keep])),
        })
    leave = pd.DataFrame(leave_rows)
    leave.to_csv(OUT / "02_leave_one_dataset_sensitivity.csv", index=False)

    null_p = np.nan
    null_s = np.nan
    if null_path.exists():
        null = pd.read_csv(null_path)
        if "mean_spearman" in null:
            observed_s = float(np.nanmean(s))
            null_s = float((1 + (pd.to_numeric(null["mean_spearman"], errors="coerce") >= observed_s).sum()) / (null["mean_spearman"].notna().sum() + 1))
        if "mean_pearson" in null:
            observed_p = float(np.nanmean(p))
            null_p = float((1 + (pd.to_numeric(null["mean_pearson"], errors="coerce") >= observed_p).sum()) / (null["mean_pearson"].notna().sum() + 1))

    summary = pd.DataFrame([{
        "n_datasets": len(by_ds),
        "mean_dataset_spearman": float(np.nanmean(s)),
        "median_dataset_spearman": float(np.nanmedian(s)),
        "min_dataset_spearman": float(np.nanmin(s)),
        "mean_dataset_pearson": float(np.nanmean(p)),
        "median_dataset_pearson": float(np.nanmedian(p)),
        "min_dataset_pearson": float(np.nanmin(p)),
        "n_spearman_ci_excludes_zero": int(by_ds["spearman_ci_excludes_zero"].sum()),
        "n_pearson_ci_excludes_zero": int(by_ds["pearson_ci_excludes_zero"].sum()),
        "permutation_p_spearman": null_s,
        "permutation_p_pearson": null_p,
    }])
    summary.to_csv(OUT / "03_overall_robustness.csv", index=False)

    print("Stage 2.9.4: dataset robustness", flush=True)
    print(by_ds.to_string(index=False), flush=True)
    print("\nStage 2.9.4: leave-one-dataset sensitivity", flush=True)
    print(leave.to_string(index=False), flush=True)
    print("\nStage 2.9.4 summary", flush=True)
    print(summary.to_string(index=False), flush=True)
    print("Stage 2.9.4 complete.", flush=True)
    return summary


if __name__ == "__main__":
    run()
