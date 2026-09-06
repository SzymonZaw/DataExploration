"""Stage 2.9.3: stability and null validation of latent reprogramming progress.

This stage does not fit an ODE. It asks whether the Stage 2.9.2 latent progress
coordinate is stable to gene resampling and resistant to a within-held-out
permutation null. Bootstrap refits are leakage-free: genes/programs/PCA are
relearned from training datasets only for every replicate.
"""
from pathlib import Path
import numpy as np
import pandas as pd

from .stage291 import _load_common_space, _trajectories, _profile_cube, _stability, _select, _fit_programs, _state, _state_traj
from .stage292 import _program_grid, _fit_progress_with_scores, _project_progress

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_3"
OUT.mkdir(parents=True, exist_ok=True)


def _fit_bootstrap_progress(matrix, meta, held_out, n_programs, max_genes, seed, gene_fraction=0.8):
    """Fit one leakage-free bootstrap model and return held-out latent progress."""
    datasets = sorted(meta.dataset.astype(str).unique())
    train = [d for d in datasets if d != held_out]
    tr = _trajectories(matrix, meta, train)
    te = _trajectories(matrix, meta, [held_out])
    if held_out not in te or len(tr) < 2:
        return None
    cube, used = _profile_cube(tr)
    stab = _stability(cube)
    ids, _, _ = _select(stab, max_genes)
    if len(ids) < max(2, n_programs):
        return None
    rng = np.random.default_rng(seed)
    n_keep = max(n_programs, int(np.ceil(len(ids) * gene_fraction)))
    ids = rng.choice(ids, size=min(n_keep, len(ids)), replace=False)
    programs = _fit_programs(cube, ids, n_programs, seed)
    state = _state(matrix, meta, programs, train)
    trajs = _state_traj(state, meta)
    train_trajs = {d: trajs[d] for d in train if d in trajs}
    test = trajs.get(held_out)
    model = _fit_progress_with_scores(train_trajs, seed, 17)
    if model is None or test is None or len(test) < 3:
        return None
    pred = _project_progress(test, model)
    actual = (test.index.to_numpy(float) - test.index.min()) / (test.index.max() - test.index.min())
    if pred is None or len(pred) != len(actual):
        return None
    return pd.DataFrame({"time_hours": test.index.to_numpy(float), "normalized_time": actual, "latent_progress": pred})


def _corr(a, b, method="pearson"):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) < 1e-12 or np.std(b[ok]) < 1e-12:
        return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]), method=method))


def _observed():
    """Reuse Stage 2.9.2 observed trajectories, avoiding a second model fit."""
    path = ROOT / "results" / "Dynamics" / "stage2_9_2" / "02_latent_progress_trajectories.csv"
    summary = ROOT / "results" / "Dynamics" / "stage2_9_2" / "01_latent_progress_summary.csv"
    if not path.exists() or not summary.exists():
        raise RuntimeError("Stage 2.9.2 outputs are missing. Run: python validate_pipeline.py --stage292")
    return pd.read_csv(path), pd.read_csv(summary)


def _bootstrap_stability(matrix, meta, held_out, reference, n_bootstrap, n_programs, max_genes, seed):
    vals = []
    ref = reference.sort_values("time_hours").latent_progress.to_numpy(float)
    for b in range(n_bootstrap):
        pred = _fit_bootstrap_progress(matrix, meta, held_out, n_programs, max_genes, seed + b)
        if pred is None:
            continue
        pred = pred.sort_values("time_hours").latent_progress.to_numpy(float)
        if len(pred) == len(ref):
            vals.append({
                "held_out_dataset": held_out,
                "bootstrap": b,
                "progress_pearson_vs_reference": _corr(ref, pred, "pearson"),
                "progress_spearman_vs_reference": _corr(ref, pred, "spearman"),
                "n_timepoints": len(pred),
            })
    return vals


def _permutation_null(observed, permutations, seed):
    """Conditional null: shuffle time labels within each held-out trajectory."""
    rows = []
    rng = np.random.default_rng(seed)
    for _, g in observed.groupby("held_out_dataset", sort=True):
        p = g.latent_progress.to_numpy(float)
        t = g.normalized_time.to_numpy(float)
        if len(p) < 3:
            continue
        for i in range(permutations):
            tp = rng.permutation(t)
            rows.append({
                "held_out_dataset": str(g.held_out_dataset.iloc[0]),
                "permutation": i,
                "spearman": _corr(tp, p, "spearman"),
                "pearson": _corr(tp, p, "pearson"),
            })
    return pd.DataFrame(rows)


def stage2_9_3(n_programs=10, max_genes=2500, bootstrap=200, permutations=1000, seed=42):
    print("Stage 2.9.3: loading Stage 2.9.2 observed latent progress...", flush=True)
    observed, summary292 = _observed()
    matrix, meta = _load_common_space()

    boot_rows = []
    for ds in sorted(observed.held_out_dataset.astype(str).unique()):
        ref = observed[observed.held_out_dataset.astype(str) == ds]
        print(f"Stage 2.9.3: bootstrap stability, held out={ds} ({bootstrap} replicates)...", flush=True)
        vals = _bootstrap_stability(matrix, meta, ds, ref, bootstrap, n_programs, max_genes, seed + 10000)
        boot_rows.extend(vals)

    boot = pd.DataFrame(boot_rows)
    boot.to_csv(OUT / "02_bootstrap_stability.csv", index=False)

    if boot.empty:
        raise RuntimeError("No bootstrap replicates succeeded; inspect Stage 2.9.2 fold status and common-space quality.")

    boot_summary = boot.groupby("held_out_dataset").agg(
        n_bootstrap=("bootstrap", "size"),
        mean_pearson=("progress_pearson_vs_reference", "mean"),
        median_pearson=("progress_pearson_vs_reference", "median"),
        p05_pearson=("progress_pearson_vs_reference", lambda x: x.quantile(.05)),
        p95_pearson=("progress_pearson_vs_reference", lambda x: x.quantile(.95)),
        mean_spearman=("progress_spearman_vs_reference", "mean"),
        median_spearman=("progress_spearman_vs_reference", "median"),
        p05_spearman=("progress_spearman_vs_reference", lambda x: x.quantile(.05)),
        p95_spearman=("progress_spearman_vs_reference", lambda x: x.quantile(.95)),
    ).reset_index()
    boot_summary.to_csv(OUT / "03_bootstrap_stability_summary.csv", index=False)

    print(f"Stage 2.9.3: conditional time-label permutation null ({permutations} per held-out dataset)...", flush=True)
    null = _permutation_null(observed, permutations, seed + 20000)
    null.to_csv(OUT / "04_permutation_null.csv", index=False)

    obs_rows = []
    for ds, g in observed.groupby("held_out_dataset", sort=True):
        obs_rows.append({
            "held_out_dataset": ds,
            "n_timepoints": len(g),
            "observed_spearman": _corr(g.normalized_time, g.latent_progress, "spearman"),
            "observed_pearson": _corr(g.normalized_time, g.latent_progress, "pearson"),
        })
    obs = pd.DataFrame(obs_rows)
    records = []
    for _, r in obs.iterrows():
        ds = r.held_out_dataset
        n = null[null.held_out_dataset == ds]
        ps = (1 + int((n.spearman >= r.observed_spearman).sum())) / (len(n) + 1) if len(n) else np.nan
        pp = (1 + int((n.pearson >= r.observed_pearson).sum())) / (len(n) + 1) if len(n) else np.nan
        b = boot_summary[boot_summary.held_out_dataset == ds]
        records.append({**r.to_dict(), "permutation_p_spearman": ps, "permutation_p_pearson": pp,
                        "bootstrap_mean_spearman": float(b.mean_spearman.iloc[0]) if len(b) else np.nan,
                        "bootstrap_p05_spearman": float(b.p05_spearman.iloc[0]) if len(b) else np.nan,
                        "bootstrap_p95_spearman": float(b.p95_spearman.iloc[0]) if len(b) else np.nan})
    final = pd.DataFrame(records)
    final.to_csv(OUT / "01_fold_stability_summary.csv", index=False)

    mean_rho = float(final.observed_spearman.mean())
    mean_pear = float(final.observed_pearson.mean())
    mean_p_rho = float(final.permutation_p_spearman.mean())
    mean_p_pear = float(final.permutation_p_pearson.mean())
    print("\nStage 2.9.3 fold stability:")
    print(final.to_string(index=False))
    print(f"Observed mean Spearman: {mean_rho:.3f}; Pearson: {mean_pear:.3f}")
    print(f"Mean conditional permutation p: Spearman={mean_p_rho:.4f}; Pearson={mean_p_pear:.4f}")
    print("Bootstrap tests latent-axis stability; permutation tests temporal ordering conditional on the learned axis.")
    print("Do not fit ODE/state-space yet unless latent progress is stable across bootstrap resampling and resistant to the permutation null.")
    return final


if __name__ == "__main__":
    stage2_9_3()
