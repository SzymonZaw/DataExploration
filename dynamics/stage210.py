"""Stage 2.10: leakage-free conserved biological transition modules.

Goal: test whether small, biologically defined programs preserve their temporal
ordering across heterogeneous reprogramming datasets even when a global state
coordinate is not transferable.  This stage is deliberately descriptive:
there is no ODE/state-space fitting and Stage 3 remains closed.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_10"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
GRID = np.linspace(0.0, 1.0, 40)
N_BOOT = 500
N_PERM = 1000

PROGRAMS = {
    "P01_PLURIPOTENCY": ["POU5F1","SOX2","NANOG","LIN28A","LIN28B","DPPA4","UTF1","ESRRB","KLF4","MYC"],
    "P02_PROLIFERATION": ["MKI67","PCNA","TOP2A","CCNB1","CCNB2","CCNE1","CDK1","CDC20","UBE2C","TYMS","MCM2","MCM3","MCM5","MCM6","MCM7"],
    "P03_EMT_MESENCHYMAL": ["VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","FN1","ITGA5","COL1A1","COL1A2","CDH2","CDH1","EPCAM","KRT8","KRT18","KRT19"],
    "P04_STRESS_RESPONSE": ["DDIT3","ATF4","HSPA1A","HSPA1B","HMOX1","XBP1","JUN","FOS","DUSP1","PPP1R15A","DNAJB1"],
    "P05_GLYCOLYTIC_METABOLISM": ["SLC2A1","HK2","PFKP","ALDOA","GAPDH","ENO1","PKM","LDHA","PGK1","TPI1","PDK1"],
    "P06_FGFR_PI3K_MAPK": ["FGFR1","FGFR2","FGFR3","FGFR4","FRS2","PLCG1","PIK3CA","PIK3CB","AKT1","AKT2","MAPK1","MAPK3","RAF1","SOS1"],
    "P07_CHROMATIN_EPIGENETIC": ["KMT2A","KMT2B","EZH2","DNMT1","DNMT3A","DNMT3B","TET1","TET2","HDAC1","HDAC2","SMARCA4","ARID1A","CHD4","SUZ12"],
    "P08_ECM_ADHESION": ["FN1","ITGA5","ITGB1","COL1A1","COL1A2","COL3A1","SPARC","VCAN","THBS1","LAMC1","LAMA4"],
}


def log(x):
    print(f"Stage 2.10: {x}", flush=True)


def finite(X, fill=None):
    X = np.asarray(X, float).copy()
    if X.ndim != 2:
        raise ValueError("Stage 2.10 expects a 2D matrix")
    if fill is None:
        med = np.zeros(X.shape[1], float)
        for j in range(X.shape[1]):
            v = X[np.isfinite(X[:, j]), j]
            med[j] = float(np.median(v)) if len(v) else 0.0
    else:
        med = np.asarray(fill, float)
        if len(med) != X.shape[1]:
            raise ValueError("Stage 2.10 imputation vector has wrong length")
        med = np.where(np.isfinite(med), med, 0.0)
    bad = ~np.isfinite(X)
    if bad.any():
        ii = np.where(bad)
        X[ii] = med[ii[1]]
    return X, med


def normalize_time(t):
    t = np.asarray(t, float)
    lo, hi = np.min(t), np.max(t)
    return np.zeros_like(t) if hi <= lo else (t - lo) / (hi - lo)


def resample_train(t, X, fill=None):
    tn = normalize_time(t)
    X, med = finite(X, fill)
    Y = np.empty((len(GRID), X.shape[1]))
    for j in range(X.shape[1]):
        Y[:, j] = np.interp(GRID, tn, X[:, j])
    return Y, med


def load_trajectories():
    # Reuse the already validated Stage 2.7 sample/time mapping. This is
    # important for GSE28688, where row order resolves ambiguous GSM labels.
    from dynamics.validation import _load_common_space
    matrix, meta = _load_common_space()
    meta = meta[meta["dataset"].astype(str).isin(TARGET)].copy()
    out = {}
    for ds in TARGET:
        g = meta[meta["dataset"].astype(str).eq(ds)].copy()
        g["time_hours"] = pd.to_numeric(g["time_hours"], errors="coerce")
        g = g[g["time_hours"].notna() & g["matrix_column"].notna()]
        if g["time_hours"].nunique() < 3:
            continue
        cols = g["matrix_column"].astype(str).tolist()
        X = matrix.loc[:, cols].T.copy()
        X.index = g["time_hours"].to_numpy(float)
        X = X.groupby(level=0, sort=True).mean()
        out[ds] = (X.index.to_numpy(float), X.to_numpy(float), [str(v) for v in matrix.index])
        log(f"{ds}: {len(g)} timed samples -> {len(X)} unique timepoints")
    return out


def training_imputation(train, n_genes):
    pieces = [np.asarray(v[1], float) for v in train.values()]
    stacked = np.vstack(pieces)
    return finite(stacked)[1]


def program_activity(X, genes, fill):
    X, _ = finite(X, fill)
    frame = pd.DataFrame(X.T, index=[str(g).upper() for g in genes])
    frame = frame.rank(axis=0, pct=True)
    result = {}
    for pid, members in PROGRAMS.items():
        idx = [g for g in members if g in frame.index]
        if len(idx) >= 3:
            result[pid] = frame.loc[idx].mean(axis=0).to_numpy(float)
    return result


def smooth_peak_time(activity):
    # Transition time is defined from the maximum absolute first derivative.
    # This is robust to monotone programs and avoids assuming a common peak.
    y = np.asarray(activity, float)
    d = np.gradient(y, GRID)
    k = int(np.argmax(np.abs(d)))
    return float(GRID[k]), float(d[k])


def pairwise_ordering(times):
    ids = list(times)
    rows = []
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            ta, tb = times[a], times[b]
            if not (np.isfinite(ta) and np.isfinite(tb)) or abs(ta - tb) < 1e-9:
                continue
            rows.append({
                "program_a": a,
                "program_b": b,
                "delta_transition": float(ta - tb),
                "ordering": int(np.sign(ta - tb)),
            })
    return pd.DataFrame(rows)


def train_module_consensus(train):
    # Strict training-only definition: each fold gets its own gene activities,
    # transition times, and consensus order. No held-out values enter selection.
    ds_features = {}
    for ds, (t, X, genes) in train.items():
        fill = training_imputation(train, X.shape[1])
        Y, _ = resample_train(t, X, fill)
        acts = program_activity(Y, genes, fill)
        feats = {}
        for pid, a in acts.items():
            tt, slope = smooth_peak_time(a)
            feats[pid] = {"transition_time": tt, "transition_slope": slope, "activity": a}
        ds_features[ds] = feats

    common = sorted(set.intersection(*(set(v) for v in ds_features.values()))) if ds_features else []
    if len(common) < 3:
        return ds_features, pd.DataFrame(), common

    rows = []
    for pid in common:
        vals = [ds_features[ds][pid]["transition_time"] for ds in ds_features]
        slopes = [ds_features[ds][pid]["transition_slope"] for ds in ds_features]
        rows.append({
            "program_id": pid,
            "n_training_datasets": len(vals),
            "mean_transition_time": float(np.mean(vals)),
            "transition_time_sd": float(np.std(vals)),
            "mean_transition_slope": float(np.mean(slopes)),
        })
    consensus = pd.DataFrame(rows).sort_values(["mean_transition_time", "program_id"]).reset_index(drop=True)
    consensus["consensus_rank"] = np.arange(1, len(consensus) + 1)
    return ds_features, consensus, common


def rank_corr(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]), method="spearman"))


def evaluate_heldout(heldout, train_names, traj):
    train = {d: traj[d] for d in train_names}
    ht, hX, hgenes = traj[heldout]
    ds_features, consensus, common = train_module_consensus(train)
    if len(common) < 3:
        return None, None, None

    # Held-out projection uses the training-defined programs only. Gene-level
    # medians are computed from training data and then applied to heldout.
    fill = training_imputation(train, hX.shape[1])
    hY, _ = resample_train(ht, hX, fill)
    hActs = program_activity(hY, hgenes, fill)
    held_times = {}
    held_slopes = {}
    for pid in common:
        if pid not in hActs:
            continue
        tt, slope = smooth_peak_time(hActs[pid])
        held_times[pid] = tt
        held_slopes[pid] = slope

    ids = [p for p in consensus.program_id if p in held_times]
    if len(ids) < 3:
        return None, None, None
    train_times = consensus.set_index("program_id").loc[ids, "mean_transition_time"].to_numpy(float)
    test_times = np.asarray([held_times[p] for p in ids], float)
    time_spearman = rank_corr(train_times, test_times)
    time_pearson = float(pd.Series(train_times).corr(pd.Series(test_times))) if np.std(train_times) > 1e-12 and np.std(test_times) > 1e-12 else np.nan

    order = pairwise_ordering({p: held_times[p] for p in ids})
    train_order = pairwise_ordering({p: consensus.set_index("program_id").loc[p, "mean_transition_time"] for p in ids})
    merged = order.merge(train_order, on=["program_a", "program_b"], suffixes=("_heldout", "_train"))
    ordering_agreement = float(np.mean(merged["ordering_heldout"].to_numpy() == merged["ordering_train"].to_numpy())) if len(merged) else np.nan

    # A compact biological interpretation: transition-time ranks should agree
    # while program trajectories remain related to the module activity itself.
    slopes = np.asarray([held_slopes[p] for p in ids], float)
    slope_consistency = float(np.mean(np.sign(slopes) == np.sign(
        consensus.set_index("program_id").loc[ids, "mean_transition_slope"].to_numpy(float)
    ))) if len(slopes) else np.nan

    row = {
        "heldout_dataset": heldout,
        "n_training_datasets": len(train_names),
        "n_common_programs": len(ids),
        "transition_rank_spearman": time_spearman,
        "transition_time_pearson": time_pearson,
        "pairwise_ordering_agreement": ordering_agreement,
        "slope_direction_agreement": slope_consistency,
    }
    detail = pd.DataFrame({
        "heldout_dataset": heldout,
        "program_id": ids,
        "train_consensus_transition_time": train_times,
        "heldout_transition_time": test_times,
        "train_consensus_rank": consensus.set_index("program_id").loc[ids, "consensus_rank"].to_numpy(int),
        "heldout_rank": pd.Series(test_times).rank(method="average").to_numpy(float),
        "heldout_transition_slope": slopes,
    })
    return row, detail, consensus


def bootstrap_fold(heldout, train_names, traj, n=N_BOOT):
    # Bootstrap module membership within each fixed biological program. The
    # program definitions remain fixed; only genes inside each program vary.
    train = {d: traj[d] for d in train_names}
    ht, hX, hgenes = traj[heldout]
    fill = training_imputation(train, hX.shape[1])
    hY, _ = resample_train(ht, hX, fill)
    gene_pos = {str(g).upper(): i for i, g in enumerate(hgenes)}
    rows = []
    rng = np.random.default_rng(210000 + sum(ord(c) for c in heldout))
    for b in range(n):
        times = {}
        for pid, members in PROGRAMS.items():
            available = [gene_pos[g] for g in members if g in gene_pos]
            if len(available) < 3:
                continue
            sample = rng.choice(available, size=len(available), replace=True)
            activity = np.nanmean(pd.DataFrame(hY[:, sample]).rank(axis=0, pct=True).to_numpy(float), axis=1)
            times[pid] = smooth_peak_time(activity)[0]
        if len(times) >= 3:
            vals = np.asarray([times[p] for p in sorted(times)], float)
            rows.append({"bootstrap": b + 1, "n_programs": len(times), "mean_transition_time": float(np.mean(vals)), "transition_time_sd": float(np.std(vals))})
    return pd.DataFrame(rows)


def permutation_null(observed_rows, traj, n=N_PERM):
    # Null: independently permute normalized time positions within each held-out
    # program trajectory, destroying temporal ordering while preserving values.
    rows = []
    rng = np.random.default_rng(210031)
    for b in range(n):
        fold_agreements = []
        for heldout in TARGET:
            if heldout not in traj:
                continue
            train_names = [d for d in TARGET if d != heldout and d in traj]
            train = {d: traj[d] for d in train_names}
            ht, hX, hgenes = traj[heldout]
            _, consensus, common = train_module_consensus(train)
            if len(common) < 3:
                continue
            fill = training_imputation(train, hX.shape[1])
            hY, _ = resample_train(ht, hX, fill)
            acts = program_activity(hY, hgenes, fill)
            perm_times = {}
            for pid in common:
                if pid not in acts:
                    continue
                perm = acts[pid][rng.permutation(len(acts[pid]))]
                perm_times[pid] = smooth_peak_time(perm)[0]
            ids = [p for p in consensus.program_id if p in perm_times]
            if len(ids) >= 3:
                train_times = consensus.set_index("program_id").loc[ids, "mean_transition_time"].to_numpy(float)
                test_times = np.asarray([perm_times[p] for p in ids], float)
                fold_agreements.append(rank_corr(train_times, test_times))
        rows.append({"permutation": b + 1, "mean_transition_rank_spearman": float(np.nanmean(fold_agreements)) if fold_agreements else np.nan})
    return pd.DataFrame(rows)


def run():
    traj = load_trajectories()
    names = [d for d in TARGET if d in traj]
    log(f"trajectory datasets: {', '.join(names)}")
    if len(names) < 3:
        raise RuntimeError("Stage 2.10 requires at least three trajectory datasets.")

    fold_rows = []
    detail_frames = []
    consensus_frames = []
    bootstrap_frames = []
    for heldout in names:
        train_names = [d for d in names if d != heldout]
        log(f"LODO held out {heldout}; training on {', '.join(train_names)}")
        row, detail, consensus = evaluate_heldout(heldout, train_names, traj)
        if row is None:
            continue
        fold_rows.append(row)
        detail_frames.append(detail)
        c = consensus.copy(); c.insert(0, "heldout_dataset", heldout); consensus_frames.append(c)
        boot = bootstrap_fold(heldout, train_names, traj)
        if not boot.empty:
            boot.insert(0, "heldout_dataset", heldout); bootstrap_frames.append(boot)

    R = pd.DataFrame(fold_rows)
    if R.empty:
        raise RuntimeError("No valid Stage 2.10 LODO folds.")
    D = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    C = pd.concat(consensus_frames, ignore_index=True) if consensus_frames else pd.DataFrame()
    B = pd.concat(bootstrap_frames, ignore_index=True) if bootstrap_frames else pd.DataFrame()
    R.to_csv(OUT / "01_lodo_module_conservation.csv", index=False)
    D.to_csv(OUT / "02_module_transition_order.csv", index=False)
    C.to_csv(OUT / "03_training_consensus_modules.csv", index=False)
    B.to_csv(OUT / "04_bootstrap_transition_stability.csv", index=False)

    P = permutation_null(R)
    P.to_csv(OUT / "05_time_permutation_null.csv", index=False)
    obs = float(R.transition_rank_spearman.mean())
    pv = P.mean_transition_rank_spearman.dropna().to_numpy(float)
    # Two-sided rank-correlation null.
    p = float((1 + np.sum(np.abs(pv) >= abs(obs))) / (len(pv) + 1)) if len(pv) else np.nan
    mean_order = float(R.pairwise_ordering_agreement.mean())
    mean_slope = float(R.slope_direction_agreement.mean())
    supported = bool(len(R) == len(names) and obs > 0.5 and mean_order > 0.65 and np.isfinite(p) and p < 0.05)

    summary = pd.DataFrame([{
        "n_trajectory_datasets": len(names),
        "n_valid_lodo_folds": len(R),
        "mean_transition_rank_spearman": obs,
        "mean_transition_time_pearson": float(R.transition_time_pearson.mean()),
        "mean_pairwise_ordering_agreement": mean_order,
        "mean_slope_direction_agreement": mean_slope,
        "mean_common_programs": float(R.n_common_programs.mean()),
        "bootstrap_folds_with_results": int(B.heldout_dataset.nunique()) if not B.empty else 0,
        "time_permutation_p": p,
        "conserved_transition_modules_supported": supported,
        "stage3_readiness": False,
        "interpretation": "LODO test of fixed biological programs as conserved transition modules; module definitions and gene imputation are training-only; bootstrap tests membership stability and time permutation tests temporal ordering; no ODE/state-space model.",
    }])
    summary.to_csv(OUT / "06_stage210_summary.csv", index=False)
    log("overall:")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    run()
