"""Stage 2.9.29: trajectory family discovery.

Tests whether the trajectory datasets are better described by one universal
trajectory family or by multiple reproducible families.  The analysis is
strictly representation-level: no ODE/state-space model is fitted and no
family is interpreted as a biological state without independent validation.

The unit of clustering is the dataset trajectory, not individual samples.
With the current three eligible trajectory datasets this is deliberately a
small-sample diagnostic.  Bootstrap gene resampling assesses robustness and
time-permutation provides a null for temporal shape structure.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_29"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
GRID = np.linspace(0, 1, 20)
N_BOOT = 500
N_PERM = 1000
MIN_GENES = 100
BOOT_GENE_FRACTION = 0.7


def log(x):
    print(f"Stage 2.9.29: {x}", flush=True)


def finite_matrix(X):
    X = np.asarray(X, float).copy()
    if np.isfinite(X).all():
        return X
    for j in range(X.shape[1]):
        ok = np.isfinite(X[:, j])
        fill = float(np.median(X[ok, j])) if ok.any() else 0.0
        X[~ok, j] = fill
    return X


def corr(a, b):
    a = np.asarray(a, float).ravel()
    b = np.asarray(b, float).ravel()
    n = min(a.size, b.size)
    a, b = a[:n], b[:n]
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok])))


def load():
    from dynamics.validation import _load_common_space
    m, meta = _load_common_space()
    meta = meta[meta["dataset"].astype(str).isin(TARGET)].copy()
    meta["matrix_column"] = meta["matrix_column"].astype(str)
    cols = [c for c in meta["matrix_column"] if c in m.columns]
    meta = meta[meta["matrix_column"].isin(cols)].copy()
    if "time_hours" not in meta.columns:
        from dynamics.validation import _time_hours_for_validation, _strip_dataset_prefix
        meta["time_hours"] = [
            _time_hours_for_validation(
                str(r["dataset"]),
                _strip_dataset_prefix(str(r["sample"])),
                i if str(r["dataset"]) == "GSE28688" else None,
            )
            for i, (_, r) in enumerate(meta.iterrows())
        ]
    meta["time_hours"] = pd.to_numeric(meta["time_hours"], errors="coerce")
    return m.loc[:, cols], meta


def dataset_trajectories(m, meta):
    out = {}
    for ds in TARGET:
        g = meta[meta["dataset"].astype(str).eq(ds)].copy()
        g = g[np.isfinite(g["time_hours"])].copy()
        if g["time_hours"].nunique() < 3:
            continue
        cols = list(g["matrix_column"])
        X = m.loc[:, cols].T.copy()
        X["time_hours"] = g["time_hours"].to_numpy(float)
        X = X.groupby("time_hours", sort=True).mean()
        times = X.index.to_numpy(float)
        values = finite_matrix(X.to_numpy(float))
        if len(times) >= 3:
            out[ds] = (times, values, list(m.index))
    return out


def normalize_trajectory(t, X, gene_idx=None):
    if gene_idx is not None:
        X = X[:, gene_idx]
    X = finite_matrix(X)
    tn = (t - t.min()) / (t.max() - t.min())
    Y = np.empty((len(GRID), X.shape[1]), float)
    for j in range(X.shape[1]):
        Y[:, j] = np.interp(GRID, tn, X[:, j])
    Y = Y - Y[0:1, :]
    amp = np.sqrt(np.mean(Y * Y, axis=0))
    amp = np.where(amp > 1e-8, amp, 1.0)
    return Y / amp


def gene_variance_scores(traj, gene_idx=None):
    blocks = []
    for t, X, _ in traj.values():
        blocks.append(normalize_trajectory(t, X, gene_idx))
    A = np.stack(blocks)
    mean = A.mean(axis=0)
    hetero = np.mean((A - mean[None, :, :]) ** 2, axis=(0, 1))
    shared_var = np.var(mean, axis=0)
    score = shared_var / (shared_var + hetero + 1e-12)
    return score


def select_genes(traj, gene_idx=None):
    score = gene_variance_scores(traj, gene_idx)
    all_genes = np.asarray(next(iter(traj.values()))[2])
    genes = all_genes if gene_idx is None else all_genes[gene_idx]
    order = np.argsort(score)[::-1]
    chosen = order[score[order] >= 0.5]
    if len(chosen) < MIN_GENES:
        chosen = order[:min(MIN_GENES, len(order))]
    return genes[chosen]


def trajectory_vectors(traj, genes, gene_idx=None):
    allg = list(next(iter(traj.values()))[2])
    idx = [allg.index(str(g)) for g in genes]
    return np.vstack([normalize_trajectory(traj[ds][0], traj[ds][1], idx).ravel() for ds in TARGET])


def distance_matrix(vectors):
    V = np.asarray(vectors, float)
    norms = np.linalg.norm(V, axis=1, keepdims=True)
    V = V / np.where(norms > 1e-12, norms, 1.0)
    D = np.sqrt(np.maximum(0.0, 2.0 - 2.0 * (V @ V.T)))
    np.fill_diagonal(D, 0.0)
    return D


def labels_for_k(D, k):
    n = D.shape[0]
    if k < 2 or k >= n:
        return None
    model = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage="average")
    return model.fit_predict(D)


def cluster_quality(D, labels):
    n = len(labels)
    if labels is None or len(np.unique(labels)) < 2 or len(np.unique(labels)) >= n:
        return np.nan
    try:
        return float(silhouette_score(D, labels, metric="precomputed"))
    except Exception:
        return np.nan


def best_k(D):
    n = D.shape[0]
    rows = []
    for k in range(2, min(3, n - 1) + 1):
        labels = labels_for_k(D, k)
        rows.append({"n_families": k, "silhouette": cluster_quality(D, labels), "labels": labels})
    if not rows:
        return 1, np.nan, np.zeros(n, dtype=int)
    valid = [r for r in rows if np.isfinite(r["silhouette"])]
    chosen = max(valid, key=lambda r: r["silhouette"]) if valid else rows[0]
    return int(chosen["n_families"]), float(chosen["silhouette"]) if np.isfinite(chosen["silhouette"]) else np.nan, chosen["labels"]


def partition_agreement(reference, candidate):
    R = np.asarray(reference)
    C = np.asarray(candidate)
    A = R[:, None] == R[None, :]
    B = C[:, None] == C[None, :]
    return float(np.mean(A == B))


def run():
    m, meta = load()
    traj = dataset_trajectories(m, meta)
    names = [d for d in TARGET if d in traj]
    log(f"trajectory datasets: {', '.join(names)}")
    if len(names) < 3:
        raise RuntimeError("Stage 2.9.29 requires at least three trajectory datasets.")
    traj = {d: traj[d] for d in names}

    genes = select_genes(traj)
    V = trajectory_vectors(traj, genes)
    D = distance_matrix(V)
    k, sil, labels = best_k(D)

    pair_rows = []
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            pair_rows.append({"dataset_a": a, "dataset_b": names[j], "shape_distance": float(D[i, j]), "shape_similarity": float(1.0 - D[i, j] / 2.0)})
    pd.DataFrame(pair_rows).to_csv(OUT / "01_trajectory_pairwise_distances.csv", index=False)

    observed = pd.DataFrame([{
        "n_trajectory_datasets": len(names),
        "n_selected_genes": len(genes),
        "n_families": k,
        "silhouette_score": sil,
        "family_partition": ";".join(f"{names[i]}=F{int(labels[i]) + 1}" for i in range(len(names))),
        "mean_shape_distance": float(D[np.triu_indices_from(D, 1)].mean()),
    }])
    observed.to_csv(OUT / "02_observed_family_structure.csv", index=False)
    pd.DataFrame({"dataset": names, "family": labels + 1}).to_csv(OUT / "03_observed_family_assignments.csv", index=False)

    rng = np.random.default_rng(292900)
    boot_rows = []
    all_gene_count = len(next(iter(traj.values()))[2])
    for b in range(N_BOOT):
        n_genes = max(MIN_GENES, int(all_gene_count * BOOT_GENE_FRACTION))
        idx = rng.choice(all_gene_count, size=n_genes, replace=False)
        bg = select_genes(traj, idx)
        BV = trajectory_vectors(traj, bg)
        BD = distance_matrix(BV)
        bk, bsil, blabels = best_k(BD)
        agreement = partition_agreement(labels, blabels) if blabels is not None else np.nan
        boot_rows.append({"bootstrap": b + 1, "n_selected_genes": len(bg), "n_families": bk, "silhouette_score": bsil, "partition_agreement": agreement})
    B = pd.DataFrame(boot_rows)
    B.to_csv(OUT / "04_bootstrap_family_stability.csv", index=False)

    perm_rows = []
    for p in range(N_PERM):
        PV = []
        for di, ds in enumerate(names):
            t, X, allg = traj[ds]
            perm = np.random.default_rng(400000 + p * len(names) + di).permutation(len(t))
            # Permute expression rows while keeping the chronological grid fixed.
            # This destroys temporal ordering without creating non-monotonic xp.
            XP = X[perm, :]
            idx = [allg.index(str(g)) for g in genes]
            PV.append(normalize_trajectory(t, XP, idx).ravel())
        PD = distance_matrix(np.vstack(PV))
        pk, psil, plabels = best_k(PD)
        perm_rows.append({"permutation": p + 1, "n_families": pk, "silhouette_score": psil, "partition_agreement_to_observed": partition_agreement(labels, plabels) if plabels is not None else np.nan})
    P = pd.DataFrame(perm_rows)
    P.to_csv(OUT / "05_time_permutation_null.csv", index=False)

    obs_sil = float(sil) if np.isfinite(sil) else np.nan
    null_sil = P["silhouette_score"].to_numpy(float)
    null_sil = null_sil[np.isfinite(null_sil)]
    p_sil = float((1 + np.sum(null_sil >= obs_sil)) / (len(null_sil) + 1)) if np.isfinite(obs_sil) and len(null_sil) else np.nan
    boot_family2 = float(np.mean(B["n_families"].to_numpy() == 2)) if len(B) else np.nan
    boot_agreement = float(B["partition_agreement"].mean()) if len(B) else np.nan
    boot_p05 = float(B["partition_agreement"].quantile(0.05)) if len(B) else np.nan
    multiple_supported = bool(k >= 2 and np.isfinite(obs_sil) and obs_sil > 0.3 and boot_family2 >= 0.8 and boot_p05 >= 0.7 and np.isfinite(p_sil) and p_sil < 0.05)
    universal_supported = bool(k == 1 or (np.isfinite(obs_sil) and obs_sil <= 0.3) or (np.isfinite(p_sil) and p_sil >= 0.05))

    summary = pd.DataFrame([{
        "n_trajectory_datasets": len(names),
        "observed_n_families": k,
        "observed_silhouette": obs_sil,
        "bootstrap_fraction_two_family_solutions": boot_family2,
        "mean_bootstrap_partition_agreement": boot_agreement,
        "bootstrap_partition_agreement_p05": boot_p05,
        "time_permutation_p_silhouette": p_sil,
        "multiple_trajectory_families_supported": multiple_supported,
        "single_universal_trajectory_supported": universal_supported,
        "stage3_readiness": False,
        "interpretation": "Dataset-level trajectory-shape clustering on a common normalized time grid; gene bootstrap and time-permutation null; representation diagnostic only; no ODE/state-space model",
    }])
    summary.to_csv(OUT / "06_stage2929_summary.csv", index=False)
    log("overall:")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    run()
