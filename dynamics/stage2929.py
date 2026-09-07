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
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_29"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
GRID = np.linspace(0, 1, 20)
N_COMPONENTS = 5
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
    # Remove the dataset-specific starting level and scale each gene by its
    # trajectory amplitude. This compares temporal shape rather than absolute
    # expression level.
    Y = Y - Y[0:1, :]
    amp = np.sqrt(np.mean(Y * Y, axis=0))
    amp = np.where(amp > 1e-8, amp, 1.0)
    return Y / amp


def gene_variance_scores(traj, gene_idx=None):
    blocks = []
    for t, X, _ in traj.values():
        blocks.append(normalize_trajectory(t, X, gene_idx))
    A = np.stack(blocks)  # dataset x grid x gene
    mean = A.mean(axis=0)
    hetero = np.mean((A - mean[None, :, :]) ** 2, axis=(0, 1))
    shared_var = np.var(mean, axis=0)
    score = shared_var / (shared_var + hetero + 1e-12)
    return score, mean, A


def select_genes(traj, gene_idx=None):
    score, mean, A = gene_variance_scores(traj, gene_idx)
    all_genes = list(next(iter(traj.values()))[2])
    if gene_idx is None:
        genes = np.asarray(all_genes)
    else:
        genes = np.asarray(all_genes)[gene_idx]
    order = np.argsort(score)[::-1]
    chosen_local = order[score[order] >= 0.5]
    if len(chosen_local) < MIN_GENES:
        chosen_local = order[:min(MIN_GENES, len(order))]
    return genes[chosen_local], score[chosen_local], mean[:, chosen_local], A[:, :, chosen_local]


def trajectory_vectors(traj, genes):
    vectors = []
    allg = list(next(iter(traj.values()))[2])
    idx = [allg.index(str(g)) for g in genes]
    for ds in TARGET:
        t, X, _ = traj[ds]
        Y = normalize_trajectory(t, X, idx)
        # PCA is fitted across the pooled training trajectories only. Here the
        # caller supplies the complete set of trajectories used for clustering.
        vectors.append(Y.ravel())
    return np.vstack(vectors)


def distance_matrix(vectors):
    V = np.asarray(vectors, float)
    # Euclidean shape distance after per-trajectory L2 normalization. This is
    # invariant to overall amplitude and preserves the common normalized time grid.
    norms = np.linalg.norm(V, axis=1, keepdims=True)
    V = V / np.where(norms > 1e-12, norms, 1.0)
    D = np.sqrt(np.maximum(0.0, 2.0 - 2.0 * (V @ V.T)))
    np.fill_diagonal(D, 0.0)
    return D


def labels_for_k(D, k):
    n = D.shape[0]
    if k < 2 or k >= n + 1:
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
    # With n=3, k=2 is the only meaningful non-degenerate partition.
    valid = [r for r in rows if np.isfinite(r["silhouette"])]
    if valid:
        chosen = max(valid, key=lambda r: r["silhouette"])
    else:
        chosen = rows[0]
    return int(chosen["n_families"]), float(chosen["silhouette"]) if np.isfinite(chosen["silhouette"]) else np.nan, chosen["labels"]


def family_partition(labels):
    return {int(k): tuple(np.where(np.asarray(labels) == k)[0].tolist()) for k in np.unique(labels)}


def partition_agreement(reference, candidate):
    # Pairwise co-clustering agreement, insensitive to arbitrary cluster labels.
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

    # Global descriptive representation is used only to describe the observed
    # family structure. Bootstrap and permutation analyses repeat the complete
    # representation construction from resampled genes.
    genes, score, _, _ = select_genes(traj)
    V = trajectory_vectors(traj, genes)
    D = distance_matrix(V)
    k, sil, labels = best_k(D)
    families = family_partition(labels)

    pair_rows = []
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            pair_rows.append({"dataset_a": a, "dataset_b": names[j], "shape_distance": float(D[i, j]), "shape_similarity": float(1.0 - D[i, j] / 2.0)})
    pd.DataFrame(pair_rows).to_csv(OUT / "01_trajectory_distance_matrix.csv", index=False)

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
    for b in range(N_BOOT):
        n_genes = max(MIN_GENES, int(len(next(iter(traj.values()))[2]) * BOOT_GENE_FRACTION))
        idx = rng.choice(len(next(iter(traj.values()))[2]), size=n_genes, replace=False)
        bg, _, _, _ = select_genes(traj, idx)
        BV = trajectory_vectors(traj, bg)
        BD = distance_matrix(BV)
        bk, bsil, blabels = best_k(BD)
        if blabels is None:
            agreement = np.nan
        else:
            agreement = partition_agreement(labels, blabels)
        boot_rows.append({"bootstrap": b + 1, "n_selected_genes": len(bg), "n_families": bk, "silhouette_score": bsil, "partition_agreement": agreement})
    B = pd.DataFrame(boot_rows)
    B.to_csv(OUT / "04_bootstrap_family_stability.csv", index=False)

    # Null: independently permute the temporal ordering within each dataset,
    # then rebuild the trajectory representation and clustering. A family is
    # considered more structured than null if its silhouette is unusually high.
    perm_rows = []
    for p in range(N_PERM):
        PV = []
        for ds in names:
            t, X, allg = traj[ds]
            # Preserve each dataset's number of observed timepoints but destroy
            # temporal ordering before interpolation.
            perm_t = t[np.random.default_rng(400000 + p * 10 + names.index(ds)).permutation(len(t))]
            PV.append(normalize_trajectory(perm_t, X[:, [allg.index(str(g)) for g in genes]]).ravel())
        PD = distance_matrix(np.vstack(PV))
        pk, psil, plabels = best_k(PD)
        perm_rows.append({"permutation": p + 1, "n_families": pk, "silhouette_score": psil, "partition_agreement_to_observed": partition_agreement(labels, plabels) if plabels is not None else np.nan})
    P = pd.DataFrame(perm_rows)
    P.to_csv(OUT / "05_time_permutation_null.csv", index=False)

    obs_sil = float(sil) if np.isfinite(sil) else np.nan
    null_sil = P["silhouette_score"].to_numpy(float)
    null_sil = null_sil[np.isfinite(null_sil)]
    if np.isfinite(obs_sil) and len(null_sil):
        p_sil = float((1 + np.sum(null_sil >= obs_sil)) / (len(null_sil) + 1))
    else:
        p_sil = np.nan
    boot_family2 = float(np.mean(B["n_families"].to_numpy() == 2)) if len(B) else np.nan
    boot_agreement = float(B["partition_agreement"].mean()) if len(B) else np.nan
    boot_p05 = float(B["partition_agreement"].quantile(0.05)) if len(B) else np.nan
    multiple_supported = bool(k >= 2 and np.isfinite(obs_sil) and obs_sil > 0.3 and boot_family2 >= 0.8 and boot_p05 >= 0.7 and np.isfinite(p_sil) and p_sil < 0.05)

    summary = pd.DataFrame([{
        "n_trajectory_datasets": len(names),
        "observed_n_families": k,
        "observed_silhouette": obs_sil,
        "bootstrap_fraction_two_family_solutions": boot_family2,
        "mean_bootstrap_partition_agreement": boot_agreement,
        "bootstrap_partition_agreement_p05": boot_p05,
        "time_permutation_p_silhouette": p_sil,
        "multiple_trajectory_families_supported": multiple_supported,
        "single_universal_trajectory_supported": bool(not multiple_supported and k == 1),
        "stage3_readiness": False,
        "interpretation": "Dataset-level trajectory-shape clustering on a common normalized time grid; gene bootstrap and time-permutation null; representation diagnostic only; no ODE/state-space model",
    }])
    summary.to_csv(OUT / "06_stage2929_summary.csv", index=False)
    log("overall:")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    run()
