"""Stage 2.9.19: diagnose common-state quality and cross-dataset harmonization.

This stage does not fit an ODE or claim biological equivalence. It asks whether
GSE67462, GSE28688 and GSE297234 occupy a reasonably comparable state space.
Diagnostics include matched-time agreement, PCA dataset mixing, dataset
predictability, and gene-level dataset-vs-time variance decomposition.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_19"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]

def log(x):
    print(f"Stage 2.9.19: {x}", flush=True)

def _load():
    from dynamics.validation import _load_common_space
    matrix, metadata = _load_common_space()
    metadata = metadata[metadata.dataset.isin(TARGET)].copy()
    cols = [c for c in metadata.matrix_column if c in matrix.columns]
    metadata = metadata[metadata.matrix_column.isin(cols)].copy()
    return matrix.loc[:, cols], metadata

def _sample_table(matrix, meta):
    rows = []
    for ds, g in meta.groupby("dataset"):
        for t, gt in g.groupby("time_hours"):
            cols = gt.matrix_column.tolist()
            if not cols: continue
            v = matrix[cols].mean(axis=1).to_numpy(float)
            rows.append({"dataset": ds, "time_hours": float(t), "n_samples": len(cols), "vector": v})
    return rows

def matched_time_agreement(matrix, meta):
    rows = []
    by = {(r["dataset"], r["time_hours"]): r["vector"] for r in _sample_table(matrix, meta)}
    pairs = [(TARGET[i], TARGET[j]) for i in range(len(TARGET)) for j in range(i + 1, len(TARGET))]
    for a, b in pairs:
        common = sorted(set(t for (ds, t) in by if ds == a) & set(t for (ds, t) in by if ds == b))
        for t in common:
            x, y = by[(a, t)], by[(b, t)]
            ok = np.isfinite(x) & np.isfinite(y)
            r = np.corrcoef(x[ok], y[ok])[0, 1] if ok.sum() > 2 and np.std(x[ok]) > 0 and np.std(y[ok]) > 0 else np.nan
            rmse = float(np.sqrt(np.mean((x[ok] - y[ok]) ** 2))) if ok.any() else np.nan
            rows.append({"dataset_a": a, "dataset_b": b, "time_hours": t, "n_genes": int(ok.sum()), "pearson": r, "rmse": rmse})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "01_matched_time_agreement.csv", index=False)
    return out

def pca_mixing(matrix, meta, n_components=10):
    X = matrix.T.to_numpy(float)
    X = np.nan_to_num(X, nan=np.nanmedian(X, axis=0))
    X = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(n_components, X.shape[0]-1, X.shape[1])).fit_transform(X)
    ds = meta.dataset.to_numpy()
    # kNN mixing: fraction of neighbours from the same dataset. Lower is better.
    k = min(5, len(ds) - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(pca)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    same = np.mean(np.array([[ds[j] == ds[i] for j in row] for i, row in enumerate(idx)]))
    rows = [{"metric": "mean_same_dataset_knn_fraction", "value": float(same), "k": k},
            {"metric": "expected_same_dataset_fraction", "value": float(np.mean(pd.Series(ds).value_counts(normalize=True).to_numpy() ** 1)), "k": k}]
    coords = pd.DataFrame(pca, columns=[f"PC{i+1}" for i in range(pca.shape[1])])
    coords["dataset"] = ds
    coords["sample"] = meta["sample"].to_numpy()
    coords["time_hours"] = meta["time_hours"].to_numpy()
    coords.to_csv(OUT / "02_pca_coordinates.csv", index=False)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "03_pca_mixing.csv", index=False)
    return out, coords

def dataset_predictability(coords):
    X = coords.filter(regex=r"^PC\\d+$").to_numpy(float)
    y = coords.dataset.to_numpy()
    result = {"accuracy": np.nan, "n_samples": len(y), "n_datasets": len(np.unique(y)), "status": "insufficient"}
    counts = pd.Series(y).value_counts()
    if len(np.unique(y)) >= 2 and counts.min() >= 2:
        cv = StratifiedKFold(n_splits=int(min(3, counts.min())), shuffle=True, random_state=2919)
        clf = LogisticRegression(max_iter=2000, multi_class="auto")
        scores = cross_val_score(clf, X[:, :min(10, X.shape[1])], y, cv=cv)
        result.update(accuracy=float(scores.mean()), status="ok")
    out = pd.DataFrame([result])
    out.to_csv(OUT / "04_dataset_predictability.csv", index=False)
    return out

def variance_decomposition(matrix, meta):
    # Gene-level ANOVA-like R2: dataset-only, time-only, and additive dataset+time.
    rows = []
    d = pd.get_dummies(meta.dataset, dtype=float).to_numpy()
    t = meta.time_hours.to_numpy(float)
    finite_t = np.isfinite(t)
    tn = np.zeros_like(t)
    if finite_t.any() and np.ptp(t[finite_t]) > 0: tn[finite_t] = (t[finite_t] - np.mean(t[finite_t])) / np.std(t[finite_t])
    Xd = np.column_stack([np.ones(len(meta)), d[:, 1:]])
    Xt = np.column_stack([np.ones(len(meta)), tn])
    Xdt = np.column_stack([Xd, tn])
    for gene, y in matrix.loc[:, meta.matrix_column].iterrows():
        y = y.to_numpy(float); ok = np.isfinite(y) & finite_t
        if ok.sum() < 6 or np.var(y[ok]) <= 1e-12: continue
        yy = y[ok]
        def r2(X):
            xx = X[ok]; beta = np.linalg.lstsq(xx, yy, rcond=None)[0]; pred = xx @ beta
            return max(0.0, 1.0 - np.sum((yy-pred)**2) / np.sum((yy-yy.mean())**2))
        rows.append({"gene": gene, "r2_dataset": r2(Xd), "r2_time": r2(Xt), "r2_dataset_time": r2(Xdt)})
    out = pd.DataFrame(rows)
    if len(out):
        out["dataset_excess_r2"] = np.maximum(0, out.r2_dataset_time - out.r2_time)
        out = out.sort_values("dataset_excess_r2", ascending=False)
    out.to_csv(OUT / "05_gene_variance_decomposition.csv", index=False)
    return out

def pairwise_time_trajectory_similarity(matrix, meta):
    rows = []
    trajectories = {}
    for ds, g in meta.groupby("dataset"):
        tg = g.groupby("time_hours").matrix_column.apply(list)
        trajectories[ds] = {float(t): matrix[cols].mean(axis=1).to_numpy(float) for t, cols in tg.items()}
    for i, a in enumerate(TARGET):
        for b in TARGET[i+1:]:
            common = sorted(set(trajectories.get(a, {})) & set(trajectories.get(b, {})))
            if len(common) < 3: continue
            cors = []
            for t in common:
                x, y = trajectories[a][t], trajectories[b][t]; ok=np.isfinite(x)&np.isfinite(y)
                if ok.sum()>2 and np.std(x[ok])>0 and np.std(y[ok])>0: cors.append(np.corrcoef(x[ok],y[ok])[0,1])
            rows.append({"dataset_a":a,"dataset_b":b,"n_common_times":len(common),"mean_matched_time_gene_correlation":float(np.mean(cors)) if cors else np.nan,"median_matched_time_gene_correlation":float(np.median(cors)) if cors else np.nan})
    out=pd.DataFrame(rows); out.to_csv(OUT/"06_trajectory_similarity.csv",index=False); return out

def run():
    log("loading common gene space")
    matrix, meta = _load()
    log(f"datasets={','.join(TARGET)}; genes={matrix.shape[0]}; samples={matrix.shape[1]}")
    agreement = matched_time_agreement(matrix, meta)
    mix, coords = pca_mixing(matrix, meta)
    pred = dataset_predictability(coords)
    var = variance_decomposition(matrix, meta)
    traj = pairwise_time_trajectory_similarity(matrix, meta)
    summary = {
        "n_genes": matrix.shape[0], "n_samples": matrix.shape[1], "n_datasets": meta.dataset.nunique(),
        "n_matched_time_cases": len(agreement),
        "mean_matched_time_pearson": float(agreement.pearson.mean()) if len(agreement) else np.nan,
        "median_matched_time_pearson": float(agreement.pearson.median()) if len(agreement) else np.nan,
        "mean_same_dataset_knn_fraction": float(mix.loc[mix.metric=="mean_same_dataset_knn_fraction","value"].iloc[0]),
        "dataset_predictability_accuracy": float(pred.accuracy.iloc[0]) if len(pred) else np.nan,
        "median_gene_r2_dataset": float(var.r2_dataset.median()) if len(var) else np.nan,
        "median_gene_r2_time": float(var.r2_time.median()) if len(var) else np.nan,
        "median_dataset_excess_r2": float(var.dataset_excess_r2.median()) if len(var) else np.nan,
        "n_genes_dataset_excess_r2_gt_0_10": int((var.dataset_excess_r2 > .10).sum()) if len(var) else 0,
        "n_pairs_with_trajectory_similarity": len(traj),
        "harmonization_interpretation": "inspect matched-time agreement, kNN mixing, dataset predictability and gene-level excess dataset variance; no automatic ODE gate",
        "stage3_readiness": False,
    }
    out=pd.DataFrame([summary]); out.to_csv(OUT/"07_STAGE2_9_19_SUMMARY.csv",index=False)
    log("complete")
    print("\nStage 2.9.19 summary:")
    print(out.to_string(index=False))
    if len(agreement): print("\nMatched-time agreement:"); print(agreement.groupby(["dataset_a","dataset_b"]).pearson.agg(["count","mean","median"]).to_string())
    if len(traj): print("\nTrajectory similarity:"); print(traj.to_string(index=False))
    return out

if __name__ == "__main__": run()
