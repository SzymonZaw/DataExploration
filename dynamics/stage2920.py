"""Stage 2.9.20: controlled repair of the common biological state space.

Compares the current common gene space with three harmonization variants:
1) baseline all genes,
2) genes with low dataset-specific excess variance,
3) genes whose time signal dominates dataset signal,
4) fixed biological program activities from Stage 2.9.14.

The stage is diagnostic only. It does not fit an ODE and never enables Stage 3.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_20"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]


def log(x): print(f"Stage 2.9.20: {x}", flush=True)


def load_space():
    from dynamics.validation import _load_common_space
    matrix, meta = _load_common_space()
    meta = meta[meta.dataset.isin(TARGET)].copy()
    cols = [c for c in meta.matrix_column if c in matrix.columns]
    meta = meta[meta.matrix_column.isin(cols)].copy()
    meta = meta.set_index("matrix_column").loc[cols].reset_index()
    return matrix.loc[:, cols], meta


def variance_table(matrix, meta):
    dummies = pd.get_dummies(meta.dataset, dtype=float).to_numpy()
    t = meta.time_hours.to_numpy(float)
    ok_t = np.isfinite(t)
    tn = np.zeros_like(t)
    if ok_t.any() and np.ptp(t[ok_t]) > 0:
        tn[ok_t] = (t[ok_t] - np.mean(t[ok_t])) / np.std(t[ok_t])
    Xd = np.column_stack([np.ones(len(meta)), dummies[:, 1:]])
    Xt = np.column_stack([np.ones(len(meta)), tn])
    Xdt = np.column_stack([Xd, tn])
    rows = []
    for gene, s in matrix.loc[:, meta.matrix_column].iterrows():
        y = s.to_numpy(float)
        good = np.isfinite(y) & ok_t
        if good.sum() < 6 or np.var(y[good]) <= 1e-12:
            continue
        yy = y[good]
        def r2(X):
            xx = X[good]
            beta = np.linalg.lstsq(xx, yy, rcond=None)[0]
            pred = xx @ beta
            den = np.sum((yy - yy.mean()) ** 2)
            return max(0.0, 1.0 - np.sum((yy - pred) ** 2) / den) if den > 0 else 0.0
        rd, rt, rdt = r2(Xd), r2(Xt), r2(Xdt)
        rows.append({"gene": gene, "r2_dataset": rd, "r2_time": rt,
                     "r2_dataset_time": rdt,
                     "dataset_excess_r2": max(0.0, rdt - rt),
                     "time_minus_dataset_r2": rt - rd})
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "01_gene_variance_diagnostics.csv", index=False)
    return out


def feature_sets(matrix, var):
    genes = set(matrix.index.astype(str))
    sets = {"baseline_all": sorted(genes)}
    if len(var):
        v = var.set_index("gene")
        sets["low_dataset_excess_q80"] = sorted(set(v.index[v.dataset_excess_r2 <= v.dataset_excess_r2.quantile(.80)]).intersection(genes))
        sets["low_dataset_excess_q70"] = sorted(set(v.index[v.dataset_excess_r2 <= v.dataset_excess_r2.quantile(.70)]).intersection(genes))
        sets["time_dominant"] = sorted(set(v.index[(v.r2_time >= v.r2_dataset) & (v.r2_time >= v.r2_time.quantile(.50))]).intersection(genes))
    for k, g in list(sets.items()):
        if len(g) < 100: sets[k] = sorted(genes)
    rows = [{"variant": k, "n_genes": len(g)} for k, g in sets.items()]
    pd.DataFrame(rows).to_csv(OUT / "02_feature_set_summary.csv", index=False)
    return sets


def aggregate_vectors(matrix, meta, genes):
    m = matrix.loc[genes, meta.matrix_column]
    rows = []
    for ds, g in meta.groupby("dataset"):
        for t, gt in g.groupby("time_hours"):
            cols = gt.matrix_column.tolist()
            if not cols: continue
            x = m[cols].mean(axis=1).to_numpy(float)
            rows.append({"dataset": ds, "time_hours": float(t), "vector": x})
    return rows


def evaluate_variant(matrix, meta, genes, variant):
    m = matrix.loc[genes, meta.matrix_column]
    X = m.T.to_numpy(float)
    X = np.nan_to_num(X, nan=np.nanmedian(X, axis=0))
    X = StandardScaler().fit_transform(X)
    ncomp = min(10, X.shape[0] - 1, X.shape[1])
    coords = PCA(n_components=ncomp).fit_transform(X)
    ds = meta.dataset.to_numpy()
    k = min(5, len(ds) - 1)
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    same = np.mean(np.array([[ds[j] == ds[i] for j in row] for i, row in enumerate(idx)]))
    props = pd.Series(ds).value_counts(normalize=True).to_numpy()
    expected = float(np.sum(props ** 2))

    # Dataset-vs-time variance after restricting to the candidate state genes.
    agg = aggregate_vectors(matrix, meta, genes)
    pairs = []
    by = {(r["dataset"], r["time_hours"]): r["vector"] for r in agg}
    for i, a in enumerate(TARGET):
        for b in TARGET[i + 1:]:
            common = sorted(set(t for d, t in by if d == a) & set(t for d, t in by if d == b))
            cors = []
            for t in common:
                x, y = by[(a, t)], by[(b, t)]
                good = np.isfinite(x) & np.isfinite(y)
                if good.sum() > 2 and np.std(x[good]) > 0 and np.std(y[good]) > 0:
                    cors.append(float(np.corrcoef(x[good], y[good])[0, 1]))
            if cors:
                pairs.append({"dataset_a": a, "dataset_b": b, "n_common_times": len(cors), "mean_pearson": np.mean(cors)})
    return {"variant": variant, "n_genes": len(genes), "same_dataset_knn": same,
            "expected_same_dataset_knn": expected, "knn_excess": same - expected,
            "mean_matched_time_pearson": np.mean([p["mean_pearson"] for p in pairs]) if pairs else np.nan,
            "n_pairs": len(pairs)}, pairs


def program_variant(matrix, meta):
    """Evaluate fixed Stage 2.9.14 programs without using discovered programs."""
    from dynamics.stage2914 import PROGRAMS, _activity
    rows = []
    for _, r in meta.iterrows():
        vals = matrix[r.matrix_column]
        activities = _activity(vals, PROGRAMS)
        rows.append({"dataset": r.dataset, "sample": r["sample"], "time_hours": r.time_hours,
                     **activities})
    a = pd.DataFrame(rows)
    a.to_csv(OUT / "05_fixed_program_activity.csv", index=False)
    cors = []
    for p in PROGRAMS:
        vals = []
        times = []
        for ds, g in a.groupby("dataset"):
            g = g[np.isfinite(g.time_hours)]
            if g.time_hours.nunique() < 3: continue
            vals.append(g[p].to_numpy(float)); times.append(g.time_hours.to_numpy(float))
        for i in range(len(vals)):
            for j in range(i + 1, len(vals)):
                # Only descriptive cross-dataset program profile correlation.
                n = min(len(vals[i]), len(vals[j]))
                if n >= 3 and np.std(vals[i][:n]) > 0 and np.std(vals[j][:n]) > 0:
                    cors.append({"program_id": p, "dataset_a": list(a.dataset.unique())[i], "dataset_b": list(a.dataset.unique())[j],
                                 "pearson": np.corrcoef(vals[i][:n], vals[j][:n])[0, 1]})
    out = pd.DataFrame(cors)
    out.to_csv(OUT / "06_fixed_program_cross_dataset_similarity.csv", index=False)
    return out


def run():
    log("loading common gene space")
    matrix, meta = load_space()
    log(f"datasets={','.join(TARGET)}; genes={matrix.shape[0]}; samples={matrix.shape[1]}")
    var = variance_table(matrix, meta)
    sets = feature_sets(matrix, var)
    results, pair_rows = [], []
    for variant, genes in sets.items():
        log(f"evaluating {variant}: {len(genes)} genes")
        r, pairs = evaluate_variant(matrix, meta, genes, variant)
        results.append(r)
        pair_rows.extend([{**p, "variant": variant} for p in pairs])
    summary = pd.DataFrame(results)
    summary.to_csv(OUT / "03_harmonization_variant_comparison.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(OUT / "04_matched_time_by_variant.csv", index=False)
    try:
        program_variant(matrix, meta)
    except Exception as exc:
        log(f"fixed-program diagnostic skipped: {exc}")
    best = summary.sort_values(["mean_matched_time_pearson", "knn_excess"], ascending=[False, True]).iloc[0] if len(summary) else None
    overall = pd.DataFrame([{
        "n_variants": len(summary), "baseline_matched_time_pearson": float(summary.loc[summary.variant=="baseline_all", "mean_matched_time_pearson"].iloc[0]),
        "best_variant": str(best.variant) if best is not None else "none",
        "best_matched_time_pearson": float(best.mean_matched_time_pearson) if best is not None else np.nan,
        "baseline_knn_excess": float(summary.loc[summary.variant=="baseline_all", "knn_excess"].iloc[0]),
        "best_knn_excess": float(best.knn_excess) if best is not None else np.nan,
        "repair_improves_matched_time": bool(best is not None and best.mean_matched_time_pearson > summary.loc[summary.variant=="baseline_all", "mean_matched_time_pearson"].iloc[0] + .05),
        "repair_reduces_dataset_knn_enrichment": bool(best is not None and best.knn_excess < summary.loc[summary.variant=="baseline_all", "knn_excess"].iloc[0] - .02),
        "stage3_readiness": False,
        "interpretation": "compare variants; no automatic acceptance of a repaired state"
    }])
    overall.to_csv(OUT / "07_STAGE2_9_20_SUMMARY.csv", index=False)
    log("complete")
    print("\nStage 2.9.20 variant comparison:")
    print(summary.to_string(index=False))
    print("\nStage 2.9.20 summary:")
    print(overall.to_string(index=False))
    return overall


if __name__ == "__main__": run()
