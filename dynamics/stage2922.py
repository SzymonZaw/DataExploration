"""Stage 2.9.22: leakage-free predictive-state validation.

Tests whether a repaired common-state coordinate at time t predicts the
future state at t+delta across datasets, rather than merely reconstructing
current time. All gene filtering, scaling and PCA are fitted on training
datasets only. Held-out datasets are never used to choose a repair variant.

This is a pre-Stage-3 diagnostic: no ODE/state-space model is fitted.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_22"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
VARIANTS = ("baseline_all", "low_dataset_excess_q80", "low_dataset_excess_q70", "time_dominant")
HORIZONS = (24.0, 48.0, 72.0)
N_PERM = 500
PRIMARY_VARIANT = "time_dominant"


def log(x):
    print(f"Stage 2.9.22: {x}", flush=True)


def _safe_corr(x, y, method="pearson"):
    x = np.asarray(x, float); y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 3 or np.std(x[good]) < 1e-12 or np.std(y[good]) < 1e-12:
        return np.nan
    return float(pd.Series(x[good]).corr(pd.Series(y[good]), method=method))


def load_space():
    from dynamics.validation import _load_common_space
    matrix, meta = _load_common_space()
    meta = meta[meta["dataset"].isin(TARGET)].copy()
    cols = [str(c) for c in meta["matrix_column"] if str(c) in matrix.columns]
    meta = meta[meta["matrix_column"].astype(str).isin(cols)].copy()
    meta["matrix_column"] = meta["matrix_column"].astype(str)
    meta = meta.set_index("matrix_column").loc[cols].reset_index()
    return matrix.loc[:, cols], meta


def _times(meta):
    if "time_hours" in meta.columns:
        return pd.to_numeric(meta["time_hours"], errors="coerce").to_numpy(float)
    from dynamics.validation import _time_hours_for_validation, _strip_dataset_prefix
    vals = []
    for i, r in meta.iterrows():
        idx = i if str(r["dataset"]) == "GSE28688" else None
        vals.append(_time_hours_for_validation(str(r["dataset"]), _strip_dataset_prefix(str(r["sample"])), idx))
    return np.asarray(vals, float)


def _design(meta):
    ds = pd.get_dummies(meta["dataset"].astype(str), dtype=float).to_numpy()
    t = _times(meta); ok = np.isfinite(t)
    tn = np.zeros_like(t)
    if ok.any() and np.std(t[ok]) > 1e-12:
        tn[ok] = (t[ok] - np.mean(t[ok])) / np.std(t[ok])
    Xd = np.column_stack([np.ones(len(meta)), ds[:, 1:]])
    Xt = np.column_stack([np.ones(len(meta)), tn])
    Xdt = np.column_stack([Xd, tn])
    return Xd, Xt, Xdt, t, ok


def _r2(y, X):
    good = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if good.sum() < 6 or np.var(y[good]) <= 1e-12:
        return 0.0
    yy = y[good]; xx = X[good]
    beta = np.linalg.lstsq(xx, yy, rcond=None)[0]
    den = np.sum((yy - yy.mean()) ** 2)
    return max(0.0, 1.0 - np.sum((yy - xx @ beta) ** 2) / den)


def training_feature_sets(matrix, meta):
    """Fit all gene filters on training datasets only."""
    Xd, Xt, Xdt, t, ok = _design(meta)
    rows = []
    for gene, s in matrix.loc[:, meta["matrix_column"]].iterrows():
        y = s.to_numpy(float); good = np.isfinite(y) & ok
        if good.sum() < 6 or np.var(y[good]) <= 1e-12:
            continue
        rd = _r2(y, Xd); rt = _r2(y, Xt); rdt = _r2(y, Xdt)
        rows.append({"gene": str(gene), "r2_dataset": rd, "r2_time": rt,
                     "r2_dataset_time": rdt, "dataset_excess_r2": max(0.0, rdt - rt),
                     "time_minus_dataset_r2": rt - rd})
    var = pd.DataFrame(rows)
    genes = set(matrix.index.astype(str))
    sets = {"baseline_all": sorted(genes)}
    if len(var):
        v = var.set_index("gene")
        sets["low_dataset_excess_q80"] = sorted(set(v.index[v.dataset_excess_r2 <= v.dataset_excess_r2.quantile(.80)]) & genes)
        sets["low_dataset_excess_q70"] = sorted(set(v.index[v.dataset_excess_r2 <= v.dataset_excess_r2.quantile(.70)]) & genes)
        sets["time_dominant"] = sorted(set(v.index[(v.r2_time >= v.r2_dataset) & (v.r2_time >= v.r2_time.quantile(.50))]) & genes)
    for k in list(sets):
        if len(sets[k]) < 100:
            sets[k] = sets["baseline_all"]
    return var, sets


def _impute_scale(train, test):
    tr = np.array(train, float, copy=True); te = np.array(test, float, copy=True)
    med = np.nanmedian(tr, axis=0); med = np.where(np.isfinite(med), med, 0.0)
    for X in (tr, te):
        bad = ~np.isfinite(X)
        if bad.any():
            X[bad] = np.take(med, np.where(bad)[1])
    scaler = StandardScaler().fit(tr)
    return scaler.transform(tr), scaler.transform(te)


def _state_trajectories(matrix, meta, train_datasets, test_dataset, genes):
    """Fit state axis on training samples and return per-dataset time/state means."""
    trm = meta[meta.dataset.astype(str).isin(train_datasets)].copy()
    allm = meta[meta.dataset.astype(str).isin(list(train_datasets) + [test_dataset])].copy()
    trm["time"] = _times(trm); allm["time"] = _times(allm)
    trm = trm[np.isfinite(trm.time)].copy(); allm = allm[np.isfinite(allm.time)].copy()
    genes = [g for g in genes if g in matrix.index]
    if len(genes) < 100 or trm.dataset.nunique() < 2:
        return None
    tr_cols = trm.matrix_column.tolist()
    Xtr, _ = _impute_scale(matrix.loc[genes, tr_cols].T.to_numpy(float), matrix.loc[genes, tr_cols].T.to_numpy(float))
    ncomp = min(10, Xtr.shape[0] - 1, Xtr.shape[1])
    if ncomp < 1:
        return None
    pca = PCA(n_components=ncomp, random_state=0).fit(Xtr)
    # Project every sample with training-fitted preprocessing.
    all_cols = allm.matrix_column.tolist()
    _, Xall = _impute_scale(matrix.loc[genes, tr_cols].T.to_numpy(float), matrix.loc[genes, all_cols].T.to_numpy(float))
    z = pca.transform(Xall)[:, 0]
    allm["state"] = z
    out = {}
    for ds, g in allm.groupby("dataset"):
        rec = g.groupby("time", as_index=False).state.mean().sort_values("time")
        if len(rec) >= 2:
            out[str(ds)] = (rec.time.to_numpy(float), rec.state.to_numpy(float))
    return out, pca


def _transition_pairs(trajectories, datasets, delta):
    X=[]; Y=[]
    for ds in datasets:
        if ds not in trajectories: continue
        t,z = trajectories[ds]
        for i, ti in enumerate(t):
            j = np.where(np.isclose(t, ti + delta, atol=1e-6))[0]
            if len(j):
                X.append(z[i]); Y.append(z[int(j[0])])
    return np.asarray(X, float), np.asarray(Y, float)


def _nearest_time_prediction(train_traj, target_time):
    vals=[]
    for t,z in train_traj:
        if len(t):
            j=int(np.argmin(np.abs(t-target_time))); vals.append(z[j])
    return float(np.mean(vals)) if vals else np.nan


def _linear_time_model(train_traj, target_time):
    tt=[]; zz=[]
    for t,z in train_traj:
        tt.extend(t.tolist()); zz.extend(z.tolist())
    if len(tt)<3 or np.std(tt)<1e-12:return np.nan
    model=LinearRegression().fit(np.asarray(tt).reshape(-1,1),np.asarray(zz))
    return float(model.predict([[target_time]])[0])


def evaluate_fold(matrix, meta, heldout, variant, genes):
    train_datasets=[d for d in TARGET if d != heldout]
    test_meta=meta[meta.dataset.astype(str)==heldout].copy(); test_meta["time"]=_times(test_meta)
    test_meta=test_meta[np.isfinite(test_meta.time)].copy()
    if len(test_meta)<3 or test_meta.time.nunique()<3:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"reason":"insufficient_test_trajectory"}, []
    fitted=_state_trajectories(matrix, meta, train_datasets, heldout, genes)
    if fitted is None:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"reason":"state_axis_fit_failed"}, []
    trajectories,pca=fitted
    if heldout not in trajectories:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"reason":"no_test_state_trajectory"}, []
    train_traj=[trajectories[d] for d in train_datasets if d in trajectories]
    results=[]; pred_rows=[]
    for delta in HORIZONS:
        X,Y=_transition_pairs(trajectories, train_datasets, delta)
        test_t,test_z=trajectories[heldout]
        pairs=[]
        for i,ti in enumerate(test_t):
            j=np.where(np.isclose(test_t,ti+delta,atol=1e-6))[0]
            if len(j):pairs.append((ti,float(test_z[i]),float(test_z[int(j[0])]),float(test_t[int(j[0])])))
        if len(pairs)<2 or len(X)<3:
            continue
        model=LinearRegression().fit(X.reshape(-1,1),Y)
        pred=[]; true=[]; persistence=[]; nearest=[]; linear=[]
        for ti,zi,y,tf in pairs:
            pred.append(float(model.predict([[zi]])[0])); true.append(y); persistence.append(zi)
            nearest.append(_nearest_time_prediction(train_traj,tf)); linear.append(_linear_time_model(train_traj,tf))
            pred_rows.append({"held_out_dataset":heldout,"variant":variant,"horizon_hours":delta,"current_time_hours":ti,"future_time_hours":tf,"current_state":zi,"true_future_state":y,"predicted_future_state":pred[-1],"persistence_prediction":persistence[-1],"nearest_training_time_prediction":nearest[-1],"linear_time_prediction":linear[-1]})
        pred=np.asarray(pred); true=np.asarray(true); persistence=np.asarray(persistence); nearest=np.asarray(nearest); linear=np.asarray(linear)
        good=np.isfinite(pred)&np.isfinite(true)&np.isfinite(nearest)&np.isfinite(linear)
        if good.sum()<2: continue
        rmse=float(np.sqrt(np.mean((pred[good]-true[good])**2)))
        prmse=float(np.sqrt(np.mean((persistence[good]-true[good])**2)))
        nrmse=float(np.sqrt(np.mean((nearest[good]-true[good])**2)))
        lrmse=float(np.sqrt(np.mean((linear[good]-true[good])**2)))
        results.append({"status":"ok","held_out_dataset":heldout,"variant":variant,"horizon_hours":delta,"n_training_transition_pairs":len(X),"n_test_pairs":int(good.sum()),"transition_slope":float(model.coef_[0]),"transition_intercept":float(model.intercept_),"future_state_rmse":rmse,"persistence_rmse":prmse,"nearest_time_rmse":nrmse,"linear_time_rmse":lrmse,"improvement_vs_persistence":prmse-rmse,"improvement_vs_nearest_time":nrmse-rmse,"improvement_vs_linear_time":lrmse-rmse,"future_state_pearson":_safe_corr(pred[good],true[good]),"future_state_spearman":_safe_corr(pred[good],true[good],"spearman")})
    if not results:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"reason":"no_supported_horizons"}, pred_rows
    agg={"status":"ok","held_out_dataset":heldout,"variant":variant,"n_genes":len([g for g in genes if g in matrix.index]),"n_horizons":len(results),"n_test_pairs":int(sum(r["n_test_pairs"] for r in results))}
    for key in ("future_state_rmse","persistence_rmse","nearest_time_rmse","linear_time_rmse","improvement_vs_persistence","improvement_vs_nearest_time","improvement_vs_linear_time","future_state_pearson","future_state_spearman"):
        agg["mean_"+key]=float(np.nanmean([r[key] for r in results]))
    return agg,pred_rows


def permutation_fold(matrix, meta, heldout, variant, genes, seed):
    train_datasets=[d for d in TARGET if d != heldout]
    fitted=_state_trajectories(matrix,meta,train_datasets,heldout,genes)
    if fitted is None or heldout not in fitted:return []
    trajectories,_=fitted; test_t,test_z=trajectories[heldout]; train_traj=[trajectories[d] for d in train_datasets if d in trajectories]
    rng=np.random.default_rng(seed); rows=[]
    for delta in HORIZONS:
        X,Y=_transition_pairs(trajectories,train_datasets,delta)
        pairs=[]
        for i,ti in enumerate(test_t):
            j=np.where(np.isclose(test_t,ti+delta,atol=1e-6))[0]
            if len(j):pairs.append((float(test_z[i]),float(test_z[int(j[0])]),float(test_t[int(j[0])])))
        if len(pairs)<2 or len(X)<3:continue
        model=LinearRegression().fit(X.reshape(-1,1),Y)
        pred=np.asarray([model.predict([[zi]])[0] for zi,_,_ in pairs])
        true=np.asarray([y for _,y,_ in pairs])
        observed_rmse=float(np.sqrt(np.mean((pred-true)**2)))
        null=[]
        for b in range(N_PERM):
            perm=rng.permutation(true); null.append(float(np.sqrt(np.mean((pred-perm)**2))))
        null=np.asarray(null); p=float((1+np.sum(null<=observed_rmse))/(1+len(null)))
        rows.append({"held_out_dataset":heldout,"variant":variant,"horizon_hours":delta,"observed_future_state_rmse":observed_rmse,"permutation_p_rmse":p,"null_mean_rmse":float(null.mean()),"null_p05_rmse":float(np.quantile(null,.05)),"null_p95_rmse":float(np.quantile(null,.95))})
    return rows


def run():
    log("loading common gene space")
    matrix,meta=load_space(); meta["time_hours"]=_times(meta)
    fold_rows=[]; pred_rows=[]; selection_rows=[]; perm_rows=[]
    for fold,heldout in enumerate(TARGET,1):
        log(f"fold {fold}/3: held out {heldout}")
        train_meta=meta[meta.dataset.astype(str)!=heldout].copy()
        var,sets=training_feature_sets(matrix.loc[:,train_meta.matrix_column],train_meta)
        var["held_out_dataset"]=heldout; selection_rows.append(var)
        for variant in VARIANTS:
            log(f"  {variant}: {len(sets[variant])} training-selected genes")
            r,rows=evaluate_fold(matrix,meta,heldout,variant,sets[variant]); fold_rows.append(r); pred_rows.extend(rows)
            perm_rows.extend(permutation_fold(matrix,meta,heldout,variant,sets[variant],1000+fold))
    folds=pd.DataFrame(fold_rows); preds=pd.DataFrame(pred_rows); perms=pd.DataFrame(perm_rows); sel=pd.concat(selection_rows,ignore_index=True) if selection_rows else pd.DataFrame()
    folds.to_csv(OUT/"01_lodo_future_state_results.csv",index=False); preds.to_csv(OUT/"02_future_state_predictions.csv",index=False); sel.to_csv(OUT/"03_training_gene_selection_by_fold.csv",index=False); perms.to_csv(OUT/"04_future_state_permutation_null.csv",index=False)
    ok=folds[folds.status=="ok"].copy(); summaries=[]
    for v in VARIANTS:
        g=ok[ok.variant==v]
        if len(g):
            pm=perms[perms.variant==v]
            summaries.append({"variant":v,"n_valid_folds":len(g),"mean_future_state_rmse":float(g.mean_future_state_rmse.mean()),"mean_persistence_rmse":float(g.mean_persistence_rmse.mean()),"mean_nearest_time_rmse":float(g.mean_nearest_time_rmse.mean()),"mean_linear_time_rmse":float(g.mean_linear_time_rmse.mean()),"mean_improvement_vs_persistence":float(g.mean_improvement_vs_persistence.mean()),"mean_improvement_vs_nearest_time":float(g.mean_improvement_vs_nearest_time.mean()),"mean_improvement_vs_linear_time":float(g.mean_improvement_vs_linear_time.mean()),"mean_future_state_pearson":float(g.mean_future_state_pearson.mean()),"mean_future_state_spearman":float(g.mean_future_state_spearman.mean()),"min_permutation_p_rmse":float(pm.permutation_p_rmse.min()) if len(pm) else np.nan,"mean_permutation_p_rmse":float(pm.permutation_p_rmse.mean()) if len(pm) else np.nan})
    summary=pd.DataFrame(summaries); summary.to_csv(OUT/"05_variant_summary.csv",index=False)
    primary=summary[summary.variant==PRIMARY_VARIANT]
    gate=False
    if len(primary):
        p=primary.iloc[0]
        # Primary gate is intentionally pre-specified. It must beat persistence,
        # nearest-time and linear-time baselines and have at least two LODO folds.
        gate=bool(p.n_valid_folds>=2 and p.mean_improvement_vs_persistence>0 and p.mean_improvement_vs_nearest_time>0 and p.mean_improvement_vs_linear_time>0 and np.isfinite(p.min_permutation_p_rmse) and p.min_permutation_p_rmse<0.05)
    overall=pd.DataFrame([{"n_heldout_datasets":len(TARGET),"primary_variant":PRIMARY_VARIANT,"primary_n_valid_folds":int(primary.iloc[0].n_valid_folds) if len(primary) else 0,"primary_mean_improvement_vs_persistence":float(primary.iloc[0].mean_improvement_vs_persistence) if len(primary) else np.nan,"primary_mean_improvement_vs_nearest_time":float(primary.iloc[0].mean_improvement_vs_nearest_time) if len(primary) else np.nan,"primary_mean_improvement_vs_linear_time":float(primary.iloc[0].mean_improvement_vs_linear_time) if len(primary) else np.nan,"primary_min_permutation_p_rmse":float(primary.iloc[0].min_permutation_p_rmse) if len(primary) else np.nan,"predictive_state_supported":gate,"stage3_readiness":False,"interpretation":"LODO future-state prediction; repair/PCA are training-only; no ODE/state-space model"}])
    overall.to_csv(OUT/"06_STAGE2_9_22_SUMMARY.csv",index=False)
    log("complete")
    print("\nStage 2.9.22 variant summary:"); print(summary.to_string(index=False)); print("\nStage 2.9.22 summary:"); print(overall.to_string(index=False))
    return overall


if __name__=="__main__": run()
