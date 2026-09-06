"""Stage 2.9.23: leakage-free one-step-ahead predictive-state validation.

Uses each dataset's next actually observed timepoint instead of forcing fixed
24/48/72 h horizons. Gene selection, scaling and PCA are fitted on training
datasets only. The primary question is whether the current state predicts the
next observed state better than persistence, nearest-training-time and linear
calendar-time baselines. No ODE/state-space model is fitted.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_23"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
VARIANTS = ("baseline_all", "low_dataset_excess_q80", "low_dataset_excess_q70", "time_dominant")
PRIMARY_VARIANT = "time_dominant"
N_PERM = 1000


def log(x):
    print(f"Stage 2.9.23: {x}", flush=True)


def corr(x, y, method="pearson"):
    x = np.asarray(x, float); y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 3 or np.std(x[good]) < 1e-12 or np.std(y[good]) < 1e-12:
        return np.nan
    return float(pd.Series(x[good]).corr(pd.Series(y[good]), method=method))


def load_space():
    from dynamics.validation import _load_common_space
    matrix, meta = _load_common_space()
    meta = meta[meta["dataset"].isin(TARGET)].copy()
    meta["matrix_column"] = meta["matrix_column"].astype(str)
    cols = [c for c in meta["matrix_column"] if c in matrix.columns]
    meta = meta[meta["matrix_column"].isin(cols)].copy()
    meta = meta.set_index("matrix_column").loc[cols].reset_index()
    return matrix.loc[:, cols], meta


def times(meta):
    if "time_hours" in meta.columns:
        return pd.to_numeric(meta["time_hours"], errors="coerce").to_numpy(float)
    from dynamics.validation import _time_hours_for_validation, _strip_dataset_prefix
    vals = []
    for i, r in meta.iterrows():
        idx = i if str(r["dataset"]) == "GSE28688" else None
        vals.append(_time_hours_for_validation(str(r["dataset"]), _strip_dataset_prefix(str(r["sample"])), idx))
    return np.asarray(vals, float)


def design(meta):
    ds = pd.get_dummies(meta["dataset"].astype(str), dtype=float).to_numpy()
    t = times(meta); ok = np.isfinite(t)
    tn = np.zeros_like(t)
    if ok.any() and np.std(t[ok]) > 1e-12:
        tn[ok] = (t[ok] - np.mean(t[ok])) / np.std(t[ok])
    Xd = np.column_stack([np.ones(len(meta)), ds[:, 1:]])
    Xt = np.column_stack([np.ones(len(meta)), tn])
    Xdt = np.column_stack([Xd, tn])
    return Xd, Xt, Xdt, t, ok


def r2(y, X):
    good = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if good.sum() < 6 or np.var(y[good]) <= 1e-12: return 0.0
    yy, xx = y[good], X[good]
    beta = np.linalg.lstsq(xx, yy, rcond=None)[0]
    den = np.sum((yy - yy.mean()) ** 2)
    return max(0.0, 1.0 - np.sum((yy - xx @ beta) ** 2) / den)


def feature_sets(matrix, meta):
    Xd, Xt, Xdt, t, ok = design(meta)
    rows=[]
    for gene, s in matrix.loc[:, meta.matrix_column].iterrows():
        y=s.to_numpy(float); good=np.isfinite(y)&ok
        if good.sum()<6 or np.var(y[good])<=1e-12: continue
        rd,rt,rdt=r2(y,Xd),r2(y,Xt),r2(y,Xdt)
        rows.append({"gene":str(gene),"r2_dataset":rd,"r2_time":rt,"r2_dataset_time":rdt,
                     "dataset_excess_r2":max(0.0,rdt-rt),"time_minus_dataset_r2":rt-rd})
    var=pd.DataFrame(rows); genes=set(matrix.index.astype(str))
    sets={"baseline_all":sorted(genes)}
    if len(var):
        v=var.set_index("gene")
        sets["low_dataset_excess_q80"]=sorted(set(v.index[v.dataset_excess_r2<=v.dataset_excess_r2.quantile(.80)])&genes)
        sets["low_dataset_excess_q70"]=sorted(set(v.index[v.dataset_excess_r2<=v.dataset_excess_r2.quantile(.70)])&genes)
        sets["time_dominant"]=sorted(set(v.index[(v.r2_time>=v.r2_dataset)&(v.r2_time>=v.r2_time.quantile(.50))])&genes)
    for k in list(sets):
        if len(sets[k])<100: sets[k]=sets["baseline_all"]
    return var,sets


def projection(matrix, meta, train_ds, test_ds, genes):
    tr=meta[meta.dataset.astype(str).isin(train_ds)].copy(); allm=meta[meta.dataset.astype(str).isin(list(train_ds)+[test_ds])].copy()
    tr["time"]=times(tr); allm["time"]=times(allm)
    tr=tr[np.isfinite(tr.time)].copy(); allm=allm[np.isfinite(allm.time)].copy()
    genes=[g for g in genes if g in matrix.index]
    if len(genes)<100 or tr.dataset.nunique()<2: return None
    trc,ac=tr.matrix_column.tolist(),allm.matrix_column.tolist()
    Xtr=np.array(matrix.loc[genes,trc].T.to_numpy(float),copy=True); Xa=np.array(matrix.loc[genes,ac].T.to_numpy(float),copy=True)
    med=np.nanmedian(Xtr,axis=0); med=np.where(np.isfinite(med),med,0.0)
    for X in (Xtr,Xa):
        bad=~np.isfinite(X)
        if bad.any(): X[bad]=np.take(med,np.where(bad)[1])
    scaler=StandardScaler().fit(Xtr); pca=PCA(n_components=min(10,Xtr.shape[0]-1,Xtr.shape[1]),random_state=0).fit(scaler.transform(Xtr))
    allm["state"]=pca.transform(scaler.transform(Xa))[:,0]
    traj={}
    for ds,g in allm.groupby("dataset"):
        rec=g.groupby("time",as_index=False).state.mean().sort_values("time")
        if len(rec)>=2: traj[str(ds)]=(rec.time.to_numpy(float),rec.state.to_numpy(float))
    return traj


def one_step(t,z):
    t=np.asarray(t,float); z=np.asarray(z,float); out=[]
    for i in range(len(t)-1):
        if np.isfinite(t[i]) and np.isfinite(z[i]) and np.isfinite(t[i+1]) and np.isfinite(z[i+1]) and t[i+1]>t[i]:
            out.append((float(t[i]),float(z[i]),float(t[i+1]),float(z[i+1]),float(t[i+1]-t[i])))
    return out


def training_pairs(traj, datasets):
    rows=[]
    for ds in datasets:
        if ds in traj: rows.extend(one_step(*traj[ds]))
    return rows


def nearest(train_traj, target):
    vals=[]
    for t,z in train_traj:
        if len(t): vals.append(float(z[int(np.argmin(np.abs(t-target)))]))
    return float(np.mean(vals)) if vals else np.nan


def linear(train_traj, target):
    tt=[]; zz=[]
    for t,z in train_traj: tt.extend(t.tolist()); zz.extend(z.tolist())
    if len(tt)<3 or np.std(tt)<1e-12: return np.nan
    return float(LinearRegression().fit(np.asarray(tt).reshape(-1,1),np.asarray(zz)).predict([[target]])[0])


def evaluate(matrix,meta,heldout,variant,genes):
    train_ds=[d for d in TARGET if d!=heldout]; traj=projection(matrix,meta,train_ds,heldout,genes)
    if not traj or heldout not in traj:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"n_genes":len(genes),"reason":"state_axis_fit_failed"},[]
    train_pairs=training_pairs(traj,train_ds); test_pairs=one_step(*traj[heldout])
    if len(train_pairs)<3 or len(test_pairs)<2:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"n_genes":len(genes),"n_training_transitions":len(train_pairs),"n_test_transitions":len(test_pairs),"reason":"insufficient_one_step_transitions"},[]
    X=np.asarray([r[1] for r in train_pairs]); Y=np.asarray([r[3] for r in train_pairs])
    model=LinearRegression().fit(X.reshape(-1,1),Y)
    pred=np.asarray([model.predict([[r[1]]])[0] for r in test_pairs]); true=np.asarray([r[3] for r in test_pairs])
    persistence=np.asarray([r[1] for r in test_pairs]); train_traj=[traj[d] for d in train_ds if d in traj]
    nearest=np.asarray([nearest(train_traj,r[2]) for r in test_pairs]); linear_pred=np.asarray([linear(train_traj,r[2]) for r in test_pairs])
    good=np.isfinite(pred)&np.isfinite(true)&np.isfinite(nearest)&np.isfinite(linear_pred)
    if good.sum()<2:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"n_genes":len(genes),"reason":"insufficient_finite_predictions"},[]
    rm=lambda a:float(np.sqrt(np.mean((a[good]-true[good])**2)))
    row={"status":"ok","held_out_dataset":heldout,"variant":variant,"n_genes":len(genes),"n_training_transitions":len(train_pairs),"n_test_transitions":int(good.sum()),
         "mean_test_delta_hours":float(np.mean([r[4] for r in test_pairs])),"future_state_rmse":rm(pred),"persistence_rmse":rm(persistence),"nearest_time_rmse":rm(nearest),"linear_time_rmse":rm(linear_pred),
         "improvement_vs_persistence":rm(persistence)-rm(pred),"improvement_vs_nearest_time":rm(nearest)-rm(pred),"improvement_vs_linear_time":rm(linear_pred)-rm(pred),
         "future_state_pearson":corr(pred[good],true[good]),"future_state_spearman":corr(pred[good],true[good],"spearman"),"transition_slope":float(model.coef_[0]),"transition_intercept":float(model.intercept_)}
    preds=[]
    for r,p in zip(test_pairs,pred): preds.append({"held_out_dataset":heldout,"variant":variant,"current_time_hours":r[0],"future_time_hours":r[2],"delta_hours":r[4],"current_state":r[1],"true_future_state":r[3],"predicted_future_state":p,"persistence_prediction":r[1],"nearest_training_time_prediction":nearest(train_traj,r[2]),"linear_time_prediction":linear(train_traj,r[2])})
    return row,preds


def permutation(matrix,meta,heldout,variant,genes,seed):
    train_ds=[d for d in TARGET if d!=heldout]; traj=projection(matrix,meta,train_ds,heldout,genes)
    if not traj or heldout not in traj: return []
    train_pairs=training_pairs(traj,train_ds); test_pairs=one_step(*traj[heldout])
    if len(train_pairs)<3 or len(test_pairs)<2: return []
    X=np.asarray([r[1] for r in train_pairs]); Y=np.asarray([r[3] for r in train_pairs]); model=LinearRegression().fit(X.reshape(-1,1),Y)
    pred=np.asarray([model.predict([[r[1]]])[0] for r in test_pairs]); true=np.asarray([r[3] for r in test_pairs]); obs=float(np.sqrt(np.mean((pred-true)**2)))
    rng=np.random.default_rng(seed); null=np.empty(N_PERM)
    for b in range(N_PERM): null[b]=np.sqrt(np.mean((pred-rng.permutation(true))**2))
    return [{"held_out_dataset":heldout,"variant":variant,"observed_future_state_rmse":obs,"permutation_p_rmse":float((1+np.sum(null<=obs))/(N_PERM+1)),"null_mean_rmse":float(null.mean()),"null_p05_rmse":float(np.quantile(null,.05)),"null_p95_rmse":float(np.quantile(null,.95))}]


def run():
    log("loading common gene space")
    matrix,meta=load_space(); meta["time_hours"]=times(meta)
    folds=[]; preds=[]; perms=[]; selections=[]
    for fold,heldout in enumerate(TARGET,1):
        log(f"fold {fold}/3: held out {heldout}")
        train_meta=meta[meta.dataset.astype(str)!=heldout].copy(); var,sets=feature_sets(matrix.loc[:,train_meta.matrix_column],train_meta)
        if len(var): var["held_out_dataset"]=heldout; selections.append(var)
        for v in VARIANTS:
            log(f"  {v}: {len(sets[v])} training-selected genes")
            r,pr=evaluate(matrix,meta,heldout,v,sets[v]); folds.append(r); preds.extend(pr); perms.extend(permutation(matrix,meta,heldout,v,sets[v],7000+fold))
    f=pd.DataFrame(folds); p=pd.DataFrame(preds); s=pd.concat(selections,ignore_index=True) if selections else pd.DataFrame()
    pc=["held_out_dataset","variant","observed_future_state_rmse","permutation_p_rmse","null_mean_rmse","null_p05_rmse","null_p95_rmse"]; pm=pd.DataFrame(perms,columns=pc)
    summaries=[]
    for v in VARIANTS:
        g=f[(f.variant==v)&(f.status=="ok")]
        row={"variant":v,"n_valid_folds":len(g)}
        for k in ("future_state_rmse","persistence_rmse","nearest_time_rmse","linear_time_rmse","improvement_vs_persistence","improvement_vs_nearest_time","improvement_vs_linear_time","future_state_pearson","future_state_spearman"):
            vals=g[k].replace([np.inf,-np.inf],np.nan).dropna() if len(g) and k in g else pd.Series(dtype=float); row["mean_"+k]=float(vals.mean()) if len(vals) else np.nan
        q=pm[pm["variant"]==v] if "variant" in pm.columns else pd.DataFrame(columns=pc); row["min_permutation_p_rmse"]=float(q.permutation_p_rmse.min()) if len(q) else np.nan; row["mean_permutation_p_rmse"]=float(q.permutation_p_rmse.mean()) if len(q) else np.nan
        summaries.append(row)
    sm=pd.DataFrame(summaries); pri=sm[sm.variant==PRIMARY_VARIANT].iloc[0]
    supported=bool(pri.n_valid_folds>=2 and pri.mean_improvement_vs_persistence>0 and pri.mean_improvement_vs_nearest_time>0 and pri.mean_improvement_vs_linear_time>0 and np.isfinite(pri.min_permutation_p_rmse) and pri.min_permutation_p_rmse<0.05)
    overall=pd.DataFrame([{"n_heldout_datasets":len(TARGET),"primary_variant":PRIMARY_VARIANT,"primary_n_valid_folds":int(pri.n_valid_folds),"primary_mean_improvement_vs_persistence":float(pri.mean_improvement_vs_persistence) if pd.notna(pri.mean_improvement_vs_persistence) else np.nan,"primary_mean_improvement_vs_nearest_time":float(pri.mean_improvement_vs_nearest_time) if pd.notna(pri.mean_improvement_vs_nearest_time) else np.nan,"primary_mean_improvement_vs_linear_time":float(pri.mean_improvement_vs_linear_time) if pd.notna(pri.mean_improvement_vs_linear_time) else np.nan,"primary_min_permutation_p_rmse":float(pri.min_permutation_p_rmse) if pd.notna(pri.min_permutation_p_rmse) else np.nan,"one_step_predictive_state_supported":supported,"stage3_readiness":False,"interpretation":"LODO one-step-ahead prediction using each dataset's next observed timepoint; training-only gene repair/scaling/PCA; no fixed-horizon interpolation and no ODE/state-space model"}])
    f.to_csv(OUT/"01_fold_results.csv",index=False); p.to_csv(OUT/"02_one_step_predictions.csv",index=False); s.to_csv(OUT/"03_training_gene_selection.csv",index=False); pm.to_csv(OUT/"04_permutation_null.csv",index=False); sm.to_csv(OUT/"05_variant_summary.csv",index=False); overall.to_csv(OUT/"06_stage2923_summary.csv",index=False)
    log("variant summary:"); print(sm.to_string(index=False)); log("overall:"); print(overall.to_string(index=False)); return overall

if __name__=="__main__": run()
