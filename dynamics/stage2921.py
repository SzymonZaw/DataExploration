"""Stage 2.9.21: leakage-free validation of common-state repair.

For each trajectory dataset held out from GSE67462/GSE28688/GSE297234,
feature selection and PCA are fitted on training datasets only. Candidate
repair variants are pre-specified rather than selected on held-out data.
No ODE/state-space model is fitted here.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_21"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
VARIANTS = ("baseline_all", "low_dataset_excess_q80", "low_dataset_excess_q70", "time_dominant")
N_PERM = 500


def log(x):
    print(f"Stage 2.9.21: {x}", flush=True)


def load_space():
    from dynamics.validation import _load_common_space
    matrix, meta = _load_common_space()
    meta = meta[meta["dataset"].isin(TARGET)].copy()
    cols = [c for c in meta["matrix_column"].astype(str) if c in matrix.columns]
    meta = meta[meta["matrix_column"].isin(cols)].copy()
    meta = meta.set_index("matrix_column").loc[cols].reset_index()
    return matrix.loc[:, cols], meta


def _time(meta):
    if "time_hours" in meta.columns:
        t = pd.to_numeric(meta["time_hours"], errors="coerce").to_numpy(float)
    else:
        from dynamics.validation import _time_hours_for_validation, _strip_dataset_prefix
        vals=[]
        for i, r in meta.iterrows():
            vals.append(_time_hours_for_validation(str(r["dataset"]), _strip_dataset_prefix(str(r["sample"])), i if str(r["dataset"]) == "GSE28688" else None))
        t=np.asarray(vals,float)
    return t


def _design(meta):
    ds = pd.get_dummies(meta["dataset"].astype(str), dtype=float).to_numpy()
    t = _time(meta); ok=np.isfinite(t)
    tn=np.zeros_like(t)
    if ok.any() and np.ptp(t[ok]) > 0:
        tn[ok]=(t[ok]-np.mean(t[ok]))/np.std(t[ok])
    Xd=np.column_stack([np.ones(len(meta)), ds[:,1:]])
    Xt=np.column_stack([np.ones(len(meta)), tn])
    Xdt=np.column_stack([Xd, tn])
    return Xd, Xt, Xdt, t, ok


def _r2(y, X):
    good=np.isfinite(y) & np.all(np.isfinite(X),axis=1)
    if good.sum()<6 or np.var(y[good])<=1e-12: return 0.0
    xx=X[good]; yy=y[good]
    beta=np.linalg.lstsq(xx,yy,rcond=None)[0]
    den=np.sum((yy-yy.mean())**2)
    return max(0.0,1.0-np.sum((yy-xx@beta)**2)/den)


def training_feature_sets(matrix, meta):
    """Compute repair filters using training datasets only."""
    Xd,Xt,Xdt,t,ok=_design(meta)
    rows=[]
    for gene, s in matrix.loc[:,meta["matrix_column"]].iterrows():
        y=s.to_numpy(float); good=np.isfinite(y)&ok
        if good.sum()<6 or np.var(y[good])<=1e-12: continue
        rd=_r2(y,Xd); rt=_r2(y,Xt); rdt=_r2(y,Xdt)
        rows.append({"gene":str(gene),"r2_dataset":rd,"r2_time":rt,"r2_dataset_time":rdt,
                     "dataset_excess_r2":max(0.0,rdt-rt),"time_minus_dataset_r2":rt-rd})
    var=pd.DataFrame(rows)
    genes=set(matrix.index.astype(str))
    sets={"baseline_all":sorted(genes)}
    if len(var):
        v=var.set_index("gene")
        sets["low_dataset_excess_q80"]=sorted(set(v.index[v.dataset_excess_r2<=v.dataset_excess_r2.quantile(.80)]) & genes)
        sets["low_dataset_excess_q70"]=sorted(set(v.index[v.dataset_excess_r2<=v.dataset_excess_r2.quantile(.70)]) & genes)
        sets["time_dominant"]=sorted(set(v.index[(v.r2_time>=v.r2_dataset)&(v.r2_time>=v.r2_time.quantile(.50))]) & genes)
    for k in list(sets):
        if len(sets[k])<100: sets[k]=sets["baseline_all"]
    return var,sets


def _impute_scale(train, test):
    tr=np.array(train, dtype=float, copy=True); te=np.array(test, dtype=float, copy=True)
    med=np.nanmedian(tr,axis=0)
    med=np.where(np.isfinite(med),med,0.0)
    for X in (tr,te):
        bad=~np.isfinite(X)
        if bad.any(): X[bad]=np.take(med,np.where(bad)[1])
    scaler=StandardScaler().fit(tr)
    return scaler.transform(tr),scaler.transform(te),scaler


def _aggregate(matrix, meta, genes):
    m=matrix.loc[genes,meta["matrix_column"]]
    rows=[]
    for ds,g in meta.groupby("dataset"):
        for t,gt in g.groupby("time_hours"):
            if not np.isfinite(float(t)): continue
            cols=gt["matrix_column"].tolist()
            rows.append({"dataset":str(ds),"time_hours":float(t),"vector":m[cols].mean(axis=1).to_numpy(float)})
    return rows


def _safe_corr(x,y,method="pearson"):
    good=np.isfinite(x)&np.isfinite(y)
    if good.sum()<3 or np.std(x[good])==0 or np.std(y[good])==0:return np.nan
    if method=="spearman":
        return float(pd.Series(x[good]).corr(pd.Series(y[good]),method="spearman"))
    return float(np.corrcoef(x[good],y[good])[0,1])


def _nearest_time_baseline(train_times, train_progress, test_times):
    pred=[]
    for t in test_times:
        j=int(np.argmin(np.abs(train_times-t)));pred.append(train_progress[j])
    return np.asarray(pred,float)


def evaluate_fold(matrix, meta, heldout, variant, genes, rng):
    train_meta=meta[meta["dataset"]!=heldout].copy(); test_meta=meta[meta["dataset"]==heldout].copy()
    train_meta=train_meta[np.isfinite(pd.to_numeric(train_meta["time_hours"],errors="coerce"))].copy()
    test_meta=test_meta[np.isfinite(pd.to_numeric(test_meta["time_hours"],errors="coerce"))].copy()
    if train_meta["dataset"].nunique()<2 or test_meta["time_hours"].nunique()<2:
        return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"reason":"insufficient_training_or_test_trajectory"}, []
    genes=[g for g in genes if g in matrix.index]
    if len(genes)<100:return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"reason":"too_few_training_selected_genes"},[]
    tr_cols=train_meta["matrix_column"].tolist(); te_cols=test_meta["matrix_column"].tolist()
    tr,test,_=_impute_scale(matrix.loc[genes,tr_cols].T.to_numpy(float),matrix.loc[genes,te_cols].T.to_numpy(float))
    ncomp=min(10,tr.shape[0]-1,tr.shape[1])
    if ncomp<1:return {"status":"skipped","held_out_dataset":heldout,"variant":variant,"reason":"no_pca_component"},[]
    pca=PCA(n_components=ncomp,random_state=0).fit(tr); ztr=pca.transform(tr)[:,0]; zte=pca.transform(test)[:,0]
    tt=pd.to_numeric(train_meta["time_hours"],errors="coerce").to_numpy(float); tv=pd.to_numeric(test_meta["time_hours"],errors="coerce").to_numpy(float)
    scale=max(np.ptp(tt),1e-12); t0=np.min(tt); trn=(tt-t0)/scale; ten=(tv-t0)/scale
    orient=1.0 if _safe_corr(ztr,trn,"pearson")>=0 else -1.0; ztr*=orient;zte*=orient
    model=LinearRegression().fit(ztr.reshape(-1,1),trn)
    pred=np.clip(model.predict(zte.reshape(-1,1)),0,1)
    nearest=_nearest_time_baseline(tt,trn,tv)
    rmse=float(np.sqrt(np.mean((pred-ten)**2))); nrmse=float(np.sqrt(np.mean((nearest-ten)**2)))
    spearman=_safe_corr(pred,ten,"spearman"); pearson=_safe_corr(pred,ten,"pearson")
    monot=float(np.mean(np.diff(pred[np.argsort(tv)])>=0)) if len(pred)>1 else np.nan
    rows=[]
    for t,p,n in zip(tv,pred,nearest): rows.append({"held_out_dataset":heldout,"variant":variant,"time_hours":float(t),"predicted_normalized_time":float(p),"actual_normalized_time":float((t-t0)/scale),"nearest_time_prediction":float(n)})
    return {"status":"ok","held_out_dataset":heldout,"variant":variant,"n_genes":len(genes),"n_train_samples":len(tr_cols),"n_test_samples":len(te_cols),"n_train_datasets":train_meta["dataset"].nunique(),"n_test_timepoints":test_meta["time_hours"].nunique(),"explained_variance_pc1":float(pca.explained_variance_ratio_[0]),"train_pc1_time_pearson":float(_safe_corr(ztr,trn,"pearson")),"test_progress_spearman":spearman,"test_progress_pearson":pearson,"test_progress_rmse":rmse,"nearest_time_rmse":nrmse,"rmse_improvement":nrmse-rmse,"monotonicity":monot,"axis_orientation":orient},rows


def permutation_fold(matrix, meta, heldout, variant, genes, seed):
    """Permutation null for held-out state/time association; repair and PCA stay training-only."""
    train_meta=meta[meta["dataset"]!=heldout].copy(); test_meta=meta[meta["dataset"]==heldout].copy()
    train_meta=train_meta[np.isfinite(pd.to_numeric(train_meta["time_hours"],errors="coerce"))].copy(); test_meta=test_meta[np.isfinite(pd.to_numeric(test_meta["time_hours"],errors="coerce"))].copy()
    if train_meta["dataset"].nunique()<2 or test_meta["time_hours"].nunique()<2:return []
    genes=[g for g in genes if g in matrix.index]
    tr_cols=train_meta["matrix_column"].tolist(); te_cols=test_meta["matrix_column"].tolist()
    tr,test,_=_impute_scale(matrix.loc[genes,tr_cols].T.to_numpy(float),matrix.loc[genes,te_cols].T.to_numpy(float))
    pca=PCA(n_components=min(10,tr.shape[0]-1,tr.shape[1]),random_state=0).fit(tr); ztr=pca.transform(tr)[:,0]; zte=pca.transform(test)[:,0]
    tt=pd.to_numeric(train_meta["time_hours"],errors="coerce").to_numpy(float); tv=pd.to_numeric(test_meta["time_hours"],errors="coerce").to_numpy(float); scale=max(np.ptp(tt),1e-12);tn=(tt-np.min(tt))/scale;te=(tv-np.min(tt))/scale
    if _safe_corr(ztr,tn)>=0:ztr=-ztr* -1
    else:ztr=-ztr* -1
    model=LinearRegression().fit(ztr.reshape(-1,1),tn); pred=model.predict(zte.reshape(-1,1))
    rng=np.random.default_rng(seed); out=[]
    for p in range(N_PERM):
        perm=rng.permutation(te);out.append({"held_out_dataset":heldout,"variant":variant,"permutation":p,"spearman":_safe_corr(pred,perm,"spearman"),"pearson":_safe_corr(pred,perm,"pearson")})
    return out


def run():
    log("loading common gene space");matrix,meta=load_space();meta["time_hours"]=_time(meta)
    fold_results=[];trajectory_rows=[];var_rows=[];perm_rows=[]
    for fold,heldout in enumerate(TARGET,1):
        train_meta=meta[meta["dataset"]!=heldout].copy(); log(f"fold {fold}/3: held out {heldout}")
        var,sets=training_feature_sets(matrix.loc[:,train_meta["matrix_column"]],train_meta)
        var["held_out_dataset"]=heldout;var_rows.append(var)
        for variant in VARIANTS:
            log(f"  {variant}: {len(sets[variant])} training-selected genes")
            r,rows=evaluate_fold(matrix,meta,heldout,variant,sets[variant],np.random.default_rng(fold))
            fold_results.append(r);trajectory_rows.extend(rows)
        if fold==1: log("  permutation nulls: 500 per fold/variant")
        for variant in VARIANTS:
            perm_rows.extend(permutation_fold(matrix,meta,heldout,variant,sets[variant],1000+fold))
    folds=pd.DataFrame(fold_results);traj=pd.DataFrame(trajectory_rows);perms=pd.DataFrame(perm_rows);vard=pd.concat(var_rows,ignore_index=True) if var_rows else pd.DataFrame()
    folds.to_csv(OUT/"01_lodo_variant_results.csv",index=False);traj.to_csv(OUT/"02_lodo_state_predictions.csv",index=False);vard.to_csv(OUT/"03_training_gene_selection_by_fold.csv",index=False);perms.to_csv(OUT/"04_time_permutation_null.csv",index=False)
    ok=folds[folds.status=="ok"].copy()
    summary=[]
    for v in VARIANTS:
        g=ok[ok.variant==v]
        if len(g):
            obs=float(g.test_progress_spearman.mean()); null=perms[perms.variant==v]["spearman"].dropna().to_numpy(float)
            p=float((1+np.sum(np.abs(null)>=abs(obs)))/(1+len(null))) if len(null) else np.nan
            summary.append({"variant":v,"n_valid_folds":len(g),"mean_test_progress_spearman":obs,"mean_test_progress_pearson":float(g.test_progress_pearson.mean()),"mean_test_rmse":float(g.test_progress_rmse.mean()),"mean_nearest_time_rmse":float(g.nearest_time_rmse.mean()),"mean_rmse_improvement":float(g.rmse_improvement.mean()),"mean_monotonicity":float(g.monotonicity.mean()),"permutation_p_spearman":p})
    s=pd.DataFrame(summary);s.to_csv(OUT/"05_variant_summary.csv",index=False)
    base=s[s.variant=="baseline_all"].iloc[0] if len(s[s.variant=="baseline_all"]) else None
    repair=s[s.variant!="baseline_all"].copy();best=repair.sort_values("mean_rmse_improvement",ascending=False).iloc[0] if len(repair) else None
    gate=bool(best is not None and best.mean_rmse_improvement>0 and best.permutation_p_spearman<0.05 and best.n_valid_folds>=2)
    overall=pd.DataFrame([{"n_heldout_datasets":len(TARGET),"n_valid_folds_baseline":int(base.n_valid_folds) if base is not None else 0,"best_repair_variant":str(best.variant) if best is not None else "","best_repair_mean_rmse_improvement":float(best.mean_rmse_improvement) if best is not None else np.nan,"best_repair_permutation_p_spearman":float(best.permutation_p_spearman) if best is not None else np.nan,"leakage_free_repair_supported":gate,"stage3_readiness":False,"interpretation":"held-out validation only; repair variant is not selected from held-out performance; no ODE/state-space model"}]);overall.to_csv(OUT/"06_STAGE2_9_21_SUMMARY.csv",index=False)
    log("complete");print("\nStage 2.9.21 variant summary:");print(s.to_string(index=False));print("\nStage 2.9.21 summary:");print(overall.to_string(index=False));return overall

if __name__=="__main__":run()
