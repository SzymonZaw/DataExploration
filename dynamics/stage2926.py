"""Stage 2.9.26: leakage-free multivariate biological state validation.

Tests whether a fixed biologically anchored multi-dimensional state built from
Stage 2.9.14 programs predicts the next observed state better than simple
baselines. No ODE/state-space model is fitted.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"stage2_9_26"
OUT.mkdir(parents=True,exist_ok=True)
TARGET=["GSE67462","GSE28688","GSE297234"]
N_PERM=1000
PROGRAMS=(
"P01_PLURIPOTENCY","P02_PROLIFERATION","P03_EMT_MESENCHYMAL",
"P04_STRESS_RESPONSE","P05_GLYCOLYTIC_METABOLISM","P06_FGFR_PI3K_MAPK",
"P07_CHROMATIN_EPIGENETIC","P08_ECM_ADHESION")

def log(x): print(f"Stage 2.9.26: {x}",flush=True)

def rmse(a,b): return float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)))

def corr(a,b,method="pearson"):
    a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method=method))

def load():
    from dynamics.validation import _load_common_space
    m,meta=_load_common_space();meta=meta[meta.dataset.astype(str).isin(TARGET)].copy()
    meta["matrix_column"]=meta.matrix_column.astype(str)
    cols=[c for c in meta.matrix_column if c in m.columns]
    meta=meta[meta.matrix_column.isin(cols)].copy().set_index("matrix_column").loc[cols].reset_index()
    if "time_hours" not in meta.columns:
        from dynamics.validation import _time_hours_for_validation,_strip_dataset_prefix
        meta["time_hours"]=[_time_hours_for_validation(str(r.dataset),_strip_dataset_prefix(str(r.sample)),i if str(r.dataset)=="GSE28688" else None) for i,r in meta.iterrows()]
    meta["time_hours"]=pd.to_numeric(meta.time_hours,errors="coerce")
    return m.loc[:,cols],meta

def activity(matrix):
    from dynamics.stage2914 import PROGRAMS as P
    ranks=matrix.rank(axis=0,method="average",pct=True)
    lookup={str(g).upper():g for g in matrix.index};frames=[]
    for pid,spec in P.items():
        pos=[lookup[g] for g in spec["positive"] if g in lookup];neg=[lookup[g] for g in spec["negative"] if g in lookup]
        if len(pos)<3:continue
        s=ranks.loc[pos].mean(axis=0)
        if neg:s=s-ranks.loc[neg].mean(axis=0)
        frames.append(pd.DataFrame({"program_id":pid,"activity":s.values},index=matrix.columns))
    return pd.concat(frames).reset_index(names="matrix_column")

def trajectories(matrix,meta):
    a=activity(matrix);out={}
    for ds in TARGET:
        g=meta[meta.dataset.astype(str).eq(ds)&np.isfinite(meta.time_hours)].copy()
        g=g[g.matrix_column.isin(set(a.matrix_column.astype(str)))]
        if g.time_hours.nunique()<3:continue
        q=a[a.matrix_column.isin(g.matrix_column)].merge(g[["matrix_column","time_hours"]],on="matrix_column")
        q=q.groupby(["time_hours","program_id"],as_index=False).activity.mean()
        p=q.pivot(index="time_hours",columns="program_id",values="activity").sort_index()
        if all(pid in p.columns for pid in PROGRAMS):out[ds]=(p.index.to_numpy(float),p[list(PROGRAMS)].to_numpy(float))
    return out

def project(train_traj):
    """Fit a training-only program scaling; state is eight standardized programs."""
    vals=np.vstack([v for _,v in train_traj.values()])
    sc=StandardScaler().fit(vals)
    return sc

def transitions(t,z):
    return [(i,float(t[i]),float(t[i+1]),z[i].copy(),z[i+1].copy()) for i in range(len(t)-1) if t[i]<t[i+1] and np.all(np.isfinite(z[i])) and np.all(np.isfinite(z[i+1]))]

def run():
    log("loading common gene space and fixed biological programs")
    m,meta=load();trajs=trajectories(m,meta)
    log(f"constructed multivariate trajectories: {', '.join(sorted(trajs))}")
    fold_rows=[];pred_rows=[];perm_rows=[];selection=[]
    for fold,held in enumerate(TARGET,1):
        train={d:v for d,v in trajs.items() if d!=held}
        log(f"fold {fold}/3: held out {held}")
        if held not in trajs or len(train)<2:
            continue
        sc=project(train)
        train_t={d:(t,sc.transform(z)) for d,(t,z) in train.items()}
        test_t, test_z=trajs[held];test_z=sc.transform(test_z)
        test=transitions(test_t,test_z);trs=[]
        for d,(t,z) in train_t.items():trs.extend(transitions(t,z))
        if len(trs)<3 or len(test)<2:continue
        X=np.vstack([r[3] for r in trs]);Y=np.vstack([r[4] for r in trs])
        model=LinearRegression().fit(X,Y)
        pred=model.predict(np.vstack([r[3] for r in test]));true=np.vstack([r[4] for r in test]);cur=np.vstack([r[3] for r in test])
        # nearest training-time vector
        nearest=[];linear=[]
        flat_t=np.concatenate([t for t,_ in train_t.values()]);flat_z=np.vstack([z for _,z in train_t.values()])
        lt=LinearRegression().fit(flat_t.reshape(-1,1),flat_z)
        for r in test:
            j=int(np.argmin(np.abs(flat_t-r[2])));nearest.append(flat_z[j]);linear.append(lt.predict([[r[2]]])[0])
        nearest=np.vstack(nearest);linear=np.vstack(linear)
        row={"status":"ok","held_out_dataset":held,"n_training_datasets":len(train),"n_training_transitions":len(trs),"n_test_transitions":len(test),"n_state_dimensions":len(PROGRAMS),"state_model_rmse":rmse(pred,true),"persistence_rmse":rmse(cur,true),"nearest_time_rmse":rmse(nearest,true),"linear_time_rmse":rmse(linear,true),"improvement_vs_persistence":rmse(cur,true)-rmse(pred,true),"improvement_vs_nearest_time":rmse(nearest,true)-rmse(pred,true),"improvement_vs_linear_time":rmse(linear,true)-rmse(pred,true),"mean_dimension_pearson":float(np.nanmean([corr(pred[:,k],true[:,k]) for k in range(true.shape[1])])),"mean_dimension_spearman":float(np.nanmean([corr(pred[:,k],true[:,k],"spearman") for k in range(true.shape[1])]))}
        fold_rows.append(row)
        for r,a,b,n,l in zip(test,pred,true,nearest,linear):
            pred_rows.append({"held_out_dataset":held,"previous_time_hours":r[1],"future_time_hours":r[2],"state_model_rmse":rmse(a,b),"persistence_rmse":rmse(r[3],b),"nearest_time_rmse":rmse(n,b),"linear_time_rmse":rmse(l,b),"true_state":";".join(map(str,b)),"predicted_state":";".join(map(str,a))})
        selection.append({"held_out_dataset":held,"n_training_datasets":len(train),"n_state_dimensions":len(PROGRAMS),"state_scaling":"training-only StandardScaler on 8 fixed biological program activities"})
        rng=np.random.default_rng(26000+fold);obs=rmse(pred,true);null=np.empty(N_PERM)
        for b in range(N_PERM):null[b]=rmse(pred,rng.permutation(true))
        perm_rows.append({"held_out_dataset":held,"observed_state_model_rmse":obs,"permutation_p_rmse":float((1+np.sum(null<=obs))/(N_PERM+1)),"null_mean_rmse":float(null.mean()),"null_p05_rmse":float(np.quantile(null,.05)),"null_p95_rmse":float(np.quantile(null,.95))})
    f=pd.DataFrame(fold_rows);p=pd.DataFrame(pred_rows);s=pd.DataFrame(selection);pm=pd.DataFrame(perm_rows)
    sums=[]
    for col in ["state_model_rmse","persistence_rmse","nearest_time_rmse","linear_time_rmse","improvement_vs_persistence","improvement_vs_nearest_time","improvement_vs_linear_time","mean_dimension_pearson","mean_dimension_spearman"]:
        pass
    if len(f):
        summary=pd.DataFrame([{"variant":"fixed_8_program_state","n_valid_folds":len(f),"mean_state_model_rmse":f.state_model_rmse.mean(),"mean_persistence_rmse":f.persistence_rmse.mean(),"mean_nearest_time_rmse":f.nearest_time_rmse.mean(),"mean_linear_time_rmse":f.linear_time_rmse.mean(),"mean_improvement_vs_persistence":f.improvement_vs_persistence.mean(),"mean_improvement_vs_nearest_time":f.improvement_vs_nearest_time.mean(),"mean_improvement_vs_linear_time":f.improvement_vs_linear_time.mean(),"mean_dimension_pearson":f.mean_dimension_pearson.mean(),"mean_dimension_spearman":f.mean_dimension_spearman.mean(),"min_permutation_p_rmse":pm.permutation_p_rmse.min() if len(pm) else np.nan}])
    else:summary=pd.DataFrame()
    if len(summary):
        r=summary.iloc[0];supported=bool(r.n_valid_folds>=2 and r.mean_improvement_vs_persistence>0 and r.mean_improvement_vs_nearest_time>0 and r.mean_improvement_vs_linear_time>0 and np.isfinite(r.min_permutation_p_rmse) and r.min_permutation_p_rmse<0.05)
        overall=pd.DataFrame([{"n_heldout_datasets":len(TARGET),"primary_state":"8_fixed_biological_programs","n_valid_folds":int(r.n_valid_folds),"mean_improvement_vs_persistence":r.mean_improvement_vs_persistence,"mean_improvement_vs_nearest_time":r.mean_improvement_vs_nearest_time,"mean_improvement_vs_linear_time":r.mean_improvement_vs_linear_time,"min_permutation_p_rmse":r.min_permutation_p_rmse,"multivariate_predictive_state_supported":supported,"stage3_readiness":False,"interpretation":"LODO one-step prediction in an 8-dimensional biologically anchored state; training-only scaling; compared with persistence, nearest-time and linear-time; no ODE/state-space model"}])
    else:overall=pd.DataFrame([{"n_heldout_datasets":len(TARGET),"primary_state":"8_fixed_biological_programs","n_valid_folds":0,"multivariate_predictive_state_supported":False,"stage3_readiness":False}])
    s.to_csv(OUT/"01_training_state_scaling.csv",index=False);f.to_csv(OUT/"02_fold_results.csv",index=False);p.to_csv(OUT/"03_one_step_predictions.csv",index=False);pm.to_csv(OUT/"04_permutation_null.csv",index=False);summary.to_csv(OUT/"05_state_summary.csv",index=False);overall.to_csv(OUT/"06_stage2926_summary.csv",index=False)
    log("state summary:");print(summary.to_string(index=False));log("overall:");print(overall.to_string(index=False));return overall

if __name__=="__main__":run()
