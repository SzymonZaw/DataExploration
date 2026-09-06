"""Stage 2.9.27: leakage-free invariant biological coordinate discovery.

Discovers low-dimensional coordinates from cross-dataset temporal concordance,
using training datasets only in each LODO fold. This is a representation-
identifiability diagnostic, not an ODE/state-space fit.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"stage2_9_27"
OUT.mkdir(parents=True,exist_ok=True)
TARGET=["GSE67462","GSE28688","GSE297234"]
N_BOOT=200
N_PERM=500
MIN_TRAIN_GENES=100


def log(x): print(f"Stage 2.9.27: {x}",flush=True)

def corr(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok])))

def rmse(a,b): return float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)))

def load():
    from dynamics.validation import _load_common_space
    m,meta=_load_common_space()
    meta=meta[meta.dataset.astype(str).isin(TARGET)].copy()
    meta["matrix_column"]=meta.matrix_column.astype(str)
    cols=[c for c in meta.matrix_column if c in m.columns]
    meta=meta[meta.matrix_column.isin(cols)].copy().set_index("matrix_column").loc[cols].reset_index()
    if "time_hours" not in meta.columns:
        from dynamics.validation import _time_hours_for_validation,_strip_dataset_prefix
        meta["time_hours"]=[_time_hours_for_validation(str(r.dataset),_strip_dataset_prefix(str(r.sample)),i if str(r.dataset)=="GSE28688" else None) for i,r in meta.iterrows()]
    meta["time_hours"]=pd.to_numeric(meta.time_hours,errors="coerce")
    return m.loc[:,cols],meta

def dataset_trajectories(m,meta):
    out={}
    for ds in TARGET:
        g=meta[meta.dataset.astype(str).eq(ds)&np.isfinite(meta.time_hours)].copy()
        if g.time_hours.nunique()<3:continue
        # Dataset mean at each observed time; genes are the candidate coordinates.
        X=m[g.matrix_column].T.copy()
        X["time_hours"]=g.time_hours.to_numpy()
        X=X.groupby("time_hours").mean().sort_index()
        if len(X)>=3: out[ds]=(X.index.to_numpy(float),X.drop(columns="time_hours").to_numpy(float),list(m.index))
    return out

def invariant_scores(train):
    # Fit a low-order temporal shape independently per training dataset.
    # The score is cross-dataset agreement of normalized trajectory coefficients;
    # it never uses the held-out dataset.
    ds=list(train); genes=train[ds[0]][2]; coefs=[]
    for d in ds:
        t,X,_=train[d]; tn=(t-t.min())/(t.max()-t.min())
        B=np.column_stack([np.ones(len(tn)),tn,tn*tn])
        coef=np.linalg.lstsq(B,X,rcond=None)[0]
        scale=np.linalg.norm(coef[1:],axis=0)
        scale[scale==0]=1
        coefs.append(coef[1:]/scale)
    A=np.stack(coefs)
    mean=A.mean(axis=0)
    score=np.nanmean(np.sum(A*mean[None,:,:],axis=1),axis=0)
    magnitude=np.nanmean([np.linalg.norm(c,axis=0) for c in coefs],axis=0)
    return pd.DataFrame({"gene":genes,"invariance_score":score,"temporal_magnitude":magnitude,"selection_score":score*magnitude})

def fit_axis(train,genes,k=1):
    rows=[]
    for d,(t,X,all_genes) in train.items():
        idx=[all_genes.index(g) for g in genes]
        Y=X[:,idx]
        Y=StandardScaler().fit_transform(Y)
        tn=(t-t.min())/(t.max()-t.min())
        for i in range(len(t)):
            rows.append((d,float(tn[i]),Y[i]))
    A=np.vstack([r[2] for r in rows])
    pca=PCA(n_components=k).fit(A)
    return pca

def project_axis(pca,train,genes,dsdata):
    all_genes=dsdata[2]; idx=[all_genes.index(g) for g in genes]
    # Scale held-out samples using the training pooled gene means/stds.
    trainY=np.vstack([v[1][:,[v[2].index(g) for g in genes]] for v in train.values()])
    sc=StandardScaler().fit(trainY)
    def proj(X): return pca.transform(sc.transform(X[:,idx]))
    return proj

def one_step(t,z):
    return [(i,float(t[i]),float(t[i+1]),z[i].copy(),z[i+1].copy()) for i in range(len(t)-1)]

def run():
    m,meta=load(); traj=dataset_trajectories(m,meta)
    log(f"trajectory datasets: {', '.join(sorted(traj))}")
    fold_rows=[]; axis_rows=[]; gene_rows=[]; pred_rows=[]; perm_rows=[]
    for fold,held in enumerate(TARGET,1):
        train={d:v for d,v in traj.items() if d!=held}
        if held not in traj or len(train)<2: continue
        log(f"fold {fold}/3: held out {held}")
        scores=invariant_scores(train)
        scores=scores.sort_values("selection_score",ascending=False)
        genes=scores.loc[scores.invariance_score>=0.75,"gene"].astype(str).head(2000).tolist()
        if len(genes)<MIN_TRAIN_GENES: genes=scores.head(MIN_TRAIN_GENES)["gene"].astype(str).tolist()
        gene_rows.append({"held_out_dataset":held,"n_training_datasets":len(train),"n_selected_genes":len(genes),"selection_rule":"cross-dataset normalized temporal-shape concordance; training only","median_invariance_score":float(scores.head(len(genes)).invariance_score.median())})
        for k in [1,2,3]:
            pca=fit_axis(train,genes,k=k)
            axis_rows.append({"held_out_dataset":held,"n_coordinates":k,"explained_variance":float(pca.explained_variance_ratio_.sum())})
        pca=fit_axis(train,genes,k=1)
        # Refit scaler explicitly on training data for projection.
        trainY=np.vstack([v[1][:,[v[2].index(g) for g in genes]] for v in train.values()])
        sc=StandardScaler().fit(trainY)
        train_axes={d:(t,pca.transform(sc.transform(X[:,[allg.index(g) for g in genes]]))[:,0]) for d,(t,X,allg) in train.items()}
        test_t,test_X,allg=traj[held]
        test_axis=pca.transform(sc.transform(test_X[:,[allg.index(g) for g in genes]]))[:,0]
        # Orient only with training temporal direction.
        tt=np.concatenate([t for t,z in train_axes.values()]); zz=np.concatenate([z for t,z in train_axes.values()])
        if corr(tt,zz)<0:
            for d,(t,z) in train_axes.items():train_axes[d]=(t,-z)
            test_axis=-test_axis
        tr=[]
        for d,(t,z) in train_axes.items():tr.extend(one_step(t,z))
        te=one_step(test_t,test_axis)
        if len(tr)<3 or len(te)<2: continue
        X=np.array([r[3] for r in tr]).reshape(-1,1);Y=np.array([r[4] for r in tr])
        model=LinearRegression().fit(X,Y)
        pred=model.predict(np.array([r[3] for r in te]).reshape(-1,1))
        true=np.array([r[4] for r in te]);cur=np.array([r[3] for r in te])
        flat_t=np.concatenate([t for t,z in train_axes.values()]);flat_z=np.concatenate([z for t,z in train_axes.values()])
        lt=LinearRegression().fit(flat_t.reshape(-1,1),flat_z)
        nearest=np.array([flat_z[int(np.argmin(np.abs(flat_t-r[2])))] for r in te])
        linear=lt.predict(np.array([r[2] for r in te]).reshape(-1,1))
        obs=rmse(pred,true); persistence=rmse(cur,true); near=rmse(nearest,true); lin=rmse(linear,true)
        fold_rows.append({"held_out_dataset":held,"n_training_datasets":len(train),"n_selected_genes":len(genes),"n_training_transitions":len(tr),"n_test_transitions":len(te),"train_axis_time_correlation":corr(tt,zz),"test_progress_spearman":corr(np.array([r[2] for r in te]),true),"state_model_rmse":obs,"persistence_rmse":persistence,"nearest_time_rmse":near,"linear_time_rmse":lin,"improvement_vs_persistence":persistence-obs,"improvement_vs_nearest_time":near-obs,"improvement_vs_linear_time":lin-obs})
        for r,a,b in zip(te,pred,true): pred_rows.append({"held_out_dataset":held,"previous_time_hours":r[1],"future_time_hours":r[2],"true_state":float(b),"predicted_state":float(a),"persistence_state":float(r[3])})
        rng=np.random.default_rng(27000+fold); null=np.empty(N_PERM)
        for j in range(N_PERM): null[j]=rmse(pred,rng.permutation(true))
        perm_rows.append({"held_out_dataset":held,"observed_rmse":obs,"permutation_p_rmse":float((1+np.sum(null<=obs))/(N_PERM+1)),"null_mean_rmse":float(null.mean()),"null_p05_rmse":float(np.quantile(null,.05)),"null_p95_rmse":float(np.quantile(null,.95))})
        # Bootstrap genes to test whether the discovered axis is identifiable.
        rng=np.random.default_rng(28000+fold); base=pca.components_[0]; cos=[]
        gene_idx=np.arange(len(genes))
        for _ in range(N_BOOT):
            sel=rng.choice(gene_idx,size=len(gene_idx),replace=True)
            A=np.vstack([v[1][:,[v[2].index(genes[i]) for i in sel]] for v in train.values()])
            A=StandardScaler().fit_transform(A); pc=PCA(n_components=1).fit(A).components_[0]
            # Compare with the original direction after matching duplicated bootstrap columns.
            ref=base[sel]
            den=np.linalg.norm(pc)*np.linalg.norm(ref)
            if den: cos.append(abs(float(np.dot(pc,ref)/den)))
        axis_rows.append({"held_out_dataset":held,"n_coordinates":"bootstrap_axis","bootstrap_cosine_mean":float(np.mean(cos)) if cos else np.nan,"bootstrap_cosine_p05":float(np.quantile(cos,.05)) if cos else np.nan,"bootstrap_cosine_p95":float(np.quantile(cos,.95)) if cos else np.nan})
    f=pd.DataFrame(fold_rows);g=pd.DataFrame(gene_rows);a=pd.DataFrame(axis_rows);p=pd.DataFrame(pred_rows);pm=pd.DataFrame(perm_rows)
    if len(f):
        summary=pd.DataFrame([{"n_valid_folds":len(f),"mean_state_model_rmse":f.state_model_rmse.mean(),"mean_persistence_rmse":f.persistence_rmse.mean(),"mean_nearest_time_rmse":f.nearest_time_rmse.mean(),"mean_linear_time_rmse":f.linear_time_rmse.mean(),"mean_improvement_vs_persistence":f.improvement_vs_persistence.mean(),"mean_improvement_vs_nearest_time":f.improvement_vs_nearest_time.mean(),"mean_improvement_vs_linear_time":f.improvement_vs_linear_time.mean(),"mean_test_progress_spearman":f.test_progress_spearman.mean(),"min_permutation_p_rmse":pm.permutation_p_rmse.min() if len(pm) else np.nan}])
    else: summary=pd.DataFrame()
    supported=False
    if len(summary):
        r=summary.iloc[0]; supported=bool(r.n_valid_folds>=2 and r.mean_improvement_vs_persistence>0 and r.mean_improvement_vs_nearest_time>0 and r.mean_improvement_vs_linear_time>0 and r.min_permutation_p_rmse<0.05)
    overall=pd.DataFrame([{"n_trajectory_datasets":len(traj),"n_valid_lodo_folds":len(f),"mean_selected_genes":float(g.n_selected_genes.mean()) if len(g) else np.nan,"mean_state_improvement_vs_persistence":summary.iloc[0].mean_improvement_vs_persistence if len(summary) else np.nan,"mean_state_improvement_vs_nearest_time":summary.iloc[0].mean_improvement_vs_nearest_time if len(summary) else np.nan,"min_permutation_p_rmse":summary.iloc[0].min_permutation_p_rmse if len(summary) else np.nan,"invariant_coordinate_predictive_support":supported,"stage3_readiness":False,"interpretation":"Training-only cross-dataset temporal-shape concordance; one-dimensional invariant coordinate; bootstrap axis stability; held-out one-step prediction; no ODE/state-space model"}])
    g.to_csv(OUT/"01_gene_selection.csv",index=False);f.to_csv(OUT/"02_lodo_results.csv",index=False);a.to_csv(OUT/"03_axis_stability.csv",index=False);p.to_csv(OUT/"04_one_step_predictions.csv",index=False);pm.to_csv(OUT/"05_permutation_null.csv",index=False);summary.to_csv(OUT/"06_summary.csv",index=False);overall.to_csv(OUT/"07_stage2927_summary.csv",index=False)
    log("overall:");print(overall.to_string(index=False));return overall

if __name__=="__main__": run()
