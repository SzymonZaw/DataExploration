"""Stage 2.9.25: leakage-free state identifiability / Markov sufficiency.

Tests whether the current latent state is sufficient for predicting the next
observed state, or whether adding the previous state materially improves
one-step prediction. No ODE/state-space model is fitted.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"stage2_9_25"; OUT.mkdir(parents=True,exist_ok=True)
TARGET=["GSE67462","GSE28688","GSE297234"]
VARIANTS=("baseline_all","time_dominant")
N_PERM=1000
MIN_TRAIN_TRIPLES=4
MIN_TEST_TRIPLES=2

def log(x): print(f"Stage 2.9.25: {x}",flush=True)

def corr(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<2 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok])))

def load():
    from dynamics.validation import _load_common_space
    m,meta=_load_common_space(); meta=meta[meta.dataset.astype(str).isin(TARGET)].copy()
    meta["matrix_column"]=meta.matrix_column.astype(str); cols=[c for c in meta.matrix_column if c in m.columns]
    meta=meta[meta.matrix_column.isin(cols)].copy().set_index("matrix_column").loc[cols].reset_index()
    if "time_hours" in meta: meta["time_hours"]=pd.to_numeric(meta.time_hours,errors="coerce")
    else:
        from dynamics.validation import _time_hours_for_validation,_strip_dataset_prefix
        vals=[]
        for i,r in meta.iterrows(): vals.append(_time_hours_for_validation(str(r.dataset),_strip_dataset_prefix(str(r.sample)),i if str(r.dataset)=="GSE28688" else None))
        meta["time_hours"]=vals
    return m.loc[:,cols],meta

def feature_sets(m,meta):
    from dynamics.stage2924 import feature_sets
    var,sets=feature_sets(m,meta)
    return sets

def project(m,meta,train_ds,test_ds,genes):
    tr=meta[meta.dataset.astype(str).isin(train_ds)&np.isfinite(meta.time_hours)].copy(); te=meta[meta.dataset.astype(str).eq(test_ds)&np.isfinite(meta.time_hours)].copy()
    if len(tr)<4 or len(te)<3 or tr.dataset.nunique()<2:return None
    genes=[g for g in genes if g in m.index]
    if len(genes)<100:return None
    Xtr=np.array(m.loc[genes,tr.matrix_column].T.to_numpy(float),copy=True); Xte=np.array(m.loc[genes,te.matrix_column].T.to_numpy(float),copy=True)
    med=np.nanmedian(Xtr,axis=0); med=np.where(np.isfinite(med),med,0.)
    for X in (Xtr,Xte):
        bad=~np.isfinite(X)
        if bad.any():X[bad]=np.take(med,np.where(bad)[1])
    sc=StandardScaler().fit(Xtr); pca=PCA(n_components=min(10,Xtr.shape[0]-1,Xtr.shape[1]),random_state=0).fit(sc.transform(Xtr))
    tr["state"]=pca.transform(sc.transform(Xtr))[:,0]; te["state"]=pca.transform(sc.transform(Xte))[:,0]
    out={}
    for ds,g in pd.concat([tr,te]).groupby("dataset"):
        q=g.groupby("time_hours",as_index=False).state.mean().sort_values("time_hours")
        if len(q)>=3:out[str(ds)]=(q.time_hours.to_numpy(float),q.state.to_numpy(float))
    return out

def triples(t,z):
    return [(float(t[i-2]),float(z[i-2]),float(t[i-1]),float(z[i-1]),float(t[i]),float(z[i])) for i in range(2,len(t)) if np.all(np.isfinite([t[i-2],z[i-2],t[i-1],z[i-1],t[i],z[i]])) and t[i-2]<t[i-1]<t[i]]

def eval_fold(m,meta,heldout,variant,genes):
    train=[d for d in TARGET if d!=heldout]; trj=project(m,meta,train,heldout,genes)
    if not trj or heldout not in trj:return None,[]
    train_rows=[]; test_rows=triples(*trj[heldout])
    for d in train:
        if d in trj:train_rows.extend(triples(*trj[d]))
    if len(train_rows)<MIN_TRAIN_TRIPLES or len(test_rows)<MIN_TEST_TRIPLES:return None,[]
    Xcur=np.asarray([[r[3]] for r in train_rows]); Xhist=np.asarray([[r[3],r[1]] for r in train_rows]); y=np.asarray([r[5] for r in train_rows])
    cur=LinearRegression().fit(Xcur,y); hist=LinearRegression().fit(Xhist,y)
    Xtcur=np.asarray([[r[3]] for r in test_rows]); Xth=np.asarray([[r[3],r[1]] for r in test_rows]); yt=np.asarray([r[5] for r in test_rows])
    pc=cur.predict(Xtcur); ph=hist.predict(Xth)
    rm=lambda a:float(np.sqrt(np.mean((a-yt)**2)))
    row={"status":"ok","held_out_dataset":heldout,"variant":variant,"n_genes":len(genes),"n_training_triples":len(train_rows),"n_test_triples":len(test_rows),"current_only_rmse":rm(pc),"history_augmented_rmse":rm(ph),"rmse_improvement_history":rm(pc)-rm(ph),"current_only_pearson":corr(pc,yt),"history_augmented_pearson":corr(ph,yt),"current_only_spearman":corr(pd.Series(pc).rank(),pd.Series(yt).rank()),"history_augmented_spearman":corr(pd.Series(ph).rank(),pd.Series(yt).rank()),"history_coef_previous_state":float(hist.coef_[1]),"history_intercept":float(hist.intercept_)}
    preds=[]
    for r,a,b in zip(test_rows,pc,ph):preds.append({"held_out_dataset":heldout,"variant":variant,"previous_time_hours":r[0],"previous_state":r[1],"current_time_hours":r[2],"current_state":r[3],"future_time_hours":r[4],"true_future_state":r[5],"current_only_prediction":a,"history_augmented_prediction":b})
    return row,preds

def perm_fold(m,meta,heldout,variant,genes,seed):
    train=[d for d in TARGET if d!=heldout]; trj=project(m,meta,train,heldout,genes)
    if not trj or heldout not in trj:return None
    rows=[]
    for d in train:
        if d in trj:rows.extend(triples(*trj[d]))
    test=triples(*trj[heldout])
    if len(rows)<MIN_TRAIN_TRIPLES or len(test)<MIN_TEST_TRIPLES:return None
    model=LinearRegression().fit(np.asarray([[r[3],r[1]] for r in rows]),np.asarray([r[5] for r in rows])); pred=model.predict(np.asarray([[r[3],r[1]] for r in test])); y=np.asarray([r[5] for r in test]); obs=float(np.sqrt(np.mean((pred-y)**2)))
    rng=np.random.default_rng(seed); null=np.empty(N_PERM)
    for b in range(N_PERM):null[b]=np.sqrt(np.mean((pred-rng.permutation(y))**2))
    return {"held_out_dataset":heldout,"variant":variant,"observed_history_rmse":obs,"permutation_p_rmse":float((1+np.sum(null<=obs))/(N_PERM+1)),"null_mean_rmse":float(null.mean()),"null_p05_rmse":float(np.quantile(null,.05)),"null_p95_rmse":float(np.quantile(null,.95))}

def run():
    log("loading common gene space"); m,meta=load(); folds=[]; preds=[]; perms=[]; selections=[]
    for i,h in enumerate(TARGET,1):
        log(f"fold {i}/3: held out {h}"); trmeta=meta[meta.dataset.astype(str)!=h].copy(); sets=feature_sets(m.loc[:,trmeta.matrix_column],trmeta)
        for v in VARIANTS:
            log(f"  {v}: {len(sets[v])} training-selected genes")
            r,p=eval_fold(m,meta,h,v,sets[v]);
            if r is not None:folds.append(r);preds.extend(p)
            q=perm_fold(m,meta,h,v,sets[v],25000+i); 
            if q is not None:perms.append(q)
            selections.append({"held_out_dataset":h,"variant":v,"n_genes":len(sets[v])})
    f=pd.DataFrame(folds); p=pd.DataFrame(preds); pm=pd.DataFrame(perms,columns=["held_out_dataset","variant","observed_history_rmse","permutation_p_rmse","null_mean_rmse","null_p05_rmse","null_p95_rmse"]); sel=pd.DataFrame(selections)
    sums=[]
    for v in VARIANTS:
        g=f[(f.variant==v)&(f.status=="ok")]; q=pm[pm.variant==v] if len(pm) else pd.DataFrame()
        sums.append({"variant":v,"n_valid_folds":len(g),"mean_current_only_rmse":g.current_only_rmse.mean() if len(g) else np.nan,"mean_history_augmented_rmse":g.history_augmented_rmse.mean() if len(g) else np.nan,"mean_rmse_improvement_history":g.rmse_improvement_history.mean() if len(g) else np.nan,"mean_current_only_pearson":g.current_only_pearson.mean() if len(g) else np.nan,"mean_history_augmented_pearson":g.history_augmented_pearson.mean() if len(g) else np.nan,"min_permutation_p_rmse":q.permutation_p_rmse.min() if len(q) else np.nan})
    sm=pd.DataFrame(sums); pri=sm[sm.variant=="time_dominant"].iloc[0]
    supported=bool(pri.n_valid_folds>=2 and pri.mean_rmse_improvement_history>0 and np.isfinite(pri.min_permutation_p_rmse) and pri.min_permutation_p_rmse<0.05)
    summary=pd.DataFrame([{"n_heldout_datasets":len(TARGET),"primary_variant":"time_dominant","primary_n_valid_folds":int(pri.n_valid_folds),"primary_mean_rmse_improvement_history":pri.mean_rmse_improvement_history,"primary_min_permutation_p_rmse":pri.min_permutation_p_rmse,"markov_sufficiency_supported":supported,"stage3_readiness":False,"interpretation":"LODO test of Markov sufficiency: compare z(t)->z(t+1) with history-augmented z(t),z(t-1); train-only gene selection/scaling/PCA; no ODE/state-space model"}])
    sel.to_csv(OUT/"01_training_gene_selection.csv",index=False); f.to_csv(OUT/"02_fold_results.csv",index=False); p.to_csv(OUT/"03_one_step_predictions.csv",index=False); pm.to_csv(OUT/"04_permutation_null.csv",index=False); sm.to_csv(OUT/"05_variant_summary.csv",index=False); summary.to_csv(OUT/"06_stage2925_summary.csv",index=False)
    log("variant summary:"); print(sm.to_string(index=False)); log("overall:"); print(summary.to_string(index=False)); return summary

if __name__=="__main__":run()
