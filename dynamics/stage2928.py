"""Stage 2.9.28: leakage-free shared temporal biological component discovery.

Decomposes training trajectories into dataset-centered temporal shapes, extracts
features that reproduce across independent experiments, and evaluates whether
the resulting shared component is reproducible in a held-out dataset.
No ODE/state-space model is fitted.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"stage2_9_28"
OUT.mkdir(parents=True,exist_ok=True)
TARGET=["GSE67462","GSE28688","GSE297234"]
GRID=np.linspace(0,1,9)
MAX_GENES=3000
N_BOOT=200
N_PERM=1000


def log(x): print(f"Stage 2.9.28: {x}",flush=True)

def finite(X):
    X=np.array(X,float,copy=True)
    if np.isfinite(X).all(): return X
    med=np.zeros(X.shape[1],float)
    for j in range(X.shape[1]):
        vals=X[np.isfinite(X[:,j]),j]
        med[j]=float(np.median(vals)) if len(vals) else 0.0
    r,c=np.where(~np.isfinite(X))
    if len(r): X[r,c]=med[c]
    return X

def corr(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if a.shape!=b.shape:
        n=min(a.size,b.size);a=a.reshape(-1)[:n];b=b.reshape(-1)[:n]
    ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok])))

def interp_shape(t,X):
    t=np.asarray(t,float); X=finite(X)
    if len(t)<2:return np.repeat(X[:1],len(GRID),axis=0)
    tn=(t-t.min())/(t.max()-t.min())
    out=np.empty((len(GRID),X.shape[1]))
    for g in range(X.shape[1]): out[:,g]=np.interp(GRID,tn,X[:,g])
    return out

def load():
    from dynamics.validation import _load_common_space
    m,meta=_load_common_space()
    meta=meta[meta["dataset"].astype(str).isin(TARGET)].copy()
    meta["matrix_column"]=meta["matrix_column"].astype(str)
    cols=[c for c in meta["matrix_column"] if c in m.columns]
    meta=meta[meta["matrix_column"].isin(cols)].copy()
    return m.loc[:,cols],meta

def trajectories(m,meta):
    out={}
    for ds in TARGET:
        g=meta[meta["dataset"].astype(str).eq(ds)].copy()
        if "time_hours" not in g.columns:
            from dynamics.validation import _time_hours_for_validation,_strip_dataset_prefix
            g["time_hours"]=[_time_hours_for_validation(ds,_strip_dataset_prefix(str(s)),i if ds=="GSE28688" else None) for i,s in enumerate(g["sample"])]
        g["time_hours"]=pd.to_numeric(g["time_hours"],errors="coerce")
        g=g[np.isfinite(g["time_hours"])].copy()
        if g["time_hours"].nunique()<3:continue
        cols=list(g["matrix_column"])
        X=m.loc[:,cols].T.copy();X["time_hours"]=g["time_hours"].to_numpy(float)
        X=X.groupby("time_hours",sort=True).mean()
        out[ds]=(X.index.to_numpy(float),finite(X.to_numpy(float)),list(m.index))
    return out

def shared_shapes(train):
    ds=list(train); genes=train[ds[0]][2]; shapes=[]
    for d in ds:
        t,X,_=train[d]
        Y=interp_shape(t,X);Y=Y-Y[0:1,:]
        amp=np.sqrt(np.mean(Y*Y,axis=0));amp=np.where(amp>1e-8,amp,1.0)
        shapes.append(Y/amp)
    A=np.stack(shapes)
    shared=A.mean(axis=0)
    hetero=np.mean((A-shared[None,:,:])**2,axis=(0,1))
    shared_var=np.var(shared,axis=0)
    score=shared_var/(shared_var+hetero+1e-12)
    magnitude=np.sqrt(np.mean(shared**2,axis=0))
    return genes,shared,score,magnitude,A

def select(train):
    genes,shared,score,mag,A=shared_shapes(train)
    tab=pd.DataFrame({"gene":genes,"shared_temporal_variance":np.var(shared,axis=0),"dataset_heterogeneity":np.mean((A-shared[None,:,:])**2,axis=(0,1)),"shared_fraction":score,"shared_magnitude":mag})
    tab["selection_score"]=tab["shared_fraction"]*tab["shared_magnitude"]
    tab=tab.sort_values(["shared_fraction","selection_score"],ascending=False)
    chosen=tab.loc[tab["shared_fraction"]>=0.5,"gene"].head(MAX_GENES).tolist()
    if len(chosen)<100: chosen=tab.head(min(100,MAX_GENES))["gene"].tolist()
    return tab,chosen,shared

def fit_template(train,genes):
    blocks=[]
    for _,(t,X,allg) in train.items():
        idx=[allg.index(g) for g in genes]
        Y=interp_shape(t,X[:,idx]);Y=Y-Y[0:1,:]
        blocks.append(Y)
    A=np.vstack(blocks)
    mean=A.mean(axis=0);sd=A.std(axis=0);sd=np.where(sd>1e-8,sd,1.0)
    Z=(A-mean)/sd
    pca=PCA(n_components=1).fit(Z)
    return mean,sd,pca,pca.transform(Z)[:,0]

def project(t,X,genes,allg,mean,sd,pca):
    idx=[allg.index(g) for g in genes]
    Y=interp_shape(t,X[:,idx]);Y=Y-Y[0:1,:]
    return pca.transform((Y-mean)/sd)[:,0]

def run():
    m,meta=load();traj=trajectories(m,meta)
    log(f"trajectory datasets: {', '.join(sorted(traj))}")
    fold=[];gene_all=[];boot=[];perm=[];template_rows=[]
    for fi,held in enumerate(TARGET,1):
        if held not in traj:continue
        train={d:v for d,v in traj.items() if d!=held}
        if len(train)<2:continue
        log(f"fold {fi}/3: held out {held}")
        tab,genes,_=select(train)
        tab.insert(0,"held_out_dataset",held);gene_all.append(tab)
        mean,sd,pca,_=fit_template(train,genes)
        test_t,test_X,allg=traj[held]
        zte=project(test_t,test_X,genes,allg,mean,sd,pca)
        # project() already returns the held-out trajectory on the common GRID.
        # The original held-out time vector can have a different number of points.
        zte_grid=zte
        time_corr=corr(GRID,zte_grid)
        rng=np.random.default_rng(29000+fi);cos=[]
        base_sign=np.sign(time_corr) if np.isfinite(time_corr) and time_corr!=0 else 1.0
        for _ in range(N_BOOT):
            sel=rng.choice(len(genes),len(genes),replace=True)
            bg=[genes[i] for i in sel]
            bmean,bsd,bpca,_=fit_template(train,bg)
            bz=project(test_t,test_X,bg,allg,bmean,bsd,bpca)
            bz_grid=bz
            c=corr(zte_grid*base_sign,bz_grid)
            if np.isfinite(c):cos.append(abs(c))
        boot.append({"held_out_dataset":held,"n_selected_genes":len(genes),"axis_bootstrap_correlation_mean":np.mean(cos) if cos else np.nan,"axis_bootstrap_correlation_p05":np.quantile(cos,.05) if cos else np.nan,"axis_bootstrap_correlation_p95":np.quantile(cos,.95) if cos else np.nan})
        obs=np.abs(time_corr)
        rng=np.random.default_rng(29100+fi);null=[]
        for _ in range(N_PERM): null.append(np.abs(corr(GRID,rng.permutation(zte_grid))))
        perm.append({"held_out_dataset":held,"observed_abs_time_correlation":obs,"permutation_p":(1+np.sum(np.asarray(null)>=obs))/(N_PERM+1),"null_mean":np.mean(null)})
        fold.append({"held_out_dataset":held,"n_training_datasets":len(train),"n_selected_genes":len(genes),"mean_shared_fraction":float(tab.head(len(genes))["shared_fraction"].mean()),"test_time_correlation":time_corr,"test_abs_time_correlation":obs})
        template_rows.append(pd.DataFrame({"held_out_dataset":held,"normalized_time":GRID,"shared_coordinate":zte_grid}))
    F=pd.DataFrame(fold);G=pd.concat(gene_all,ignore_index=True) if gene_all else pd.DataFrame();B=pd.DataFrame(boot);P=pd.DataFrame(perm);T=pd.concat(template_rows,ignore_index=True) if template_rows else pd.DataFrame()
    summary=pd.DataFrame([{"n_trajectory_datasets":len(traj),"n_valid_lodo_folds":len(F),"mean_selected_genes":F["n_selected_genes"].mean() if len(F) else np.nan,"mean_test_abs_time_correlation":F["test_abs_time_correlation"].mean() if len(F) else np.nan,"mean_bootstrap_axis_correlation":B["axis_bootstrap_correlation_mean"].mean() if len(B) else np.nan,"min_bootstrap_p05":B["axis_bootstrap_correlation_p05"].min() if len(B) else np.nan,"max_permutation_p":P["permutation_p"].max() if len(P) else np.nan,"shared_temporal_component_supported":bool(len(F)>=2 and B["axis_bootstrap_correlation_p05"].min()>0.7 and P["permutation_p"].max()<0.05) if len(B) and len(P) else False,"stage3_readiness":False,"interpretation":"Training-only dataset-centered temporal decomposition; shared component extracted before held-out evaluation; bootstrap gene stability; held-out temporal reproducibility; NaN-safe projection; no ODE/state-space model"}])
    G.to_csv(OUT/"01_gene_shared_variance.csv",index=False);F.to_csv(OUT/"02_lodo_shared_component.csv",index=False);B.to_csv(OUT/"03_bootstrap_axis_stability.csv",index=False);P.to_csv(OUT/"04_time_permutation_null.csv",index=False);T.to_csv(OUT/"05_heldout_shared_coordinate.csv",index=False);summary.to_csv(OUT/"06_stage2928_summary.csv",index=False)
    log("overall:");print(summary.to_string(index=False));return summary

if __name__=="__main__":run()
