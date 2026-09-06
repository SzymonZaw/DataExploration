"""Stage 2.9.6: consensus biological-program stabilization.

Builds a consensus gene set from independently discovered Stage 2.9.1
training folds, then evaluates a leakage-free PCA progress axis using only
consensus genes available in each training fold. This reduces dependence on
arbitrary KMeans program labels and provides a stricter gate before dynamics.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT=Path(__file__).resolve().parents[1]
COMMON=ROOT/"results"/"Dynamics"/"stage2_6"
IN=ROOT/"results"/"Dynamics"/"stage2_9_1"
OUT=ROOT/"results"/"Dynamics"/"stage2_9_6"
OUT.mkdir(parents=True,exist_ok=True)

def _corr(a,b,method="spearman"):
    a=np.asarray(a,float); b=np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method=method))

def _load():
    m=pd.read_csv(COMMON/"06_common_human_gene_matrix.csv",index_col=0)
    meta=pd.read_csv(COMMON/"07_common_gene_sample_metadata.csv")
    mem=pd.read_csv(IN/"02_program_membership.csv")
    return m,meta,mem

def _time(meta,ds,c):
    q=meta[(meta.dataset.astype(str)==str(ds))&(meta.matrix_column.astype(str)==str(c))]
    if q.empty:return np.nan
    return pd.to_numeric(q.iloc[0].time_hours,errors="coerce")

def _trajectory(matrix,meta,ds,genes):
    g=meta[(meta.dataset.astype(str)==str(ds))&meta.matrix_column.notna()].copy()
    g["t"]=[_time(meta,ds,c) for c in g.matrix_column]
    g=g[np.isfinite(g.t.to_numpy(float))]
    if len(g)<3 or g.t.nunique()<3:return None
    f=matrix.loc[genes,g.matrix_column.astype(str).tolist()].T
    f.index=g.t.to_numpy(float); f=f.groupby(level=0).mean().sort_index()
    if len(f)<3:return None
    t=f.index.to_numpy(float); span=t[-1]-t[0]
    return (t-t[0])/span,f.to_numpy(float)

def _axis(train,genes):
    cubes=[]
    for ds in train:
        tr=_trajectory(MATRIX,META,ds,genes)
        if tr is None:continue
        t,y=tr; grid=np.linspace(0,1,9); out=np.full((9,len(genes)),np.nan)
        for j in range(len(genes)):
            ok=np.isfinite(y[:,j])
            if ok.sum()>=3:out[:,j]=np.interp(grid,t[ok],y[ok,j])
        if np.isfinite(out).mean()>=.7:cubes.append(out)
    if len(cubes)<2:return None
    X=np.concatenate(cubes); med=np.nanmedian(X,axis=0); X=np.where(np.isfinite(X),X,med); X=np.where(np.isfinite(X),X,0)
    sc=StandardScaler().fit(X); pc=PCA(n_components=1,random_state=42).fit(sc.transform(X)); scores=pc.transform(sc.transform(X))[:,0]
    target=np.tile(np.linspace(0,1,9),len(cubes)); r=_corr(scores,target,"pearson"); sign=1 if np.isfinite(r) and r>=0 else -1
    return sc,pc,sign,scores

def _project(sc,pc,sign,ref,y):
    med=np.nanmedian(y,axis=0); y=np.where(np.isfinite(y),y,med); y=np.where(np.isfinite(y),y,0)
    z=sign*pc.transform(sc.transform(y))[:,0]; lo,hi=np.nanpercentile(ref,[2.5,97.5])
    if hi<=lo:lo,hi=np.nanmin(ref),np.nanmax(ref)
    return np.full(len(z),.5) if hi<=lo else np.clip((z-lo)/(hi-lo),0,1)

def run(min_fold_recurrence=2,seed=42):
    global MATRIX,META
    MATRIX,META,MEM=_load()
    folds=sorted(MEM.held_out_dataset.astype(str).unique())
    # Gene identity is stable across folds; program labels are not. Recurrence
    # therefore acts as a consensus-program proxy without forcing label matching.
    rec=MEM.groupby("gene").held_out_dataset.nunique().sort_values(ascending=False)
    genes2=rec[rec>=min_fold_recurrence].index.astype(str).tolist()
    genes3=rec[rec>=3].index.astype(str).tolist()
    pd.DataFrame({"gene":rec.index,"n_discovery_folds":rec.values,"consensus_ge_2":rec.values>=2,"consensus_ge_3":rec.values>=3}).to_csv(OUT/"02_gene_recurrence.csv",index=False)
    rows=[]; sens=[]
    for ds in folds:
        train=[x for x in folds if x!=ds]
        # use only genes discovered in at least min_fold_recurrence of the
        # training discovery folds, not including the held-out fold
        train_mem=MEM[MEM.held_out_dataset.astype(str).isin(train)]
        rr=train_mem.groupby("gene").held_out_dataset.nunique(); genes=rr[rr>=min_fold_recurrence].index.astype(str).tolist()
        genes=[g for g in genes if g in MATRIX.index]
        if len(genes)<20:
            rows.append({"held_out_dataset":ds,"status":"insufficient_consensus_genes","n_consensus_genes":len(genes)});continue
        ax=_axis(train,genes)
        test=_trajectory(MATRIX,META,ds,genes)
        if ax is None or test is None:
            rows.append({"held_out_dataset":ds,"status":"insufficient_training_or_test_trajectory","n_consensus_genes":len(genes)});continue
        sc,pc,sign,ref=ax; tt,y=test; z=_project(sc,pc,sign,ref,y); mono=float(np.mean(np.diff(z)>=-.02));
        rows.append({"held_out_dataset":ds,"status":"ok","n_consensus_genes":len(genes),"train_axis_time_pearson":_corr(ref,np.tile(np.linspace(0,1,9),len(ref)//9),"pearson"),"progress_spearman":_corr(tt,z,"spearman"),"progress_pearson":_corr(tt,z,"pearson"),"progress_rmse":float(np.sqrt(np.mean((tt-z)**2))),"progress_monotonicity":mono})
        # sensitivity to consensus threshold
        if len(genes3)>=20:
            g3=[g for g in genes3 if g in MATRIX.index]; ax3=_axis(train,g3); te3=_trajectory(MATRIX,META,ds,g3)
            if ax3 and te3:
                sc3,pc3,sg3,r3=ax3; t3,y3=te3; z3=_project(sc3,pc3,sg3,r3,y3); sens.append({"held_out_dataset":ds,"threshold":"3_folds","n_genes":len(g3),"spearman":_corr(t3,z3),"pearson":_corr(t3,z3,"pearson")})
    result=pd.DataFrame(rows); result.to_csv(OUT/"01_consensus_fold_results.csv",index=False); pd.DataFrame(sens).to_csv(OUT/"03_consensus_threshold_sensitivity.csv",index=False)
    ok=result[result.status=="ok"]
    summary=pd.DataFrame([{"n_folds":len(result),"n_valid_folds":len(ok),"consensus_min_fold_recurrence":min_fold_recurrence,"mean_progress_spearman":ok.progress_spearman.mean() if len(ok) else np.nan,"mean_progress_pearson":ok.progress_pearson.mean() if len(ok) else np.nan,"mean_progress_rmse":ok.progress_rmse.mean() if len(ok) else np.nan,"mean_monotonicity":ok.progress_monotonicity.mean() if len(ok) else np.nan,"median_consensus_genes":ok.n_consensus_genes.median() if len(ok) else np.nan,"interpretation":"consensus-gene stabilization; program discovery is represented by recurrence across independent Stage 2.9.1 folds"}]); summary.to_csv(OUT/"04_overall_summary.csv",index=False)
    print("Stage 2.9.6: consensus program stabilization",flush=True); print(result.to_string(index=False),flush=True); print("\nStage 2.9.6 summary",flush=True); print(summary.to_string(index=False),flush=True); print("Stage 2.9.6 complete. ODE/state-space remains gated on reproducible consensus progress.",flush=True); return summary

if __name__=="__main__":run()
