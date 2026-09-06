"""Stage 2.8 diagnostics for cross-dataset generalisation."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
STAGE26=ROOT/"results"/"Dynamics"/"stage2_6"
STAGE27=ROOT/"results"/"Dynamics"/"stage2_7"
OUT=ROOT/"results"/"Dynamics"/"stage2_8"
OUT.mkdir(parents=True,exist_ok=True)

def metrics(a,b):
    a=np.asarray(a,float).ravel(); b=np.asarray(b,float).ravel(); ok=np.isfinite(a)&np.isfinite(b)
    if not ok.any(): return (np.nan,np.nan,np.nan)
    a,b=a[ok],b[ok]; c=np.corrcoef(a,b)[0,1] if len(a)>1 and np.std(a)>0 and np.std(b)>0 else np.nan
    return (float(np.sqrt(np.mean((a-b)**2))),float(np.mean(np.abs(a-b))),float(c))

def traj(ds,matrix,meta):
    g=meta[(meta.dataset==ds)&meta.matrix_column.notna()&meta.time_hours.notna()].copy()
    if len(g)<2:return None
    f=pd.DataFrame(matrix[g.matrix_column].T.to_numpy(),index=g.time_hours.to_numpy(float),columns=matrix.index)
    return f.groupby(level=0).mean().sort_index()

def interp(frame,t):
    if t<frame.index.min() or t>frame.index.max():return None
    x=frame.index.to_numpy(float); y=frame.to_numpy()
    return np.asarray([np.interp(t,x,y[:,j]) for j in range(y.shape[1])])

def run_diagnostics(matrix,meta):
    ts={d:traj(d,matrix,meta) for d in sorted(meta.dataset.unique())}
    ts={d:f for d,f in ts.items() if f is not None}
    rows=[]
    for ds,test in ts.items():
        train={d:f for d,f in ts.items() if d!=ds}
        for t in test.index:
            actual=test.loc[t].to_numpy(float)
            pred=[interp(f,t) for f in train.values()]; pred=[p for p in pred if p is not None]
            if not pred:continue
            cross=np.mean(pred,axis=0)
            nearest=[f.iloc[int(np.argmin(abs(f.index.to_numpy(float)-t)))].to_numpy(float) for f in train.values()]
            base=np.mean(nearest,axis=0); rm,ma,co=metrics(actual,cross); br,bm,bc=metrics(actual,base)
            rows.append(dict(dataset=ds,time_hours=float(t),n_training_datasets=len(pred),cross_dataset_rmse=rm,cross_dataset_mae=ma,cross_dataset_correlation=co,nearest_time_baseline_rmse=br,nearest_time_baseline_mae=bm,nearest_time_baseline_correlation=bc,rmse_improvement_vs_nearest=br-rm))
    return pd.DataFrame(rows)

def time_coverage(meta):
    rows=[]
    for ds,g in meta.groupby("dataset"):
        t=g[g.time_hours.notna()].time_hours
        n=t.nunique(); rows.append(dict(dataset=ds,n_samples=len(g),n_timed=len(t),n_unique_times=n,min_time_hours=t.min() if len(t) else np.nan,max_time_hours=t.max() if len(t) else np.nan,coverage_span_hours=t.max()-t.min() if len(t) else np.nan,role="trajectory" if n>=2 else ("single_timepoint" if n==1 else "context_only")))
    return pd.DataFrame(rows)

def gene_temporal_concordance(matrix,meta):
    """Score each gene by cross-dataset agreement of its temporal profile.

    For every pair of trajectory datasets, only exact overlapping timepoints are
    used. Scores are computed per gene and then aggregated across pairs. This is
    a diagnostic/selection score, not a claim of biological equivalence.
    """
    ts={d:traj(d,matrix,meta) for d in sorted(meta.dataset.unique())}
    ts={d:f for d,f in ts.items() if f is not None and len(f.index)>=2}
    genes=list(matrix.index)
    pairs=[]
    for i,a in enumerate(sorted(ts)):
        for b in sorted(ts)[i+1:]:
            common=sorted(set(ts[a].index).intersection(ts[b].index))
            if len(common)>=3:pairs.append((a,b,common))
    print(f"Stage 2.8: gene concordance across {len(pairs)} dataset pairs and {len(genes)} genes...")
    accum={g:[] for g in genes}
    for pi,(a,b,common) in enumerate(pairs,1):
        xa=ts[a].loc[common].to_numpy(float); xb=ts[b].loc[common].to_numpy(float)
        for j,g in enumerate(genes):
            va=xa[:,j]; vb=xb[:,j]; ok=np.isfinite(va)&np.isfinite(vb)
            if ok.sum()<3 or np.std(va[ok])==0 or np.std(vb[ok])==0: continue
            r=float(np.corrcoef(va[ok],vb[ok])[0,1])
            sign=float(np.sign(np.sum(np.diff(va[ok])*np.diff(vb[ok]))))
            accum[g].append((r,sign))
        print(f"  pair {pi}/{len(pairs)}: {a} vs {b}, overlap={len(common)}",flush=True)
    rows=[]
    for g,vals in accum.items():
        if not vals: continue
        corr=np.asarray([v[0] for v in vals],float); signs=np.asarray([v[1] for v in vals],float)
        rows.append(dict(gene=g,n_pairs=len(vals),mean_correlation=float(np.mean(corr)),median_correlation=float(np.median(corr)),positive_pair_fraction=float(np.mean(corr>0)),directional_agreement=float(np.mean(signs>0)),concordance_score=float(np.mean(corr)*np.mean(signs>0))) )
    out=pd.DataFrame(rows)
    if not out.empty: out=out.sort_values(["concordance_score","median_correlation","n_pairs"],ascending=[False,False,False]).reset_index(drop=True)
    out.to_csv(OUT/"04_gene_concordance_ranked.csv",index=False)
    return out

def stage2_8():
    from .validation import _load_common_space
    print("Stage 2.8: loading validated common gene space...")
    matrix,meta=_load_common_space(); print(f"Stage 2.8: matrix={matrix.shape[0]} genes x {matrix.shape[1]} samples")
    print("Stage 2.8: per-dataset cross-validation and nearest-time baseline...")
    per=run_diagnostics(matrix,meta); per.to_csv(OUT/"01_per_dataset_diagnostics.csv",index=False); print(f"Stage 2.8: wrote {len(per)} timepoint diagnostics.")
    cov=time_coverage(meta); cov.to_csv(OUT/"02_time_coverage.csv",index=False); print("Stage 2.8: time coverage written.")
    if not per.empty:
        s=per.groupby("dataset").agg(n_cases=("time_hours","size"),mean_rmse=("cross_dataset_rmse","mean"),mean_mae=("cross_dataset_mae","mean"),mean_correlation=("cross_dataset_correlation","mean"),mean_nearest_rmse=("nearest_time_baseline_rmse","mean"),mean_rmse_improvement=("rmse_improvement_vs_nearest","mean")).reset_index()
    else:s=pd.DataFrame()
    s.to_csv(OUT/"03_per_dataset_summary.csv",index=False)
    gene=gene_temporal_concordance(matrix,meta)
    print(f"Stage 2.8: wrote {len(gene)} gene concordance scores.")
    print("Stage 2.8 complete.")
    print("Per-dataset summary:"); print(s.to_string(index=False) if not s.empty else "No diagnostics available.")
    if not gene.empty: print("Top gene concordance scores:"); print(gene.head(20).to_string(index=False))
    return s

if __name__=="__main__": stage2_8()
