"""Stage 2.9.17: program-transition and lag diagnostics.

Uses fixed biological program activities from Stage 2.9.14. The goal is to
ask whether programs show reproducible lead/lag structure rather than merely
tracking absolute time. This is descriptive and does not fit an ODE.
"""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
SRC=ROOT/"results"/"Dynamics"/"stage2_9_14"
OUT=ROOT/"results"/"Dynamics"/"stage2_9_17"; OUT.mkdir(parents=True,exist_ok=True)

def log(x): print(f"Stage 2.9.17: {x}",flush=True)
def _load():
 p=SRC/"03_program_trajectories.csv"
 if not p.exists(): raise RuntimeError("Run Stage 2.9.14 first")
 return pd.read_csv(p)
def _curve(g):
 g=g.sort_values("time_hours"); t=g.time_hours.to_numpy(float); y=g.activity.to_numpy(float)
 if len(np.unique(t))<3:return None
 tn=(t-t.min())/np.ptp(t); grid=np.linspace(0,1,9); return grid,np.interp(grid,tn,y)
def _lag(a,b,max_lag=3):
 vals=[]
 for lag in range(-max_lag,max_lag+1):
  if lag<0:x,y=a[:lag],b[-lag:]
  elif lag>0:x,y=a[lag:],b[:-lag]
  else:x,y=a,b
  if len(x)>=4 and np.std(x)>1e-12 and np.std(y)>1e-12: vals.append((lag,float(pd.Series(x).corr(pd.Series(y))),len(x)))
 if not vals:return (np.nan,np.nan,0)
 return max(vals,key=lambda z:abs(z[1]))
def run():
 log("starting program lead/lag diagnostics")
 df=_load(); curves={}
 for (ds,pid),g in df.groupby(["dataset","program_id"]):
  c=_curve(g)
  if c is not None:curves[(ds,pid)]=c[1]
 rows=[]
 pids=sorted(df.program_id.unique()); datasets=sorted(df.dataset.unique())
 for ds in datasets:
  for i,a in enumerate(pids):
   if (ds,a) not in curves:continue
   for b in pids[i+1:]:
    if (ds,b) not in curves:continue
    lag,r,n=_lag(curves[(ds,a)],curves[(ds,b)])
    rows.append({"dataset":ds,"program_a":a,"program_b":b,"best_lag_grid_steps":lag,"max_abs_pearson":r,"n_grid_points":n})
 lagdf=pd.DataFrame(rows);lagdf.to_csv(OUT/"01_within_dataset_program_lags.csv",index=False)
 if len(lagdf):
  summ=lagdf.groupby(["program_a","program_b"],as_index=False).agg(n_datasets=("dataset","nunique"),mean_abs_pearson=("max_abs_pearson",lambda x:np.nanmean(np.abs(x))),median_abs_pearson=("max_abs_pearson",lambda x:np.nanmedian(np.abs(x))),median_best_lag=("best_lag_grid_steps","median")); summ.to_csv(OUT/"02_cross_dataset_lag_summary.csv",index=False)
 else:summ=pd.DataFrame()
 # Directional reproducibility: same non-zero lag sign in >=2 datasets.
 if len(lagdf):
  dirs=lagdf[lagdf.best_lag_grid_steps!=0].groupby(["program_a","program_b"]).agg(n_nonzero=("best_lag_grid_steps","count"),n_positive=("best_lag_grid_steps",lambda x:int(np.sum(np.asarray(x)>0))),n_negative=("best_lag_grid_steps",lambda x:int(np.sum(np.asarray(x)<0)))).reset_index();dirs["directional_consistency"]=dirs.apply(lambda r:max(r.n_positive,r.n_negative)/r.n_nonzero if r.n_nonzero else np.nan,axis=1)
 else:dirs=pd.DataFrame()
 dirs.to_csv(OUT/"03_directional_consistency.csv",index=False)
 overall=pd.DataFrame([{"n_datasets":len(datasets),"n_program_pairs":len(summ),"n_pairs_with_reproducible_direction":int(np.sum((dirs.directional_consistency>=0.67)&(dirs.n_nonzero>=2))) if len(dirs) else 0,"mean_cross_dataset_abs_lag_correlation":float(summ.mean_abs_pearson.mean()) if len(summ) else np.nan,"transition_signal_supported":bool(len(dirs) and np.any((dirs.directional_consistency>=0.67)&(dirs.n_nonzero>=2)))}])
 overall.to_csv(OUT/"04_STAGE2_9_17_SUMMARY.csv",index=False);log("complete");print("\nStage 2.9.17 lag summary:");print(summ.sort_values("mean_abs_pearson",ascending=False).head(20).to_string(index=False) if len(summ) else "none");print("\nStage 2.9.17 gate:");print(overall.to_string(index=False));return overall
if __name__=="__main__":run()
