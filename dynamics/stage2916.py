"""Stage 2.9.16: residual biological-state validation.

Removes explicit time trends using training-only LODO models, then asks whether
held-out residual structure is reproducible across datasets. Diagnostic only.
"""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; SRC=ROOT/"results"/"Dynamics"/"stage2_9_14"; OUT=ROOT/"results"/"Dynamics"/"stage2_9_16"; OUT.mkdir(parents=True,exist_ok=True)
def log(x): print(f"Stage 2.9.16: {x}",flush=True)
def corr(a,b,method="pearson"):
 a=np.asarray(a,float);b=np.asarray(b,float);ok=np.isfinite(a)&np.isfinite(b)
 if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
 return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method="spearman" if method=="spearman" else "pearson"))
def _load():
 p=SRC/"03_program_trajectories.csv"
 if not p.exists():raise RuntimeError("Run Stage 2.9.14 first")
 x=pd.read_csv(p);x["dataset"]=x.dataset.astype(str);x["program_id"]=x.program_id.astype(str);return x
def _norm(g):
 g=g.sort_values("time_hours").copy();t=g.time_hours.to_numpy(float)
 if len(g)<3 or np.ptp(t)<=0:return None
 g["normalized_time"]=(t-t.min())/np.ptp(t);return g
def _lodo(df):
 rows=[]
 for held in sorted(df.dataset.unique()):
  train=df[df.dataset!=held];test=df[df.dataset==held];rr=[]
  for pid in sorted(set(train.program_id)&set(test.program_id)):
   tr=_norm(train[train.program_id==pid]);te=_norm(test[test.program_id==pid])
   if tr is None or te is None:continue
   X=np.column_stack([np.ones(len(tr)),tr.normalized_time]);b=np.linalg.lstsq(X,tr.activity.to_numpy(float),rcond=None)[0]
   pred=b[0]+b[1]*te.normalized_time.to_numpy(float);res=te.activity.to_numpy(float)-pred
   for t,r in zip(te.time_hours,res):rr.append({"held_out_dataset":held,"program_id":pid,"time_hours":t,"residual_activity":r})
   rs=corr(res,te.normalized_time.to_numpy(float),"spearman")
   rows.append({"held_out_dataset":held,"program_id":pid,"n_test_timepoints":len(te),"residual_rmse":float(np.sqrt(np.mean(res**2))),"residual_mae":float(np.mean(np.abs(res))),"residual_time_spearman":rs,"residual_time_abs_spearman":abs(rs) if np.isfinite(rs) else np.nan})
  if rr:pd.DataFrame(rr).to_csv(OUT/f"residuals_{held}.csv",index=False)
 out=pd.DataFrame(rows);out.to_csv(OUT/"01_residual_lodo_by_program.csv",index=False)
 s=out.groupby("held_out_dataset").agg(n_programs=("program_id","count"),mean_abs_residual_time_spearman=("residual_time_abs_spearman","mean"),mean_residual_rmse=("residual_rmse","mean"),mean_residual_mae=("residual_mae","mean")).reset_index() if len(out) else pd.DataFrame()
 s.to_csv(OUT/"02_residual_lodo_summary.csv",index=False);return out,s
def _cross(df):
 rows=[];grid=np.linspace(0,1,9)
 for pid,g in df.groupby("program_id"):
  curves=[]
  for ds,h in g.groupby("dataset"):
   h=_norm(h)
   if h is None:continue
   t=h.normalized_time.to_numpy(float);y=h.activity.to_numpy(float);X=np.column_stack([np.ones(len(t)),t]);b=np.linalg.lstsq(X,y,rcond=None)[0]
   curves.append((ds,np.interp(grid,t,y-X@b)))
  for i in range(len(curves)):
   for j in range(i+1,len(curves)):
    rows.append({"program_id":pid,"dataset_a":curves[i][0],"dataset_b":curves[j][0],"residual_pearson":corr(curves[i][1],curves[j][1])})
 out=pd.DataFrame(rows);out.to_csv(OUT/"03_cross_dataset_residual_signal.csv",index=False);return out
def _perm(df,n=500,seed=2916):
 rng=np.random.default_rng(seed);obs=[];null=[]
 for _,g in df.groupby(["dataset","program_id"]):
  h=_norm(g)
  if h is None:continue
  t=h.normalized_time.to_numpy(float);y=h.activity.to_numpy(float);X=np.column_stack([np.ones(len(t)),t]);b=np.linalg.lstsq(X,y,rcond=None)[0];obs.append(abs(corr(y-X@b,t,"spearman")))
 observed=float(np.nanmean(obs)) if obs else np.nan
 for _ in range(n):
  vals=[]
  for _,g in df.groupby(["dataset","program_id"]):
   h=_norm(g)
   if h is None:continue
   t=h.normalized_time.to_numpy(float);y=rng.permutation(h.activity.to_numpy(float));X=np.column_stack([np.ones(len(t)),t]);b=np.linalg.lstsq(X,y,rcond=None)[0];v=abs(corr(y-X@b,t,"spearman"));
   if np.isfinite(v):vals.append(v)
  null.append(float(np.mean(vals)) if vals else np.nan)
 null=np.asarray(null);p=float((1+np.sum(null>=observed))/(1+np.isfinite(null).sum())) if np.isfinite(observed) else np.nan
 return pd.DataFrame([{"observed_mean_abs_residual_time_spearman":observed,"permutation_n":n,"empirical_p":p,"null_mean":float(np.nanmean(null))}]),pd.DataFrame({"permutation":np.arange(n),"mean_abs_residual_time_spearman":null})
def run(permutations=500):
 log("starting residual-state analysis; LODO training-only time models")
 df=_load();_,s=_lodo(df);cross=_cross(df);perm,null=_perm(df,permutations);perm.to_csv(OUT/"04_residual_time_permutation_summary.csv",index=False);null.to_csv(OUT/"05_residual_time_permutation_null.csv",index=False)
 # Gate asks whether residual structure is reproducible across datasets. Time
 # permutation is secondary: residuals should not be significant against time.
 p=float(perm.iloc[0].empirical_p) if len(perm) else np.nan; cp=float(cross.residual_pearson.mean()) if len(cross) else np.nan
 overall=pd.DataFrame([{"n_datasets":int(df.dataset.nunique()),"n_programs":int(df.program_id.nunique()),"n_lodo_datasets":int(s.held_out_dataset.nunique()) if len(s) else 0,"mean_lodo_abs_residual_time_spearman":float(s.mean_abs_residual_time_spearman.mean()) if len(s) else np.nan,"mean_cross_dataset_residual_pearson":cp,"residual_time_permutation_p":p,"residual_state_supported":bool(len(cross)>=3 and np.isfinite(cp) and cp>0.3 and np.isfinite(p) and p>0.05)}])
 overall.to_csv(OUT/"06_STAGE2_9_16_SUMMARY.csv",index=False);log("complete");print("\nStage 2.9.16 residual LODO:");print(s.to_string(index=False) if len(s) else "none");print("\nStage 2.9.16 cross-dataset residual signal:");print(cross.groupby("program_id").residual_pearson.mean().to_string() if len(cross) else "none");print("\nStage 2.9.16 gate:");print(overall.to_string(index=False));return overall
if __name__=="__main__":run()
