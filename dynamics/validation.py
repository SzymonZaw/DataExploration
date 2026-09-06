"""Stage 2.7 validation of the biological common gene space."""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
STAGE26=ROOT/"results"/"Dynamics"/"stage2_6"
OUT=ROOT/"results"/"Dynamics"/"stage2_7"
OUT.mkdir(parents=True,exist_ok=True)


def _metrics(y_true,y_pred):
    a=np.asarray(y_true,dtype=float).ravel(); b=np.asarray(y_pred,dtype=float).ravel(); ok=np.isfinite(a)&np.isfinite(b)
    if not ok.any(): return {"rmse":np.nan,"mae":np.nan,"correlation":np.nan}
    a,b=a[ok],b[ok]; corr=np.corrcoef(a,b)[0,1] if len(a)>1 and np.std(a)>0 and np.std(b)>0 else np.nan
    return {"rmse":float(np.sqrt(np.mean((a-b)**2))),"mae":float(np.mean(np.abs(a-b))),"correlation":float(corr)}


def _interpolate(points,times,target):
    times=np.asarray(times,dtype=float); points=np.asarray(points,dtype=float)
    if len(times)<2 or target<times.min() or target>times.max(): return None
    return np.asarray([np.interp(target,times,points[:,j]) for j in range(points.shape[1])])


def _normalise_name(value):
    return re.sub(r"\s+","",str(value).strip().strip('"').replace("\\","/")).lower()


def _strip_dataset_prefix(sample):
    s=str(sample).strip().strip('"')
    while "__" in s and re.match(r"^GSE\d+__",s,re.I): s=s.split("__",1)[1]
    return s


def _candidate_column_names(dataset,sample):
    ds,sm=str(dataset).strip().strip('"'),str(sample).strip().strip('"')
    return [sm,f"{ds}__{sm}",f"{ds}_{sm}",f"{ds}/{sm}"]


def _build_matrix_column_map(matrix,metadata):
    actual=list(matrix.columns); normalised={}; ambiguous=set()
    for col in actual:
        key=_normalise_name(col)
        if key in normalised and normalised[key]!=col: ambiguous.add(key)
        else: normalised[key]=col
    resolved=[]
    for _,row in metadata.iterrows():
        found=None
        for candidate in _candidate_column_names(row["dataset"],row["sample"]):
            key=_normalise_name(candidate)
            if key in normalised and key not in ambiguous: found=normalised[key]; break
        resolved.append(found)
    metadata=metadata.copy(); metadata["matrix_column"]=resolved
    return metadata


def _time_hours_for_validation(dataset,sample,row_index=None):
    """Resolve times deterministically from the dataset/GSM labels used by Stage 2.6."""
    raw=_strip_dataset_prefix(sample); gsm=re.search(r"GSM(\d+)",raw,re.I)
    n=int(gsm.group(1)) if gsm else None
    if dataset=="GSE148158":
        return {4455240:48.,4455241:48.,4455242:72.,4455243:72.,4455244:48.,4455245:72.}.get(n,np.nan)
    if dataset=="GSE28688":
        times=[0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
        return times[row_index] if row_index is not None and 0<=row_index<len(times) else np.nan
    if dataset=="GSE52052":
        return 264. if n in {1258008,1258009,1258010,1258011,1258012,1258013} else np.nan
    if dataset=="GSE67462":
        groups={0.:range(1647454,1647456),24.:range(1647456,1647458),72.:range(1647458,1647460),120.:range(1647460,1647462),168.:range(1647462,1647464),264.:range(1647464,1647466),360.:range(1647466,1647468),432.:range(1647468,1647470)}
        for t,ids in groups.items():
            if n in ids:return t
        return np.nan
    if dataset=="GSE297234":
        return {8986586:0.,8986587:72.,8986588:168.,8986589:240.,8986590:0.,8986591:72.,8986592:168.,8986593:240.}.get(n,np.nan)
    return np.nan


def _print_mapping_diagnostics(matrix,metadata):
    rows=[]
    for ds,g in metadata.groupby("dataset",sort=True):
        timed=g[g.time_hours.notna()]; matched=g[g.matrix_column.notna()]; tm=timed[timed.matrix_column.notna()]
        rows.append({"dataset":ds,"matrix_columns":int(tm.matrix_column.nunique()),"metadata_samples":len(g),"matched_samples":len(matched),"timed_samples":len(timed),"timed_matched":len(tm),"unique_times":int(tm.time_hours.nunique()),"time_values":",".join(map(str,sorted(tm.time_hours.unique()))),"replicates":",".join(sorted(map(str,tm.replicate.dropna().unique())))})
    diag=pd.DataFrame(rows); diag.to_csv(OUT/"00_mapping_diagnostics.csv",index=False); print("\nStage 2.7 sample-to-matrix mapping:"); print(diag.to_string(index=False)); return diag


def _load_common_space():
    mp=STAGE26/"06_common_human_gene_matrix.csv"; pp=STAGE26/"07_common_gene_sample_metadata.csv"
    if not mp.exists() or not pp.exists(): raise FileNotFoundError("Stage 2.6 outputs are missing; run Dynamics.py first.")
    matrix=pd.read_csv(mp,index_col=0).apply(pd.to_numeric,errors="coerce"); metadata=pd.read_csv(pp)
    metadata["dataset"]=metadata["dataset"].astype(str); metadata["sample"]=metadata["sample"].astype(str)
    from Dynamics import condition,replicate,GSE28688_ROW_SAMPLE
    times=[]; conditions=[]; replicates=[]
    for ds,sample in zip(metadata.dataset,metadata['sample']):
        raw=_strip_dataset_prefix(sample); idx=None
        if ds=="GSE28688":
            m=re.fullmatch(r"GSM(\d+)",raw,re.I)
            if m: idx=int(m.group(1))-710513
        times.append(_time_hours_for_validation(ds,raw,idx)); conditions.append(condition(ds,raw)); replicates.append(replicate(raw))
    metadata["time_hours"]=times; metadata["condition"]=conditions; metadata["replicate"]=replicates
    metadata=_build_matrix_column_map(matrix,metadata); _print_mapping_diagnostics(matrix,metadata)
    if metadata.matrix_column.notna().sum()==0: raise RuntimeError("Stage 2.7 could not match metadata to Stage 2.6 matrix columns.")
    return matrix,metadata[metadata.matrix_column.notna()].copy()


def _trajectories(matrix,metadata,time_override=None):
    out={}
    for ds,g in metadata.groupby("dataset"):
        g=g[g.matrix_column.notna()&g.time_hours.notna()].copy()
        if len(g)<2: continue
        times=g.time_hours.astype(float).to_numpy()
        if time_override: times=np.asarray([time_override.get((ds,c),t) for c,t in zip(g.matrix_column,times)],dtype=float)
        frame=pd.DataFrame(matrix[g.matrix_column].T.to_numpy(),index=times,columns=matrix.index).groupby(level=0).mean().sort_index()
        if len(frame)>=2: out[ds]=(frame.index.to_numpy(float),frame.to_numpy(float))
    return out


def leave_one_dataset_out(matrix,metadata):
    traj=_trajectories(matrix,metadata); rows=[]
    for test_ds,(test_times,_) in sorted(traj.items()):
        for target in test_times:
            preds=[]
            for train_ds,(times,values) in traj.items():
                if train_ds==test_ds: continue
                p=_interpolate(values,times,target)
                if p is not None: preds.append(p)
            if not preds: continue
            cols=metadata[(metadata.dataset==test_ds)&(metadata.time_hours==target)].matrix_column.tolist()
            if cols: rows.append({"validation":"leave_one_dataset_out","test_dataset":test_ds,"time_hours":float(target),"n_training_datasets":len(preds),**_metrics(matrix[cols].mean(axis=1).to_numpy(float),np.mean(preds,axis=0))})
    return pd.DataFrame(rows)


def leave_one_replicate_out(matrix,metadata):
    rows=[]
    for ds,g in metadata.groupby("dataset"):
        reps=sorted(r for r in g.replicate.dropna().unique() if str(r).lower()!="unknown")
        if len(reps)<2: continue
        for held in reps:
            train=g[(g.replicate!=held)&g.matrix_column.notna()&g.time_hours.notna()]; test=g[(g.replicate==held)&g.matrix_column.notna()&g.time_hours.notna()]
            frame=pd.DataFrame(matrix[train.matrix_column].T.to_numpy(),index=train.time_hours.to_numpy(),columns=matrix.index).groupby(level=0).mean().sort_index()
            for target in sorted(test.time_hours.unique()):
                pred=_interpolate(frame.to_numpy(),frame.index.to_numpy(float),float(target))
                if pred is not None: rows.append({"validation":"leave_one_replicate_out","dataset":ds,"held_out_replicate":str(held),"time_hours":float(target),**_metrics(matrix[test[test.time_hours==target].matrix_column].mean(axis=1).to_numpy(float),pred)})
    return pd.DataFrame(rows)


def permutation_null(matrix,metadata,n_permutations=25,seed=42):
    rng=np.random.default_rng(seed); rows=[]
    for permutation in range(n_permutations):
        override={}
        for ds,g in metadata.groupby("dataset"):
            g=g[g.time_hours.notna()&g.matrix_column.notna()]; vals=g.time_hours.to_numpy(float)
            for col,value in zip(g.matrix_column,rng.permutation(vals)): override[(ds,col)]=float(value)
        traj=_trajectories(matrix,metadata,override)
        for test_ds in sorted(traj):
            test=metadata[(metadata.dataset==test_ds)&metadata.time_hours.notna()]
            for target in sorted(test.time_hours.unique()):
                preds=[_interpolate(values,times,float(target)) for train_ds,(times,values) in traj.items() if train_ds!=test_ds and _interpolate(values,times,float(target)) is not None]
                if preds:
                    cols=test[test.time_hours==target].matrix_column.tolist(); rows.append({"validation":"time_permutation_null","permutation":permutation,"test_dataset":test_ds,"time_hours":float(target),**_metrics(matrix[cols].mean(axis=1).to_numpy(float),np.mean(preds,axis=0))})
    return pd.DataFrame(rows)


def stage2_7(n_permutations=25,seed=42):
    matrix,metadata=_load_common_space(); dataset_df=leave_one_dataset_out(matrix,metadata); replicate_df=leave_one_replicate_out(matrix,metadata); null_df=permutation_null(matrix,metadata,n_permutations,seed)
    dataset_df.to_csv(OUT/"01_leave_one_dataset_out.csv",index=False); replicate_df.to_csv(OUT/"02_leave_one_replicate_out.csv",index=False); null_df.to_csv(OUT/"03_time_permutation_null.csv",index=False)
    summary=[]
    for name,frame in (("leave_one_dataset_out",dataset_df),("leave_one_replicate_out",replicate_df),("time_permutation_null",null_df)):
        summary.append({"validation":name,"n_cases":len(frame),"mean_rmse":frame.rmse.mean() if not frame.empty else np.nan,"median_rmse":frame.rmse.median() if not frame.empty else np.nan,"mean_mae":frame.mae.mean() if not frame.empty else np.nan,"mean_correlation":frame.correlation.mean() if not frame.empty else np.nan})
    summary_df=pd.DataFrame(summary); summary_df.to_csv(OUT/"04_validation_summary.csv",index=False); return summary_df

__all__=["stage2_7","leave_one_dataset_out","leave_one_replicate_out","permutation_null"]
