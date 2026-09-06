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
    while re.match(r"^GSE\d+__",s,re.I): s=re.sub(r"^GSE\d+__","",s,count=1,flags=re.I)
    return s

def _candidate_column_names(dataset,sample):
    ds,sm=str(dataset).strip().strip('"'),str(sample).strip().strip('"'); raw=_strip_dataset_prefix(sm)
    return [sm,raw,f"{ds}__{raw}",f"{ds}_{raw}",f"{ds}/{raw}"]

def _build_matrix_column_map(matrix,metadata):
    actual=list(matrix.columns); normalised={}; ambiguous=set()
    for col in actual:
        key=_normalise_name(col)
        if key in normalised and normalised[key]!=col: ambiguous.add(key)
        else: normalised[key]=col
    metadata=metadata.copy(); resolved=[]
    for _,row in metadata.iterrows():
        found=None
        for candidate in _candidate_column_names(row["dataset"],row["sample"]):
            key=_normalise_name(candidate)
            if key in normalised and key not in ambiguous: found=normalised[key]; break
        resolved.append(found)
    metadata["matrix_column"]=resolved
    matched=metadata["matrix_column"].notna().sum()
    # Legacy Stage 2.6 fallback: an older row-wise z-score operation could
    # serialize sample columns as 0,1,2,... while metadata retained the same
    # sample order. If every metadata row corresponds to one matrix column,
    # positional matching is deterministic and does not infer biology.
    if matched==0 and len(actual)==len(metadata):
        numeric_like=all(re.fullmatch(r"\d+(?:\.0+)?",str(c).strip()) for c in actual)
        if numeric_like:
            metadata["matrix_column"]=actual
            print(f"Stage 2.7: recovered {len(actual)} legacy matrix columns by metadata order.")
    return metadata

def _time_from_text(dataset,sample):
    s=_strip_dataset_prefix(sample).lower().replace("_"," ").replace("-"," ")
    patterns={
        "GSE148158":[(r"48\s*h|48h|day\s*2",48.),(r"72\s*h|72h|day\s*3",72.)],
        "GSE52052":[(r"day\s*11|11\s*d|11d",264.)],
        "GSE67462":[(r"day\s*0\b|d\s*0\b|0\s*h",0.),(r"day\s*1\b|d\s*1\b|24\s*h|24h",24.),(r"day\s*3\b|d\s*3\b|72\s*h|72h",72.),(r"day\s*5\b|d\s*5\b|120\s*h|120h",120.),(r"day\s*7\b|d\s*7\b|168\s*h|168h",168.),(r"day\s*11\b|d\s*11\b|264\s*h|264h",264.),(r"day\s*15\b|d\s*15\b|360\s*h|360h",360.),(r"day\s*18\b|d\s*18\b|432\s*h|432h",432.)],
        "GSE297234":[(r"d\s*0\b|day\s*0\b",0.),(r"d\s*3\b|day\s*3\b",72.),(r"d\s*7\b|day\s*7\b",168.),(r"d\s*10\b|day\s*10\b",240.)],
        "GSE28688":[(r"24\s*h|24h",24.),(r"48\s*h|48h",48.),(r"72\s*h|72h",72.)],
    }
    for p,v in patterns.get(dataset,[]):
        if re.search(p,s): return v
    return np.nan

def _time_hours_for_validation(dataset,sample,row_index=None,gsm_time=None):
    raw=_strip_dataset_prefix(sample); gsm=re.search(r"GSM(\d+)",raw,re.I)
    if gsm_time is not None and gsm:
        key=f"GSM{gsm.group(1)}"
        if key in gsm_time: return float(gsm_time[key])
    if dataset=="GSE28688" and row_index is not None:
        times=[0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
        if 0<=row_index<len(times) and np.isfinite(times[row_index]): return times[row_index]
    return _time_from_text(dataset,raw)

def _recover_ordered_sample_labels(metadata):
    try:
        from Dynamics import PCA_FILES,GSE28688_ROW_SAMPLE
    except Exception:
        return metadata
    metadata=metadata.copy(); recovered=0
    for ds,idxs in metadata.groupby("dataset",sort=False).groups.items():
        if ds not in PCA_FILES or not PCA_FILES[ds].exists(): continue
        try:
            pca=pd.read_csv(PCA_FILES[ds],index_col=0); labels=[str(x) for x in pca.index]
            if ds=="GSE28688" and len(labels)==14: labels=list(GSE28688_ROW_SAMPLE)
        except Exception: continue
        idxs=list(idxs)
        for pos,idx in enumerate(idxs):
            sample=str(metadata.at[idx,"sample"]); raw=_strip_dataset_prefix(sample)
            if re.fullmatch(r"\d+",raw) and pos<len(labels): metadata.at[idx,"sample"]=labels[pos]; recovered+=1
    if recovered: print(f"Stage 2.7: recovered {recovered} legacy sample labels from PCA/GEO order.")
    return metadata

def _print_mapping_diagnostics(matrix,metadata):
    rows=[]
    for ds,g in metadata.groupby("dataset",sort=True):
        timed=g[g.time_hours.notna()]; matched=g[g.matrix_column.notna()]; tm=timed[timed.matrix_column.notna()]
        unresolved=g[g.time_hours.isna()]["sample"].astype(str).tolist()[:5]
        rows.append({"dataset":ds,"matrix_columns":int(g.matrix_column.nunique()),"metadata_samples":len(g),"matched_samples":len(matched),"timed_samples":len(timed),"timed_matched":len(tm),"unique_times":int(tm.time_hours.nunique()),"time_values":",".join(map(str,sorted(tm.time_hours.unique()))),"unresolved_time_examples":";".join(unresolved),"replicates":",".join(sorted(map(str,tm.replicate.dropna().unique())))})
    diag=pd.DataFrame(rows); diag.to_csv(OUT/"00_mapping_diagnostics.csv",index=False); print("\nStage 2.7 sample-to-matrix mapping:"); print(diag.to_string(index=False)); return diag

def _load_common_space():
    mp=STAGE26/"06_common_human_gene_matrix.csv"; pp=STAGE26/"07_common_gene_sample_metadata.csv"
    if not mp.exists() or not pp.exists(): raise FileNotFoundError("Stage 2.6 outputs are missing; run Dynamics.py first.")
    matrix=pd.read_csv(mp,index_col=0).apply(pd.to_numeric,errors="coerce"); metadata=pd.read_csv(pp)
    metadata["dataset"]=metadata["dataset"].astype(str); metadata["sample"]=metadata["sample"].astype(str)
    metadata=_recover_ordered_sample_labels(metadata)
    from Dynamics import condition,replicate,GSM_TIME
    times=[]; conditions=[]; replicates=[]; row_counts={}
    for ds,sample in zip(metadata.dataset,metadata["sample"]):
        raw=_strip_dataset_prefix(sample); idx=row_counts.get(ds,0); row_counts[ds]=idx+1
        times.append(_time_hours_for_validation(ds,raw,idx if ds=="GSE28688" else None,GSM_TIME)); conditions.append(condition(ds,raw)); replicates.append(replicate(raw))
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
