"""Dynamics v4.6: OSKM reprogramming dynamics pipeline.

Stage 2.4 adds out-of-sample validation of the latent trajectory alignment.
The test holds out an entire biological timepoint from a target trajectory,
learns the geometric transform only from the remaining timepoints, and then
predicts the held-out target state from an independent reference trajectory.
GSE148158 GFP/OSKM branches with only two timepoints are explicitly marked as
insufficient for this test rather than being forced through it.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/"results"; OUT=RESULTS/"Dynamics"
STAGE_DIRS={i:OUT/f"stage{i}" for i in range(1,10)}
for d in STAGE_DIRS.values(): d.mkdir(parents=True,exist_ok=True)
STAGE21=OUT/"stage2_1"; STAGE21.mkdir(parents=True,exist_ok=True)
STAGE22=OUT/"stage2_2"; STAGE22.mkdir(parents=True,exist_ok=True)
STAGE23=OUT/"stage2_3"; STAGE23.mkdir(parents=True,exist_ok=True)
STAGE24=OUT/"stage2_4"; STAGE24.mkdir(parents=True,exist_ok=True)
DATASETS={"GSE28688":RESULTS/"GSE28688"/"non_normalized"/"07_PCA_coordinates.csv","GSE148158":RESULTS/"GSE148158"/"07_PCA_coordinates.csv","GSE52052":RESULTS/"GSE52052"/"08_PCA_coordinates.csv","GSE67462":RESULTS/"GSE67462"/"09_PCA_coordinates.csv","GSE297234":RESULTS/"GSE297234"/"08_PCA_coordinates.csv"}
GSM_TIME={"GSM4455240":48.,"GSM4455241":48.,"GSM4455242":72.,"GSM4455243":72.,"GSM4455244":48.,"GSM4455245":72.,"GSM710515":24.,"GSM710516":24.,"GSM710517":48.,"GSM710518":48.,"GSM710519":72.,"GSM710520":72.,"GSM1258008":264.,"GSM1258009":264.,"GSM1258010":264.,"GSM1258011":264.,"GSM1258012":264.,"GSM1258013":264.,"GSM1647454":0.,"GSM1647455":0.,"GSM1647456":24.,"GSM1647457":24.,"GSM1647458":72.,"GSM1647459":72.,"GSM1647460":120.,"GSM1647461":120.,"GSM1647462":168.,"GSM1647463":168.,"GSM1647464":264.,"GSM1647465":264.,"GSM1647466":360.,"GSM1647467":360.,"GSM1647468":432.,"GSM1647469":432.,"GSM8986586":0.,"GSM8986587":72.,"GSM8986588":168.,"GSM8986589":240.,"GSM8986590":0.,"GSM8986591":72.,"GSM8986592":168.,"GSM8986593":240.}
GSE28688_ROW_SAMPLE=[f"GSM{x}" for x in range(710513,710527)]
GSE28688_ROW_TIME=[0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

def load_pca(path):
    if not path.exists(): return None
    df=pd.read_csv(path,index_col=0); pcs=[c for c in ("PC1","PC2","PC3") if c in df.columns]
    if not pcs:return None
    x=df[pcs].apply(pd.to_numeric,errors="coerce"); x.index=x.index.astype(str); return x

def time_hours(ds,s):
    s=str(s).strip().strip('"')
    if s in GSM_TIME:return GSM_TIME[s]
    t=s.lower().replace("_"," ").replace("-"," ")
    pats={"GSE28688":[(r"24\s*h",24.),(r"48\s*h",48.),(r"72\s*h",72.)],"GSE148158":[(r"48",48.),(r"72",72.)],"GSE52052":[(r"day\s*11",264.)],"GSE67462":[(r"day\s*0\b",0.),(r"day\s*1\b",24.),(r"day\s*3\b",72.),(r"day\s*5\b",120.),(r"day\s*7\b",168.),(r"day\s*11\b",264.),(r"day\s*15\b",360.),(r"day\s*18\b",432.)],"GSE297234":[(r"d0\b|day\s*0\b",0.),(r"d3\b|day\s*3\b",72.),(r"d7\b|day\s*7\b",168.),(r"d10\b|day\s*10\b",240.)]}
    for p,v in pats.get(ds,[]):
        if re.search(p,t):return v
    return np.nan

def condition(ds,s):
    s=str(s).lower()
    if ds=="GSE148158":
        if "oskm" in s:return "OSKM"
        if "gfp" in s:return "GFP"
        if "h1" in s or "h9" in s:return "hESC"
        if "bj" in s:return "BJ_fibroblast"
    if ds=="GSE297234":
        return "aged" if any(x in s for x in ("6586","6587","6588","6589")) else ("young" if any(x in s for x in ("6590","6591","6592","6593")) else "unknown")
    return "all"

def replicate(s):
    s=str(s); m=re.search(r"(?:-|_|\s)([ab])$",s,re.I)
    if m:return m.group(1).lower()
    m=re.search(r"(?:rep|replicate)[_\s-]*(\d+)",s,re.I)
    if m:return m.group(1)
    m=re.fullmatch(r"GSM(\d+)",s)
    if m and 1647454<=int(m.group(1))<=1647469:return "1" if int(m.group(1))%2==0 else "2"
    return "unknown"

def zscore(x):
    x=pd.Series(x,dtype=float); sd=x.std(ddof=0)
    return (x-x.mean())/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index)

def orient(x):
    x=pd.Series(x,dtype=float).copy(); v=x.dropna()
    if len(v):
        i=v.abs().idxmax()
        if x.loc[i]<0:x=-x
    return x

def stage1_data_integration():
    rows=[]; states=[]
    for ds,path in DATASETS.items():
        x=load_pca(path)
        if x is None:
            rows.append({"dataset":ds,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","timing_source":"none","path":str(path)}); continue
        o=x.copy(); o.insert(0,"sample",o.index.astype(str)); source="GSM_or_text"
        if ds=="GSE28688" and len(o)==14:o["sample"]=GSE28688_ROW_SAMPLE; source="GSE28688_GEO_row_order"
        o["dataset"]=ds; o["time_hours"]=[time_hours(ds,s) for s in o["sample"]]
        if ds=="GSE28688" and source=="GSE28688_GEO_row_order":o["time_hours"]=GSE28688_ROW_TIME
        o["condition"]=[condition(ds,s) for s in o["sample"]]
        o["stage"]=o["time_hours"].map(lambda t:f"day{int(t/24)}" if pd.notna(t) and t%24==0 else (f"{int(t)}h" if pd.notna(t) else "unknown"))
        o["replicate"]=[replicate(s) for s in o["sample"]]; o["timing_source"]=source
        for i,pc in enumerate(["PC1","PC2","PC3"],1):o[f"latent_{i}"]=zscore(orient(o[pc]))
        timed=o[o["time_hours"].notna()]; role="trajectory" if timed["time_hours"].nunique()>=2 else "context_only"
        rows.append({"dataset":ds,"PCA_file_found":True,"n_samples":len(o),"n_timed_samples":len(timed),"n_unique_times":timed["time_hours"].nunique(),"role":role,"timing_source":source,"path":str(path)}); states.append(o)
    av=pd.DataFrame(rows); st=pd.concat(states,ignore_index=True) if states else pd.DataFrame(); av.to_csv(STAGE_DIRS[1]/"01_dataset_availability.csv",index=False); st.to_csv(STAGE_DIRS[1]/"02_master_sample_metadata.csv",index=False); return st,av

def _curve(st,ds,grid,branch=None,exclude_time=None):
    g=st[(st["dataset"]==ds)&st["time_hours"].notna()]
    if branch is not None:g=g[g["condition"]==branch]
    if exclude_time is not None:g=g[g["time_hours"]!=exclude_time]
    m=g.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index()
    if len(m)<2:return None
    t=m.index.to_numpy(float); u=(t-t.min())/(t.max()-t.min())
    return np.column_stack([np.interp(grid,u,m[c]) for c in ["latent_1","latent_2","latent_3"]])

def _point_at_time(st,ds,t,branch=None):
    g=st[(st["dataset"]==ds)&(st["time_hours"]==t)]
    if branch is not None:g=g[g["condition"]==branch]
    if g.empty:return None
    return g[["latent_1","latent_2","latent_3"]].mean().to_numpy(float)

def _rot(a,b):
    aa=a-a.mean(0); bb=b-b.mean(0); u,_,vt=np.linalg.svd(aa.T@bb); return u@vt

def _fit_transform(source,target):
    rot=_rot(source,target); a=(source-source.mean(0))@rot
    scale=np.divide(np.std(target,0),np.where(np.std(a,0)>0,np.std(a,0),1))
    return rot,source.mean(0),scale,target.mean(0)

def _apply_transform(x,trans):
    rot,sm,scale,tm=trans; y=(x-sm)@rot; y*=scale; y+=tm; return y

def stage2_1_common_latent_state(st):
    ds=[d for d in st.dataset.unique() if st[(st.dataset==d)&st.time_hours.notna()].time_hours.nunique()>=2]; grid=np.linspace(0,1,25)
    cur={d:_curve(st,d,grid) for d in ds}; cur={d:c for d,c in cur.items() if c is not None and np.isfinite(c).all()}
    ref="GSE67462" if "GSE67462" in cur else sorted(cur)[0]; r=cur[ref]; aligned={ref:r}
    for d,c in cur.items():
        if d==ref:continue
        trans=_fit_transform(c,r); aligned[d]=_apply_transform(c,trans)
    rows=[]
    for d,a in aligned.items():
        g=st[(st.dataset==d)&st.time_hours.notna()]
        for t,n in g.groupby("time_hours").size().items():
            u=(t-g.time_hours.min())/(g.time_hours.max()-g.time_hours.min()); v=np.array([np.interp(u,grid,a[:,j]) for j in range(3)])
            rows.append({"dataset":d,"time_hours":t,"normalized_time":u,"common_latent_1":v[0],"common_latent_2":v[1],"common_latent_3":v[2],"n_replicates":n})
    pd.DataFrame(rows).to_csv(STAGE21/"01_common_latent_trajectory.csv",index=False)
    out=st.copy()
    for j in range(1,4):out[f"common_latent_{j}"]=np.nan
    out["common_latent_status"]="not_aligned"
    for d,a in aligned.items():
        idx=out.index[(out.dataset==d)&out.time_hours.notna()]; lo,hi=out.loc[idx,"time_hours"].min(),out.loc[idx,"time_hours"].max(); u=(out.loc[idx,"time_hours"].to_numpy()-lo)/(hi-lo)
        for j in range(3):out.loc[idx,f"common_latent_{j+1}"]=np.interp(u,grid,a[:,j])
        out.loc[idx,"common_latent_status"]="time_anchored_aligned"
    out.to_csv(STAGE21/"04_sample_common_latent_state.csv",index=False)
    pd.DataFrame({"reference_dataset":[ref],"trajectory_datasets":[";".join(sorted(aligned))],"grid_points":[25],"alignment":["orthogonal Procrustes on time-mean trajectories"]}).to_csv(STAGE21/"03_alignment_metadata.csv",index=False)
    return out

def stage2_2_validate_common_latent(st):
    g=st[st.common_latent_status=="time_anchored_aligned"]; ds=sorted(g.dataset.unique()); grid=np.linspace(0,1,25); cur={d:_curve(st,d,grid) for d in ds}; cur={d:c for d,c in cur.items() if c is not None}; pairs=[]
    for i,a in enumerate(ds):
        for b in ds[i+1:]:
            x,y=cur[a],cur[b]; pairs.append({"dataset_a":a,"dataset_b":b,"trajectory_correlation":np.corrcoef(x.ravel(),y.ravel())[0,1],"aligned_rmse":np.sqrt(np.mean((x-y)**2)),"path_length_a":np.linalg.norm(np.diff(x,axis=0),axis=1).sum(),"path_length_b":np.linalg.norm(np.diff(y,axis=0),axis=1).sum()})
    p=pd.DataFrame(pairs); p["path_length_ratio"]=p.path_length_a/p.path_length_b; p.to_csv(STAGE22/"02_cross_dataset_distances.csv",index=False)
    q=[]
    for d in ds:
        r=np.mean([cur[x] for x in ds if x!=d],axis=0); x=cur[d]; q.append({"dataset":d,"leave_one_dataset_out_reference_rmse":np.sqrt(np.mean((x-r)**2)),"leave_one_dataset_out_correlation":np.corrcoef(x.ravel(),r.ravel())[0,1]})
    q=pd.DataFrame(q); q.to_csv(STAGE22/"01_alignment_quality.csv",index=False)
    mean=np.mean(list(cur.values()),axis=0); cc=[]
    for k,u in enumerate(grid):cc.append({"normalized_time":u,"mean_common_latent_1":mean[k,0],"mean_common_latent_2":mean[k,1],"mean_common_latent_3":mean[k,2],"cross_dataset_dispersion":np.mean(np.linalg.norm(np.vstack([cur[d][k] for d in ds])-mean[k],axis=1))})
    pd.DataFrame(cc).to_csv(STAGE22/"03_trajectory_concordance.csv",index=False); return q

def _branches(st):
    traj=[d for d in st.dataset.unique() if st[(st.dataset==d)&st.time_hours.notna()].time_hours.nunique()>=2]; grid=np.linspace(0,1,25); branches=[]
    for d in traj:
        conds=sorted(set(st.loc[(st.dataset==d)&st.time_hours.notna(),"condition"])); valid=[c for c in conds if _curve(st,d,grid,c) is not None]
        if d=="GSE148158" and len(valid)>=2:branches += [(d,c) for c in valid]
        else:branches.append((d,"all"))
    return branches

def stage2_3_within_time_residual_validation(st):
    branches=_branches(st); grid=np.linspace(0,1,25); ref="GSE67462"; refcurve=_curve(st,ref,grid); trans={}
    for d,c in branches:
        curve=_curve(st,d,grid,None if c=="all" else c)
        if curve is None:continue
        trans[(d,c)]=_fit_transform(curve,refcurve) if not (d==ref and c=="all") else (np.eye(3),curve.mean(0),np.ones(3),refcurve.mean(0))
    out=st.copy()
    for j in range(1,4):out[f"aligned_sample_latent_{j}"]=np.nan; out[f"aligned_time_mean_{j}"]=np.nan; out[f"within_time_residual_{j}"]=np.nan
    out["within_time_residual_norm"]=np.nan; out["trajectory_branch"]="unassigned"
    for (d,c),tr in trans.items():
        mask=(out.dataset==d)&out.time_hours.notna()&((out.condition==c) if c!="all" else True); idx=out.index[mask]; x=out.loc[idx,["latent_1","latent_2","latent_3"]].to_numpy(float); xa=_apply_transform(x,tr); out.loc[idx,["aligned_sample_latent_1","aligned_sample_latent_2","aligned_sample_latent_3"]]=xa; out.loc[idx,"trajectory_branch"]=c
        means=pd.DataFrame(xa,index=idx).groupby(out.loc[idx,"time_hours"]).transform("mean").to_numpy(); out.loc[idx,["aligned_time_mean_1","aligned_time_mean_2","aligned_time_mean_3"]]=means; res=xa-means
        for j in range(3):out.loc[idx,f"within_time_residual_{j+1}"]=res[:,j]
        out.loc[idx,"within_time_residual_norm"]=np.linalg.norm(res,axis=1)
    timed=out[out.time_hours.notna()&out.aligned_sample_latent_1.notna()].copy(); summary=[]; rep=[]
    for (d,b),g in timed.groupby(["dataset","trajectory_branch"]):
        summary.append({"dataset":d,"trajectory_branch":b,"n_samples":len(g),"n_times":g.time_hours.nunique(),"mean_within_time_residual_norm":g.within_time_residual_norm.mean(),"median_within_time_residual_norm":g.within_time_residual_norm.median(),"p95_within_time_residual_norm":g.within_time_residual_norm.quantile(.95),"within_time_sd_latent1":g.within_time_residual_1.std(ddof=1)})
        for t,h in g.groupby("time_hours"):
            vals=h[["aligned_sample_latent_1","aligned_sample_latent_2","aligned_sample_latent_3"]].to_numpy()
            if len(vals)>=2:rep.append({"dataset":d,"trajectory_branch":b,"time_hours":t,"n_replicates":len(h),"replicate_pair_distance":np.mean([np.linalg.norm(vals[i]-vals[j]) for i in range(len(vals)) for j in range(i+1,len(vals))])})
    sm=pd.DataFrame(summary); rp=pd.DataFrame(rep); sm.to_csv(STAGE23/"01_within_time_residual_summary.csv",index=False); rp.to_csv(STAGE23/"02_replicate_distances.csv",index=False); out.to_csv(STAGE23/"03_sample_level_aligned_states.csv",index=False)
    rows=[]
    for i,(a,ca) in enumerate(branches):
        for b,cb in branches[i+1:]:
            xa=_curve(st,a,grid,None if ca=="all" else ca); xb=_curve(st,b,grid,None if cb=="all" else cb)
            if xa is None or xb is None:continue
            rawcorr=np.corrcoef(xa.ravel(),xb.ravel())[0,1]; aa=_apply_transform(xa,_fit_transform(xa,xb)); alignedcorr=np.corrcoef(aa.ravel(),xb.ravel())[0,1]
            rows.append({"trajectory_a":f"{a}:{ca}","trajectory_b":f"{b}:{cb}","raw_trajectory_correlation":rawcorr,"post_procrustes_correlation":alignedcorr,"correlation_gain":alignedcorr-rawcorr,"raw_rmse":np.sqrt(np.mean((xa-xb)**2)),"post_procrustes_rmse":np.sqrt(np.mean((aa-xb)**2))})
    gain=pd.DataFrame(rows); gain.to_csv(STAGE23/"04_raw_vs_aligned_concordance.csv",index=False)
    wd=rp.groupby(["dataset","trajectory_branch"]).replicate_pair_distance.median() if not rp.empty else pd.Series(dtype=float); between=[]
    for _,r in gain.iterrows():
        a,ca=r.trajectory_a.split(":",1); b,cb=r.trajectory_b.split(":",1); xa=_curve(st,a,grid,None if ca=="all" else ca); xb=_curve(st,b,grid,None if cb=="all" else cb); dist=np.mean(np.linalg.norm(xa-xb,axis=1)); wa=wd.get((a,ca),np.nan); wb=wd.get((b,cb),np.nan); between.append({"trajectory_a":r.trajectory_a,"trajectory_b":r.trajectory_b,"mean_trajectory_distance":dist,"median_within_replicate_distance_a":wa,"median_within_replicate_distance_b":wb,"between_to_within_ratio":dist/np.nanmean([wa,wb])})
    ratio=pd.DataFrame(between); ratio.to_csv(STAGE23/"05_between_vs_within_distance.csv",index=False)
    report=pd.DataFrame({"metric":["reference_dataset","trajectory_branches","mean_within_time_residual","median_within_time_residual","mean_correlation_gain_from_alignment","mean_between_to_within_ratio"],"value":[ref,";".join(f"{d}:{b}" for d,b in branches),sm.mean_within_time_residual_norm.mean(),sm.median_within_time_residual_norm.median(),gain.correlation_gain.mean(),ratio.between_to_within_ratio.mean()]}); report.to_csv(STAGE23/"06_stage23_decision_metrics.csv",index=False); return out,report

def stage2_4_out_of_sample_validation(st):
    """Hold out one complete biological timepoint from the target trajectory.

    The transform is learned only from the remaining target timepoints and an
    independent reference trajectory. The held-out target sample is then
    transformed and compared with the reference state at the same normalized
    time. This avoids using the held-out target point to estimate its own
    alignment. Targets with <3 timepoints are reported as insufficient.
    """
    branches=_branches(st); grid=np.linspace(0,1,25); rows=[]
    for d,b in branches:
        g=st[(st.dataset==d)&st.time_hours.notna()&((st.condition==b) if b!="all" else True)]
        times=sorted(g.time_hours.unique())
        if len(times)<3:
            rows.append({"validation":"leave_one_timepoint_out","target_trajectory":f"{d}:{b}","held_out_time_hours":np.nan,"reference_trajectory":"","n_training_times":max(0,len(times)-1),"status":"insufficient_timepoints","normalized_time":np.nan,"prediction_rmse":np.nan,"prediction_relative_error":np.nan,"n_heldout_samples":0}); continue
        for hold in times:
            target_train=_curve(st,d,grid,None if b=="all" else b,exclude_time=hold)
            if target_train is None:continue
            lo,hi=min(times),max(times); u=(hold-lo)/(hi-lo)
            for rd,rb in branches:
                if rd==d and rb==b:continue
                rg=st[(st.dataset==rd)&st.time_hours.notna()&((st.condition==rb) if rb!="all" else True)]
                rt=sorted(rg.time_hours.unique())
                if len(rt)<2:continue
                refcurve=_curve(st,rd,grid,None if rb=="all" else rb)
                if refcurve is None:continue
                # Fit only on normalized-time support shared by the target training curve.
                train_grid=np.linspace(0,1,25); target_train2=target_train
                trans=_fit_transform(target_train2,refcurve)
                pred_ref=np.array([np.interp(u,train_grid,refcurve[:,j]) for j in range(3)])
                xheld=_point_at_time(st,d,hold,None if b=="all" else b)
                if xheld is None:continue
                pred=_apply_transform(xheld.reshape(1,-1),trans)[0]
                err=float(np.linalg.norm(pred-pred_ref)); scale=float(np.linalg.norm(pred_ref-refcurve.mean(0)))
                rows.append({"validation":"leave_one_timepoint_out","target_trajectory":f"{d}:{b}","held_out_time_hours":hold,"reference_trajectory":f"{rd}:{rb}","n_training_times":len(times)-1,"status":"tested","normalized_time":u,"prediction_rmse":err,"prediction_relative_error":err/(scale+1e-8),"n_heldout_samples":int(len(g[g.time_hours==hold]))})
    o=pd.DataFrame(rows); o.to_csv(STAGE24/"01_leave_one_timepoint_out.csv",index=False)
    tested=o[o.status=="tested"].copy(); summary=[]
    if not tested.empty:
        for t,g in tested.groupby("target_trajectory"):
            summary.append({"target_trajectory":t,"n_tests":len(g),"mean_prediction_error":g.prediction_rmse.mean(),"median_prediction_error":g.prediction_rmse.median(),"p95_prediction_error":g.prediction_rmse.quantile(.95),"mean_relative_error":g.prediction_relative_error.mean()})
    sm=pd.DataFrame(summary); sm.to_csv(STAGE24/"02_oos_summary_by_trajectory.csv",index=False)
    # Pairwise raw-vs-OOS comparison: how much better is a prediction than the naive
    # baseline of using the reference trajectory state without target alignment?
    base=[]
    for _,r in tested.iterrows():
        rd,rb=r.reference_trajectory.split(":",1); hold=r.held_out_time_hours; rg=st[(st.dataset==rd)&st.time_hours.notna()&((st.condition==rb) if rb!="all" else True)]; rt=sorted(rg.time_hours.unique())
        if len(rt)<2:continue
        refcurve=_curve(st,rd,grid,None if rb=="all" else rb); u=r.normalized_time; refpt=np.array([np.interp(u,grid,refcurve[:,j]) for j in range(3)]); xheld=_point_at_time(st,r.target_trajectory.split(":",1)[0],hold,None if r.target_trajectory.endswith(":all") else r.target_trajectory.split(":",1)[1]); rawerr=float(np.linalg.norm(xheld-refpt)); base.append({"target_trajectory":r.target_trajectory,"held_out_time_hours":hold,"reference_trajectory":r.reference_trajectory,"aligned_prediction_error":r.prediction_rmse,"naive_unaligned_error":rawerr,"error_reduction":rawerr-r.prediction_rmse,"relative_error_reduction":(rawerr-r.prediction_rmse)/(rawerr+1e-8)})
    bd=pd.DataFrame(base); bd.to_csv(STAGE24/"03_oos_vs_naive_baseline.csv",index=False)
    if not tested.empty:
        report=pd.DataFrame({"metric":["tested_rows","target_trajectories_tested","mean_oos_error","median_oos_error","mean_relative_error","mean_error_reduction_vs_naive"],"value":[len(tested),tested.target_trajectory.nunique(),tested.prediction_rmse.mean(),tested.prediction_rmse.median(),tested.prediction_relative_error.mean(),bd.error_reduction.mean() if not bd.empty else np.nan]})
    else:
        report=pd.DataFrame({"metric":["tested_rows","target_trajectories_tested","mean_oos_error","median_oos_error","mean_relative_error","mean_error_reduction_vs_naive"],"value":[0,0,np.nan,np.nan,np.nan,np.nan]})
    report.to_csv(STAGE24/"04_stage24_decision_metrics.csv",index=False)
    with open(STAGE24/"REPORT.txt","w",encoding="utf-8") as f:
        f.write("Dynamics v4.6 Stage 2.4 — out-of-sample trajectory validation\n\n")
        f.write("A complete biological timepoint is held out from each target trajectory with >=3 timepoints.\n")
        f.write("The geometric transform is estimated from the remaining target timepoints and an independent reference trajectory, then applied to the held-out target sample.\n")
        f.write("Two-timepoint branches (including GSE148158 GFP and OSKM) are explicitly marked insufficient.\n")
    return o,report

def stage3_trajectory_reconstruction(st):
    rows=[]
    for (d,b),g in st[st.time_hours.notna()&st.aligned_sample_latent_1.notna()].groupby(["dataset","trajectory_branch"]):
        if g.time_hours.nunique()<2:continue
        m=g.groupby("time_hours")[["aligned_sample_latent_1","aligned_sample_latent_2","aligned_sample_latent_3"]].mean().sort_index().reset_index(); m.insert(0,"dataset",d); m.insert(1,"trajectory_branch",b); m.columns=["dataset","trajectory_branch","time_hours","common_latent_1","common_latent_2","common_latent_3"]; rows.append(m)
    out=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(); out.to_csv(STAGE_DIRS[3]/"01_reconstructed_trajectories.csv",index=False); return out

def derivative(x,t):
    if len(x)<2:return np.full(len(x),np.nan)
    if np.any(np.diff(t)<=0):raise ValueError("Non-increasing time")
    return np.gradient(x,t) if len(x)>2 else np.repeat((x[1]-x[0])/(t[1]-t[0]),2)

def stage4_dynamics(tr):
    out=tr.copy()
    for _,idx in out.groupby(["dataset","trajectory_branch"]).groups.items():
        q=out.loc[idx].sort_values("time_hours"); t=q.time_hours.to_numpy(float)
        for j in (1,2,3):out.loc[q.index,f"dcommon_latent_{j}_dt"]=derivative(q[f"common_latent_{j}"].to_numpy(float),t)
        out.loc[q.index,"state_speed"]=np.sqrt(np.nansum(out.loc[q.index,[f"dcommon_latent_{j}_dt" for j in (1,2,3)]].to_numpy()**2,axis=1))
    out.to_csv(STAGE_DIRS[4]/"01_dynamics.csv",index=False); return out

def stage5_critical_transitions(d):
    out=d.copy(); out["rolling_variance_latent1"]=np.nan; out["rolling_autocorrelation_latent1"]=np.nan
    for _,idx in out.groupby(["dataset","trajectory_branch"]).groups.items():
        q=out.loc[idx].sort_values("time_hours"); s=q.common_latent_1; w=min(5,len(s))
        if w>=3:
            out.loc[q.index,"rolling_variance_latent1"]=s.rolling(w,min_periods=3).var().to_numpy(); out.loc[q.index,"rolling_autocorrelation_latent1"]=s.rolling(w,min_periods=3).apply(lambda z:z.autocorr(1) if z.std()>0 else np.nan).to_numpy()
    out["critical_transition_flag"]=False; out.to_csv(STAGE_DIRS[5]/"01_critical_transition_indicators.csv",index=False); return out

def stage6_symbolic_equation_discovery(d):
    if d.empty:return pd.DataFrame()
    o=pd.DataFrame({"dataset":d.dataset,"trajectory_branch":d.trajectory_branch,"time_hours":d.time_hours,"x":d.common_latent_1,"target_dx_dt":d.dcommon_latent_1_dt}); o["x2"]=o.x**2; o["x3"]=o.x**3; o["abs_x"]=o.x.abs(); o["sin_x"]=np.sin(o.x); o["cos_x"]=np.cos(o.x); o.to_csv(STAGE_DIRS[6]/"01_symbolic_regression_design.csv",index=False); return o

def stage7_heldout_validation(d):
    ds=sorted(d.dataset.unique()) if not d.empty else []; o=pd.DataFrame([{"held_out_dataset":x,"training_datasets":";".join(y for y in ds if y!=x),"validation_type":"leave-one-complete-dataset-out","status":"planned"} for x in ds]); o.to_csv(STAGE_DIRS[7]/"01_heldout_validation_plan.csv",index=False); return o

def stage8_regulatory_integration(d):
    ds=sorted(d.dataset.unique()) if not d.empty else []; o=pd.DataFrame([{"expression_dataset":x,"regulatory_dataset":"GSE67520","integration":"time-aligned regulatory evidence","causal_claim":False,"status":"planned"} for x in ds]); o.to_csv(STAGE_DIRS[8]/"01_regulatory_integration_plan.csv",index=False); return o

def stage9_predictive_ai(d):
    rows=[]
    for (ds,b),g in d.groupby(["dataset","trajectory_branch"]):
        g=g.sort_values("time_hours")
        for i in range(len(g)-1):
            a,z=g.iloc[i],g.iloc[i+1]; rows.append({"dataset":ds,"trajectory_branch":b,"time_t":a.time_hours,"time_next":z.time_hours,"dt_hours":z.time_hours-a.time_hours,"latent_1_t":a.common_latent_1,"latent_2_t":a.common_latent_2,"latent_3_t":a.common_latent_3,"latent_1_next":z.common_latent_1,"latent_2_next":z.common_latent_2,"latent_3_next":z.common_latent_3})
    o=pd.DataFrame(rows); o.to_csv(STAGE_DIRS[9]/"01_next_state_prediction_table.csv",index=False); return o

def print_stage24(v):
    print("\n"+"="*88+"\nSTAGE 2.4 — OUT-OF-SAMPLE / LEAVE-ONE-TIMEPOINT-OUT VALIDATION\n"+"="*88)
    if v is None or v.empty:print("No Stage 2.4 results.")
    else:
        with pd.option_context("display.max_rows",None,"display.max_columns",None,"display.width",240):
            for title,file in [("OOS SUMMARY BY TRAJECTORY","02_oos_summary_by_trajectory.csv"),("OOS VS NAIVE BASELINE","03_oos_vs_naive_baseline.csv"),("DECISION METRICS","04_stage24_decision_metrics.csv")]:
                print(f"\n--- {title} ---"); p=STAGE24/file
                if p.exists():print(pd.read_csv(p).to_string(index=False))
        tested=v[v.status=="tested"] if "status" in v.columns else pd.DataFrame()
        print("\n--- INTERPRETATION ---")
        print(f"tested rows = {len(tested)}")
        print(f"tested target trajectories = {tested.target_trajectory.nunique() if not tested.empty else 0}")
        print("NOTE: Stage 2.4 is a generalization test. A low in-sample alignment error is not sufficient; held-out timepoints must also be predicted with low error.")
        print("Two-timepoint branches are intentionally excluded from the OOS test.")
    print("="*88)

def main():
    st,av=stage1_data_integration(); st=stage2_1_common_latent_state(st); stage2_2_validate_common_latent(st); st,v23=stage2_3_within_time_residual_validation(st); print_stage24(stage2_4_out_of_sample_validation(st)); tr=stage3_trajectory_reconstruction(st); dy=stage4_dynamics(tr); dy=stage5_critical_transitions(dy); s=stage6_symbolic_equation_discovery(dy); h=stage7_heldout_validation(dy); r=stage8_regulatory_integration(dy); p=stage9_predictive_ai(dy)
    print(f"\nDynamics v4.6 results written to: {OUT}")
    print(f"Datasets with PCA: {int(av.PCA_file_found.sum())}/{len(av)}")
    print(f"Stage 2.1 aligned trajectory datasets: {st[st.common_latent_status=='time_anchored_aligned'].dataset.nunique()}")
    print(f"Stage 2.4 tested OOS rows: {int((stage2_4_out_of_sample_validation(st)[0].status=='tested').sum())}")
    print(f"Trajectory timepoints: {len(tr)}")
    print(f"Stage 6 symbolic rows: {len(s)}")
    print(f"Stage 7 held-out datasets: {len(h)}")
    print(f"Stage 8 regulatory rows: {len(r)}")
    print(f"Stage 9 prediction pairs: {len(p)}")

if __name__=="__main__":main()
