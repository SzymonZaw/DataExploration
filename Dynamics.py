"""Dynamics v4.8: OSKM reprogramming dynamics pipeline.

Stage 2.5 adds a feature-level biological harmonization audit. The existing
PCA/Procrustes representation remains an experimental baseline. Stage 2.5
searches existing lightweight expression outputs for comparable feature/sample
matrices, reports feature overlap and data modality, and constructs a
provisional shared feature matrix only when a sufficiently compatible set is
available. It deliberately does not treat time-anchored PCA coordinates as
biological equivalence.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent; RESULTS=ROOT/"results"; OUT=RESULTS/"Dynamics"
for i in range(1,10): (OUT/f"stage{i}").mkdir(parents=True,exist_ok=True)
STAGE21=OUT/"stage2_1"; STAGE22=OUT/"stage2_2"; STAGE23=OUT/"stage2_3"; STAGE24=OUT/"stage2_4"; STAGE25=OUT/"stage2_5"
for d in (STAGE21,STAGE22,STAGE23,STAGE24,STAGE25): d.mkdir(parents=True,exist_ok=True)
DATASETS={"GSE28688":RESULTS/"GSE28688"/"non_normalized"/"07_PCA_coordinates.csv","GSE148158":RESULTS/"GSE148158"/"07_PCA_coordinates.csv","GSE52052":RESULTS/"GSE52052"/"08_PCA_coordinates.csv","GSE67462":RESULTS/"GSE67462"/"09_PCA_coordinates.csv","GSE297234":RESULTS/"GSE297234"/"08_PCA_coordinates.csv"}
GSM_TIME={"GSM4455240":48.,"GSM4455241":48.,"GSM4455242":72.,"GSM4455243":72.,"GSM4455244":48.,"GSM4455245":72.,"GSM710515":24.,"GSM710516":24.,"GSM710517":48.,"GSM710518":48.,"GSM710519":72.,"GSM710520":72.,"GSM1258008":264.,"GSM1258009":264.,"GSM1258010":264.,"GSM1258011":264.,"GSM1258012":264.,"GSM1258013":264.,"GSM1647454":0.,"GSM1647455":0.,"GSM1647456":24.,"GSM1647457":24.,"GSM1647458":72.,"GSM1647459":72.,"GSM1647460":120.,"GSM1647461":120.,"GSM1647462":168.,"GSM1647463":168.,"GSM1647464":264.,"GSM1647465":264.,"GSM1647466":360.,"GSM1647467":360.,"GSM1647468":432.,"GSM1647469":432.,"GSM8986586":0.,"GSM8986587":72.,"GSM8986588":168.,"GSM8986589":240.,"GSM8986590":0.,"GSM8986591":72.,"GSM8986592":168.,"GSM8986593":240.}
GSE28688_ROW_SAMPLE=[f"GSM{x}" for x in range(710513,710527)]; GSE28688_ROW_TIME=[0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
def load_pca(path):
    if not path.exists(): return None
    df=pd.read_csv(path,index_col=0); pcs=[c for c in ("PC1","PC2","PC3") if c in df.columns]
    if len(pcs)<3:return None
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
    if ds=="GSE297234":return "aged" if any(x in s for x in ("6586","6587","6588","6589")) else ("young" if any(x in s for x in ("6590","6591","6592","6593")) else "unknown")
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
    x=pd.Series(x,dtype=float); sd=x.std(ddof=0); return (x-x.mean())/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index)
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
        if x is None:rows.append({"dataset":ds,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","timing_source":"none","path":str(path)}); continue
        o=x.copy(); o.insert(0,"sample",o.index.astype(str)); source="GSM_or_text"
        if ds=="GSE28688" and len(o)==14:o["sample"]=GSE28688_ROW_SAMPLE; source="GSE28688_GEO_row_order"
        o["dataset"]=ds; o["time_hours"]=[time_hours(ds,s) for s in o["sample"]]
        if ds=="GSE28688" and source=="GSE28688_GEO_row_order":o["time_hours"]=GSE28688_ROW_TIME
        o["condition"]=[condition(ds,s) for s in o["sample"]]; o["stage"]=o["time_hours"].map(lambda t:f"day{int(t/24)}" if pd.notna(t) and t%24==0 else (f"{int(t)}h" if pd.notna(t) else "unknown")); o["replicate"]=[replicate(s) for s in o["sample"]]; o["timing_source"]=source
        for i,pc in enumerate(["PC1","PC2","PC3"],1):o[f"latent_{i}"]=zscore(orient(o[pc]))
        timed=o[o.time_hours.notna()]; role="trajectory" if timed.time_hours.nunique()>=2 else "context_only"; rows.append({"dataset":ds,"PCA_file_found":True,"n_samples":len(o),"n_timed_samples":len(timed),"n_unique_times":timed.time_hours.nunique(),"role":role,"timing_source":source,"path":str(path)}); states.append(o)
    av=pd.DataFrame(rows); st=pd.concat(states,ignore_index=True) if states else pd.DataFrame(); av.to_csv(OUT/"stage1"/"01_dataset_availability.csv",index=False); st.to_csv(OUT/"stage1"/"02_master_sample_metadata.csv",index=False); return st,av
def _curve(st,ds,grid,branch=None,exclude_time=None):
    g=st[(st.dataset==ds)&st.time_hours.notna()]
    if branch is not None:g=g[g.condition==branch]
    if exclude_time is not None:g=g[g.time_hours!=exclude_time]
    m=g.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index()
    if len(m)<2:return None
    t=m.index.to_numpy(float); u=(t-t.min())/(t.max()-t.min()); return np.column_stack([np.interp(grid,u,m[c]) for c in ["latent_1","latent_2","latent_3"]])
def _point_at_time(st,ds,t,branch=None):
    g=st[(st.dataset==ds)&(st.time_hours==t)]
    if branch is not None:g=g[g.condition==branch]
    if g.empty:return None
    return g[["latent_1","latent_2","latent_3"]].mean().to_numpy(float)
def _rot(a,b):
    aa=a-a.mean(0); bb=b-b.mean(0); u,_,vt=np.linalg.svd(aa.T@bb); return u@vt
def _fit_transform(source,target):
    rot=_rot(source,target); a=(source-source.mean(0))@rot; scale=np.divide(np.std(target,0),np.where(np.std(a,0)>0,np.std(a,0),1)); return rot,source.mean(0),scale,target.mean(0)
def _apply_transform(x,trans):
    rot,sm,scale,tm=trans; y=(x-sm)@rot; y*=scale; y+=tm; return y
def stage2_1_common_latent_state(st):
    ds=[d for d in st.dataset.unique() if st[(st.dataset==d)&st.time_hours.notna()].time_hours.nunique()>=2]; grid=np.linspace(0,1,25); cur={d:_curve(st,d,grid) for d in ds}; cur={d:c for d,c in cur.items() if c is not None and np.isfinite(c).all()}; ref="GSE67462" if "GSE67462" in cur else sorted(cur)[0]; r=cur[ref]; aligned={ref:r}
    for d,c in cur.items():
        if d!=ref:aligned[d]=_apply_transform(c,_fit_transform(c,r))
    out=st.copy(); out["common_latent_1"]=np.nan; out["common_latent_2"]=np.nan; out["common_latent_3"]=np.nan; out["common_latent_status"]="not_aligned"
    for d,a in aligned.items():
        idx=out.index[(out.dataset==d)&out.time_hours.notna()]; lo,hi=out.loc[idx,"time_hours"].min(),out.loc[idx,"time_hours"].max(); u=(out.loc[idx,"time_hours"].to_numpy()-lo)/(hi-lo)
        for j in range(3):out.loc[idx,f"common_latent_{j+1}"]=np.interp(u,grid,a[:,j])
        out.loc[idx,"common_latent_status"]="time_anchored_aligned"
    out.to_csv(STAGE21/"04_sample_common_latent_state.csv",index=False); return out
def stage2_2_validate_common_latent(st):
    g=st[st.common_latent_status=="time_anchored_aligned"]; ds=sorted(g.dataset.unique()); grid=np.linspace(0,1,25); cur={d:_curve(st,d,grid) for d in ds}; cur={d:c for d,c in cur.items() if c is not None}; rows=[]
    for i,a in enumerate(ds):
        for b in ds[i+1:]:
            x,y=cur[a],cur[b]; rows.append({"dataset_a":a,"dataset_b":b,"trajectory_correlation":np.corrcoef(x.ravel(),y.ravel())[0,1],"aligned_rmse":np.sqrt(np.mean((x-y)**2)),"path_length_a":np.linalg.norm(np.diff(x,axis=0),axis=1).sum(),"path_length_b":np.linalg.norm(np.diff(y,axis=0),axis=1).sum()})
    p=pd.DataFrame(rows); p["path_length_ratio"]=p.path_length_a/p.path_length_b; p.to_csv(STAGE22/"02_cross_dataset_distances.csv",index=False); return p
def _branches(st):
    traj=[d for d in st.dataset.unique() if st[(st.dataset==d)&st.time_hours.notna()].time_hours.nunique()>=2]; grid=np.linspace(0,1,25); out=[]
    for d in traj:
        conds=sorted(set(st.loc[(st.dataset==d)&st.time_hours.notna(),"condition"])); valid=[c for c in conds if _curve(st,d,grid,c) is not None]
        if d=="GSE148158" and len(valid)>=2:out += [(d,c) for c in valid]
        else:out.append((d,"all"))
    return out
def stage2_3_within_time_residual_validation(st):
    branches=_branches(st); grid=np.linspace(0,1,25); refcurve=_curve(st,"GSE67462",grid); trans={}
    for d,c in branches:
        curve=_curve(st,d,grid,None if c=="all" else c)
        if curve is not None:trans[(d,c)]=_fit_transform(curve,refcurve) if not (d=="GSE67462" and c=="all") else (np.eye(3),curve.mean(0),np.ones(3),refcurve.mean(0))
    out=st.copy(); out["aligned_sample_latent_1"]=np.nan; out["aligned_sample_latent_2"]=np.nan; out["aligned_sample_latent_3"]=np.nan; out["within_time_residual_norm"]=np.nan; out["trajectory_branch"]="unassigned"
    for (d,c),tr in trans.items():
        mask=(out.dataset==d)&out.time_hours.notna()&((out.condition==c) if c!="all" else True); idx=out.index[mask]; x=out.loc[idx,["latent_1","latent_2","latent_3"]].to_numpy(float); xa=_apply_transform(x,tr); out.loc[idx,["aligned_sample_latent_1","aligned_sample_latent_2","aligned_sample_latent_3"]]=xa; out.loc[idx,"trajectory_branch"]=c; means=pd.DataFrame(xa,index=idx).groupby(out.loc[idx,"time_hours"]).transform("mean").to_numpy(); out.loc[idx,"within_time_residual_norm"]=np.linalg.norm(xa-means,axis=1)
    out.to_csv(STAGE23/"03_sample_level_aligned_states.csv",index=False); sm=out[out.time_hours.notna()&out.aligned_sample_latent_1.notna()].groupby(["dataset","trajectory_branch"]).agg(n_samples=("sample","size"),n_times=("time_hours","nunique"),mean_within_time_residual_norm=("within_time_residual_norm","mean"),median_within_time_residual_norm=("within_time_residual_norm","median"),p95_within_time_residual_norm=("within_time_residual_norm",lambda x:x.quantile(.95))).reset_index(); sm.to_csv(STAGE23/"01_within_time_residual_summary.csv",index=False); return out,sm
def stage2_4_out_of_sample_validation(st):
    branches=_branches(st); grid=np.linspace(0,1,25); rows=[]
    for d,b in branches:
        g=st[(st.dataset==d)&st.time_hours.notna()&((st.condition==b) if b!="all" else True)]; times=sorted(g.time_hours.unique())
        if len(times)<3:continue
        lo,hi=min(times),max(times)
        for hold in times:
            train_times=[x for x in times if x!=hold]
            if len(train_times)<2:continue
            m=g[g.time_hours.isin(train_times)].groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index(); tu=(m.index.to_numpy(float)-lo)/(hi-lo); target_train=m.to_numpy(float); u=(hold-lo)/(hi-lo)
            for rd,rb in branches:
                if (rd,rb)==(d,b):continue
                rg=st[(st.dataset==rd)&st.time_hours.notna()&((st.condition==rb) if rb!="all" else True)]; rm=rg.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index();
                if len(rm)<2:continue
                ru=(rm.index.to_numpy(float)-rm.index.min())/(rm.index.max()-rm.index.min()); refcurve=np.column_stack([np.interp(np.linspace(0,1,25),ru,rm[c]) for c in ["latent_1","latent_2","latent_3"]]); ref_at=np.array([np.interp(u,np.linspace(0,1,25),refcurve[:,j]) for j in range(3)]); target_curve=np.column_stack([np.interp(np.linspace(0,1,25),tu,target_train[:,j]) for j in range(3)]); tr=_fit_transform(target_curve,refcurve); xheld=_point_at_time(st,d,hold,None if b=="all" else b); pred=_apply_transform(xheld.reshape(1,-1),tr)[0]; err=float(np.linalg.norm(pred-ref_at)); naive=float(np.linalg.norm(xheld-ref_at)); rows.append({"target_trajectory":f"{d}:{b}","held_out_time_hours":hold,"reference_trajectory":f"{rd}:{rb}","n_training_times":len(train_times),"status":"tested","prediction_error":err,"naive_unaligned_error":naive,"error_reduction":naive-err,"relative_error_reduction":(naive-err)/(naive+1e-8)})
    o=pd.DataFrame(rows); o.to_csv(STAGE24/"01_leave_one_timepoint_out.csv",index=False); sm=o.groupby("target_trajectory").agg(n_tests=("prediction_error","size"),mean_prediction_error=("prediction_error","mean"),median_prediction_error=("prediction_error","median"),p95_prediction_error=("prediction_error",lambda x:x.quantile(.95)),mean_relative_error=("prediction_error","mean")).reset_index() if not o.empty else pd.DataFrame(); sm.to_csv(STAGE24/"02_oos_summary_by_trajectory.csv",index=False); o.to_csv(STAGE24/"03_oos_vs_naive_baseline.csv",index=False); pd.DataFrame({"metric":["tested_rows","target_trajectories_tested","mean_oos_error","median_oos_error","mean_error_reduction_vs_naive"],"value":[len(o),o.target_trajectory.nunique() if not o.empty else 0,o.prediction_error.mean() if not o.empty else np.nan,o.prediction_error.median() if not o.empty else np.nan,o.error_reduction.mean() if not o.empty else np.nan]}).to_csv(STAGE24/"04_stage24_decision_metrics.csv",index=False); return o
def _find_feature_matrix(ds):
    root=RESULTS/ds; files=[p for p in root.rglob("*.csv") if "PCA" not in p.name and "correlation" not in p.name.lower() and "variance" not in p.name.lower() and "metadata" not in p.name.lower() and p.stat().st_size<250_000_000]
    best=None
    for p in files:
        try:
            df=pd.read_csv(p,index_col=0,nrows=5)
            if df.shape[1]>=2 and df.shape[0]>=2 and any(str(c).startswith("GSM") for c in df.columns):best=p; break
        except Exception:pass
    return best
def _read_feature_matrix(path):
    df=pd.read_csv(path,index_col=0); df=df.apply(pd.to_numeric,errors="coerce"); df=df.dropna(how="all").dropna(axis=1,how="all"); df.index=df.index.astype(str); df.columns=df.columns.astype(str); return df
def stage2_5_biological_feature_harmonization(st):
    out=[]; mats={}
    for ds in DATASETS:
        p=_find_feature_matrix(ds)
        if p is None:out.append({"dataset":ds,"candidate_matrix":False,"path":"","n_features":0,"n_samples":0,"feature_id_type":"unknown","status":"no_compatible_csv_found"}); continue
        try:
            df=_read_feature_matrix(p); mats[ds]=df; out.append({"dataset":ds,"candidate_matrix":True,"path":str(p),"n_features":len(df),"n_samples":df.shape[1],"feature_id_type":"row_index","status":"candidate_feature_matrix"})
        except Exception as e:out.append({"dataset":ds,"candidate_matrix":True,"path":str(p),"n_features":0,"n_samples":0,"feature_id_type":"unknown","status":f"read_error:{type(e).__name__}"})
    audit=pd.DataFrame(out); audit.to_csv(STAGE25/"01_feature_harmonization_audit.csv",index=False)
    common=set.intersection(*(set(m.index) for m in mats.values())) if len(mats)>=2 else set(); overlap=[]
    for a in mats:
        for b in mats:
            if a<b:overlap.append({"dataset_a":a,"dataset_b":b,"common_feature_count":len(set(mats[a].index)&set(mats[b].index))})
    pd.DataFrame(overlap).to_csv(STAGE25/"02_pairwise_feature_overlap.csv",index=False)
    status="insufficient_compatible_feature_matrices"
    if len(common)>=50:
        parts=[]
        for ds,m in mats.items():
            q=m.loc[sorted(common)].copy(); q=q.groupby(level=0).mean(); q=(q.sub(q.mean(axis=1),axis=0)).div(q.std(axis=1).replace(0,np.nan),axis=0); q.columns=[f"{ds}:{c}" for c in q.columns]; parts.append(q)
        shared=pd.concat(parts,axis=1); shared.to_csv(STAGE25/"03_shared_feature_matrix.csv"); status="shared_feature_matrix_constructed"
    report=pd.DataFrame({"metric":["datasets_with_candidate_matrices","common_feature_count","status","note"],"value":[len(mats),len(common),status,"No cross-species ortholog mapping, scRNA pseudobulk or platform correction is applied yet; this stage is an audit/provisional feature space, not final biological validation."]}); report.to_csv(STAGE25/"04_stage25_decision.csv",index=False)
    print("\n"+"="*88+"\nSTAGE 2.5 — BIOLOGICAL FEATURE HARMONIZATION AUDIT\n"+"="*88)
    print(audit.to_string(index=False)); print(f"\ncommon exact feature IDs across compatible matrices = {len(common)}"); print(f"status = {status}"); print("NOTE: Stage 2.5 does not infer biological equivalence from PCA/time alignment. Human–mouse orthology, scRNA pseudobulk and platform-aware integration remain explicit next steps."); print("="*88)
    return audit

def stage3_trajectory_reconstruction(st):
    rows=[]
    for (d,b),g in st[st.time_hours.notna()&st.aligned_sample_latent_1.notna()].groupby(["dataset","trajectory_branch"]):
        if g.time_hours.nunique()<2:continue
        m=g.groupby("time_hours")[["aligned_sample_latent_1","aligned_sample_latent_2","aligned_sample_latent_3"]].mean().sort_index().reset_index(); m.insert(0,"dataset",d); m.insert(1,"trajectory_branch",b); m.columns=["dataset","trajectory_branch","time_hours","common_latent_1","common_latent_2","common_latent_3"]; rows.append(m)
    o=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(); o.to_csv(OUT/"stage3"/"01_reconstructed_trajectories.csv",index=False); return o
def derivative(x,t):
    if len(x)<2:return np.full(len(x),np.nan)
    return np.gradient(x,t) if len(x)>2 else np.repeat((x[1]-x[0])/(t[1]-t[0]),2)
def stage4_dynamics(tr):
    o=tr.copy()
    for _,idx in o.groupby(["dataset","trajectory_branch"]).groups.items():
        q=o.loc[idx].sort_values("time_hours"); t=q.time_hours.to_numpy(float)
        for j in (1,2,3):o.loc[q.index,f"dcommon_latent_{j}_dt"]=derivative(q[f"common_latent_{j}"].to_numpy(float),t)
        o.loc[q.index,"state_speed"]=np.sqrt(np.nansum(o.loc[q.index,[f"dcommon_latent_{j}_dt" for j in (1,2,3)]].to_numpy()**2,axis=1))
    o.to_csv(OUT/"stage4"/"01_dynamics.csv",index=False); return o
def stage5_critical_transitions(d):
    o=d.copy(); o["rolling_variance_latent1"]=np.nan; o["rolling_autocorrelation_latent1"]=np.nan
    for _,idx in o.groupby(["dataset","trajectory_branch"]).groups.items():
        q=o.loc[idx].sort_values("time_hours"); w=min(5,len(q))
        if w>=3:o.loc[q.index,"rolling_variance_latent1"]=q.common_latent_1.rolling(w,min_periods=3).var().to_numpy(); o.loc[q.index,"rolling_autocorrelation_latent1"]=q.common_latent_1.rolling(w,min_periods=3).apply(lambda z:z.autocorr(1) if z.std()>0 else np.nan).to_numpy()
    o["critical_transition_flag"]=False; o.to_csv(OUT/"stage5"/"01_critical_transition_indicators.csv",index=False); return o
def stage6_symbolic_equation_discovery(d):
    if d.empty:return pd.DataFrame()
    o=pd.DataFrame({"dataset":d.dataset,"trajectory_branch":d.trajectory_branch,"time_hours":d.time_hours,"x":d.common_latent_1,"target_dx_dt":d.dcommon_latent_1_dt}); o["x2"]=o.x**2; o["x3"]=o.x**3; o["abs_x"]=o.x.abs(); o["sin_x"]=np.sin(o.x); o["cos_x"]=np.cos(o.x); o.to_csv(OUT/"stage6"/"01_symbolic_regression_design.csv",index=False); return o
def stage7_heldout_validation(d):
    ds=sorted(d.dataset.unique()) if not d.empty else []; o=pd.DataFrame([{"held_out_dataset":x,"training_datasets":";".join(y for y in ds if y!=x),"validation_type":"leave-one-complete-dataset-out","status":"planned"} for x in ds]); o.to_csv(OUT/"stage7"/"01_heldout_validation_plan.csv",index=False); return o
def stage8_regulatory_integration(d):
    ds=sorted(d.dataset.unique()) if not d.empty else []; o=pd.DataFrame([{"expression_dataset":x,"regulatory_dataset":"GSE67520","integration":"time-aligned regulatory evidence","causal_claim":False,"status":"planned"} for x in ds]); o.to_csv(OUT/"stage8"/"01_regulatory_integration_plan.csv",index=False); return o
def stage9_predictive_ai(d):
    rows=[]
    for (ds,b),g in d.groupby(["dataset","trajectory_branch"]):
        g=g.sort_values("time_hours")
        for i in range(len(g)-1):
            a,z=g.iloc[i],g.iloc[i+1]; rows.append({"dataset":ds,"trajectory_branch":b,"time_t":a.time_hours,"time_next":z.time_hours,"dt_hours":z.time_hours-a.time_hours,"latent_1_t":a.common_latent_1,"latent_2_t":a.common_latent_2,"latent_3_t":a.common_latent_3,"latent_1_next":z.common_latent_1,"latent_2_next":z.common_latent_2,"latent_3_next":z.common_latent_3})
    o=pd.DataFrame(rows); o.to_csv(OUT/"stage9"/"01_next_state_prediction_table.csv",index=False); return o
def main():
    st,av=stage1_data_integration(); st=stage2_1_common_latent_state(st); stage2_2_validate_common_latent(st); st,_=stage2_3_within_time_residual_validation(st); o24=stage2_4_out_of_sample_validation(st); stage2_5_biological_feature_harmonization(st); tr=stage3_trajectory_reconstruction(st); dy=stage4_dynamics(tr); dy=stage5_critical_transitions(dy); s=stage6_symbolic_equation_discovery(dy); h=stage7_heldout_validation(dy); r=stage8_regulatory_integration(dy); p=stage9_predictive_ai(dy)
    print(f"\nDynamics v4.8 results written to: {OUT}"); print(f"Datasets with PCA: {int(av.PCA_file_found.sum())}/{len(av)}"); print(f"Stage 2.1 aligned trajectory datasets: {st[st.common_latent_status=='time_anchored_aligned'].dataset.nunique()}"); print(f"Stage 2.4 tested OOS rows: {len(o24)}"); print(f"Stage 2.5 candidate matrices: {len(pd.read_csv(STAGE25/'01_feature_harmonization_audit.csv').query('candidate_matrix == True'))}"); print(f"Trajectory timepoints: {len(tr)}"); print(f"Stage 6 symbolic rows: {len(s)}"); print(f"Stage 7 held-out datasets: {len(h)}"); print(f"Stage 8 regulatory rows: {len(r)}"); print(f"Stage 9 prediction pairs: {len(p)}")
if __name__=="__main__":main()
