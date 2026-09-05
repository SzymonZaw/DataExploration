"""Dynamics v4.4: OSKM reprogramming dynamics pipeline.

Stages 1-5: integration, latent representation, trajectories, derivatives,
and critical-transition indicators.
Stage 2.1: time-anchored cross-study trajectory geometry.
Stage 2.2: alignment/concordance diagnostics.
Stage 2.3: sample-level residual validation of the common latent space.
Stages 6-9: symbolic discovery, held-out validation, regulatory integration,
and predictive-AI preparation.

Scientific boundary: heterogeneous datasets are not assumed to measure an
identical molecular state. Stage 2.3 explicitly tests whether the apparent
cross-study geometry survives when actual sample-level variation is retained.
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
DATASETS={
"GSE28688":RESULTS/"GSE28688"/"non_normalized"/"07_PCA_coordinates.csv",
"GSE148158":RESULTS/"GSE148158"/"07_PCA_coordinates.csv",
"GSE52052":RESULTS/"GSE52052"/"08_PCA_coordinates.csv",
"GSE67462":RESULTS/"GSE67462"/"09_PCA_coordinates.csv",
"GSE297234":RESULTS/"GSE297234"/"08_PCA_coordinates.csv"}
GSM_TIME={
"GSM4455240":48.,"GSM4455241":48.,"GSM4455242":72.,"GSM4455243":72.,"GSM4455244":48.,"GSM4455245":72.,
"GSM710515":24.,"GSM710516":24.,"GSM710517":48.,"GSM710518":48.,"GSM710519":72.,"GSM710520":72.,
"GSM1258008":264.,"GSM1258009":264.,"GSM1258010":264.,"GSM1258011":264.,"GSM1258012":264.,"GSM1258013":264.,
"GSM1647454":0.,"GSM1647455":0.,"GSM1647456":24.,"GSM1647457":24.,"GSM1647458":72.,"GSM1647459":72.,
"GSM1647460":120.,"GSM1647461":120.,"GSM1647462":168.,"GSM1647463":168.,"GSM1647464":264.,"GSM1647465":264.,
"GSM1647466":360.,"GSM1647467":360.,"GSM1647468":432.,"GSM1647469":432.,
"GSM8986586":0.,"GSM8986587":72.,"GSM8986588":168.,"GSM8986589":240.,"GSM8986590":0.,"GSM8986591":72.,"GSM8986592":168.,"GSM8986593":240.}
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
            rows.append({"dataset":ds,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","timing_source":"none","path":str(path)});continue
        o=x.copy();o.insert(0,"sample",o.index.astype(str)); source="GSM_or_text"
        if ds=="GSE28688" and len(o)==14:o["sample"]=GSE28688_ROW_SAMPLE;source="GSE28688_GEO_row_order"
        o["dataset"]=ds;o["time_hours"]=[time_hours(ds,s) for s in o["sample"]]
        if ds=="GSE28688" and source=="GSE28688_GEO_row_order":o["time_hours"]=GSE28688_ROW_TIME
        o["stage"]=o["time_hours"].map(lambda t:f"day{int(t/24)}" if pd.notna(t) and t%24==0 else (f"{int(t)}h" if pd.notna(t) else "unknown"))
        o["replicate"]=[replicate(s) for s in o["sample"]];o["timing_source"]=source
        for i,pc in enumerate(["PC1","PC2","PC3"],1):o[f"latent_{i}"]=zscore(orient(o[pc]))
        timed=o[o.time_hours.notna()];role="trajectory" if timed.time_hours.nunique()>=2 else "context_only"
        rows.append({"dataset":ds,"PCA_file_found":True,"n_samples":len(o),"n_timed_samples":len(timed),"n_unique_times":timed.time_hours.nunique(),"role":role,"timing_source":source,"path":str(path)});states.append(o)
    av=pd.DataFrame(rows);st=pd.concat(states,ignore_index=True) if states else pd.DataFrame();av.to_csv(STAGE_DIRS[1]/"01_dataset_availability.csv",index=False);st.to_csv(STAGE_DIRS[1]/"02_master_sample_metadata.csv",index=False);return st,av

def _curve(st,ds,grid):
    g=st[(st.dataset==ds)&st.time_hours.notna()];m=g.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index()
    if len(m)<2:return None
    t=m.index.to_numpy(float);u=(t-t.min())/(t.max()-t.min());return np.column_stack([np.interp(grid,u,m[c]) for c in ["latent_1","latent_2","latent_3"]])

def _rot(a,b):
    aa=a-a.mean(0);bb=b-b.mean(0);u,_,vt=np.linalg.svd(aa.T@bb);return u@vt

def stage2_1_common_latent_state(st):
    ds=[d for d in st.dataset.unique() if st[(st.dataset==d)&st.time_hours.notna()].time_hours.nunique()>=2];grid=np.linspace(0,1,25);cur={d:_curve(st,d,grid) for d in ds};cur={d:c for d,c in cur.items() if c is not None and np.isfinite(c).all()}
    if not cur:return st
    ref="GSE67462" if "GSE67462" in cur else sorted(cur)[0];aligned={ref:cur[ref]};r=cur[ref]
    for d,c in cur.items():
        if d==ref:continue
        a=(c-c.mean(0))@_rot(c,r);a*=np.divide(np.std(r,0),np.where(np.std(a,0)>0,np.std(a,0),1));a+=r.mean(0);aligned[d]=a
    rows=[]
    for d,a in aligned.items():
        g=st[(st.dataset==d)&st.time_hours.notna()];m=g.groupby("time_hours").size();lo,hi=g.time_hours.min(),g.time_hours.max()
        for t,n in m.items():
            u=(t-lo)/(hi-lo);v=np.array([np.interp(u,grid,a[:,j]) for j in range(3)]);rows.append({"dataset":d,"time_hours":t,"normalized_time":u,"common_latent_1":v[0],"common_latent_2":v[1],"common_latent_3":v[2],"n_replicates":n})
    pd.DataFrame(rows).to_csv(STAGE21/"01_common_latent_trajectory.csv",index=False)
    out=st.copy()
    for j in range(1,4):out[f"common_latent_{j}"]=np.nan
    out["common_latent_status"]="not_aligned"
    # IMPORTANT: these columns remain the historical time-anchored trajectory coordinates;
    # Stage 2.3 recomputes sample-level aligned coordinates from the actual latent samples.
    for d,a in aligned.items():
        idx=out.index[(out.dataset==d)&out.time_hours.notna()];lo,hi=out.loc[idx,"time_hours"].min(),out.loc[idx,"time_hours"].max();u=(out.loc[idx,"time_hours"].to_numpy()-lo)/(hi-lo)
        for j in range(3):out.loc[idx,f"common_latent_{j+1}"]=np.interp(u,grid,a[:,j])
        out.loc[idx,"common_latent_status"]="time_anchored_aligned"
    out.to_csv(STAGE21/"04_sample_common_latent_state.csv",index=False);pd.DataFrame({"reference_dataset":[ref],"trajectory_datasets":[";".join(sorted(aligned))],"grid_points":[25],"alignment":["orthogonal Procrustes on time-mean trajectories"]}).to_csv(STAGE21/"03_alignment_metadata.csv",index=False);return out

def stage2_2_validate_common_latent(st):
    g=st[st.common_latent_status=="time_anchored_aligned"];ds=sorted(g.dataset.unique());grid=np.linspace(0,1,25);cur={d:_curve(st,d,grid) for d in ds};cur={d:c for d,c in cur.items() if c is not None}
    pairs=[]
    for i,a in enumerate(ds):
        for b in ds[i+1:]:
            x,y=cur[a],cur[b];pairs.append({"dataset_a":a,"dataset_b":b,"trajectory_correlation":np.corrcoef(x.ravel(),y.ravel())[0,1],"aligned_rmse":np.sqrt(np.mean((x-y)**2)),"path_length_a":np.linalg.norm(np.diff(x,axis=0),axis=1).sum(),"path_length_b":np.linalg.norm(np.diff(y,axis=0),axis=1).sum()})
    p=pd.DataFrame(pairs);p["path_length_ratio"]=p.path_length_a/p.path_length_b;p.to_csv(STAGE22/"02_cross_dataset_distances.csv",index=False)
    q=[]
    for d in ds:
        r=np.mean([cur[x] for x in ds if x!=d],axis=0);x=cur[d];q.append({"dataset":d,"leave_one_dataset_out_reference_rmse":np.sqrt(np.mean((x-r)**2)),"leave_one_dataset_out_correlation":np.corrcoef(x.ravel(),r.ravel())[0,1]})
    q=pd.DataFrame(q);q.to_csv(STAGE22/"01_alignment_quality.csv",index=False)
    mean=np.mean(list(cur.values()),axis=0);cc=[]
    for k,u in enumerate(grid):cc.append({"normalized_time":u,"mean_common_latent_1":mean[k,0],"mean_common_latent_2":mean[k,1],"mean_common_latent_3":mean[k,2],"cross_dataset_dispersion":np.mean(np.linalg.norm(np.vstack([cur[d][k] for d in ds])-mean[k],axis=1))})
    cc=pd.DataFrame(cc);cc.to_csv(STAGE22/"03_trajectory_concordance.csv",index=False)
    sens=[]
    for ref in ds:
        for test in ds:
            if ref==test:continue
            c=_curve(st,test,grid);a=(c-c.mean(0))@_rot(c,cur[ref]);a*=np.divide(np.std(cur[ref],0),np.where(np.std(a,0)>0,np.std(a,0),1));a+=cur[ref].mean(0);sens.append({"reference_dataset":ref,"test_dataset":test,"reference_sensitivity_rmse":np.sqrt(np.mean((a-cur[ref])**2))})
    s=pd.DataFrame(sens);s.to_csv(STAGE22/"04_reference_sensitivity.csv",index=False)
    summary=pd.DataFrame({"metric":["n_trajectory_datasets","n_pairwise_comparisons","mean_pairwise_correlation","median_pairwise_rmse","mean_cross_dataset_dispersion","mean_reference_sensitivity_rmse"],"value":[len(ds),len(p),p.trajectory_correlation.mean(),p.aligned_rmse.median(),cc.cross_dataset_dispersion.mean(),s.reference_sensitivity_rmse.mean()]});summary.to_csv(STAGE22/"05_common_latent_validation.csv",index=False);return summary

def stage2_3_within_time_residual_validation(st):
    """Apply trajectory-derived Procrustes transforms to REAL samples.

    The transform is fitted only from time-mean trajectories. Individual PCA
    samples are then transformed without replacing them by their time mean.
    This preserves replicate/sample variation and allows a genuine residual test.
    """
    traj=[d for d in st.dataset.unique() if st[(st.dataset==d)&st.time_hours.notna()].time_hours.nunique()>=2];grid=np.linspace(0,1,25);raw={d:_curve(st,d,grid) for d in traj};raw={d:c for d,c in raw.items() if c is not None};ref="GSE67462" if "GSE67462" in raw else sorted(raw)[0];trans={ref:(np.eye(3),np.zeros(3),np.ones(3))};R=raw[ref]
    # Store rotation, source/target centering and axis scale; fit only on means.
    for d,c in raw.items():
        if d==ref:continue
        rot=_rot(c,R);a=(c-c.mean(0))@rot;scale=np.divide(np.std(R,0),np.where(np.std(a,0)>0,np.std(a,0),1));trans[d]=(rot,c.mean(0),scale)
    out=st.copy();
    for j in range(1,4):out[f"aligned_sample_latent_{j}"]=np.nan;out[f"aligned_time_mean_{j}"]=np.nan;out[f"within_time_residual_{j}"]=np.nan
    out["within_time_residual_norm"]=np.nan
    for d in traj:
        if d not in trans:continue
        rot,cm,scale=trans[d];idx=out.index[(out.dataset==d)&out.time_hours.notna()];x=out.loc[idx,["latent_1","latent_2","latent_3"]].to_numpy(float);xa=(x-cm)@rot;xa*=scale;xa+=R.mean(0)
        out.loc[idx,["aligned_sample_latent_1","aligned_sample_latent_2","aligned_sample_latent_3"]]=xa
        means=pd.DataFrame(xa,index=idx).groupby(out.loc[idx,"time_hours"]).transform("mean").to_numpy();out.loc[idx,["aligned_time_mean_1","aligned_time_mean_2","aligned_time_mean_3"]]=means
        res=xa-means
        for j in range(3):out.loc[idx,f"within_time_residual_{j+1}"]=res[:,j]
        out.loc[idx,"within_time_residual_norm"]=np.linalg.norm(res,axis=1)
    timed=out[out.time_hours.notna() & out.aligned_sample_latent_1.notna()].copy();summary=[];rep=[]
    for d,g in timed.groupby("dataset"):
        summary.append({"dataset":d,"n_samples":len(g),"n_times":g.time_hours.nunique(),"mean_within_time_residual_norm":g.within_time_residual_norm.mean(),"median_within_time_residual_norm":g.within_time_residual_norm.median(),"p95_within_time_residual_norm":g.within_time_residual_norm.quantile(.95),"within_time_sd_latent1":g.within_time_residual_1.std(ddof=1)})
        for t,h in g.groupby("time_hours"):
            if len(h)>=2:
                vals=h[["aligned_sample_latent_1","aligned_sample_latent_2","aligned_sample_latent_3"]].to_numpy();dist=np.linalg.norm(vals[0]-vals[1]) if len(vals)==2 else np.mean([np.linalg.norm(vals[i]-vals[j]) for i in range(len(vals)) for j in range(i+1,len(vals))]);rep.append({"dataset":d,"time_hours":t,"n_replicates":len(h),"replicate_pair_distance":dist})
    sm=pd.DataFrame(summary);rp=pd.DataFrame(rep);sm.to_csv(STAGE23/"01_within_time_residual_summary.csv",index=False);rp.to_csv(STAGE23/"02_replicate_distances.csv",index=False);out.to_csv(STAGE23/"03_sample_level_aligned_states.csv",index=False)
    # Compare raw vs aligned mean-trajectory concordance to detect alignment-created agreement.
    rows=[]
    for i,a in enumerate(traj):
        for b in traj[i+1:]:
            ca,cb=raw[a],raw[b];aa=(ca-ca.mean(0))@_rot(ca,cb);aa*=np.divide(np.std(cb,0),np.where(np.std(aa,0)>0,np.std(aa,0),1));aa+=cb.mean(0)
            rawcorr=np.corrcoef(ca.ravel(),cb.ravel())[0,1];alignedcorr=np.corrcoef(aa.ravel(),cb.ravel())[0,1]
            rows.append({"dataset_a":a,"dataset_b":b,"raw_trajectory_correlation":rawcorr,"post_procrustes_correlation":alignedcorr,"correlation_gain":alignedcorr-rawcorr,"raw_rmse":np.sqrt(np.mean((ca-cb)**2)),"post_procrustes_rmse":np.sqrt(np.mean((aa-cb)**2))})
    gain=pd.DataFrame(rows);gain.to_csv(STAGE23/"04_raw_vs_aligned_concordance.csv",index=False)
    # Between-study distance versus typical within-study replicate distance.
    between=[]
    for i,a in enumerate(traj):
        for b in traj[i+1:]:
            between.append({"dataset_a":a,"dataset_b":b,"mean_trajectory_distance":np.mean(np.linalg.norm(raw[a]-raw[b],axis=1))})
    bd=pd.DataFrame(between);wd=rp.groupby("dataset").replicate_pair_distance.median().rename("median_replicate_distance") if not rp.empty else pd.Series(dtype=float);ratio=[]
    for _,r in bd.iterrows():ratio.append({**r.to_dict(),"median_within_replicate_distance_a":wd.get(r.dataset_a,np.nan),"median_within_replicate_distance_b":wd.get(r.dataset_b,np.nan),"between_to_within_ratio":r.mean_trajectory_distance/np.nanmean([wd.get(r.dataset_a,np.nan),wd.get(r.dataset_b,np.nan)])})
    ratio=pd.DataFrame(ratio);ratio.to_csv(STAGE23/"05_between_vs_within_distance.csv",index=False)
    report=pd.DataFrame({"metric":["reference_dataset","trajectory_datasets","mean_within_time_residual","median_within_time_residual","mean_correlation_gain_from_alignment","mean_between_to_within_ratio"],"value":[ref,";".join(traj),sm.mean_within_time_residual_norm.mean(),sm.median_within_time_residual_norm.median(),gain.correlation_gain.mean(),ratio.between_to_within_ratio.mean()]});report.to_csv(STAGE23/"06_stage23_decision_metrics.csv",index=False)
    return out,report

def stage3_trajectory_reconstruction(st):
    rows=[]
    for d,g in st[st.time_hours.notna()].groupby("dataset"):
        if g.time_hours.nunique()<2:continue
        cols=[f"aligned_sample_latent_{i}" for i in (1,2,3)]
        if not all(c in g and g[c].notna().any() for c in cols):cols=[f"common_latent_{i}" for i in (1,2,3)]
        m=g.groupby("time_hours")[cols].mean().sort_index().reset_index();m.insert(0,"dataset",d);m.columns=["dataset","time_hours","common_latent_1","common_latent_2","common_latent_3"];rows.append(m)
    out=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame();out.to_csv(STAGE_DIRS[3]/"01_reconstructed_trajectories.csv",index=False);return out

def derivative(x,t):
    if len(x)<2:return np.full(len(x),np.nan)
    if np.any(np.diff(t)<=0):raise ValueError("Non-increasing time")
    return np.gradient(x,t) if len(x)>2 else np.repeat((x[1]-x[0])/(t[1]-t[0]),2)

def stage4_dynamics(tr):
    out=tr.copy()
    for d,idx in out.groupby("dataset").groups.items():
        q=out.loc[idx].sort_values("time_hours");t=q.time_hours.to_numpy(float)
        for j in (1,2,3):
            v=derivative(q[f"common_latent_{j}"].to_numpy(float),t);out.loc[q.index,f"dcommon_latent_{j}_dt"]=v
        out.loc[q.index,"state_speed"]=np.sqrt(np.nansum(out.loc[q.index,[f"dcommon_latent_{j}_dt" for j in (1,2,3)]].to_numpy()**2,axis=1))
    out.to_csv(STAGE_DIRS[4]/"01_dynamics.csv",index=False);return out

def stage5_critical_transitions(d):
    out=d.copy();out["rolling_variance_latent1"]=np.nan;out["rolling_autocorrelation_latent1"]=np.nan
    for _,idx in out.groupby("dataset").groups.items():
        q=out.loc[idx].sort_values("time_hours");s=q.common_latent_1;w=min(5,len(s));out.loc[q.index,"rolling_variance_latent1"]=s.rolling(w,min_periods=3).var().to_numpy();out.loc[q.index,"rolling_autocorrelation_latent1"]=s.rolling(w,min_periods=3).apply(lambda z:z.autocorr(1) if z.std()>0 else np.nan).to_numpy()
    out["critical_transition_flag"]=False;out.to_csv(STAGE_DIRS[5]/"01_critical_transition_indicators.csv",index=False);return out

def stage6_symbolic_equation_discovery(d):
    if d.empty:return pd.DataFrame()
    o=pd.DataFrame({"dataset":d.dataset,"time_hours":d.time_hours,"x":d.common_latent_1,"target_dx_dt":d.dcommon_latent_1_dt});o["x2"]=o.x**2;o["x3"]=o.x**3;o["abs_x"]=o.x.abs();o["sin_x"]=np.sin(o.x);o["cos_x"]=np.cos(o.x);o.to_csv(STAGE_DIRS[6]/"01_symbolic_regression_design.csv",index=False);return o

def stage7_heldout_validation(d):
    ds=sorted(d.dataset.unique()) if not d.empty else [];o=pd.DataFrame([{"held_out_dataset":x,"training_datasets":";".join(y for y in ds if y!=x),"validation_type":"leave-one-complete-dataset-out","status":"planned"} for x in ds]);o.to_csv(STAGE_DIRS[7]/"01_heldout_validation_plan.csv",index=False);return o

def stage8_regulatory_integration(d):
    ds=sorted(d.dataset.unique()) if not d.empty else [];o=pd.DataFrame([{"expression_dataset":x,"regulatory_dataset":"GSE67520","integration":"time-aligned regulatory evidence","causal_claim":False,"status":"planned"} for x in ds]);o.to_csv(STAGE_DIRS[8]/"01_regulatory_integration_plan.csv",index=False);return o

def stage9_predictive_ai(d):
    rows=[]
    for ds,g in d.groupby("dataset"):
        g=g.sort_values("time_hours")
        for i in range(len(g)-1):
            a,b=g.iloc[i],g.iloc[i+1];rows.append({"dataset":ds,"time_t":a.time_hours,"time_next":b.time_hours,"dt_hours":b.time_hours-a.time_hours,"latent_1_t":a.common_latent_1,"latent_2_t":a.common_latent_2,"latent_3_t":a.common_latent_3,"latent_1_next":b.common_latent_1,"latent_2_next":b.common_latent_2,"latent_3_next":b.common_latent_3})
    o=pd.DataFrame(rows);o.to_csv(STAGE_DIRS[9]/"01_next_state_prediction_table.csv",index=False);return o

def print_stage23(out,report):
    print("\n"+"="*88+"\nSTAGE 2.3 — WITHIN-TIME RESIDUAL VALIDATION\n"+"="*88)
    for title,file in [("RESIDUAL SUMMARY","01_within_time_residual_summary.csv"),("REPLICATE DISTANCES","02_replicate_distances.csv"),("RAW VS ALIGNED","04_raw_vs_aligned_concordance.csv"),("BETWEEN VS WITHIN","05_between_vs_within_distance.csv"),("DECISION METRICS","06_stage23_decision_metrics.csv")]:
        print(f"\n--- {title} ---");p=STAGE23/file
        if p.exists():
            with pd.option_context("display.max_rows",None,"display.max_columns",None,"display.width",220):print(pd.read_csv(p).to_string(index=False))
    print("\n--- INTERPRETATION ---")
    if report is not None and not report.empty:
        v=dict(zip(report.metric,report.value));print(f"mean_within_time_residual = {float(v['mean_within_time_residual']):.6f}");print(f"median_within_time_residual = {float(v['median_within_time_residual']):.6f}");print(f"mean_correlation_gain_from_alignment = {float(v['mean_correlation_gain_from_alignment']):.6f}");print(f"mean_between_to_within_ratio = {float(v['mean_between_to_within_ratio']):.6f}")
    print("NOTE: Stage 2.3 is the first check that retains actual replicate/sample variation. A low residual alone is not sufficient; alignment gain and between/within distance must also be considered.");print("="*88+"\n")

def main():
    st,av=stage1_data_integration();st=stage2_1_common_latent_state(st);v22=stage2_2_validate_common_latent(st);st,v23=stage2_3_within_time_residual_validation(st);print_stage23(st,v23);tr=stage3_trajectory_reconstruction(st);dy=stage4_dynamics(tr);dy=stage5_critical_transitions(dy);s=stage6_symbolic_equation_discovery(dy);h=stage7_heldout_validation(dy);r=stage8_regulatory_integration(dy);p=stage9_predictive_ai(dy)
    (OUT/"REPORT.txt").write_text("Dynamics v4.4\n\nStage 2.3 now validates sample-level residuals after applying transforms fitted only on time-mean trajectories. Stages 6-9 remain preparation layers; symbolic equations and AI models are not fitted.\n",encoding="utf-8")
    print(f"Dynamics v4.4 results written to: {OUT}");print(f"Datasets with PCA: {int(av.PCA_file_found.sum())}/{len(DATASETS)}");print(f"Stage 2.1 aligned trajectory datasets: {st[st.common_latent_status=='time_anchored_aligned'].dataset.nunique()}");print(f"Stage 2.2 validation rows: {len(v22)}");print(f"Stage 2.3 residual datasets: {st[st.aligned_sample_latent_1.notna()].dataset.nunique()}");print(f"Trajectory timepoints: {len(tr)}");print(f"Stage 6 symbolic rows: {len(s)}");print(f"Stage 7 held-out datasets: {len(h)}");print(f"Stage 8 regulatory rows: {len(r)}");print(f"Stage 9 prediction pairs: {len(p)}")

if __name__=="__main__":main()
