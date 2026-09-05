"""Dynamics v4.3: staged OSKM reprogramming dynamics pipeline.

Stages 1-5: integration, latent representation, trajectories, derivatives,
and critical-transition indicators.
Stage 2.1: time-anchored cross-study common latent representation.
Stage 2.2: validation of common latent geometry and leakage-safe alignment.
Stages 6-9: symbolic equation discovery, held-out validation, regulatory
integration, and predictive AI preparation.

Scientific boundary: Stage 2.1 is a geometry/time-alignment layer, not a
claim that heterogeneous platforms measure an identical molecular state.
Stage 2.2 quantifies alignment/concordance; it does not establish biological
universality. GSE67520 remains regulatory evidence and GSE52052 is context-only.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "Dynamics"
STAGE_DIRS = {i: OUT / f"stage{i}" for i in range(1, 10)}
for d in STAGE_DIRS.values(): d.mkdir(parents=True, exist_ok=True)
STAGE21 = OUT / "stage2_1"; STAGE21.mkdir(parents=True, exist_ok=True)
STAGE22 = OUT / "stage2_2"; STAGE22.mkdir(parents=True, exist_ok=True)

DATASETS = {
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "07_PCA_coordinates.csv",
    "GSE148158": RESULTS / "GSE148158" / "07_PCA_coordinates.csv",
    "GSE52052": RESULTS / "GSE52052" / "08_PCA_coordinates.csv",
    "GSE67462": RESULTS / "GSE67462" / "09_PCA_coordinates.csv",
    "GSE297234": RESULTS / "GSE297234" / "08_PCA_coordinates.csv",
}
GSM_TIME = {
    "GSM4455240":48.,"GSM4455241":48.,"GSM4455242":72.,"GSM4455243":72.,"GSM4455244":48.,"GSM4455245":72.,
    "GSM710515":24.,"GSM710516":24.,"GSM710517":48.,"GSM710518":48.,"GSM710519":72.,"GSM710520":72.,
    "GSM1258008":264.,"GSM1258009":264.,"GSM1258010":264.,"GSM1258011":264.,"GSM1258012":264.,"GSM1258013":264.,
    "GSM1647454":0.,"GSM1647455":0.,"GSM1647456":24.,"GSM1647457":24.,"GSM1647458":72.,"GSM1647459":72.,
    "GSM1647460":120.,"GSM1647461":168.,"GSM1647462":168.,"GSM1647463":168.,"GSM1647464":264.,"GSM1647465":264.,
    "GSM1647466":360.,"GSM1647467":360.,"GSM1647468":432.,"GSM1647469":432.,
    "GSM8986586":0.,"GSM8986587":72.,"GSM8986588":168.,"GSM8986589":240.,"GSM8986590":0.,"GSM8986591":72.,"GSM8986592":168.,"GSM8986593":240.,
}
# Correct the two remaining GSE67462 replicate times explicitly.
GSM_TIME.update({"GSM1647460":120.,"GSM1647461":120.})
GSE28688_ROW_SAMPLE = [f"GSM{x}" for x in range(710513,710527)]
GSE28688_ROW_TIME = [0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
TIME_PATTERNS = {
    "GSE28688":[(r"24\s*h",24.),(r"48\s*h",48.),(r"72\s*h",72.)],
    "GSE148158":[(r"48",48.),(r"72",72.)], "GSE52052":[(r"day\s*11",264.)],
    "GSE67462":[(r"day\s*0\b",0.),(r"day\s*1\b",24.),(r"day\s*3\b",72.),(r"day\s*5\b",120.),(r"day\s*7\b",168.),(r"day\s*11\b",264.),(r"day\s*15\b",360.),(r"day\s*18\b",432.)],
    "GSE297234":[(r"d0\b|day\s*0\b",0.),(r"d3\b|day\s*3\b",72.),(r"d7\b|day\s*7\b",168.),(r"d10\b|day\s*10\b",240.)],
}

def load_pca(path):
    if not path.exists(): return None
    df=pd.read_csv(path,index_col=0); pcs=[c for c in ("PC1","PC2","PC3") if c in df.columns]
    if not pcs: return None
    out=df[pcs].apply(pd.to_numeric,errors="coerce"); out.index=out.index.astype(str); return out

def time_hours(dataset,sample):
    s=str(sample).strip().strip('"')
    if s in GSM_TIME: return GSM_TIME[s]
    text=s.lower().replace("_"," ").replace("-"," ")
    for pat,val in TIME_PATTERNS.get(dataset,[]):
        if re.search(pat,text): return val
    return np.nan

def stage_label(sample,t):
    text=str(sample).lower()
    if re.search(r"ipsc|\bips\b",text): return "iPSC"
    if re.search(r"hesc",text): return "hESC"
    if pd.notna(t): return f"day{int(t/24)}" if t%24==0 else f"{int(t)}h"
    return "unknown"

def replicate(sample):
    s=str(sample); m=re.search(r"(?:-|_|\s)([ab])$",s,re.I)
    if m: return m.group(1).lower()
    m=re.search(r"(?:rep|replicate)[_\s-]*(\d+)",s,re.I)
    if m: return m.group(1)
    m=re.fullmatch(r"GSM(\d+)",s)
    if m and 1647454<=int(m.group(1))<=1647469: return "1" if int(m.group(1))%2==0 else "2"
    return "unknown"

def zscore(x):
    x=pd.Series(x,dtype=float); sd=x.std(ddof=0)
    return (x-x.mean())/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index)

def orient(x):
    x=pd.Series(x,dtype=float).copy(); v=x.dropna()
    if v.empty: return x
    i=v.abs().idxmax(); return -x if x.loc[i]<0 else x

def derivative(x,t):
    x=np.asarray(x,float); t=np.asarray(t,float)
    if len(x)<2: return np.full(len(x),np.nan)
    if np.any(np.diff(t)<=0): raise ValueError("Duplicate/non-increasing time reached derivative().")
    if len(x)==2:
        v=(x[1]-x[0])/(t[1]-t[0]); return np.array([v,v])
    return np.gradient(x,t)

def stage1_data_integration():
    availability=[]; parts=[]
    for dataset,path in DATASETS.items():
        coords=load_pca(path)
        if coords is None:
            availability.append({"dataset":dataset,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","timing_source":"none","path":str(path)}); continue
        pcs=[c for c in ("PC1","PC2","PC3") if c in coords.columns]; out=coords[pcs].copy(); out.insert(0,"sample",out.index.astype(str)); source="GSM_or_text"
        if dataset=="GSE28688" and len(out)==14: out["sample"]=GSE28688_ROW_SAMPLE; source="GSE28688_GEO_row_order"
        out["dataset"]=dataset; out["time_hours"]=[time_hours(dataset,s) for s in out["sample"]]
        if dataset=="GSE28688" and source=="GSE28688_GEO_row_order": out["time_hours"]=GSE28688_ROW_TIME
        out["stage"]=[stage_label(s,t) for s,t in zip(out["sample"],out["time_hours"])]
        out["replicate"]=[replicate(s) for s in out["sample"]]; out["timing_source"]=source
        for pc in pcs: out[f"{pc}_z"]=zscore(orient(out[pc]))
        timed=out[out["time_hours"].notna()]; role="trajectory" if timed["time_hours"].nunique()>=2 else "context_only"
        availability.append({"dataset":dataset,"PCA_file_found":True,"n_samples":len(out),"n_timed_samples":len(timed),"n_unique_times":timed["time_hours"].nunique(),"role":role,"timing_source":source,"path":str(path)}); parts.append(out)
    av=pd.DataFrame(availability); states=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame()
    av.to_csv(STAGE_DIRS[1]/"01_dataset_availability.csv",index=False); states.to_csv(STAGE_DIRS[1]/"02_master_sample_metadata.csv",index=False); return states,av

def stage2_latent_state(states):
    if states.empty: return states.copy()
    out=states.copy()
    for i,pc in enumerate(("PC1_z","PC2_z","PC3_z"),1): out[f"latent_{i}"]=out[pc] if pc in out else np.nan
    out["latent_space_type"]="study_normalized_PCA"; return out

def _trajectory_grid(states,dataset,grid):
    g=states[(states["dataset"]==dataset)&states["time_hours"].notna()].copy()
    if g["time_hours"].nunique()<2: return None
    m=g.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index(); t=m.index.to_numpy(float); u=(t-t.min())/(t.max()-t.min())
    return np.column_stack([np.interp(grid,u,m[c].to_numpy(float)) for c in ("latent_1","latent_2","latent_3")])

def _procrustes_rotation(source,target):
    a=source-source.mean(axis=0); b=target-target.mean(axis=0); ua,_,vta=np.linalg.svd(a.T@b); return ua@vta

def stage2_1_common_latent_state(states):
    if states.empty: return states.copy()
    trajectory_datasets=[d for d in states["dataset"].unique() if states.loc[(states["dataset"]==d)&states["time_hours"].notna(),"time_hours"].nunique()>=2]; grid=np.linspace(0,1,25)
    curves={d:_trajectory_grid(states,d,grid) for d in trajectory_datasets}; curves={d:c for d,c in curves.items() if c is not None and np.isfinite(c).all()}
    if not curves: return states.copy()
    reference="GSE67462" if "GSE67462" in curves else sorted(curves)[0]; ref=curves[reference]; aligned={reference:ref}
    for d,c in curves.items():
        if d==reference: continue
        aligned[d]=(c-c.mean(axis=0))@_procrustes_rotation(c,ref); s=np.std(aligned[d],axis=0,ddof=0); r=np.std(ref,axis=0,ddof=0); aligned[d]=aligned[d]*np.divide(r,np.where(s>0,s,1.0)); aligned[d]+=ref.mean(axis=0)
    rows=[]
    for dataset in trajectory_datasets:
        if dataset not in aligned: continue
        g=states[(states["dataset"]==dataset)&states["time_hours"].notna()].copy(); m=g.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index(); t=m.index.to_numpy(float); u=(t-t.min())/(t.max()-t.min()); a=aligned[dataset]
        for tt,uu in zip(t,u): rows.append({"dataset":dataset,"time_hours":tt,"normalized_time":uu,"common_latent_1":np.interp(uu,grid,a[:,0]),"common_latent_2":np.interp(uu,grid,a[:,1]),"common_latent_3":np.interp(uu,grid,a[:,2]),"n_replicates":int((g["time_hours"]==tt).sum())})
    common=pd.DataFrame(rows); common.to_csv(STAGE21/"01_common_latent_trajectory.csv",index=False)
    out=states.copy(); out["common_latent_1"]=np.nan; out["common_latent_2"]=np.nan; out["common_latent_3"]=np.nan; out["common_latent_status"]="not_aligned"
    for dataset in trajectory_datasets:
        if dataset not in aligned: continue
        idx=out.index[(out["dataset"]==dataset)&out["time_hours"].notna()]; g=out.loc[idx]; all_t=g["time_hours"].to_numpy(float); lo,hi=all_t.min(),all_t.max(); u=(g["time_hours"].to_numpy(float)-lo)/(hi-lo); a=aligned[dataset]
        for j in range(3): out.loc[idx,f"common_latent_{j+1}"]=np.interp(u,grid,a[:,j])
        out.loc[idx,"common_latent_status"]="time_anchored_aligned"
    context=out[out["common_latent_status"]!="time_anchored_aligned"][["dataset","sample","time_hours","stage"]]; context.to_csv(STAGE21/"02_context_only_observations.csv",index=False)
    pd.DataFrame({"reference_dataset":[reference],"trajectory_datasets":[";".join(sorted(curves))],"grid_points":[len(grid)],"alignment":["orthogonal Procrustes on normalized-time trajectories"],"fit_uses_context_only_data":[False],"interpretation":["trajectory geometry; not raw-expression equivalence"]}).to_csv(STAGE21/"03_alignment_metadata.csv",index=False); out.to_csv(STAGE21/"04_sample_common_latent_state.csv",index=False); return out

def _pairwise_metrics(a,b):
    xa=a[["common_latent_1","common_latent_2","common_latent_3"]].to_numpy(float); xb=b[["common_latent_1","common_latent_2","common_latent_3"]].to_numpy(float)
    corr=np.corrcoef(xa.ravel(),xb.ravel())[0,1] if np.std(xa)>0 and np.std(xb)>0 else np.nan
    rmse=float(np.sqrt(np.mean((xa-xb)**2))); path_a=float(np.sum(np.linalg.norm(np.diff(xa,axis=0),axis=1))); path_b=float(np.sum(np.linalg.norm(np.diff(xb,axis=0),axis=1)))
    return corr,rmse,path_a,path_b

def stage2_2_validate_common_latent(states):
    """Stage 2.2: quantify concordance and test whether alignment is reference-dependent."""
    if states.empty: return pd.DataFrame()
    traj=states[(states["common_latent_status"]=="time_anchored_aligned")&states["time_hours"].notna()].copy(); datasets=sorted(traj["dataset"].unique()); grid=np.linspace(0,1,25); curves={}
    for d in datasets:
        g=traj[traj["dataset"]==d].groupby("time_hours")[["common_latent_1","common_latent_2","common_latent_3"]].mean().sort_index(); t=g.index.to_numpy(float); u=(t-t.min())/(t.max()-t.min()); curves[d]=np.column_stack([np.interp(grid,u,g[c].to_numpy(float)) for c in ("common_latent_1","common_latent_2","common_latent_3")])
    pair_rows=[]
    for i,a in enumerate(datasets):
        for b in datasets[i+1:]:
            aa=pd.DataFrame(curves[a],columns=["common_latent_1","common_latent_2","common_latent_3"]); bb=pd.DataFrame(curves[b],columns=aa.columns); corr,rmse,pa,pb=_pairwise_metrics(aa,bb); pair_rows.append({"dataset_a":a,"dataset_b":b,"trajectory_correlation":corr,"aligned_rmse":rmse,"path_length_a":pa,"path_length_b":pb,"path_length_ratio":pa/pb if pb>0 else np.nan})
    pairwise=pd.DataFrame(pair_rows); pairwise.to_csv(STAGE22/"02_cross_dataset_distances.csv",index=False)
    quality=[]
    for d in datasets:
        c=curves[d]; ref=np.mean([curves[x] for x in datasets if x!=d],axis=0) if len(datasets)>1 else c; rmse=float(np.sqrt(np.mean((c-ref)**2))); corr=float(np.corrcoef(c.ravel(),ref.ravel())[0,1]) if np.std(c)>0 and np.std(ref)>0 else np.nan; quality.append({"dataset":d,"leave_one_dataset_out_reference_rmse":rmse,"leave_one_dataset_out_correlation":corr})
    q=pd.DataFrame(quality); q.to_csv(STAGE22/"01_alignment_quality.csv",index=False)
    mean_curve=np.mean(list(curves.values()),axis=0); rows=[]
    for k,u in enumerate(grid):
        vals=np.vstack([curves[d][k] for d in datasets]); rows.append({"normalized_time":u,"mean_common_latent_1":mean_curve[k,0],"mean_common_latent_2":mean_curve[k,1],"mean_common_latent_3":mean_curve[k,2],"cross_dataset_dispersion":float(np.mean(np.linalg.norm(vals-mean_curve[k],axis=1)))})
    concord=pd.DataFrame(rows); concord.to_csv(STAGE22/"03_trajectory_concordance.csv",index=False)
    sens=[]; raw_states=states.copy()
    for ref in datasets:
        ref_curve=curves[ref]
        for d in datasets:
            if d==ref: continue
            c=_trajectory_grid(raw_states,d,grid); aligned=(c-c.mean(axis=0))@_procrustes_rotation(c,ref_curve); s=np.std(aligned,axis=0); r=np.std(ref_curve,axis=0); aligned=aligned*np.divide(r,np.where(s>0,s,1.0)); aligned+=ref_curve.mean(axis=0); rmse=float(np.sqrt(np.mean((aligned-ref_curve)**2))); sens.append({"reference_dataset":ref,"test_dataset":d,"reference_sensitivity_rmse":rmse})
    sensitivity=pd.DataFrame(sens); sensitivity.to_csv(STAGE22/"04_reference_sensitivity.csv",index=False)
    summary=pd.DataFrame({"metric":["n_trajectory_datasets","n_pairwise_comparisons","mean_pairwise_correlation","median_pairwise_rmse","mean_cross_dataset_dispersion","mean_reference_sensitivity_rmse"],"value":[len(datasets),len(pairwise),pairwise["trajectory_correlation"].mean() if len(pairwise) else np.nan,pairwise["aligned_rmse"].median() if len(pairwise) else np.nan,concord["cross_dataset_dispersion"].mean() if len(concord) else np.nan,sensitivity["reference_sensitivity_rmse"].mean() if len(sensitivity) else np.nan]}); summary.to_csv(STAGE22/"05_common_latent_validation.csv",index=False)
    report=["Dynamics v4.3 — Stage 2.2 validation","","Purpose: quantify whether the common latent geometry is stable across independent datasets.","","Metrics: pairwise trajectory correlation, aligned RMSE, path-length ratio, leave-one-dataset-out reference comparison, normalized-time dispersion, and reference sensitivity.","","Interpretation boundary: high concordance supports reproducibility of trajectory geometry, but does not prove molecular equivalence, causality, or universality across species/platforms."]
    (STAGE22/"REPORT.txt").write_text("\n".join(report),encoding="utf-8"); return summary

def print_stage22_validation(summary):
    """Print the Stage 2.2 evidence needed for the go/no-go decision."""
    files = {
        "ALIGNMENT QUALITY": STAGE22 / "01_alignment_quality.csv",
        "PAIRWISE DISTANCES": STAGE22 / "02_cross_dataset_distances.csv",
        "TRAJECTORY CONCORDANCE": STAGE22 / "03_trajectory_concordance.csv",
        "REFERENCE SENSITIVITY": STAGE22 / "04_reference_sensitivity.csv",
        "SUMMARY": STAGE22 / "05_common_latent_validation.csv",
    }
    print("\n" + "=" * 88)
    print("STAGE 2.2 VALIDATION — FULL CONSOLE INSPECTION")
    print("=" * 88)
    for title, path in files.items():
        print(f"\n--- {title} ---")
        if not path.exists():
            print(f"MISSING: {path}")
            continue
        try:
            df = pd.read_csv(path)
            if df.empty:
                print("EMPTY")
            else:
                with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 220, "display.max_colwidth", 80):
                    print(df.to_string(index=False))
        except Exception as exc:
            print(f"ERROR READING {path}: {exc}")
    if summary is not None and not summary.empty:
        print("\n--- AUTOMATIC GO/NO-GO FLAGS ---")
        vals=dict(zip(summary["metric"],summary["value"]))
        corr=vals.get("mean_pairwise_correlation",np.nan); rmse=vals.get("median_pairwise_rmse",np.nan); disp=vals.get("mean_cross_dataset_dispersion",np.nan); sens=vals.get("mean_reference_sensitivity_rmse",np.nan)
        print(f"mean_pairwise_correlation = {corr:.6f}")
        print(f"median_pairwise_rmse = {rmse:.6f}")
        print(f"mean_cross_dataset_dispersion = {disp:.6f}")
        print(f"mean_reference_sensitivity_rmse = {sens:.6f}")
        print("NOTE: these are diagnostics, not statistical proof. High concordance can be produced by the time-anchored alignment itself.")
        print("NEXT DECISION: inspect within-time residuals before treating the common latent space as biologically informative.")
    print("=" * 88 + "\n")

def stage3_trajectory_reconstruction(states):
    if states.empty: return pd.DataFrame()
    rows=[]; timed=states[states["time_hours"].notna()]
    for dataset,g in timed.groupby("dataset"):
        if g["time_hours"].nunique()<2: continue
        cols=[c for c in ("common_latent_1","common_latent_2","common_latent_3") if c in g and g[c].notna().any()]
        if len(cols)<3: continue
        m=g.groupby("time_hours")[cols].mean().sort_index().reset_index(); m.insert(0,"dataset",dataset); m["n_replicates"]=g.groupby("time_hours").size().reindex(m["time_hours"]).to_numpy(); rows.append(m)
    traj=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(); traj.to_csv(STAGE_DIRS[3]/"01_reconstructed_trajectories.csv",index=False); return traj

def stage4_dynamics(traj):
    if traj.empty: return traj.copy()
    out=traj.copy()
    for dataset,idx in out.groupby("dataset").groups.items():
        sub=out.loc[idx].sort_values("time_hours"); t=sub["time_hours"].to_numpy(float)
        for axis in ("common_latent_1","common_latent_2","common_latent_3"):
            x=sub[axis].to_numpy(float); v=derivative(x,t) if len(x)>=2 else np.full(len(x),np.nan); a=derivative(v,t) if len(x)>=3 else np.full(len(x),np.nan); out.loc[sub.index,f"d{axis}_dt"]=v; out.loc[sub.index,f"d2{axis}_dt2"]=a
    vcols=[f"dcommon_latent_{i}_dt" for i in (1,2,3)]; acols=[f"d2common_latent_{i}_dt2" for i in (1,2,3)]; out["state_speed"]=np.sqrt(np.nansum(out[vcols].to_numpy(float)**2,axis=1)); out["state_acceleration"]=np.sqrt(np.nansum(out[acols].to_numpy(float)**2,axis=1)); out.to_csv(STAGE_DIRS[4]/"01_dynamics.csv",index=False); return out

def stage5_critical_transitions(dynamics):
    if dynamics.empty: return dynamics.copy()
    out=dynamics.copy(); out["rolling_variance_latent1"]=np.nan; out["rolling_autocorrelation_latent1"]=np.nan
    for dataset,idx in out.groupby("dataset").groups.items():
        sub=out.loc[idx].sort_values("time_hours"); s=sub["common_latent_1"]
        if len(s)>=3:
            w=min(5,len(s)); out.loc[sub.index,"rolling_variance_latent1"]=s.rolling(w,min_periods=3).var().to_numpy(); out.loc[sub.index,"rolling_autocorrelation_latent1"]=s.rolling(w,min_periods=3).apply(lambda q:q.autocorr(1) if q.std()>0 else np.nan).to_numpy()
    out["critical_transition_flag"]=False; out.to_csv(STAGE_DIRS[5]/"01_critical_transition_indicators.csv",index=False); return out

def stage6_symbolic_equation_discovery(dynamics):
    if dynamics.empty: return pd.DataFrame()
    rows=[]
    for dataset,g in dynamics.groupby("dataset"):
        x=g["common_latent_1"].to_numpy(float); target=g["dcommon_latent_1_dt"].to_numpy(float); rows.append(pd.DataFrame({"dataset":dataset,"time_hours":g["time_hours"].to_numpy(float),"x":x,"target_dx_dt":target,"x2":x*x,"x3":x*x*x,"abs_x":np.abs(x),"sin_x":np.sin(x),"cos_x":np.cos(x),"log1p_abs_x":np.log1p(np.abs(x))}))
    table=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(); table.to_csv(STAGE_DIRS[6]/"01_symbolic_regression_design.csv",index=False); pd.DataFrame({"candidate_form":["dx/dt = f(x)","dx/dt = f(x,y)","dx/dt = f(x,y,z)"],"method":["symbolic regression"]*3,"validation":["complete-dataset holdout"]*3,"status":["planned"]*3}).to_csv(STAGE_DIRS[6]/"02_symbolic_plan.csv",index=False); return table

def stage7_heldout_validation(dynamics):
    datasets=sorted(dynamics["dataset"].unique()) if not dynamics.empty else []; rows=[{"held_out_dataset":test,"training_datasets":";".join(d for d in datasets if d!=test),"validation_type":"leave-one-complete-dataset-out","status":"planned"} for test in datasets]; plan=pd.DataFrame(rows); plan.to_csv(STAGE_DIRS[7]/"01_heldout_validation_plan.csv",index=False); return plan

def stage8_regulatory_integration(dynamics):
    if dynamics.empty: return pd.DataFrame()
    plan=pd.DataFrame([{"expression_dataset":d,"regulatory_dataset":"GSE67520","integration":"time-aligned regulatory evidence","causal_claim":False,"status":"planned"} for d in sorted(dynamics["dataset"].unique())]); plan.to_csv(STAGE_DIRS[8]/"01_regulatory_integration_plan.csv",index=False); return plan

def stage9_predictive_ai(dynamics):
    if dynamics.empty: return pd.DataFrame()
    rows=[]
    for dataset,g in dynamics.groupby("dataset"):
        g=g.sort_values("time_hours")
        for i in range(len(g)-1):
            a,b=g.iloc[i],g.iloc[i+1]; rows.append({"dataset":dataset,"time_t":a["time_hours"],"time_next":b["time_hours"],"dt_hours":b["time_hours"]-a["time_hours"],"latent_1_t":a["common_latent_1"],"latent_2_t":a["common_latent_2"],"latent_3_t":a["common_latent_3"],"latent_1_next":b["common_latent_1"],"latent_2_next":b["common_latent_2"],"latent_3_next":b["common_latent_3"]})
    table=pd.DataFrame(rows); table.to_csv(STAGE_DIRS[9]/"01_next_state_prediction_table.csv",index=False); pd.DataFrame({"model_target":["z(t+dt) from z(t), dt"],"validation_unit":["complete_dataset"],"status":["preparation_only"]}).to_csv(STAGE_DIRS[9]/"02_ai_prediction_plan.csv",index=False); return table

def main():
    states,av=stage1_data_integration(); states=stage2_latent_state(states); states=stage2_1_common_latent_state(states); validation=stage2_2_validate_common_latent(states); print_stage22_validation(validation); traj=stage3_trajectory_reconstruction(states); dyn=stage4_dynamics(traj); dyn=stage5_critical_transitions(dyn); sym=stage6_symbolic_equation_discovery(dyn); hold=stage7_heldout_validation(dyn); reg=stage8_regulatory_integration(dyn); pred=stage9_predictive_ai(dyn)
    report=["Dynamics v4.3 — stages 1-9","","Stage 1: data integration and metadata harmonisation.","Stage 2: study-normalized latent state representation.","Stage 2.1: time-anchored common trajectory geometry.","Stage 2.2: common latent-space validation with full console inspection.","Stage 3: time-aware trajectory reconstruction.","Stage 4: derivative-based dynamics.","Stage 5: critical-transition/stability indicators.","Stage 6: symbolic equation-discovery design; no equation fitted.","Stage 7: complete-dataset held-out validation plan.","Stage 8: regulatory integration plan using GSE67520.","Stage 9: next-state predictive AI design; no model fitted.","","Scientific boundary: Stage 2.1/2.2 test trajectory geometry, not raw-expression equivalence or biological universality. Stages 6-9 remain preparation layers until validation is satisfactory."]
    (OUT/"REPORT.txt").write_text("\n".join(report),encoding="utf-8")
    print(f"Dynamics v4.3 results written to: {OUT}"); print(f"Datasets with PCA: {int(av['PCA_file_found'].sum())}/{len(DATASETS)}"); print(f"Stage 2.1 aligned trajectory datasets: {states['common_latent_status'].eq('time_anchored_aligned').groupby(states['dataset']).any().sum() if not states.empty else 0}"); print(f"Stage 2.2 validation rows: {len(validation)}"); print(f"Trajectory timepoints: {len(traj)}"); print(f"Stage 6 symbolic rows: {len(sym)}"); print(f"Stage 7 held-out datasets: {len(hold)}"); print(f"Stage 8 regulatory rows: {len(reg)}"); print(f"Stage 9 prediction pairs: {len(pred)}")

if __name__=="__main__": main()
