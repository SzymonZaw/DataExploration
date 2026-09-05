"""Dynamics v4.0: staged OSKM reprogramming dynamics pipeline.

Stages 1-5: integration, latent representation, trajectories, derivatives,
and critical-transition indicators.
Stages 6-9: symbolic equation discovery, held-out validation, regulatory
integration, and predictive AI preparation.

Scientific boundary: Stage 2 remains study-normalized PCA rather than a proven
shared biological latent space. Stages 6-9 prepare analyses; they do not claim
a fitted biological law or causality.
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
    "GSM1647460":120.,"GSM1647461":120.,"GSM1647462":168.,"GSM1647463":168.,"GSM1647464":264.,"GSM1647465":264.,
    "GSM1647466":360.,"GSM1647467":360.,"GSM1647468":432.,"GSM1647469":432.,
    "GSM8986586":0.,"GSM8986587":72.,"GSM8986588":168.,"GSM8986589":240.,"GSM8986590":0.,"GSM8986591":72.,"GSM8986592":168.,"GSM8986593":240.,
}
GSE28688_ROW_SAMPLE = [f"GSM{x}" for x in range(710513,710527)]
GSE28688_ROW_TIME = [0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

TIME_PATTERNS = {
    "GSE28688":[(r"24\s*h",24.),(r"48\s*h",48.),(r"72\s*h",72.)],
    "GSE148158":[(r"48",48.),(r"72",72.)],
    "GSE52052":[(r"day\s*11",264.)],
    "GSE67462":[(r"day\s*0\b",0.),(r"day\s*1\b",24.),(r"day\s*3\b",72.),(r"day\s*5\b",120.),(r"day\s*7\b",168.),(r"day\s*11\b",264.),(r"day\s*15\b",360.),(r"day\s*18\b",432.)],
    "GSE297234":[(r"d0\b|day\s*0\b",0.),(r"d3\b|day\s*3\b",72.),(r"d7\b|day\s*7\b",168.),(r"d10\b|day\s*10\b",240.)],
}


def load_pca(path):
    if not path.exists(): return None
    df=pd.read_csv(path,index_col=0)
    pcs=[c for c in ("PC1","PC2","PC3") if c in df.columns]
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
    s=str(sample)
    m=re.search(r"(?:-|_|\s)([ab])$",s,re.I)
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
    """Stage 1: integrate PCA outputs and harmonise metadata."""
    availability=[]; parts=[]
    for dataset,path in DATASETS.items():
        coords=load_pca(path)
        if coords is None:
            availability.append({"dataset":dataset,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","timing_source":"none","path":str(path)}); continue
        pcs=[c for c in ("PC1","PC2","PC3") if c in coords.columns]
        out=coords[pcs].copy(); out.insert(0,"sample",out.index.astype(str)); source="GSM_or_text"
        if dataset=="GSE28688" and len(out)==14: out["sample"]=GSE28688_ROW_SAMPLE; source="GSE28688_GEO_row_order"
        out["dataset"]=dataset; out["time_hours"]=[time_hours(dataset,s) for s in out.sample]
        if dataset=="GSE28688" and source=="GSE28688_GEO_row_order": out["time_hours"]=GSE28688_ROW_TIME
        out["stage"]=[stage_label(s,t) for s,t in zip(out.sample,out.time_hours)]; out["replicate"]=[replicate(s) for s in out.sample]; out["timing_source"]=source
        for pc in pcs: out[f"{pc}_z"]=zscore(orient(out[pc]))
        timed=out[out.time_hours.notna()]; role="trajectory" if timed.time_hours.nunique()>=2 else "context_only"
        availability.append({"dataset":dataset,"PCA_file_found":True,"n_samples":len(out),"n_timed_samples":len(timed),"n_unique_times":timed.time_hours.nunique(),"role":role,"timing_source":source,"path":str(path)}); parts.append(out)
    av=pd.DataFrame(availability); states=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame()
    av.to_csv(STAGE_DIRS[1]/"01_dataset_availability.csv",index=False); states.to_csv(STAGE_DIRS[1]/"02_master_sample_metadata.csv",index=False)
    return states,av


def stage2_latent_state(states):
    """Stage 2: study-normalized low-dimensional state representation."""
    if states.empty: return states.copy()
    out=states.copy()
    for i,pc in enumerate(("PC1_z","PC2_z","PC3_z"),1): out[f"latent_{i}"]=out[pc] if pc in out else np.nan
    out["latent_space_type"]="study_normalized_PCA"
    out[["dataset","sample","time_hours","stage","replicate","latent_1","latent_2","latent_3","latent_space_type"]].to_csv(STAGE_DIRS[2]/"01_latent_state_coordinates.csv",index=False)
    return out


def stage3_trajectory_reconstruction(states):
    """Stage 3: reconstruct trajectories after replicate averaging."""
    if states.empty: return pd.DataFrame()
    rows=[]
    for dataset,g in states[states.time_hours.notna()].groupby("dataset"):
        if g.time_hours.nunique()<2: continue
        cols=[c for c in ("latent_1","latent_2","latent_3") if c in g]
        m=g.groupby("time_hours",as_index=False)[cols].mean().sort_values("time_hours"); m.insert(0,"dataset",dataset)
        m["n_replicates"]=g.groupby("time_hours").size().reindex(m.time_hours).to_numpy(); rows.append(m)
    traj=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame()
    traj.to_csv(STAGE_DIRS[3]/"01_reconstructed_trajectories.csv",index=False)
    return traj


def stage4_dynamics(traj):
    """Stage 4: first/second derivatives, speed and acceleration."""
    if traj.empty: return traj.copy()
    out=traj.copy()
    for dataset,idx in out.groupby("dataset").groups.items():
        sub=out.loc[idx].sort_values("time_hours"); t=sub.time_hours.to_numpy(float)
        for axis in ("latent_1","latent_2","latent_3"):
            x=sub[axis].to_numpy(float); v=derivative(x,t) if len(x)>=2 else np.full(len(x),np.nan); a=derivative(v,t) if len(x)>=3 else np.full(len(x),np.nan)
            out.loc[sub.index,f"d{axis}_dt"]=v; out.loc[sub.index,f"d2{axis}_dt2"]=a
    vcols=[f"dlatent_{i}_dt" for i in (1,2,3)]; acols=[f"d2latent_{i}_dt2" for i in (1,2,3)]
    out["state_speed"]=np.sqrt(np.nansum(out[vcols].to_numpy(float)**2,axis=1)); out["state_acceleration"]=np.sqrt(np.nansum(out[acols].to_numpy(float)**2,axis=1))
    out.to_csv(STAGE_DIRS[4]/"01_dynamics.csv",index=False); return out


def stage5_critical_transitions(dynamics):
    """Stage 5: exploratory variance/autocorrelation stability indicators."""
    if dynamics.empty: return dynamics.copy()
    out=dynamics.copy(); out["rolling_variance_latent1"]=np.nan; out["rolling_autocorrelation_latent1"]=np.nan
    for dataset,idx in out.groupby("dataset").groups.items():
        sub=out.loc[idx].sort_values("time_hours"); s=sub.latent_1
        if len(s)>=3:
            w=min(5,len(s)); var=s.rolling(w,min_periods=3).var(); ac=s.rolling(w,min_periods=3).apply(lambda q:q.autocorr(1) if q.std()>0 else np.nan)
            out.loc[sub.index,"rolling_variance_latent1"]=var.to_numpy(); out.loc[sub.index,"rolling_autocorrelation_latent1"]=ac.to_numpy()
    out["critical_transition_flag"]=False; out.to_csv(STAGE_DIRS[5]/"01_critical_transition_indicators.csv",index=False); return out


def stage6_symbolic_equation_discovery(dynamics):
    """Stage 6: build symbolic-regression design; do not fit an equation yet."""
    if dynamics.empty: return pd.DataFrame()
    rows=[]
    for dataset,g in dynamics.groupby("dataset"):
        x=g.latent_1.to_numpy(float); target=g.dlatent_1_dt.to_numpy(float)
        rows.append(pd.DataFrame({"dataset":dataset,"time_hours":g.time_hours.to_numpy(float),"x":x,"target_dx_dt":target,"x2":x*x,"x3":x*x*x,"abs_x":np.abs(x),"sin_x":np.sin(x),"cos_x":np.cos(x),"log1p_abs_x":np.log1p(np.abs(x))}))
    table=pd.concat(rows,ignore_index=True) if rows else pd.DataFrame(); table.to_csv(STAGE_DIRS[6]/"01_symbolic_regression_design.csv",index=False)
    pd.DataFrame({"candidate_form":["dx/dt = f(x)","dx/dt = f(x,y)","dx/dt = f(x,y,z)"],"method":["symbolic regression","symbolic regression","symbolic regression"],"validation":["complete-dataset holdout"]*3,"status":["planned"]*3}).to_csv(STAGE_DIRS[6]/"02_symbolic_plan.csv",index=False)
    return table


def stage7_heldout_validation(dynamics):
    """Stage 7: leave-one-complete-dataset-out validation plan."""
    datasets=sorted(dynamics.dataset.unique()) if not dynamics.empty else []; rows=[]
    for test in datasets: rows.append({"held_out_dataset":test,"training_datasets":";".join(d for d in datasets if d!=test),"validation_type":"leave-one-complete-dataset-out","status":"planned"})
    plan=pd.DataFrame(rows); plan.to_csv(STAGE_DIRS[7]/"01_heldout_validation_plan.csv",index=False); return plan


def stage8_regulatory_integration(dynamics):
    """Stage 8: prepare time-aligned integration with GSE67520 ChIP-seq."""
    if dynamics.empty: return pd.DataFrame()
    rows=[{"expression_dataset":d,"regulatory_dataset":"GSE67520","integration":"time-aligned regulatory evidence","causal_claim":False,"status":"planned"} for d in sorted(dynamics.dataset.unique())]
    plan=pd.DataFrame(rows); plan.to_csv(STAGE_DIRS[8]/"01_regulatory_integration_plan.csv",index=False); return plan


def stage9_predictive_ai(dynamics):
    """Stage 9: prepare next-state prediction pairs; do not train AI yet."""
    if dynamics.empty: return pd.DataFrame()
    rows=[]
    for dataset,g in dynamics.groupby("dataset"):
        g=g.sort_values("time_hours")
        for i in range(len(g)-1):
            a,b=g.iloc[i],g.iloc[i+1]; rows.append({"dataset":dataset,"time_t":a.time_hours,"time_next":b.time_hours,"dt_hours":b.time_hours-a.time_hours,"latent_1_t":a.latent_1,"latent_2_t":a.latent_2,"latent_3_t":a.latent_3,"latent_1_next":b.latent_1,"latent_2_next":b.latent_2,"latent_3_next":b.latent_3})
    table=pd.DataFrame(rows); table.to_csv(STAGE_DIRS[9]/"01_next_state_prediction_table.csv",index=False)
    pd.DataFrame({"model_target":["z(t+dt) from z(t), dt"],"validation_unit":["complete_dataset"],"status":["preparation_only"]}).to_csv(STAGE_DIRS[9]/"02_ai_prediction_plan.csv",index=False); return table


def main():
    states,av=stage1_data_integration(); states=stage2_latent_state(states); traj=stage3_trajectory_reconstruction(states); dyn=stage4_dynamics(traj); dyn=stage5_critical_transitions(dyn)
    sym=stage6_symbolic_equation_discovery(dyn); hold=stage7_heldout_validation(dyn); reg=stage8_regulatory_integration(dyn); pred=stage9_predictive_ai(dyn)
    report=["Dynamics v4.0 — stages 1-9","","Stage 1: data integration and metadata harmonisation.","Stage 2: study-normalized latent state representation.","Stage 3: time-aware trajectory reconstruction.","Stage 4: derivative-based dynamics.","Stage 5: critical-transition/stability indicators.","Stage 6: symbolic equation-discovery design; no equation fitted.","Stage 7: complete-dataset held-out validation plan.","Stage 8: regulatory integration plan using GSE67520.","Stage 9: next-state predictive AI design; no model fitted.","","Scientific boundary: Stage 2 is not yet a proven cross-study biological latent space. Stages 6-9 remain preparation layers until the common representation and leakage-safe validation are established."]
    (OUT/"REPORT.txt").write_text("\n".join(report),encoding="utf-8")
    print(f"Dynamics v4.0 results written to: {OUT}"); print(f"Datasets with PCA: {int(av.PCA_file_found.sum())}/{len(DATASETS)}"); print(f"Trajectory timepoints: {len(traj)}"); print(f"Stage 6 symbolic rows: {len(sym)}"); print(f"Stage 7 held-out datasets: {len(hold)}"); print(f"Stage 8 regulatory rows: {len(reg)}"); print(f"Stage 9 prediction pairs: {len(pred)}")


if __name__=="__main__": main()
