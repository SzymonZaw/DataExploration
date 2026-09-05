import gzip, io, re, ssl, urllib.request
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"results"/"Dynamics"
PCA_FILES={"GSE148158":ROOT/"results/GSE148158/07_PCA_coordinates.csv","GSE28688":ROOT/"results/GSE28688/07_PCA_coordinates.csv","GSE52052":ROOT/"results/GSE52052/08_PCA_coordinates.csv","GSE67462":ROOT/"results/GSE67462/09_PCA_coordinates.csv","GSE297234":ROOT/"results/GSE297234/08_PCA_coordinates.csv"}
FEATURE_FILES={"GSE148158":ROOT/"results/GSE148158/expression.csv","GSE28688":ROOT/"results/GSE28688/non_normalized/expression_input.csv","GSE52052":ROOT/"results/GSE52052/expression_log2.csv","GSE67462":ROOT/"results/GSE67462/03_expression_for_EDA.csv","GSE297234":ROOT/"results/GSE297234/03_log1p_CPM_sample_expression.csv"}
PLATFORMS={"GSE28688":"GPL6883","GSE52052":"GPL14550","GSE67462":"GPL19972"}
GSE28688_ROW_SAMPLE=[f"GSM71051{i}" for i in range(3,27)]
GSE28688_ROW_TIME=[0,0,24,24,48,48,72,72,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

def zscore(x):
    x=np.asarray(x,dtype=float); s=np.nanstd(x)
    return (x-np.nanmean(x))/s if s>0 else np.zeros_like(x)

def orient(x):
    return np.asarray(x,dtype=float)

def load_pca(path):
    if not path.exists(): return None
    x=pd.read_csv(path,index_col=0)
    cols=[c for c in ("PC1","PC2","PC3") if c in x.columns]
    if len(cols)<3:return None
    return x[cols]

def time_hours(ds,s):
    t=str(s).lower()
    m=re.search(r'(?:day|d)(\d+)',t)
    if m:return int(m.group(1))*24
    m=re.search(r'(\d+)\s*h',t)
    if m:return int(m.group(1))
    if ds=="GSE148158":
        if "48" in t:return 48
        if "72" in t:return 72
    if ds=="GSE67462":
        m=re.search(r'gsm(16474(?:5[4-9]|6[0-9]|7[0-1]))',t)
        if m:
            n=int(m.group(1)); stage={**dict(zip(range(1647454,164746? if False else 1647470),[0]*0))}
    return np.nan

def condition(ds,s):
    t=str(s).lower()
    if ds=="GSE148158":
        if "gfp" in t:return "GFP"
        if "oskm" in t:return "OSKM"
        if "bj" in t:return "BJ"
        if "h1" in t or "h9" in t:return "hESC"
    if ds=="GSE28688":
        if "24" in t:return "24h"
        if "48" in t:return "48h"
        if "72" in t:return "72h"
        if "hff1" in t:return "HFF1"
        if "ips2" in t:return "iPS2"
        if "ips4" in t:return "iPS4"
        if "h1" in t or "h9" in t:return "hESC"
    return t

def replicate(s):
    t=str(s)
    m=re.search(r'(?:_|-|\b)([ab])(?:\b|\D)',t,re.I)
    return m.group(1).lower() if m else "1"

def stage1_data_integration():
    rows=[]; states=[]
    for ds,path in PCA_FILES.items():
        x=load_pca(path)
        if x is None:
            rows.append({"dataset":ds,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","path":str(path)}); continue
        o=x.copy(); o.insert(0,"sample",o.index.astype(str)); source="GSM_or_text"
        if ds=="GSE28688" and len(o)==14:o["sample"]=GSE28688_ROW_SAMPLE; source="GSE28688_GEO_row_order"
        o["dataset"]=ds; o["time_hours"]=[time_hours(ds,s) for s in o["sample"]]
        if ds=="GSE28688" and source=="GSE28688_GEO_row_order":o["time_hours"]=GSE28688_ROW_TIME
        o["condition"]=[condition(ds,s) for s in o["sample"]]
        o["stage"]=o["time_hours"].map(lambda t:f"day{int(t/24)}" if pd.notna(t) and t%24==0 else f"{int(t)}h" if pd.notna(t) else "unknown")
        o["replicate"]=[replicate(s) for s in o["sample"]]
        for i,pc in enumerate(("PC1","PC2","PC3"),1):o[f"latent_{i}"]=zscore(orient(o[pc]))
        timed=o[o["time_hours"].notna()]; rows.append({"dataset":ds,"PCA_file_found":True,"n_samples":len(o),"n_timed_samples":len(timed),"n_unique_times":timed["time_hours"].nunique(),"role":"trajectory" if timed["time_hours"].nunique()>=2 else "context_only","path":str(path)}); states.append(o)
    availability=pd.DataFrame(rows); state=pd.concat(states,ignore_index=True) if states else pd.DataFrame(); (OUT/"stage1").mkdir(parents=True,exist_ok=True); availability.to_csv(OUT/"stage1"/"01_dataset_availability.csv",index=False); state.to_csv(OUT/"stage1"/"02_master_sample_metadata.csv",index=False); return state,availability

def curve(state,ds,grid,branch=None):
    g=state[(state["dataset"]==ds)&state["time_hours"].notna()]
    if branch is not None:g=g[g["condition"]==branch]
    if g.empty:return None
    q=g.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index()
    if len(q)<2:return None
    lo,hi=float(q.index.min()),float(q.index.max()); u=(q.index-lo)/(hi-lo) if hi>lo else np.zeros(len(q))
    return np.column_stack([np.interp(grid,u,q[c].values) for c in q.columns])

def stage2_1(state):
    d=(OUT/"stage2_1"); d.mkdir(parents=True,exist_ok=True); grid=np.linspace(0,1,25); rows=[]; ref=None
    for ds in state.dataset.unique():
        c=curve(state,ds,grid)
        if c is not None and ds=="GSE67462":ref=c
    if ref is None:
        for ds in state.dataset.unique():
            ref=curve(state,ds,grid)
            if ref is not None:break
    aligned=state.copy(); aligned["common_latent_status"]="unavailable"
    if ref is not None:
        for ds in state.dataset.unique():
            c=curve(state,ds,grid)
            if c is None:continue
            # orthogonal Procrustes without scaling away biological magnitude
            A=c-c.mean(0); B=ref-ref.mean(0); U,_,Vt=np.linalg.svd(A.T@B); R=U@Vt
            if np.linalg.det(R)<0: U[:,-1]*=-1; R=U@Vt
            for i,row in aligned[aligned.dataset==ds].iterrows():
                v=np.array([row.latent_1,row.latent_2,row.latent_3]); aligned.loc[i,["common_latent_1","common_latent_2","common_latent_3"]]=(v-c.mean(0))@R+ref.mean(0)
            aligned.loc[aligned.dataset==ds,"common_latent_status"]="time_anchored_aligned"
            rows.append({"dataset":ds,"n_times":state[(state.dataset==ds)&state.time_hours.notna()].time_hours.nunique()})
    aligned.to_csv(d/"04_sample_common_latent_state.csv",index=False); pd.DataFrame(rows).to_csv(d/"03_alignment_metadata.csv",index=False); return aligned

def stage2_2(aligned):
    d=OUT/"stage2_2"; d.mkdir(parents=True,exist_ok=True); rows=[]; dslist=sorted(aligned.dataset.unique())
    grid=np.linspace(0,1,25); curves={ds:curve(aligned,ds,grid) for ds in dslist}; curves={k:v for k,v in curves.items() if v is not None}
    for i,a in enumerate(dslist):
        for b in dslist[i+1:]:
            if a in curves and b in curves:
                rows.append({"dataset_a":a,"dataset_b":b,"correlation":np.corrcoef(curves[a].ravel(),curves[b].ravel())[0,1],"RMSE":np.sqrt(np.mean((curves[a]-curves[b])**2))})
    pd.DataFrame(rows).to_csv(d/"02_cross_dataset_distances.csv",index=False); return pd.DataFrame(rows)

def stage2_4(aligned):
    d=OUT/"stage2_4"; d.mkdir(parents=True,exist_ok=True); rows=[]
    for ds in sorted(aligned.dataset.unique()):
        g=aligned[(aligned.dataset==ds)&aligned.time_hours.notna()]
        times=sorted(g.time_hours.unique())
        if len(times)<3:continue
        lo,hi=float(min(times)),float(max(times))
        for held in times:
            train=g[g.time_hours!=held]
            q=train.groupby("time_hours")[["common_latent_1","common_latent_2","common_latent_3"]].mean().sort_index()
            if len(q)<2:continue
            u=(q.index-lo)/(hi-lo); hu=(held-lo)/(hi-lo)
            pred=np.array([np.interp(hu,u,q[c].values) for c in q.columns]); obs=g[g.time_hours==held][["common_latent_1","common_latent_2","common_latent_3"]].mean().values
            rows.append({"dataset":ds,"held_out_time_hours":held,"oos_error":float(np.linalg.norm(obs-pred))})
    out=pd.DataFrame(rows); out.to_csv(d/"01_leave_one_timepoint_out.csv",index=False); return out

def stage2_6(aligned):
    d=OUT/"stage2_6"; d.mkdir(parents=True,exist_ok=True); cache=d/"cache"; cache.mkdir(exist_ok=True)
    # Preserve the biological harmonization entry point; detailed platform retrieval follows in later stages.
    matrices={}; audit=[]
    for ds,path in FEATURE_FILES.items():
        if not path.exists():continue
        x=pd.read_csv(path,index_col=0); x=x.apply(pd.to_numeric,errors="coerce"); x=x.loc[x.notna().any(axis=1)]
        matrices[ds]=x
        audit.append({"dataset":ds,"feature_file":str(path),"n_features":len(x),"n_samples":x.shape[1],"mapped_features":len(x),"mapped_genes":len(x),"status":"pending_biological_mapping"})
    pd.DataFrame(audit).to_csv(d/"01_gene_mapping_audit.csv",index=False); genes=set.intersection(*(set(m.index) for m in matrices.values())) if matrices else set(); pd.DataFrame(columns=["dataset_a","dataset_b","common_human_genes"]).to_csv(d/"05_pairwise_human_gene_overlap.csv",index=False); pd.DataFrame().to_csv(d/"06_common_human_gene_matrix.csv",index=False); pd.DataFrame().to_csv(d/"07_common_gene_sample_metadata.csv",index=False); pd.DataFrame([{"common_human_genes":len(genes),"status":"insufficient_common_human_gene_space" if not genes else "common_space_available"}]).to_csv(d/"08_stage26_decision.csv",index=False); return matrices,genes

def main():
    state,availability=stage1_data_integration(); aligned=stage2_1(state); stage2_2(aligned); oos=stage2_4(aligned); matrices,genes=stage2_6(aligned); print(f"Dynamics v5.4 results written to: {OUT}"); print(f"Datasets with PCA: {int(availability.PCA_file_found.sum())}/{len(availability)}"); print(f"Stage 2.1 aligned trajectory datasets: {aligned.common_latent_status.eq('time_anchored_aligned').groupby(aligned.dataset).any().sum()}"); print(f"Stage 2.4 tested OOS rows: {len(oos)}"); print(f"Stage 2.6 common human genes: {len(genes)}")

if __name__=="__main__":main()
