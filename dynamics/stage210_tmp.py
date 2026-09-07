"""Stage 2.10: leakage-free conserved biological transition modules."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"stage2_10"
OUT.mkdir(parents=True,exist_ok=True)
TARGET=["GSE67462","GSE28688","GSE297234"]
GRID=np.linspace(0,1,40)
N_BOOT=500
N_PERM=1000
PROGRAMS={
"P01_PLURIPOTENCY":["POU5F1","SOX2","NANOG","LIN28A","LIN28B","DPPA4","UTF1","ESRRB","KLF4","MYC"],
"P02_PROLIFERATION":["MKI67","PCNA","TOP2A","CCNB1","CCNB2","CCNE1","CDK1","CDC20","UBE2C","TYMS","MCM2","MCM3","MCM5","MCM6","MCM7"],
"P03_EMT_MESENCHYMAL":["VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","FN1","ITGA5","COL1A1","COL1A2","CDH2","CDH1","EPCAM","KRT8","KRT18","KRT19"],
"P04_STRESS_RESPONSE":["DDIT3","ATF4","HSPA1A","HSPA1B","HMOX1","XBP1","JUN","FOS","DUSP1","PPP1R15A","DNAJB1"],
"P05_GLYCOLYTIC_METABOLISM":["SLC2A1","HK2","PFKP","ALDOA","GAPDH","ENO1","PKM","LDHA","PGK1","TPI1","PDK1"],
"P06_FGFR_PI3K_MAPK":["FGFR1","FGFR2","FGFR3","FGFR4","FRS2","PLCG1","PIK3CA","PIK3CB","AKT1","AKT2","MAPK1","MAPK3","RAF1","SOS1"],
"P07_CHROMATIN_EPIGENETIC":["KMT2A","KMT2B","EZH2","DNMT1","DNMT3A","DNMT3B","TET1","TET2","HDAC1","HDAC2","SMARCA4","ARID1A","CHD4","SUZ12"],
"P08_ECM_ADHESION":["FN1","ITGA5","ITGB1","COL1A1","COL1A2","COL3A1","SPARC","VCAN","THBS1","LAMC1","LAMA4"]}

def log(x): print(f"Stage 2.10: {x}",flush=True)
def finite(X,fill=None):
    X=np.asarray(X,float).copy()
    if fill is None:
        fill=np.array([np.median(v[np.isfinite(v)]) if np.isfinite(v).any() else 0.0 for v in X.T])
    else: fill=np.where(np.isfinite(np.asarray(fill,float)),np.asarray(fill,float),0.0)
    bad=~np.isfinite(X)
    if bad.any():
        ii=np.where(bad); X[ii]=fill[ii[1]]
    return X,fill

def normtime(t):
    t=np.asarray(t,float); lo,hi=np.min(t),np.max(t)
    return np.zeros_like(t) if hi<=lo else (t-lo)/(hi-lo)

def resample(t,X,fill):
    X,_=finite(X,fill); tn=normtime(t); Y=np.empty((len(GRID),X.shape[1]))
    for j in range(X.shape[1]): Y[:,j]=np.interp(GRID,tn,X[:,j])
    return Y

def load():
    from dynamics.validation import _load_common_space
    matrix,meta=_load_common_space(); meta=meta[meta.dataset.astype(str).isin(TARGET)].copy(); out={}
    for ds in TARGET:
        g=meta[meta.dataset.astype(str).eq(ds)].copy(); g["time_hours"]=pd.to_numeric(g["time_hours"],errors="coerce"); g=g[g.time_hours.notna()&g.matrix_column.notna()]
        if g.time_hours.nunique()<3: continue
        X=matrix.loc[:,g.matrix_column.astype(str)].T.copy(); X.index=g.time_hours.to_numpy(float); X=X.groupby(level=0,sort=True).mean()
        out[ds]=(X.index.to_numpy(float),X.to_numpy(float),[str(x).upper() for x in matrix.index]); log(f"{ds}: {len(g)} timed samples -> {len(X)} unique timepoints")
    return out

def train_fill(train): return finite(np.vstack([x[1] for x in train.values()]))[1]
def activity(Y,genes,fill):
    Y,_=finite(Y,fill); frame=pd.DataFrame(Y.T,index=genes).rank(axis=0,pct=True); out={}
    for pid,members in PROGRAMS.items():
        ix=[g for g in members if g in frame.index]
        if len(ix)>=3: out[pid]=frame.loc[ix].mean(axis=0).to_numpy(float)
    return out

def transition(y):
    d=np.gradient(np.asarray(y,float),GRID); k=int(np.argmax(np.abs(d))); return float(GRID[k]),float(d[k])
def spearman(a,b):
    return float(pd.Series(a).corr(pd.Series(b),method="spearman")) if len(a)>=3 else np.nan

def fold(held,train_names,traj):
    train={d:traj[d] for d in train_names}; ht,hX,hgenes=traj[held]; fill=train_fill(train)
    train_times={}
    for ds,(t,X,genes) in train.items():
        A=activity(resample(t,X,fill),genes,fill); train_times[ds]={p:transition(a)[0] for p,a in A.items()}
    common=sorted(set.intersection(*(set(v) for v in train_times.values())))
    if len(common)<3:return None,None
    consensus={p:float(np.mean([train_times[d][p] for d in train_names])) for p in common}
    hA=activity(resample(ht,hX,fill),hgenes,fill); common=[p for p in common if p in hA]
    if len(common)<3:return None,None
    held={p:transition(hA[p])[0] for p in common}; a=np.array([consensus[p] for p in common]); b=np.array([held[p] for p in common])
    pairs=[]
    for i in range(len(common)):
        for j in range(i+1,len(common)):
            if a[i]!=a[j] and b[i]!=b[j]: pairs.append(np.sign(a[i]-a[j])==np.sign(b[i]-b[j]))
    row={"heldout_dataset":held,"n_training_datasets":len(train_names),"n_common_programs":len(common),"transition_rank_spearman":spearman(a,b),"transition_time_pearson":float(pd.Series(a).corr(pd.Series(b))) if np.std(a)>1e-12 and np.std(b)>1e-12 else np.nan,"pairwise_ordering_agreement":float(np.mean(pairs)) if pairs else np.nan}
    detail=pd.DataFrame({"heldout_dataset":held,"program_id":common,"train_consensus_transition_time":a,"heldout_transition_time":b,"train_consensus_rank":pd.Series(a).rank(method="average").to_numpy(),"heldout_rank":pd.Series(b).rank(method="average").to_numpy()})
    return row,detail

def bootstrap(held,train_names,traj):
    train={d:traj[d] for d in train_names}; ht,hX,hgenes=traj[held]; fill=train_fill(train); Y=resample(ht,hX,fill); pos={g:i for i,g in enumerate(hgenes)}; rng=np.random.default_rng(210000+sum(map(ord,held))); rows=[]
    for k in range(N_BOOT):
        times={}
        for pid,members in PROGRAMS.items():
            ix=[pos[g] for g in members if g in pos]
            if len(ix)<3: continue
            sample=rng.choice(ix,len(ix),replace=True); A=pd.DataFrame(Y[:,sample]).rank(axis=0,pct=True).mean(axis=1).to_numpy(); times[pid]=transition(A)[0]
        if len(times)>=3: rows.append({"bootstrap":k+1,"n_programs":len(times),"transition_time_sd":float(np.std(list(times.values())))})
    return pd.DataFrame(rows)

def permutation(traj):
    rng=np.random.default_rng(210031); rows=[]
    for k in range(N_PERM):
        vals=[]
        for held in TARGET:
            if held not in traj: continue
            train_names=[d for d in TARGET if d!=held and d in traj]; train={d:traj[d] for d in train_names}; fill=train_fill(train); _,consensus_detail=fold(held,train_names,traj)
            # Rebuild training consensus without touching heldout values.
            train_times={}
            for ds,(t,X,genes) in train.items(): train_times[ds]={p:transition(a)[0] for p,a in activity(resample(t,X,fill),genes,fill).items()}
            common=sorted(set.intersection(*(set(v) for v in train_times.values()))) if train_times else []
            if len(common)<3: continue
            ht,hX,hgenes=traj[held]; hA=activity(resample(ht,hX,fill),hgenes,fill); common=[p for p in common if p in hA]
            if len(common)<3: continue
            a=np.array([np.mean([train_times[d][p] for d in train_names]) for p in common]); b=np.array([transition(hA[p][rng.permutation(len(hA[p]))])[0] for p in common]); vals.append(spearman(a,b))
        rows.append({"permutation":k+1,"mean_transition_rank_spearman":float(np.nanmean(vals)) if vals else np.nan})
    return pd.DataFrame(rows)

def run():
    traj=load(); names=[d for d in TARGET if d in traj]; log("trajectory datasets: "+", ".join(names))
    if len(names)<3: raise RuntimeError("Stage 2.10 requires at least three trajectory datasets.")
    rows=[]; details=[]; boots=[]
    for held in names:
        train=[d for d in names if d!=held]; log(f"LODO held out {held}; training on {', '.join(train)}"); row,detail=fold(held,train,traj)
        if row is not None: rows.append(row); details.append(detail); bt=bootstrap(held,train,traj); bt.insert(0,"heldout_dataset",held); boots.append(bt)
    R=pd.DataFrame(rows)
    if len(R)!=len(names): raise RuntimeError("No valid Stage 2.10 LODO folds.")
    D=pd.concat(details,ignore_index=True); B=pd.concat(boots,ignore_index=True); P=permutation(traj); R.to_csv(OUT/"01_lodo_module_conservation.csv",index=False); D.to_csv(OUT/"02_module_transition_order.csv",index=False); B.to_csv(OUT/"03_bootstrap_transition_stability.csv",index=False); P.to_csv(OUT/"04_time_permutation_null.csv",index=False)
    obs=float(R.transition_rank_spearman.mean()); null=P.mean_transition_rank_spearman.dropna().to_numpy(); p=float((1+np.sum(np.abs(null)>=abs(obs)))/(len(null)+1)); order=float(R.pairwise_ordering_agreement.mean()); support=bool(obs>0.5 and order>0.65 and p<0.05)
    S=pd.DataFrame([{"n_trajectory_datasets":len(names),"n_valid_lodo_folds":len(R),"mean_transition_rank_spearman":obs,"mean_transition_time_pearson":float(R.transition_time_pearson.mean()),"mean_pairwise_ordering_agreement":order,"mean_common_programs":float(R.n_common_programs.mean()),"time_permutation_p":p,"conserved_transition_modules_supported":support,"stage3_readiness":False,"interpretation":"LODO validation of fixed biological programs as conserved transition modules; training-only imputation and consensus; bootstrap tests module stability; time permutation tests ordering; no ODE/state-space model."}]); S.to_csv(OUT/"05_stage210_summary.csv",index=False); log("overall:"); print(S.to_string(index=False)); return S

if __name__=="__main__": run()
