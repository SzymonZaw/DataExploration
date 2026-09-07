"""Stage 2.9.32: shared-vs-dataset decomposition with biological invariance tests.

This stage asks whether the cross-dataset temporal component contains a
biologically interpretable invariant signal rather than merely a time clock.
All feature selection, scaling and PCA fitting are training-only under LODO.
The held-out dataset is evaluated for (1) cross-dataset agreement, (2)
biological-program association, (3) dataset-specific residual magnitude, and
(4) one-step/final-point predictive value against persistence. No ODE or
state-space model is fit.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_32"
OUT.mkdir(parents=True, exist_ok=True)
TARGET = ["GSE67462", "GSE28688", "GSE297234"]
GRID = np.linspace(0.0, 1.0, 40)
MIN_GENES, MAX_GENES = 200, 3000
N_PERM = 1000

PROGRAMS = {
    "P01_PLURIPOTENCY": ["POU5F1","SOX2","NANOG","LIN28A","LIN28B","DPPA4","UTF1","ESRRB","KLF4","MYC"],
    "P02_PROLIFERATION": ["MKI67","PCNA","TOP2A","CCNB1","CCNB2","CCNE1","CDK1","CDC20","UBE2C","TYMS","MCM2","MCM3","MCM5","MCM6","MCM7"],
    "P03_EMT_MESENCHYMAL": ["VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","FN1","ITGA5","COL1A1","COL1A2","CDH2","CDH1","EPCAM","KRT8","KRT18","KRT19"],
    "P04_STRESS_RESPONSE": ["DDIT3","ATF4","HSPA1A","HSPA1B","HMOX1","XBP1","JUN","FOS","DUSP1","PPP1R15A","DNAJB1"],
    "P05_GLYCOLYTIC_METABOLISM": ["SLC2A1","HK2","PFKP","ALDOA","GAPDH","ENO1","PKM","LDHA","PGK1","TPI1","PDK1"],
    "P06_FGFR_PI3K_MAPK": ["FGFR1","FGFR2","FGFR3","FGFR4","FRS2","PLCG1","PIK3CA","PIK3CB","AKT1","AKT2","MAPK1","MAPK3","RAF1","SOS1"],
    "P07_CHROMATIN_EPIGENETIC": ["KMT2A","KMT2B","EZH2","DNMT1","DNMT3A","DNMT3B","TET1","TET2","HDAC1","HDAC2","SMARCA4","ARID1A","CHD4","SUZ12"],
    "P08_ECM_ADHESION": ["FN1","ITGA5","ITGB1","COL1A1","COL1A2","COL3A1","SPARC","VCAN","THBS1","LAMC1","LAMA4"],
}

def log(x): print(f"Stage 2.9.32: {x}", flush=True)

def finite(X, fill=None):
    X = np.asarray(X, float).copy()
    med = np.nanmedian(X, axis=0) if fill is None else np.asarray(fill, float)
    med = np.where(np.isfinite(med), med, 0.0)
    bad = ~np.isfinite(X)
    if bad.any():
        ii = np.where(bad)
        X[ii] = med[ii[1]]
    return X, med

def normalize_time(t):
    t = np.asarray(t, float)
    lo, hi = np.min(t), np.max(t)
    return np.zeros_like(t) if hi <= lo else (t-lo)/(hi-lo)

def resample(t, X, fill=None):
    tn = normalize_time(t)
    X, med = finite(X, fill)
    Y = np.empty((len(GRID), X.shape[1]))
    for j in range(X.shape[1]): Y[:,j] = np.interp(GRID, tn, X[:,j])
    return Y, med

def baseline_scale(Y, baseline=None, amp=None):
    Y, _ = finite(Y)
    baseline = Y[0].copy() if baseline is None else np.asarray(baseline,float)
    Z = Y-baseline
    if amp is None:
        amp=np.sqrt(np.mean(Z*Z,axis=0)); amp=np.where(amp>1e-8,amp,1.0)
    return Z/amp, baseline, amp

def corr(a,b,method="pearson"):
    a,b=np.asarray(a,float),np.asarray(b,float); ok=np.isfinite(a)&np.isfinite(b)
    if ok.sum()<3 or np.std(a[ok])<1e-12 or np.std(b[ok])<1e-12:return np.nan
    return float(pd.Series(a[ok]).corr(pd.Series(b[ok]),method=method))

def load_trajectories():
    from dynamics.validation import _load_common_space
    matrix, meta = _load_common_space()
    meta=meta[meta["dataset"].astype(str).isin(TARGET)].copy(); out={}
    for ds in TARGET:
        g=meta[meta["dataset"].astype(str).eq(ds)].copy()
        g["time_hours"]=pd.to_numeric(g["time_hours"],errors="coerce")
        g=g[g.time_hours.notna()&g.matrix_column.notna()]
        if g.time_hours.nunique()<3: continue
        cols=g.matrix_column.astype(str).tolist(); X=matrix.loc[:,cols].T.copy()
        X.index=g.time_hours.to_numpy(float)
        X=X.groupby(level=0,sort=True).mean()
        out[ds]=(X.index.to_numpy(float),X.to_numpy(float),list(matrix.index))
        log(f"{ds}: {len(g)} timed samples -> {len(X)} unique timepoints")
    return out

def select_shared_genes(train, all_genes):
    states=[]
    for t,X,_ in train.values(): states.append(baseline_scale(resample(t,X)[0])[0])
    A=np.stack(states); mean=A.mean(axis=0)
    shared=np.var(mean,axis=0); hetero=np.mean((A-mean[None,:,:])**2,axis=(0,1))
    score=shared/(shared+hetero+1e-12); order=np.argsort(score)[::-1]
    chosen=order[score[order]>=0.55]
    if len(chosen)<MIN_GENES: chosen=order[:min(MIN_GENES,len(order))]
    return np.asarray(all_genes)[chosen[:MAX_GENES]], score[chosen[:MAX_GENES]]

def fit_shared(train, all_genes):
    genes,scores=select_shared_genes(train,all_genes); pos={str(g):i for i,g in enumerate(all_genes)}
    idx=[pos[str(g)] for g in genes]; states=[]
    for t,X,_ in train.values(): states.append(baseline_scale(resample(t,X[:,idx])[0])[0])
    stack=np.vstack(states); pca=PCA(n_components=1,random_state=2932).fit(stack)
    pcs=[pca.transform(s).ravel() for s in states]
    # Orient the shared axis so its training mean has positive endpoint change.
    shared=np.mean(np.stack(pcs),axis=0)
    if shared[-1]<shared[0]: shared=-shared; pcs=[-p for p in pcs]; pca.components_=-pca.components_
    return genes,idx,scores,pca,shared,pcs

def program_activity(Y, genes):
    frame=pd.DataFrame(Y,index=[str(g).upper() for g in genes]).rank(axis=0,pct=True)
    out={}
    for pid, members in PROGRAMS.items():
        ix=[m for m in members if m in frame.index]
        if len(ix)>=3: out[pid]=frame.loc[ix].mean(axis=0).to_numpy(float)
    return out

def prepare_fold(heldout, train_names, traj):
    train={d:traj[d] for d in train_names}; t,X,all_genes=traj[heldout]
    if len(t)<4:return None
    genes,idx,sel_scores,pca,shared,_=fit_shared(train,all_genes)
    prefix_raw, med=resample(t[:-1],X[:-1,idx])
    prefix_state,base,amp=baseline_scale(prefix_raw)
    full_raw,_=resample(t,X[:,idx],med)
    full_state,_,_=baseline_scale(full_raw,base,amp)
    prefix_pc=pca.transform(prefix_state).ravel(); true_pc=float(pca.transform(full_state)[-1,0])
    # Cross-dataset agreement: held-out observed prefix vs shared training curve.
    prefix_shared=np.interp(normalize_time(t[:-1]),GRID,shared)
    agreement_p=corr(prefix_pc,prefix_shared); agreement_s=corr(prefix_pc,prefix_shared,"spearman")
    # Dataset-specific residual and a simple residual-persistence predictor.
    residual=prefix_pc-shared
    residual_final=float(residual[-1])
    pred_shared=float(shared[-1]); pred_dataset=pred_shared+residual_final
    persistence=float(prefix_pc[-1])
    nearest=float(np.interp(normalize_time(t[:-1])[-1],GRID,shared))
    # Biological association is measured on the held-out prefix only, with
    # fixed gene panels defined independently of this fold.
    acts=program_activity(prefix_state,genes)
    biological=[]
    for pid,a in acts.items():
        biological.append({"program_id":pid,"component_program_pearson":corr(prefix_pc,a),"component_program_spearman":corr(prefix_pc,a,"spearman")})
    return {"heldout_dataset":heldout,"n_selected_genes":len(genes),"selected_gene_score_mean":float(np.mean(sel_scores)),"pca":pca,"shared":shared,"prefix_pc":prefix_pc,"true_pc":true_pc,"agreement_pearson":agreement_p,"agreement_spearman":agreement_s,"residual_rms":float(np.sqrt(np.mean(residual**2))),"residual_endpoint_abs":abs(residual_final),"pred_shared":pred_shared,"pred_dataset":pred_dataset,"persistence":persistence,"nearest":nearest,"biological":biological}

def evaluate(a):
    true=a["true_pc"]
    rmse_d=abs(a["pred_dataset"]-true); rmse_s=abs(a["pred_shared"]-true); rmse_p=abs(a["persistence"]-true); rmse_n=abs(a["nearest"]-true)
    return {"heldout_dataset":a["heldout_dataset"],"n_selected_genes":a["n_selected_genes"],"selected_gene_score_mean":a["selected_gene_score_mean"],"agreement_pearson":a["agreement_pearson"],"agreement_spearman":a["agreement_spearman"],"residual_rms":a["residual_rms"],"residual_endpoint_abs":a["residual_endpoint_abs"],"abs_error_shared_plus_dataset":rmse_d,"abs_error_shared_only":rmse_s,"abs_error_persistence":rmse_p,"abs_error_nearest_shared":rmse_n,"improvement_vs_persistence":rmse_p-rmse_d,"improvement_vs_shared_only":rmse_s-rmse_d,"improvement_vs_nearest":rmse_n-rmse_d}

def permutation_null(artifacts,n=N_PERM):
    rows=[]
    for b in range(n):
        vals=[]
        for di,(ds,a) in enumerate(artifacts.items()):
            rng=np.random.default_rng(932000+b*17+di); perm=a["prefix_pc"][rng.permutation(len(a["prefix_pc"]))]
            pred=float(a["shared"][-1]+(np.mean(perm-a["shared"])))
            err=abs(pred-a["true_pc"]); per=abs(a["persistence"]-a["true_pc"]); vals.append(per-err)
        rows.append({"permutation":b+1,"mean_improvement_vs_persistence":float(np.mean(vals)) if vals else np.nan})
    return pd.DataFrame(rows)

def run():
    traj=load_trajectories(); names=[d for d in TARGET if d in traj]
    log(f"trajectory datasets: {', '.join(names)}")
    if len(names)<3: raise RuntimeError("Stage 2.9.32 requires at least three trajectory datasets.")
    artifacts={}; rows=[]; bio_rows=[]
    for held in names:
        train=[d for d in names if d!=held]; log(f"LODO held out {held}; training on {', '.join(train)}")
        a=prepare_fold(held,train,traj)
        if a is None: continue
        artifacts[held]=a; rows.append(evaluate(a))
        for r in a["biological"]: bio_rows.append({"heldout_dataset":held,**r})
    R=pd.DataFrame(rows)
    if R.empty: raise RuntimeError("No valid Stage 2.9.32 LODO folds.")
    R.to_csv(OUT/"01_lodo_shared_biological_invariance.csv",index=False)
    B=pd.DataFrame(bio_rows); B.to_csv(OUT/"02_biological_program_association.csv",index=False)
    P=permutation_null(artifacts);P.to_csv(OUT/"03_time_permutation_null.csv",index=False)
    obs=float(R.improvement_vs_persistence.mean()); nv=P.mean_improvement_vs_persistence.dropna().to_numpy(float)
    p=float((1+np.sum(nv>=obs))/(len(nv)+1)) if len(nv) else np.nan
    # Invariance requires positive cross-dataset agreement, small residual
    # relative to shared signal, and reproducible biological association.
    mean_agree=float(R.agreement_spearman.mean()); mean_resid=float(R.residual_rms.mean())
    bio_strength=float(B.component_program_spearman.abs().mean()) if len(B) else np.nan
    support=bool(mean_agree>0.5 and obs>0 and np.isfinite(p) and p<0.05 and bio_strength>0.3)
    summary=pd.DataFrame([{"n_trajectory_datasets":len(names),"n_valid_lodo_folds":len(R),"mean_selected_genes":float(R.n_selected_genes.mean()),"mean_agreement_spearman":mean_agree,"mean_agreement_pearson":float(R.agreement_pearson.mean()),"mean_residual_rms":mean_resid,"mean_residual_endpoint_abs":float(R.residual_endpoint_abs.mean()),"mean_biological_component_spearman_abs":bio_strength,"mean_abs_error_shared_plus_dataset":float(R.abs_error_shared_plus_dataset.mean()),"mean_abs_error_persistence":float(R.abs_error_persistence.mean()),"mean_improvement_vs_persistence":obs,"mean_improvement_vs_shared_only":float(R.improvement_vs_shared_only.mean()),"time_permutation_p":p,"shared_biological_invariance_supported":support,"stage3_readiness":False,"interpretation":"LODO decomposition into training-derived shared temporal component and held-out dataset-specific residual; fixed biological programs test whether the shared component is more than a temporal clock; predictive gate remains closed unless invariance, biology, prediction and permutation criteria all pass."}])
    summary.to_csv(OUT/"04_stage2932_summary.csv",index=False)
    log("overall:"); print(summary.to_string(index=False)); return summary

if __name__=="__main__": run()
