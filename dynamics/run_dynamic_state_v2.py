"""Validated DynamicStateModel v2 experiment.

LODO prefix-to-future forecasting with training-only preprocessing, explicit delta-t,
multiple seeds, strong baselines and a time-permutation null. This module is a
predictive benchmark, not a biological or causal claim.
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
from .dynamic_state_model import DynamicStateModel
from .validation import _load_common_space

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"dynamic_state_model_v2"
DATASETS=["GSE67462","GSE28688","GSE297234"]
DEFAULT_SEEDS=(211,212,213,214,215)


def set_seed(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)


def load_data():
    matrix,meta=_load_common_space(); data={}
    for ds in DATASETS:
        g=meta[(meta["dataset"]==ds)&meta["time_hours"].notna()&meta["matrix_column"].notna()].copy()
        if g["time_hours"].nunique()<3: continue
        g=g.sort_values("time_hours")
        X=matrix.loc[:,g["matrix_column"].astype(str).tolist()].T
        X.index=g["time_hours"].to_numpy(float)
        X=X.groupby(level=0,sort=True).mean()
        data[ds]=(X.index.to_numpy(float),X.to_numpy(float))
    if len(data)<2: raise RuntimeError("At least two longitudinal datasets are required.")
    return data


def training_statistics(train,max_genes):
    arrays=[X for _,X in train.values()]
    raw=np.vstack(arrays)
    med=np.nanmedian(np.where(np.isfinite(raw),raw,np.nan),axis=0)
    med=np.where(np.isfinite(med),med,0.0)
    pooled=np.vstack([np.where(np.isfinite(X),X,med) for X in arrays])
    var=np.var(pooled,axis=0)
    keep=np.argsort(var)[::-1][:min(max_genes,pooled.shape[1])]
    mean=pooled[:,keep].mean(axis=0); std=pooled[:,keep].std(axis=0)
    return keep,mean,np.where(std>1e-8,std,1.0)


def transform(X,keep,mean,std):
    X=X[:,keep]; X=np.where(np.isfinite(X),X,mean); return (X-mean)/std


def transition_pairs(train,keep,mean,std):
    xs,ys,dts=[],[],[]
    for times,X in train.values():
        Z=transform(X,keep,mean,std); scale=max(float(times[-1]-times[0]),1.0)
        tn=(times-times[0])/scale
        for i in range(len(times)-1):
            xs.append(Z[i]); ys.append(Z[i+1]); dts.append([float(tn[i+1]-tn[i])])
    return tuple(torch.tensor(np.asarray(v),dtype=torch.float32) for v in (xs,ys,dts))


def train_model(train_data,keep,mean,std,state_dim,hidden_dim,epochs,lr,seed):
    set_seed(seed); x,y,dt=transition_pairs(train_data,keep,mean,std)
    model=DynamicStateModel(input_dim=x.shape[1],state_dim=state_dim,hidden_dim=hidden_dim,context_dim=1,dropout=.05)
    opt=torch.optim.AdamW(model.parameters(),lr=lr,weight_decay=1e-4); best=float("inf"); stale=0
    patience=max(20,epochs//8)
    for _ in range(epochs):
        model.train(); opt.zero_grad(set_to_none=True)
        z=model.encode(x); z_next=model.encode(y).detach(); rec=model.decode(z)
        pred=model.transition(z,context=dt); pred_obs=model.decode(pred)
        loss=((rec-x)**2).mean()+.25*((pred-z_next)**2).mean()+((pred_obs-y)**2).mean()
        loss.backward(); nn.utils.clip_grad_norm_(model.parameters(),1.0); opt.step()
        value=float(loss.detach())
        if value<best-1e-6: best=value; stale=0
        else:
            stale+=1
            if stale>=patience: break
    return model,best


def predict_step(model,current,dt):
    x=torch.tensor(current[None,:],dtype=torch.float32); d=torch.tensor([[float(dt)]],dtype=torch.float32)
    model.eval()
    with torch.no_grad(): return model.predict_observation(x,context=d).numpy()[0]


def nearest_prediction(train_data,target_t,keep,mean,std):
    candidates=[]
    for times,X in train_data.values():
        idx=int(np.argmin(np.abs(times-target_t))); candidates.append(transform(X,keep,mean,std)[idx])
    return np.mean(candidates,axis=0)


def linear_prediction(prefix_t,prefix_X,target_t):
    if len(prefix_t)<2: return prefix_X[-1].copy()
    t0,t1=float(prefix_t[-2]),float(prefix_t[-1])
    if t1==t0: return prefix_X[-1].copy()
    return prefix_X[-1]+((float(target_t)-t1)/(t1-t0))*(prefix_X[-1]-prefix_X[-2])


def metrics(y,p):
    e=np.asarray(y)-np.asarray(p); return float(np.sqrt(np.mean(e**2))),float(np.mean(np.abs(e)))


def evaluate_prefix(model,train_data,heldout,keep,mean,std,prefix_fraction):
    times,Xraw=heldout; X=transform(Xraw,keep,mean,std)
    n_prefix=min(max(2,int(np.ceil(len(times)*prefix_fraction))),len(times)-1)
    prefix_t=times[:n_prefix]; prefix_X=X[:n_prefix]
    true=[]; pm=[]; pp=[]; pn=[]; pl=[]
    scale=max(float(times[-1]-times[0]),1.0)
    # Genuine multi-step rollout: future observations are never fed back into the model.
    current=prefix_X[-1].copy(); current_t=float(prefix_t[-1])
    for j in range(n_prefix,len(times)):
        dt=(float(times[j])-current_t)/scale
        pm.append(predict_step(model,current,dt)); pp.append(prefix_X[-1])
        pn.append(nearest_prediction(train_data,float(times[j]),keep,mean,std)); pl.append(linear_prediction(prefix_t,prefix_X,float(times[j])))
        true.append(X[j]); current=pm[-1]; current_t=float(times[j])
    if not true: return None
    out={"n_future_points":len(true)}
    for name,p in (("model",pm),("persistence",pp),("nearest",pn),("linear",pl)):
        out[f"rmse_{name}"],out[f"mae_{name}"]=metrics(true,p)
    out["improvement_vs_persistence"]=out["rmse_persistence"]-out["rmse_model"]
    out["improvement_vs_nearest"]=out["rmse_nearest"]-out["rmse_model"]
    out["improvement_vs_linear"]=out["rmse_linear"]-out["rmse_model"]
    return out


def permutation_null(data,max_genes,state_dim,hidden_dim,epochs,lr,seed,prefix_fraction,n_perm):
    rng=np.random.default_rng(seed); vals=[]; perm_epochs=min(50,max(20,epochs//5))
    for k in range(n_perm):
        permuted={ds:(t,X[rng.permutation(len(t))]) for ds,(t,X) in data.items()}; fold_vals=[]
        for held_name in sorted(permuted):
            train={k:v for k,v in permuted.items() if k!=held_name}
            keep,mean,std=training_statistics(train,max_genes)
            model,_=train_model(train,keep,mean,std,state_dim,hidden_dim,perm_epochs,lr,seed+1000+k)
            r=evaluate_prefix(model,train,permuted[held_name],keep,mean,std,prefix_fraction)
            if r: fold_vals.append(r["improvement_vs_persistence"])
        if fold_vals: vals.append(float(np.mean(fold_vals)))
        if (k+1)%10==0: print(f"DynamicStateModel v2 permutation {k+1}/{n_perm}",flush=True)
    return np.asarray(vals)


def run(max_genes=2000,state_dim=8,hidden_dim=128,epochs=250,lr=1e-3,seeds=DEFAULT_SEEDS,prefix_fraction=.6,n_perm=100):
    OUT.mkdir(parents=True,exist_ok=True); data=load_data(); rows=[]
    for seed in seeds:
        for held_name in sorted(data):
            train={k:v for k,v in data.items() if k!=held_name}; keep,mean,std=training_statistics(train,max_genes)
            model,loss=train_model(train,keep,mean,std,state_dim,hidden_dim,epochs,lr,seed)
            r=evaluate_prefix(model,train,data[held_name],keep,mean,std,prefix_fraction)
            if r:
                rows.append({"seed":seed,"heldout_dataset":held_name,"n_training_datasets":len(train),"n_selected_genes":len(keep),"train_loss":loss,**r})
                print(f"DynamicStateModel v2: seed={seed} heldout={held_name} RMSE={r['rmse_model']:.4f} persistence={r['rmse_persistence']:.4f} nearest={r['rmse_nearest']:.4f} linear={r['rmse_linear']:.4f}",flush=True)
    df=pd.DataFrame(rows); df.to_csv(OUT/"01_lodo_prefix_metrics.csv",index=False)
    if df.empty: raise RuntimeError("No valid LODO prefix-forecast folds.")
    obs=float(df.improvement_vs_persistence.mean())
    summary=pd.DataFrame([{"n_valid_folds":len(df),"n_seeds":df.seed.nunique(),"mean_rmse_model":df.rmse_model.mean(),"mean_rmse_persistence":df.rmse_persistence.mean(),"mean_rmse_nearest":df.rmse_nearest.mean(),"mean_rmse_linear":df.rmse_linear.mean(),"mean_improvement_vs_persistence":obs,"mean_improvement_vs_nearest":df.improvement_vs_nearest.mean(),"mean_improvement_vs_linear":df.improvement_vs_linear.mean(),"median_improvement_vs_persistence":df.improvement_vs_persistence.median(),"q05_improvement_vs_persistence":df.improvement_vs_persistence.quantile(.05),"q95_improvement_vs_persistence":df.improvement_vs_persistence.quantile(.95)}])
    null=permutation_null(data,max_genes,state_dim,hidden_dim,epochs,lr,seeds[0],prefix_fraction,n_perm)
    p=float((1+np.sum(null>=obs))/(1+len(null))) if len(null) else np.nan
    summary["permutation_p_improvement_vs_persistence"]=p
    summary["dynamic_state_predictive_support"]=bool(obs>0 and summary.loc[0,"mean_improvement_vs_nearest"]>0 and summary.loc[0,"mean_improvement_vs_linear"]>0 and p<.05)
    summary.to_csv(OUT/"02_summary.csv",index=False); pd.DataFrame({"null_improvement_vs_persistence":null}).to_csv(OUT/"03_permutation_null.csv",index=False)
    print(summary.to_string(index=False),flush=True); return summary

if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--max-genes",type=int,default=2000); p.add_argument("--state-dim",type=int,default=8); p.add_argument("--hidden-dim",type=int,default=128); p.add_argument("--epochs",type=int,default=250); p.add_argument("--lr",type=float,default=1e-3); p.add_argument("--prefix-fraction",type=float,default=.6); p.add_argument("--permutations",type=int,default=100); p.add_argument("--seeds",type=int,nargs="+",default=list(DEFAULT_SEEDS)); a=p.parse_args(); run(a.max_genes,a.state_dim,a.hidden_dim,a.epochs,a.lr,tuple(a.seeds),a.prefix_fraction,a.permutations)
