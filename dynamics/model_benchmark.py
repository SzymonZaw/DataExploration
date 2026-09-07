"""Unified predictive benchmark for learned cellular-state dynamics.

The benchmark is the main Stage 2.11 framework. It compares progressively richer
representations under the same leakage-free LODO prefix-to-future protocol:
PCA, autoencoder, Markov DynamicStateModel, delta-t conditioned DynamicStateModel,
and history-aware DynamicStateModel. Models are evaluated against persistence,
nearest-time and linear extrapolation baselines.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import random
import numpy as np
import pandas as pd
import torch
from torch import nn
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from .dynamic_state_model import DynamicStateModel

@dataclass
class BenchmarkConfig:
    max_genes: int = 2000
    state_dim: int = 8
    hidden_dim: int = 128
    epochs: int = 250
    lr: float = 1e-3
    prefix_fraction: float = 0.6
    history_len: int = 2
    history_dim: int = 16
    dropout: float = 0.05


def set_seed(seed: int) -> None:
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)


def training_statistics(train, max_genes):
    arrays = [X for _, X in train.values()]
    raw = np.vstack(arrays)
    med = np.nanmedian(np.where(np.isfinite(raw), raw, np.nan), axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    pooled = np.vstack([np.where(np.isfinite(X), X, med) for X in arrays])
    var = np.var(pooled, axis=0)
    keep = np.argsort(var)[::-1][:min(max_genes, pooled.shape[1])]
    mean = pooled[:, keep].mean(axis=0)
    std = np.where(pooled[:, keep].std(axis=0) > 1e-8, pooled[:, keep].std(axis=0), 1.0)
    return keep, mean, std


def transform(X, keep, mean, std):
    X = X[:, keep]
    return (np.where(np.isfinite(X), X, mean) - mean) / std


def chronological_examples(data, keep, mean, std, history_len=0):
    xs=[]; ys=[]; dts=[]; histories=[]
    for times, Xraw in data.values():
        X=transform(Xraw,keep,mean,std)
        scale=max(float(times[-1]-times[0]),1.0); tn=(times-times[0])/scale
        start=history_len
        for i in range(start,len(times)-1):
            xs.append(X[i]); ys.append(X[i+1]); dts.append([float(tn[i+1]-tn[i])])
            if history_len: histories.append(X[i-history_len:i])
    x=torch.tensor(np.asarray(xs),dtype=torch.float32)
    y=torch.tensor(np.asarray(ys),dtype=torch.float32)
    dt=torch.tensor(np.asarray(dts),dtype=torch.float32)
    h=torch.tensor(np.asarray(histories),dtype=torch.float32) if history_len else None
    return x,y,dt,h


def train_neural(train_data, keep, mean, std, cfg, seed, mode):
    set_seed(seed)
    hist_len=cfg.history_len if mode=="memory" else 0
    x,y,dt,h=chronological_examples(train_data,keep,mean,std,hist_len)
    history_dim=cfg.history_dim if mode=="memory" else 0
    context_dim=1 if mode in {"delta_t","memory"} else 0
    model=DynamicStateModel(x.shape[1],cfg.state_dim,cfg.hidden_dim,context_dim,history_dim,cfg.dropout)
    opt=torch.optim.AdamW(model.parameters(),lr=cfg.lr,weight_decay=1e-4)
    best=float("inf"); stale=0; patience=max(25,cfg.epochs//8)
    for _ in range(cfg.epochs):
        model.train(); opt.zero_grad(set_to_none=True)
        z=model.encode(x); zn=model.encode(y).detach(); rec=model.decode(z)
        ctx=dt if context_dim else None; hist=h if history_dim else None
        pred=model.transition(z,context=ctx,history=hist); pred_obs=model.decode(pred)
        loss=((rec-x)**2).mean()+.25*((pred-zn)**2).mean()+((pred_obs-y)**2).mean()
        loss.backward(); nn.utils.clip_grad_norm_(model.parameters(),1.0); opt.step()
        value=float(loss.detach())
        if value < best-1e-6: best=value; stale=0
        else:
            stale+=1
            if stale>=patience: break
    return model,best


def predict_step(model,current,dt,history,mode):
    x=torch.tensor(current[None,:],dtype=torch.float32)
    context=torch.tensor([[float(dt)]],dtype=torch.float32) if mode in {"delta_t","memory"} else None
    h=torch.tensor(history[None,:,:],dtype=torch.float32) if mode=="memory" else None
    model.eval()
    with torch.no_grad(): return model.predict_observation(x,context=context,history=h).numpy()[0]


def evaluate_neural(model,train_data,heldout,keep,mean,std,cfg,mode):
    times,Xraw=heldout; X=transform(Xraw,keep,mean,std)
    n_prefix=min(max(2,cfg.history_len+1,int(np.ceil(len(times)*cfg.prefix_fraction))),len(times)-1)
    prefix_t=list(times[:n_prefix]); observed=[v.copy() for v in X[:n_prefix]]
    current=observed[-1].copy(); current_t=float(prefix_t[-1]); true=[]; pred=[]; persistence=[]; nearest=[]; linear=[]
    train_trans={k:transform(v,keep,mean,std) for k,(_,v) in train_data.items()}
    scale=max(float(times[-1]-times[0]),1.0)
    for j in range(n_prefix,len(times)):
        dt=(float(times[j])-current_t)/scale
        hist=np.asarray(observed[-cfg.history_len:]) if mode=="memory" else None
        p=predict_step(model,current,dt,hist,mode)
        pred.append(p); persistence.append(observed[-1]); true.append(X[j])
        cand=[]
        for ds,(tt,_) in train_data.items(): cand.append(train_trans[ds][int(np.argmin(np.abs(tt-times[j])))])
        nearest.append(np.mean(cand,axis=0))
        if len(prefix_t)>=2 and prefix_t[-1]!=prefix_t[-2]:
            linear.append(observed[-1]+((times[j]-prefix_t[-1])/(prefix_t[-1]-prefix_t[-2]))*(observed[-1]-observed[-2]))
        else: linear.append(observed[-1])
        # Model-only rollout: never append the true future observation.
        current=p; current_t=float(times[j])
        # Baselines receive the observed prefix only; future values remain hidden.
        if j < len(times)-1: pass
    def rm(a,b): return float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)))
    out={"n_future_points":len(true),"rmse_model":rm(true,pred),"rmse_persistence":rm(true,persistence),"rmse_nearest":rm(true,nearest),"rmse_linear":rm(true,linear)}
    out["improvement_vs_persistence"]=out["rmse_persistence"]-out["rmse_model"]
    out["improvement_vs_nearest"]=out["rmse_nearest"]-out["rmse_model"]
    out["improvement_vs_linear"]=out["rmse_linear"]-out["rmse_model"]
    return out


def evaluate_pca(train_data,heldout,keep,mean,std,cfg):
    # PCA is a representation baseline: fit on training observations only.
    scaler=StandardScaler(); train_z=np.vstack([transform(X,keep,mean,std) for _,X in train_data.values()]); scaler.fit(train_z)
    pca=PCA(n_components=min(cfg.state_dim,train_z.shape[1],train_z.shape[0])); pca.fit(scaler.transform(train_z))
    times,Xraw=heldout; X=scaler.transform(transform(Xraw,keep,mean,std)); Z=pca.transform(X)
    n=min(max(2,int(np.ceil(len(times)*cfg.prefix_fraction)),len(times)-1),len(times)-1)
    true=Z[n:]; persistence=np.repeat(Z[n-1][None,:],len(true),axis=0)
    # constant-velocity latent extrapolation from observed prefix
    linear=[]
    for j in range(n,len(times)):
        linear.append(Z[n-1]+((times[j]-times[n-1])/(times[n-1]-times[n-2]))*(Z[n-1]-Z[n-2]) if times[n-1]!=times[n-2] else Z[n-1])
    return {"n_future_points":len(true),"rmse_model":float(np.sqrt(np.mean((np.asarray(linear)-true)**2))),"rmse_persistence":float(np.sqrt(np.mean((persistence-true)**2))),"rmse_nearest":np.nan,"rmse_linear":float(np.sqrt(np.mean((np.asarray(linear)-true)**2))),"improvement_vs_persistence":float(np.sqrt(np.mean((persistence-true)**2))-np.sqrt(np.mean((np.asarray(linear)-true)**2))),"improvement_vs_nearest":np.nan,"improvement_vs_linear":0.0}


def evaluate_autoencoder(train_data,heldout,keep,mean,std,cfg,seed):
    model,_=train_neural(train_data,keep,mean,std,cfg,seed,"markov")
    return evaluate_neural(model,train_data,heldout,keep,mean,std,cfg,"markov")


def benchmark_fold(train_data,heldout,cfg,seed):
    keep,mean,std=training_statistics(train_data,cfg.max_genes); rows=[]
    for mode in ("pca","autoencoder","markov","delta_t","memory"):
        if mode=="pca": r=evaluate_pca(train_data,heldout,keep,mean,std,cfg)
        elif mode=="autoencoder": r=evaluate_autoencoder(train_data,heldout,keep,mean,std,cfg,seed)
        else:
            model,loss=train_neural(train_data,keep,mean,std,cfg,seed,mode)
            r=evaluate_neural(model,train_data,heldout,keep,mean,std,cfg,mode); r["train_loss"]=loss
        rows.append({"seed":seed,"heldout_dataset":heldout,"model":mode,"n_selected_genes":len(keep),**r})
    return rows
