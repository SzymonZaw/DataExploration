"""Unified predictive benchmark for learned cellular-state dynamics."""
from __future__ import annotations
from dataclasses import dataclass
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

def set_seed(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)

def training_statistics(train, max_genes):
    arrays=[X for _,X in train.values()]; raw=np.vstack(arrays)
    med=np.nanmedian(np.where(np.isfinite(raw),raw,np.nan),axis=0); med=np.where(np.isfinite(med),med,0.0)
    pooled=np.vstack([np.where(np.isfinite(X),X,med) for X in arrays]); var=np.var(pooled,axis=0)
    keep=np.argsort(var)[::-1][:min(max_genes,pooled.shape[1])]
    mean=pooled[:,keep].mean(axis=0); std=np.where(pooled[:,keep].std(axis=0)>1e-8,pooled[:,keep].std(axis=0),1.0)
    return keep,mean,std

def transform(X,keep,mean,std):
    X=X[:,keep]; return (np.where(np.isfinite(X),X,mean)-mean)/std

def chronological_examples(data,keep,mean,std,history_len=0):
    xs=[];ys=[];dts=[];histories=[]
    for times,Xraw in data.values():
        X=transform(Xraw,keep,mean,std); scale=max(float(times[-1]-times[0]),1.0); tn=(times-times[0])/scale
        for i in range(history_len,len(times)-1):
            xs.append(X[i]);ys.append(X[i+1]);dts.append([float(tn[i+1]-tn[i])])
            if history_len: histories.append(X[i-history_len:i])
    return (torch.tensor(np.asarray(xs),dtype=torch.float32),torch.tensor(np.asarray(ys),dtype=torch.float32),torch.tensor(np.asarray(dts),dtype=torch.float32),torch.tensor(np.asarray(histories),dtype=torch.float32) if history_len else None)

def train_neural(train_data,keep,mean,std,cfg,seed,mode):
    set_seed(seed); hist_len=cfg.history_len if mode=="memory" else 0
    x,y,dt,h_raw=chronological_examples(train_data,keep,mean,std,hist_len); context_dim=1 if mode in {"delta_t","memory"} else 0; history_dim=cfg.history_dim if mode=="memory" else 0
    model=DynamicStateModel(x.shape[1],cfg.state_dim,cfg.hidden_dim,context_dim,history_dim,cfg.dropout)
    opt=torch.optim.AdamW(model.parameters(),lr=cfg.lr,weight_decay=1e-4); best=float("inf");stale=0;patience=max(25,cfg.epochs//8)
    for _ in range(cfg.epochs):
        model.train();opt.zero_grad(set_to_none=True);z=model.encode(x);zn=model.encode(y).detach();rec=model.decode(z)
        if history_dim:
            b,t,d=h_raw.shape;hist_lat=model.encode(h_raw.reshape(b*t,d)).reshape(b,t,cfg.state_dim);pred=model.transition(z,context=dt,history=hist_lat)
        else: pred=model.transition(z,context=dt if context_dim else None)
        loss=((rec-x)**2).mean()+.25*((pred-zn)**2).mean()+((model.decode(pred)-y)**2).mean();loss.backward();nn.utils.clip_grad_norm_(model.parameters(),1.0);opt.step();value=float(loss.detach())
        if value<best-1e-6:best=value;stale=0
        else:stale+=1
        if stale>=patience:break
    return model,best

def predict_step(model,current,dt,history,mode):
    x=torch.tensor(current[None,:],dtype=torch.float32);ctx=torch.tensor([[float(dt)]],dtype=torch.float32) if mode in {"delta_t","memory"} else None;h=None
    if mode=="memory":
        raw=torch.tensor(history,dtype=torch.float32)
        if raw.ndim!=3: raise ValueError(f"Memory history must be [batch,time,features], got shape {tuple(raw.shape)}")
        b,t,d=raw.shape;model.eval()
        with torch.no_grad():h=model.encode(raw.reshape(b*t,d)).reshape(b,t,-1)
    model.eval()
    with torch.no_grad():return model.predict_observation(x,context=ctx,history=h).numpy()[0]

def evaluate_neural(model,train_data,heldout_data,keep,mean,std,cfg,mode):
    times,Xraw=heldout_data;X=transform(Xraw,keep,mean,std);n_prefix=min(max(2,cfg.history_len+1,int(np.ceil(len(times)*cfg.prefix_fraction))),len(times)-1);prefix_t=list(times[:n_prefix]);observed=[v.copy() for v in X[:n_prefix]];current=observed[-1].copy();current_t=float(prefix_t[-1]);true=[];pred=[];persistence=[];nearest=[];linear=[];train_trans={k:transform(v,keep,mean,std) for k,(_,v) in train_data.items()};scale=max(float(times[-1]-times[0]),1.0)
    for j in range(n_prefix,len(times)):
        dt=(float(times[j])-current_t)/scale
        hist=np.asarray(observed[-cfg.history_len:])[None,:,:] if mode=="memory" else None
        p=predict_step(model,current,dt,hist,mode);pred.append(p);persistence.append(observed[-1]);true.append(X[j]);nearest.append(np.mean([train_trans[ds][int(np.argmin(np.abs(tt-times[j])))] for ds,(tt,_) in train_data.items()],axis=0));linear.append(observed[-1]+((times[j]-prefix_t[-1])/(prefix_t[-1]-prefix_t[-2]))*(observed[-1]-observed[-2]) if len(prefix_t)>=2 and prefix_t[-1]!=prefix_t[-2] else observed[-1]);current=p;current_t=float(times[j]);observed.append(p)
    rm=lambda a,b:float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)));out={"n_future_points":len(true),"rmse_model":rm(true,pred),"rmse_persistence":rm(true,persistence),"rmse_nearest":rm(true,nearest),"rmse_linear":rm(true,linear)};out["improvement_vs_persistence"]=out["rmse_persistence"]-out["rmse_model"];out["improvement_vs_nearest"]=out["rmse_nearest"]-out["rmse_model"];out["improvement_vs_linear"]=out["rmse_linear"]-out["rmse_model"];return out

def evaluate_pca(train_data,heldout_data,keep,mean,std,cfg):
    train_z=np.vstack([transform(X,keep,mean,std) for _,X in train_data.values()]);scaler=StandardScaler().fit(train_z);pca=PCA(n_components=min(cfg.state_dim,train_z.shape[0],train_z.shape[1])).fit(scaler.transform(train_z));times,Xraw=heldout_data;Z=pca.transform(scaler.transform(transform(Xraw,keep,mean,std)));n_prefix=min(max(2,int(np.ceil(len(times)*cfg.prefix_fraction)),len(times)-1),len(times)-1);true=Z[n_prefix:];pers=np.repeat(Z[n_prefix-1][None,:],len(true),axis=0);linear=[Z[n_prefix-1]+((times[j]-times[n_prefix-1])/(times[n_prefix-1]-times[n_prefix-2]))*(Z[n_prefix-1]-Z[n_prefix-2]) if times[n_prefix-1]!=times[n_prefix-2] else Z[n_prefix-1] for j in range(n_prefix,len(times))];linear=np.asarray(linear);rm=lambda a,b:float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)));a=rm(true,pers);b=rm(true,linear);return {"n_future_points":len(true),"rmse_model":b,"rmse_persistence":a,"rmse_nearest":np.nan,"rmse_linear":b,"improvement_vs_persistence":a-b,"improvement_vs_nearest":np.nan,"improvement_vs_linear":0.0}

def evaluate_autoencoder(train_data,heldout_data,keep,mean,std,cfg,seed): return evaluate_neural(train_neural(train_data,keep,mean,std,cfg,seed,"markov")[0],train_data,heldout_data,keep,mean,std,cfg,"markov")
def benchmark_fold(train_data,heldout_data,cfg,seed,heldout_name):
    keep,mean,std=training_statistics(train_data,cfg.max_genes);rows=[]
    for mode in ("pca","autoencoder","markov","delta_t","memory"):
        if mode=="pca":r=evaluate_pca(train_data,heldout_data,keep,mean,std,cfg)
        elif mode=="autoencoder":r=evaluate_autoencoder(train_data,heldout_data,keep,mean,std,cfg,seed)
        else:
            model,loss=train_neural(train_data,keep,mean,std,cfg,seed,mode);r=evaluate_neural(model,train_data,heldout_data,keep,mean,std,cfg,mode);r["train_loss"]=loss
        rows.append({"seed":seed,"heldout_dataset":heldout_name,"model":mode,"n_selected_genes":len(keep),**r})
    return rows
