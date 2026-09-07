"""Unified predictive benchmark for cellular-state dynamics.

All methods share the same leakage-free LODO folds, training-only preprocessing,
held-out prefix, recursive rollout, baselines and permutation protocol.
"""
from __future__ import annotations
from dataclasses import dataclass
import random
import numpy as np
import pandas as pd
import torch
from torch import nn
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
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
    permutation_n: int = 1000

class StaticAutoencoder(nn.Module):
    """Non-dynamical AE baseline; latent dynamics are fitted separately."""
    def __init__(self, input_dim, state_dim, hidden_dim, dropout):
        super().__init__()
        self.encoder=nn.Sequential(nn.Linear(input_dim,hidden_dim),nn.GELU(),nn.Dropout(dropout),nn.Linear(hidden_dim,state_dim))
        self.decoder=nn.Sequential(nn.Linear(state_dim,hidden_dim),nn.GELU(),nn.Linear(hidden_dim,input_dim))
    def encode(self,x): return self.encoder(x)
    def decode(self,z): return self.decoder(z)
    def forward(self,x):
        z=self.encode(x);return z,self.decode(z)

def set_seed(seed):
    random.seed(seed);np.random.seed(seed);torch.manual_seed(seed)

def training_statistics(train,max_genes):
    arrays=[X for _,X in train.values()];raw=np.vstack(arrays)
    safe=np.where(np.isfinite(raw),raw,np.nan);med=np.nanmedian(safe,axis=0);med=np.where(np.isfinite(med),med,0.0)
    pooled=np.vstack([np.where(np.isfinite(X),X,med) for X in arrays]);var=np.var(pooled,axis=0)
    keep=np.argsort(var)[::-1][:min(max_genes,pooled.shape[1])];mean=pooled[:,keep].mean(axis=0);std=np.where(pooled[:,keep].std(axis=0)>1e-8,pooled[:,keep].std(axis=0),1.0)
    return keep,mean,std

def transform(X,keep,mean,std):
    X=X[:,keep];return (np.where(np.isfinite(X),X,mean)-mean)/std

def chronological_examples(data,keep,mean,std,history_len=0):
    xs=[];ys=[];dts=[];histories=[]
    for times,Xraw in data.values():
        X=transform(Xraw,keep,mean,std);scale=max(float(times[-1]-times[0]),1.0);tn=(times-times[0])/scale
        for i in range(history_len,len(times)-1):
            xs.append(X[i]);ys.append(X[i+1]);dts.append([float(tn[i+1]-tn[i])])
            if history_len:histories.append(X[i-history_len:i])
    return np.asarray(xs),np.asarray(ys),np.asarray(dts),np.asarray(histories) if history_len else None

def train_dynamic(train_data,keep,mean,std,cfg,seed,mode):
    set_seed(seed);hist_len=cfg.history_len if mode=="memory" else 0
    xa,ya,dta,ha=chronological_examples(train_data,keep,mean,std,hist_len)
    x=torch.tensor(xa,dtype=torch.float32);y=torch.tensor(ya,dtype=torch.float32);dt=torch.tensor(dta,dtype=torch.float32)
    context_dim=1 if mode in {"delta_t","memory"} else 0;history_dim=cfg.history_dim if mode=="memory" else 0
    model=DynamicStateModel(x.shape[1],cfg.state_dim,cfg.hidden_dim,context_dim,history_dim,cfg.dropout)
    opt=torch.optim.AdamW(model.parameters(),lr=cfg.lr,weight_decay=1e-4);best=float("inf");stale=0;patience=max(25,cfg.epochs//8)
    for _ in range(cfg.epochs):
        model.train();opt.zero_grad(set_to_none=True);z=model.encode(x);zn=model.encode(y).detach();rec=model.decode(z)
        if history_dim:
            b,t,d=ha.shape;hraw=torch.tensor(ha,dtype=torch.float32);hl=model.encode(hraw.reshape(b*t,d)).reshape(b,t,cfg.state_dim);pred=model.transition(z,context=dt,history=hl)
        else:pred=model.transition(z,context=dt if context_dim else None)
        loss=((rec-x)**2).mean()+.25*((pred-zn)**2).mean()+((model.decode(pred)-y)**2).mean();loss.backward();nn.utils.clip_grad_norm_(model.parameters(),1.0);opt.step();value=float(loss.detach())
        if value<best-1e-6:best=value;stale=0
        else:stale+=1
        if stale>=patience:break
    return model,best

def train_autoencoder(train_data,keep,mean,std,cfg,seed):
    set_seed(seed);xa,_,_,_=chronological_examples(train_data,keep,mean,std,0);x=torch.tensor(xa,dtype=torch.float32)
    model=StaticAutoencoder(x.shape[1],cfg.state_dim,cfg.hidden_dim,cfg.dropout);opt=torch.optim.AdamW(model.parameters(),lr=cfg.lr,weight_decay=1e-4)
    for _ in range(cfg.epochs):
        opt.zero_grad(set_to_none=True);z,rec=model(x);loss=((rec-x)**2).mean();loss.backward();nn.utils.clip_grad_norm_(model.parameters(),1.0);opt.step()
    # Fit a simple latent-space one-step transition only after AE training.
    z=[];zn=[]
    with torch.no_grad():
        z=model.encode(xa_tensor:=x).numpy();zn=model.encode(torch.tensor(chronological_examples(train_data,keep,mean,std,0)[1],dtype=torch.float32)).numpy()
    transition=Ridge(alpha=1.0).fit(z,zn)
    return model,transition

def predict_dynamic(model,current,dt,history,mode):
    x=torch.tensor(current[None,:],dtype=torch.float32);ctx=torch.tensor([[float(dt)]],dtype=torch.float32) if mode in {"delta_t","memory"} else None;h=None
    if mode=="memory":
        raw=torch.tensor(history,dtype=torch.float32)
        if raw.ndim!=3:raise ValueError(f"Memory history must be [batch,time,features], got {tuple(raw.shape)}")
        b,t,d=raw.shape
        with torch.no_grad():h=model.encode(raw.reshape(b*t,d)).reshape(b,t,-1)
    with torch.no_grad():return model.predict_observation(x,context=ctx,history=h).numpy()[0]

def evaluate_dynamic(model,train_data,heldout_data,keep,mean,std,cfg,mode):
    times,Xraw=heldout_data;X=transform(Xraw,keep,mean,std);n_prefix=min(max(2,cfg.history_len+1,int(np.ceil(len(times)*cfg.prefix_fraction))),len(times)-1)
    prefix_t=list(times[:n_prefix]);observed=[v.copy() for v in X[:n_prefix]];current=observed[-1].copy();current_t=float(prefix_t[-1]);true=[];pred=[];persistence=[];nearest=[];linear=[]
    train_trans={k:transform(v,keep,mean,std) for k,(_,v) in train_data.items()};scale=max(float(times[-1]-times[0]),1.0)
    for j in range(n_prefix,len(times)):
        dt=(float(times[j])-current_t)/scale;hist=np.asarray(observed[-cfg.history_len:])[None,:,:] if mode=="memory" else None;p=predict_dynamic(model,current,dt,hist,mode)
        pred.append(p);persistence.append(observed[-1]);true.append(X[j]);nearest.append(np.mean([train_trans[ds][int(np.argmin(np.abs(tt-times[j])))] for ds,(tt,_) in train_data.items()],axis=0))
        linear.append(observed[-1]+((times[j]-prefix_t[-1])/(prefix_t[-1]-prefix_t[-2]))*(observed[-1]-observed[-2]) if len(prefix_t)>=2 and prefix_t[-1]!=prefix_t[-2] else observed[-1]);current=p;current_t=float(times[j]);observed.append(p)
    return _metrics(true,pred,persistence,nearest,linear)

def evaluate_autoencoder(model,transition,train_data,heldout_data,keep,mean,std,cfg):
    times,Xraw=heldout_data;X=transform(Xraw,keep,mean,std);n_prefix=min(max(2,int(np.ceil(len(times)*cfg.prefix_fraction)),len(times)-1),len(times)-1);Z=model.encode(torch.tensor(X,dtype=torch.float32)).detach().numpy();true=Z[n_prefix:];pred=[];current=Z[n_prefix-1].copy()
    for _ in range(n_prefix,len(times)):current=transition.predict(current[None,:])[0];pred.append(current.copy())
    persistence=np.repeat(Z[n_prefix-1][None,:],len(true),axis=0);return _metrics(true,np.asarray(pred),persistence,np.full_like(true,np.nan),np.full_like(true,np.nan))

def _metrics(true,pred,persistence,nearest,linear):
    rm=lambda a,b:float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)))
    out={"n_future_points":len(true),"rmse_model":rm(true,pred),"rmse_persistence":rm(true,persistence),"rmse_nearest":rm(true,nearest) if np.isfinite(nearest).all() else np.nan,"rmse_linear":rm(true,linear) if np.isfinite(linear).all() else np.nan,"_true":np.asarray(true),"_pred":np.asarray(pred),"_persistence":np.asarray(persistence)}
    out["improvement_vs_persistence"]=out["rmse_persistence"]-out["rmse_model"];out["improvement_vs_nearest"]=out["rmse_nearest"]-out["rmse_model"] if np.isfinite(out["rmse_nearest"]) else np.nan;out["improvement_vs_linear"]=out["rmse_linear"]-out["rmse_model"] if np.isfinite(out["rmse_linear"]) else np.nan;return out

def evaluate_pca(train_data,heldout_data,keep,mean,std,cfg):
    train_z=np.vstack([transform(X,keep,mean,std) for _,X in train_data.values()]);scaler=StandardScaler().fit(train_z);pca=PCA(n_components=min(cfg.state_dim,train_z.shape[0],train_z.shape[1])).fit(scaler.transform(train_z));times,Xraw=heldout_data;Z=pca.transform(scaler.transform(transform(Xraw,keep,mean,std)));n_prefix=min(max(2,int(np.ceil(len(times)*cfg.prefix_fraction)),len(times)-1),len(times)-1);true=Z[n_prefix:];pers=np.repeat(Z[n_prefix-1][None,:],len(true),axis=0);linear=[]
    for j in range(n_prefix,len(times)):linear.append(Z[n_prefix-1]+((times[j]-times[n_prefix-1])/(times[n_prefix-1]-times[n_prefix-2]))*(Z[n_prefix-1]-Z[n_prefix-2]) if times[n_prefix-1]!=times[n_prefix-2] else Z[n_prefix-1])
    return _metrics(true,np.asarray(linear),pers,np.full_like(true,np.nan),np.asarray(linear))

def permutation_p(runs,n_perm,seed=12345):
    rng=np.random.default_rng(seed);obs=float(np.mean([r["improvement_vs_persistence"] for r in runs]));null=np.empty(n_perm)
    for k in range(n_perm):
        vals=[]
        for r in runs:
            y=r["_true"];p=r["_pred"];b=r["_persistence"];perm=y[rng.permutation(len(y))];rm=lambda a,c:np.sqrt(np.mean((a-c)**2));vals.append(float(rm(perm,b)-rm(perm,p)))
        null[k]=np.mean(vals)
    return obs,float((np.sum(null>=obs)+1)/(n_perm+1)),null

def benchmark_fold(train_data,heldout_data,cfg,seed,heldout_name):
    keep,mean,std=training_statistics(train_data,cfg.max_genes);rows=[]
    results={}
    for mode in ("pca","autoencoder","markov","delta_t","memory"):
        if mode=="pca":r=evaluate_pca(train_data,heldout_data,keep,mean,std,cfg)
        elif mode=="autoencoder":
            model,transition=train_autoencoder(train_data,keep,mean,std,cfg,seed);r=evaluate_autoencoder(model,transition,train_data,heldout_data,keep,mean,std,cfg)
        else:
            model,loss=train_dynamic(train_data,keep,mean,std,cfg,seed,mode);r=evaluate_dynamic(model,train_data,heldout_data,keep,mean,std,cfg,mode);r["train_loss"]=loss
        results[mode]=r
        rows.append({"seed":seed,"heldout_dataset":heldout_name,"model":mode,"n_selected_genes":len(keep),**{k:v for k,v in r.items() if not k.startswith("_")}})
    return rows,results

def benchmark(data,cfg,seeds):
    all_rows=[];all_runs={m:[] for m in ("pca","autoencoder","markov","delta_t","memory")}
    for seed in seeds:
        for held in sorted(data):
            train={k:v for k,v in data.items() if k!=held};rows,runs=benchmark_fold(train,data[held],cfg,seed,held);all_rows.extend(rows)
            for m,r in runs.items():all_runs[m].append(r)
    df=pd.DataFrame(all_rows);summary=df.groupby("model",as_index=False).agg(n_runs=("rmse_model","size"),mean_rmse=("rmse_model","mean"),median_rmse=("rmse_model","median"),mean_improvement_vs_persistence=("improvement_vs_persistence","mean"),q05_improvement_vs_persistence=("improvement_vs_persistence",lambda x:x.quantile(.05)),q95_improvement_vs_persistence=("improvement_vs_persistence",lambda x:x.quantile(.95)))
    pvals=[]
    for m in summary.model:
        obs,p,null=permutation_p(all_runs[m],cfg.permutation_n,seed=9000+list(summary.model).index(m));pvals.append({"model":m,"observed_mean_improvement":obs,"permutation_p":p,"null_mean":float(null.mean()),"null_q95":float(np.quantile(null,.95))})
    perm=pd.DataFrame(pvals);summary=summary.merge(perm,on="model",how="left");summary["predictive_support"]=(summary.mean_improvement_vs_persistence>0)&(summary.q05_improvement_vs_persistence>0)&(summary.permutation_p<0.05)
    return df,summary
