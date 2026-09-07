"""Unified leakage-free predictive benchmark for cellular-state dynamics.

All model scores are computed in the same transformed observation space. The
held-out trajectory is used only as an observed prefix plus future ground truth;
model rollouts never consume future held-out observations.
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

MODELS = ("pca", "autoencoder", "markov", "delta_t", "memory")

@dataclass
class BenchmarkConfig:
    max_genes: int = 2000
    state_dim: int = 8
    hidden_dim: int = 128
    epochs: int = 250
    lr: float = 1e-3
    prefix_fraction: float = .6
    history_len: int = 2
    history_dim: int = 16
    dropout: float = .05
    permutation_n: int = 1000
    time_scale_hours: float = 168.0

class StaticAutoencoder(nn.Module):
    def __init__(self, input_dim, state_dim, hidden_dim, dropout):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(input_dim, hidden_dim), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden_dim, state_dim))
        self.decoder = nn.Sequential(nn.Linear(state_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, input_dim))
    def encode(self, x): return self.encoder(x)
    def decode(self, z): return self.decoder(z)
    def forward(self, x):
        z = self.encode(x); return z, self.decode(z)

def set_seed(seed):
    random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)

def training_statistics(train, max_genes):
    arrays = [np.asarray(X, dtype=float) for _, X in train.values()]
    raw = np.vstack(arrays)
    med = np.nanmedian(np.where(np.isfinite(raw), raw, np.nan), axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    pooled = np.vstack([np.where(np.isfinite(X), X, med) for X in arrays])
    var = np.var(pooled, axis=0)
    keep = np.argsort(var)[::-1][:min(max_genes, pooled.shape[1])]
    selected = pooled[:, keep]
    mean = selected.mean(0)
    std = np.where(selected.std(0) > 1e-8, selected.std(0), 1.0)
    return keep, mean, std

def transform(X, keep, mean, std):
    X = np.asarray(X, dtype=float)[:, keep]
    return (np.where(np.isfinite(X), X, mean) - mean) / std

def chronological_examples(data, keep, mean, std, history_len=0, time_scale_hours=168.0):
    xs, ys, dts, hist = [], [], [], []
    for times, Xraw in data.values():
        order = np.argsort(np.asarray(times, dtype=float))
        times = np.asarray(times, dtype=float)[order]
        X = transform(Xraw, keep, mean, std)[order]
        for i in range(history_len, len(times) - 1):
            xs.append(X[i]); ys.append(X[i + 1])
            dts.append([(times[i + 1] - times[i]) / max(time_scale_hours, 1e-8)])
            if history_len: hist.append(X[i-history_len:i])
    return np.asarray(xs), np.asarray(ys), np.asarray(dts), np.asarray(hist) if history_len else None

def train_dynamic(train_data, keep, mean, std, cfg, seed, mode):
    set_seed(seed)
    hl = cfg.history_len if mode == "memory" else 0
    xa, ya, dta, ha = chronological_examples(train_data, keep, mean, std, hl, cfg.time_scale_hours)
    if len(xa) == 0: raise RuntimeError(f"No training transitions available for mode={mode}")
    x = torch.tensor(xa, dtype=torch.float32); y = torch.tensor(ya, dtype=torch.float32); dt = torch.tensor(dta, dtype=torch.float32)
    context_dim = 1 if mode in {"delta_t", "memory"} else 0
    history_dim = cfg.history_dim if mode == "memory" else 0
    model = DynamicStateModel(x.shape[1], cfg.state_dim, cfg.hidden_dim, context_dim, history_dim, cfg.dropout)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg.lr, weight_decay=1e-4)
    best = float("inf"); stale = 0; patience = max(25, cfg.epochs // 8)
    for _ in range(cfg.epochs):
        model.train(); opt.zero_grad(set_to_none=True)
        z = model.encode(x); zn = model.encode(y).detach(); rec = model.decode(z)
        if history_dim:
            b, t, d = ha.shape; hraw = torch.tensor(ha, dtype=torch.float32)
            h = model.encode(hraw.reshape(b*t, d)).reshape(b, t, cfg.state_dim)
            pred = model.transition(z, context=dt, history=h)
        else: pred = model.transition(z, context=dt if context_dim else None)
        loss = ((rec-x)**2).mean() + .25*((pred-zn)**2).mean() + ((model.decode(pred)-y)**2).mean()
        loss.backward(); nn.utils.clip_grad_norm_(model.parameters(), 1.0); opt.step(); value=float(loss.detach())
        if value < best - 1e-6: best, stale = value, 0
        else: stale += 1
        if stale >= patience: break
    return model, best

def train_autoencoder(train_data, keep, mean, std, cfg, seed):
    set_seed(seed)
    all_x = np.vstack([transform(X, keep, mean, std) for _, X in train_data.values()])
    xa, ya, _, _ = chronological_examples(train_data, keep, mean, std, 0, cfg.time_scale_hours)
    x_all = torch.tensor(all_x, dtype=torch.float32); x_pairs = torch.tensor(xa, dtype=torch.float32); y_pairs = torch.tensor(ya, dtype=torch.float32)
    model = StaticAutoencoder(x_all.shape[1], cfg.state_dim, cfg.hidden_dim, cfg.dropout)
    opt = torch.optim.AdamW(model.parameters(), lr=cfg.lr, weight_decay=1e-4)
    best = float("inf"); stale = 0; patience = max(25, cfg.epochs // 8)
    for _ in range(cfg.epochs):
        opt.zero_grad(set_to_none=True); _, rec = model(x_all); loss=((rec-x_all)**2).mean(); loss.backward(); nn.utils.clip_grad_norm_(model.parameters(),1.0); opt.step(); value=float(loss.detach())
        if value < best - 1e-6: best, stale = value, 0
        else: stale += 1
        if stale >= patience: break
    with torch.no_grad(): z=model.encode(x_pairs).numpy(); zn=model.encode(y_pairs).numpy()
    return model, Ridge(alpha=1.0).fit(z, zn), best

def fit_pca_transition(train_data, keep, mean, std, cfg):
    X = np.vstack([transform(v, keep, mean, std) for _, v in train_data.values()])
    scaler = StandardScaler().fit(X); pca = PCA(n_components=min(cfg.state_dim, X.shape[0], X.shape[1])).fit(scaler.transform(X))
    xa, ya, _, _ = chronological_examples(train_data, keep, mean, std, 0, cfg.time_scale_hours)
    transition = Ridge(alpha=1.0).fit(pca.transform(scaler.transform(xa)), pca.transform(scaler.transform(ya)))
    return scaler, pca, transition

def _prepare_holdout(heldout_data, keep, mean, std, cfg, order=None):
    times, Xraw = heldout_data; times=np.asarray(times,dtype=float); order=np.argsort(times) if order is None else np.asarray(order)
    times=times[order]; X=transform(np.asarray(Xraw,dtype=float),keep,mean,std)[order]
    n_prefix=min(max(2,cfg.history_len+1,int(np.ceil(len(times)*cfg.prefix_fraction))),len(times)-1)
    return times,X,n_prefix

def _dynamic_history_latent(model, history_observations):
    if history_observations is None or len(history_observations)==0: return None
    raw=torch.tensor(np.asarray(history_observations,dtype=np.float32)[None,:,:]); b,t,d=raw.shape
    with torch.no_grad(): return model.encode(raw.reshape(b*t,d)).reshape(b,t,-1)

def predict_dynamic(model,current,dt,history_observations,mode):
    x=torch.tensor(np.asarray(current,dtype=np.float32)[None,:]); ctx=torch.tensor([[float(dt)]],dtype=torch.float32) if mode in {"delta_t","memory"} else None
    h=_dynamic_history_latent(model,history_observations) if mode=="memory" else None
    with torch.no_grad():
        z_next=model.transition(model.encode(x),context=ctx,history=h)
        return model.decode(z_next).numpy()[0], z_next.numpy()[0]

def _metrics(true,pred,persistence,nearest,linear):
    rm=lambda a,b:float(np.sqrt(np.mean((np.asarray(a)-np.asarray(b))**2)))
    true=np.asarray(true); pred=np.asarray(pred); persistence=np.asarray(persistence); nearest=np.asarray(nearest); linear=np.asarray(linear)
    out={"n_future_points":len(true),"rmse_model":rm(true,pred),"rmse_persistence":rm(true,persistence),"rmse_nearest":rm(true,nearest) if np.isfinite(nearest).all() else np.nan,"rmse_linear":rm(true,linear) if np.isfinite(linear).all() else np.nan,"_true":true,"_pred":pred,"_persistence":persistence}
    out["improvement_vs_persistence"]=out["rmse_persistence"]-out["rmse_model"]
    out["improvement_vs_nearest"]=out["rmse_nearest"]-out["rmse_model"] if np.isfinite(out["rmse_nearest"]) else np.nan
    out["improvement_vs_linear"]=out["rmse_linear"]-out["rmse_model"] if np.isfinite(out["rmse_linear"]) else np.nan
    return out

def evaluate_dynamic(model,train_data,heldout_data,keep,mean,std,cfg,mode,order=None):
    times,X,n_prefix=_prepare_holdout(heldout_data,keep,mean,std,cfg,order); scale=max(cfg.time_scale_hours,1e-8)
    current=X[n_prefix-1].copy(); current_t=float(times[n_prefix-1]); history=list(X[max(0,n_prefix-cfg.history_len):n_prefix]) if mode=="memory" else []
    true=[]; pred=[]; persistence=[]; nearest=[]; linear=[]
    train_trans={k:(np.asarray(t,dtype=float),transform(v,keep,mean,std)) for k,(t,v) in train_data.items()}
    persist_value=X[n_prefix-1].copy()
    if n_prefix>=2 and times[n_prefix-1]!=times[n_prefix-2]: slope=(X[n_prefix-1]-X[n_prefix-2])/(times[n_prefix-1]-times[n_prefix-2])
    else: slope=np.zeros_like(persist_value)
    for j in range(n_prefix,len(times)):
        p,_=predict_dynamic(model,current,(float(times[j])-current_t)/scale,history if mode=="memory" else None,mode)
        pred.append(p); true.append(X[j]); persistence.append(persist_value)
        nearest.append(np.mean([xx[int(np.argmin(np.abs(tt-times[j])))] for tt,xx in train_trans.values()],axis=0))
        linear.append(X[n_prefix-1]+slope*(times[j]-times[n_prefix-1]))
        current=p; current_t=float(times[j])
        if mode=="memory": history=(history+[p])[-cfg.history_len:]
    return _metrics(true,pred,persistence,np.asarray(nearest),np.asarray(linear))

def evaluate_autoencoder(model,transition,heldout_data,keep,mean,std,cfg,order=None):
    times,X,n_prefix=_prepare_holdout(heldout_data,keep,mean,std,cfg,order); current_z=model.encode(torch.tensor(X[n_prefix-1][None,:],dtype=torch.float32)).detach().numpy()[0]
    pred=[]; true=[]
    for j in range(n_prefix,len(times)):
        current_z=transition.predict(current_z[None,:])[0]
        pred.append(model.decode(torch.tensor(current_z[None,:],dtype=torch.float32)).detach().numpy()[0]); true.append(X[j])
    persistence=np.repeat(X[n_prefix-1][None,:],len(true),axis=0)
    return _metrics(true,pred,persistence,np.full_like(np.asarray(true),np.nan),np.full_like(np.asarray(true),np.nan))

def evaluate_pca(fitted,heldout_data,keep,mean,std,cfg,order=None):
    scaler,pca,transition=fitted; times,X,n_prefix=_prepare_holdout(heldout_data,keep,mean,std,cfg,order); current_z=pca.transform(scaler.transform(X[n_prefix-1][None,:]))[0]
    pred=[]; true=[]
    for j in range(n_prefix,len(times)):
        current_z=transition.predict(current_z[None,:])[0]
        pred.append(scaler.inverse_transform(pca.inverse_transform(current_z[None,:]))[0]); true.append(X[j])
    persistence=np.repeat(X[n_prefix-1][None,:],len(true),axis=0)
    return _metrics(true,pred,persistence,np.full_like(np.asarray(true),np.nan),np.full_like(np.asarray(true),np.nan))

def benchmark_fold(train_data,heldout_data,cfg,seed,heldout_name):
    keep,mean,std=training_statistics(train_data,cfg.max_genes); rows=[]; runs={}
    for mode in MODELS:
        if mode=="pca": fitted=fit_pca_transition(train_data,keep,mean,std,cfg); r=evaluate_pca(fitted,heldout_data,keep,mean,std,cfg); train_loss=np.nan
        elif mode=="autoencoder":
            model,transition,train_loss=train_autoencoder(train_data,keep,mean,std,cfg,seed); fitted=(model,transition); r=evaluate_autoencoder(model,transition,heldout_data,keep,mean,std,cfg)
        else:
            model,train_loss=train_dynamic(train_data,keep,mean,std,cfg,seed,mode); fitted=model; r=evaluate_dynamic(model,train_data,heldout_data,keep,mean,std,cfg,mode)
        runs[mode]={"model":mode,"seed":seed,"heldout_dataset":heldout_name,"train_data":train_data,"heldout_data":heldout_data,"keep":keep,"mean":mean,"std":std,"cfg":cfg,"fitted":fitted,"observed":r}
        rows.append({"seed":seed,"heldout_dataset":heldout_name,"model":mode,"n_selected_genes":len(keep),**{k:v for k,v in r.items() if not k.startswith("_")},"train_loss":train_loss})
    return rows,runs

def permutation_null(runs,n_perm,seed=12345):
    rng=np.random.default_rng(seed); obs=float(np.mean([r["observed"]["improvement_vs_persistence"] for r in runs])); null=np.empty(n_perm)
    for k in range(n_perm):
        vals=[]
        for r in runs:
            times,X=r["heldout_data"]; shuffled=np.asarray(X)[rng.permutation(len(X))]; pseudo=(np.asarray(times,dtype=float),shuffled); mode=r["model"]
            if mode=="pca": q=evaluate_pca(r["fitted"],pseudo,r["keep"],r["mean"],r["std"],r["cfg"])
            elif mode=="autoencoder": q=evaluate_autoencoder(r["fitted"][0],r["fitted"][1],pseudo,r["keep"],r["mean"],r["std"],r["cfg"])
            else: q=evaluate_dynamic(r["fitted"],r["train_data"],pseudo,r["keep"],r["mean"],r["std"],r["cfg"],mode)
            vals.append(q["improvement_vs_persistence"])
        null[k]=np.mean(vals)
    p=float((np.sum(null>=obs)+1)/(n_perm+1)); return obs,p,null

def benchmark(data,cfg,seeds):
    all_rows=[]; all_runs={m:[] for m in MODELS}
    for seed in seeds:
        for held in sorted(data):
            train={k:v for k,v in data.items() if k!=held}; rows,runs=benchmark_fold(train,data[held],cfg,seed,held); all_rows.extend(rows)
            for m,r in runs.items(): all_runs[m].append(r)
    df=pd.DataFrame(all_rows)
    summary=df.groupby("model",as_index=False).agg(
        n_runs=("rmse_model","size"),
        mean_rmse=("rmse_model","mean"),
        median_rmse=("rmse_model","median"),
        mean_improvement_vs_persistence=("improvement_vs_persistence","mean"),
        q05_improvement_vs_persistence=("improvement_vs_persistence",lambda x:x.quantile(.05)),
        q95_improvement_vs_persistence=("improvement_vs_persistence",lambda x:x.quantile(.95)),
        mean_improvement_vs_nearest=("improvement_vs_nearest","mean"),
        q05_improvement_vs_nearest=("improvement_vs_nearest",lambda x:x.quantile(.05)),
        q95_improvement_vs_nearest=("improvement_vs_nearest",lambda x:x.quantile(.95)),
    )
    pvals=[]
    for idx,m in enumerate(summary.model):
        obs,p,null=permutation_null(all_runs[m],cfg.permutation_n,seed=9000+idx)
        pvals.append({"model":m,"observed_mean_improvement":obs,"permutation_p":p,"null_mean":float(null.mean()),"null_q95":float(np.quantile(null,.95))})
    summary=summary.merge(pd.DataFrame(pvals),on="model",how="left")
    summary["predictive_support"]=(summary.mean_improvement_vs_persistence>0)&(summary.q05_improvement_vs_persistence>0)&(summary.mean_improvement_vs_nearest>0)&(summary.q05_improvement_vs_nearest>0)&(summary.permutation_p<.05)
    return df,summary
