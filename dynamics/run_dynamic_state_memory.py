"""Memory-aware DynamicStateModel experiment.

Compares a Markov latent dynamics model z(t+dt)=F(z(t),dt) against a
history-aware model z(t+dt)=F(z(t),h(t),dt), using strict LODO prefix-to-future
forecasting, multiple seeds and a time-permutation null. The purpose is to test
whether current latent state is sufficient or whether observed history adds
predictive information.
"""
from __future__ import annotations
import argparse, random
from pathlib import Path
import numpy as np, pandas as pd, torch
from torch import nn
from .dynamic_state_model import DynamicStateModel
from .validation import _load_common_space
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results"/"Dynamics"/"dynamic_state_memory"
DATASETS=["GSE67462","GSE28688","GSE297234"]; DEFAULT_SEEDS=(311,312,313,314,315); HISTORY_LEN=2

def set_seed(seed): random.seed(seed); np.random.seed(seed); torch.manual_seed(seed)

def load_data():
    matrix,meta=_load_common_space(); data={}
    for ds in DATASETS:
        g=meta[(meta.dataset==ds)&meta.time_hours.notna()&meta.matrix_column.notna()].copy().sort_values("time_hours")
        if g.time_hours.nunique()<3: continue
        X=matrix.loc[:,g.matrix_column.astype(str).tolist()].T.copy(); X.index=g.time_hours.to_numpy(float); X=X.groupby(level=0,sort=True).mean(); data[ds]=(X.index.to_numpy(float),X.to_numpy(float))
    if len(data)<2: raise RuntimeError("At least two longitudinal datasets are required.")
    return data

def stats(train,max_genes):
    arrays=[X for _,X in train.values()]; pooled=np.vstack(arrays)
    med=np.nanmedian(np.where(np.isfinite(pooled),pooled,np.nan),axis=0); med=np.where(np.isfinite(med),med,0.)
    pooled=np.vstack([np.where(np.isfinite(X),X,med) for X in arrays]); var=np.var(pooled,axis=0)
    keep=np.argsort(var)[::-1][:min(max_genes,len(var))]; mean=pooled[:,keep].mean(0); std=pooled[:,keep].std(0); std=np.where(std>1e-8,std,1.)
    return keep,mean,std

def transform(X,keep,mean,std):
    X=X[:,keep]; X=np.where(np.isfinite(X),X,mean); return (X-mean)/std

def sequence_examples(train,keep,mean,std):
    xs=[]; hs=[]; ys=[]; dts=[]
    for times,Xraw in train.values():
        X=transform(Xraw,keep,mean,std); scale=max(float(times[-1]-times[0]),1.); tn=(times-times[0])/scale
        for i in range(HISTORY_LEN,len(times)-1):
            xs.append(X[i]); hs.append(X[i-HISTORY_LEN:i]); ys.append(X[i+1]); dts.append(tn[i+1]-tn[i])
    if not xs: raise RuntimeError("Training trajectories are too short for the requested history length.")
    return (torch.tensor(np.asarray(xs),dtype=torch.float32),torch.tensor(np.asarray(hs),dtype=torch.float32),torch.tensor(np.asarray(ys),dtype=torch.float32),torch.tensor(np.asarray(dts)[:,None],dtype=torch.float32))

def latent_history(model,h):
    b,k,d=h.shape; return model.encode(h.reshape(b*k,d)).reshape(b,k,model.state_dim)

def train_model(train_data,keep,mean,std,state_dim,hidden_dim,epochs,lr,seed,use_history):
    set_seed(seed); x,h_obs,y,dt=sequence_examples(train_data,keep,mean,std)
    model=DynamicStateModel(input_dim=x.shape[1],state_dim=state_dim,hidden_dim=hidden_dim,context_dim=1,history_dim=(16 if use_history else 0),dropout=.05)
    opt=torch.optim.AdamW(model.parameters(),lr=lr,weight_decay=1e-4); best=float("inf"); patience=max(25,epochs//8); stale=0
    for _ in range(epochs):
        model.train(); opt.zero_grad(set_to_none=True); z=model.encode(x); zn=model.encode(y).detach(); rec=model.decode(z)
        h_lat=latent_history(model,h_obs) if use_history else None; pred=model.transition(z,context=dt,history=h_lat); pred_obs=model.decode(pred)
        loss=((rec-x)**2).mean()+.25*((pred-zn)**2).mean()+((pred_obs-y)**2).mean(); loss.backward(); nn.utils.clip_grad_norm_(model.parameters(),1.); opt.step(); val=float(loss.detach())
        if val<best-1e-6: best=val; stale=0
        else: stale+=1
        if stale>=patience: break
    return model,best

def predict(model,current,history,dt,use_history):
    x=torch.tensor(current[None,:],dtype=torch.float32); d=torch.tensor([[float(dt)]],dtype=torch.float32); h=None
    if use_history:
        ho=torch.tensor(history[None,:,:],dtype=torch.float32); h=latent_history(model,ho)
    model.eval()
    with torch.no_grad(): return model.predict_observation(x,context=d,history=h).numpy()[0]

def metrics(y,p):
    e=np.asarray(y)-np.asarray(p); return float(np.sqrt(np.mean(e*e))),float(np.mean(np.abs(e)))

def evaluate(model,train_data,held,keep,mean,std,prefix_fraction,use_history):
    times,Xraw=held; X=transform(Xraw,keep,mean,std); scale=max(float(times[-1]-times[0]),1.); n=max(HISTORY_LEN+1,int(np.ceil(len(times)*prefix_fraction))); n=min(n,len(times)-1)
    prefix=list(X[:n]); true=[]; pm=[]; pp=[]
    for j in range(n,len(times)):
        hist=np.asarray(prefix[-HISTORY_LEN:]); dt=(times[j]-times[j-1])/scale
        pm.append(predict(model,prefix[-1],hist,dt,use_history)); pp.append(prefix[-1]); true.append(X[j]); prefix.append(X[j])
    if not true: return None
    rm,ma=metrics(true,pm); rp,_=metrics(true,pp)
    return {"n_future_points":len(true),"rmse":rm,"mae":ma,"rmse_persistence":rp,"improvement_vs_persistence":rp-rm}

def permutation_null(data,max_genes,state_dim,hidden_dim,epochs,lr,seed,prefix_fraction,n_perm):
    rng=np.random.default_rng(seed); vals=[]
    for k in range(n_perm):
        perm={ds:(t,X[rng.permutation(len(t))]) for ds,(t,X) in data.items()}; fold=[]
        for held_name in sorted(perm):
            train={q:v for q,v in perm.items() if q!=held_name}
            if len(train)<2: continue
            keep,mean,std=stats(train,max_genes); m,_=train_model(train,keep,mean,std,state_dim,hidden_dim,max(20,epochs//3),lr,seed+1000+k,True); r=evaluate(m,train,perm[held_name],keep,mean,std,prefix_fraction,True)
            if r: fold.append(r["improvement_vs_persistence"])
        if fold: vals.append(float(np.mean(fold)))
        if (k+1)%10==0: print(f"DynamicStateModel memory permutation {k+1}/{n_perm}",flush=True)
    return np.asarray(vals)

def run(max_genes=2000,state_dim=8,hidden_dim=128,epochs=250,lr=1e-3,seeds=DEFAULT_SEEDS,prefix_fraction=.6,n_perm=100):
    OUT.mkdir(parents=True,exist_ok=True); data=load_data(); rows=[]
    for seed in seeds:
        for held_name in sorted(data):
            train={k:v for k,v in data.items() if k!=held_name}
            if len(train)<2: continue
            keep,mean,std=stats(train,max_genes)
            for kind,use_history in (("markov",False),("memory",True)):
                model,loss=train_model(train,keep,mean,std,state_dim,hidden_dim,epochs,lr,seed,use_history); r=evaluate(model,train,data[held_name],keep,mean,std,prefix_fraction,use_history)
                rows.append({"seed":seed,"heldout_dataset":held_name,"model":kind,"n_selected_genes":len(keep),"train_loss":loss,**r})
                print(f"DynamicStateModel memory: seed={seed} heldout={held_name} model={kind} RMSE={r['rmse']:.4f} persistence={r['rmse_persistence']:.4f}",flush=True)
    df=pd.DataFrame(rows); df.to_csv(OUT/"01_lodo_memory_metrics.csv",index=False)
    piv=df.pivot_table(index=["seed","heldout_dataset"],columns="model",values="rmse").reset_index(); piv["memory_gain_vs_markov"]=piv.markov-piv.memory; piv.to_csv(OUT/"04_memory_vs_markov.csv",index=False)
    m=df[df.model=="memory"]; summary=pd.DataFrame([{"n_valid_fold_models":len(df),"n_seeds":df.seed.nunique(),"mean_rmse_markov":df[df.model=="markov"].rmse.mean(),"mean_rmse_memory":m.rmse.mean(),"mean_memory_gain_vs_markov":piv.memory_gain_vs_markov.mean(),"median_memory_gain_vs_markov":piv.memory_gain_vs_markov.median(),"q05_memory_gain_vs_markov":piv.memory_gain_vs_markov.quantile(.05),"q95_memory_gain_vs_markov":piv.memory_gain_vs_markov.quantile(.95),"mean_improvement_memory_vs_persistence":m.improvement_vs_persistence.mean()}])
    null=permutation_null(data,max_genes,state_dim,hidden_dim,epochs,lr,seeds[0],prefix_fraction,n_perm); obs=float(summary.loc[0,"mean_memory_gain_vs_markov"]); p=float((1+np.sum(null>=obs))/(1+len(null))) if len(null) else np.nan
    summary["permutation_p_memory_gain"]=p; summary["memory_predictive_support"]=bool(summary.loc[0,"mean_memory_gain_vs_markov"]>0 and p<.05); summary["stage3_readiness"]=False
    summary.to_csv(OUT/"02_summary.csv",index=False); pd.DataFrame({"null_memory_gain":null}).to_csv(OUT/"03_permutation_null.csv",index=False); print(summary.to_string(index=False),flush=True); return summary

if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--max-genes",type=int,default=2000); p.add_argument("--state-dim",type=int,default=8); p.add_argument("--hidden-dim",type=int,default=128); p.add_argument("--epochs",type=int,default=250); p.add_argument("--lr",type=float,default=1e-3); p.add_argument("--prefix-fraction",type=float,default=.6); p.add_argument("--permutations",type=int,default=100); p.add_argument("--seeds",type=int,nargs="+",default=list(DEFAULT_SEEDS)); a=p.parse_args(); run(a.max_genes,a.state_dim,a.hidden_dim,a.epochs,a.lr,tuple(a.seeds),a.prefix_fraction,a.permutations)
