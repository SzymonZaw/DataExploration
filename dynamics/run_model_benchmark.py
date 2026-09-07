"""Run the unified DynamicStateModel comparison framework."""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
from .validation import _load_common_space
from .model_benchmark import BenchmarkConfig, benchmark

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results"/"Dynamics"/"dynamic_state_benchmark"
DATASETS=("GSE67462","GSE28688","GSE297234")
SEEDS=(411,412,413,414,415)

def load_data():
    matrix,meta=_load_common_space();data={}
    for ds in DATASETS:
        g=meta[(meta.dataset==ds)&meta.time_hours.notna()&meta.matrix_column.notna()].copy().sort_values("time_hours")
        if g.time_hours.nunique()<3:continue
        X=matrix.loc[:,g.matrix_column.astype(str).tolist()].T.copy();X.index=g.time_hours.to_numpy(float);X=X.groupby(level=0,sort=True).mean();data[ds]=(X.index.to_numpy(float),X.to_numpy(float))
    if len(data)<3:raise RuntimeError(f"Expected 3 trajectory datasets, got {len(data)}")
    return data

def run(max_genes=2000,state_dim=8,hidden_dim=128,epochs=250,lr=1e-3,prefix_fraction=.6,seeds=SEEDS,permutation_n=1000):
    OUT.mkdir(parents=True,exist_ok=True);data=load_data();cfg=BenchmarkConfig(max_genes=max_genes,state_dim=state_dim,hidden_dim=hidden_dim,epochs=epochs,lr=lr,prefix_fraction=prefix_fraction,permutation_n=permutation_n)
    df,summary=benchmark(data,cfg,seeds);df.to_csv(OUT/"01_model_metrics.csv",index=False);summary.to_csv(OUT/"02_model_summary.csv",index=False)
    pivot=df.pivot_table(index=["seed","heldout_dataset"],columns="model",values="rmse_model");pivot.to_csv(OUT/"03_fold_comparison.csv")
    print("\nUnified DynamicStateModel benchmark:");print(summary.to_string(index=False),flush=True);return summary

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--max-genes",type=int,default=2000);p.add_argument("--state-dim",type=int,default=8);p.add_argument("--hidden-dim",type=int,default=128);p.add_argument("--epochs",type=int,default=250);p.add_argument("--lr",type=float,default=1e-3);p.add_argument("--prefix-fraction",type=float,default=.6);p.add_argument("--permutation-n",type=int,default=1000);p.add_argument("--seeds",type=int,nargs="+",default=list(SEEDS));a=p.parse_args();run(a.max_genes,a.state_dim,a.hidden_dim,a.epochs,a.lr,a.prefix_fraction,tuple(a.seeds),a.permutation_n)
