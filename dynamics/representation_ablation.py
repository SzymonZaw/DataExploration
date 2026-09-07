"""Phase 1: ablate the biological state representation.

The forecasting benchmark is kept fixed while the input representation changes:
raw common genes, PROGENy pathway activity and DoRothEA TF activity.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .validation import _load_common_space
from .model_benchmark import BenchmarkConfig, benchmark
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'Dynamics' / 'phase1_representation_ablation'
DATASETS=('GSE67462','GSE28688','GSE297234')
SEEDS=(511,512,513,514,515)

def _trajectory_data(matrix,metadata):
 data={}
 for ds in DATASETS:
  g=metadata[(metadata.dataset==ds)&metadata.time_hours.notna()&metadata.matrix_column.notna()].copy().sort_values('time_hours')
  if g.time_hours.nunique()<3: continue
  X=matrix.loc[:,g.matrix_column.astype(str).tolist()].T.copy(); X.index=g.time_hours.to_numpy(float); X=X.groupby(level=0,sort=True).mean()
  data[ds]=(X.index.to_numpy(float),X.to_numpy(float),list(X.columns))
 return data

def _score_network(data,net):
 import decoupler as dc
 scored={}
 for ds,(times,X,genes) in data.items():
  samples_by_gene=pd.DataFrame(X,columns=genes)
  samples_by_gene.index=[f'{ds}__{i}' for i in range(len(samples_by_gene))]
  arr=samples_by_gene.to_numpy(dtype=float)
  arr[~np.isfinite(arr)]=0.0
  samples_by_gene.iloc[:,:]=arr
  acts,_=dc.mt.ulm(data=samples_by_gene,net=net)
  acts=acts.apply(pd.to_numeric,errors='coerce').replace([np.inf,-np.inf],np.nan).fillna(0.0)
  scored[ds]=(times,acts.to_numpy(float),list(acts.columns))
 return scored

# remainder unchanged intentionally minimal
from pathlib import Path as _P
exec((ROOT/'dynamics'/'representation_ablation.py').read_text())