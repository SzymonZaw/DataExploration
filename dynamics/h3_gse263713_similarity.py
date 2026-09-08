from __future__ import annotations

import ast, hashlib, json, re
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; DATA=ROOT/'Data'; OUT=ROOT/'results'/'Dynamics'/'stage2_11c_h3_control_redesign'; TARGET=ROOT/'results'/'Dynamics'/'stage2_11c_control_null'; H3=DATA/'GSE263713_raw_counts.tsv.gz'; EXPECTED_SHA='8ce1e16a0039d93449e70aa6e827bbc8510a9b15dfbb78dc920cc2a2de5615a6'; SEED=20260911; N_PERMUTATIONS=10000; SAMPLE_RE=re.compile(r'^([^-]+)-TP-(\d+(?:\.\d+)?)$')
def sha256(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for c in iter(lambda:f.read(1024*1024),b''): h.update(c)
 return h.hexdigest()
def read_counts():
 d=pd.read_csv(H3,sep='\t',compression='gzip',comment='#',low_memory=False); cols=[c for c in d.columns if SAMPLE_RE.fullmatch(str(c))]
 if len(cols)!=78: raise ValueError(f'Expected 78 samples, found {len(cols)}')
 x=d[cols].apply(pd.to_numeric,errors='coerce').fillna(0.0); x.index=d.iloc[:,0].astype(str).str.strip(); return x.groupby(level=0).sum()
def map_hgnc(x):
 ids=pd.Index(x.index.astype(str)); frac=sum(bool(re.fullmatch(r'ENSG\d+(?:\.\d+)?',i.upper())) for i in ids)/max(len(ids),1)
 if frac<.8: x.index=ids.str.upper(); return x.groupby(level=0).mean()
 import mygene
 mg=mygene.MyGeneInfo(); mp={}; clean=sorted({i.split('.',1)[0].upper() for i in ids})
 for s in range(0,len(clean),1000):
  for r in mg.querymany(clean[s:s+1000],scopes='ensembl.gene',fields='symbol',species='human',as_dataframe=False,verbose=False):
   if not r.get('notfound') and r.get('query') and r.get('symbol'): mp[str(r['query']).upper()]=str(r['symbol']).upper()
 syms=[mp.get(i.split('.',1)[0].upper(),'') for i in ids]; keep=[bool(s) for s in syms]; y=x.loc[keep].copy(); y.index=np.asarray(syms)[keep]; return y.groupby(level=0).mean()
def networks():
 import decoupler as dc
 return {'PROGENy':dc.op.progeny(organism='human',top=100),'DoRothEA':dc.op.dorothea(organism='human',levels=['A','B','C'])}
def activity(expr,net):
 import decoupler as dc
 n=net.copy(); n['target']=n['target'].astype(str).str.upper(); common=expr.index.intersection(n['target'])
 if len(common)<50: raise ValueError(f'Too few shared genes: {len(common)}')
 r=dc.mt.ulm(data=expr.loc[common].T,net=n[n.target.isin(common)],tmin=5)
 if isinstance(r,tuple): r=r[0]
 if isinstance(r,pd.DataFrame) and {'source','score'}.issubset(r.columns): return r.pivot(index=r.index,columns='source',values='score')
 return r
def control_traj(a):
 rows=[]
 for sample in a.index:
  m=SAMPLE_RE.fullmatch(str(sample))
  if m:
   for feature,value in a.loc[sample].items(): rows.append((m.group(1),float(m.group(2)),feature,float(value)))
 d=pd.DataFrame(rows,columns=['person','time','feature','value']); return d.groupby(['feature','time'],as_index=False).value.mean().pivot(index='feature',columns='time',values='value')
def parse_list(v):
 if isinstance(v,(list,tuple,np.ndarray)): return list(v)
 if pd.isna(v): return []
 try: return list(ast.literal_eval(str(v)))
 except Exception: return [float(v)]
def target_traj():
 p=TARGET/'01_candidate_trajectories.csv'; d=pd.read_csv(p)
 fc=next((c for c in ('feature','candidate','source') if c in d.columns),None)
 if fc is None: raise ValueError('No feature column in target trajectories')
 if 'a_values' not in d.columns or 'b_values' not in d.columns: raise ValueError('Target trajectories must expose a_values and b_values in the locked historical format')
 rows=[]
 for _,r in d.iterrows():
  av=parse_list(r['a_values']); bv=parse_list(r['b_values'])
  if len(av)!=4 or len(bv)!=4: raise ValueError(f'Unexpected a_values/b_values lengths for {r[fc]}: {len(av)}/{len(bv)}')
  # a_values and b_values are two cohort trajectories for the same feature.
  # They represent the same four time points, so average them pointwise.
  rows.append([r[fc],*np.nanmean(np.asarray([av,bv],dtype=float),axis=0)])
 return pd.DataFrame(rows,columns=[fc,0.,3.,7.,10.]).set_index(fc).apply(pd.to_numeric,errors='coerce')
def stat(a,b):
 common=a.index.intersection(b.index); aa=a.loc[common].to_numpy(float); bb=b.loc[common].to_numpy(float); vals=[]
 for x,y in zip(aa,bb):
  m=np.isfinite(x)&np.isfinite(y)
  if m.sum()>=4 and np.std(x[m])>0 and np.std(y[m])>0: vals.append(np.corrcoef(x[m],y[m])[0,1])
 if not vals: raise ValueError('No valid feature correlations')
 return float(np.median(vals)),len(common)
def null(a,b,rng):
 common=a.index.intersection(b.index); aa=a.loc[common].to_numpy(float); bb=b.loc[common].to_numpy(float); out=np.empty(N_PERMUTATIONS)
 for k in range(N_PERMUTATIONS):
  vals=[]
  for x,y in zip(aa,bb[rng.permutation(len(common))]):
   m=np.isfinite(x)&np.isfinite(y)
   if m.sum()>=4 and np.std(x[m])>0 and np.std(y[m])>0: vals.append(np.corrcoef(x[m],y[m])[0,1])
  out[k]=np.median(vals) if vals else np.nan
 return out
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if sha256(H3)!=EXPECTED_SHA: raise RuntimeError('GSE263713 SHA mismatch')
 gate=OUT/'H3_GSE263713_REPRESENTATION_GATE.json'
 if not gate.exists() or not json.loads(gate.read_text())['decision']['eligible_for_h3_similarity']: raise RuntimeError('Representation gate missing or failed')
 expr=map_hgnc(read_counts()); target=target_traj(); rng=np.random.default_rng(SEED); rows=[]
 for family,net in networks().items():
  ctl=control_traj(activity(expr,net)); obs,n=stat(target,ctl); nd=null(target,ctl,rng); p=(1+int(np.nansum(nd>=obs)))/(N_PERMUTATIONS+1)
  rows.append({'control':'GSE263713','feature_family':family,'n_common_features':n,'observed_median_trajectory_corr':obs,'null_median':float(np.nanmedian(nd)),'null_q95':float(np.nanquantile(nd,.95)),'empirical_one_sided_p':float(p)})
  pd.DataFrame({'null_statistic':nd}).to_csv(OUT/f'H3_GSE263713_{family}_NULL.csv',index=False)
 result={'status':'ok','accession':'GSE263713','seed':SEED,'n_permutations':N_PERMUTATIONS,'target_trajectory_format':'a_values and b_values are two cohort trajectories averaged pointwise at the four Yamanaka time points','temporal_alignment':[0.0,3.0,7.0,10.0],'controls':rows,'decision':{'status':'H3_UNRESOLVED','reason':'One orthogonal temporal control is insufficient for the precommitted final H3 decision.'},'guardrail':'Similarity is a falsification diagnostic, not evidence of causality or mechanism.'}
 (OUT/'H3_GSE263713_SIMILARITY_RESULT.json').write_text(json.dumps(result,indent=2)); pd.DataFrame(rows).to_csv(OUT/'H3_GSE263713_SIMILARITY_RESULT.csv',index=False); print(json.dumps(result,indent=2))
if __name__=='__main__': main()
