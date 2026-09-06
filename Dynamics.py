"""Dynamics v5.8 — biological gene-level harmonization for OSKM dynamics."""
from pathlib import Path
import json
import re
import ssl
import urllib.request
import urllib.parse
from io import StringIO
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/"results"; OUT=RESULTS/"Dynamics"
for n in range(1,10):(OUT/f"stage{n}").mkdir(parents=True,exist_ok=True)
for n in ("stage2_1","stage2_2","stage2_3","stage2_4","stage2_5","stage2_6"):(OUT/n).mkdir(parents=True,exist_ok=True)
STAGE21=OUT/"stage2_1"; STAGE22=OUT/"stage2_2"; STAGE24=OUT/"stage2_4"; STAGE26=OUT/"stage2_6"; CACHE=STAGE26/"cache"; CACHE.mkdir(parents=True,exist_ok=True)
PCA_FILES={"GSE28688":RESULTS/"GSE28688/non_normalized/07_PCA_coordinates.csv","GSE148158":RESULTS/"GSE148158/07_PCA_coordinates.csv","GSE52052":RESULTS/"GSE52052/08_PCA_coordinates.csv","GSE67462":RESULTS/"GSE67462/09_PCA_coordinates.csv","GSE297234":RESULTS/"GSE297234/08_PCA_coordinates.csv"}
FEATURE_FILES={"GSE148158":RESULTS/"GSE148158/expression.csv","GSE28688":RESULTS/"GSE28688/non_normalized/expression_input.csv","GSE52052":RESULTS/"GSE52052/expression_log2.csv","GSE67462":RESULTS/"GSE67462/03_expression_for_EDA.csv","GSE297234":RESULTS/"GSE297234/03_log1p_CPM_sample_expression.csv"}
PLATFORMS={"GSE28688":("GPL6883","human"),"GSE52052":("GPL14550","human")}
GSM_TIME={"GSM4455240":48.,"GSM4455241":48.,"GSM4455242":72.,"GSM4455243":72.,"GSM4455244":48.,"GSM4455245":72.,"GSM710515":24.,"GSM710516":24.,"GSM710517":48.,"GSM710518":48.,"GSM710519":72.,"GSM710520":72.,"GSM1258008":264.,"GSM1258009":264.,"GSM1258010":264.,"GSM1258011":264.,"GSM1258012":264.,"GSM1258013":264.,"GSM1647454":0.,"GSM1647455":0.,"GSM1647456":24.,"GSM1647457":24.,"GSM1647458":72.,"GSM1647459":72.,"GSM1647460":120.,"GSM1647461":120.,"GSM1647462":168.,"GSM1647463":168.,"GSM1647464":264.,"GSM1647465":264.,"GSM1647466":360.,"GSM1647467":360.,"GSM1647468":432.,"GSM1647469":432.,"GSM8986586":0.,"GSM8986587":72.,"GSM8986588":168.,"GSM8986589":240.,"GSM8986590":0.,"GSM8986591":72.,"GSM8986592":168.,"GSM8986593":240.}
GSE28688_ROW_SAMPLE=[f"GSM{x}" for x in range(710513,710527)]; GSE28688_ROW_TIME=[0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

def load_pca(path):
    if not path.exists():return None
    df=pd.read_csv(path,index_col=0); pcs=[c for c in ("PC1","PC2","PC3") if c in df.columns]
    if len(pcs)<3:return None
    x=df[pcs].apply(pd.to_numeric,errors="coerce"); x.index=x.index.astype(str); return x

def time_hours(ds,sample):
    s=str(sample).strip().strip('"')
    if s in GSM_TIME:return GSM_TIME[s]
    t=s.lower().replace("_"," ").replace("-"," ")
    patterns={"GSE28688":[(r"24\s*h",24.),(r"48\s*h",48.),(r"72\s*h",72.)],"GSE148158":[(r"48",48.),(r"72",72.)],"GSE52052":[(r"day\s*11",264.)],"GSE67462":[(r"day\s*0\b",0.),(r"day\s*1\b",24.),(r"day\s*3\b",72.),(r"day\s*5\b",120.),(r"day\s*7\b",168.),(r"day\s*11\b",264.),(r"day\s*15\b",360.),(r"day\s*18\b",432.)],"GSE297234":[(r"d0\b|day\s*0\b",0.),(r"d3\b|day\s*3\b",72.),(r"d7\b|day\s*7\b",168.),(r"d10\b|day\s*10\b",240.)]}
    for p,v in patterns.get(ds,[]):
        if re.search(p,t):return v
    return np.nan

def condition(ds,sample):
    s=str(sample).lower()
    if ds=="GSE148158":
        if "oskm" in s:return "OSKM"
        if "gfp" in s:return "GFP"
        if "h1" in s or "h9" in s:return "hESC"
        if "bj" in s:return "BJ_fibroblast"
    if ds=="GSE297234":return "aged" if any(x in s for x in ("6586","6587","6588","6589")) else "young" if any(x in s for x in ("6590","6591","6592","6593")) else "all"
    return "all"

def replicate(sample):
    s=str(sample); m=re.search(r"(?:-|_|\s)([ab])$",s,re.I)
    if m:return m.group(1).lower()
    m=re.search(r"(?:rep|replicate)[_\s-]*(\d+)",s,re.I)
    if m:return m.group(1)
    m=re.fullmatch(r"GSM(\d+)",s)
    if m and 1647454<=int(m.group(1))<=1647469:return "1" if int(m.group(1))%2==0 else "2"
    return "unknown"

def zscore(x):
    x=pd.Series(x,dtype=float); sd=x.std(ddof=0); return (x-x.mean())/sd if np.isfinite(sd) and sd>0 else pd.Series(np.nan,index=x.index)

def orient(x):
    x=pd.Series(x,dtype=float).copy(); v=x.dropna()
    if len(v) and x.loc[v.abs().idxmax()]<0:x=-x
    return x

def stage1_data_integration():
    rows=[]; states=[]
    for ds,path in PCA_FILES.items():
        x=load_pca(path)
        if x is None:rows.append({"dataset":ds,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","path":str(path)});continue
        o=x.copy(); o.insert(0,"sample",o.index.astype(str)); source="GSM_or_text"
        if ds=="GSE28688" and len(o)==14:o["sample"]=GSE28688_ROW_SAMPLE; source="GSE28688_GEO_row_order"
        o["dataset"]=ds; o["time_hours"]=[time_hours(ds,s) for s in o["sample"]]
        if ds=="GSE28688" and source=="GSE28688_GEO_row_order":o["time_hours"]=GSE28688_ROW_TIME
        o["condition"]=[condition(ds,s) for s in o["sample"]]; o["stage"]=o["time_hours"].map(lambda t:f"day{int(t/24)}" if pd.notna(t) and t%24==0 else f"{int(t)}h" if pd.notna(t) else "unknown"); o["replicate"]= [replicate(s) for s in o["sample"]]
        for i,pc in enumerate(("PC1","PC2","PC3"),1):o[f"latent_{i}"]=zscore(orient(o[pc]))
        timed=o[o.time_hours.notna()]; rows.append({"dataset":ds,"PCA_file_found":True,"n_samples":len(o),"n_timed_samples":len(timed),"n_unique_times":timed.time_hours.nunique(),"role":"trajectory" if timed.time_hours.nunique()>=2 else "context_only","path":str(path)}); states.append(o)
    availability=pd.DataFrame(rows); state=pd.concat(states,ignore_index=True) if states else pd.DataFrame(); availability.to_csv(OUT/"stage1/01_dataset_availability.csv",index=False); state.to_csv(OUT/"stage1/02_master_sample_metadata.csv",index=False); return state,availability

def curve(state,ds,grid,branch=None):
    g=state[(state.dataset==ds)&state.time_hours.notna()]
    if branch is not None:g=g[g.condition==branch]
    m=g.groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index()
    if len(m)<2:return None
    t=m.index.to_numpy(float); u=(t-t.min())/(t.max()-t.min()); return np.column_stack([np.interp(grid,u,m[c]) for c in ("latent_1","latent_2","latent_3")])

def rotation(a,b):
    aa,bb=a-a.mean(0),b-b.mean(0); u,_,vt=np.linalg.svd(aa.T@bb); return u@vt

def fit_transform(source,target):
    rot=rotation(source,target); a=(source-source.mean(0))@rot; scale=np.divide(np.std(target,0),np.where(np.std(a,0)>0,np.std(a,0),1)); return rot,source.mean(0),scale,target.mean(0)

def apply_transform(x,tr):
    rot,sm,scale,tm=tr; y=(x-sm)@rot; y*=scale; return y+tm

def stage2_1(state):
    datasets=[d for d in state.dataset.unique() if state[(state.dataset==d)&state.time_hours.notna()].time_hours.nunique()>=2]; grid=np.linspace(0,1,25); curves={d:curve(state,d,grid) for d in datasets}; curves={d:c for d,c in curves.items() if c is not None and np.isfinite(c).all()}; ref="GSE67462" if "GSE67462" in curves else sorted(curves)[0]; aligned={ref:curves[ref]}
    for d,c in curves.items():
        if d!=ref:aligned[d]=apply_transform(c,fit_transform(c,curves[ref]))
    out=state.copy()
    for j in range(1,4):out[f"common_latent_{j}"]=np.nan
    out["common_latent_status"]="not_aligned"
    for d,a in aligned.items():
        idx=out.index[(out.dataset==d)&out.time_hours.notna()]; lo,hi=out.loc[idx,"time_hours"].min(),out.loc[idx,"time_hours"].max(); u=(out.loc[idx,"time_hours"].to_numpy()-lo)/(hi-lo)
        for j in range(3):out.loc[idx,f"common_latent_{j+1}"]=np.interp(u,grid,a[:,j])
        out.loc[idx,"common_latent_status"]="time_anchored_aligned"
    out.to_csv(STAGE21/"04_sample_common_latent_state.csv",index=False); return out

def stage2_2(state):
    datasets=sorted(state.loc[state.common_latent_status=="time_anchored_aligned","dataset"].unique()); grid=np.linspace(0,1,25); curves={d:curve(state,d,grid) for d in datasets}; rows=[]
    for i,a in enumerate(datasets):
        for b in datasets[i+1:]:
            x,y=curves[a],curves[b]; rows.append({"dataset_a":a,"dataset_b":b,"trajectory_correlation":np.corrcoef(x.ravel(),y.ravel())[0,1],"aligned_rmse":np.sqrt(np.mean((x-y)**2)),"path_length_a":np.linalg.norm(np.diff(x,axis=0),axis=1).sum(),"path_length_b":np.linalg.norm(np.diff(y,axis=0),axis=1).sum()})
    df=pd.DataFrame(rows)
    if not df.empty:df["path_length_ratio"]=df.path_length_a/df.path_length_b
    df.to_csv(STAGE22/"02_cross_dataset_distances.csv",index=False); return df

def branches(state):
    result=[]
    for ds in state.dataset.unique():
        g=state[(state.dataset==ds)&state.time_hours.notna()]
        if g.time_hours.nunique()<2:continue
        valid=[c for c in sorted(g.condition.unique()) if curve(state,ds,np.linspace(0,1,25),None if c=="all" else c) is not None]
        result.extend((ds,c) for c in valid) if ds=="GSE148158" and len(valid)>=2 else result.append((ds,"all"))
    return result

def stage2_4(state):
    rows=[]; grid=np.linspace(0,1,25); ref="GSE67462"
    for ds,branch in branches(state):
        mask=(state.dataset==ds)&state.time_hours.notna()
        if branch!="all":mask&=state.condition==branch
        g=state[mask]; times=sorted(g.time_hours.unique())
        if len(times)<3:continue
        lo,hi=min(times),max(times); ref_curve=curve(state,ref,grid)
        if ref_curve is None:continue
        for held in times:
            train=[t for t in times if t!=held]
            if len(train)<2:continue
            m=g[g.time_hours.isin(train)].groupby("time_hours")[["latent_1","latent_2","latent_3"]].mean().sort_index()
            if len(m)<2:continue
            tu=(m.index.to_numpy(float)-lo)/(hi-lo); target=np.column_stack([np.interp(grid,tu,m[c]) for c in ("latent_1","latent_2","latent_3")]); tr=fit_transform(target,ref_curve); hu=(held-lo)/(hi-lo); observed=g[g.time_hours==held][["latent_1","latent_2","latent_3"]].mean().to_numpy(float); observed_common=apply_transform(observed.reshape(1,-1),tr)[0]; ref_pred=np.array([np.interp(hu,grid,ref_curve[:,j]) for j in range(3)]); naive=float(np.linalg.norm(observed_common-ref_pred)); pred=np.array([np.interp(hu,grid,target[:,j]) for j in range(3)]); err=float(np.linalg.norm(observed_common-pred)); rows.append({"dataset":ds,"trajectory_branch":branch,"held_out_time_hours":held,"n_train_times":len(train),"oos_error":err,"naive_error":naive,"error_reduction_vs_naive":naive-err,"relative_error":err/max(naive,1e-12)})
    df=pd.DataFrame(rows); df.to_csv(STAGE24/"01_leave_one_timepoint_out.csv",index=False); return df

def read_feature_matrix(path):
    df=pd.read_csv(path,index_col=0); df=df.apply(pd.to_numeric,errors="coerce"); df=df.dropna(axis=0,how="all").dropna(axis=1,how="all"); df.index=df.index.astype(str).str.strip().str.strip('"'); return df.groupby(level=0).mean()

def normalize_gene_symbol(x):
    s=str(x).strip().strip('"').upper(); return s if s and s not in {"NAN","NA","NONE","NULL","-","."} else None

def clean_refseq(x):
    s=str(x).strip().strip('"'); s=re.sub(r"_at$","",s,flags=re.I); s=re.sub(r"\.\d+$","",s); return s.upper()

def extract_refseq(x):
    s=str(x).strip().strip('"'); m=re.search(r"((?:NM|NR|XM|XR)_\d+(?:\.\d+)?)",s,re.I); return clean_refseq(m.group(1)) if m else None

def extract_ensembl(x):
    s=str(x).strip().strip('"'); m=re.search(r"(ENS[A-Z]*G\d+)",s,re.I); return m.group(1).upper() if m else None

def read_mapping_cache(name):
    p=CACHE/name
    if not p.exists():return {}
    try:return pd.read_csv(p).set_index(0)[1].dropna().astype(str).to_dict()
    except Exception:return {}

def write_mapping_cache(name,mapping):
    pd.Series(mapping,dtype=str).rename_axis(0).reset_index(name=1).to_csv(CACHE/name,index=False)

def _fetch_url(url,timeout=30):
    req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 DataExploration"}); return urllib.request.urlopen(req,timeout=timeout).read().decode("utf-8","replace")

def biomart_ensembl_to_human(ids):
    ids=list(dict.fromkeys(ids)); cache=read_mapping_cache("ensembl_to_human.csv"); missing=[x for x in ids if x not in cache]
    if not missing:return cache
    for host in ("https://www.ensembl.org","https://useast.ensembl.org"):
        try:
            chunks=[]
            for i in range(0,len(missing),2000):
                vals="\n".join(missing[i:i+2000]); query='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Filter name="ensembl_gene_id" value="'+",".join(missing[i:i+2000])+'"/><Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/></Dataset></Query>'
                q=urllib.parse.quote(query); chunks.append(_fetch_url(host+"/biomart/martservice?query="+q))
            for txt in chunks:
                for line in txt.splitlines():
                    p=line.split("\t");
                    if len(p)>=2 and p[0] and p[1]:cache[p[0].upper()]=p[1].upper()
            if cache:write_mapping_cache("ensembl_to_human.csv",cache); return cache
        except Exception as e: print(f"  BioMart human Ensembl batch 1 failed: {e}")
    try:
        import mygene
        mg=mygene.MyGeneInfo(); res=mg.querymany(missing,scopes="ensembl.gene",fields="symbol",species="human",as_dataframe=False,returnall=False,verbose=False)
        for r in res:
            if r.get("query") and r.get("symbol"):cache[str(r["query"]).upper()]=str(r["symbol"]).upper()
        write_mapping_cache("ensembl_to_human.csv",cache)
    except Exception:pass
    return cache

def biomart_refseq_to_human(ids):
    ids=list(dict.fromkeys(ids)); cache=read_mapping_cache("refseq_to_human.csv"); missing=[x for x in ids if x not in cache]
    if not missing:return cache
    for host in ("https://www.ensembl.org","https://useast.ensembl.org"):
        try:
            for i in range(0,len(missing),2000):
                vals=",".join(missing[i:i+2000]); query='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default"><Filter name="refseq_mrna" value="'+vals+'"/><Attribute name="refseq_mrna"/><Attribute name="external_gene_name"/></Dataset></Query>'; q=urllib.parse.quote(query); txt=_fetch_url(host+"/biomart/martservice?query="+q)
                for line in txt.splitlines():
                    p=line.split("\t");
                    if len(p)>=2 and p[0] and p[1]:cache[clean_refseq(p[0])]=p[1].upper()
            if cache:write_mapping_cache("refseq_to_human.csv",cache); return cache
        except Exception as e: print(f"  BioMart RefSeq batch failed: {e}")
    try:
        import mygene
        mg=mygene.MyGeneInfo(); res=mg.querymany(missing,scopes="refseq.rna",fields="symbol",species="human",as_dataframe=False,returnall=False,verbose=False)
        for r in res:
            if r.get("query") and r.get("symbol"):cache[clean_refseq(r["query"])]=str(r["symbol"]).upper()
        write_mapping_cache("refseq_to_human.csv",cache)
    except Exception:pass
    return cache

def _feature_gene_mapping(ds,df):
    direct={}
    for feat in df.index:
        s=normalize_gene_symbol(feat)
        if s and re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}",s) and not re.match(r"^(ENSG|ENSMUSG|NM_|NR_|XM_|XR_)",s):direct[feat]=s
    ens={feat:extract_ensembl(feat) for feat in df.index}; ens_ids=sorted({x for x in ens.values() if x}); ref={feat:extract_refseq(feat) for feat in df.index}; ref_ids=sorted({x for x in ref.values() if x}); return direct,ens,ens_ids,ref,ref_ids

def merge_direct_and_mapping(direct,df,mapping,extractor):
    rows={}; mapped=0
    for feat in df.index:
        gene=direct.get(feat)
        if gene is None:
            key=extractor(feat); gene=mapping.get(key) if key else None
        if gene:
            mapped+=1; rows.setdefault(gene,[]).append(feat)
    out=pd.DataFrame(index=sorted(rows),columns=df.columns,dtype=float)
    for gene,features in rows.items():out.loc[gene]=df.loc[features].mean(axis=0).to_numpy()
    return out,mapped

def stage2_6():
    matrices={}; audits=[]
    for ds,path in FEATURE_FILES.items():
        if not path.exists():audits.append({"dataset":ds,"feature_file":str(path),"n_features":0,"n_samples":0,"mapped_features":0,"mapped_genes":0,"status":"missing"});continue
        df=read_feature_matrix(path); n_features,n_samples=df.shape; direct,ens,ens_ids,ref,ref_ids=_feature_gene_mapping(ds,df); mapped_direct=len(direct)
        if ds=="GSE67462":
            bm=biomart_refseq_to_human(ref_ids); mat,mapped=merge_direct_and_mapping(direct,df,bm,extract_refseq); matrices[ds]=mat; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped+mapped_direct,"mapped_genes":len(mat),"status":"mapped_mouse_refseq_to_human" if len(mat) else "no_gene_mapping"})
        elif ds in ("GSE148158",):
            if ens_ids:
                bm=biomart_ensembl_to_human(ens_ids); mat,mapped=merge_direct_and_mapping(direct,df,bm,extract_ensembl); matrices[ds]=mat; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped+mapped_direct,"mapped_genes":len(mat),"status":"mapped_human_ensembl_and_symbols" if len(mat) else "no_gene_mapping"})
            else:mat=direct; matrices[ds]=mat; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped_direct,"mapped_genes":len(mat),"status":"direct_gene_ids" if len(mat) else "no_gene_mapping"})
        elif ds=="GSE28688":
            matrices[ds]=direct; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped_direct,"mapped_genes":len(direct),"status":"mapped_to_human_gene_space" if len(direct) else "no_gene_mapping"})
        elif ds=="GSE52052":
            matrices[ds]=direct; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped_direct,"mapped_genes":len(direct),"status":"mapped_to_human_gene_space" if len(direct) else "no_gene_mapping"})
        else:
            if ens_ids:
                bm=biomart_ensembl_to_human(ens_ids); mat,mapped=merge_direct_and_mapping(direct,df,bm,extract_ensembl); matrices[ds]=mat; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped+mapped_direct,"mapped_genes":len(mat),"status":"mapped_human_ensembl_and_symbols" if len(mat) else "no_gene_mapping"})
            elif ref_ids:
                bm=biomart_refseq_to_human(ref_ids); mat,mapped=merge_direct_and_mapping(direct,df,bm,extract_refseq); matrices[ds]=mat; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped+mapped_direct,"mapped_genes":len(mat),"status":"mapped_human_refseq_and_symbols" if len(mat) else "no_gene_mapping"})
            else:
                matrices[ds]=direct; audits.append({"dataset":ds,"feature_file":str(path),"n_features":n_features,"n_samples":n_samples,"mapped_features":mapped_direct,"mapped_genes":len(direct),"status":"direct_gene_ids" if len(direct) else "no_gene_mapping"})
    audit=pd.DataFrame(audits); audit.to_csv(STAGE26/"01_gene_mapping_audit.csv",index=False); print(audit.to_string(index=False))
    contributing=[d for d,m in matrices.items() if len(m)>0]; common=set.intersection(*(set(m.index) for m in matrices.values() if len(m)>0)) if matrices else set(); print(f"\ncommon human genes across all contributing datasets = {len(common)}"); print(f"datasets contributing = {len(contributing)}")
    if common:
        genes=sorted(common); blocks=[]; meta=[]
        for ds in contributing:
            m=matrices[ds].loc[genes].copy(); original_samples=m.columns.astype(str).tolist(); z=m.apply(lambda x:zscore(x).to_numpy(),axis=1,result_type="expand"); z.index=genes; z.columns=[f"{ds}__{s}" for s in original_samples]; blocks.append(z); meta.extend({"dataset":ds,"sample":s,"matrix_column":f"{ds}__{s}"} for s in original_samples)
        pd.concat(blocks,axis=1).to_csv(STAGE26/"06_common_human_gene_matrix.csv"); pd.DataFrame(meta).to_csv(STAGE26/"07_common_gene_sample_metadata.csv",index=False)
    status="sufficient_common_human_gene_space" if len(common)>=1000 and len(contributing)>=4 else "insufficient_common_human_gene_space"
    pd.DataFrame([{"common_human_genes":len(common),"datasets_contributing":len(contributing),"status":status,"time_used_for_feature_space":False}]).to_csv(STAGE26/"08_stage26_decision.csv",index=False)
    print(f"status = {status}"); print("NOTE: time is not used to construct the feature space. Human datasets are harmonized to human gene symbols; GSE67462 uses mouse RefSeq → human ortholog mapping with BioMart and MyGene fallback."); print("="*88); return common

def main():
    state,availability=stage1_data_integration(); aligned=stage2_1(state); stage2_2(aligned); oos=stage2_4(state); common=stage2_6(); print(f"\nDynamics v5.8 results written to: {OUT}"); print(f"Datasets with PCA: {availability.PCA_file_found.sum() if not availability.empty else 0}/5"); print(f"Stage 2.1 aligned trajectory datasets: {aligned.loc[aligned.common_latent_status=='time_anchored_aligned','dataset'].nunique() if not aligned.empty else 0}"); print(f"Stage 2.4 tested OOS rows: {len(oos)}"); print(f"Stage 2.6 common human genes: {len(common)}")

if __name__=="__main__":main()
