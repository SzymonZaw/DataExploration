"""Dynamics v5.4 — biological gene-level harmonization for OSKM dynamics."""
from pathlib import Path
import json
import re
import ssl
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
OUT = RESULTS / "Dynamics"
for n in range(1, 10):
    (OUT / f"stage{n}").mkdir(parents=True, exist_ok=True)
for name in ("stage2_1", "stage2_2", "stage2_3", "stage2_4", "stage2_5", "stage2_6"):
    (OUT / name).mkdir(parents=True, exist_ok=True)
STAGE21 = OUT / "stage2_1"; STAGE22 = OUT / "stage2_2"; STAGE24 = OUT / "stage2_4"; STAGE26 = OUT / "stage2_6"
CACHE = STAGE26 / "cache"; CACHE.mkdir(parents=True, exist_ok=True)

PCA_FILES = {
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "07_PCA_coordinates.csv",
    "GSE148158": RESULTS / "GSE148158" / "07_PCA_coordinates.csv",
    "GSE52052": RESULTS / "GSE52052" / "08_PCA_coordinates.csv",
    "GSE67462": RESULTS / "GSE67462" / "09_PCA_coordinates.csv",
    "GSE297234": RESULTS / "GSE297234" / "08_PCA_coordinates.csv",
}
FEATURE_FILES = {
    "GSE148158": RESULTS / "GSE148158" / "expression.csv",
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "expression_input.csv",
    "GSE52052": RESULTS / "GSE52052" / "expression_log2.csv",
    "GSE67462": RESULTS / "GSE67462" / "03_expression_for_EDA.csv",
    "GSE297234": RESULTS / "GSE297234" / "03_log1p_CPM_sample_expression.csv",
}
PLATFORMS = {"GSE28688": ("GPL6883", "human"), "GSE52052": ("GPL14550", "human"), "GSE67462": ("GPL19972", "mouse")}
GSM_TIME = {
    "GSM4455240":48.,"GSM4455241":48.,"GSM4455242":72.,"GSM4455243":72.,"GSM4455244":48.,"GSM4455245":72.,
    "GSM710515":24.,"GSM710516":24.,"GSM710517":48.,"GSM710518":48.,"GSM710519":72.,"GSM710520":72.,
    "GSM1258008":264.,"GSM1258009":264.,"GSM1258010":264.,"GSM1258011":264.,"GSM1258012":264.,"GSM1258013":264.,
    "GSM1647454":0.,"GSM1647455":0.,"GSM1647456":24.,"GSM1647457":24.,"GSM1647458":72.,"GSM1647459":72.,"GSM1647460":120.,"GSM1647461":120.,"GSM1647462":168.,"GSM1647463":168.,"GSM1647464":264.,"GSM1647465":264.,"GSM1647466":360.,"GSM1647467":360.,"GSM1647468":432.,"GSM1647469":432.,
    "GSM8986586":0.,"GSM8986587":72.,"GSM8986588":168.,"GSM8986589":240.,"GSM8986590":0.,"GSM8986591":72.,"GSM8986592":168.,"GSM8986593":240.,
}
GSE28688_ROW_SAMPLE=[f"GSM{x}" for x in range(710513,710527)]
GSE28688_ROW_TIME=[0.,0.,24.,24.,48.,48.,72.,72.,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]

def load_pca(path):
    if not path.exists(): return None
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
    if ds=="GSE297234":
        if any(x in s for x in ("6586","6587","6588","6589")):return "aged"
        if any(x in s for x in ("6590","6591","6592","6593")):return "young"
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
    if len(v):
        i=v.abs().idxmax()
        if x.loc[i]<0:x=-x
    return x

def stage1_data_integration():
    rows=[]; states=[]
    for ds,path in PCA_FILES.items():
        x=load_pca(path)
        if x is None:
            rows.append({"dataset":ds,"PCA_file_found":False,"n_samples":0,"n_timed_samples":0,"n_unique_times":0,"role":"unavailable","path":str(path)}); continue
        o=x.copy(); o.insert(0,"sample",o.index.astype(str)); source="GSM_or_text"
        if ds=="GSE28688" and len(o)==14:o["sample"]=GSE28688_ROW_SAMPLE; source="GSE28688_GEO_row_order"
        o["dataset"]=ds; o["time_hours"]=[time_hours(ds,s) for s in o["sample"]]
        if ds=="GSE28688" and source=="GSE28688_GEO_row_order":o["time_hours"]=GSE28688_ROW_TIME
        o["condition"]=[condition(ds,s) for s in o["sample"]]; o["stage"]=o.time_hours.map(lambda t:f"day{int(t/24)}" if pd.notna(t) and t%24==0 else f"{int(t)}h" if pd.notna(t) else "unknown"); o["replicate"]=[replicate(s) for s in o.sample]
        for i,pc in enumerate(("PC1","PC2","PC3"),1):o[f"latent_{i}"]=zscore(orient(o[pc]))
        timed=o[o.time_hours.notna()]; rows.append({"dataset":ds,"PCA_file_found":True,"n_samples":len(o),"n_timed_samples":len(timed),"n_unique_times":timed.time_hours.nunique(),"role":"trajectory" if timed.time_hours.nunique()>=2 else "context_only","path":str(path)}); states.append(o)
    availability=pd.DataFrame(rows); state=pd.concat(states,ignore_index=True) if states else pd.DataFrame(); availability.to_csv(OUT/"stage1"/"01_dataset_availability.csv",index=False); state.to_csv(OUT/"stage1"/"02_master_sample_metadata.csv",index=False); return state,availability

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
        valid=[]
        for c in sorted(g.condition.unique()):
            if curve(state,ds,np.linspace(0,1,25),None if c=="all" else c) is not None:valid.append(c)
        if ds=="GSE148158" and len(valid)>=2:result.extend((ds,c) for c in valid)
        else:result.append((ds,"all"))
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
            tu=(m.index.to_numpy(float)-lo)/(hi-lo); target=np.column_stack([np.interp(grid,tu,m[c]) for c in ("latent_1","latent_2","latent_3")]); tr=fit_transform(target,ref_curve); hu=(held-lo)/(hi-lo); src_pred=np.array([np.interp(hu,grid,target[:,j]) for j in range(3)]); pred_common=apply_transform(src_pred.reshape(1,-1),tr)[0]; observed=g[g.time_hours==held][["latent_1","latent_2","latent_3"]].mean().to_numpy(float); observed_common=apply_transform(observed.reshape(1,-1),tr)[0]; ref_pred=np.array([np.interp(hu,grid,ref_curve[:,j]) for j in range(3)]); err=float(np.linalg.norm(pred_common-observed_common)); naive=float(np.linalg.norm(observed_common-ref_pred)); rows.append({"dataset":ds,"trajectory_branch":branch,"held_out_time_hours":held,"n_train_times":len(train),"oos_error":err,"naive_error":naive,"error_reduction_vs_naive":naive-err,"relative_error":err/(naive+1e-12)})
    df=pd.DataFrame(rows); df.to_csv(STAGE24/"01_leave_one_timepoint_out.csv",index=False)
    if not df.empty:
        df.groupby(["dataset","trajectory_branch"]).agg(n_tests=("oos_error","size"),mean_oos_error=("oos_error","mean"),median_oos_error=("oos_error","median"),p95_oos_error=("oos_error",lambda x:x.quantile(.95)),mean_relative_error=("relative_error","mean"),mean_error_reduction_vs_naive=("error_reduction_vs_naive","mean")).reset_index().to_csv(STAGE24/"02_oos_summary_by_trajectory.csv",index=False); df.to_csv(STAGE24/"03_oos_vs_naive_baseline.csv",index=False)
    return df

def download_text(url,target):
    if target.exists() and target.stat().st_size>0:return target
    target.parent.mkdir(parents=True,exist_ok=True); part=target.with_suffix(target.suffix+".part"); req=urllib.request.Request(url,headers={"User-Agent":"DataExploration-Dynamics/5.4"}); errors=[]
    for label,context in (("verified",None),("unverified",ssl._create_unverified_context())):
        try:
            if part.exists():part.unlink()
            kwargs={"timeout":120}
            if context is not None:kwargs["context"]=context
            with urllib.request.urlopen(req,**kwargs) as r,open(part,"wb") as f:
                while True:
                    chunk=r.read(1024*1024)
                    if not chunk:break
                    f.write(chunk)
            if part.stat().st_size==0:raise RuntimeError("empty response")
            part.replace(target); print(f"  downloaded {url} ({label} TLS)"); return target
        except Exception as exc:errors.append(f"{label}: {exc}")
    if part.exists():part.unlink()
    raise RuntimeError("; ".join(errors))

def read_geo_platform(gpl):
    cache=CACHE/f"{gpl}_full.txt"; url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&view=full&form=text"
    try:download_text(url,cache)
    except Exception as exc:print(f"  platform download failed for {gpl}: {exc}");return None
    return cache.read_text(encoding="utf-8",errors="replace")

def normalize_gene_symbol(x):
    if x is None:return None
    s=str(x).strip().strip('"').upper()
    if not s or s in {"NA","N/A","NAN","NULL","-","NONE"}:return None
    return re.sub(r"\s+","",s)

def parse_platform_table(text):
    if not text:return pd.DataFrame()
    lines=text.splitlines(); begin=end=None
    for i,line in enumerate(lines):
        if line.startswith("!platform_table_begin"):begin=i+1
        elif line.startswith("!platform_table_end"):end=i;break
    if begin is None:return pd.DataFrame()
    try:return pd.read_csv(pd.io.common.StringIO("\n".join(lines[begin:end])),sep="\t",dtype=str,comment="#")
    except Exception:return pd.DataFrame()

def pick_column(df,candidates):
    lower={str(c).strip().lower():c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:return lower[c.lower()]
    for c in df.columns:
        if any(k in str(c).lower() for k in candidates):return c
    return None

def fetch_platform_mapping(gpl):
    table=parse_platform_table(read_geo_platform(gpl))
    if table.empty:return pd.DataFrame(columns=["feature_id","gene_symbol"])
    id_col=pick_column(table,["id","id_ref","probeid","probe_id","probename","array_address_id"]); gene_col=pick_column(table,["gene symbol","gene_symbol","genesymbol","gene_assignment","symbol","gene"])
    if id_col is None or gene_col is None:
        print(f"  no usable ID/gene columns found for {gpl}; columns={list(table.columns)[:12]}"); return pd.DataFrame(columns=["feature_id","gene_symbol"])
    out=pd.DataFrame({"feature_id":table[id_col].astype(str),"gene_symbol":table[gene_col].map(normalize_gene_symbol)}).dropna(subset=["gene_symbol"]); return out[out.feature_id.str.strip().ne("")].drop_duplicates()

def read_feature_matrix(path):
    if not path.exists():return None
    try:df=pd.read_csv(path,index_col=0,low_memory=False)
    except Exception as exc:print(f"  cannot read {path}: {exc}");return None
    x=df.apply(pd.to_numeric,errors="coerce"); x=x.loc[x.notna().sum(axis=1)>0]; x=x.loc[:,x.notna().sum(axis=0)>0]; return x.groupby(level=0).mean()

def map_matrix_to_genes(expr,mapping,dataset):
    if expr is None:return None,pd.DataFrame()
    m=mapping.copy(); m.feature_id=m.feature_id.astype(str).str.strip(); expr=expr.copy(); expr.index=expr.index.astype(str).str.strip(); merged=expr.reset_index(names="feature_id").merge(m,on="feature_id",how="left"); audit=merged[["feature_id","gene_symbol"]].copy(); audit["dataset"]=dataset; audit["mapped"]=audit.gene_symbol.notna(); mapped=merged.dropna(subset=["gene_symbol"]); sample_cols=[c for c in expr.columns if c in mapped.columns]; return mapped.groupby("gene_symbol")[sample_cols].mean(),audit

def direct_gene_matrix(expr,dataset):
    if expr is None:return None,pd.DataFrame()
    ids=pd.Series(expr.index.astype(str),index=expr.index); mapped=ids.map(normalize_gene_symbol); valid=mapped.notna()&mapped.str.match(r"^[A-Z0-9][A-Z0-9._-]*$",na=False); audit=pd.DataFrame({"feature_id":ids.values,"gene_symbol":mapped.values,"dataset":dataset,"mapped":valid.values}); gene=expr.loc[valid].copy(); gene.index=mapped[valid].values; return gene.groupby(level=0).mean(),audit

def http_text(url,data=None,content_type="text/plain",timeout=180):
    req=urllib.request.Request(url,data=data,headers={"Content-Type":content_type,"User-Agent":"DataExploration-Dynamics/5.4"},method="POST" if data is not None else "GET")
    errors=[]
    for label,context in (("verified",None),("unverified",ssl._create_unverified_context())):
        try:
            kwargs={"timeout":timeout}
            if context is not None:kwargs["context"]=context
            with urllib.request.urlopen(req,**kwargs) as r:return r.read().decode("utf-8",errors="replace")
        except Exception as exc:errors.append(f"{label}: {exc}")
    raise RuntimeError("; ".join(errors))

def biomart_refseq_to_human(accessions):
    """Bulk-map mouse RefSeq mRNA accessions directly to human symbols via Ensembl BioMart."""
    accessions=sorted({str(x).strip() for x in accessions if str(x).strip()})
    cache=CACHE/"GPL19972_refseq_to_human.tsv"
    if cache.exists() and cache.stat().st_size>0:
        try:return pd.read_csv(cache,sep="\t",dtype=str)
        except Exception:pass
    rows=[]; endpoint="https://www.ensembl.org/biomart/martservice"
    for start in range(0,len(accessions),500):
        batch=accessions[start:start+500]
        xml='''<?xml version="1.0" encoding="UTF-8"?>\n<!DOCTYPE Query>\n<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="0" datasetConfigVersion="0.6">\n<Dataset name="mmusculus_gene_ensembl" interface="default">\n<Filter name="refseq_mrna" value="%s"/>\n<Attribute name="refseq_mrna"/>\n<Attribute name="external_gene_name"/>\n<Attribute name="hsapiens_homolog_associated_gene_name"/>\n</Dataset>\n</Query>''' % ','.join(batch)
        try:
            text=http_text(endpoint,data=xml.encode("utf-8"),content_type="application/xml",timeout=180)
            if text.startswith("Query ERROR") or "ERROR" in text[:500]:
                raise RuntimeError(text[:500])
            from io import StringIO
            tab=pd.read_csv(StringIO(text),sep="\t",dtype=str)
            if not tab.empty:
                tab.columns=["refseq_mrna","mouse_gene","human_gene"][:len(tab.columns)]; rows.append(tab)
            print(f"  BioMart RefSeq mapping: {min(start+500,len(accessions))}/{len(accessions)}")
        except Exception as exc:
            print(f"  BioMart batch {start//500+1} failed: {exc}")
    if rows:
        out=pd.concat(rows,ignore_index=True).fillna(""); out["mouse_gene"]=out.mouse_gene.map(normalize_gene_symbol); out["human_gene"]=out.human_gene.map(normalize_gene_symbol); out=out[(out.refseq_mrna!="")&(out.human_gene.notna())].drop_duplicates(); out.to_csv(cache,sep="\t",index=False); return out
    return pd.DataFrame(columns=["refseq_mrna","mouse_gene","human_gene"])

def stage2_6(state):
    audit_rows=[]; matrices={}
    for ds,path in FEATURE_FILES.items():
        expr=read_feature_matrix(path)
        if expr is None:
            audit_rows.append({"dataset":ds,"feature_file":str(path),"n_features":0,"n_samples":0,"mapped_features":0,"mapped_genes":0,"status":"feature_matrix_unavailable"});continue
        if ds in PLATFORMS:
            gpl,species=PLATFORMS[ds]
            if ds=="GSE67462":
                table=parse_platform_table(read_geo_platform(gpl)); id_col=pick_column(table,["id"]); ref_col=pick_column(table,["gb_acc","gb acc","refseq"])
                if id_col is None or ref_col is None:
                    print(f"  GPL19972 requires ID + GB_ACC; columns={list(table.columns)[:12]}"); gene=pd.DataFrame(); audit=pd.DataFrame({"feature_id":expr.index.astype(str),"gene_symbol":np.nan,"dataset":ds,"mapped":False})
                else:
                    platform=pd.DataFrame({"feature_id":table[id_col].astype(str).str.strip(),"refseq_mrna":table[ref_col].astype(str).str.strip()}); platform=platform[platform.refseq_mrna.ne("")&platform.refseq_mrna.ne("nan")].drop_duplicates()
                    merged=expr.reset_index(names="feature_id").merge(platform,on="feature_id",how="left"); bm=biomart_refseq_to_human(merged.refseq_mrna.dropna().astype(str).tolist()); merged=merged.merge(bm[["refseq_mrna","human_gene"]],on="refseq_mrna",how="left"); audit=merged[["feature_id","human_gene"]].rename(columns={"human_gene":"gene_symbol"}); audit["dataset"]=ds; audit["mapped"]=audit.gene_symbol.notna(); mapped=merged.dropna(subset=["human_gene"]); sample_cols=[c for c in expr.columns if c in mapped.columns]; gene=mapped.groupby("human_gene")[sample_cols].mean()
                audit.to_csv(STAGE26/f"{ds}_feature_to_gene.csv",index=False); matrices[ds]=gene; audit_rows.append({"dataset":ds,"feature_file":str(path),"n_features":expr.shape[0],"n_samples":expr.shape[1],"mapped_features":int(audit.mapped.sum()),"mapped_genes":gene.shape[0],"status":"mapped_to_human_gene_space" if not gene.empty else "no_gene_mapping"}); continue
            mapping=fetch_platform_mapping(gpl); gene,audit=map_matrix_to_genes(expr,mapping,ds); audit.to_csv(STAGE26/f"{ds}_feature_to_gene.csv",index=False); matrices[ds]=gene; audit_rows.append({"dataset":ds,"feature_file":str(path),"n_features":expr.shape[0],"n_samples":expr.shape[1],"mapped_features":int(audit.mapped.sum()),"mapped_genes":gene.shape[0],"status":"mapped_to_human_gene_space" if not gene.empty else "no_gene_mapping"})
        else:
            gene,audit=direct_gene_matrix(expr,ds); audit.to_csv(STAGE26/f"{ds}_feature_to_gene.csv",index=False); matrices[ds]=gene; audit_rows.append({"dataset":ds,"feature_file":str(path),"n_features":expr.shape[0],"n_samples":expr.shape[1],"mapped_features":int(audit.mapped.sum()),"mapped_genes":gene.shape[0],"status":"direct_gene_ids" if not gene.empty else "no_gene_mapping"})
    audit_df=pd.DataFrame(audit_rows); audit_df.to_csv(STAGE26/"01_gene_mapping_audit.csv",index=False); nonempty=[x for x in matrices.values() if x is not None and not x.empty]; genes=sorted(set.intersection(*(set(x.index) for x in nonempty))) if len(nonempty)==len(matrices) and matrices else []
    overlap_rows=[]; dslist=sorted(matrices)
    for i,a in enumerate(dslist):
        for b in dslist[i+1:]:overlap_rows.append({"dataset_a":a,"dataset_b":b,"common_human_genes":len(set(matrices[a].index)&set(matrices[b].index))})
    pd.DataFrame(overlap_rows).to_csv(STAGE26/"05_pairwise_human_gene_overlap.csv",index=False)
    if genes:
        parts=[]
        for ds in dslist:
            x=matrices[ds].loc[genes].copy(); x=x.sub(x.mean(axis=1),axis=0).div(x.std(axis=1,ddof=0).replace(0,np.nan),axis=0).fillna(0); x.columns=[f"{ds}__{c}" for c in x.columns]; parts.append(x)
        common=pd.concat(parts,axis=1); common.to_csv(STAGE26/"06_common_human_gene_matrix.csv"); meta=pd.DataFrame({"sample":common.columns}); meta["dataset"]=meta.sample.str.split("__",n=1).str[0]; meta["original_sample"]=meta.sample.str.split("__",n=1).str[1]; meta.to_csv(STAGE26/"07_common_gene_sample_metadata.csv",index=False)
    status="success" if genes else "insufficient_common_human_gene_space"; pd.DataFrame([{"common_human_genes":len(genes),"n_datasets":len(matrices),"status":status,"time_used_for_feature_construction":False}]).to_csv(STAGE26/"08_stage26_decision.csv",index=False)
    print("\n"+"="*88);print("STAGE 2.6 — TIME-INDEPENDENT BIOLOGICAL GENE-LEVEL HARMONIZATION");print("="*88);print(audit_df.to_string(index=False));print(f"\ncommon human genes across all contributing datasets = {len(genes)}");print(f"datasets contributing = {len(matrices)}");print(f"status = {status}");print("NOTE: time is not used to construct the feature space. GPL19972 RefSeq accessions are bulk-mapped through Ensembl BioMart directly to human ortholog gene symbols.");print("="*88+"\n");return matrices,genes

def main():
    state,availability=stage1_data_integration();aligned=stage2_1(state);stage2_2(aligned);oos=stage2_4(aligned);matrices,genes=stage2_6(aligned);print(f"Dynamics v5.4 results written to: {OUT}");print(f"Datasets with PCA: {availability.PCA_file_found.sum()}/{len(availability)}");print(f"Stage 2.1 aligned trajectory datasets: {aligned.common_latent_status.eq('time_anchored_aligned').groupby(aligned.dataset).any().sum()}");print(f"Stage 2.4 tested OOS rows: {len(oos)}");print(f"Stage 2.6 common human genes: {len(genes)}")

if __name__=="__main__":main()
