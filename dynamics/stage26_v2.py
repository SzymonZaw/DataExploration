"""Robust Stage 2.6 gene harmonisation without giant BioMart requests.

The refresh is deliberately conservative: mapping is performed in bounded HTTP batches,
failures are isolated per batch, and validated common-space files are replaced only when
the new intersection satisfies the minimum quality gate.
"""
from pathlib import Path
from io import StringIO
import json
import re
import time
import urllib.parse
import urllib.request
import ssl
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
STAGE26 = RESULTS / "Dynamics" / "stage2_6"
CACHE = STAGE26 / "cache"
STAGE26.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

FEATURE_FILES = {
    "GSE148158": RESULTS / "GSE148158" / "expression.csv",
    "GSE28688": RESULTS / "GSE28688" / "non_normalized" / "expression_input.csv",
    "GSE52052": RESULTS / "GSE52052" / "expression_log2.csv",
    "GSE67462": RESULTS / "GSE67462" / "03_expression_for_EDA.csv",
    "GSE297234": RESULTS / "GSE297234" / "03_log1p_CPM_sample_expression.csv",
}
PLATFORMS = {"GSE28688": "GPL6883", "GSE52052": "GPL14550"}


def _request(url, data=None, timeout=90):
    headers = {"User-Agent": "DataExploration-Dynamics/6.1", "Accept": "application/json,text/plain"}
    body = None if data is None else (data if isinstance(data, bytes) else data.encode())
    req = urllib.request.Request(url, data=body, headers=headers, method="POST" if body is not None else "GET")
    last = None
    for context in (None, ssl._create_unverified_context()):
        try:
            kw = {"timeout": timeout}
            if context is not None: kw["context"] = context
            with urllib.request.urlopen(req, **kw) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:
            last = exc
    raise RuntimeError(str(last))


def normalize_symbol(x):
    s = str(x).strip().strip('"').upper()
    return None if s in {"", "NAN", "NA", "NONE", "NULL", "-", "."} else s


def extract_ensembl(x):
    m = re.search(r"(ENS[A-Z]*G\d+)", str(x), re.I)
    return m.group(1).upper() if m else None


def clean_refseq(x):
    s = str(x).strip().strip('"')
    s = re.sub(r"_at$", "", s, flags=re.I)
    return re.sub(r"\.\d+$", "", s).upper()


def extract_refseq(x):
    m = re.search(r"((?:NM|NR|XM|XR)_\d+(?:\.\d+)?)", str(x), re.I)
    return clean_refseq(m.group(1)) if m else None


def read_feature_matrix(path):
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    df.index = df.index.astype(str).str.strip().str.strip('"')
    return df.groupby(level=0).mean()


def _cache_read(name):
    p = CACHE / name
    if not p.exists(): return {}
    try:
        d = pd.read_csv(p, dtype=str)
        return dict(zip(d.iloc[:,0], d.iloc[:,1]))
    except Exception:
        return {}


def _cache_write(name, mapping):
    pd.DataFrame({"query": list(mapping), "symbol": [mapping[x] for x in mapping]}).to_csv(CACHE / name, index=False)


def mygene_map(ids, scopes, species, fields, cache_name):
    ids = sorted(set(str(x) for x in ids if str(x)))
    out = _cache_read(cache_name)
    missing = [x for x in ids if x not in out]
    for start in range(0, len(missing), 200):
        batch = missing[start:start+200]
        form = urllib.parse.urlencode({"q": ",".join(batch), "scopes": scopes, "fields": fields, "species": species, "size": len(batch)})
        ok = False
        for attempt in range(3):
            try:
                obj = json.loads(_request("https://mygene.info/v3/query", data=form))
                hits = obj if isinstance(obj, list) else obj.get("hits", obj.get("out", []))
                for hit in hits:
                    q = str(hit.get("query", ""))
                    sym = normalize_symbol(hit.get("symbol"))
                    if q and sym: out[q] = sym
                ok = True
                break
            except Exception as exc:
                if attempt == 2: print(f"  MyGene {cache_name} batch {start//200+1} failed: {exc}")
                else: time.sleep(1.0 + attempt)
        if ok and start + 200 < len(missing): time.sleep(0.05)
    _cache_write(cache_name, out)
    return out


def mouse_refseq_human_map(ids):
    ids = sorted(set(clean_refseq(x) for x in ids))
    out = _cache_read("mouse_refseq_human.csv")
    missing = [x for x in ids if x not in out]
    # MyGene supplies HomoloGene/orthology information for mouse RefSeq.
    for start in range(0, len(missing), 200):
        batch = missing[start:start+200]
        form = urllib.parse.urlencode({"q": ",".join(batch), "scopes": "refseq.rna", "fields": "homologene", "species": "mouse", "size": len(batch)})
        try:
            obj = json.loads(_request("https://mygene.info/v3/query", data=form))
            hits = obj if isinstance(obj, list) else obj.get("hits", [])
            for hit in hits:
                q = clean_refseq(hit.get("query", "")); hom = hit.get("homologene", {})
                genes = hom.get("genes", []) if isinstance(hom, dict) else []
                human = [str(p[1]) for p in genes if isinstance(p,(list,tuple)) and len(p)>=2 and str(p[0])=="9606"]
                if human: out[q] = "ENTREZ:" + human[0]
        except Exception as exc:
            print(f"  MyGene mouse orthology batch {start//200+1} failed: {exc}")
        time.sleep(0.05)
    entrez = sorted({v.split(":",1)[1] for v in out.values() if str(v).startswith("ENTREZ:")})
    symbols = mygene_map(entrez, "entrezgene", "human", "symbol", "human_entrez_symbol.csv")
    for q, v in list(out.items()):
        if str(v).startswith("ENTREZ:"):
            s = symbols.get(v.split(":",1)[1])
            if s: out[q] = s
    _cache_write("mouse_refseq_human.csv", out)
    return out


def platform_mapping(gpl):
    p = CACHE / f"{gpl}_platform_mapping.tsv"
    if p.exists(): return pd.read_csv(p, sep="\t", dtype=str).fillna("")
    text = _request(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&view=full&form=text", timeout=120)
    m = re.search(r"!platform_table_begin\n(.*?)\n!platform_table_end", text, re.S)
    if not m: raise RuntimeError(f"No platform table for {gpl}")
    tab = pd.read_csv(StringIO(m.group(1)), sep="\t", dtype=str)
    cols = {str(c).lower().strip(): c for c in tab.columns}
    idcol = next((cols[k] for k in ("id","probe id","probeid") if k in cols), tab.columns[0])
    symcol = next((c for k,c in cols.items() if "gene symbol" in k or k in {"symbol","gene_symbol"}), None)
    if symcol is None: raise RuntimeError(f"No gene-symbol column for {gpl}")
    out = pd.DataFrame({"feature":tab[idcol].astype(str),"human_gene":tab[symcol].map(normalize_symbol)})
    out = out.dropna(subset=["human_gene"]).drop_duplicates("feature")
    out.to_csv(p, sep="\t", index=False)
    return out


def map_dataset(ds, path):
    raw = read_feature_matrix(path)
    acc = {}
    mapped = 0
    if ds in PLATFORMS:
        ann = platform_mapping(PLATFORMS[ds]).set_index("feature")["human_gene"].to_dict()
        for feat, row in raw.iterrows():
            gene = normalize_symbol(ann.get(str(feat)))
            if gene: acc.setdefault(gene, []).append(row.to_numpy(float)); mapped += 1
        source = "GEO_platform"
    elif ds == "GSE67462":
        refs = [extract_refseq(x) for x in raw.index]; refs = [x for x in refs if x]
        mp = mouse_refseq_human_map(refs)
        for feat,row in raw.iterrows():
            gene = mp.get(extract_refseq(feat))
            if gene: acc.setdefault(gene, []).append(row.to_numpy(float)); mapped += 1
        source = "mouse_RefSeq_to_human_ortholog_MyGene"
    else:
        ens = [extract_ensembl(x) for x in raw.index]; ens = [x for x in ens if x]
        mp = mygene_map(ens, "ensembl.gene", "human", "symbol", "human_ensembl_symbol.csv") if ens else {}
        for feat,row in raw.iterrows():
            gene = mp.get(extract_ensembl(feat))
            if gene is None and extract_ensembl(feat) is None: gene = normalize_symbol(feat)
            if gene: acc.setdefault(gene, []).append(row.to_numpy(float)); mapped += 1
        source = "Ensembl_MyGene+direct_symbol"
    mat = pd.DataFrame({g:np.mean(v,axis=0) for g,v in acc.items()}, index=raw.columns).T if acc else pd.DataFrame(index=[],columns=raw.columns)
    return raw,mat,mapped,len(raw)-mapped,len(acc),source


def _atomic(df,path,index=True):
    tmp = path.with_suffix(path.suffix+".tmp")
    df.to_csv(tmp,index=index)
    tmp.replace(path)


def stage2_6_robust(min_genes=1000,min_datasets=4):
    audits=[]; matrices={}
    for ds,path in FEATURE_FILES.items():
        try:
            raw,mat,mapped,unmapped,genes,source=map_dataset(ds,path)
            matrices[ds]=mat
            audits.append({"dataset":ds,"n_features":len(raw),"n_samples":raw.shape[1],"mapped_features":mapped,"unmapped_features":unmapped,"mapped_genes":genes,"mapping_coverage":mapped/max(len(raw),1),"mapping_source":source,"status":"mapped" if mapped else "no_gene_mapping"})
        except Exception as exc:
            matrices[ds]=pd.DataFrame()
            audits.append({"dataset":ds,"n_features":0,"n_samples":0,"mapped_features":0,"unmapped_features":0,"mapped_genes":0,"mapping_coverage":0.0,"mapping_source":"failed","status":f"failed:{type(exc).__name__}"})
            print(f"Stage 2.6 {ds} failed: {type(exc).__name__}: {exc}")
    audit=pd.DataFrame(audits); _atomic(audit,STAGE26/"05_mapping_audit.csv",False)
    valid={d:m for d,m in matrices.items() if isinstance(m,pd.DataFrame) and not m.empty and len(m.index)>0}
    common=set.intersection(*(set(m.index) for m in valid.values())) if valid else set()
    print(f"Stage 2.6: common human genes={len(common)}, contributing datasets={len(valid)}")
    if len(common)<min_genes or len(valid)<min_datasets:
        print("WARNING: insufficient new common space; previous validated files preserved.")
        return {"status":"insufficient_common_human_gene_space","common_genes":len(common),"contributing_datasets":len(valid)}
    genes=sorted(common); blocks=[]; meta=[]
    for ds,m in valid.items():
        b=m.loc[genes].copy()
        b=b.apply(lambda r:(r-r.mean())/r.std(ddof=0) if r.std(ddof=0)>0 else np.nan,axis=1)
        cols=[f"{ds}__{i}" for i in range(b.shape[1])]; b.columns=cols; blocks.append(b)
        for i,c in enumerate(cols): meta.append({"dataset":ds,"sample":str(m.columns[i]),"matrix_column":c})
    _atomic(pd.concat(blocks,axis=1),STAGE26/"06_common_human_gene_matrix.csv",True)
    _atomic(pd.DataFrame(meta),STAGE26/"07_common_gene_sample_metadata.csv",False)
    decision=pd.DataFrame([{"common_human_genes":len(genes),"datasets_contributing":len(valid),"status":"sufficient_common_human_gene_space","time_used_for_feature_space":False}])
    _atomic(decision,STAGE26/"08_stage26_decision.csv",False)
    return {"status":"sufficient_common_human_gene_space","common_genes":len(genes),"contributing_datasets":len(valid)}
