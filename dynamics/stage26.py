"""Robust Stage 2.6 common human-gene space construction.

This module deliberately keeps the validated previous common-space files intact when
refreshing fails or produces an insufficient intersection. BioMart requests use POST
first and recursively split on HTTP 414; MyGene is the fallback mapping service.
"""
from pathlib import Path
import json
import re
import ssl
import urllib.parse
import urllib.request
from io import StringIO
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


def normalize_symbol(x):
    s = str(x).strip().strip('"').upper()
    return None if s in {"", "NAN", "NA", "NONE", "NULL", "-", "."} else s


def clean_refseq(x):
    s = str(x).strip().strip('"')
    s = re.sub(r"_at$", "", s, flags=re.I)
    return re.sub(r"\.\d+$", "", s).upper()


def extract_ensembl(x):
    m = re.search(r"(ENS[A-Z]*G\d+)", str(x).strip().strip('"'), re.I)
    return m.group(1).upper() if m else None


def extract_refseq(x):
    m = re.search(r"((?:NM|NR|XM|XR)_\d+(?:\.\d+)?)", str(x), re.I)
    return clean_refseq(m.group(1)) if m else None


def read_feature_matrix(path):
    df = pd.read_csv(path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=0, how="all").dropna(axis=1, how="all")
    df.index = df.index.astype(str).str.strip().str.strip('"')
    return df.groupby(level=0).mean()


def _request(url, data=None, timeout=180):
    headers = {"User-Agent": "DataExploration-Dynamics/6.0", "Accept": "text/plain,application/json"}
    method = "GET"
    body = None
    if data is not None:
        body = data if isinstance(data, bytes) else data.encode("utf-8")
        headers["Content-Type"] = "application/x-www-form-urlencoded"
        method = "POST"
    req = urllib.request.Request(url, data=body, headers=headers, method=method)
    for context in (None, ssl._create_unverified_context()):
        try:
            kw = {"timeout": timeout}
            if context is not None:
                kw["context"] = context
            with urllib.request.urlopen(req, **kw) as r:
                return r.read().decode("utf-8", "replace")
        except Exception:
            continue
    raise RuntimeError(f"request failed: {url}")


def _biomart_xml(dataset, filter_name, ids, attributes):
    vals = ",".join(ids).replace("&", "&amp;").replace('"', "&quot;")
    attrs = "".join(f'<Attribute name="{a}"/>' for a in attributes)
    return (f'<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query '
            f'virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" '
            f'count="" datasetConfigVersion="0.6"><Dataset name="{dataset}" interface="default">'
            f'<Filter name="{filter_name}" value="{vals}"/>{attrs}</Dataset></Query>')


def biomart_query(dataset, filter_name, ids, attributes, label, batch_no=1):
    """BioMart POST with recursive splitting for URI-size failures."""
    if not ids:
        return ""
    xml = _biomart_xml(dataset, filter_name, ids, attributes)
    endpoints = ["https://www.ensembl.org/biomart/martservice", "https://useast.ensembl.org/biomart/martservice"]
    errors = []
    for endpoint in endpoints:
        try:
            txt = _request(endpoint, data=xml)
            if txt.strip() and not txt.lstrip().lower().startswith(("query error", "error")):
                return txt
        except Exception as exc:
            errors.append(str(exc))
    if len(ids) > 25:
        mid = len(ids) // 2
        print(f"  BioMart {label} batch {batch_no}: request failed; splitting {len(ids)} -> {mid}+{len(ids)-mid}")
        return biomart_query(dataset, filter_name, ids[:mid], attributes, label, batch_no * 2) + biomart_query(dataset, filter_name, ids[mid:], attributes, label, batch_no * 2 + 1)
    print(f"  BioMart {label} batch {batch_no} failed for {len(ids)} IDs; fallback will be used.")
    return ""


def mygene_query(ids, scopes, species, fields):
    rows = []
    ids = list(dict.fromkeys(str(x) for x in ids if str(x)))
    for start in range(0, len(ids), 500):
        batch = ids[start:start + 500]
        form = urllib.parse.urlencode({"q": ",".join(batch), "scopes": scopes, "fields": fields, "species": species, "size": 500})
        result = None
        try:
            text = _request("https://mygene.info/v3/query", data=form)
            result = json.loads(text)
        except Exception as exc:
            print(f"  MyGene batch {start // 500 + 1} failed: {exc}")
        hits = result if isinstance(result, list) else (result.get("out", result.get("hits", [])) if isinstance(result, dict) else [])
        rows.extend(hits)
    return rows


def platform_mapping(gpl):
    cache = CACHE / f"{gpl}_platform_mapping.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t", dtype=str).fillna("")
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&view=full&form=text"
    text = _request(url, timeout=120)
    m = re.search(r"!platform_table_begin\n(.*?)\n!platform_table_end", text, re.S)
    if not m:
        raise RuntimeError(f"no platform table found for {gpl}")
    tab = pd.read_csv(StringIO(m.group(1)), sep="\t", dtype=str)
    cols = {str(c).lower().strip(): c for c in tab.columns}
    idcol = next((cols[k] for k in ("id", "probe id", "probeid") if k in cols), tab.columns[0])
    sym = next((c for k, c in cols.items() if any(q in k for q in ("gene symbol", "gene_symbol", "symbol", "gene assignment"))), None)
    if sym is None:
        raise RuntimeError(f"no gene-symbol column for {gpl}")
    out = pd.DataFrame({"feature": tab[idcol].astype(str), "human_gene": tab[sym].map(normalize_symbol)})
    out = out.dropna(subset=["human_gene"]).drop_duplicates("feature")
    out.to_csv(cache, sep="\t", index=False)
    return out


def map_human_ensembl(ids):
    ids = sorted({extract_ensembl(x) for x in ids if extract_ensembl(x)})
    rows = []
    txt = biomart_query("hsapiens_gene_ensembl", "ensembl_gene_id", ids, ["ensembl_gene_id", "external_gene_name"], "human Ensembl")
    for line in txt.splitlines():
        p = line.split("\t")
        if len(p) >= 2 and p[0] and normalize_symbol(p[1]):
            rows.append((extract_ensembl(p[0]), normalize_symbol(p[1])))
    out = dict(rows)
    missing = [x for x in ids if x not in out]
    for hit in mygene_query(missing, "ensembl.gene", "human", "symbol"):
        q, s = extract_ensembl(hit.get("query", "")), normalize_symbol(hit.get("symbol"))
        if q and s:
            out[q] = s
    return out


def map_mouse_refseq(ids):
    ids = sorted({clean_refseq(x) for x in ids if extract_refseq(x)})
    rows = []
    txt = biomart_query("mmusculus_gene_ensembl", "refseq_mrna", ids, ["refseq_mrna", "external_gene_name", "hsapiens_homolog_associated_gene_name"], "mouse RefSeq")
    for line in txt.splitlines():
        p = line.split("\t")
        if len(p) >= 3 and p[0] and normalize_symbol(p[2]):
            rows.append((clean_refseq(p[0]), normalize_symbol(p[2])))
    out = dict(rows)
    missing = [x for x in ids if x not in out]
    hits = mygene_query(missing, "refseq", "mouse", "symbol,homologene")
    human_entrez = set()
    temp = []
    for hit in hits:
        q = clean_refseq(hit.get("query", "")); hom = hit.get("homologene", {})
        genes = hom.get("genes", []) if isinstance(hom, dict) else []
        ent = [str(p[1]) for p in genes if isinstance(p, (list, tuple)) and len(p) >= 2 and str(p[0]) == "9606"]
        human_entrez.update(ent)
        temp.append((q, ent))
    symbols = {}
    for hit in mygene_query(sorted(human_entrez), "entrezgene", "human", "symbol"):
        s = normalize_symbol(hit.get("symbol"))
        if s:
            symbols[str(hit.get("query", ""))] = s
    for q, ents in temp:
        for ent in ents:
            if ent in symbols:
                out[q] = symbols[ent]
                break
    return out


def direct_map(df, dataset):
    acc = {}
    source = {}
    ambiguous = 0
    unmapped = 0
    for idx, row in df.iterrows():
        s = str(idx).strip().strip('"')
        gene = None
        ens = extract_ensembl(s)
        ref = extract_refseq(s)
        if ens:
            source.setdefault("ensembl_pending", []).append((s, ens))
        elif ref:
            source.setdefault("refseq_pending", []).append((s, ref))
        else:
            gene = normalize_symbol(s)
            if gene:
                source.setdefault("direct", []).append(s)
        if gene:
            acc.setdefault(gene, []).append(row.to_numpy(float))
    ens_map = map_human_ensembl([x[1] for x in source.get("ensembl_pending", [])]) if source.get("ensembl_pending") else {}
    ref_map = map_mouse_refseq([x[1] for x in source.get("refseq_pending", [])]) if source.get("refseq_pending") else {}
    mapped = 0
    for s in source.get("direct", []):
        acc.setdefault(normalize_symbol(s), []).append(df.loc[s].to_numpy(float)); mapped += 1
    for s, ens in source.get("ensembl_pending", []):
        gene = ens_map.get(ens)
        if gene:
            acc.setdefault(gene, []).append(df.loc[s].to_numpy(float)); mapped += 1
        else:
            unmapped += 1
    for s, ref in source.get("refseq_pending", []):
        gene = ref_map.get(ref)
        if gene:
            acc.setdefault(gene, []).append(df.loc[s].to_numpy(float)); mapped += 1
        else:
            unmapped += 1
    mat = pd.DataFrame({g: np.mean(v, axis=0) for g, v in acc.items()}, index=df.columns).T if acc else pd.DataFrame(index=[], columns=df.columns)
    return mat, mapped, unmapped, len(acc)


def map_dataset(ds, path):
    df = read_feature_matrix(path)
    if ds in PLATFORMS:
        ann = platform_mapping(PLATFORMS[ds]).set_index("feature")["human_gene"].to_dict()
        acc = {}
        mapped = 0
        for idx, row in df.iterrows():
            gene = normalize_symbol(ann.get(str(idx)))
            if gene:
                acc.setdefault(gene, []).append(row.to_numpy(float)); mapped += 1
        mat = pd.DataFrame({g: np.mean(v, axis=0) for g, v in acc.items()}, index=df.columns).T if acc else pd.DataFrame(index=[], columns=df.columns)
        return df, mat, mapped, len(df) - mapped, len(acc), "GEO_platform"
    mat, mapped, unmapped, genes = direct_map(df, ds)
    return df, mat, mapped, unmapped, genes, "direct+BioMart/MyGene"


def _write_if_atomic(df, path):
    tmp = path.with_suffix(path.suffix + ".tmp")
    df.to_csv(tmp, index=True)
    tmp.replace(path)


def stage2_6_robust(min_genes=1000, min_datasets=4):
    rows = []
    matrices = {}
    for ds, path in FEATURE_FILES.items():
        try:
            raw, mat, mapped, unmapped, genes, source = map_dataset(ds, path)
            coverage = mapped / max(len(raw), 1)
            status = "mapped" if mapped else "no_gene_mapping"
            matrices[ds] = mat
            rows.append({"dataset": ds, "n_features": len(raw), "n_samples": raw.shape[1], "mapped_features": mapped, "unmapped_features": unmapped, "mapped_genes": genes, "mapping_coverage": coverage, "mapping_source": source, "status": status})
        except Exception as exc:
            print(f"Stage 2.6 {ds} mapping failed: {type(exc).__name__}: {exc}")
            matrices[ds] = pd.DataFrame()
            rows.append({"dataset": ds, "n_features": 0, "n_samples": 0, "mapped_features": 0, "unmapped_features": 0, "mapped_genes": 0, "mapping_coverage": 0.0, "mapping_source": "failed", "status": f"failed:{type(exc).__name__}"})
    audit = pd.DataFrame(rows)
    _write_if_atomic(audit, STAGE26 / "05_mapping_audit.csv")
    valid = {ds: m for ds, m in matrices.items() if isinstance(m, pd.DataFrame) and not m.empty and len(m.index) > 0}
    common = set.intersection(*(set(m.index) for m in valid.values())) if valid else set()
    contributing = len(valid)
    print(f"Stage 2.6: common human genes={len(common)}, contributing datasets={contributing}")
    if len(common) < min_genes or contributing < min_datasets:
        msg = f"insufficient common human gene space ({len(common)} genes across {contributing} datasets); previous valid files preserved"
        print("WARNING: " + msg)
        return {"status": "insufficient_common_human_gene_space", "common_genes": len(common), "contributing_datasets": contributing, "audit": audit}
    genes = sorted(common)
    blocks = []
    metadata = []
    for ds, mat in valid.items():
        block = mat.loc[genes].copy()
        block = block.apply(lambda r: (r - r.mean()) / r.std(ddof=0) if r.std(ddof=0) > 0 else np.nan, axis=1)
        block.columns = [f"{ds}__{i}" for i in range(block.shape[1])]
        blocks.append(block)
        for i, col in enumerate(mat.columns):
            metadata.append({"dataset": ds, "sample": str(col), "matrix_column": f"{ds}__{i}", "sample_index": i})
    common_matrix = pd.concat(blocks, axis=1)
    meta = pd.DataFrame(metadata)
    _write_if_atomic(common_matrix, STAGE26 / "06_common_human_gene_matrix.csv")
    _write_if_atomic(meta, STAGE26 / "07_common_gene_sample_metadata.csv")
    return {"status": "sufficient_common_human_gene_space", "common_genes": len(genes), "contributing_datasets": contributing, "audit": audit}


__all__ = ["stage2_6_robust"]
