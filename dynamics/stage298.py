"""Stage 2.9.8: biological annotation of consensus genes.

This stage deliberately does not fit an ODE. It asks whether genes recurring
across independent Stage 2.9.1 discovery folds have coherent biology.

Primary enrichment uses g:Profiler with the 11,899-gene common human space as
background. GO Biological Process, Reactome and KEGG are requested directly.
Hallmark and transcription-factor enrichment are optional and fail gracefully
when network access is unavailable. Results are cached locally so repeated
runs do not need to repeat requests.
"""
from pathlib import Path
import json
import os
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

try:
    import certifi
except ImportError:  # pragma: no cover - standard-library fallback
    certifi = None

ROOT = Path(__file__).resolve().parents[1]
COMMON = ROOT / "results" / "Dynamics" / "stage2_6"
IN = ROOT / "results" / "Dynamics" / "stage2_9_6"
OUT = ROOT / "results" / "Dynamics" / "stage2_9_8"
CACHE = OUT / "cache"
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)


def _log(msg):
    print(f"Stage 2.9.8: {msg}", flush=True)


def _ssl_context(verify=True):
    """Build an HTTPS context, preferring certifi on Windows."""
    if not verify:
        return ssl._create_unverified_context()
    if certifi is not None:
        return ssl.create_default_context(cafile=certifi.where())
    return ssl.create_default_context()


def _allow_insecure_fallback():
    """Permit a last-resort TLS fallback for broken local CA stores.

    This is intentionally opt-in via STAGE298_INSECURE_SSL=1. The normal path
    always verifies certificates. The fallback is useful on isolated Windows
    installations where Python's CA bundle is intercepted by enterprise TLS.
    """
    return os.environ.get("STAGE298_INSECURE_SSL", "0").strip() == "1"


def _request_json(url, payload, timeout=90, retries=2):
    body = json.dumps(payload).encode("utf-8")
    contexts = [_ssl_context(True)]
    if _allow_insecure_fallback():
        contexts.append(_ssl_context(False))
    for context_i, context in enumerate(contexts):
        for attempt in range(retries + 1):
            try:
                req = urllib.request.Request(
                    url,
                    data=body,
                    headers={
                        "Content-Type": "application/json",
                        "User-Agent": "DataExploration-stage2.9.8/1.2",
                    },
                    method="POST",
                )
                with urllib.request.urlopen(req, timeout=timeout, context=context) as r:
                    return json.loads(r.read().decode("utf-8"))
            except Exception as exc:
                if attempt >= retries:
                    if context_i == 0 and len(contexts) > 1:
                        _log(f"verified TLS failed; retrying once with certificate verification disabled: {exc}")
                    else:
                        _log(f"request failed after {retries + 1} attempts: {exc}")
                else:
                    time.sleep(2 ** attempt)
    return None


def _request_text(url, data=None, timeout=90, method="GET"):
    contexts = [_ssl_context(True)]
    if _allow_insecure_fallback():
        contexts.append(_ssl_context(False))
    last_exc = None
    for context_i, context in enumerate(contexts):
        try:
            req = urllib.request.Request(
                url,
                data=data,
                headers={"User-Agent": "DataExploration-stage2.9.8/1.2"},
                method=method,
            )
            with urllib.request.urlopen(req, timeout=timeout, context=context) as r:
                return r.read().decode("utf-8")
        except Exception as exc:
            last_exc = exc
            if context_i == 0 and len(contexts) > 1:
                _log(f"verified TLS failed for {url}; retrying with certificate verification disabled")
    raise last_exc


def _load_genes():
    rec = pd.read_csv(IN / "02_gene_recurrence.csv")
    rec["gene"] = rec["gene"].astype(str).str.upper().str.strip()
    rec["n_discovery_folds"] = pd.to_numeric(rec["n_discovery_folds"], errors="coerce").fillna(0).astype(int)
    background = pd.read_csv(COMMON / "06_common_human_gene_matrix.csv", index_col=0).index.astype(str).str.upper().tolist()
    background = sorted(set(background))
    background_set = set(background)
    sets = {
        "recurrence_ge_2": sorted(set(rec.loc[(rec.n_discovery_folds >= 2) & rec.gene.isin(background_set), "gene"])),
        "recurrence_ge_3": sorted(set(rec.loc[(rec.n_discovery_folds >= 3) & rec.gene.isin(background_set), "gene"])),
    }
    _log(f"loaded {len(background):,} background genes")
    _log(f"consensus sets: >=2 folds = {len(sets['recurrence_ge_2']):,}; >=3 folds = {len(sets['recurrence_ge_3']):,}")
    return sets, background, rec


def _gprofiler(gene_list, background, label):
    cache = CACHE / f"gprofiler_{label}.json"
    if cache.exists():
        try:
            return json.loads(cache.read_text(encoding="utf-8"))
        except Exception:
            pass
    payload = {
        "organism": "hsapiens",
        "query": gene_list,
        "background": background,
        "sources": ["GO:BP", "REAC", "KEGG"],
        "user_threshold": 0.05,
        "significance_threshold_method": "g_SCS",
        "no_evidences": False,
    }
    _log(f"g:Profiler enrichment for {label} ({len(gene_list):,} genes)...")
    data = _request_json("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", payload)
    if data is not None:
        cache.write_text(json.dumps(data), encoding="utf-8")
    return data


def _flatten_gprofiler(data, label):
    rows = []
    if not data:
        return pd.DataFrame()
    for r in data.get("result", []):
        rows.append({
            "gene_set": label,
            "source": r.get("source"),
            "term_id": r.get("native"),
            "term_name": r.get("name"),
            "intersection_size": r.get("intersection_size"),
            "term_size": r.get("term_size"),
            "query_size": r.get("query_size"),
            "effective_domain_size": r.get("effective_domain_size"),
            "p_value": r.get("p_value"),
            "precision": r.get("precision"),
            "recall": r.get("recall"),
            "intersection_genes": ",".join(map(str, r.get("intersections", []))),
        })
    return pd.DataFrame(rows)


def _enrichr_libraries():
    """Return currently available Enrichr library names, if the API permits it."""
    cache = CACHE / "enrichr_libraries.json"
    if cache.exists():
        try:
            return json.loads(cache.read_text(encoding="utf-8"))
        except Exception:
            pass
    try:
        raw = _request_text("https://maayanlab.cloud/Enrichr/datasetStatistics")
        data = json.loads(raw)
        names = []
        if isinstance(data, list):
            for row in data:
                if isinstance(row, dict):
                    name = row.get("libraryName") or row.get("library")
                    if name:
                        names.append(str(name))
        elif isinstance(data, dict):
            for key, value in data.items():
                if isinstance(value, dict) and (value.get("libraryName") or value.get("library")):
                    names.append(str(value.get("libraryName") or value.get("library")))
                elif isinstance(key, str):
                    names.append(key)
        names = sorted(set(names))
        if names:
            cache.write_text(json.dumps(names), encoding="utf-8")
        return names
    except Exception as exc:
        _log(f"could not query Enrichr library list: {exc}")
        return []


def _resolve_library(preferred):
    available = _enrichr_libraries()
    if not available:
        return None
    if preferred in available:
        return preferred
    preferred_upper = preferred.upper()
    for name in available:
        if str(name).upper() == preferred_upper:
            return name
    tokens = {
        "MSigDB_Hallmark_2020": ("MSIGDB", "HALLMARK"),
        "ChEA_2022": ("CHEA",),
    }.get(preferred, ())
    for name in available:
        upper = str(name).upper()
        if tokens and all(token in upper for token in tokens):
            return name
    return None


def _enrichr(gene_list, library, label):
    """Optional Enrichr enrichment; returns empty frame on unavailable network/API."""
    resolved = _resolve_library(library)
    if not resolved:
        _log(f"Enrichr {library}: no compatible current library name found; skipping")
        return []
    cache = CACHE / f"enrichr_{label}_{resolved}.json"
    if cache.exists():
        try:
            return json.loads(cache.read_text(encoding="utf-8"))
        except Exception:
            pass
    genes = "\n".join(gene_list)
    try:
        data = urllib.parse.urlencode({"list": genes, "description": f"DataExploration {label}"}).encode("utf-8")
        raw = _request_text("https://maayanlab.cloud/Enrichr/addList", data=data, method="POST")
        added = json.loads(raw)
        user_list_id = added.get("userListId")
        if not user_list_id:
            _log(f"Enrichr {resolved}: addList returned no userListId")
            return []
        url = "https://maayanlab.cloud/Enrichr/enrich?" + urllib.parse.urlencode({"userListId": user_list_id, "backgroundType": resolved})
        raw = _request_text(url)
        data = json.loads(raw)
        rows = data.get(resolved, [])
        cache.write_text(json.dumps(rows), encoding="utf-8")
        return rows
    except urllib.error.HTTPError as exc:
        _log(f"Enrichr {resolved} unavailable: HTTP {exc.code}")
        return []
    except Exception as exc:
        _log(f"Enrichr {resolved} unavailable: {exc}")
        return []


def _flatten_enrichr(rows, label, library):
    out = []
    for r in rows:
        if not isinstance(r, list) or len(r) < 8:
            continue
        out.append({
            "gene_set": label,
            "library": library,
            "rank": r[0],
            "term_name": r[1],
            "p_value": r[2],
            "z_score": r[3],
            "combined_score": r[4],
            "overlap_genes": ",".join(map(str, r[5])) if isinstance(r[5], list) else str(r[5]),
            "adjusted_p_value": r[6],
            "old_p_value": r[7],
        })
    return pd.DataFrame(out)


def _program_summary(go, reactome, kegg, hallmark, tf):
    frames = []
    if len(go):
        x = go[go.source.eq("GO:BP")].copy(); x["category"] = "GO_BP"; frames.append(x[["gene_set","category","term_id","term_name","p_value","intersection_size","intersection_genes"]])
    if len(reactome):
        x = reactome[reactome.source.eq("REAC")].copy(); x["category"] = "Reactome"; frames.append(x[["gene_set","category","term_id","term_name","p_value","intersection_size","intersection_genes"]])
    if len(kegg):
        x = kegg[kegg.source.eq("KEGG")].copy(); x["category"] = "KEGG"; frames.append(x[["gene_set","category","term_id","term_name","p_value","intersection_size","intersection_genes"]])
    if len(hallmark):
        x = hallmark.copy(); x["category"] = "Hallmark"; frames.append(x[["gene_set","category","term_name","p_value","adjusted_p_value","overlap_genes"]].rename(columns={"overlap_genes":"intersection_genes"}))
    if len(tf):
        x = tf.copy(); x["category"] = "TF"; frames.append(x[["gene_set","category","term_name","p_value","adjusted_p_value","overlap_genes"]].rename(columns={"overlap_genes":"intersection_genes"}))
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["gene_set", "p_value"], na_position="last")


def run():
    _log("starting biological annotation; no ODE/state-space model is fitted")
    if _allow_insecure_fallback():
        _log("WARNING: STAGE298_INSECURE_SSL=1; HTTPS certificate verification may be bypassed if normal TLS fails")
    sets, background, rec = _load_genes()
    all_go, all_reactome, all_kegg = [], [], []
    for label, genes in sets.items():
        gp = _gprofiler(genes, background, label)
        flat = _flatten_gprofiler(gp, label)
        if len(flat):
            all_go.append(flat[flat.source.eq("GO:BP")].copy())
            all_reactome.append(flat[flat.source.eq("REAC")].copy())
            all_kegg.append(flat[flat.source.eq("KEGG")].copy())
    go = pd.concat([x for x in all_go if len(x)], ignore_index=True) if any(len(x) for x in all_go) else pd.DataFrame()
    reactome = pd.concat([x for x in all_reactome if len(x)], ignore_index=True) if any(len(x) for x in all_reactome) else pd.DataFrame()
    kegg = pd.concat([x for x in all_kegg if len(x)], ignore_index=True) if any(len(x) for x in all_kegg) else pd.DataFrame()
    hallmark_rows, tf_rows = [], []
    for label, genes in sets.items():
        hallmark_rows.extend(_flatten_enrichr(_enrichr(genes, "MSigDB_Hallmark_2020", label), label, "MSigDB_Hallmark_2020").to_dict("records"))
        tf_rows.extend(_flatten_enrichr(_enrichr(genes, "ChEA_2022", label), label, "ChEA_2022").to_dict("records"))
    hallmark = pd.DataFrame(hallmark_rows); tf = pd.DataFrame(tf_rows)
    rec.to_csv(OUT / "01_gene_recurrence.csv", index=False)
    go.to_csv(OUT / "02_go_biological_process.csv", index=False)
    reactome.to_csv(OUT / "03_reactome.csv", index=False)
    kegg.to_csv(OUT / "04_kegg.csv", index=False)
    hallmark.to_csv(OUT / "05_hallmark.csv", index=False)
    tf.to_csv(OUT / "06_tf_enrichment.csv", index=False)
    summary = _program_summary(go, reactome, kegg, hallmark, tf)
    summary.to_csv(OUT / "07_program_summary.csv", index=False)
    meta = pd.DataFrame([{
        "gene_set": label,
        "n_genes": len(genes),
        "n_enriched_go_bp_fdr05": int(((go.gene_set == label) & (pd.to_numeric(go.p_value, errors="coerce") < 0.05)).sum()) if len(go) else 0,
        "n_enriched_reactome_fdr05": int(((reactome.gene_set == label) & (pd.to_numeric(reactome.p_value, errors="coerce") < 0.05)).sum()) if len(reactome) else 0,
        "n_enriched_kegg_fdr05": int(((kegg.gene_set == label) & (pd.to_numeric(kegg.p_value, errors="coerce") < 0.05)).sum()) if len(kegg) else 0,
        "n_enriched_hallmark_fdr05": int(((hallmark.gene_set == label) & (pd.to_numeric(hallmark.adjusted_p_value, errors="coerce") < 0.05)).sum()) if len(hallmark) else 0,
        "n_enriched_tf_fdr05": int(((tf.gene_set == label) & (pd.to_numeric(tf.adjusted_p_value, errors="coerce") < 0.05)).sum()) if len(tf) else 0,
    } for label, genes in sets.items()])
    meta.to_csv(OUT / "08_enrichment_summary.csv", index=False)
    _log("complete.")
    print("\nStage 2.9.8 enrichment summary", flush=True)
    print(meta.to_string(index=False), flush=True)
    if len(summary):
        print("\nTop biological terms:", flush=True)
        print(summary.head(30).to_string(index=False), flush=True)
    else:
        print("No enrichment results were returned; inspect network/cache diagnostics.", flush=True)
    return meta


if __name__ == "__main__":
    run()
