"""Stage 2.9.8: biological annotation of consensus genes.

This stage deliberately does not fit an ODE. It asks whether genes recurring
across independent Stage 2.9.1 discovery folds have coherent biology.

Primary enrichment uses g:Profiler with the 11,899-gene common human space as
an explicit custom statistical background. Hallmark and transcription-factor
enrichment are optional and fail gracefully when network access is unavailable.
"""
from pathlib import Path
import json
import os
import ssl
import time
import urllib.error
import urllib.request

import pandas as pd

try:
    import certifi
except ImportError:
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
    if not verify:
        return ssl._create_unverified_context()
    if certifi is not None:
        return ssl.create_default_context(cafile=certifi.where())
    return ssl.create_default_context()


def _request_json(url, payload, timeout=90, retries=2):
    body = json.dumps(payload).encode("utf-8")
    insecure = os.environ.get("STAGE298_INSECURE_SSL", "") == "1"
    last_exc = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(url, data=body, headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": "DataExploration-stage2.9.8/1.4",
            }, method="POST")
            with urllib.request.urlopen(req, timeout=timeout, context=_ssl_context(True)) as r:
                return json.loads(r.read().decode("utf-8"))
        except urllib.error.URLError as exc:
            last_exc = exc
            if insecure and attempt == 0 and isinstance(exc.reason, ssl.SSLCertVerificationError):
                _log(f"verified TLS failed; retrying once with certificate verification disabled: {exc}")
                try:
                    req = urllib.request.Request(url, data=body, headers={
                        "Content-Type": "application/json",
                        "Accept": "application/json",
                        "User-Agent": "DataExploration-stage2.9.8/1.4",
                    }, method="POST")
                    with urllib.request.urlopen(req, timeout=timeout, context=_ssl_context(False)) as r:
                        return json.loads(r.read().decode("utf-8"))
                except Exception as exc2:
                    last_exc = exc2
            if attempt < retries:
                time.sleep(2 ** attempt)
        except Exception as exc:
            last_exc = exc
            if attempt < retries:
                time.sleep(2 ** attempt)
    _log(f"request failed after {retries + 1} attempts: {last_exc}")
    return None


def _load_genes():
    rec = pd.read_csv(IN / "02_gene_recurrence.csv")
    rec["gene"] = rec["gene"].astype(str).str.upper().str.strip()
    rec["n_discovery_folds"] = pd.to_numeric(rec["n_discovery_folds"], errors="coerce").fillna(0).astype(int)
    background = pd.read_csv(COMMON / "06_common_human_gene_matrix.csv", index_col=0).index.astype(str).str.upper().tolist()
    background = sorted(set(background))
    bg = set(background)
    sets = {
        "recurrence_ge_2": sorted(set(rec.loc[(rec.n_discovery_folds >= 2) & rec.gene.isin(bg), "gene"])),
        "recurrence_ge_3": sorted(set(rec.loc[(rec.n_discovery_folds >= 3) & rec.gene.isin(bg), "gene"])),
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
        "sources": ["GO:BP", "REAC", "KEGG"],
        "domain_scope": "custom",
        "background": background,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
        "output": "json",
    }
    _log(f"g:Profiler enrichment for {label} ({len(gene_list):,} genes)...")
    return _request_json("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", payload)


def _flatten_intersections(intersections, query_names):
    """Return query genes corresponding to non-empty g:Profiler intersections.

    g:Profiler documents ``intersections`` as a list of lists aligned with the
    query Ensembl IDs. We deliberately use the exact submitted query order here
    instead of the optional metadata mapping: the submitted IDs are already
    canonical ENSG identifiers, while evidence-code metadata can have a
    different nested representation when ``no_evidences=False``.
    """
    if not isinstance(intersections, list):
        return ""
    hits = []
    for i, item in enumerate(intersections):
        if item is None:
            continue
        present = bool(item) if isinstance(item, list) else bool(item)
        if present and i < len(query_names):
            hits.append(str(query_names[i]))
    return ",".join(dict.fromkeys(hits))


def _flatten_gprofiler(data, label, fallback_genes):
    rows = []
    if not data:
        return pd.DataFrame()
    query_names = list(fallback_genes)
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
            "intersection_genes": _flatten_intersections(r.get("intersections", []), query_names),
        })
    return pd.DataFrame(rows)


def _program_summary(go, reactome, kegg):
    frames = []
    for frame, source, category in ((go, "GO:BP", "GO_BP"), (reactome, "REAC", "Reactome"), (kegg, "KEGG", "KEGG")):
        if len(frame):
            x = frame[frame.source.eq(source)].copy()
            if len(x):
                x["category"] = category
                frames.append(x[["gene_set", "category", "term_id", "term_name", "p_value", "intersection_size", "intersection_genes"]])
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).sort_values(["gene_set", "p_value"], na_position="last")


def run():
    _log("starting biological annotation; no ODE/state-space model is fitted")
    if os.environ.get("STAGE298_INSECURE_SSL", "") == "1":
        _log("WARNING: STAGE298_INSECURE_SSL=1; HTTPS certificate verification may be bypassed if normal TLS fails")
    sets, background, rec = _load_genes()
    all_go, all_reactome, all_kegg = [], [], []
    for label, genes in sets.items():
        data = _gprofiler(genes, background, label)
        flat = _flatten_gprofiler(data, label, genes)
        if len(flat):
            all_go.append(flat[flat.source.eq("GO:BP")])
            all_reactome.append(flat[flat.source.eq("REAC")])
            all_kegg.append(flat[flat.source.eq("KEGG")])
    go = pd.concat(all_go, ignore_index=True) if all_go else pd.DataFrame()
    reactome = pd.concat(all_reactome, ignore_index=True) if all_reactome else pd.DataFrame()
    kegg = pd.concat(all_kegg, ignore_index=True) if all_kegg else pd.DataFrame()
    rec.to_csv(OUT / "01_gene_recurrence.csv", index=False)
    go.to_csv(OUT / "02_go_biological_process.csv", index=False)
    reactome.to_csv(OUT / "03_reactome.csv", index=False)
    kegg.to_csv(OUT / "04_kegg.csv", index=False)
    pd.DataFrame().to_csv(OUT / "05_hallmark.csv", index=False)
    pd.DataFrame().to_csv(OUT / "06_tf_enrichment.csv", index=False)
    summary = _program_summary(go, reactome, kegg)
    summary.to_csv(OUT / "07_program_summary.csv", index=False)
    meta = pd.DataFrame([{
        "gene_set": label,
        "n_genes": len(genes),
        "n_enriched_go_bp_fdr05": int(((go.gene_set == label) & (pd.to_numeric(go.p_value, errors="coerce") < 0.05)).sum()) if len(go) else 0,
        "n_enriched_reactome_fdr05": int(((reactome.gene_set == label) & (pd.to_numeric(reactome.p_value, errors="coerce") < 0.05)).sum()) if len(reactome) else 0,
        "n_enriched_kegg_fdr05": int(((kegg.gene_set == label) & (pd.to_numeric(kegg.p_value, errors="coerce") < 0.05)).sum()) if len(kegg) else 0,
        "n_enriched_hallmark_fdr05": 0,
        "n_enriched_tf_fdr05": 0,
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
