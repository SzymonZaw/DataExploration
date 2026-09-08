"""Audit GSE129486 for predeclared H2 nuisance-axis activity.

No Yamanaka comparison, H3 similarity, or biological decision is performed.
The audit asks whether the inflammatory control actually engages the nuisance
axes that make it non-orthogonal to H2: JAK-STAT, NF-kB, and interferon/STAT.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
EXPR = DATA / "GSE129486_rnaseq-data-1_gene-tpm.tsv.gz"
META = DATA / "GSE129486_rnaseq-data-1_metadata.tsv.gz"
EXPECTED_EXPR_SHA = "6e3d7860f4f38d95830a15b8dd570d58f9226170df6002343e432a05c96fcbf0"
EXPECTED_META_SHA = "78a60e353461ab819672e478a9f38b20822fdfc493a1773677ddef3eb167ff04"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_inputs():
    expr = pd.read_csv(EXPR, sep="\t", compression="gzip", low_memory=False)
    meta = pd.read_csv(META, sep="\t", compression="gzip", low_memory=False)
    gene_col = expr.columns[0]
    sample_cols = [c for c in expr.columns[1:]]
    missing = [c for c in sample_cols if c not in set(meta["sample"].astype(str))]
    extra = [c for c in meta["sample"].astype(str) if c not in set(sample_cols)]
    return expr, meta, gene_col, sample_cols, missing, extra


def normalize_ensembl(ids):
    return [str(x).split(".", 1)[0].upper() for x in ids]


def map_to_hgnc(ids):
    import mygene
    info = mygene.MyGeneInfo()
    clean = sorted(set(normalize_ensembl(ids)))
    out = {}
    for i in range(0, len(clean), 1000):
        result = info.querymany(clean[i:i + 1000], scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=False, returnall=False, step=1000, verbose=False)
        for row in result:
            if row.get("notfound"):
                continue
            q, symbol = row.get("query"), row.get("symbol")
            if q and symbol:
                out[str(q).upper()] = str(symbol).upper()
    return out


def pathway_scores(expr_hgnc: pd.DataFrame):
    """Score predeclared PROGENy nuisance pathways across samples.

    Current decoupler releases accept ``data=`` and ``net=`` for ULM;
    source/target/weight are network columns, not method keyword arguments.
    """
    import decoupler as dc

    progeny = dc.op.progeny(organism="human", top=100)
    result = dc.mt.ulm(data=expr_hgnc, net=progeny, tmin=5, verbose=False)
    scores = result[0] if isinstance(result, tuple) else result
    if not isinstance(scores, pd.DataFrame):
        scores = pd.DataFrame(scores)
    scores.index = scores.index.astype(str)

    normalized = {str(c).upper().replace("_", "-"): c for c in scores.columns}
    wanted = {}
    for alias in ("JAK-STAT", "JAK_STAT", "JAKSTAT"):
        key = alias.upper().replace("_", "-")
        if key in normalized:
            wanted["JAK-STAT"] = normalized[key]
            break
    for alias in ("NFKB", "NF-KB", "NF_KB"):
        key = alias.upper().replace("_", "-")
        if key in normalized:
            wanted["NF-kB"] = normalized[key]
            break
    if not wanted:
        return pd.DataFrame(index=expr_hgnc.index)
    return scores.loc[:, list(wanted.values())].rename(columns={v: k for k, v in wanted.items()})


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    expr_sha = sha256(EXPR)
    meta_sha = sha256(META)
    report = {
        "status": "ok",
        "scope": "GSE129486 nuisance-axis audit only; no Yamanaka comparison, H3 similarity, trajectory-agreement statistic, or H3 decision",
        "input_hashes": {
            "GSE129486_gene_tpm": {"sha256": expr_sha, "expected": EXPECTED_EXPR_SHA, "matches": expr_sha == EXPECTED_EXPR_SHA},
            "GSE129486_metadata": {"sha256": meta_sha, "expected": EXPECTED_META_SHA, "matches": meta_sha == EXPECTED_META_SHA},
        },
        "predeclared_h2_axes": ["JAK-STAT", "NF-kB", "interferon/STAT-like inflammatory response"],
    }
    expr, meta, gene_col, sample_cols, missing, extra = load_inputs()
    report["sample_alignment"] = {"expression_samples": len(sample_cols), "metadata_rows": len(meta), "missing_metadata_for_expression_samples": missing, "metadata_samples_not_in_expression": extra}

    mapping = map_to_hgnc(expr[gene_col].tolist())
    symbols = pd.Series(normalize_ensembl(expr[gene_col])).map(mapping)
    valid = symbols.notna()
    mapped = expr.loc[valid, sample_cols].copy()
    mapped.index = symbols.loc[valid].values
    mapped = mapped.groupby(mapped.index).mean().apply(pd.to_numeric, errors="coerce").T
    mapped.index = mapped.index.astype(str)
    report["mapping"] = {
        "input_ensembl_ids": int(len(expr)),
        "mapped_ids": int(len(mapping)),
        "unique_hgnc_symbols": int(mapped.shape[1]),
        "mapping_rate": float(len(mapping) / len(set(normalize_ensembl(expr[gene_col])))),
        "duplicate_symbols_aggregated": int(len(set(symbols.dropna())) - mapped.shape[1]),
    }

    scores = pathway_scores(mapped)
    if scores.empty:
        report["pathway_activity"] = {"status": "NOT_AVAILABLE", "reason": "PROGENy output did not expose the predeclared JAK-STAT/NFkB columns in this decoupler version."}
    else:
        score_meta = meta.copy()
        score_meta["sample"] = score_meta["sample"].astype(str)
        score_meta = score_meta.set_index("sample").reindex(scores.index)
        records = []
        for pathway in scores.columns:
            s = scores[pathway].astype(float)
            rec = {"pathway": pathway, "mean": float(s.mean()), "sd": float(s.std(ddof=1)), "min": float(s.min()), "max": float(s.max()), "abs_mean": float(np.abs(s).mean())}
            for col in ["time", "cell_line", "stimulation"]:
                if col in score_meta.columns:
                    grouped = pd.DataFrame({"score": s, col: score_meta[col].values}).dropna().groupby(col)["score"]
                    rec[col + "_means"] = {str(k): float(v) for k, v in grouped.mean().items()}
            records.append(rec)
        report["pathway_activity"] = {"status": "DESCRIPTIVE", "progeny_pathways": records}

    report["interpretation"] = {
        "current_role": "GSE129486 is an inflammatory temporal control with a predeclared H2 risk; this audit quantifies pathway engagement without using Yamanaka similarity.",
        "decision_status": "NO_H3_DECISION",
        "next_step": "If pathway scoring is available, use the persisted JAK-STAT/NFkB activity distributions to predefine an orthogonality rule; then audit interferon/STAT-like DoRothEA activity separately before any H3 similarity.",
    }
    out = OUT / "H3_GSE129486_H2_NUISANCE_AUDIT.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"WROTE {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
