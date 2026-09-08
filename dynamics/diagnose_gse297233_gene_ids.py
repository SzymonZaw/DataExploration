"""Lightweight diagnostic for the GSE297233 gene identifier namespace.

This script deliberately does NOT run the Stage 2.11C pipeline. It reads only
GSE297233 and compares its gene identifiers with targets from decoupler's
PROGENy and DoRothEA networks.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .yamanaka_poc import load_expression

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_control_null"
PATH = DATA / "GSE297233_raw_counts_matrix.csv.gz"


def _gene_id_kind(values):
    values = [str(v) for v in values if str(v) and str(v) != "nan"]
    if not values:
        return "empty"
    n = len(values)
    ensembl = sum(bool(re.fullmatch(r"ENSG\\d+(?:\\.\\d+)?", v, flags=re.I)) for v in values)
    numeric = sum(bool(re.fullmatch(r"\\d+", v)) for v in values)
    if ensembl / n >= 0.8:
        return "ensembl_gene"
    if numeric / n >= 0.8:
        return "numeric_gene_id"
    return "symbol_or_other"


def main():
    import decoupler as dc

    X, _ = load_expression(PATH)
    ids = pd.Index(X.index.astype(str))
    progeny = dc.op.progeny(organism="human", top=100)
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"])

    def network_diag(name, net):
        targets = pd.Index(net["target"].astype(str).unique())
        raw_shared = ids.intersection(targets)
        norm_ids = pd.Index([re.sub(r"\\.\\d+$", "", x).upper() for x in ids])
        norm_targets = pd.Index([re.sub(r"\\.\\d+$", "", x).upper() for x in targets])
        normalized_shared = norm_ids.intersection(norm_targets)
        return {
            "network": name,
            "n_network_targets": int(len(targets)),
            "n_shared_raw": int(len(raw_shared)),
            "n_shared_after_ensembl_version_strip_upper": int(len(normalized_shared)),
            "example_shared_raw": raw_shared[:20].tolist(),
            "example_shared_normalized": normalized_shared[:20].tolist(),
        }

    summary = {
        "dataset": "GSE297233",
        "path": str(PATH.relative_to(ROOT)),
        "n_genes": int(len(ids)),
        "gene_id_kind": _gene_id_kind(ids[: min(1000, len(ids))]),
        "first_30_gene_ids": ids[:30].tolist(),
        "ensembl_like_fraction_first_1000": float(sum(bool(re.fullmatch(r"ENSG\\d+(?:\\.\\d+)?", x, flags=re.I)) for x in ids[:1000]) / max(1, min(1000, len(ids)))),
        "progeny": network_diag("PROGENy", progeny),
        "dorothea": network_diag("DoRothEA", dorothea),
        "interpretation": "This diagnostic tests identifier namespace mismatch only; it does not perform gene annotation or biological inference.",
    }

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "08_gse297233_gene_id_diagnostic.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
