"""Lightweight diagnostic for the GSE297233 gene identifier namespace.

This script deliberately does NOT run the Stage 2.11C pipeline. It reads only
GSE297233 and compares its gene identifiers with targets from decoupler's
PROGENy and DoRothEA networks. It also tests the Ensembl -> HGNC mapping used
by ``dynamics.yamanaka_poc._signature``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

from .yamanaka_poc import _gene_id_kind, _map_ensembl_to_symbols, load_expression

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_control_null"
PATH = DATA / "GSE297233_raw_counts_matrix.csv.gz"


def _canonical(value: object) -> str:
    s = str(value).strip()
    if re.fullmatch(r"ENSG\d+(?:\.\d+)?", s, flags=re.IGNORECASE):
        s = s.split(".", 1)[0]
    return s.upper()


def _network_diag(name, ids, net):
    targets = pd.Index(net["target"].astype(str).unique())
    raw_shared = ids.intersection(targets)
    normalized_shared = pd.Index([_canonical(x) for x in ids]).intersection(
        pd.Index([_canonical(x) for x in targets])
    )
    return {
        "network": name,
        "n_network_targets": int(len(targets)),
        "n_shared_raw": int(len(raw_shared)),
        "n_shared_after_format_normalization": int(len(normalized_shared)),
        "example_shared_raw": raw_shared[:20].tolist(),
        "example_shared_normalized": normalized_shared[:20].tolist(),
    }


def main():
    import decoupler as dc

    X, _ = load_expression(PATH)
    ids = pd.Index(X.index.astype(str))
    mapped, mapping_diag = _map_ensembl_to_symbols(pd.Series(0.0, index=ids))
    mapped_ids = pd.Index(mapped.index.astype(str))

    progeny = dc.op.progeny(organism="human", top=100)
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"])

    summary = {
        "dataset": "GSE297233",
        "path": str(PATH.relative_to(ROOT)),
        "n_genes": int(len(ids)),
        "gene_id_kind": _gene_id_kind(ids),
        "first_30_gene_ids": ids[:30].tolist(),
        "ensembl_like_fraction_first_1000": float(
            sum(bool(re.fullmatch(r"ENSG\d+(?:\.\d+)?", x, flags=re.IGNORECASE)) for x in ids[:1000])
            / max(1, min(1000, len(ids)))
        ),
        "mapping": mapping_diag,
        "n_mapped_symbol_ids": int(len(mapped_ids)),
        "first_30_mapped_symbols": mapped_ids[:30].tolist(),
        "progeny_before_mapping": _network_diag("PROGENy", ids, progeny),
        "progeny_after_mapping": _network_diag("PROGENy", mapped_ids, progeny),
        "dorothea_before_mapping": _network_diag("DoRothEA", ids, dorothea),
        "dorothea_after_mapping": _network_diag("DoRothEA", mapped_ids, dorothea),
        "interpretation": "This diagnostic tests identifier namespace compatibility; it does not perform biological inference.",
    }

    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "08_gse297233_gene_id_diagnostic.json"
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
