"""Audit identifier compatibility between GSE67462 expression and GSE67520 mm9 TSS.

Diagnostic only. This audit does not modify frozen Z6 support or any scientific
threshold. It determines whether the unexpectedly tiny 2-11 gene overlap is a
true biological intersection or an identifier/namespace mismatch.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import pandas as pd

from dynamics import validation
from dynamics.run_z6_gse67462_multimodal_validation_analysis import _load_tss


def _norm(value: str) -> str:
    value = str(value).strip().strip('"')
    value = re.sub(r"\.[0-9]+$", "", value)
    return value.lower()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", default="Data/GSE67520/mm9.refGene.gtf.gz")
    parser.add_argument("--output", default="results/Dynamics/z6_gse67462_identifier_audit")
    args = parser.parse_args()

    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    matrix, metadata = validation._load_common_space()
    tss = _load_tss(Path(args.gtf))

    expr_ids = pd.Index([str(x) for x in matrix.index])
    tss_ids = pd.Index([str(x) for x in tss["gene"]])

    exact = set(expr_ids) & set(tss_ids)
    norm_expr = {_norm(x): x for x in expr_ids}
    norm_tss = {_norm(x): x for x in tss_ids}
    norm_overlap = sorted(set(norm_expr) & set(norm_tss))

    expr_upper = {str(x).upper(): x for x in expr_ids}
    tss_upper = {str(x).upper(): x for x in tss_ids}
    upper_overlap = sorted(set(expr_upper) & set(tss_upper))

    report = {
        "expression_matrix_shape": [int(matrix.shape[0]), int(matrix.shape[1])],
        "expression_identifier_examples": expr_ids[:30].tolist(),
        "expression_identifier_examples_with_type": [type(x).__name__ for x in matrix.index[:10]],
        "tss_count": int(len(tss_ids)),
        "tss_identifier_examples": tss_ids[:30].tolist(),
        "exact_overlap": int(len(exact)),
        "case_insensitive_overlap": int(len(upper_overlap)),
        "version_stripped_case_insensitive_overlap": int(len(norm_overlap)),
        "exact_overlap_examples": sorted(exact)[:50],
        "normalized_overlap_examples": [norm_expr[x] for x in norm_overlap[:50]],
    }

    pd.DataFrame({
        "expression_id": [norm_expr[x] for x in norm_overlap],
        "tss_id": [norm_tss[x] for x in norm_overlap],
        "normalized_id": norm_overlap,
    }).to_csv(out / "01_identifier_overlap.csv", index=False)

    (out / "02_identifier_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    print("GSE67462/GSE67520 identifier audit complete.")
    print(f"Expression matrix: {matrix.shape[0]} genes x {matrix.shape[1]} samples")
    print(f"Expression ID examples: {expr_ids[:15].tolist()}")
    print(f"GSE67520 TSS genes: {len(tss_ids)}")
    print(f"TSS ID examples: {tss_ids[:15].tolist()}")
    print(f"Exact overlap: {len(exact)}")
    print(f"Case-insensitive overlap: {len(upper_overlap)}")
    print(f"Version-stripped overlap: {len(norm_overlap)}")
    if norm_overlap:
        print(f"Overlap examples: {[norm_expr[x] for x in norm_overlap[:20]]}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
