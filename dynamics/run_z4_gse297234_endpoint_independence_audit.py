"""Audit endpoint independence for frozen Z4 validation on GSE297234.

This script does not create a Z4-positive result. It classifies candidate
variables according to provenance and leakage risk, using the already-produced
sample metadata and frozen transfer scores. Expression-derived annotations are
not treated as independent biological endpoints.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

DEFAULT_OUT = "results/Dynamics/z4_gse297234_frozen_transfer_audit"
EXPECTED_SAMPLES = [
    "GM00731_D0", "GM00731_D3", "GM00731_D7", "GM00731_D10",
    "GM23815_D0", "GM23815_D3", "GM23815_D7", "GM23815_D10",
]

# Provenance is deliberately conservative. These variables are known from the
# current RDS inspection to coexist with the expression-derived object. Unless
# an external measurement is documented separately, they cannot establish
# endpoint independence.
CANDIDATES = {
    "cell_state": ("derived_from_same_expression_object", "state annotation stored in Seurat metadata"),
    "PartialReprog1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "NonReprog1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "EarlyPluripotency1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "Pluripotency1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "Fibroblast1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "PI16_fibtype1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "LRRC15_fibtype1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "COL3A1_fibtype1": ("derived_from_same_expression_object", "expression-derived score/annotation"),
    "HALLMARK_EMT1": ("derived_from_same_expression_object", "expression-derived pathway score"),
    "HALLMARK_TGFB1": ("derived_from_same_expression_object", "expression-derived pathway score"),
    "ordered_clusters": ("derived_from_same_expression_object", "cluster/state assignment in Seurat object"),
    "age_group": ("raw_metadata", "donor/age metadata"),
    "age_ident": ("raw_metadata", "donor identity metadata"),
    "orig.ident": ("raw_metadata", "sample identity metadata"),
}


def _constantness(df: pd.DataFrame, col: str) -> tuple[int, int, float]:
    if col not in df.columns or "sample" not in df.columns:
        return 0, 0, float("nan")
    per_sample = df.groupby("sample")[col].nunique(dropna=True)
    return int(per_sample.size), int((per_sample <= 1).sum()), float((per_sample <= 1).mean()) if len(per_sample) else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default=DEFAULT_OUT)
    args = ap.parse_args()
    out = Path(args.output)
    meta_path = out / "02_gse297234_sample_metadata.csv"
    score_path = out / "04_frozen_module_scores_with_metadata.csv"
    transfer_summary = out / "05_transfer_summary.csv"
    structure_path = Path("results/Dynamics/z4_gse297234_rds_structure_audit/structure.json")
    geo_path = Path("results/Dynamics/z4_external_metadata_audit/GSE297234/metadata.csv")

    if not meta_path.exists():
        raise SystemExit(f"Missing required input: {meta_path}")
    meta = pd.read_csv(meta_path)
    scores = pd.read_csv(score_path) if score_path.exists() else pd.DataFrame()
    transfer = pd.read_csv(transfer_summary) if transfer_summary.exists() else pd.DataFrame()

    observed = sorted(meta["sample"].astype(str).unique()) if "sample" in meta.columns else []
    missing = [s for s in EXPECTED_SAMPLES if s not in observed]
    unexpected = [s for s in observed if s not in EXPECTED_SAMPLES]
    sample_gate = not missing and not unexpected and len(observed) == len(EXPECTED_SAMPLES)

    rows = []
    for col, (provenance, reason) in CANDIDATES.items():
        if col not in meta.columns:
            available = False
            n_samples = n_constant = 0
            frac_constant = float("nan")
        else:
            available = True
            n_samples, n_constant, frac_constant = _constantness(meta, col)
        independent = provenance == "independent_external_measurement"
        leakage = provenance == "derived_from_same_expression_object"
        usable = bool(available and independent and not leakage and n_samples >= 3)
        rows.append({
            "candidate": col,
            "available_in_sample_metadata": available,
            "provenance_class": provenance,
            "provenance_reason": reason,
            "sample_count": n_samples,
            "samples_constant_within_sample": n_constant,
            "fraction_samples_constant": frac_constant,
            "same_expression_leakage_risk": leakage,
            "independent_endpoint_eligible": usable,
            "target_vs_nontarget_testable": False,
            "decision": "NOT_ELIGIBLE_UNTIL_INDEPENDENT_PROVENANCE" if not usable else "ELIGIBLE_FOR_Z4_REVIEW",
        })

    endpoint_df = pd.DataFrame(rows)
    endpoint_df.to_csv(out / "07_endpoint_candidates.csv", index=False)

    availability = []
    for sample in EXPECTED_SAMPLES:
        row = {"sample": sample, "present_in_transfer_metadata": sample in observed}
        if not scores.empty and "sample" in scores.columns:
            row["module_score_rows"] = int((scores["sample"] == sample).sum())
        else:
            row["module_score_rows"] = 0
        availability.append(row)
    pd.DataFrame(availability).to_csv(out / "08_endpoint_sample_availability.csv", index=False)

    result = {
        "candidate": "GSE297234",
        "sample_resolution": "PASS" if sample_gate else "FAIL",
        "expected_samples": EXPECTED_SAMPLES,
        "observed_samples": observed,
        "missing_samples": missing,
        "unexpected_samples": unexpected,
        "candidate_endpoint_count": int(len(endpoint_df)),
        "independent_external_measurements_found": int((endpoint_df["provenance_class"] == "independent_external_measurement").sum()),
        "expression_derived_candidates": int(endpoint_df["same_expression_leakage_risk"].sum()),
        "target_vs_nontarget_independent_endpoint": False,
        "endpoint_independence": "NOT_ESTABLISHED",
        "decision": "Z4_UNRESOLVED",
        "interpretation": "TRANSFER_ASSOCIATION_CAN_BE_REPORTED_BUT_CURRENT_RDS_ANNOTATIONS_CANNOT_ESTABLISH_INDEPENDENT_ENDPOINT",
        "frozen_z4_rule_unchanged": True,
        "inputs": {
            "sample_metadata": str(meta_path),
            "transfer_scores": str(score_path) if score_path.exists() else None,
            "transfer_summary": str(transfer_summary) if transfer_summary.exists() else None,
            "rds_structure": str(structure_path) if structure_path.exists() else None,
            "geo_metadata": str(geo_path) if geo_path.exists() else None,
        },
    }
    (out / "09_endpoint_independence_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("GSE297234 Z4 ENDPOINT-INDEPENDENCE AUDIT")
    print(f"sample resolution: {result['sample_resolution']}")
    print(f"independent external measurements found: {result['independent_external_measurements_found']}")
    print(f"expression-derived candidates: {result['expression_derived_candidates']}")
    print("endpoint independence: NOT ESTABLISHED")
    print("decision: Z4_UNRESOLVED")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
