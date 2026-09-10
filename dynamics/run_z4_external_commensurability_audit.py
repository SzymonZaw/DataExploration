"""Frozen pre-transfer commensurability audit for Z4 external specificity.

This script intentionally does NOT compute frozen transfer scores.
It audits only dataset design metadata, temporal structure and feature-space
mapping coverage. Any validation-result-driven adaptation is prohibited.

Expected inputs are simple metadata/sample tables and a deterministic
external-to-discovery gene mapping table supplied before validation results
are inspected.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REQUIRED_GATES = [
    "C1_temporal_commensurability",
    "C2_cell_state_commensurability",
    "C3_intervention_commensurability",
    "C4_feature_space_commensurability",
    "C5_endpoint_independence",
    "C6_nuisance_observability",
]


def _read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() == ".csv":
        return pd.read_csv(p)
    if p.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(p, sep="\t")
    raise ValueError(f"Unsupported table format: {p}")


def _bool_arg(value: str) -> bool:
    value = value.strip().lower()
    if value in {"1", "true", "yes", "y"}:
        return True
    if value in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Expected boolean, got {value!r}")


def _gate(name: str, passed: bool, reason: str) -> dict:
    return {"gate": name, "passed": bool(passed), "reason": reason}


def audit_candidate(
    candidate: str,
    metadata: pd.DataFrame,
    mapping: pd.DataFrame,
    discovery_features: pd.Series,
    *,
    intervention_ok: bool,
    cell_state_ok: bool,
    endpoint_independent: bool,
    nuisance_columns: list[str],
    specificity_control: bool = False,
) -> dict:
    required_cols = {"sample_id", "time"}
    missing = required_cols - set(metadata.columns)
    if missing:
        raise ValueError(f"{candidate}: metadata missing columns {sorted(missing)}")

    n_samples = len(metadata)
    n_times = metadata["time"].nunique(dropna=True)
    donors = metadata["donor"].nunique(dropna=True) if "donor" in metadata else 0

    c1 = n_times >= 3
    c4 = False
    mapping_cols = set(mapping.columns)
    if {"discovery_feature", "external_feature"}.issubset(mapping_cols):
        mapped = set(mapping["discovery_feature"].dropna().astype(str))
        frozen = set(discovery_features.dropna().astype(str))
        coverage = len(mapped & frozen) / max(1, len(frozen))
        c4 = coverage >= 0.50
    else:
        coverage = None

    nuisance_ok = all(col in metadata.columns for col in nuisance_columns)
    c6 = nuisance_ok and n_samples > 0

    gates = [
        _gate("C1_temporal_commensurability", c1, f"n_timepoints={n_times}; require >=3"),
        _gate("C2_cell_state_commensurability", cell_state_ok, "pre-specified biological cell-state relationship"),
        _gate("C3_intervention_commensurability", intervention_ok, "pre-specified OSKM/homologous intervention"),
        _gate("C4_feature_space_commensurability", c4, f"frozen-feature mapping coverage={coverage}"),
        _gate("C5_endpoint_independence", endpoint_independent, "endpoint is independent of transferred score"),
        _gate("C6_nuisance_observability", c6, f"required nuisance columns={nuisance_columns}"),
    ]
    if specificity_control:
        gates.append(_gate("C7_specificity_control_compatibility", True, "candidate is predefined context challenge"))

    passed = all(g["passed"] for g in gates if g["gate"] in REQUIRED_GATES)
    decision = "ELIGIBLE_FOR_FROZEN_TRANSFER" if passed else "INELIGIBLE"

    return {
        "candidate": candidate,
        "decision": decision,
        "n_samples": int(n_samples),
        "n_timepoints": int(n_times),
        "n_donors": int(donors),
        "mapping_coverage": coverage,
        "gates": gates,
        "transfer_score_evaluated": False,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--mapping", required=True)
    ap.add_argument("--discovery-features", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--intervention-ok", required=True)
    ap.add_argument("--cell-state-ok", required=True)
    ap.add_argument("--endpoint-independent", required=True)
    ap.add_argument("--nuisance-columns", nargs="+", required=True)
    ap.add_argument("--specificity-control", action="store_true")
    args = ap.parse_args()

    metadata = _read_table(args.metadata)
    mapping = _read_table(args.mapping)
    features = pd.read_csv(args.discovery_features, header=None).iloc[:, 0]

    result = audit_candidate(
        args.candidate,
        metadata,
        mapping,
        features,
        intervention_ok=_bool_arg(args.intervention_ok),
        cell_state_ok=_bool_arg(args.cell_state_ok),
        endpoint_independent=_bool_arg(args.endpoint_independent),
        nuisance_columns=args.nuisance_columns,
        specificity_control=args.specificity_control,
    )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
