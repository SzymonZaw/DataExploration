"""Lightweight audit of the existing Stage 2.11C outputs.

This module deliberately does NOT rerun Stage 2.11C, load raw expression matrices,
or recompute pathway/TF activities. It audits the files already produced by the
pipeline plus the Git history and source/protocol consistency.

Audit order:
1. mapping integrity (D): inspect the recorded Ensembl -> HGNC mapping evidence;
2. preregistration provenance (A): inspect Git history and output timestamps;
3. statistical structure (B): distinguish deterministic Tier-B filtering from
   inferential tests and calculate BH q-values for the existing Null-1 p-values;
4. candidate audit (C): identify the surviving candidate only after D/A/B.

The audit is intentionally conservative. It never upgrades a candidate or
changes the historical Stage 2.11C result in-place.
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_control_null"
PROTOCOL = ROOT / "dynamics" / "STAGE_2_11C_CONTROL_AND_NULL_PROTOCOL.md"
CONTROL_SOURCE = ROOT / "dynamics" / "yamanaka_control_null.py"
SUMMARY = OUT / "stage2_11c_summary.json"
NULL1 = OUT / "02_null1_temporal_order.csv"
CONTEXT = OUT / "04_context_control_OSK.csv"
STABILITY = OUT / "05_candidate_stability.csv"
TIERS = OUT / "06_candidate_tiers.csv"
PERT_AUDIT = OUT / "07_perturbation_activity_audit.csv"
MAPPING_DIAG = OUT / "08_gse297233_gene_id_diagnostic.json"

PROTOCOL_COMMIT = "f90972f5f312b0e4ba39ec64693ccd63b7807c79"
EXPECTED_DECISIONS = {"PROCEED", "MIXED", "SPECIFICITY_UNSUPPORTED", "INCONCLUSIVE"}


def _sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def _git(*args: str) -> str:
    try:
        return subprocess.check_output(
            ["git", *args], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return ""


def _git_commit_meta(commit: str) -> dict:
    raw = _git("show", "-s", "--format=%H%x1f%aI%x1f%cI%x1f%s", commit)
    parts = raw.split("\x1f", 3)
    if len(parts) != 4:
        return {"commit": commit, "available": False}
    return {
        "commit": parts[0],
        "author_date": parts[1],
        "committer_date": parts[2],
        "message": parts[3],
        "available": True,
    }


def _file_mtime(path: Path) -> str | None:
    if not path.exists():
        return None
    return datetime.fromtimestamp(path.stat().st_mtime, tz=timezone.utc).isoformat()


def _bh(values: pd.Series) -> pd.Series:
    """Benjamini-Hochberg adjustment, preserving the input index."""
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if not valid.any():
        return out
    pv = p.loc[valid].to_numpy(float)
    n = len(pv)
    order = np.argsort(pv, kind="mergesort")
    ranked = pv[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)
    restored = np.empty(n, dtype=float)
    restored[order] = q
    out.loc[valid] = restored
    return out


def _mapping_audit() -> dict:
    if not MAPPING_DIAG.exists():
        return {
            "status": "missing",
            "independent_mapping_validation": "not_available",
            "interpretation": "No saved GSE297233 identifier diagnostic is available.",
        }
    data = json.loads(MAPPING_DIAG.read_text(encoding="utf-8"))
    mapping = data.get("mapping", {})
    n_input = int(data.get("n_genes", 0))
    n_mapped = int(mapping.get("n_mapped_genes", 0))
    n_unmapped = int(mapping.get("n_unmapped_genes", 0))
    n_unique_symbols = int(data.get("n_mapped_symbol_ids", 0))
    collisions = max(0, n_mapped - n_unique_symbols)
    progeny = data.get("progeny_after_mapping", {})
    dorothea = data.get("dorothea_after_mapping", {})
    return {
        "status": "recorded",
        "gene_id_kind": data.get("gene_id_kind"),
        "n_input_gene_ids": n_input,
        "n_mapped_gene_ids": n_mapped,
        "n_unmapped_gene_ids": n_unmapped,
        "n_unique_mapped_symbols": n_unique_symbols,
        "mapping_collision_rows": collisions,
        "mapping_collision_present": collisions > 0,
        "mapping_rate": float(n_mapped / n_input) if n_input else None,
        "progeny_shared_after_mapping": int(progeny.get("n_shared_raw", 0)),
        "dorothea_shared_after_mapping": int(dorothea.get("n_shared_raw", 0)),
        "aggregation_in_code": "duplicate canonical HGNC symbols are averaged",
        "independent_mapping_validation": "partial",
        "limitations": [
            "The saved diagnostic establishes namespace compatibility and recorded mapping counts, but does not independently validate every Ensembl->HGNC mapping against a second annotation source.",
            "A non-zero collision count means multiple mapped input rows collapse onto the same HGNC symbol; the current code averages those rows before network scoring.",
            "MyGene mapping is an external annotation service, so a future reproducibility lock should persist the exact mapping table or annotation release/version.",
        ],
    }


def _provenance_audit() -> dict:
    protocol_meta = _git_commit_meta(PROTOCOL_COMMIT)
    head = _git("rev-parse", "HEAD")
    head_meta = _git_commit_meta(head) if head else {"available": False}
    protocol_hash = _sha256(PROTOCOL)
    summary_hash = None
    if SUMMARY.exists():
        try:
            summary_hash = json.loads(SUMMARY.read_text(encoding="utf-8")).get("protocol_hash")
        except Exception:
            pass
    output_mtime = _file_mtime(SUMMARY)
    protocol_before_output = None
    if protocol_meta.get("committer_date") and output_mtime:
        try:
            pdt = datetime.fromisoformat(protocol_meta["committer_date"].replace("Z", "+00:00"))
            odt = datetime.fromisoformat(output_mtime)
            protocol_before_output = pdt <= odt
        except Exception:
            protocol_before_output = None
    history_ok = bool(protocol_meta.get("available")) and (
        not head or bool(_git("merge-base", "--is-ancestor", PROTOCOL_COMMIT, head))
    )
    return {
        "protocol_commit": protocol_meta,
        "current_head": head_meta,
        "protocol_hash_current": protocol_hash,
        "protocol_hash_recorded_in_summary": summary_hash,
        "summary_output_mtime_utc": output_mtime,
        "protocol_commit_precedes_output_file_mtime": protocol_before_output,
        "protocol_commit_reachable_from_current_head": history_ok,
        "provenance_conclusion": (
            "The decision thresholds are present in a Git commit dated before the existing output file timestamp. This supports pre-specification of the thresholds for this local result, although the output file timestamp is not a cryptographic record of execution time."
            if protocol_before_output is True
            else "Pre-specification timing could not be established from the available Git/output metadata."
        ),
    }


def _implementation_audit() -> dict:
    protocol = PROTOCOL.read_text(encoding="utf-8") if PROTOCOL.exists() else ""
    source = CONTROL_SOURCE.read_text(encoding="utf-8") if CONTROL_SOURCE.exists() else ""
    protocol_p01 = bool(re.search(r"p-value is \*\*< 0\.01\*\*", protocol))
    source_p005 = "p<=.05" in source or "p <= .05" in source
    source_confound_label = "CONFOUNDED_UNSUPPORTED" in source
    return {
        "status": "checked",
        "protocol_requires_null1_p_lt_001": protocol_p01,
        "implementation_uses_null1_p_le_005": source_p005,
        "null1_threshold_mismatch": protocol_p01 and source_p005,
        "implementation_uses_confound_label": source_confound_label,
        "context_control_is_heterologous": "heterologous" in protocol.lower(),
        "implementation_conclusion": "The current protocol and implementation are internally inconsistent on the Null-1 p-value threshold; this must be resolved before treating the run as a strict protocol-compliant confirmatory result." if protocol_p01 and source_p005 else "No Null-1 threshold mismatch detected.",
    }


def _statistics_audit() -> dict:
    if not NULL1.exists() or not TIERS.exists():
        return {"status": "missing", "interpretation": "Required existing result tables are unavailable."}
    null1 = pd.read_csv(NULL1)
    tiers = pd.read_csv(TIERS)
    pcol = "empirical_two_sided_p"
    if pcol in null1.columns:
        null1["bh_q_value_global"] = _bh(null1[pcol])
        null1["bh_significant_global"] = null1["bh_q_value_global"] < 0.05
        q_values = null1[["feature", pcol, "bh_q_value_global", "bh_significant_global"]]
    else:
        q_values = pd.DataFrame(columns=["feature", pcol, "bh_q_value_global", "bh_significant_global"])
    merged = tiers.merge(q_values, on="feature", how="left")
    tier_a = merged[merged.get("tier", "") == "A"]
    tier_b = merged[merged.get("tier", "") == "B"]
    return {
        "status": "checked",
        "n_null1_tests": int(len(null1)),
        "n_null1_bh_significant_global": int(null1.get("bh_significant_global", pd.Series(dtype=bool)).sum()),
        "n_tier_a": int(len(tier_a)),
        "n_tier_b": int(len(tier_b)),
        "tier_b_is_deterministic_filter": True,
        "tier_b_inferential_p_value_present": False,
        "tier_b_multiple_testing_correction_applicable": False,
        "bh_scope": "all existing Null-1 candidate p-values in the saved 02_null1_temporal_order.csv",
        "interpretation": "Tier-B survival itself is a deterministic conjunction of direction/rank/confound filters, not an alpha-level statistical test. Multiple-testing correction is therefore applied here to the existing Null-1 p-values, not to the 1/28 survival count.",
        "candidate_q_values": q_values.to_dict(orient="records"),
    }


def _candidate_audit() -> dict:
    if not TIERS.exists():
        return {"status": "missing"}
    tiers = pd.read_csv(TIERS)
    b = tiers[tiers["tier"].astype(str) == "B"].copy()
    records = b.to_dict(orient="records")
    for r in records:
        r["independent_literature_validation"] = "not performed by this audit"
        r["candidate_status"] = "candidate_only_pending_external_validation"
    return {
        "status": "checked",
        "n_tier_b": int(len(b)),
        "surviving_candidates": records,
        "promotion_rule": "Do not promote a Tier-B candidate to a mechanistic claim from Stage 2.11C alone.",
    }


def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    mapping = _mapping_audit()
    provenance = _provenance_audit()
    implementation = _implementation_audit()
    statistics = _statistics_audit()
    candidate = _candidate_audit()

    historical_decision = None
    if SUMMARY.exists():
        try:
            historical_decision = json.loads(SUMMARY.read_text(encoding="utf-8")).get("decision")
        except Exception:
            historical_decision = None

    if historical_decision and historical_decision.get("decision") == "CONFOUNDED_UNSUPPORTED":
        interpreted_decision = "SPECIFICITY_UNSUPPORTED"
    else:
        interpreted_decision = historical_decision.get("decision") if isinstance(historical_decision, dict) else "INCONCLUSIVE"

    audit = {
        "status": "ok",
        "audit_scope": "existing Stage 2.11C outputs only; no full pipeline rerun",
        "historical_decision": historical_decision,
        "audited_decision_label": interpreted_decision,
        "confound_source": "unresolved",
        "evidence_status": "exploratory_until_protocol_mismatch_is_resolved" if implementation.get("null1_threshold_mismatch") else "exploratory",
        "mapping_integrity": mapping,
        "preregistration_provenance": provenance,
        "protocol_implementation_consistency": implementation,
        "statistical_multiplicity": statistics,
        "surviving_candidate_audit": candidate,
        "context_identifiability": {
            "status": "unresolved",
            "statement": "The GSE304042 comparison cannot isolate delivery effects from cell-type, perturbation (OSK vs OSKM), and broader experimental-design differences.",
            "not_supported_claims": [
                "Sendai is the identified confounder",
                "interferon/STAT is proven to cause the temporal signal",
                "the single Tier-B candidate is a true positive",
            ],
        },
        "recommended_next_gate": "Resolve the Null-1 protocol/implementation mismatch, lock the exact gene mapping table or annotation version, then reassess the saved candidate p-values with BH correction before any biological interpretation.",
    }

    (OUT / "09_stage2_11c_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    if statistics.get("candidate_q_values"):
        pd.DataFrame(statistics["candidate_q_values"]).to_csv(
            OUT / "09_stage2_11c_null1_bh_audit.csv", index=False
        )
    if candidate.get("surviving_candidates"):
        pd.DataFrame(candidate["surviving_candidates"]).to_csv(
            OUT / "09_stage2_11c_surviving_candidate_audit.csv", index=False
        )
    print(json.dumps(audit, indent=2, ensure_ascii=False, default=str))
    return audit


if __name__ == "__main__":
    run()
