"""Stage 2.11C: direct attribution audit for process-specific state information.

This module audits saved outputs only. It deliberately does not rerun the expensive
Yamanaka pipeline and reports unavailable design blocks as unresolved.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_control_null"
PROTOCOL = ROOT / "dynamics" / "STAGE_2_11C_SIGNAL_ATTRIBUTION_PROTOCOL.md"
SEED = 20260911
N_PERMUTATIONS = 10_000


def _read_csv(name: str) -> pd.DataFrame:
    path = OUT / name
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def _protocol_hash() -> str:
    return hashlib.sha256(PROTOCOL.read_bytes()).hexdigest() if PROTOCOL.exists() else "missing"


def _parse_array(value: object) -> np.ndarray:
    if isinstance(value, (list, tuple, np.ndarray)):
        return np.asarray(value, dtype=float)
    return np.asarray(json.loads(str(value)), dtype=float)


def _feature_family_counts(rows: pd.DataFrame) -> dict:
    if rows.empty or "representation" not in rows.columns:
        return {}
    return {str(k): int(v) for k, v in rows["representation"].value_counts().sort_index().items()}


def _trajectory_source(null1: pd.DataFrame, tiers: pd.DataFrame) -> pd.DataFrame:
    """Return the saved table that actually contains observed correlations and trajectories.

    02_null1_temporal_order.csv contains the null-test results only; the observed
    trajectory vectors live in the candidate trajectory/tier table. Earlier code
    incorrectly assumed both schemas were identical, causing KeyError: observed_corr.
    """
    required = {"observed_corr", "abs_temporal_signal", "conserved_direction", "a_values", "b_values"}
    if required.issubset(tiers.columns):
        return tiers
    if required.issubset(null1.columns):
        return null1
    return pd.DataFrame()


def _trajectory_signal_audit(rows: pd.DataFrame) -> dict:
    if rows.empty:
        return {"status": "missing", "n_features": 0}
    observed = pd.to_numeric(rows["observed_corr"], errors="coerce").dropna().to_numpy(float)
    temporal = pd.to_numeric(rows["abs_temporal_signal"], errors="coerce").dropna().to_numpy(float)
    conserved = rows["conserved_direction"].fillna(False).astype(bool)
    return {
        "status": "descriptive_only",
        "n_features": int(len(rows)),
        "n_features_corr_ge_0_9": int(np.sum(observed >= 0.9)),
        "median_feature_correlation": float(np.median(observed)) if len(observed) else np.nan,
        "mean_feature_correlation": float(np.mean(observed)) if len(observed) else np.nan,
        "median_abs_temporal_signal": float(np.median(temporal)) if len(temporal) else np.nan,
        "conserved_direction_fraction": float(conserved.mean()) if len(conserved) else np.nan,
        "warning": "These quantities measure reproducible trajectory structure, not process specificity.",
    }


def _whole_representation_permutation(rows: pd.DataFrame) -> dict:
    if rows.empty or not {"a_values", "b_values"}.issubset(rows.columns):
        return {"status": "missing"}
    pairs = []
    for r in rows.itertuples(index=False):
        try:
            a, b = _parse_array(r.a_values), _parse_array(r.b_values)
        except (ValueError, TypeError, json.JSONDecodeError):
            continue
        if len(a) == len(b) and len(a) >= 3 and np.std(a) > 0 and np.std(b) > 0:
            pairs.append((a, b))
    if len(pairs) < 10:
        return {"status": "insufficient_features", "n_features": len(pairs)}
    A = np.vstack([a for a, _ in pairs])
    B = np.vstack([b for _, b in pairs])
    def z(x):
        s = x.std(axis=1, keepdims=True)
        return (x - x.mean(axis=1, keepdims=True)) / np.where(s == 0, 1, s)
    ZA, ZB = z(A), z(B)
    observed = float(np.mean(np.sum(ZA * ZB, axis=1) / ZA.shape[1]))
    rng = np.random.default_rng(SEED)
    null = np.empty(N_PERMUTATIONS, dtype=float)
    for i in range(N_PERMUTATIONS):
        perm = rng.permutation(len(pairs))
        null[i] = float(np.mean(np.sum(ZA * ZB[perm], axis=1) / ZA.shape[1]))
    ge = int(np.sum(null >= observed - 1e-12))
    return {
        "status": "secondary_diagnostic",
        "n_features": int(len(pairs)),
        "observed_mean_feature_similarity": observed,
        "null_median": float(np.median(null)),
        "null_q95": float(np.quantile(null, 0.95)),
        "empirical_one_sided_p": float((ge + 1) / (N_PERMUTATIONS + 1)),
        "interpretation": "Feature identity contributes to cross-dataset trajectory alignment under this diagnostic; this is not proof of process specificity.",
    }


def _context_audit(context: pd.DataFrame) -> dict:
    if context.empty:
        return {"status": "missing"}
    if "known_delivery_confounded" not in context.columns:
        return {"status": "incomplete"}
    return {
        "status": "heterologous_only",
        "n_rows": int(len(context)),
        "n_survives": int(context.get("survives", pd.Series(False, index=context.index)).fillna(False).astype(bool).sum()),
        "n_known_delivery_confounded": int(context["known_delivery_confounded"].fillna(False).astype(bool).sum()),
        "limitation": "GSE304042 changes multiple biological/experimental factors and cannot identify a unique confounder.",
    }


def _negative_control_audit() -> dict:
    manifest = OUT / "signal_attribution_negative_controls.json"
    if not manifest.exists():
        return {"status": "missing", "required": True, "reason": "No pre-registered unrelated temporal negative-control manifest is present."}
    try:
        spec = json.loads(manifest.read_text(encoding="utf-8"))
    except Exception as exc:
        return {"status": "invalid_manifest", "error": type(exc).__name__}
    datasets = spec.get("datasets", [])
    if len(datasets) < 2:
        return {"status": "insufficient", "n_declared": len(datasets), "required_minimum": 2}
    missing = [str(x.get("path", "")) for x in datasets if not (ROOT / str(x.get("path", ""))).exists()]
    return {"status": "ready" if not missing else "missing_files", "n_declared": len(datasets), "missing_files": missing, "datasets": datasets}


def _decision(blocks: dict) -> tuple[str, str]:
    neg = blocks["unrelated_temporal_controls"]
    process = blocks["process_discrimination"]
    nuisance = blocks["nuisance_residualization"]
    context = blocks["matched_context_controls"]
    if neg["status"] not in {"ready", "analyzed"}:
        return "SIGNAL_UNRESOLVED", "Unrelated temporal negative controls are not available; H3 cannot be tested."
    if process["status"] != "analyzed":
        return "SIGNAL_UNRESOLVED", "Sample-level process discrimination was not analyzed from a locked independent validation design."
    if nuisance["status"] != "analyzed":
        return "SIGNAL_UNRESOLVED", "Nuisance attribution/residualization was not analyzed."
    if context["status"] == "dominant":
        return "CONFOUNDER_DOMINANT", "A measurable nuisance/context axis dominates the reproducible signal."
    if neg.get("artifact_plausible"):
        return "AGGREGATION_ARTIFACT_PLAUSIBLE", "Unrelated temporal controls produce comparable representation agreement."
    if process.get("passes") and nuisance.get("passes") and neg.get("passes"):
        return "SPECIFIC_SIGNAL_SUPPORTED", "All predeclared specificity blocks passed."
    return "SIGNAL_UNRESOLVED", "The available evidence does not distinguish process-specific signal from nuisance or generic temporal structure."


def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    null1 = _read_csv("02_null1_temporal_order.csv")
    context = _read_csv("04_context_control_OSK.csv")
    tiers = _read_csv("06_candidate_tiers.csv")
    trajectory_rows = _trajectory_source(null1, tiers)

    process_discrimination = {"status": "not_analyzed", "reason": "Saved Stage 2.11C outputs contain feature trajectories, not independent sample-level state vectors/outcomes required for cross-validated process discrimination."}
    nuisance = {"status": "not_analyzed", "reason": "No locked sample-level nuisance matrix is persisted in the historical output."}
    negative = _negative_control_audit()
    blocks = {
        "process_discrimination": process_discrimination,
        "matched_context_controls": _context_audit(context),
        "unrelated_temporal_controls": negative,
        "nuisance_residualization": nuisance,
    }
    decision, reason = _decision(blocks)
    trajectory = _trajectory_signal_audit(trajectory_rows)
    alignment = _whole_representation_permutation(trajectory_rows)
    summary = {
        "status": "ok",
        "audit_scope": "saved Stage 2.11C outputs only; no expensive pipeline rerun",
        "protocol_hash": _protocol_hash(),
        "seed": SEED,
        "feature_family_counts": _feature_family_counts(tiers if not tiers.empty else null1),
        "historical_trajectory_structure": trajectory,
        "whole_representation_alignment_diagnostic": alignment,
        "blocks": blocks,
        "decision": {"primary": decision, "reason": reason},
        "guardrail": "High reproducibility of S(t) is not interpreted as process specificity until independent process, matched nuisance, unrelated temporal, and residualization blocks are available.",
    }
    (OUT / "10_signal_attribution_audit.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    return summary


if __name__ == "__main__":
    run()
