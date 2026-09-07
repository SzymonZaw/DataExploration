"""Stage 2.11C: preregistered control and null analysis for Yamanaka signals.

This is a falsification gate, not a causal mechanism-discovery model.
It evaluates candidate pathway/TF activity trajectories with:
1) exact temporal-order permutations;
2) matched synthetic monotonic-shape controls;
3) heterologous OSK perturbation comparison (GSE297233 vs GSE304042);
4) candidate stability under sample/group bootstrap.

Large RDS inputs remain local and are never written to Git.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from .yamanaka_poc import condition_table, load_expression, _signature
from .yamanaka_trajectory_poc import RDS_FILES, _convert_rds, _log_cpm

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_control_null"
SEED = 20260911
N_MONOTONIC = 10_000
N_BOOTSTRAP = 1_000
EXPECTED_DAYS = (0.0, 3.0, 7.0, 10.0)
CONFOUNDER_PROGENY = {"JAK-STAT"}
CONFOUNDER_TF = {"IRF1", "IRF2", "IRF9", "STAT1", "STAT2"}

PERTURBATION_FILES = {
    "GSE297233": DATA / "GSE297233_raw_counts_matrix.csv.gz",
    "GSE304042": DATA / "GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz",
}


def _protocol_hash() -> str:
    path = ROOT / "dynamics" / "STAGE_2_11C_CONTROL_AND_NULL_PROTOCOL.md"
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else "missing"


def _activity_networks():
    try:
        import decoupler as dc
    except ImportError as exc:
        raise RuntimeError("Stage 2.11C requires decoupler. Install requirements-ai.txt.") from exc
    return dc.op.progeny(organism="human", top=100), dc.op.dorothea(organism="human", levels=["A", "B", "C"])


def _score_time_series(X: pd.DataFrame, meta: pd.DataFrame, net: pd.DataFrame) -> pd.DataFrame:
    import decoupler as dc

    meta = meta.copy()
    meta["day"] = pd.to_numeric(meta["day"], errors="coerce")
    meta["group"] = meta["group"].astype(str)
    meta = meta.dropna(subset=["day"])
    meta = meta[meta["day"].isin(EXPECTED_DAYS)]
    meta = meta[meta.group.isin(X.columns)].sort_values(["day", "group"])
    if meta.empty:
        return pd.DataFrame()
    samples = _log_cpm(X.loc[:, meta.group.tolist()]).T
    samples.index = meta.group.tolist()
    acts, _ = dc.mt.ulm(data=samples, net=net)
    acts = acts.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    acts = acts.fillna(0.0)
    acts["day"] = meta.set_index("group").loc[acts.index, "day"].astype(float)
    return acts


def _day_centroids(frame: pd.DataFrame) -> pd.DataFrame:
    features = [c for c in frame.columns if c != "day"]
    return frame.groupby("day")[features].mean().sort_index()


def _trajectory_rows(scored: dict[str, pd.DataFrame]) -> pd.DataFrame:
    names = sorted(scored)
    if len(names) < 2:
        return pd.DataFrame()
    A, B = scored[names[0]], scored[names[1]]
    ca, cb = _day_centroids(A), _day_centroids(B)
    common = sorted(set(ca.index) & set(cb.index))
    if len(common) != len(EXPECTED_DAYS):
        return pd.DataFrame()
    rows = []
    for feature in sorted(set(ca.columns) & set(cb.columns)):
        va = ca.loc[common, feature].to_numpy(float)
        vb = cb.loc[common, feature].to_numpy(float)
        if np.std(va) == 0 or np.std(vb) == 0:
            continue
        rho_a = pd.Series(common).corr(pd.Series(va), method="spearman")
        rho_b = pd.Series(common).corr(pd.Series(vb), method="spearman")
        rows.append({
            "feature": feature,
            "observed_corr": float(np.corrcoef(va, vb)[0, 1]),
            "rho_day_a": float(rho_a),
            "rho_day_b": float(rho_b),
            "conserved_direction": bool(np.sign(va[-1] - va[0]) == np.sign(vb[-1] - vb[0]) and va[-1] != va[0] and vb[-1] != vb[0]),
            "abs_temporal_signal": float((abs(rho_a) + abs(rho_b)) / 2),
            "a_values": va.tolist(),
            "b_values": vb.tolist(),
        })
    return pd.DataFrame(rows)


def _temporal_null(rows: pd.DataFrame) -> pd.DataFrame:
    perms = list(itertools.permutations(range(len(EXPECTED_DAYS))))
    out = []
    for r in rows.itertuples(index=False):
        va, vb = np.asarray(r.a_values, float), np.asarray(r.b_values, float)
        null = np.array([np.corrcoef(va, vb[list(p)])[0, 1] for p in perms], dtype=float)
        obs = float(r.observed_corr)
        ge = int(np.sum(np.abs(null) >= abs(obs) - 1e-12))
        p = (ge + 1) / (len(null) + 1)
        out.append({
            "feature": r.feature,
            "n_exact_permutations": len(perms),
            "null_q95": float(np.quantile(null, 0.95)),
            "null_q99": float(np.quantile(null, 0.99)),
            "empirical_two_sided_p": float(p),
            "passes_null1": bool(abs(obs) > np.quantile(np.abs(null), 0.99) and p <= 0.05),
        })
    return pd.DataFrame(out)


def _monotonic_series(rng: np.random.Generator, start: float, end: float) -> np.ndarray:
    if start < end:
        interior = np.sort(rng.uniform(start, end, size=2))
    else:
        interior = np.sort(rng.uniform(end, start, size=2))[::-1]
    return np.array([start, interior[0], interior[1], end], dtype=float)


def _monotonic_null(rows: pd.DataFrame, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    out = []
    for r in rows.itertuples(index=False):
        va, vb = np.asarray(r.a_values, float), np.asarray(r.b_values, float)
        null = np.empty(N_MONOTONIC, dtype=float)
        for i in range(N_MONOTONIC):
            a = _monotonic_series(rng, va[0], va[-1])
            b = _monotonic_series(rng, vb[0], vb[-1])
            null[i] = np.corrcoef(a, b)[0, 1]
        obs = float(r.observed_corr)
        out.append({
            "feature": r.feature,
            "n_monotonic_controls": N_MONOTONIC,
            "control_q95": float(np.quantile(null, 0.95)),
            "control_q99": float(np.quantile(null, 0.99)),
            "passes_null2": bool(obs > np.quantile(null, 0.95)),
        })
    return pd.DataFrame(out)


def _perturbation_activity(network: pd.DataFrame, condition: str) -> dict[str, pd.Series]:
    result = {}
    for ds, path in PERTURBATION_FILES.items():
        X, _ = load_expression(path)
        labels = condition_table(X, ds).set_index("sample")["condition"]
        sig = _signature(X, labels, condition)
        if sig is None:
            continue
        import decoupler as dc
        sample = pd.DataFrame([sig.to_numpy(float)], columns=sig.index, index=[f"{ds}__{condition}"])
        acts, _ = dc.mt.ulm(data=sample, net=network)
        result[ds] = acts.iloc[0].apply(float)
    return result


def _context_control(progeny: pd.DataFrame, dorothea: pd.DataFrame, candidate_features: dict[str, list[str]]) -> tuple[pd.DataFrame, dict]:
    p = _perturbation_activity(progeny, "OSK")
    t = _perturbation_activity(dorothea, "OSK")
    rows = []
    detail = {}
    for representation, acts in (("PROGENy", p), ("DoRothEA", t)):
        features = candidate_features.get(representation, [])
        for feature in features:
            if len(acts) < 2 or any(feature not in a.index for a in acts.values()):
                rows.append({"representation": representation, "feature": feature, "status": "inconclusive", "survives": False})
                continue
            names = sorted(acts)
            a, b = acts[names[0]][feature], acts[names[1]][feature]
            direction = np.sign(a) == np.sign(b) and a != 0 and b != 0
            rank_lists = [s.abs().sort_values(ascending=False).index.tolist() for s in acts.values()]
            ranks = [rl.index(feature) + 1 for rl in rank_lists]
            percentile = max(ranks) / max(len(rank_lists[0]), 1)
            survives = bool(direction and percentile <= 0.25)
            confound = bool((representation == "PROGENy" and feature in CONFOUNDER_PROGENY) or (representation == "DoRothEA" and feature in CONFOUNDER_TF))
            rows.append({"representation": representation, "feature": feature, "activity_a": float(a), "activity_b": float(b), "concordant_direction": bool(direction), "worst_rank_fraction": float(percentile), "known_delivery_confounded": confound, "survives": survives and not confound, "status": "ok"})
        detail[representation] = sorted(acts)
    return pd.DataFrame(rows), detail


def _bootstrap_stability(scored: dict[str, pd.DataFrame], candidate_features: dict[str, list[str]], seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    for representation, features in candidate_features.items():
        for feature in features:
            same = 0
            total = 0
            for _ in range(N_BOOTSTRAP):
                vals = []
                for ds in sorted(scored):
                    f = scored[ds]
                    sub = f[f["day"].isin(EXPECTED_DAYS)]
                    if feature not in sub.columns:
                        continue
                    days = sorted(sub.day.unique())
                    if len(days) != 4:
                        continue
                    # Resample groups within each day when replicates exist.
                    parts = []
                    for day in days:
                        g = sub[sub.day == day][feature].to_numpy(float)
                        if len(g):
                            parts.append(float(rng.choice(g)))
                    if len(parts) == 4:
                        vals.append(parts)
                if len(vals) == 2:
                    a, b = map(np.asarray, vals)
                    if np.std(a) and np.std(b):
                        corr = np.corrcoef(a, b)[0, 1]
                        same += int(np.sign(a[-1] - a[0]) == np.sign(b[-1] - b[0]))
                        total += 1
            rows.append({"representation": representation, "feature": feature, "bootstrap_replicates": N_BOOTSTRAP, "stable_direction_fraction": float(same / total) if total else np.nan, "stable": bool(total and same / total >= 0.8)})
    return pd.DataFrame(rows)


def _assign_tiers(detail: pd.DataFrame, context: pd.DataFrame, stability: pd.DataFrame) -> pd.DataFrame:
    out = detail.merge(context[["representation", "feature", "survives", "known_delivery_confounded"]], on=["representation", "feature"], how="left")
    out = out.merge(stability[["representation", "feature", "stable"]], on=["representation", "feature"], how="left")
    out["tier"] = "NONE"
    tier_a = out.passes_null1 & out.passes_null2 & out.conserved_direction & out.stable
    out.loc[tier_a, "tier"] = "A"
    tier_b = tier_a & out.survives.fillna(False) & ~out.known_delivery_confounded.fillna(False)
    out.loc[tier_b, "tier"] = "B"
    return out


def _decision(tiered: pd.DataFrame) -> dict:
    a = tiered[tiered.tier == "A"]
    if a.empty:
        return {"decision": "CONFOUNDED_UNSUPPORTED", "tier_a_n": 0, "tier_b_n": 0, "survival_rate": 0.0, "reason": "No candidate passed the temporal-order, monotonic-shape and stability gates."}
    b = tiered[tiered.tier == "B"]
    rate = len(b) / len(a)
    if rate >= 0.60:
        decision = "PROCEED"
    elif rate >= 0.30:
        decision = "MIXED"
    else:
        decision = "CONFOUNDED_UNSUPPORTED"
    return {"decision": decision, "tier_a_n": int(len(a)), "tier_b_n": int(len(b)), "survival_rate": float(rate), "reason": "Decision follows the fixed Stage 2.11C thresholds."}


def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    progeny, dorothea = _activity_networks()
    scored = {"PROGENy": {}, "DoRothEA": {}}
    audits = []
    with tempfile.TemporaryDirectory(prefix="stage2_11c_", dir=str(ROOT)) as td:
        tmp = Path(td)
        for ds, rds in RDS_FILES.items():
            X, meta = _convert_rds(ds, rds, tmp)
            audits.append({"dataset": ds, "n_genes": int(X.shape[0]), "n_groups": int(X.shape[1]), "days": sorted(pd.to_numeric(meta.day, errors="coerce").dropna().unique().tolist())})
            scored["PROGENy"][ds] = _score_time_series(X, meta, progeny)
            scored["DoRothEA"][ds] = _score_time_series(X, meta, dorothea)

    candidate_features = {}
    trajectory_frames = []
    for rep, values in scored.items():
        rows = _trajectory_rows(values)
        if not rows.empty:
            rows["representation"] = rep
            trajectory_frames.append(rows)
            candidate_features[rep] = rows.feature.tolist()
    trajectory = pd.concat(trajectory_frames, ignore_index=True) if trajectory_frames else pd.DataFrame()
    null1 = _temporal_null(trajectory) if not trajectory.empty else pd.DataFrame()
    null2 = _monotonic_null(trajectory, SEED + 1) if not trajectory.empty else pd.DataFrame()
    detail = trajectory.merge(null1, on="feature", how="left").merge(null2, on="feature", how="left")
    context, context_detail = _context_control(progeny, dorothea, candidate_features)
    stability = _bootstrap_stability(scored["PROGENy"], {"PROGENy": candidate_features.get("PROGENy", [])}, SEED + 2)
    stability_tf = _bootstrap_stability(scored["DoRothEA"], {"DoRothEA": candidate_features.get("DoRothEA", [])}, SEED + 3)
    stability = pd.concat([stability, stability_tf], ignore_index=True)
    if not detail.empty:
        tiered = _assign_tiers(detail, context, stability)
    else:
        tiered = pd.DataFrame()
    decision = _decision(tiered) if not tiered.empty else {"decision": "INCONCLUSIVE", "reason": "No complete two-dataset temporal activity matrix was available."}

    OUT.mkdir(parents=True, exist_ok=True)
    trajectory.to_csv(OUT / "01_candidate_trajectories.csv", index=False)
    null1.to_csv(OUT / "02_null1_temporal_order.csv", index=False)
    null2.to_csv(OUT / "03_null2_monotonic_shape.csv", index=False)
    context.to_csv(OUT / "04_context_control_OSK.csv", index=False)
    stability.to_csv(OUT / "05_candidate_stability.csv", index=False)
    tiered.to_csv(OUT / "06_candidate_tiers.csv", index=False)

    summary = {
        "status": "ok",
        "protocol_hash": _protocol_hash(),
        "seed": SEED,
        "n_monotonic_controls": N_MONOTONIC,
        "n_bootstrap": N_BOOTSTRAP,
        "dataset_audit": audits,
        "context_control_datasets": context_detail,
        "decision": decision,
        "interpretation_guardrail": "Stage 2.11C is a falsification gate. Passing supports candidate temporal regulators; it does not establish causality, lineage, or a molecular mechanism.",
        "important_limitations": [
            "The temporal-order permutation has only 24 exact permutations for a four-day trajectory, so the minimum attainable +1-corrected empirical p-value is 1/25 = 0.04.",
            "GSE304042 is a heterologous biological/context control, not a clean Sendai-only control because cell type and experimental design differ.",
            "Pseudobulk cannot establish single-cell lineage or cell-level transition dynamics.",
        ],
    }
    (OUT / "stage2_11c_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=str), encoding="utf-8")
    return summary


def main() -> None:
    print(json.dumps(run(), indent=2, ensure_ascii=False, default=str))


if __name__ == "__main__":
    main()
