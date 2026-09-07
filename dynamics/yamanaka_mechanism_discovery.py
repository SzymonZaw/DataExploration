"""Stage 2.11: perturbation-conditioned mechanism discovery.

Consumes the real Yamanaka POC inputs. The two GSE297234 RDS objects provide
0/3/7/10-day trajectories; GSE297233/GSE304042 provide independent intervention
signatures. The stage generates interpretable pathway/TF candidates and reports
whether perturbation-conditioned forecasting is identifiable.

No causal or lineage claim is made.
"""
from __future__ import annotations

import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from .yamanaka_poc import condition_table, load_expression, _signature
from .yamanaka_trajectory_poc import RDS_FILES, _convert_rds, _log_cpm

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11_mechanism_discovery"
PERTURBATION_FILES = {
    "GSE297233": DATA / "GSE297233_raw_counts_matrix.csv.gz",
    "GSE304042": DATA / "GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz",
}


def _activity_networks():
    try:
        import decoupler as dc
    except ImportError as exc:
        raise RuntimeError("Stage 2.11 requires decoupler. Install requirements-ai.txt.") from exc
    return (
        dc.op.progeny(organism="human", top=100),
        dc.op.dorothea(organism="human", levels=["A", "B", "C"]),
    )


def _score_time_series(X: pd.DataFrame, meta: pd.DataFrame, net: pd.DataFrame) -> pd.DataFrame:
    import decoupler as dc

    meta = meta.copy()
    meta["day"] = pd.to_numeric(meta["day"], errors="coerce")
    meta = meta.dropna(subset=["day"])
    meta["group"] = meta["group"].astype(str)
    groups = [g for g in meta.group.tolist() if g in X.columns]
    meta = meta[meta.group.isin(groups)].sort_values(["day", "group"])
    X = _log_cpm(X.loc[:, meta.group.tolist()])
    samples = X.T.copy()
    samples.index = meta.group.tolist()
    acts, _ = dc.mt.ulm(data=samples, net=net)
    acts = acts.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    acts["day"] = meta.set_index("group").loc[acts.index, "day"].astype(float)
    return acts


def _trajectory_activity_summary(scored: dict[str, pd.DataFrame]) -> dict:
    names = sorted(scored)
    if len(names) < 2:
        return {"status": "insufficient_datasets"}
    a, b = names[:2]
    A, B = scored[a], scored[b]
    features = sorted((set(A.columns) & set(B.columns)) - {"day"})
    rows = []
    for feature in features:
        aa = A.groupby("day")[feature].mean()
        bb = B.groupby("day")[feature].mean()
        common = sorted(set(aa.index) & set(bb.index))
        if len(common) < 3:
            continue
        va, vb = aa.loc[common].to_numpy(float), bb.loc[common].to_numpy(float)
        corr = np.corrcoef(va, vb)[0, 1] if np.std(va) > 0 and np.std(vb) > 0 else np.nan
        ta = pd.Series(common).corr(pd.Series(va), method="spearman")
        tb = pd.Series(common).corr(pd.Series(vb), method="spearman")
        rows.append({
            "feature": feature,
            "trajectory_correlation": float(corr) if np.isfinite(corr) else np.nan,
            "spearman_day_dataset_a": float(ta) if pd.notna(ta) else np.nan,
            "spearman_day_dataset_b": float(tb) if pd.notna(tb) else np.nan,
            "conserved_direction": bool(np.sign(va[-1] - va[0]) == np.sign(vb[-1] - vb[0])),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return {"status": "no_common_activity_features"}
    out["temporal_signal"] = (out.spearman_day_dataset_a.abs() + out.spearman_day_dataset_b.abs()) / 2
    out = out.sort_values(["conserved_direction", "temporal_signal", "trajectory_correlation"], ascending=[False, False, False])
    return {
        "status": "ok",
        "dataset_a": a,
        "dataset_b": b,
        "n_common_features": int(len(out)),
        "n_conserved_direction": int(out.conserved_direction.sum()),
        "median_trajectory_correlation": float(out.trajectory_correlation.median()),
        "top_candidates": out.head(20).to_dict("records"),
    }


def _perturbation_evidence() -> dict:
    exprs, labels = {}, {}
    for ds, path in PERTURBATION_FILES.items():
        X, _ = load_expression(path)
        labels[ds] = condition_table(X, ds).set_index("sample")["condition"]
        exprs[ds] = X

    signatures = {}
    for ds, X in exprs.items():
        lab = labels[ds]
        for condition in sorted(set(lab) - {"CONTROL", "OTHER"}):
            sig = _signature(X, lab, condition)
            if sig is not None:
                signatures[(ds, condition)] = sig

    rows = []
    for (ds, condition), sig in signatures.items():
        top = sig.abs().sort_values(ascending=False).head(100)
        rows.append({
            "dataset": ds,
            "perturbation": condition,
            "n_genes": int(sig.notna().sum()),
            "mean_abs_effect": float(sig.abs().mean()),
            "top100_mean_abs_effect": float(top.mean()),
        })

    osk = {ds: sig for (ds, condition), sig in signatures.items() if condition == "OSK"}
    agreement = None
    if len(osk) >= 2:
        names = sorted(osk)
        common = osk[names[0]].index.intersection(osk[names[1]].index)
        if len(common) >= 100:
            agreement = {
                "dataset_a": names[0],
                "dataset_b": names[1],
                "n_common_genes": int(len(common)),
                "pearson": float(osk[names[0]].loc[common].corr(osk[names[1]].loc[common])),
                "spearman": float(osk[names[0]].loc[common].corr(osk[names[1]].loc[common], method="spearman")),
            }
    return {"perturbation_signatures": rows, "osk_cross_dataset": agreement}


def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    progeny, dorothea = _activity_networks()
    scored = {"PROGENy": {}, "DoRothEA": {}}
    dataset_audits = []

    with tempfile.TemporaryDirectory(prefix="stage2_11_", dir=str(ROOT)) as td:
        tmp = Path(td)
        for ds, rds in RDS_FILES.items():
            X, meta = _convert_rds(ds, rds, tmp)
            dataset_audits.append({
                "dataset": ds,
                "n_genes": int(X.shape[0]),
                "n_groups": int(X.shape[1]),
                "days": sorted(pd.to_numeric(meta.day, errors="coerce").dropna().unique().tolist()),
            })
            scored["PROGENy"][ds] = _score_time_series(X, meta, progeny)
            scored["DoRothEA"][ds] = _score_time_series(X, meta, dorothea)

    activity_results = {}
    for representation, values in scored.items():
        feature_dir = OUT / representation
        feature_dir.mkdir(parents=True, exist_ok=True)
        for ds, frame in values.items():
            frame.to_csv(feature_dir / f"{ds}_activity.csv")
        activity_results[representation] = _trajectory_activity_summary(values)

    perturbation = _perturbation_evidence()
    result = {
        "status": "ok",
        "dataset_audit": dataset_audits,
        "activity_representations": activity_results,
        "independent_perturbation_evidence": perturbation,
        "identifiability": {
            "trajectory_datasets_have_multiple_perturbations": False,
            "reason": "Both GSE297234 trajectory objects represent OSKM reprogramming, so a perturbation coefficient cannot be identified from these time-resolved data alone.",
            "next_experiment": "Add timed datasets with distinct interventions and test whether candidate pathway/TF activities predict intervention-specific state changes."
        },
        "interpretation_guardrail": "Activity trajectories and perturbation signatures generate falsifiable mechanism candidates; they do not establish causality without intervention-specific temporal validation.",
    }
    (OUT / "stage2_11_mechanism_discovery_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=str), encoding="utf-8"
    )
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, ensure_ascii=False, default=str))
