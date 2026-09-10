"""Frozen donor-stratified external transfer audit for GSE297234.

Tests whether frozen GSE67462 temporal modules retain temporal direction and
shape independently in the young and aged human donors of GSE297234.
This is a transferability/replication diagnostic, not a Z4 acceptance test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DEFAULT_INPUT = "results/Dynamics/z4_gse297234_frozen_transfer_audit/04_frozen_module_scores_with_metadata.csv"
DEFAULT_OUT = "results/Dynamics/z4_gse297234_donor_transfer_audit"


def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan
    return float(spearmanr(x[m], y[m]).statistic)


def infer_donor(sample, age_group=None):
    s = str(sample).upper()
    if s.startswith("GM00731") or (isinstance(age_group, str) and "aged" in age_group.lower()):
        return "aged_GM00731"
    if s.startswith("GM23815") or (isinstance(age_group, str) and "young" in age_group.lower()):
        return "young_GM23815"
    return "unresolved"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=DEFAULT_INPUT)
    ap.add_argument("--output", default=DEFAULT_OUT)
    args = ap.parse_args()
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    if "age_group" in df.columns:
        df["donor"] = [infer_donor(s, a) for s, a in zip(df["sample"], df["age_group"])]
    else:
        df["donor"] = [infer_donor(s) for s in df["sample"]]

    df.to_csv(out / "01_sample_scores_with_donor.csv", index=False)

    rows = []
    for (donor, module), g in df.groupby(["donor", "module"]):
        g = g.dropna(subset=["inferred_time_days"]).sort_values("inferred_time_days")
        time = g["inferred_time_days"].to_numpy(float)
        score = g["score"].to_numpy(float)
        rho = spearman(time, score)
        rows.append({
            "donor": donor,
            "module": int(module),
            "n_genes": int(g["n_genes"].max()) if len(g) else 0,
            "n_samples": int(g["sample"].nunique()),
            "n_timepoints": int(g["inferred_time_days"].nunique()),
            "time_spearman": rho,
            "monotonic_increasing": bool(np.isfinite(rho) and rho > 0),
            "monotonic_decreasing": bool(np.isfinite(rho) and rho < 0),
        })
    summary = pd.DataFrame(rows).sort_values(["module", "donor"])
    summary.to_csv(out / "02_donor_module_transfer_summary.csv", index=False)

    # Compare donor-specific temporal direction for modules available in both donors.
    direction_rows = []
    for module, g in summary.groupby("module"):
        vals = dict(zip(g["donor"], g["time_spearman"]))
        a, y = vals.get("aged_GM00731", np.nan), vals.get("young_GM23815", np.nan)
        direction_rows.append({
            "module": int(module),
            "aged_time_spearman": a,
            "young_time_spearman": y,
            "same_direction": bool(np.isfinite(a) and np.isfinite(y) and np.sign(a) == np.sign(y) and a != 0 and y != 0),
            "both_positive": bool(np.isfinite(a) and np.isfinite(y) and a > 0 and y > 0),
        })
    direction = pd.DataFrame(direction_rows).sort_values("module")
    direction.to_csv(out / "03_donor_direction_concordance.csv", index=False)

    common = direction[direction["same_direction"]]
    both_pos = direction[direction["both_positive"]]
    result = {
        "candidate": "GSE297234",
        "input": str(Path(args.input)),
        "donors": ["aged_GM00731", "young_GM23815"],
        "modules_evaluated": int(summary["module"].nunique()),
        "modules_same_direction_across_donors": int(len(common)),
        "modules_positive_in_both_donors": int(len(both_pos)),
        "transfer_diagnostic": "DONOR_STRATIFIED_EXTERNAL_TEMPORAL_TRANSFER",
        "z4_D_replication": "DIAGNOSTIC_ONLY",
        "z4_A_external_target_association": "NOT_EVALUATED",
        "z4_B_target_vs_nontarget_specificity": "NOT_EVALUATED",
        "z4_C_context_robustness": "NOT_EVALUATED",
        "endpoint_independence": "NOT_ESTABLISHED",
        "decision": "Z4_UNRESOLVED",
        "interpretation": "DONOR_STRATIFIED_TRANSFER_DIAGNOSTIC_ONLY",
        "reason": "The two GSE297234 donors provide independent human biological samples, but donor agreement of frozen module trajectories remains a transfer/replication diagnostic because no independent biological endpoint has been established.",
        "frozen_z4_rule_unchanged": True,
    }
    (out / "04_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("GSE297234 Z4 DONOR-STRATIFIED TRANSFER AUDIT")
    print(f"modules evaluated: {result['modules_evaluated']}")
    print(f"same direction across donors: {result['modules_same_direction_across_donors']}")
    print(f"positive in both donors: {result['modules_positive_in_both_donors']}")
    print("endpoint independence: NOT ESTABLISHED")
    print("decision: Z4_UNRESOLVED")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
