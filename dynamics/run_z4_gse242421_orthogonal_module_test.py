"""Frozen-module orthogonal test for GSE242421 scATAC gene activity.

No module discovery or target-data fitting is performed. The script maps the
frozen GSE67462 mouse module genes to human genes with the frozen 1:1 mapping,
computes target gene-activity module scores, and evaluates their temporal
association against a time-permutation null.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DEFAULT_ACTIVITY = Path("results/Dynamics/z4_gse242421_gene_activity/01_gene_activity_log1p_cpm.tsv.gz")
DEFAULT_MAPPING = Path("results/Dynamics/z4_frozen_orthology_mapping/01_frozen_human_mouse_mapping.csv")
DEFAULT_MODULES = Path("results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv")
DEFAULT_OUT = Path("results/Dynamics/z4_gse242421_orthogonal_module_test")
TIME_ORDER = ["D0", "D2", "D4", "D6", "D8", "D10", "D12", "D14"]
ALL_ORDER = TIME_ORDER + ["iPSC"]


def find_reference_direction(explicit: Path | None, root: Path) -> tuple[pd.DataFrame | None, str | None]:
    if explicit is not None and explicit.exists():
        return pd.read_csv(explicit), str(explicit)
    candidates = [
        root / "results/Dynamics/z6_gse67462_temporal_modules/04_module_trajectories.csv",
        root / "results/Dynamics/z6_gse67462_temporal_modules/05_module_summary.csv",
        root / "results/Dynamics/z6_gse67462_temporal_modules/module_trajectories.csv",
    ]
    for p in candidates:
        if p.exists():
            return pd.read_csv(p), str(p)
    return None, None


def normalize_gene(x: object) -> str:
    return str(x).strip().upper()


def sample_day(sample: object) -> float:
    m = re.search(r"D(\d+(?:\.\d+)?)", str(sample), flags=re.I)
    return float(m.group(1)) if m else np.nan


def rho(x: np.ndarray, y: np.ndarray) -> float:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return np.nan
    return float(spearmanr(x[m], y[m]).statistic)


def direction_from_reference(ref: pd.DataFrame, module: int) -> float | None:
    if ref is None or ref.empty:
        return None
    modcol = next((c for c in ref.columns if c.lower() == "module"), None)
    if modcol is None:
        return None
    g = ref[pd.to_numeric(ref[modcol], errors="coerce") == module]
    if g.empty:
        return None
    for c in ["reference_direction", "frozen_direction", "time_spearman", "gse67462_time_spearman", "direction"]:
        if c in g.columns:
            vals = pd.to_numeric(g[c], errors="coerce").dropna()
            if len(vals):
                v = float(vals.iloc[0])
                if v > 0: return 1.0
                if v < 0: return -1.0
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--activity", type=Path, default=DEFAULT_ACTIVITY)
    ap.add_argument("--mapping", type=Path, default=DEFAULT_MAPPING)
    ap.add_argument("--modules", type=Path, default=DEFAULT_MODULES)
    ap.add_argument("--reference", type=Path, default=None)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--permutations", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=424242)
    args = ap.parse_args()
    out = args.output; out.mkdir(parents=True, exist_ok=True)
    root = Path.cwd()

    for p in [args.activity, args.mapping, args.modules]:
        if not p.exists():
            raise FileNotFoundError(f"Required input missing: {p}")

    activity = pd.read_csv(args.activity, sep="\t", compression="gzip", index_col=0)
    activity.index = activity.index.map(normalize_gene)
    activity = activity[~activity.index.duplicated(keep="first")]

    mapping = pd.read_csv(args.mapping)
    required_mapping = {"human_gene", "mouse_gene"}
    if not required_mapping.issubset(mapping.columns):
        raise ValueError(f"Frozen mapping missing {required_mapping - set(mapping.columns)}")
    mapping["human_gene_norm"] = mapping["human_gene"].map(normalize_gene)
    mapping["mouse_gene_norm"] = mapping["mouse_gene"].map(normalize_gene)
    mapping = mapping.dropna(subset=["human_gene_norm", "mouse_gene_norm"]).drop_duplicates("human_gene_norm")
    mouse_to_human = dict(zip(mapping["mouse_gene_norm"], mapping["human_gene_norm"]))

    modules = pd.read_csv(args.modules)
    if not {"module", "gene"}.issubset(modules.columns):
        raise ValueError("Frozen module table must contain module and gene columns")
    modules["gene_norm"] = modules["gene"].map(normalize_gene)

    ref, ref_path = find_reference_direction(args.reference, root)
    rng = np.random.default_rng(args.seed)
    sample_cols = [c for c in ALL_ORDER if c in activity.columns]
    if sample_cols != ALL_ORDER:
        raise ValueError(f"Gene-activity matrix missing samples: {sorted(set(ALL_ORDER)-set(sample_cols))}")
    time_days = np.array([sample_day(c) for c in TIME_ORDER], dtype=float)

    score_rows = []
    module_rows = []
    for module in sorted(pd.to_numeric(modules["module"], errors="coerce").dropna().astype(int).unique()):
        mg = modules.loc[pd.to_numeric(modules["module"], errors="coerce") == module, "gene_norm"].unique()
        human = sorted({mouse_to_human[g] for g in mg if g in mouse_to_human and mouse_to_human[g] in activity.index})
        if not human:
            continue
        score = activity.loc[human, ALL_ORDER].mean(axis=0).to_numpy(float)
        for s, v in zip(ALL_ORDER, score):
            score_rows.append({"module": module, "sample": s, "time_days": sample_day(s), "score": float(v), "n_genes": len(human)})
        obs = rho(time_days, score[:len(TIME_ORDER)])
        null = np.empty(args.permutations, dtype=float)
        for i in range(args.permutations):
            null[i] = rho(time_days, rng.permutation(score[:len(TIME_ORDER)]))
        ref_dir = direction_from_reference(ref, module)
        if ref_dir is None:
            p = np.nan
            concordant = np.nan
        else:
            signed_obs = ref_dir * obs
            signed_null = ref_dir * null
            p = float((np.sum(signed_null >= signed_obs) + 1) / (len(signed_null) + 1))
            concordant = bool(np.isfinite(obs) and np.sign(obs) == np.sign(ref_dir) and obs != 0)
        module_rows.append({
            "module": module,
            "frozen_mouse_genes": len(mg),
            "validated_human_genes": len(human),
            "coverage_fraction": len(human) / max(len(mg), 1),
            "target_time_spearman": obs,
            "frozen_reference_direction": ref_dir,
            "direction_concordant": concordant,
            "permutation_p_directional": p,
            "permutations": args.permutations,
        })

    scores = pd.DataFrame(score_rows)
    summary = pd.DataFrame(module_rows).sort_values("module")
    scores.to_csv(out / "01_orthogonal_module_scores.csv", index=False)
    summary.to_csv(out / "02_orthogonal_module_summary.csv", index=False)

    evaluable = summary[summary["frozen_reference_direction"].notna()].copy()
    valid_cov = bool(len(summary) > 0 and (summary["validated_human_genes"] >= 3).all())
    if evaluable.empty or not valid_cov:
        decision = "Z4_ORTHO_UNRESOLVED"
    else:
        passed = (evaluable["direction_concordant"] & (evaluable["permutation_p_directional"] < 0.05))
        n_pass = int(passed.sum())
        if n_pass == len(evaluable): decision = "Z4_ORTHO_SUPPORTED"
        elif n_pass > 0: decision = "Z4_ORTHO_PARTIAL"
        else: decision = "Z4_ORTHO_FAILED"

    result = {
        "candidate": "GSE242421",
        "activity_input": str(args.activity),
        "frozen_mapping": str(args.mapping),
        "frozen_modules": str(args.modules),
        "reference_direction_file": ref_path,
        "modules_evaluated": int(len(summary)),
        "modules_with_frozen_direction": int(len(evaluable)),
        "modules_passing_directional_null": int(((evaluable.get("direction_concordant", pd.Series(dtype=bool))) & (evaluable.get("permutation_p_directional", pd.Series(dtype=float)) < 0.05)).sum()) if not evaluable.empty else 0,
        "target_samples": ALL_ORDER,
        "continuous_test_samples": TIME_ORDER,
        "endpoint_iPSC_reported_separately": True,
        "permutation_n": args.permutations,
        "seed": args.seed,
        "module_fitting_on_target": False,
        "feature_selection_on_target": False,
        "frozen_rule_unchanged": True,
        "decision": decision,
        "interpretation": "ORTHOGONAL_MOLECULAR_CORROBORATION_ONLY",
    }
    (out / "03_summary.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("GSE242421 Z4 ORTHOGONAL FROZEN-MODULE TEST")
    print(f"modules evaluated: {len(summary)}")
    print(f"modules with frozen direction: {len(evaluable)}")
    print(f"modules passing directional permutation test: {result['modules_passing_directional_null']}")
    print(f"reference direction file: {ref_path or 'NOT FOUND'}")
    print(f"decision: {decision}")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
