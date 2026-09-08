"""Prospective H3 similarity/null analysis for GSE263713.

This stage is deliberately isolated from the historical Stage 2.11C run.
It uses GSE263713 as an unrelated circadian temporal-control candidate and
compares its feature trajectories against the saved Yamanaka feature
trajectories. No causal or mechanistic claim is made here.
"""
from __future__ import annotations

import ast
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
H3 = DATA / "GSE263713_raw_counts.tsv.gz"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
TARGET = ROOT / "results" / "Dynamics" / "stage2_11c_control_null"
EXPECTED_SHA = "8ce1e16a0039d93449e70aa6e827bbc8510a9b15dfbb78dc920cc2a2de5615a6"
SEED = 20260911
N_PERMUTATIONS = 10000
YAMANAKA_DAYS = (0.0, 3.0, 7.0, 10.0)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_h3() -> pd.DataFrame:
    """Read featureCounts-style matrix and return genes x samples counts."""
    df = pd.read_csv(H3, sep="\t", compression="gzip", comment="#", low_memory=False)
    sample_cols = [c for c in df.columns if re.fullmatch(r"[^-]+-TP-\d+(?:\.\d+)?", str(c))]
    if len(sample_cols) != 78:
        raise ValueError(f"Expected 78 H3 sample columns, found {len(sample_cols)}")
    genes = df.iloc[:, 0].astype(str).str.strip()
    x = df.loc[:, sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    x.index = genes
    x = x.groupby(level=0).sum()
    return x


def map_hgnc(expr: pd.DataFrame) -> pd.DataFrame:
    """Map Ensembl IDs if present; retain HGNC-like symbols directly."""
    ids = pd.Index(expr.index.astype(str))
    ensembl_fraction = np.mean([bool(re.fullmatch(r"ENSG\d+(?:\.\d+)?", x.upper())) for x in ids])
    if ensembl_fraction < 0.8:
        return expr
    import mygene
    clean = sorted(set(x.split(".", 1)[0].upper() for x in ids))
    mapping = {}
    info = mygene.MyGeneInfo()
    for i in range(0, len(clean), 1000):
        result = info.querymany(clean[i:i + 1000], scopes="ensembl.gene", fields="symbol", species="human", as_dataframe=False, returnall=False, step=1000, verbose=False)
        for row in result:
            if row.get("notfound"):
                continue
            q, symbol = row.get("query"), row.get("symbol")
            if q and symbol:
                mapping[str(q).upper()] = str(symbol).upper()
    symbols = [mapping.get(x.split(".", 1)[0].upper(), "") for x in ids]
    keep = np.array([bool(x) for x in symbols])
    out = expr.loc[keep].copy()
    out.index = np.asarray(symbols)[keep]
    return out.groupby(level=0).mean()


def activity_networks():
    import decoupler as dc
    return {
        "PROGENy": dc.op.progeny(organism="human", top=100),
        "DoRothEA": dc.op.dorothea(organism="human", levels=["A", "B", "C"]),
    }


def score_activity(expr: pd.DataFrame, net: pd.DataFrame, kind: str) -> pd.DataFrame:
    import decoupler as dc
    common = expr.index.intersection(net["target"].astype(str).str.upper())
    if len(common) < 50:
        raise ValueError(f"Too few shared genes for {kind}: {len(common)}")
    x = expr.loc[common].T
    n = net.copy()
    n["target"] = n["target"].astype(str).str.upper()
    n = n[n.target.isin(common)]
    result = dc.mt.ulm(data=x, net=n, tmin=5)
    if isinstance(result, tuple):
        result = result[0]
    if "source" in result.columns and "score" in result.columns:
        return result.pivot(index=result.index, columns="source", values="score")
    return result


def parse_sample(sample: str):
    m = re.fullmatch(r"([^-]+)-TP-(\d+(?:\.\d+)?)", str(sample))
    if not m:
        raise ValueError(f"Unexpected H3 sample name: {sample}")
    return m.group(1), float(m.group(2))


def trajectory_matrix(activity: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample in activity.index:
        individual, time = parse_sample(sample)
        for feature, value in activity.loc[sample].items():
            rows.append({"individual": individual, "time": time, "feature": feature, "value": float(value)})
    long = pd.DataFrame(rows)
    return long.groupby(["feature", "time"], as_index=False)["value"].mean().pivot(index="feature", columns="time", values="value")


def _parse_value_array(raw, column: str, feature: str) -> np.ndarray:
    """Parse one serialized four-point trajectory without concatenating variants."""
    if isinstance(raw, (list, tuple, np.ndarray)):
        values = raw
    elif pd.isna(raw):
        raise ValueError(f"Missing {column} for {feature}")
    else:
        text = str(raw).strip()
        try:
            values = json.loads(text)
        except json.JSONDecodeError:
            try:
                values = ast.literal_eval(text)
            except (ValueError, SyntaxError) as exc:
                raise ValueError(f"Cannot parse {column} for {feature}: {raw!r}") from exc
    values = np.asarray(values, dtype=float).reshape(-1)
    if len(values) != 4:
        raise ValueError(f"Unexpected {column} length for {feature}: {len(values)}; expected 4")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"Non-finite value in {column} for {feature}")
    return values


def target_traj() -> tuple[pd.DataFrame, dict]:
    """Load the saved Yamanaka trajectories using one canonical 4-point variant.

    `a_values` and `b_values` are alternative encodings of the SAME trajectory,
    not two temporal series to concatenate. Both must contain exactly the four
    Yamanaka time points. `a_values` is the locked canonical representation;
    `b_values` is parsed and checked as an integrity/consistency alternative.
    """
    path = TARGET / "01_candidate_trajectories.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing saved target trajectories: {path}")
    df = pd.read_csv(path)
    feature_col = next((c for c in ["feature", "candidate"] if c in df.columns), None)
    if feature_col is None:
        raise ValueError("Cannot identify feature column in 01_candidate_trajectories.csv")

    if not {"a_values", "b_values"}.issubset(df.columns):
        # Backward-compatible path for files that already expose four named
        # temporal columns; this does not concatenate alternative variants.
        value_cols = [c for c in df.columns if re.search(r"day|time", str(c), re.I)]
        if len(value_cols) < 4:
            numeric = df.select_dtypes(include=[np.number]).columns.tolist()
            value_cols = numeric[:4]
        if len(value_cols) < 4:
            raise ValueError("Target trajectory file does not expose four temporal values")
        out = df[[feature_col] + value_cols[:4]].copy().set_index(feature_col)
        out.columns = list(YAMANAKA_DAYS)
        return out.apply(pd.to_numeric, errors="coerce"), {
            "source_columns": value_cols[:4],
            "selection_rule": "Use the four explicit temporal columns in Yamanaka order; no concatenation.",
            "n_points": 4,
            "time_points": list(YAMANAKA_DAYS),
        }

    records = []
    max_abs_diff = 0.0
    for _, row in df.iterrows():
        feature = str(row[feature_col])
        a = _parse_value_array(row["a_values"], "a_values", feature)
        b = _parse_value_array(row["b_values"], "b_values", feature)
        diff = float(np.max(np.abs(a - b)))
        max_abs_diff = max(max_abs_diff, diff)
        if not np.allclose(a, b, rtol=1e-10, atol=1e-12):
            raise ValueError(
                f"a_values/b_values disagree for {feature}; refusing to choose silently "
                f"(max_abs_diff={diff})"
            )
        records.append((feature, a))

    out = pd.DataFrame(
        [values for _, values in records],
        index=[feature for feature, _ in records],
        columns=YAMANAKA_DAYS,
    )
    audit = {
        "source_columns": ["a_values", "b_values"],
        "selection_rule": "Treat a_values and b_values as alternative encodings of one trajectory; require both length 4 and numerically agree, then use a_values as the canonical representation.",
        "n_points": 4,
        "time_points": list(YAMANAKA_DAYS),
        "n_features": int(len(records)),
        "max_abs_difference_a_vs_b": max_abs_diff,
    }
    return out, audit


def align_and_stat(target: pd.DataFrame, control: pd.DataFrame):
    common = target.index.intersection(control.index)
    if len(common) < 20:
        raise ValueError(f"Too few common activity features: {len(common)}")
    t = target.loc[common].to_numpy(float)
    c = control.loc[common].to_numpy(float)
    corrs = []
    for a, b in zip(t, c):
        mask = np.isfinite(a) & np.isfinite(b)
        if mask.sum() >= 4 and np.std(a[mask]) > 0 and np.std(b[mask]) > 0:
            corrs.append(np.corrcoef(a[mask], b[mask])[0, 1])
    return float(np.nanmedian(corrs)), np.asarray(corrs), len(common)


def permutation_null(target: pd.DataFrame, control: pd.DataFrame, rng: np.random.Generator):
    common = target.index.intersection(control.index)
    t = target.loc[common].to_numpy(float)
    c = control.loc[common].to_numpy(float)
    vals = []
    for _ in range(N_PERMUTATIONS):
        perm = rng.permutation(len(common))
        cp = c[perm]
        corrs = []
        for a, b in zip(t, cp):
            mask = np.isfinite(a) & np.isfinite(b)
            if mask.sum() >= 4 and np.std(a[mask]) > 0 and np.std(b[mask]) > 0:
                corrs.append(np.corrcoef(a[mask], b[mask])[0, 1])
        vals.append(np.nanmedian(corrs) if corrs else np.nan)
    return np.asarray(vals, float)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    actual = sha256(H3)
    if actual != EXPECTED_SHA:
        raise RuntimeError("GSE263713 SHA mismatch; aborting H3 analysis.")
    gate = OUT / "H3_GSE263713_REPRESENTATION_GATE.json"
    if not gate.exists() or not json.loads(gate.read_text(encoding="utf-8"))["decision"]["eligible_for_h3_similarity"]:
        raise RuntimeError("Representation gate is absent or not eligible; aborting H3 analysis.")

    expr = map_hgnc(parse_h3())
    nets = activity_networks()
    target, target_audit = target_traj()
    control_features = {}
    for kind, net in nets.items():
        activity = score_activity(expr, net, kind)
        control_features[kind] = trajectory_matrix(activity)

    rng = np.random.default_rng(SEED)
    rows = []
    nulls = {}
    for kind, control in control_features.items():
        obs, corrs, n_common = align_and_stat(target, control)
        null = permutation_null(target, control, rng)
        nulls[kind] = null
        p = (1 + np.sum(null >= obs)) / (N_PERMUTATIONS + 1)
        rows.append({"control":"GSE263713","feature_family":kind,"n_common_features":n_common,"observed_median_trajectory_corr":obs,"null_median":float(np.nanmedian(null)),"null_q95":float(np.nanquantile(null,0.95)),"empirical_one_sided_p":float(p)})
        pd.DataFrame({"null_statistic":null}).to_csv(OUT / f"H3_GSE263713_{kind}_NULL.csv", index=False)
    result = {
        "status":"ok",
        "accession":"GSE263713",
        "seed":SEED,
        "n_permutations":N_PERMUTATIONS,
        "observed_target_trajectory_file":str((TARGET/"01_candidate_trajectories.csv").relative_to(ROOT)),
        "target_trajectory_mapping":target_audit,
        "controls":rows,
        "decision":{"status":"H3_UNRESOLVED","reason":"This stage establishes the prospective control comparison only; final H3 status requires the precommitted two-control set and relative-control margin."},
        "guardrail":"This analysis does not establish process specificity, causality, or mechanism. GSE263713 is one orthogonal temporal control candidate; it cannot by itself produce the final two-control H3 decision."
    }
    (OUT/"H3_GSE263713_SIMILARITY_RESULT.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
    pd.DataFrame(rows).to_csv(OUT/"H3_GSE263713_SIMILARITY_RESULT.csv",index=False)
    print(json.dumps(result,indent=2))

if __name__ == "__main__":
    main()
