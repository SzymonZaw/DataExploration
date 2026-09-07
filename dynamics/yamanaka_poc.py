"""Proof-of-feasibility for Yamanaka-factor perturbation discovery.

This module deliberately stays small.  It asks whether recent human
reprogramming datasets contain a reproducible, factor-specific expression
signal before introducing a larger mechanistic model.

Inputs
------
* Data/GSE297233_raw_counts_matrix.csv.gz (tracked, bulk RNA-seq)
* Data/GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz (tracked)
* Data/GSE297234_*_SEVOSKM.rds (local only, optional; not required for the
  first POC because Seurat/SingleCellExperiment RDS files are not reliably
  readable by pure Python).

The POC performs:
1. automatic discovery of gene/sample columns in the two CSV matrices;
2. conservative log-CPM normalization;
3. sample-condition parsing for GFP/control, OSK and individual factors;
4. factor-vs-control signatures and cosine similarity;
5. OSK cross-dataset direction concordance on overlapping genes;
6. leave-one-sample-out classification of perturbation identity using PCA +
   logistic regression, when enough labeled samples exist;
7. a compact report suitable for the PhD feasibility discussion.

The output is evidence of *perturbation signal*, not proof of mechanism or
causality.  It is intentionally independent of the Stage 2.6 common-space
benchmark so that a failure there does not contaminate this feasibility test.
"""
from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.model_selection import LeaveOneOut, cross_val_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11_yamanaka_poc"

FILES = {
    "GSE297233": DATA / "GSE297233_raw_counts_matrix.csv.gz",
    "GSE304042": DATA / "GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz",
}

FACTOR_TOKENS = {
    "OSK": ("OSK", "OSK_"),
    "OCT4": ("OCT4", "OCT4_", "OCT4YR"),
    "SOX2": ("SOX2", "SOX2_"),
    "KLF4": ("KLF4", "KLF4_"),
}
CONTROL_TOKENS = ("GFP", "CONTROL", "CTRL", "UNTRANSF", "UNTREATED", "EMPTY")


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing dataset: {path}")
    return pd.read_csv(path, compression="infer", low_memory=False)


def _gene_column(df: pd.DataFrame) -> str:
    preferred = ("gene", "genes", "gene_symbol", "symbol", "gene_id", "ensembl")
    lower = {str(c).lower(): c for c in df.columns}
    for name in preferred:
        if name in lower:
            return lower[name]
    # The first mostly-string column is a safer fallback than assuming column 0.
    scores = []
    for c in df.columns:
        s = df[c]
        numeric = pd.to_numeric(s, errors="coerce").notna().mean()
        scores.append((numeric, c))
    non_numeric = [c for numeric, c in scores if numeric < 0.5]
    if non_numeric:
        return non_numeric[0]
    return str(df.columns[0])


def _numeric_sample_columns(df: pd.DataFrame, gene_col: str) -> list[str]:
    cols = []
    for c in df.columns:
        if c == gene_col:
            continue
        numeric = pd.to_numeric(df[c], errors="coerce")
        if numeric.notna().mean() >= 0.8:
            cols.append(str(c))
    if len(cols) < 2:
        raise ValueError(f"Could not identify >=2 numeric sample columns in {gene_col} matrix.")
    return cols


def _clean_gene_names(values) -> pd.Index:
    genes = pd.Index(values.astype(str).str.strip())
    genes = genes.str.replace(r"\\.\\d+$", "", regex=True)
    genes = genes.str.upper()
    return genes


def _collapse_genes(X: pd.DataFrame) -> pd.DataFrame:
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    X = X.groupby(level=0, sort=False).sum()
    return X


def _log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    counts = counts.clip(lower=0.0)
    library = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(library, axis=1) * 1e6
    return np.log1p(cpm)


def load_expression(path: Path) -> tuple[pd.DataFrame, dict]:
    df = _read_csv(path)
    gene_col = _gene_column(df)
    sample_cols = _numeric_sample_columns(df, gene_col)
    genes = _clean_gene_names(df[gene_col])
    counts = df.loc[:, sample_cols].copy()
    counts.columns = [str(c) for c in sample_cols]
    counts.index = genes
    counts = _collapse_genes(counts)
    expr = _log_cpm(counts)
    audit = {
        "path": str(path.relative_to(ROOT)),
        "gene_column": gene_col,
        "n_genes": int(expr.shape[0]),
        "n_samples": int(expr.shape[1]),
        "sample_names": list(expr.columns),
    }
    return expr, audit


def classify_condition(sample: str) -> str:
    s = str(sample).upper()
    for factor, tokens in FACTOR_TOKENS.items():
        if any(token in s for token in tokens):
            if "OSK" in s and factor in {"OCT4", "SOX2", "KLF4"}:
                continue
            return factor
    if "OSK" in s:
        return "OSK"
    if any(token in s for token in CONTROL_TOKENS):
        return "CONTROL"
    # The 2025 dataset contains an OCT4YR+SK perturbation; retain it as a
    # distinct condition rather than silently calling it OSK.
    if "YR" in s and "SK" in s:
        return "OCT4YR_SK"
    return "OTHER"


def condition_table(expr: pd.DataFrame, dataset: str) -> pd.DataFrame:
    rows = []
    for sample in expr.columns:
        rows.append({"dataset": dataset, "sample": sample, "condition": classify_condition(sample)})
    return pd.DataFrame(rows)


def _signature(expr: pd.DataFrame, labels: pd.Series, condition: str) -> pd.Series | None:
    ctrl = expr.loc[:, labels.eq("CONTROL")]
    case = expr.loc[:, labels.eq(condition)]
    if ctrl.shape[1] == 0 or case.shape[1] == 0:
        return None
    return case.mean(axis=1) - ctrl.mean(axis=1)


def _cosine(a: pd.Series, b: pd.Series) -> float:
    common = a.index.intersection(b.index)
    if len(common) < 20:
        return float("nan")
    x = a.loc[common].to_numpy(float)
    y = b.loc[common].to_numpy(float)
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    return float(np.dot(x, y) / denom) if denom > 0 else float("nan")


def _top_genes(signature: pd.Series, n: int = 20) -> tuple[list[str], list[str]]:
    up = signature.sort_values(ascending=False).head(n).index.tolist()
    down = signature.sort_values(ascending=True).head(n).index.tolist()
    return up, down


def cross_dataset_osk(exprs: dict[str, pd.DataFrame], labels: dict[str, pd.DataFrame]) -> dict:
    sigs = {}
    for ds, expr in exprs.items():
        lab = labels[ds].set_index("sample")["condition"]
        sig = _signature(expr, lab, "OSK")
        if sig is not None:
            sigs[ds] = sig
    pairwise = []
    keys = sorted(sigs)
    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            common = sigs[a].index.intersection(sigs[b].index)
            sa, sb = sigs[a].loc[common], sigs[b].loc[common]
            pearson = float(sa.corr(sb)) if len(common) >= 20 else float("nan")
            sign_concordance = float((np.sign(sa) == np.sign(sb)).mean()) if len(common) else float("nan")
            pairwise.append({
                "dataset_a": a,
                "dataset_b": b,
                "n_common_genes": int(len(common)),
                "pearson_signature": pearson,
                "sign_concordance": sign_concordance,
                "cosine_signature": _cosine(sa, sb),
            })
    return {"n_osk_signatures": len(sigs), "pairwise": pairwise}


def perturbation_classification(exprs: dict[str, pd.DataFrame], labels: dict[str, pd.DataFrame]) -> dict:
    rows = []
    for ds, expr in exprs.items():
        lab = labels[ds].set_index("sample")["condition"]
        keep = lab.isin(["CONTROL", "OSK", "OCT4", "SOX2", "KLF4", "OCT4YR_SK"])
        expr = expr.loc[:, keep.index[keep]]
        y = lab.loc[expr.columns]
        # Need at least 3 classes and one observation per class after filtering.
        if y.nunique() < 3 or y.value_counts().min() < 1:
            continue
        X = expr.T.to_numpy(float)
        n_components = max(1, min(10, X.shape[0] - 1, X.shape[1] - 1))
        model = make_pipeline(
            StandardScaler(),
            PCA(n_components=n_components, random_state=0),
            LogisticRegression(max_iter=3000, multi_class="auto", random_state=0),
        )
        # Leave-one-out is intentionally used because these datasets are tiny.
        cv = LeaveOneOut()
        try:
            scores = cross_val_score(model, X, y.to_numpy(), cv=cv)
            accuracy = float(scores.mean())
        except Exception:
            accuracy = float("nan")
        rows.append({
            "dataset": ds,
            "n_samples": int(len(y)),
            "n_classes": int(y.nunique()),
            "accuracy_loo": accuracy,
            "class_counts": json.dumps(y.value_counts().to_dict(), sort_keys=True),
        })
    return {"datasets": rows}


def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    exprs, labels, audits = {}, {}, []
    for ds, path in FILES.items():
        expr, audit = load_expression(path)
        lab = condition_table(expr, ds)
        exprs[ds], labels[ds] = expr, lab
        audits.append({**audit, "condition_counts": lab.condition.value_counts().to_dict()})

    signature_rows = []
    for ds, expr in exprs.items():
        lab = labels[ds].set_index("sample")["condition"]
        for condition in ("OSK", "OCT4", "SOX2", "KLF4", "OCT4YR_SK"):
            sig = _signature(expr, lab, condition)
            if sig is None:
                continue
            up, down = _top_genes(sig)
            signature_rows.append({
                "dataset": ds,
                "condition": condition,
                "n_genes": int(len(sig)),
                "mean_abs_effect": float(np.abs(sig).mean()),
                "top_up_genes": ";".join(up),
                "top_down_genes": ";".join(down),
            })

    signatures = pd.DataFrame(signature_rows)
    signatures.to_csv(OUT / "02_perturbation_signatures.csv", index=False)
    pd.DataFrame(audits).to_json(OUT / "01_input_audit.json", orient="records", indent=2)

    cross = cross_dataset_osk(exprs, labels)
    classification = perturbation_classification(exprs, labels)
    summary = {
        "status": "ok",
        "datasets_loaded": sorted(exprs),
        "input_audit": audits,
        "cross_dataset_osk": cross,
        "perturbation_classification": classification,
        "interpretation_guardrail": (
            "Positive cross-dataset signature agreement is evidence of reproducible "
            "perturbation signal, not mechanistic causality. A useful next step is "
            "pathway/TF activity and time-resolved scRNA pseudobulk validation."
        ),
    }
    (OUT / "03_poc_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, default=str))
