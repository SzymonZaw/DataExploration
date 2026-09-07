from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
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


def _clean_sample_name(value: object) -> str:
    return str(value).strip().replace(" ", "_")


def load_expression(path: Path) -> tuple[pd.DataFrame, dict]:
    if not path.exists():
        raise FileNotFoundError(f"Missing dataset: {path}")
    df = pd.read_csv(path, compression="infer")
    if df.shape[1] < 2:
        raise ValueError(f"Expected genes x samples matrix in {path}")
    gene_col = df.columns[0]
    genes = df[gene_col].astype(str).str.strip()
    X = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    X.index = genes
    X = X.groupby(level=0).mean()
    X = X.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    audit = {"path": str(path.relative_to(ROOT)), "n_genes": int(X.shape[0]), "n_samples": int(X.shape[1])}
    return X, audit


def condition_table(expr: pd.DataFrame, dataset: str) -> pd.DataFrame:
    rows = []
    for sample in expr.columns:
        s = _clean_sample_name(sample).upper()
        condition = "OTHER"
        if dataset == "GSE297233":
            if "CONTROL" in s or s.startswith("CTRL") or "D0" in s:
                condition = "CONTROL"
            elif "OCT4YR" in s and "SK" in s:
                condition = "OCT4YR_SK"
            elif "OSK" in s:
                condition = "OSK"
        elif dataset == "GSE304042":
            # GSE304042 uses GFP-a/b/c as the untreated/control condition.
            if "CONTROL" in s or "CTRL" in s or s.startswith("GFP"):
                condition = "CONTROL"
            elif "OSK" in s:
                condition = "OSK"
            elif "OCT4" in s or "OCT4YR" in s:
                condition = "OCT4"
            elif "SOX2" in s:
                condition = "SOX2"
            elif "KLF4" in s:
                condition = "KLF4"
        rows.append({"sample": sample, "condition": condition})
    return pd.DataFrame(rows)


def _signature(expr: pd.DataFrame, labels: pd.Series, condition: str) -> pd.Series | None:
    keep = labels.isin(["CONTROL", condition])
    if keep.sum() < 2 or labels[keep].nunique() < 2:
        return None
    a = expr.loc[:, labels.index[keep & (labels == condition)]]
    b = expr.loc[:, labels.index[keep & (labels == "CONTROL")]]
    if a.shape[1] == 0 or b.shape[1] == 0:
        return None
    return np.log1p(a).mean(axis=1) - np.log1p(b).mean(axis=1)


def _top_genes(sig: pd.Series, n: int = 20) -> tuple[list[str], list[str]]:
    sig = sig.sort_values()
    return sig.tail(n).sort_values(ascending=False).index.tolist(), sig.head(n).index.tolist()


def _cosine(a: pd.Series, b: pd.Series) -> float:
    av, bv = a.to_numpy(float), b.to_numpy(float)
    denom = np.linalg.norm(av) * np.linalg.norm(bv)
    return float(np.dot(av, bv) / denom) if denom > 0 else float("nan")


def cross_dataset_osk(exprs: dict[str, pd.DataFrame], labels: dict[str, pd.DataFrame]) -> dict:
    signatures = {}
    for ds, expr in exprs.items():
        lab = labels[ds].set_index("sample")["condition"]
        sig = _signature(expr, lab, "OSK")
        if sig is not None:
            signatures[ds] = sig
    pairwise = []
    names = sorted(signatures)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            common = signatures[a].index.intersection(signatures[b].index)
            if len(common) < 100:
                continue
            sa, sb = signatures[a].loc[common], signatures[b].loc[common]
            pairwise.append({
                "dataset_a": a,
                "dataset_b": b,
                "n_common_genes": int(len(common)),
                "pearson": float(sa.corr(sb, method="pearson")),
                "spearman": float(sa.corr(sb, method="spearman")),
                "cosine_signature": _cosine(sa, sb),
            })
    return {"n_osk_signatures": len(signatures), "pairwise": pairwise}


def perturbation_classification(exprs: dict[str, pd.DataFrame], labels: dict[str, pd.DataFrame]) -> dict:
    rows = []
    for ds, expr in exprs.items():
        lab = labels[ds].set_index("sample")["condition"]
        keep = lab.isin(["CONTROL", "OSK", "OCT4", "SOX2", "KLF4", "OCT4YR_SK"])
        sample_names = lab.index[keep]
        expr = expr.loc[:, expr.columns.intersection(sample_names)]
        y = lab.loc[expr.columns]
        if y.nunique() < 3 or len(y) < 4:
            continue
        X = expr.T.to_numpy(float)
        n_components = max(1, min(10, X.shape[0] - 1, X.shape[1] - 1))
        model = make_pipeline(
            StandardScaler(),
            PCA(n_components=n_components, random_state=0),
            # New scikit-learn releases removed the deprecated multi_class argument.
            LogisticRegression(max_iter=3000, random_state=0),
        )
        try:
            scores = cross_val_score(model, X, y.to_numpy(), cv=LeaveOneOut())
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

    pd.DataFrame(signature_rows).to_csv(OUT / "02_perturbation_signatures.csv", index=False)
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
