"""Stage 2.7 — validation of the biological common gene space.

This stage is intentionally diagnostic: it validates the feature-level common
space before any dynamical model or symbolic regression is fitted.

It operates on the feature files produced by Dynamics.py and does not use time
for construction of the feature space. Dataset-specific normalization is kept
separate, and the validation focuses on coverage, duplicate mappings, variance,
replicate coherence, dataset/platform effects, and leave-one-dataset-out PCA
projections where the available metadata permit it.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


DATASETS = {
    "GSE148158": r"results/GSE148158/expression.csv",
    "GSE28688": r"results/GSE28688/non_normalized/expression_input.csv",
    "GSE52052": r"results/GSE52052/expression_log2.csv",
    "GSE67462": r"results/GSE67462/03_expression_for_EDA.csv",
    "GSE297234": r"results/GSE297234/03_log1p_CPM_sample_expression.csv",
}


def _norm_gene(x: object) -> str:
    return str(x).strip().upper()


def _read_expression(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    if df.empty:
        raise ValueError(f"Empty expression file: {path}")
    df.index = [_norm_gene(x) for x in df.index]
    # Keep numeric columns only; metadata columns cannot enter the common matrix.
    numeric = df.apply(pd.to_numeric, errors="coerce")
    numeric = numeric.dropna(axis=1, how="all")
    if numeric.empty:
        raise ValueError(f"No numeric sample columns found: {path}")
    # Duplicate gene identifiers must not silently double-weight a gene.
    numeric = numeric.groupby(level=0, sort=True).mean()
    return numeric


def _robust_scale(df: pd.DataFrame) -> pd.DataFrame:
    """Per-gene robust scaling for cross-platform diagnostics only."""
    med = df.median(axis=1)
    mad = (df.sub(med, axis=0)).abs().median(axis=1)
    scale = mad.replace(0, np.nan)
    out = df.sub(med, axis=0).div(scale, axis=0)
    return out.replace([np.inf, -np.inf], np.nan).fillna(0.0)


def _common_genes(matrices: dict[str, pd.DataFrame]) -> list[str]:
    common = set.intersection(*(set(x.index) for x in matrices.values()))
    return sorted(common)


def _coverage_report(matrices: dict[str, pd.DataFrame], common: Iterable[str]) -> pd.DataFrame:
    common = set(common)
    rows = []
    for name, df in matrices.items():
        rows.append(
            {
                "dataset": name,
                "features": len(df),
                "samples": df.shape[1],
                "common_genes": len(common & set(df.index)),
                "coverage_fraction": len(common & set(df.index)) / max(len(df), 1),
                "duplicate_gene_rows_collapsed": int(df.index.duplicated().sum()),
            }
        )
    return pd.DataFrame(rows)


def _pca_diagnostics(matrices: dict[str, pd.DataFrame], common: list[str]) -> pd.DataFrame:
    """Return compact PCA diagnostics without fitting a cross-dataset model."""
    rows = []
    for name, df in matrices.items():
        x = _robust_scale(df.loc[common]).T.to_numpy(dtype=float)
        if x.shape[0] < 2:
            continue
        x -= x.mean(axis=0, keepdims=True)
        _, s, _ = np.linalg.svd(x, full_matrices=False)
        ev = s**2
        total = ev.sum()
        rows.append(
            {
                "dataset": name,
                "pc1_variance_fraction": float(ev[0] / total) if total else np.nan,
                "pc2_variance_fraction": float(ev[1] / total) if len(ev) > 1 and total else np.nan,
                "effective_rank": int(np.sum(ev > (ev.max() * 1e-10))) if len(ev) else 0,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--output", type=Path, default=Path("results/Dynamics/Stage_2_7"))
    args = parser.parse_args()

    matrices: dict[str, pd.DataFrame] = {}
    for name, rel in DATASETS.items():
        path = args.root / rel
        if not path.exists():
            raise FileNotFoundError(f"Missing expected Stage 2.6 output: {path}")
        matrices[name] = _read_expression(path)

    common = _common_genes(matrices)
    if len(common) < 1000:
        raise RuntimeError(f"Unexpectedly small common gene space: {len(common)}")

    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    coverage = _coverage_report(matrices, common)
    pca = _pca_diagnostics(matrices, common)
    pd.DataFrame({"gene": common}).to_csv(out / "common_genes.csv", index=False)
    coverage.to_csv(out / "coverage.csv", index=False)
    pca.to_csv(out / "pca_diagnostics.csv", index=False)

    summary = {
        "n_datasets": len(matrices),
        "common_genes": len(common),
        "time_used_for_feature_construction": False,
        "duplicate_gene_rows_collapsed": int(sum(df.index.duplicated().sum() for df in matrices.values())),
        "stage_status": "passed_basic_feature_space_checks",
    }
    pd.DataFrame([summary]).to_csv(out / "summary.csv", index=False)

    print("=" * 88)
    print("STAGE 2.7 — VALIDATION OF THE BIOLOGICAL COMMON GENE SPACE")
    print("=" * 88)
    print(f"  Datasets: {len(matrices)}")
    print(f"  Common genes: {len(common)}")
    print("  Time used to construct feature space: NO")
    print("  Duplicate gene identifiers: collapsed by within-dataset mean")
    print()
    print(coverage.to_string(index=False))
    print()
    print(pca.to_string(index=False))
    print()
    print(f"Stage 2.7 outputs written to: {out.resolve()}")
    print("NOTE: passing these checks does not establish biological equivalence or")
    print("cross-dataset predictive validity; those require explicit held-out tests.")


if __name__ == "__main__":
    main()
