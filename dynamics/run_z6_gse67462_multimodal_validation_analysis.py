"""Analyze GSE67520 peak calls against the GSE67462 common dynamic trajectory.

Diagnostic only. It does not modify frozen Z6 predictive support.

Primary representation: promoter-associated peak burden, defined as the sum of
peak scores whose peak midpoint falls within +/-2 kb of an annotated TSS.
This deliberately avoids arbitrary nearest-gene assignment for distal peaks.
The analysis is temporal: each modality is compared with the GSE67462 common
branch trajectory over day 0,1,3,5,7,11,15,18. iPSC is retained in parsing
but is not used in the primary trajectory test because the expression audit
has no validated iPSC time point.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from dynamics.run_z6_within_system_audit import _build_system_data
from dynamics import validation

ROOT = Path(__file__).resolve().parents[1]
MODALITIES = {
    "total_oct4": ("Oct4", "broadPeak"),
    "h3k4me1": ("H3K4me1", "bed"),
    "h3k27ac": ("H3K27ac", "bed"),
    "h3k4me3": ("H3K4me3", "bed"),
    "h3k27me3": ("H3K27me3", "bed"),
    "rnapii": ("RNAPII", "bed"),
}
TIME_RE = re.compile(r"_(d(?:0|1|3|5|7|11|15|18)|ipsc)_peaks", re.I)
DEFAULT_GTF = "Data/GSE67520/mm9.refGene.gtf.gz"
GTF_URL = "https://hgdownload.cse.ucsc.edu/goldenpath/mm9/bigZips/genes/mm9.refGene.gtf.gz"
PROMOTER_HALF_WIDTH = 2000
N_PERMUTATIONS = 1000
SEED = 412


def _spearman(x, y):
    a = pd.Series(np.asarray(x, dtype=float)).rank(method="average").to_numpy()
    b = pd.Series(np.asarray(y, dtype=float)).rank(method="average").to_numpy()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _parse_attrs(text):
    out = {}
    for key, value in re.findall(r'(\S+)\s+"([^"]+)"', text):
        out[key] = value
    return out


def _ensure_gtf(path: Path, allow_download: bool):
    if path.exists():
        return
    if not allow_download:
        raise FileNotFoundError(f"Missing annotation: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading mm9 RefSeq annotation to {path} ...")
    urllib.request.urlretrieve(GTF_URL, path)


def _load_tss(gtf: Path):
    rows = []
    opener = gzip.open if str(gtf).endswith(".gz") else open
    with opener(gtf, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs = _parse_attrs(fields[8])
            name = attrs.get("gene_name") or attrs.get("gene_id")
            if not name:
                continue
            start = int(fields[3]) - 1
            end = int(fields[4])
            tss = start if fields[6] != "-" else end - 1
            rows.append((fields[0], tss, name))
    tss = pd.DataFrame(rows, columns=["chrom", "tss", "gene"])
    if tss.empty:
        raise RuntimeError("No transcript/TSS records found in GTF")
    # One TSS per gene/chromosome minimizes transcript-driven weighting.
    tss = tss.groupby(["chrom", "gene"], as_index=False)["tss"].median()
    tss["start"] = (tss["tss"] - PROMOTER_HALF_WIDTH).astype(int)
    tss["end"] = (tss["tss"] + PROMOTER_HALF_WIDTH + 1).astype(int)
    return tss


def _parse_peak_file(path: Path):
    is_broad = path.name.lower().endswith("broadpeak.gz")
    rows = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                start, end = int(f[1]), int(f[2])
                score_idx = 4
                score = float(f[score_idx]) if len(f) > score_idx else 1.0
            except ValueError:
                continue
            rows.append((f[0], start, end, score))
    return rows


def _peak_to_promoter_burden(path: Path, tss: pd.DataFrame):
    # Midpoint assignment is deterministic and prevents a long peak from
    # inflating multiple neighboring promoters.
    peaks = _parse_peak_file(path)
    by_chrom = {c: g.sort_values("tss") for c, g in tss.groupby("chrom")}
    values = {}
    assigned = 0
    for chrom, start, end, score in peaks:
        if chrom not in by_chrom:
            continue
        midpoint = (start + end) // 2
        g = by_chrom[chrom]
        pos = np.searchsorted(g["tss"].to_numpy(), midpoint)
        candidates = []
        for j in (pos - 1, pos):
            if 0 <= j < len(g):
                row = g.iloc[j]
                if row.start <= midpoint <= row.end:
                    candidates.append(row)
        if not candidates:
            continue
        row = min(candidates, key=lambda r: abs(float(r.tss) - midpoint))
        gene = str(row.gene)
        values[gene] = values.get(gene, 0.0) + max(score, 0.0)
        assigned += 1
    return values, len(peaks), assigned


def _discover_files(root: Path):
    files = list(root.glob("*.gz"))
    result = {k: [] for k in MODALITIES}
    for p in files:
        low = p.name.lower()
        if "flag_" in low:
            continue
        for key in MODALITIES:
            token = {"total_oct4": "oct4_", "h3k4me1": "k4me1_", "h3k27ac": "k27ac_", "h3k4me3": "k4me3_", "h3k27me3": "k27me3_", "rnapii": "rnapii_"}[key]
            if token in low:
                m = TIME_RE.search(p.name)
                if m:
                    tag = m.group(1).lower()
                    time = "iPSC" if tag == "ipsc" else int(tag[1:])
                    result[key].append((time, p))
                break
    for key in result:
        result[key].sort(key=lambda x: (x[0] == "iPSC", x[0] if isinstance(x[0], int) else 999))
    return result


def _expression_common(matrix, metadata):
    data, _ = _build_system_data(matrix, metadata, "GSE67462")
    keys = sorted(data)
    if len(keys) != 2:
        raise RuntimeError(f"Expected two GSE67462 branches, got {keys}")
    ta, xa = data[keys[0]]
    tb, xb = data[keys[1]]
    common_times = sorted(set(map(float, ta)) & set(map(float, tb)))
    ia = {float(t): i for i, t in enumerate(ta)}
    ib = {float(t): i for i, t in enumerate(tb)}
    xa = np.vstack([xa[ia[t]] for t in common_times])
    xb = np.vstack([xb[ib[t]] for t in common_times])
    common = ((xa - xa[0]) + (xb - xb[0])) / 2.0
    genes = pd.Index(matrix.index.astype(str))
    return np.asarray(common_times), pd.DataFrame(common, index=common_times, columns=genes), keys


def _build_modality_matrix(files, tss):
    time_rows = []
    provenance = []
    for key, entries in files.items():
        for time, path in entries:
            vals, n_peaks, n_assigned = _peak_to_promoter_burden(path, tss)
            time_rows.append((key, time, vals))
            provenance.append({"modality": key, "time": str(time), "file": str(path), "n_peaks": n_peaks, "n_promoter_assigned": n_assigned})
    return time_rows, pd.DataFrame(provenance)


def _trajectory_concordance(expr, modality_values):
    genes = expr.columns
    rows = []
    for modality in sorted(set(m for m, _, _ in modality_values)):
        entries = [(t, v) for m, t, v in modality_values if m == modality and isinstance(t, int)]
        entries = sorted(entries, key=lambda x: x[0])
        times = [t for t, _ in entries if t in expr.index]
        if len(times) < 5:
            rows.append({"modality": modality, "n_timepoints": len(times), "n_genes": 0, "median_gene_spearman": np.nan, "mean_gene_spearman": np.nan, "observed_global_spearman": np.nan, "status": "INSUFFICIENT_TIMEPOINTS"})
            continue
        modality_df = pd.DataFrame([{g: v.get(g, 0.0) for g in genes} for t, v in entries if t in expr.index], index=[t for t, _ in entries if t in expr.index])
        common_genes = genes[genes.isin(modality_df.columns)]
        gene_rhos = [_spearman(expr.loc[times, g].to_numpy(), modality_df.loc[times, g].to_numpy()) for g in common_genes]
        gene_rhos = np.asarray([r for r in gene_rhos if np.isfinite(r)])
        expr_flat = expr.loc[times, common_genes].to_numpy().ravel()
        mod_flat = modality_df.loc[times, common_genes].to_numpy().ravel()
        global_rho = _spearman(expr_flat, mod_flat)
        rows.append({"modality": modality, "n_timepoints": len(times), "n_genes": len(gene_rhos), "median_gene_spearman": float(np.nanmedian(gene_rhos)) if len(gene_rhos) else np.nan, "mean_gene_spearman": float(np.nanmean(gene_rhos)) if len(gene_rhos) else np.nan, "observed_global_spearman": global_rho, "status": "OK"})
    return pd.DataFrame(rows)


def _permutation_test(expr, modality_values, n_permutations, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for modality in sorted(set(m for m, _, _ in modality_values)):
        entries = [(t, v) for m, t, v in modality_values if m == modality and isinstance(t, int) and t in expr.index]
        entries.sort(key=lambda x: x[0])
        if len(entries) < 5:
            continue
        times = [t for t, _ in entries]
        genes = expr.columns
        mod_df = pd.DataFrame([{g: v.get(g, 0.0) for g in genes} for _, v in entries], index=times)
        observed = _spearman(expr.loc[times, genes].to_numpy().ravel(), mod_df.to_numpy().ravel())
        null = []
        for _ in range(n_permutations):
            perm = rng.permutation(times)
            null.append(_spearman(expr.loc[times, genes].to_numpy().ravel(), mod_df.loc[perm, genes].to_numpy().ravel()))
        null = np.asarray([x for x in null if np.isfinite(x)])
        p = float((1 + np.sum(null >= observed)) / (1 + len(null))) if len(null) else np.nan
        rows.append({"modality": modality, "observed_global_spearman": observed, "null_mean": float(np.mean(null)) if len(null) else np.nan, "null_q95": float(np.quantile(null, 0.95)) if len(null) else np.nan, "permutation_p": p, "n_permutations": len(null)})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--modality-root", default="Data/GSE67520")
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--output", default="results/Dynamics/z6_gse67462_multimodal_validation")
    parser.add_argument("--no-download-annotation", action="store_true")
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    args = parser.parse_args()

    root = Path(args.modality_root)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    gtf = Path(args.gtf)
    _ensure_gtf(gtf, not args.no_download_annotation)
    tss = _load_tss(gtf)

    matrix, metadata = validation._load_common_space()
    expr_times, expr, branches = _expression_common(matrix, metadata)
    files = _discover_files(root)
    if any(not files[k] for k in MODALITIES):
        missing = [k for k in MODALITIES if not files[k]]
        raise RuntimeError(f"Missing modality files: {missing}")

    values, provenance = _build_modality_matrix(files, tss)
    concordance = _trajectory_concordance(expr, values)
    perm = _permutation_test(expr, values, args.permutations, SEED)
    result = concordance.merge(perm, on="modality", how="left")

    supported = result[(result.status == "OK") & (result.permutation_p < 0.05) & (result.observed_global_spearman > result.null_q95)]
    if len(supported) >= 2:
        interpretation = "MULTIMODAL_DYNAMIC_SUPPORT"
    elif len(supported) == 1 or np.any((result.permutation_p < 0.05).fillna(False)):
        interpretation = "PARTIAL_MULTIMODAL_SUPPORT"
    else:
        interpretation = "EXPRESSION_ONLY"

    summary = {
        "dataset_expression": "GSE67462",
        "dataset_multimodal": "GSE67520",
        "branches": branches,
        "expression_timepoints": expr_times.tolist(),
        "n_expression_genes": int(expr.shape[1]),
        "n_annotated_tss": int(len(tss)),
        "promoter_definition": "peak midpoint within +/-2000 bp of one representative TSS per gene",
        "n_modalities": len(MODALITIES),
        "n_supported_modalities": int(len(supported)),
        "interpretation": interpretation,
        "biological_specificity_established": False,
        "z6_predictive_support_changed": False,
        "note": "Within-program multimodal coherence is not independent-system validation.",
    }

    provenance.to_csv(out / "01_peak_provenance.csv", index=False)
    result.to_csv(out / "02_modality_temporal_concordance.csv", index=False)
    (out / "03_annotation_summary.json").write_text(json.dumps({"gtf": str(gtf), "gtf_source": GTF_URL, "n_tss": int(len(tss)), "promoter_half_width": PROMOTER_HALF_WIDTH}, indent=2), encoding="utf-8")
    (out / "04_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.DataFrame({"expression_time": expr_times}).to_csv(out / "05_expression_time_grid.csv", index=False)

    print("\nGSE67462/GSE67520 multimodal validation complete.")
    print(result.to_string(index=False))
    print(f"\nInterpretation: {interpretation}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
