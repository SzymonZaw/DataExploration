"""Audit GSE67462 feature-identifier provenance and mapping to GSE67520 genes.

Diagnostic only. This script does not change frozen Z6 support criteria or results.
It separates three questions that were previously conflated:

1. What namespace is used by the Stage 2.6/common-space expression matrix?
2. What namespace is present in the original GSE67462/GPL19972 processed data?
3. Can the expression features be mapped unambiguously to the gene symbols used
   by the GSE67520 regulatory/TSS annotation?

The audit prefers local provenance files. It supports both a simple tabular
GPL19972 annotation and the authoritative NCBI GEO SOFT family file.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from io import StringIO
from pathlib import Path

import pandas as pd

from dynamics import validation
from dynamics.run_z6_gse67462_multimodal_validation_analysis import _load_tss


def _norm_symbol(value: object) -> str:
    s = str(value).strip().strip('"').upper()
    return "" if s in {"", "NAN", "NA", "NONE", "NULL", "-", "."} else s


def _norm_refseq(value: object) -> str:
    s = str(value).strip().strip('"')
    s = re.sub(r"_at$", "", s, flags=re.I)
    s = re.sub(r"\.\d+$", "", s)
    return s.upper()


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Platform annotation file not found: {path}")
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return pd.read_csv(fh, sep=None, engine="python", dtype=str)
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def _read_soft_platform(path: Path) -> pd.DataFrame:
    """Extract the !platform_table_begin/!platform_table_end section from GEO SOFT."""
    if not path.exists():
        raise FileNotFoundError(f"GPL19972 SOFT file not found: {path}")

    opener = gzip.open if path.suffix.lower() == ".gz" else open
    lines: list[str] = []
    in_table = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if in_table and line.startswith("!platform_table_end"):
                break
            if in_table:
                lines.append(line)

    if not lines:
        raise ValueError(
            f"No !platform_table_begin/!platform_table_end section found in SOFT file: {path}"
        )

    return pd.read_csv(StringIO("\n".join(lines)), sep="\t", dtype=str, keep_default_na=False)


def _find_feature_column(df: pd.DataFrame) -> str | None:
    preferred = ["ID_REF", "ID", "probe", "probe_id", "feature", "feature_id"]
    for name in preferred:
        if name in df.columns:
            return name
    return str(df.columns[0]) if len(df.columns) else None


def _find_symbol_column(df: pd.DataFrame) -> str | None:
    preferred = ["gene_symbol", "symbol", "Symbol", "Gene Symbol", "gene", "gene_name", "GeneName"]
    for name in preferred:
        if name in df.columns:
            return name
    return None


def _find_refseq_column(df: pd.DataFrame) -> str | None:
    preferred = ["GB_ACC", "RefSeq", "refseq", "refseq_id", "accession"]
    for name in preferred:
        if name in df.columns:
            return name
    return None


def _extract_description_symbol(description: object, tss_symbols: set[str]) -> str:
    """Recover the gene symbol from GPL19972's RefSeq-style Description field.

    GPL19972 exposes ID, GB_ACC and Description but no explicit Symbol column.
    The Description commonly contains the symbol in parentheses, e.g.
    '... (Zfp85-rs1), mRNA'. We prefer parenthesized tokens that are present in
    the independently supplied TSS gene-symbol namespace.
    """
    text = str(description).strip()
    if not text:
        return ""
    candidates = re.findall(r"\(([^()]*)\)", text)
    for candidate in candidates:
        key = _norm_symbol(candidate)
        if key in tss_symbols:
            return key
    return ""


def _candidate_platform_tables(root: Path) -> list[Path]:
    roots = [root / "Data" / "GSE67462", root / "Data", root / "results" / "GSE67462"]
    names = {
        "GPL19972.csv", "GPL19972.tsv", "GPL19972.txt", "GPL19972.txt.gz",
        "GPL19972_family.soft", "GPL19972_family.soft.gz",
        "GSE67462_GPL19972.csv", "GSE67462_GPL19972.tsv", "GSE67462_GPL19972.txt",
        "GSE67462_platform.csv", "GSE67462_platform.tsv", "GSE67462_platform.txt",
    }
    found = []
    for base in roots:
        if not base.exists():
            continue
        for name in names:
            p = base / name
            if p.exists():
                found.append(p)
    return list(dict.fromkeys(found))


def _candidate_expression_sources(root: Path) -> list[Path]:
    candidates = [
        root / "results" / "GSE67462" / "03_expression_for_EDA.csv",
        root / "results" / "GSE67462" / "01_sample_level_counts.csv",
    ]
    return [p for p in candidates if p.exists()]


def _summarise_ids(ids: pd.Index) -> dict:
    vals = [str(x) for x in ids]
    refseq_like = sum(bool(re.fullmatch(r"(?:NM|NR|XM|XR)_\d+(?:\.\d+)?(?:_at)?", x, re.I)) for x in vals)
    symbol_like = sum(bool(re.fullmatch(r"[A-Za-z][A-Za-z0-9.-]*", x)) and not re.match(r"(?:NM|NR|XM|XR)_", x, re.I) for x in vals)
    return {
        "n": len(vals),
        "n_unique": len(set(vals)),
        "n_refseq_like": refseq_like,
        "fraction_refseq_like": refseq_like / len(vals) if vals else 0.0,
        "n_symbol_like": symbol_like,
        "fraction_symbol_like": symbol_like / len(vals) if vals else 0.0,
        "examples": vals[:20],
    }


def _build_mapping_report(expression_ids: pd.Index, tss: pd.DataFrame, platform: pd.DataFrame | None) -> tuple[pd.DataFrame, dict]:
    expr = pd.DataFrame({"expression_id": [str(x) for x in expression_ids]})
    expr["symbol_norm"] = expr.expression_id.map(_norm_symbol)

    tss_gene = tss["gene"].astype(str)
    tss_by_symbol: dict[str, list[str]] = {}
    for value in tss_gene:
        key = _norm_symbol(value)
        if key:
            tss_by_symbol.setdefault(key, []).append(value)
    tss_symbols = set(tss_by_symbol)
    expr["tss_symbol_match"] = expr.symbol_norm.map(lambda x: len(tss_by_symbol.get(x, [])))
    expr["tss_symbol_examples"] = expr.symbol_norm.map(lambda x: ";".join(tss_by_symbol.get(x, [])[:5]))

    platform_info = {"available": False}
    if platform is not None and not platform.empty:
        feature_col = _find_feature_column(platform)
        symbol_col = _find_symbol_column(platform)
        refseq_col = _find_refseq_column(platform)
        platform_info = {
            "available": feature_col is not None,
            "feature_column": feature_col,
            "symbol_column": symbol_col,
            "refseq_column": refseq_col,
            "description_column": "Description" if "Description" in platform.columns else None,
            "n_rows": int(len(platform)),
        }
        p = platform.copy()
        if feature_col:
            p["feature_norm"] = p[feature_col].map(_norm_refseq)
        else:
            p["feature_norm"] = ""
        if symbol_col:
            p["symbol_norm"] = p[symbol_col].map(_norm_symbol)
        elif "Description" in p.columns:
            p["symbol_norm"] = p["Description"].map(lambda x: _extract_description_symbol(x, tss_symbols))
        else:
            p["symbol_norm"] = ""
        if refseq_col:
            p["refseq_norm"] = p[refseq_col].map(_norm_refseq)
        else:
            p["refseq_norm"] = p["feature_norm"]
        p["symbol_norm"] = p["symbol_norm"].fillna("")
        p = p.drop_duplicates("feature_norm")
        by_symbol = p[p.symbol_norm != ""].drop_duplicates("symbol_norm")
        platform_symbols = set(by_symbol.symbol_norm)
        expr["platform_symbol_match"] = expr.symbol_norm.isin(platform_symbols).astype(int)
        expr["platform_symbol"] = expr.symbol_norm.where(expr.symbol_norm.isin(platform_symbols), "")
        expr["platform_symbol_tss_match"] = expr.platform_symbol.map(lambda x: int(bool(x) and x in tss_symbols))
        expr["platform_feature_match"] = pd.NA
        platform_info["n_symbols_recovered"] = int(len(platform_symbols))
        platform_info["n_symbols_matching_tss"] = int(len(platform_symbols & tss_symbols))
    else:
        expr["platform_symbol_match"] = pd.NA
        expr["platform_symbol"] = ""
        expr["platform_symbol_tss_match"] = pd.NA
        expr["platform_feature_match"] = pd.NA

    summary = {
        "expression_id_summary": _summarise_ids(expression_ids),
        "tss_gene_count": int(len(tss_gene)),
        "tss_unique_normalized_symbols": int(len(tss_by_symbol)),
        "expression_case_insensitive_symbol_overlap": int((expr.tss_symbol_match > 0).sum()),
        "expression_platform_symbol_matches": None if platform is None else int(expr.platform_symbol_match.fillna(0).sum()),
        "expression_platform_symbol_to_tss_matches": None if platform is None else int(expr.platform_symbol_tss_match.fillna(0).sum()),
        "platform": platform_info,
    }
    return expr, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", default="Data/GSE67520/mm9.refGene.gtf.gz")
    parser.add_argument("--platform-table", default=None, help="Optional local GPL19972 annotation table")
    parser.add_argument("--platform-soft", default=None, help="Local NCBI GEO GPL19972 SOFT family file (.soft or .soft.gz)")
    parser.add_argument("--output", default="results/Dynamics/z6_gse67462_identifier_mapping_audit")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)

    matrix, metadata = validation._load_common_space()
    tss = _load_tss(Path(args.gtf))
    expression_ids = pd.Index([str(x) for x in matrix.index])

    platform_path: Path | None = None
    platform = None
    if args.platform_soft:
        platform_path = Path(args.platform_soft)
        platform = _read_soft_platform(platform_path)
    elif args.platform_table:
        platform_path = Path(args.platform_table)
        platform = _read_table(platform_path)
    else:
        candidates = _candidate_platform_tables(root)
        platform_path = candidates[0] if candidates else None
        if platform_path is not None:
            if platform_path.name.lower().endswith(".soft") or platform_path.name.lower().endswith(".soft.gz"):
                platform = _read_soft_platform(platform_path)
            else:
                platform = _read_table(platform_path)

    mapping, summary = _build_mapping_report(expression_ids, tss, platform)
    summary.update({
        "expression_matrix_shape": [int(matrix.shape[0]), int(matrix.shape[1])],
        "expression_source_candidates": [str(p) for p in _candidate_expression_sources(root)],
        "gtf_path": str(Path(args.gtf)),
        "gtf_sha256": _sha256(Path(args.gtf)) if Path(args.gtf).exists() else None,
        "platform_table_path": str(platform_path) if platform_path else None,
        "platform_table_sha256": _sha256(platform_path) if platform_path else None,
        "platform_format": "GEO_SOFT" if args.platform_soft or (platform_path and platform_path.name.lower().endswith(".soft")) else "TABLE" if platform_path else None,
        "provenance_note": "GSE67462 is GPL19972 Brainarray MoGene10stv1_Mm_REFSEQ v18; GEO processed features are RefSeq-like IDs. The Stage 2.6 common-space matrix uses gene symbols, so platform provenance must be established before regulatory concordance is interpreted.",
    })

    mapping.to_csv(out / "01_feature_mapping_audit.csv", index=False)
    report = {
        "status": "COMPLETE" if platform is not None else "PARTIAL_PLATFORM_TABLE_MISSING",
        "summary": summary,
        "metadata_datasets": sorted(metadata.dataset.unique().tolist()),
    }
    (out / "02_mapping_report.json").write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")

    print("GSE67462 identifier mapping audit complete.")
    print(f"Expression matrix: {matrix.shape[0]} genes x {matrix.shape[1]} samples")
    print(f"Expression ID examples: {expression_ids[:15].tolist()}")
    print(f"TSS genes: {len(tss)}")
    print(f"Case-insensitive expression/TSS symbol overlap: {summary['expression_case_insensitive_symbol_overlap']}")
    if platform is None:
        print("GPL19972 platform annotation: NOT FOUND LOCALLY")
        print("Provide --platform-table or --platform-soft.")
    else:
        print(f"GPL19972 platform annotation: {platform_path}")
        print(f"Platform format: {summary['platform_format']}")
        print(f"Platform rows: {len(platform)}")
        print(f"Platform symbols recovered: {summary['platform']['n_symbols_recovered']}")
        print(f"Platform symbols matching TSS: {summary['platform']['n_symbols_matching_tss']}")
        print(f"Common-space expression symbols found in platform: {summary['expression_platform_symbol_matches']}")
        print(f"Common-space expression symbols mapped to TSS via platform: {summary['expression_platform_symbol_to_tss_matches']}")
    print(f"Outputs: {out}")


if __name__ == "__main__":
    main()
