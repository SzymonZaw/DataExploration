from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"

FILES = {
    "GSE3945": DATA / "GSE3945_series_matrix.txt.gz",
    "GSE129486_gene_tpm": DATA / "GSE129486_rnaseq-data-1_gene-tpm.tsv.gz",
    "GSE129486_metadata": DATA / "GSE129486_rnaseq-data-1_metadata.tsv.gz",
}
EXPECTED_SHA256 = {
    "GSE3945": "a79ec8c7b5bb88742a43bba7495856dd6a3d3e82b802588ea58d9773ef0598eb",
    "GSE129486_gene_tpm": "6e3d7860f4f38d95830a15b8dd570d58f9226170df6002343e432a05c96fcbf0",
    "GSE129486_metadata": "78a60e353461ab819672e478a9f38b20822fdfc493a1773677ddef3eb167ff04",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def status(ok: bool, reason: str = "") -> dict:
    return {"status": "pass" if ok else "fail", "reason": reason}


def infer_time(value):
    """Infer time in hours from common GEO metadata/title formats."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip().lower()
    if not s:
        return None

    # Numeric metadata such as GSE129486's time column is already in hours.
    try:
        return float(s)
    except ValueError:
        pass

    # Explicit hours, including decimal values.
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:h|hr|hrs|hour|hours)\b", s)
    if m:
        return float(m.group(1))

    # Explicit minutes, converted to hours.
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:min|mins|minute|minutes)\b", s)
    if m:
        return float(m.group(1)) / 60.0

    # Fractional hour written as e.g. "1 1/2 hr".
    m = re.search(r"([0-9]+)\s+([0-9]+)\s*/\s*([0-9]+)\s*(?:h|hr|hrs|hour|hours)\b", s)
    if m:
        return float(m.group(1)) + float(m.group(2)) / float(m.group(3))

    # Day-based labels, converted to hours.
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:day|days|d)\b", s)
    if m:
        return float(m.group(1)) * 24.0

    # Bare zero labels such as "0 hr." or "0" embedded in text.
    if re.search(r"(?:^|[^0-9])0(?:\s*(?:h|hr|hrs|hour|hours))?(?:[^0-9]|$)", s):
        return 0.0
    return None


def parse_geo_annotation_line(line: str, key: str) -> list[str]:
    """Parse a GEO series-matrix annotation line after its key.

    GEO series matrix files use tab-separated metadata lines such as
    ``!Sample_title\t"sample 1"\t"sample 2"``; they do not use ``=`` here.
    """
    payload = line[len(key):].lstrip("\t =")
    if not payload:
        return []
    return next(csv.reader([payload], delimiter="\t"), [])


def parse_geo_series_matrix(path: Path) -> dict:
    samples = []
    titles = []
    characteristics = []
    data_lines = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                samples = parse_geo_annotation_line(line, "!Sample_geo_accession")
            elif line.startswith("!Sample_title"):
                titles = parse_geo_annotation_line(line, "!Sample_title")
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = parse_geo_annotation_line(line, "!Sample_characteristics_ch1")
                characteristics.append(vals)
            elif line.startswith("!series_matrix_table_begin"):
                data_lines = list(fh)
                break
    if not data_lines:
        return {"status": "fail", "reason": "GEO series matrix table not found"}

    # GEO series-matrix expression tables are tab-delimited.
    header = next(csv.reader([data_lines[0].rstrip("\n")], delimiter="\t"), [])
    rows = []
    for line in data_lines[1:]:
        if line.startswith("!series_matrix_table_end"):
            break
        vals = next(csv.reader([line.rstrip("\n")], delimiter="\t"), [])
        if vals:
            rows.append(vals)

    ncols = len(header)
    lengths_ok = all(len(r) == ncols for r in rows[:1000])
    sample_cols = header[1:]
    gene_col = header[0] if header else None
    id_values = [r[0] for r in rows if r]
    numeric_fraction = 0.0
    if rows and len(rows[0]) > 1:
        numeric = 0
        total = 0
        for r in rows[:200]:
            for x in r[1:]:
                total += 1
                try:
                    float(x)
                    numeric += 1
                except ValueError:
                    pass
        numeric_fraction = numeric / total if total else 0.0

    time_values = []
    for i, sid in enumerate(samples):
        text_parts = [sid]
        if i < len(titles):
            text_parts.append(titles[i])
        for ch in characteristics:
            if i < len(ch):
                text_parts.append(ch[i])
        time_values.append(infer_time(" ".join(text_parts)))

    structure_ok = bool(
        lengths_ok
        and gene_col
        and len(sample_cols) > 0
        and numeric_fraction > 0.95
    )
    return {
        "status": "pass" if structure_ok else "fail",
        "reason": "" if structure_ok else "matrix structure could not be validated",
        "matrix_orientation": "gene_x_sample" if structure_ok else "unresolved",
        "n_samples": len(sample_cols),
        "n_rows": len(rows),
        "gene_identifier_column": gene_col,
        "sample_ids": sample_cols[:20],
        "sample_geo_accessions": samples,
        "sample_titles": titles,
        "time_inference": {
            "n_inferred": sum(x is not None for x in time_values),
            "values": time_values,
            "exact_labels_available": len(samples) == len(time_values) and all(x is not None for x in time_values),
        },
        "numeric_fraction_first_200_rows": numeric_fraction,
        "row_identifier_sample": id_values[:10],
    }


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)


def infer_gene_column(df: pd.DataFrame):
    preferred = ["gene", "gene_id", "geneid", "symbol", "ensembl_gene_id", "gene_name"]
    lower = {str(c).lower(): c for c in df.columns}
    for p in preferred:
        if p in lower:
            return lower[p]
    candidates = []
    for c in df.columns:
        s = df[c].dropna().astype(str).head(100)
        if len(s) and (s.str.match(r"^ENSG\d+").mean() > 0.5 or s.str.match(r"^[A-Za-z0-9_.-]+$").mean() > 0.9):
            candidates.append(c)
    return candidates[0] if candidates else None


def validate_gse129486(expr_path: Path, meta_path: Path) -> dict:
    try:
        df = read_tsv(expr_path)
        meta = read_tsv(meta_path)
    except Exception as exc:
        return {"status": "fail", "reason": f"TSV parse failed: {type(exc).__name__}: {exc}"}
    gene_col = infer_gene_column(df)
    non_gene = [c for c in df.columns if c != gene_col]
    numeric = 0
    total = 0
    for c in non_gene[: min(100, len(non_gene))]:
        vals = pd.to_numeric(df[c].head(200), errors="coerce")
        numeric += int(vals.notna().sum())
        total += len(vals)
    numeric_fraction = numeric / total if total else 0.0

    meta_lower = {str(c).lower(): c for c in meta.columns}
    sample_key = next((meta_lower[x] for x in ["sample", "sample_id", "sampleid", "geo_accession", "gsm", "cell"] if x in meta_lower), None)
    time_cols = [c for c in meta.columns if re.search(r"time|hour|hr|day|hpi|dpi", str(c), re.I)]
    inferred_times = []
    if time_cols:
        for _, row in meta.iterrows():
            value = None
            for c in time_cols:
                value = infer_time(row[c])
                if value is not None:
                    break
            inferred_times.append(value)

    expr_samples = set(map(str, non_gene))
    meta_samples = set(map(str, meta[sample_key])) if sample_key else set()
    overlap = len(expr_samples & meta_samples)
    return {
        "status": "pass" if gene_col and numeric_fraction > 0.95 and overlap > 0 else "fail",
        "reason": "" if gene_col and numeric_fraction > 0.95 and overlap > 0 else "expression/metadata structure could not be validated",
        "matrix_orientation": "gene_x_sample" if gene_col else "unresolved",
        "n_expression_rows": len(df),
        "n_expression_columns": len(df.columns),
        "n_expression_samples": len(non_gene),
        "gene_identifier_column": gene_col,
        "gene_identifier_sample": df[gene_col].dropna().astype(str).head(10).tolist() if gene_col else [],
        "numeric_fraction_first_200_rows": numeric_fraction,
        "metadata_rows": len(meta),
        "metadata_columns": list(map(str, meta.columns)),
        "metadata_sample_key": sample_key,
        "expression_metadata_sample_overlap": overlap,
        "time_columns": list(map(str, time_cols)),
        "time_inference": {
            "n_inferred": sum(x is not None for x in inferred_times),
            "n_metadata_rows": len(inferred_times),
            "exact_labels_available": bool(inferred_times) and all(x is not None for x in inferred_times),
            "values_sample": inferred_times[:30],
        },
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "ok",
        "scope": "technical validation only; no H3 similarity, trajectory agreement, null test, or biological decision computed",
        "branch": "stage-2-11c-h3-control-redesign",
        "files": {},
        "datasets": {},
        "representation_gate": {
            "same_representation_as_yamanaka": "not_yet_computed",
            "network_overlap": "not_yet_computed",
            "reason": "This gate intentionally stops before any H3 similarity analysis; representation overlap is recorded by the subsequent H3 execution script after technical validation passes.",
        },
    }

    failures = []
    for key, path in FILES.items():
        entry = {"path": str(path.relative_to(ROOT)), "exists": path.exists()}
        if not path.exists():
            entry.update(status(False, "file missing"))
            failures.append(key)
        else:
            actual = sha256(path)
            expected = EXPECTED_SHA256.get(key)
            entry.update({"sha256": actual, "expected_sha256": expected, "sha256_matches_lock": actual == expected})
            if actual != expected:
                entry.update(status(False, "SHA-256 does not match acquisition lock"))
                failures.append(key)
            else:
                entry.update(status(True))
        report["files"][key] = entry

    if not failures:
        g = parse_geo_series_matrix(FILES["GSE3945"])
        report["datasets"]["GSE3945"] = g
        if g.get("status") != "pass":
            failures.append("GSE3945_structure")
        r = validate_gse129486(FILES["GSE129486_gene_tpm"], FILES["GSE129486_metadata"])
        report["datasets"]["GSE129486"] = r
        if r.get("status") != "pass":
            failures.append("GSE129486_structure")

    report["technical_gate"] = "PASS" if not failures else "FAIL"
    report["failures"] = failures
    if failures:
        report["next_step"] = "Repair or document the failed technical condition; do not compute H3 similarity."
    else:
        report["next_step"] = "Both controls passed the file and matrix-structure gate. Proceed to the locked representation-overlap/time/block validation stage before H3 similarity."

    out_json = OUT / "H3_TECHNICAL_VALIDATION.json"
    out_json.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"WROTE {out_json.relative_to(ROOT)}")
    print("No H3 similarity analysis was performed.")
    raise SystemExit(0 if not failures else 2)


if __name__ == "__main__":
    main()
