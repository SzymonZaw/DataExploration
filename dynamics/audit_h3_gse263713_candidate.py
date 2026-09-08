"""Technical/prospective audit of GSE263713 as an orthogonal H3 temporal control.

No Yamanaka similarity, trajectory correlation, null test, or biological decision is
computed here. The audit only verifies design properties that were available before
H3 similarity analysis: time resolution, repeated-measures structure, sample IDs,
and expression-matrix compatibility.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
MATRIX = DATA / "GSE263713_raw_counts.tsv.gz"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def infer_sample_time(text: str):
    s = str(text).strip()
    m = re.search(r"-TP-([0-9]+(?:\\.[0-9]+)?)$", s, re.I)
    return float(m.group(1)) if m else None


def infer_individual(text: str):
    s = str(text).strip()
    m = re.match(r"(.+?)-TP-[0-9]+(?:\\.[0-9]+)?$", s, re.I)
    return m.group(1) if m else None


def _is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "ok",
        "accession": "GSE263713",
        "role": "prospective H3 orthogonal temporal-control candidate",
        "scope": "design/technical audit only; no Yamanaka similarity, trajectory agreement, null test, or biological decision",
        "source": "NCBI GEO",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263713",
        "download_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263713/suppl/GSE263713_raw_counts.tsv.gz",
        "expected_design": {
            "process_class": "circadian temporal regulation in primary human skin fibroblasts",
            "n_individuals": 6,
            "n_timepoints": 13,
            "time_hours": [0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48],
            "measurement": "RNA-seq",
            "repeated_measures_unit": "individual fibroblast cell line",
        },
    }

    if not MATRIX.exists():
        report["status"] = "pending_acquisition"
        report["decision"] = {"eligible": False, "reason": "GSE263713_raw_counts.tsv.gz is not present in Data/."}
        (OUT / "H3_GSE263713_CANDIDATE_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 0

    report["input"] = {"path": str(MATRIX.relative_to(ROOT)), "sha256": sha256(MATRIX)}

    # GSE263713 raw counts contain five gene-annotation columns before the samples:
    # Chr, Start, End, Strand, Length. Samples start at column index 5.
    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        header = next(csv.reader([fh.readline().rstrip("\n")], delimiter="\t"))
        first_rows = []
        for _ in range(20):
            line = fh.readline()
            if not line:
                break
            first_rows.append(next(csv.reader([line.rstrip("\n")], delimiter="\t")))

    gene_annotation_columns = 5
    if len(header) <= gene_annotation_columns:
        report["matrix"] = {"orientation": "gene_x_sample", "n_sample_columns": 0}
        report["repeated_measures"] = {"n_individuals_inferred": 0, "individuals": [], "n_timepoints_inferred": 0, "time_hours": [], "timepoint_counts": {}, "complete_time_courses": {}}
        report["decision"] = {"eligible_for_next_technical_gate": False, "reason": "UNRESOLVED: matrix header has fewer than five annotation columns plus one sample."}
        (OUT / "H3_GSE263713_CANDIDATE_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 0

    sample_ids = header[gene_annotation_columns:]
    numeric_values = []
    for row in first_rows:
        numeric_values.extend(row[gene_annotation_columns:])
    numeric = sum(1 for x in numeric_values if str(x).strip() and _is_number(x))
    nonempty = sum(1 for x in numeric_values if str(x).strip())

    parsed = [(infer_individual(s), infer_sample_time(s)) for s in sample_ids]
    individuals = sorted({x[0] for x in parsed if x[0] is not None})
    times = sorted({x[1] for x in parsed if x[1] is not None})
    counts = {}
    for individual, time in parsed:
        if individual is not None and time is not None:
            counts.setdefault(individual, []).append(time)

    report["matrix"] = {
        "orientation": "gene_x_sample",
        "gene_annotation_columns": gene_annotation_columns,
        "gene_annotation_examples": header[:gene_annotation_columns],
        "n_sample_columns": len(sample_ids),
        "sample_id_examples": sample_ids[:10],
        "numeric_fraction_first_20_rows": numeric / nonempty if nonempty else 0.0,
    }
    report["repeated_measures"] = {
        "n_individuals_inferred": len(individuals),
        "individuals": individuals,
        "n_timepoints_inferred": len(times),
        "time_hours": times,
        "timepoint_counts": {str(t): sum(1 for _, x in parsed if x == t) for t in times},
        "complete_time_courses": {str(i): sorted(v) == times for i, v in counts.items()},
    }

    expected_times = [float(x) for x in report["expected_design"]["time_hours"]]
    complete = (
        len(individuals) == 6
        and len(times) == 13
        and times == expected_times
        and all(sorted(v) == times for v in counts.values())
        and len(counts) == len(individuals)
    )
    numeric_ok = (numeric / nonempty) >= 0.95 if nonempty else False
    report["decision"] = {
        "eligible_for_next_technical_gate": bool(complete and numeric_ok),
        "reason": "PASS: repeated 48h temporal design with six inferred individual trajectories and 13 timepoints." if complete and numeric_ok else "UNRESOLVED: sample/time structure or numeric matrix content did not match the expected design.",
    }
    report["guardrail"] = "Eligibility does not imply H3 support; candidate selection is locked before any similarity statistic."

    (OUT / "H3_GSE263713_CANDIDATE_AUDIT.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
