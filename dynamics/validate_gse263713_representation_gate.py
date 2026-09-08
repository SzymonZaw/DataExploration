"""Locked representation/mapping gate for prospective H3 candidate GSE263713.

No Yamanaka similarity, trajectory agreement, null test, or H3 decision is
computed here. This gate verifies identifier namespace, canonical mapping,
collision handling, matrix/sample structure, and compatibility with the
PROGENy/DoRothEA representation used by Stage 2.11C.
"""
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
MATRIX = DATA / "GSE263713_raw_counts.tsv.gz"
EXPECTED_SHA = "8ce1e16a0039d93449e70aa6e827bbc8510a9b15dfbb78dc920cc2a2de5615a6"

MIN_MAPPING_RATE = 0.90
MIN_PROGENY_OVERLAP = 0.80
MIN_DOROTHEA_OVERLAP = 0.80


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def is_number(value: str) -> bool:
    try:
        float(value)
        return True
    except (TypeError, ValueError):
        return False


def classify_namespace(ids: list[str]) -> str:
    clean = [str(x).strip() for x in ids if str(x).strip()]
    if not clean:
        return "unknown"
    ensembl = sum(bool(re.fullmatch(r"ENSG\d+(?:\.\d+)?", x.upper())) for x in clean)
    numeric = sum(bool(re.fullmatch(r"\d+", x)) for x in clean)
    if ensembl / len(clean) >= 0.80:
        return "ensembl_gene"
    if numeric / len(clean) >= 0.80:
        return "entrez_gene"
    return "gene_symbol_or_other"


def map_to_hgnc(ids: list[str], namespace: str) -> tuple[dict[str, str], int]:
    clean = sorted(set(str(x).split(".", 1)[0].strip().upper() for x in ids if str(x).strip()))
    if namespace == "gene_symbol_or_other":
        mapping = {x: x for x in clean if re.fullmatch(r"[A-Z][A-Z0-9-]{1,14}", x)}
        return mapping, len(clean)

    try:
        import mygene
    except ImportError as exc:
        raise RuntimeError("Install mygene from requirements-ai.txt before running this stage.") from exc

    scopes = "ensembl.gene" if namespace == "ensembl_gene" else "entrezgene"
    info = mygene.MyGeneInfo()
    mapping: dict[str, str] = {}
    for i in range(0, len(clean), 1000):
        result = info.querymany(
            clean[i:i + 1000],
            scopes=scopes,
            fields="symbol",
            species="human",
            as_dataframe=False,
            returnall=False,
            step=1000,
            verbose=False,
        )
        for row in result:
            if row.get("notfound"):
                continue
            q = row.get("query")
            symbol = row.get("symbol")
            if q and symbol:
                mapping[str(q).upper()] = str(symbol).upper()
    return mapping, len(clean)


def network_targets() -> dict[str, set[str]]:
    import decoupler as dc

    progeny = dc.op.progeny(organism="human", top=100)
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
    return {
        "PROGENy": set(progeny["target"].astype(str).str.upper()),
        "DoRothEA": set(dorothea["target"].astype(str).str.upper()),
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "ok",
        "accession": "GSE263713",
        "scope": "representation/mapping technical gate only; no H3 similarity, trajectory agreement, null test, or biological decision",
        "input": {"path": str(MATRIX.relative_to(ROOT))},
        "thresholds": {
            "min_mapping_rate": MIN_MAPPING_RATE,
            "min_progeny_overlap": MIN_PROGENY_OVERLAP,
            "min_dorothea_overlap": MIN_DOROTHEA_OVERLAP,
        },
    }

    if not MATRIX.exists():
        report["status"] = "pending_acquisition"
        report["decision"] = {"eligible_for_h3_similarity": False, "reason": "Input matrix is missing."}
        print(json.dumps(report, indent=2))
        return 0

    actual_sha = sha256(MATRIX)
    report["input"].update({"sha256": actual_sha, "expected_sha256": EXPECTED_SHA, "sha_matches": actual_sha == EXPECTED_SHA})

    with gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace") as fh:
        header = next(csv.reader([fh.readline().rstrip("\n")], delimiter="\t"))
        first_rows = []
        all_gene_ids = []
        for line_no, line in enumerate(fh):
            if not line:
                break
            row = next(csv.reader([line.rstrip("\n")], delimiter="\t"), [])
            if not row:
                continue
            all_gene_ids.append(row[0])
            if len(first_rows) < 20:
                first_rows.append(row)

    # featureCounts output contains Geneid, Chr, Start, End, Strand, Length;
    # sample columns have the locked GSE263713 -TP-<hour> naming convention.
    sample_cols = [x for x in header if re.search(r"-TP-\d+(?:\D|$)", str(x))]
    annotation_cols = [x for x in header if x not in sample_cols]
    namespace = classify_namespace(all_gene_ids)

    numeric_values = []
    pos = {name: i for i, name in enumerate(header)}
    for row in first_rows:
        numeric_values.extend(row[pos[c]] for c in sample_cols if pos[c] < len(row))
    numeric = sum(1 for x in numeric_values if str(x).strip() and is_number(x))
    nonempty = sum(1 for x in numeric_values if str(x).strip())

    mapping, n_unique_input = map_to_hgnc(all_gene_ids, namespace)
    mapped_symbols = list(mapping.values())
    collision_count = len(mapped_symbols) - len(set(mapped_symbols))
    symbols = set(mapped_symbols)

    targets = network_targets()
    overlap = {}
    for name, net in targets.items():
        shared = symbols & net
        overlap[name] = {
            "n_network_targets": len(net),
            "n_shared_genes": len(shared),
            "network_overlap_fraction": len(shared) / len(net) if net else 0.0,
        }

    parsed = []
    for sample in sample_cols:
        m = re.match(r"^([^-]+)-TP-(\d+(?:\.\d+)?)$", str(sample))
        if m:
            parsed.append((m.group(1), float(m.group(2))))
    individuals = sorted({x[0] for x in parsed})
    times = sorted({x[1] for x in parsed})
    complete = {i: sorted(t for j, t in parsed if j == i) == times for i in individuals}

    report["matrix"] = {
        "orientation": "gene_x_sample",
        "gene_annotation_columns": annotation_cols,
        "n_gene_rows": len(all_gene_ids),
        "n_sample_columns": len(sample_cols),
        "sample_id_examples": sample_cols[:10],
        "numeric_fraction_first_20_rows": numeric / nonempty if nonempty else 0.0,
    }
    report["identifier_mapping"] = {
        "input_namespace": namespace,
        "n_unique_input_ids_audited": n_unique_input,
        "n_mapped_ids": len(mapping),
        "mapping_rate": len(mapping) / n_unique_input if n_unique_input else 0.0,
        "n_unique_hgnc_symbols": len(symbols),
        "mapping_collision_count": collision_count,
        "aggregation_rule": "duplicate canonical HGNC symbols must be averaged before network scoring",
        "mapping_method": "MyGene Ensembl/Entrez -> HGNC symbol; canonical HGNC-like symbols retained directly",
    }
    report["representation_overlap"] = overlap
    report["sample_structure"] = {
        "n_individuals": len(individuals),
        "individuals": individuals,
        "n_timepoints": len(times),
        "time_hours": times,
        "complete_time_courses": complete,
    }

    mapping_ok = report["identifier_mapping"]["mapping_rate"] >= MIN_MAPPING_RATE
    overlap_ok = all(
        overlap[name]["network_overlap_fraction"] >= threshold
        for name, threshold in [("PROGENy", MIN_PROGENY_OVERLAP), ("DoRothEA", MIN_DOROTHEA_OVERLAP)]
    )
    sample_ok = len(sample_cols) == 78 and len(individuals) == 6 and len(times) == 13 and all(complete.values())
    numeric_ok = (numeric / nonempty) >= 0.95 if nonempty else False
    hash_ok = actual_sha == EXPECTED_SHA

    report["decision"] = {
        "eligible_for_h3_similarity": bool(hash_ok and mapping_ok and overlap_ok and sample_ok and numeric_ok),
        "gates": {
            "input_hash": hash_ok,
            "identifier_mapping": mapping_ok,
            "representation_overlap": overlap_ok,
            "sample_structure": sample_ok,
            "numeric_matrix": numeric_ok,
        },
        "reason": "PASS: GSE263713 is compatible with the locked HGNC/network representation and repeated-measures design." if hash_ok and mapping_ok and overlap_ok and sample_ok and numeric_ok else "FAIL: representation/mapping compatibility gate not satisfied.",
    }
    report["guardrail"] = "Eligibility permits the next H3 similarity/null gate; it is not evidence for or against H3."

    out = OUT / "H3_GSE263713_REPRESENTATION_GATE.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"WROTE {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
