"""Stage 2.11C H3: locked representation-overlap/time/block validation.

This stage does NOT compute H3 similarity, trajectory agreement, null tests,
or biological decisions. It only establishes whether the two prospective
controls can enter the same representation contract as the Yamanaka analysis.

GSE3945 is a legacy spotted cDNA array (GPL2670), so probe IDs are mapped
through the locked GEO platform annotation and then to HGNC symbols.
GSE129486 uses human Ensembl gene IDs and is mapped directly to HGNC symbols.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"

GSE3945 = DATA / "GSE3945_series_matrix.txt.gz"
GSE129486_EXPR = DATA / "GSE129486_rnaseq-data-1_gene-tpm.tsv.gz"
GSE129486_META = DATA / "GSE129486_rnaseq-data-1_metadata.tsv.gz"
GPL2670 = DATA / "GPL2670.annot.gz"
MAPPING = OUT / "H3_GSE3945_probe_to_hgnc.tsv"
MAPPING_SHA = OUT / "H3_GSE3945_probe_to_hgnc.sha256"

GPL2670_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL2nnn/GPL2670/annot/GPL2670.annot.gz"
EXPECTED_GSE3945_SHA = "a79ec8c7b5bb88742a43bba7495856dd6a3d3e82b802588ea58d9773ef0598eb"
EXPECTED_GSE129486_EXPR_SHA = "6e3d7860f4f38d95830a15b8dd570d58f9226170df6002343e432a05c96fcbf0"
EXPECTED_GSE129486_META_SHA = "78a60e353461ab819672e478a9f38b20822fdfc493a1773677ddef3eb167ff04"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def infer_time(value):
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    s = str(value).strip().lower()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        pass
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:h|hr|hrs|hour|hours)\b", s)
    if m:
        return float(m.group(1))
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*(?:min|mins|minute|minutes)\b", s)
    if m:
        return float(m.group(1)) / 60.0
    m = re.search(r"([0-9]+)\s+([0-9]+)\s*/\s*([0-9]+)\s*(?:h|hr|hrs|hour|hours)\b", s)
    if m:
        return float(m.group(1)) + float(m.group(2)) / float(m.group(3))
    return 0.0 if re.search(r"(?:^|[^0-9])0(?:\s*(?:h|hr|hrs|hour|hours))?(?:[^0-9]|$)", s) else None


def parse_geo_matrix(path: Path):
    samples, titles, data = [], [], []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                samples = next(csv.reader([line[len("!Sample_geo_accession"):].lstrip("\t =")], delimiter="\t"), [])
            elif line.startswith("!Sample_title"):
                titles = next(csv.reader([line[len("!Sample_title"):].lstrip("\t =")], delimiter="\t"), [])
            elif line.startswith("!series_matrix_table_begin"):
                for row in fh:
                    if row.startswith("!series_matrix_table_end"):
                        break
                    vals = next(csv.reader([row.rstrip("\n")], delimiter="\t"), [])
                    if vals:
                        data.append(vals)
                break
    header = data[0]
    rows = data[1:]
    expr = pd.DataFrame(rows, columns=header)
    expr = expr.set_index(header[0])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr.columns = [str(x) for x in expr.columns]
    times = [infer_time(t) for t in titles]
    return expr, samples, titles, times


def download_gpl_if_missing():
    if GPL2670.exists():
        return
    print(f"Downloading {GPL2670_URL}")
    urllib.request.urlretrieve(GPL2670_URL, GPL2670)


def parse_gpl2670():
    """Return probe -> GenBank mapping from the locked GPL2670 annotation."""
    rows = []
    with gzip.open(GPL2670, "rt", encoding="utf-8", errors="replace") as fh:
        in_table = False
        header = None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            vals = next(csv.reader([line.rstrip("\n")], delimiter="\t"), [])
            if header is None:
                header = vals
                continue
            if len(vals) != len(header):
                continue
            row = dict(zip(header, vals))
            probe = str(row.get("ID", "")).strip()
            gb = str(row.get("GenBank", "")).strip()
            if probe and gb:
                rows.append((probe, gb))
    return pd.DataFrame(rows, columns=["probe_id", "genbank"])


def map_genbank_to_hgnc(accessions: list[str]) -> dict[str, str]:
    try:
        import mygene
    except ImportError as exc:
        raise RuntimeError("Install mygene from requirements-ai.txt before running this stage.") from exc
    info = mygene.MyGeneInfo()
    mapping = {}
    unique = sorted(set(x for x in accessions if x))
    for i in range(0, len(unique), 1000):
        batch = unique[i:i + 1000]
        result = info.querymany(
            batch,
            scopes="accession",
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
            q = str(row.get("query", "")).strip()
            symbol = row.get("symbol")
            if q and symbol:
                mapping[q] = str(symbol).upper()
    return mapping


def map_ensembl_to_hgnc(ids: list[str]) -> dict[str, str]:
    try:
        import mygene
    except ImportError as exc:
        raise RuntimeError("Install mygene from requirements-ai.txt before running this stage.") from exc
    info = mygene.MyGeneInfo()
    clean = sorted(set(str(x).split(".", 1)[0].upper() for x in ids if str(x).strip()))
    mapping = {}
    for i in range(0, len(clean), 1000):
        result = info.querymany(
            clean[i:i + 1000],
            scopes="ensembl.gene",
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
            q, symbol = row.get("query"), row.get("symbol")
            if q and symbol:
                mapping[str(q).upper()] = str(symbol).upper()
    return mapping


def network_overlap(symbols: set[str]):
    import decoupler as dc
    progeny = dc.op.progeny(organism="human", top=100)
    dorothea = dc.op.dorothea(organism="human", levels=["A", "B", "C"])
    out = {}
    for name, net in [("PROGENy", progeny), ("DoRothEA", dorothea)]:
        targets = set(net["target"].astype(str).str.upper())
        out[name] = {
            "n_network_targets": len(targets),
            "n_shared_genes": len(symbols & targets),
            "network_overlap_fraction": len(symbols & targets) / len(targets) if targets else 0.0,
        }
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report = {
        "status": "ok",
        "scope": "representation-overlap/time/block validation only; no H3 similarity, trajectory agreement, null test, or biological decision",
        "controls": {},
        "mapping_artifacts": {},
        "representation_overlap": {},
        "time_block_gate": {},
    }

    # Confirm the already locked inputs before doing any mapping work.
    hashes = {
        "GSE3945": (sha256(GSE3945), EXPECTED_GSE3945_SHA),
        "GSE129486_gene_tpm": (sha256(GSE129486_EXPR), EXPECTED_GSE129486_EXPR_SHA),
        "GSE129486_metadata": (sha256(GSE129486_META), EXPECTED_GSE129486_META_SHA),
    }
    report["input_hashes"] = {k: {"sha256": a, "expected": b, "matches": a == b} for k, (a, b) in hashes.items()}

    # GSE3945: legacy spotted cDNA array. GEO identifies it as GPL2670 (SHAT).
    # Its VALUE matrix is already log2 ratios, so do not apply RNA-seq log1p/CPM.
    download_gpl_if_missing()
    g_expr, g_samples, g_titles, g_times = parse_geo_matrix(GSE3945)
    gpl_sha = sha256(GPL2670)
    gpl = parse_gpl2670()
    gpl_map = map_genbank_to_hgnc(gpl["genbank"].tolist())
    gpl["hgnc_symbol"] = gpl["genbank"].map(gpl_map)
    gpl = gpl.drop_duplicates("probe_id")
    gpl.to_csv(MAPPING, sep="\t", index=False)
    MAPPING_SHA.write_text(f"{sha256(MAPPING)}  {MAPPING.name}\n", encoding="utf-8")

    probe_to_symbol = dict(zip(gpl["probe_id"], gpl["hgnc_symbol"]))
    mapped_symbols = pd.Series(g_expr.index.astype(str)).map(probe_to_symbol).dropna().astype(str).str.upper()
    g_symbols = set(mapped_symbols)
    report["controls"]["GSE3945"] = {
        "platform": "GPL2670",
        "platform_title": "SHAT",
        "n_samples": len(g_samples),
        "n_timepoints": len(set(x for x in g_times if x is not None)),
        "time_values_hours": sorted(set(x for x in g_times if x is not None)),
        "n_probes": int(len(g_expr)),
        "n_mapped_hgnc_rows": int(mapped_symbols.size),
        "n_unique_hgnc_symbols": int(len(g_symbols)),
        "mapping_rate": float(mapped_symbols.size / len(g_expr)) if len(g_expr) else 0.0,
        "missing_value_fraction": float(g_expr.isna().to_numpy().mean()),
    }
    report["mapping_artifacts"]["GPL2670"] = {
        "url": GPL2670_URL,
        "local_path": str(GPL2670.relative_to(ROOT)),
        "sha256": gpl_sha,
        "mapping_path": str(MAPPING.relative_to(ROOT)),
        "mapping_sha256": sha256(MAPPING),
        "mapping_method": "GPL2670 GenBank -> MyGene accession -> HGNC symbol",
    }

    # GSE129486: Ensembl -> HGNC; metadata supplies time, cell line and stimulation.
    e = pd.read_csv(GSE129486_EXPR, sep="\t", compression="gzip", low_memory=False)
    m = pd.read_csv(GSE129486_META, sep="\t", compression="gzip", low_memory=False)
    gene_col = e.columns[0]
    ens = e[gene_col].astype(str).str.split(".", n=1).str[0].str.upper().tolist()
    emap = map_ensembl_to_hgnc(ens)
    e_symbols = set(emap.values())
    report["controls"]["GSE129486"] = {
        "platform": "GPL11154",
        "n_samples": int(e.shape[1] - 1),
        "n_timepoints": int(m["time"].nunique()),
        "time_values_hours": sorted(pd.to_numeric(m["time"], errors="coerce").dropna().unique().tolist()),
        "n_input_ensembl_ids": int(len(ens)),
        "n_mapped_hgnc_ids": int(len(emap)),
        "n_unique_hgnc_symbols": int(len(e_symbols)),
        "mapping_rate": float(len(emap) / len(set(ens))) if set(ens) else 0.0,
        "metadata_columns": list(m.columns),
        "n_cell_lines": int(m["cell_line"].nunique()),
        "n_stimulations": int(m["stimulation"].nunique()),
    }

    report["representation_overlap"]["GSE3945"] = network_overlap(g_symbols)
    report["representation_overlap"]["GSE129486"] = network_overlap(e_symbols)

    # Block/time audit: H3 may not proceed with a null requiring exchangeability
    # if the control has essentially one biological observation per time point.
    g_time_counts = pd.Series(g_times).value_counts().sort_index()
    g_replicated_times = int((g_time_counts >= 2).sum())
    m_time_counts = pd.to_numeric(m["time"], errors="coerce").value_counts().sort_index()
    report["time_block_gate"] = {
        "GSE3945": {
            "replicated_timepoints_n": g_replicated_times,
            "timepoint_counts": {str(k): int(v) for k, v in g_time_counts.items()},
            "exchangeability_warning": "Most time points have one sample; this is insufficient for a sample-level block permutation within time.",
        },
        "GSE129486": {
            "replicated_timepoints_n": int((m_time_counts >= 2).sum()),
            "timepoint_counts": {str(k): int(v) for k, v in m_time_counts.items()},
            "block_columns": [c for c in ["cell_line", "stimulation", "plate"] if c in m.columns],
        },
    }
    report["decision"] = {
        "representation_overlap_gate": "PASS" if all(v.get("PROGENy", {}).get("n_shared_genes", 0) > 0 and v.get("DoRothEA", {}).get("n_shared_genes", 0) > 0 for v in report["representation_overlap"].values()) else "FAIL",
        "time_block_gate": "UNRESOLVED",
        "reason": "GSE3945 has only three baseline replicates and mostly single observations at other time points; the locked H3 null must not treat these singleton time points as exchangeable biological replicates.",
        "h3_similarity_allowed": False,
    }
    report["guardrail"] = "No H3 similarity or biological conclusion is computed. The output only determines whether the controls are compatible with the locked representation and whether the declared block/null structure is executable."

    out = OUT / "H3_REPRESENTATION_OVERLAP_VALIDATION.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"WROTE {out.relative_to(ROOT)}")
    print("No H3 similarity analysis was performed.")


if __name__ == "__main__":
    main()
