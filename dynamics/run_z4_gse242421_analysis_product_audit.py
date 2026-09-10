from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

EXPECTED_SCATAC = {
    "cells.tsv",
    "peaks.bed",
    "features.tsv",
    "cell_x_peak.mtx.gz",
}
EXPECTED_INTEGRATION = {
    "peak_gene_links_fdr1e-4.tsv",
    "harmony.cca.30.feat.tsv",
    "harmony.cca.metadata.tsv",
}

SOURCE = "https://zenodo.org/records/8313962"
EXPECTED_SAMPLES = ["D0", "D2", "D4", "D6", "D8", "D10", "D12", "D14", "iPSC"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def archive_inventory(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as zf:
        return zf.namelist()


def basename_set(names: list[str]) -> set[str]:
    return {Path(n).name for n in names if not n.endswith("/")}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scatac-zip", default="Data/GSE242421/scATAC.zip")
    ap.add_argument(
        "--integration-zip",
        default="Data/GSE242421/scATAC_scRNA_integration.zip",
    )
    ap.add_argument(
        "--output",
        default="results/Dynamics/z4_gse242421_analysis_product_audit/audit.json",
    )
    args = ap.parse_args()

    scatac = Path(args.scatac_zip)
    integration = Path(args.integration_zip)
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)

    result = {
        "candidate": "GSE242421",
        "source": SOURCE,
        "assembly": "hg38",
        "expected_samples": EXPECTED_SAMPLES,
        "raw_17_3gb_tar_required": False,
        "analysis_product_strategy": "scATAC peak-by-cell counts + peaks + independent peak-gene links; derive gene/module activity without refitting frozen modules",
        "archives": {},
        "gene_activity_matrix_present": False,
        "peak_matrix_present": False,
        "peak_gene_links_present": False,
        "decision": "UNRESOLVED",
    }

    for label, path, expected in [
        ("scATAC", scatac, EXPECTED_SCATAC),
        ("scATAC_scRNA_integration", integration, EXPECTED_INTEGRATION),
    ]:
        entry = {"path": str(path), "exists": path.exists()}
        if path.exists():
            entry["size_bytes"] = path.stat().st_size
            entry["sha256"] = sha256(path)
            try:
                names = archive_inventory(path)
                bases = basename_set(names)
                entry["member_count"] = len(names)
                entry["required_members"] = sorted(expected)
                entry["missing_members"] = sorted(expected - bases)
                entry["required_members_present"] = not (expected - bases)
                entry["members_preview"] = names[:50]
                if label == "scATAC":
                    result["peak_matrix_present"] = "cell_x_peak.mtx.gz" in bases
                if label == "scATAC_scRNA_integration":
                    result["peak_gene_links_present"] = "peak_gene_links_fdr1e-4.tsv" in bases
            except (OSError, zipfile.BadZipFile) as exc:
                entry["archive_error"] = f"{type(exc).__name__}: {exc}"
        result["archives"][label] = entry

    scatac_ok = result["archives"]["scATAC"].get("required_members_present", False)
    integ_ok = result["archives"]["scATAC_scRNA_integration"].get("required_members_present", False)
    if scatac_ok and integ_ok:
        result["decision"] = "READY_FOR_GENE_ACTIVITY_DERIVATION"
    elif scatac.exists() or integration.exists():
        result["decision"] = "PARTIAL_PROVENANCE_READY"

    out.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("GSE242421 Z4 ANALYSIS-PRODUCT AUDIT")
    print(f"scATAC.zip: {'FOUND' if scatac.exists() else 'MISSING'}")
    print(f"scATAC_scRNA_integration.zip: {'FOUND' if integration.exists() else 'MISSING'}")
    print(f"peak-by-cell matrix: {result['peak_matrix_present']}")
    print(f"peak-gene links: {result['peak_gene_links_present']}")
    print("gene-activity matrix present: False (will be derived, not fitted)")
    print(f"decision: {result['decision']}")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
