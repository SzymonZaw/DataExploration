from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from urllib.request import urlopen

GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242421"
SAMPLES = {
    "GSM7763395": 0,
    "GSM7763396": 2,
    "GSM7763397": 4,
    "GSM7763398": 6,
    "GSM7763399": 8,
    "GSM7763400": 10,
    "GSM7763401": 12,
    "GSM7763402": 14,
    "GSM7763403": "iPSC",
}


def fetch_text(url: str) -> str:
    with urlopen(url, timeout=30) as r:
        return r.read().decode("utf-8", errors="replace")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", default="results/Dynamics/z4_gse242421_preflight/preflight.json")
    ap.add_argument("--skip-network", action="store_true")
    args = ap.parse_args()

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    result = {
        "candidate": "GSE242421",
        "series_url": GEO,
        "required_samples": SAMPLES,
        "raw_tar_gb": 17.3,
        "raw_download_required": False,
        "raw_data_available": False,
        "processed_supplementary_data_available": True,
        "preferred_source": "analysis_products_or_processed_fragment/gene-activity data",
        "orthogonal_modality": "scATAC",
        "assembly": "hg38",
        "timepoints_days": [0, 2, 4, 6, 8, 10, 12, 14],
        "endpoint": "iPSC",
        "multiome_companion": "GSE242419",
        "scRNA_companion": "GSE242423",
        "status": "PREFLIGHT_ONLY",
    }
    if not args.skip_network:
        try:
            html = fetch_text(GEO)
            result["geo_reachable"] = True
            result["geo_contains_processed_data"] = "Processed data provided as supplementary file" in html
            result["geo_contains_raw_data_not_provided"] = "RAW DATA not provided" in html
            result["geo_contains_hg38"] = "hg38" in html
            result["geo_contains_9_samples"] = all(x in html for x in SAMPLES)
        except Exception as e:
            result["geo_reachable"] = False
            result["geo_error"] = f"{type(e).__name__}: {e}"
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print("GSE242421 Z4 ORTHOGONAL PREFLIGHT")
    print(f"samples required: {len(SAMPLES)}")
    print("timepoints: D0,D2,D4,D6,D8,D10,D12,D14 + iPSC")
    print("assembly: hg38")
    print("raw 17.3 Gb tar: NOT REQUIRED")
    print("processed/analysis products: REQUIRED NEXT")
    print(f"output: {out}")


if __name__ == "__main__":
    main()
