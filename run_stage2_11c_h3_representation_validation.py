from __future__ import annotations

import csv
import gzip
import re

import pandas as pd

import dynamics.validate_h3_representation_overlap as validation


def _parse_gpl2670_robust():
    """Parse GEO GPL2670 .annot despite legacy column-name variants."""
    rows = []
    with gzip.open(validation.GPL2670, "rt", encoding="utf-8", errors="replace") as fh:
        in_table = False
        header = None
        probe_idx = None
        gb_idx = None
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            vals = next(csv.reader([line.rstrip("\n\r")], delimiter="\t"), [])
            if not vals:
                continue
            if header is None:
                header = [str(x).strip() for x in vals]
                normalized = [re.sub(r"[^a-z0-9]+", "", h.lower()) for h in header]
                probe_idx = next((i for i, h in enumerate(normalized) if h == "id"), None)
                gb_candidates = {"genbank", "genbankaccession", "gbacc", "gblist", "genbankaccessions"}
                gb_idx = next((i for i, h in enumerate(normalized) if h in gb_candidates), None)
                if probe_idx is None or gb_idx is None:
                    raise RuntimeError(
                        "GPL2670 annotation format changed: could not identify ID and GenBank columns. "
                        f"Columns: {header}"
                    )
                continue
            if len(vals) != len(header):
                continue
            probe = str(vals[probe_idx]).strip()
            raw_gb = str(vals[gb_idx]).strip()
            if not probe or not raw_gb:
                continue
            accessions = [x.strip() for x in re.split(r"[;|,\\s]+", raw_gb) if x.strip()]
            if accessions:
                rows.append((probe, accessions[0]))
    return pd.DataFrame(rows, columns=["probe_id", "genbank"])


# Replace only the fragile parser; the scientific gate remains unchanged.
validation.parse_gpl2670 = _parse_gpl2670_robust


if __name__ == "__main__":
    validation.main()
