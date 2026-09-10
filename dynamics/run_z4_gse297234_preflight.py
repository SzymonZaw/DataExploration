"""Preflight audit for GSE297234 before frozen Z4 transfer.

This script inspects local RDS files only for structural metadata. It does not
compute a transfer score, fit a model, or select genes from validation results.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rds", nargs="+", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    paths = [Path(p) for p in args.rds]
    missing = [str(p) for p in paths if not p.exists()]
    result = {
        "candidate": "GSE297234",
        "rds_files": [str(p) for p in paths],
        "n_rds_files": len(paths),
        "missing_files": missing,
        "all_files_present": not missing,
        "transfer_score_evaluated": False,
        "orthology_mapping_evaluated": False,
        "endpoint_evaluated": False,
        "next_step": "inspect RDS structure and construct deterministic human_mouse orthology input" if not missing else "READY_FOR_RDS_STRUCTURE_INSPECTION",
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
