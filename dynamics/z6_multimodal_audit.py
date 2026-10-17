"""Audit candidate expression/regulatory files for the Z6 multimodal test.

Run from the repository root:
    python dynamics/z6_multimodal_audit.py
Optional:
    python dynamics/z6_multimodal_audit.py --root . --out results/Dynamics/z6_audit

This is an input/schema audit, not a biological result. It deliberately does not
join datasets by approximate time labels or assume that GSM identifiers match
between separate GEO series. A valid M0/M1/M2 comparison requires an explicit,
defensible sample/time mapping and a common gene universe.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

GSE_TARGETS = ("GSE67462", "GSE67520")
TABLE_SUFFIXES = {".csv", ".tsv", ".txt", ".tab"}
GSM_RE = re.compile(r"\bGSM\d+\b", re.IGNORECASE)
TIME_PATTERNS = (
    (re.compile(r"\b(?:day|d)\s*0\b", re.I), 0.0),
    (re.compile(r"\b(?:day|d)\s*1\b", re.I), 24.0),
    (re.compile(r"\b(?:day|d)\s*2\b", re.I), 48.0),
    (re.compile(r"\b(?:day|d)\s*3\b", re.I), 72.0),
    (re.compile(r"\b(?:day|d)\s*5\b", re.I), 120.0),
    (re.compile(r"\b(?:day|d)\s*7\b", re.I), 168.0),
    (re.compile(r"\b(?:day|d)\s*10\b", re.I), 240.0),
    (re.compile(r"\b(?:day|d)\s*11\b", re.I), 264.0),
    (re.compile(r"\b(?:day|d)\s*15\b", re.I), 360.0),
    (re.compile(r"\b(?:day|d)\s*18\b", re.I), 432.0),
)


def infer_time(label: str) -> float | None:
    """Infer only explicit day labels; do not infer time from GSM numbers."""
    text = str(label).replace("_", " ").replace("-", " ")
    for pattern, hours in TIME_PATTERNS:
        if pattern.search(text):
            return hours
    return None


def read_table(path: Path) -> tuple[pd.DataFrame | None, str | None]:
    # Try common delimiters; only a small preview is needed for schema discovery.
    for sep in (",", "\t", ";"):
        try:
            frame = pd.read_csv(path, sep=sep, nrows=1000, low_memory=False)
            if frame.shape[1] > 1:
                return frame, None
        except (UnicodeDecodeError, pd.errors.ParserError, OSError, ValueError):
            continue
    return None, "Could not parse as a delimited text table"


def gsm_ids(values) -> set[str]:
    found: set[str] = set()
    for value in values:
        if pd.isna(value):
            continue
        found.update(x.upper() for x in GSM_RE.findall(str(value)))
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output directory (default: results/Dynamics/z6_audit under root)",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    out = (args.out or root / "results" / "Dynamics" / "z6_audit").resolve()
    out.mkdir(parents=True, exist_ok=True)

    results_root = root / "results"
    if not results_root.exists():
        raise SystemExit(f"Missing results directory: {results_root}")

    candidates = []
    for path in results_root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TABLE_SUFFIXES:
            continue
        name = str(path.relative_to(root)).upper()
        matched = [gse for gse in GSE_TARGETS if gse in name]
        # Keep target-study files and generically named expression/regulatory tables
        # under either target's directory.
        if not matched and not any(
            token in path.name.lower()
            for token in ("expression", "h3k27", "h3k4", "rnapii", "oct4", "chip", "peak", "regulat")
        ):
            continue
        frame, error = read_table(path)
        rel = str(path.relative_to(root))
        if frame is None:
            candidates.append({
                "path": rel, "target_gse_in_path": ";".join(matched),
                "read_status": "unreadable", "error": error,
            })
            continue

        columns = [str(c) for c in frame.columns]
        header_gsm = gsm_ids(columns)
        value_gsm = gsm_ids(frame.astype(str).to_numpy().ravel())
        all_gsm = header_gsm | value_gsm
        time_labels = [c for c in columns if infer_time(c) is not None]
        candidates.append({
            "path": rel,
            "target_gse_in_path": ";".join(matched),
            "read_status": "ok",
            "preview_rows": len(frame),
            "n_columns": len(columns),
            "first_columns": " | ".join(columns[:12]),
            "GSM_ids_in_headers": len(header_gsm),
            "GSM_ids_in_preview_cells": len(value_gsm),
            "GSM_ids_in_preview_union": ";".join(sorted(all_gsm)[:100]),
            "explicit_day_columns": ";".join(time_labels),
            "likely_expression_or_regulatory": any(
                token in path.name.lower()
                for token in ("expression", "h3k27", "h3k4", "rnapii", "oct4", "chip", "peak", "regulat")
            ),
        })

    inventory = pd.DataFrame(candidates)
    inventory_path = out / "01_file_inventory.csv"
    inventory.to_csv(inventory_path, index=False)

    # Exact GSM overlap is only a diagnostic. Different GSM IDs do not rule out
    # time-matched experiments, and matching times alone do not establish pairing.
    gsm_by_file: dict[str, set[str]] = {}
    for item in candidates:
        if item.get("read_status") != "ok":
            continue
        path = root / item["path"]
        frame, _ = read_table(path)
        if frame is None:
            continue
        gsm_by_file[item["path"]] = gsm_ids(frame.columns) | gsm_ids(frame.astype(str).to_numpy().ravel())

    overlap_rows = []
    files = list(gsm_by_file)
    for i, left in enumerate(files):
        for right in files[i + 1:]:
            a, b = gsm_by_file[left], gsm_by_file[right]
            if not a or not b:
                continue
            overlap = sorted(a & b)
            overlap_rows.append({
                "file_a": left, "file_b": right,
                "n_gsm_a": len(a), "n_gsm_b": len(b),
                "n_exact_shared_gsm": len(overlap),
                "shared_gsm_ids": ";".join(overlap[:100]),
            })
    overlap_path = out / "02_exact_sample_id_overlaps.csv"
    pd.DataFrame(overlap_rows).to_csv(overlap_path, index=False)

    target_files = inventory[
        inventory.get("target_gse_in_path", pd.Series(dtype=str)).fillna("").astype(str).str.len() > 0
    ] if not inventory.empty else pd.DataFrame()
    print(f"Root: {root}")
    print(f"Candidate tables inventoried: {len(inventory)}")
    print(f"Files whose paths mention target GEO series: {len(target_files)}")
    print(f"Inventory: {inventory_path}")
    print(f"Exact GSM overlap table: {overlap_path}")
    print()
    if target_files.empty:
        print("BLOCKER: no readable candidate table path explicitly mentions GSE67462 or GSE67520.")
    else:
        print(target_files[["path", "read_status"]].to_string(index=False))
    print()
    print("Interpretation guardrails:")
    print("- A shared time label is not proof that samples are paired.")
    print("- Exact GSM overlap across two GEO series may be zero and is not required if a documented sample map exists.")
    print("- Do not run M0/M1/M2 until the inventory identifies expression and regulatory matrices,")
    print("  a defensible sample/time mapping, and gene identifiers that can be harmonized.")
    print("- All scaling, feature selection, imputation and model tuning must be fit inside each training fold.")


if __name__ == "__main__":
    main()
