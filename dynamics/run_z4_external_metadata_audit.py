"""Pre-transfer GEO metadata audit for Z4 external specificity.

Downloads and parses GEO Series SOFT metadata only. It does not read or score
expression data and does not evaluate the frozen Z4/Z6 representation.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/soft/{accession}_family.soft.gz"
EXPECTED = {
    "GSE297234": {"min_timepoints": 3, "keywords": ["oskm", "sendai"], "role": "primary_external_target"},
    "GSE28688": {"min_timepoints": 3, "keywords": ["oskm"], "role": "specificity_context_challenge"},
}


def _url(accession: str) -> str:
    n = int(re.search(r"\d+", accession).group())
    return GEO_URL.format(prefix=f"GSE{n // 1000}nnn", accession=accession)


def _download(accession: str, cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / f"{accession}_family.soft.gz"
    if not out.exists():
        print(f"Downloading {accession} SOFT metadata...")
        urllib.request.urlretrieve(_url(accession), out)
    return out


def _parse_soft(path: Path) -> pd.DataFrame:
    rows, current = [], None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    rows.append(current)
                current = {"sample_id": line.split("=", 1)[1].strip()}
                continue
            if current is None:
                continue
            if line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_platform_id = "):
                current["platform"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_source_name_ch1 = "):
                current["source"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_characteristics_ch1 = "):
                current.setdefault("characteristics", []).append(line.split("=", 1)[1].strip())
        if current:
            rows.append(current)
    df = pd.DataFrame(rows)
    if "characteristics" in df:
        df["characteristics"] = df["characteristics"].apply(lambda x: " | ".join(x) if isinstance(x, list) else "")
    return df


def _infer_time(text: str):
    for pattern, conv in [
        (r"day\s*([0-9]+)", lambda x: float(x)),
        (r"\bd\s*([0-9]+)\b", lambda x: float(x)),
        (r"([0-9]+)\s*h(?:ours?)?", lambda x: float(x) / 24.0),
    ]:
        match = re.search(pattern, text.lower())
        if match:
            return conv(match.group(1))
    return None


def _audit(accession: str, df: pd.DataFrame) -> dict:
    cfg = EXPECTED[accession]
    text = df.fillna("").astype(str).agg(" | ".join, axis=1)
    inferred = text.map(_infer_time)
    n_times = int(inferred.dropna().nunique())
    intervention_hits = sum(any(k in t.lower() for k in cfg["keywords"]) for t in text)
    nuisance_cols = {c: c in df.columns for c in ["sample_id", "title", "source", "platform", "characteristics"]}
    return {
        "candidate": accession,
        "role": cfg["role"],
        "n_samples": int(len(df)),
        "n_inferred_timepoints": n_times,
        "temporal_gate": n_times >= cfg["min_timepoints"],
        "intervention_metadata_hits": int(intervention_hits),
        "intervention_gate": intervention_hits > 0,
        "nuisance_metadata_fields": nuisance_cols,
        "transfer_score_evaluated": False,
        "next_gate": "prepare_frozen_human_mouse_mapping_and_endpoint_audit",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidate", required=True, choices=sorted(EXPECTED))
    ap.add_argument("--cache", default="Data/GEO_Z4")
    ap.add_argument("--output", default=None)
    args = ap.parse_args()

    soft = _download(args.candidate, Path(args.cache))
    df = _parse_soft(soft)
    audit = _audit(args.candidate, df)
    text = df.fillna("").astype(str).agg(" | ".join, axis=1)
    df["inferred_time_days"] = text.map(_infer_time)

    out = Path(args.output) if args.output else Path("results/Dynamics/z4_external_metadata_audit") / args.candidate
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "metadata.csv", index=False)
    with gzip.open(soft, "rt", encoding="utf-8", errors="replace") as fh:
        (out / "metadata_raw.txt").write_text(fh.read(), encoding="utf-8")
    (out / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps(audit, indent=2))
    print(f"metadata: {out / 'metadata.csv'}")
    print(f"audit: {out / 'audit.json'}")


if __name__ == "__main__":
    main()
