from pathlib import Path
import gzip
import tarfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "Data"
OUT = ROOT / "results" / "GSE67520"
OUT.mkdir(parents=True, exist_ok=True)
ARCHIVE = DATA / "GSE67520_RAW.tar"

MARK_ALIASES = {
    "flag": "FLAG_control",
    "oct4": "Oct4",
    "k4me1": "H3K4me1",
    "h3k4me1": "H3K4me1",
    "k27ac": "H3K27ac",
    "h3k27ac": "H3K27ac",
    "k4me3": "H3K4me3",
    "h3k4me3": "H3K4me3",
    "k27me3": "H3K27me3",
    "h3k27me3": "H3K27me3",
    "rnapii": "RNAPII",
    "input_dna": "input_DNA",
}


def parse_peak_file(path):
    rows = []
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) < 3:
                continue
            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            if end <= start:
                continue
            rows.append((chrom, start, end))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def classify(name):
    n = name.lower()
    mark = "unknown"
    for token, label in MARK_ALIASES.items():
        if token in n:
            mark = label
            break
    day = None
    if "ipsc" in n:
        day = "iPSC"
    else:
        import re
        m = re.search(r"(?:^|[_-])d(?:ay)?[_-]?(\d+)", n)
        if m:
            day = f"day{m.group(1)}"
        else:
            m = re.search(r"day(\d+)", n)
            if m:
                day = f"day{m.group(1)}"
    return mark, day or "unknown"


def inspect_archive():
    if not ARCHIVE.exists():
        print("Brak GSE67520_RAW.tar")
        return
    extract = ROOT / "GSE67520_extracted"
    extract.mkdir(exist_ok=True)
    marker = extract / ".archive_extracted"
    if not marker.exists():
        with tarfile.open(ARCHIVE, "r") as tar:
            tar.extractall(extract)
        marker.write_text("ok", encoding="utf-8")

    files = [p for p in extract.rglob("*") if p.is_file() and p.name != marker.name]
    print("=== GSE67520 RAW archive ===")
    print(f"Files: {len(files)}")
    inventory = []
    for p in files:
        mark, day = classify(p.name)
        size = p.stat().st_size
        inventory.append({"file": p.name, "size_bytes": size, "mark": mark, "stage": day})
        print(f"{p.name} | {size:,} bytes | mark={mark} | stage={day}")
    inv = pd.DataFrame(inventory)
    inv.to_csv(OUT / "01_RAW_file_inventory.csv", index=False)

    summaries = []
    for _, row in inv.iterrows():
        path = next((p for p in files if p.name == row["file"]), None)
        if path is None:
            continue
        try:
            peaks = parse_peak_file(path)
            if peaks.empty:
                continue
            lengths = peaks["end"] - peaks["start"]
            summaries.append({
                "file": row["file"],
                "mark": row["mark"],
                "stage": row["stage"],
                "peaks": len(peaks),
                "median_peak_length": float(lengths.median()),
                "mean_peak_length": float(lengths.mean()),
                "max_peak_length": int(lengths.max()),
                "chromosomes": int(peaks["chrom"].nunique()),
                "total_covered_bp": int(lengths.sum()),
            })
        except Exception as e:
            print(f"Could not parse {path.name}: {e}")

    summary = pd.DataFrame(summaries)
    if summary.empty:
        raise ValueError("Nie znaleziono poprawnych BED/BROADPEAK plików w GSE67520_RAW.tar")
    summary.to_csv(OUT / "02_peak_summary.csv", index=False)

    counts = summary.pivot_table(index="stage", columns="mark", values="peaks", aggfunc="sum", fill_value=0)
    counts.to_csv(OUT / "03_peak_counts_by_stage_and_mark.csv")

    fig, ax = plt.subplots(figsize=(12, 7))
    counts.plot(kind="bar", ax=ax)
    ax.set_ylabel("Number of peaks")
    ax.set_xlabel("Reprogramming stage")
    ax.set_title("GSE67520 - peak counts by stage and chromatin mark")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(OUT / "04_peak_counts.png", dpi=250); plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 7))
    for mark, sub in summary.groupby("mark"):
        if mark in {"input_DNA", "FLAG_control", "unknown"}:
            continue
        ax.plot(sub["stage"], sub["peaks"], marker="o", label=mark)
    ax.set_ylabel("Number of peaks")
    ax.set_xlabel("Reprogramming stage")
    ax.set_title("GSE67520 - peak number across reprogramming")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "05_peak_dynamics.png", dpi=250); plt.close(fig)

    stage_order = [f"day{x}" for x in [0, 1, 3, 5, 7, 11, 15, 18]] + ["iPSC"]
    summary["stage_order"] = summary["stage"].map({s: i for i, s in enumerate(stage_order)}).fillna(999)
    summary = summary.sort_values(["stage_order", "mark"])
    summary.to_csv(OUT / "06_peak_summary_sorted.csv", index=False)

    report = [
        "Dataset: GSE67520",
        "Title: Oct4 binding and Histone modification profiling during OSKM-mediated 2nd reprogramming",
        "Organism: Mus musculus",
        "Experiment type: ChIP-seq",
        "RAW archive content: BED and broadPeak peak calls, not a gene-expression matrix.",
        "The dataset profiles Oct4 binding, histone modifications and RNAPII during reprogramming.",
        f"Parsed peak files: {len(summary)}",
        "Marks observed: " + ", ".join(sorted(summary["mark"].unique())),
        "Stages observed: " + ", ".join(stage_order),
        "Important interpretation: a peak is a genomic region with enriched ChIP-seq signal; more peaks does not automatically mean higher gene expression.",
        "FLAG files are kept as a separate control category rather than being interpreted as a histone mark.",
        "This exploration summarizes peak counts, lengths, genomic coverage and stage/mark dynamics. It does not perform differential peak calling.",
    ]
    (OUT / "REPORT.txt").write_text("\n".join(report), encoding="utf-8")
    print("=== end GSE67520 ===")


inspect_archive()
print("GSE67520 exploration complete")
