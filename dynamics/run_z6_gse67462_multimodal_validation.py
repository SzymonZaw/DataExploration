"""Run the GSE67462/GSE67520 multimodal data-availability gate.

This gate recognizes the actual GSE67520 processed peak-file naming scheme.
It does not perform statistics and does not alter frozen Z6 support.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REQUIRED_MODALITIES = {
    "total_oct4": {"label": "GSE67520 total Oct4 ChIP-seq", "tokens": ("oct4_",)},
    "h3k4me1": {"label": "GSE67520 H3K4me1", "tokens": ("k4me1_",)},
    "h3k27ac": {"label": "GSE67520 H3K27ac", "tokens": ("k27ac_",)},
    "h3k4me3": {"label": "GSE67520 H3K4me3", "tokens": ("k4me3_",)},
    "h3k27me3": {"label": "GSE67520 H3K27me3", "tokens": ("k27me3_",)},
    "rnapii": {"label": "GSE67520 RNAPII", "tokens": ("rnapii_",)},
}


def _find_candidates(root: Path, tokens):
    return sorted(p for p in root.glob("*.gz") if any(t in p.name.lower() for t in tokens))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression-results", default="results/Dynamics/z6_gse67462_signal_attribution")
    parser.add_argument("--modality-root", default="Data/GSE67520")
    parser.add_argument("--output", default="results/Dynamics/z6_gse67462_multimodal_validation")
    args = parser.parse_args()

    expression_root = Path(args.expression_results)
    root = Path(args.modality_root)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    availability = {}
    for key, spec in REQUIRED_MODALITIES.items():
        candidates = _find_candidates(root, spec["tokens"])
        availability[key] = {
            "label": spec["label"],
            "available": bool(candidates),
            "n_files": len(candidates),
            "candidates": [str(p) for p in candidates],
        }

    expression_ok = expression_root.exists() and any(expression_root.iterdir())
    ready = expression_ok and all(v["available"] for v in availability.values())
    manifest = {
        "dataset": "GSE67462",
        "paired_dataset": "GSE67520",
        "time_grid": [0, 1, 3, 5, 7, 11, 15, 18, "iPSC"],
        "expression_results_available": expression_ok,
        "modalities": availability,
        "status": "READY_FOR_ANALYSIS" if ready else "NOT_RUN_DATA_UNAVAILABLE",
        "frozen_rules": {
            "modify_z6_support": False,
            "model_tuning": False,
            "time_label_permutation": True,
        },
    }
    (output / "00_multimodal_availability.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")

    print("GSE67462 multimodal validation gate complete.")
    print(f" expression_results_available={expression_ok}")
    for key, item in availability.items():
        print(f" {key}: {'AVAILABLE' if item['available'] else 'MISSING'} ({item['n_files']} files)")
    print(f" status={manifest['status']}")
    print(f" Outputs: {output}")


if __name__ == "__main__":
    main()
