"""Run the frozen GSE67462 multimodal validation audit.

The runner deliberately performs a provenance/availability gate before any
statistical analysis. It does not download data and does not alter Z6 model
selection or support criteria.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

REQUIRED_MODALITIES = {
    "total_oct4": "GSE67520 total Oct4 ChIP-seq",
    "h3k4me1": "GSE67520 H3K4me1",
    "h3k27ac": "GSE67520 H3K27ac",
    "h3k4me3": "GSE67520 H3K4me3",
    "h3k27me3": "GSE67520 H3K27me3",
    "rnapii": "GSE67520 RNAPII",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--expression-results",
        default="results/Dynamics/z6_gse67462_signal_attribution",
    )
    parser.add_argument(
        "--modality-root",
        default="Data/GSE67520",
        help="Directory containing processed GSE67520 modality files.",
    )
    parser.add_argument(
        "--output",
        default="results/Dynamics/z6_gse67462_multimodal_validation",
    )
    args = parser.parse_args()

    expression_root = Path(args.expression_results)
    modality_root = Path(args.modality_root)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    availability = {}
    for key, label in REQUIRED_MODALITIES.items():
        candidates = sorted(modality_root.glob(f"**/*{key}*"))
        availability[key] = {
            "label": label,
            "available": bool(candidates),
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
            "endpoint_only_control": True,
            "random_size_matched_gene_set_control": True,
        },
    }

    (output / "00_multimodal_availability.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("GSE67462 multimodal validation gate complete.")
    print(f" expression_results_available={expression_ok}")
    for key, item in availability.items():
        print(f" {key}: {'AVAILABLE' if item['available'] else 'MISSING'}")
    print(f" status={manifest['status']}")
    print(f" Outputs: {output}")


if __name__ == "__main__":
    main()
