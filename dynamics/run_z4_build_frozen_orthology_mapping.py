"""Build a conservative, provenance-aware frozen human-mouse mapping for Z4.

The legacy Stage 2.6 table is treated as a candidate cross-species mapping, not
as a source of arbitrary best hits. Only one-to-one human<->mouse relationships
are retained. Ambiguous human genes and ambiguous mouse genes are excluded.

Usage:
    python -m dynamics.run_z4_build_frozen_orthology_mapping

An explicit --input may be supplied when the legacy cache is stored elsewhere.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

DEFAULT_INPUT = Path("results/Dynamics/stage2_6/cache/mouse_refseq_to_human_mygene.tsv")
DEFAULT_OUTPUT = Path("results/Dynamics/z4_frozen_orthology_mapping")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def build_mapping(input_path: Path, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(input_path, sep="\t", dtype=str, keep_default_na=False)
    required = {"mouse_gene", "human_gene"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    pairs = df[["mouse_gene", "human_gene"]].copy()
    for col in pairs.columns:
        pairs[col] = pairs[col].astype(str).str.strip()
    pairs = pairs[(pairs["mouse_gene"] != "") & (pairs["human_gene"] != "")]
    pairs = pairs.drop_duplicates().reset_index(drop=True)

    human_counts = pairs.groupby("human_gene")["mouse_gene"].nunique()
    mouse_counts = pairs.groupby("mouse_gene")["human_gene"].nunique()

    human_ambiguous = set(human_counts[human_counts > 1].index)
    mouse_ambiguous = set(mouse_counts[mouse_counts > 1].index)

    frozen = pairs[
        ~pairs["human_gene"].isin(human_ambiguous)
        & ~pairs["mouse_gene"].isin(mouse_ambiguous)
    ].copy()
    frozen = frozen.sort_values(["human_gene", "mouse_gene"]).reset_index(drop=True)
    frozen.insert(0, "mapping_status", "one_to_one")

    ambiguous_human = pd.DataFrame({"human_gene": sorted(human_ambiguous)})
    ambiguous_human["mapping_status"] = "ambiguous_excluded"
    ambiguous_human.to_csv(output_dir / "02_ambiguous_human_genes.csv", index=False)

    frozen.to_csv(output_dir / "01_frozen_human_mouse_mapping.csv", index=False)

    source_sha = sha256_file(input_path)
    stats = {
        "input": str(input_path).replace("\\", "/"),
        "input_sha256": source_sha,
        "input_rows": int(len(df)),
        "nonempty_pairs": int(len(pairs)),
        "unique_human_candidates": int(pairs["human_gene"].nunique()),
        "unique_mouse_candidates": int(pairs["mouse_gene"].nunique()),
        "ambiguous_human_genes": int(len(human_ambiguous)),
        "ambiguous_mouse_genes": int(len(mouse_ambiguous)),
        "frozen_one_to_one_pairs": int(len(frozen)),
        "frozen_unique_human_genes": int(frozen["human_gene"].nunique()),
        "frozen_unique_mouse_genes": int(frozen["mouse_gene"].nunique()),
        "policy": "retain only one-to-one human<->mouse relationships; exclude ambiguous mappings",
        "mapping_role": "candidate frozen orthology for Z4 external validation; not a causal annotation",
    }
    with (output_dir / "03_mapping_summary.json").open("w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2, sort_keys=True)

    return stats


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input mapping not found: {args.input}")

    stats = build_mapping(args.input, args.output)
    print("Z4 FROZEN ORTHOLOGY MAPPING")
    print(f"input: {stats['input']}")
    print(f"input rows: {stats['input_rows']}")
    print(f"nonempty unique pairs: {stats['nonempty_pairs']}")
    print(f"ambiguous human genes excluded: {stats['ambiguous_human_genes']}")
    print(f"ambiguous mouse genes excluded: {stats['ambiguous_mouse_genes']}")
    print(f"frozen one-to-one pairs: {stats['frozen_one_to_one_pairs']}")
    print(f"input SHA256: {stats['input_sha256']}")
    print(f"output: {args.output.as_posix()}")


if __name__ == "__main__":
    main()
