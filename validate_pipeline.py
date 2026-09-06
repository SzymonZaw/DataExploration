"""Run the independent validation milestones added after Stage 2.6.

Usage:
    python validate_pipeline.py

Stage 2.7 reads the Stage 2.6 common human gene matrix and writes its reports
under results/Dynamics/stage2_7. The leakage-free Stage 2.4 implementation is
exposed as a reusable function in dynamics.stage24.
"""
from pathlib import Path

from Dynamics import stage2_6
from dynamics.validation import stage2_7

ROOT = Path(__file__).resolve().parent
COMMON_MATRIX = ROOT / "results" / "Dynamics" / "stage2_6" / "06_common_human_gene_matrix.csv"
COMMON_METADATA = ROOT / "results" / "Dynamics" / "stage2_6" / "07_common_gene_sample_metadata.csv"


if __name__ == "__main__":
    print("Running Stage 2.6 to refresh the common gene space...")
    try:
        stage2_6()
    except Exception as exc:
        if COMMON_MATRIX.exists() and COMMON_METADATA.exists():
            print(f"WARNING: Stage 2.6 refresh failed: {type(exc).__name__}: {exc}")
            print("WARNING: continuing Stage 2.7 with the last valid Stage 2.6 common-space files.")
        else:
            raise

    print("Running Stage 2.7 independent validation...")
    summary = stage2_7()
    print("\nStage 2.7 summary:")
    print(summary.to_string(index=False))
