"""Run the independent validation milestones added after Stage 2.6.

Usage:
    python validate_pipeline.py

Stage 2.7 reads the Stage 2.6 common human gene matrix and writes its reports
under results/Dynamics/stage2_7. The leakage-free Stage 2.4 implementation is
exposed as a reusable function in dynamics.stage24.
"""
from Dynamics import stage2_6
from dynamics.validation import stage2_7


if __name__ == "__main__":
    print("Running Stage 2.6 to refresh the common gene space...")
    stage2_6()
    print("Running Stage 2.7 independent validation...")
    summary = stage2_7()
    print("\nStage 2.7 summary:")
    print(summary.to_string(index=False))
