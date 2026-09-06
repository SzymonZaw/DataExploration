"""Run the independent validation milestones added after Stage 2.6."""
from pathlib import Path
import re
import pandas as pd

from Dynamics import stage2_6, PCA_FILES, GSE28688_ROW_SAMPLE
from dynamics.validation import stage2_7

ROOT = Path(__file__).resolve().parent
STAGE26 = ROOT / "results" / "Dynamics" / "stage2_6"
COMMON_MATRIX = STAGE26 / "06_common_human_gene_matrix.csv"
COMMON_METADATA = STAGE26 / "07_common_gene_sample_metadata.csv"


def _recover_legacy_metadata():
    """Repair old Stage 2.6 metadata whose matrix columns are numeric labels."""
    if not COMMON_MATRIX.exists() or not COMMON_METADATA.exists():
        return 0
    matrix = pd.read_csv(COMMON_MATRIX, index_col=0)
    metadata = pd.read_csv(COMMON_METADATA)
    if not {"dataset", "sample"}.issubset(metadata.columns):
        return 0

    changed = 0
    for ds in metadata["dataset"].astype(str).unique():
        actual = [str(c) for c in matrix.columns if str(c).startswith(ds + "__")]
        if not actual or not all(re.fullmatch(rf"{re.escape(ds)}__\d+", c) for c in actual):
            continue

        samples = []
        pca = PCA_FILES.get(ds)
        if pca is not None and pca.exists():
            try:
                samples = [str(v) for v in pd.read_csv(pca, index_col=0).index]
            except Exception:
                samples = []
        if ds == "GSE28688" and len(actual) == len(GSE28688_ROW_SAMPLE):
            samples = list(GSE28688_ROW_SAMPLE)
        if len(samples) != len(actual):
            continue

        numeric_to_sample = {str(i): samples[i] for i in range(len(samples))}
        mask = metadata["dataset"].astype(str).eq(ds)
        for idx in metadata.index[mask]:
            sample = str(metadata.at[idx, "sample"])
            m = re.search(r"(?:__)?(\d+)$", sample)
            if not m and "matrix_column" in metadata.columns:
                m = re.search(r"__(\d+)$", str(metadata.at[idx, "matrix_column"]))
            if m and m.group(1) in numeric_to_sample:
                n = m.group(1)
                metadata.at[idx, "sample"] = numeric_to_sample[n]
                metadata.at[idx, "matrix_column"] = f"{ds}__{n}"
                changed += 1

    if changed:
        metadata.to_csv(COMMON_METADATA, index=False)
        print(f"Stage 2.7: recovered {changed} legacy sample-to-column mappings from PCA/GEO order.")
    return changed


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

    _recover_legacy_metadata()
    print("Running Stage 2.7 independent validation...")
    summary = stage2_7()
    print("\nStage 2.7 summary:")
    print(summary.to_string(index=False))
