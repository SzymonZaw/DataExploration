"""Run independent validation milestones.

Default mode is fast and reuses the existing Stage 2.6 common gene space.
Use --refresh to rebuild Stage 2.6, --stage28 for Stage 2.8, --stage29
for Stage 2.9, --stage291 for robust invariant temporal programs, and
--stage292 for latent dataset-invariant reprogramming progress validation.
"""
from pathlib import Path
import argparse
import re
import pandas as pd

from dynamics.validation import stage2_7
from dynamics.stage28 import stage2_8
from Dynamics import PCA_FILES, GSE28688_ROW_SAMPLE

ROOT = Path(__file__).resolve().parent
STAGE26 = ROOT / "results" / "Dynamics" / "stage2_6"
COMMON_MATRIX = STAGE26 / "06_common_human_gene_matrix.csv"
COMMON_METADATA = STAGE26 / "07_common_gene_sample_metadata.csv"


def _recover_legacy_metadata():
    if not COMMON_MATRIX.exists() or not COMMON_METADATA.exists(): return 0
    matrix = pd.read_csv(COMMON_MATRIX, index_col=0); metadata = pd.read_csv(COMMON_METADATA)
    if not {"dataset", "sample"}.issubset(metadata.columns): return 0
    changed = 0
    for ds in metadata["dataset"].astype(str).unique():
        actual = [str(c) for c in matrix.columns if str(c).startswith(ds + "__")]
        if not actual or not all(re.fullmatch(rf"{re.escape(ds)}__\d+", c) for c in actual): continue
        samples=[]; pca=PCA_FILES.get(ds)
        if pca is not None and pca.exists():
            try: samples=[str(v) for v in pd.read_csv(pca,index_col=0).index]
            except Exception: samples=[]
        if ds=="GSE28688" and len(actual)==len(GSE28688_ROW_SAMPLE): samples=list(GSE28688_ROW_SAMPLE)
        if len(samples)!=len(actual): continue
        mask=metadata["dataset"].astype(str).eq(ds)
        for idx in metadata.index[mask]:
            sample=str(metadata.at[idx,"sample"]); m=re.search(r"(?:__)?(\d+)$",sample)
            if not m and "matrix_column" in metadata.columns: m=re.search(r"__(\d+)$",str(metadata.at[idx,"matrix_column"]))
            if m and int(m.group(1))<len(samples):
                n=m.group(1); metadata.at[idx,"sample"]=samples[int(n)]; metadata.at[idx,"matrix_column"]=f"{ds}__{n}"; changed+=1
    if changed:
        metadata.to_csv(COMMON_METADATA,index=False); print(f"Stage 2.7: recovered {changed} legacy sample-to-column mappings from PCA/GEO order.")
    return changed


def _write_dataset_roles():
    if not COMMON_METADATA.exists(): return pd.DataFrame()
    metadata=pd.read_csv(COMMON_METADATA)
    from dynamics.validation import _time_hours_for_validation, _strip_dataset_prefix
    rows=[]
    for ds,g in metadata.groupby("dataset",sort=True):
        times=[]
        for row_index, sample in enumerate(g["sample"].astype(str)):
            idx=row_index if ds=="GSE28688" else None
            t=_time_hours_for_validation(ds,_strip_dataset_prefix(sample),idx)
            if pd.notna(t): times.append(float(t))
        unique=sorted(set(times))
        role="trajectory" if len(unique)>=2 else "single_timepoint" if len(unique)==1 else "context_only"
        rows.append({"dataset":ds,"n_samples":len(g),"n_timed_samples":len(times),"n_unique_times":len(unique),"time_hours":",".join(map(str,unique)),"role":role})
    out=pd.DataFrame(rows); out.to_csv(ROOT/"results/Dynamics/stage2_7/00_dataset_roles.csv",index=False)
    print("\nStage 2.7 dataset roles:"); print(out.to_string(index=False)); return out


def _parse_args():
    parser=argparse.ArgumentParser(description="Run DataExploration validation stages.")
    parser.add_argument("--refresh", action="store_true", help="Rebuild the Stage 2.6 common human-gene space before validation.")
    parser.add_argument("--stage28", action="store_true", help="Run only Stage 2.8 diagnostics using the existing common space.")
    parser.add_argument("--stage29", action="store_true", help="Run only Stage 2.9 leakage-free invariant module-state validation.")
    parser.add_argument("--stage291", action="store_true", help="Run only Stage 2.9.1 robust invariant temporal-program validation with permutation null.")
    parser.add_argument("--stage292", action="store_true", help="Run only Stage 2.9.2 latent dataset-invariant reprogramming progress validation with permutation null.")
    return parser.parse_args()


def _require_common_space():
    if not COMMON_MATRIX.exists() or not COMMON_METADATA.exists():
        raise RuntimeError("Stage 2.6 common-space files are missing. Run: python validate_pipeline.py --refresh")


def main():
    args = _parse_args()

    if args.refresh:
        from dynamics.stage26_v2 import stage2_6_robust
        print("Running Stage 2.6 refresh (MyGene/NCBI; no giant BioMart requests)...")
        refresh = stage2_6_robust()
        if refresh.get("status") != "sufficient_common_human_gene_space":
            if not (COMMON_MATRIX.exists() and COMMON_METADATA.exists()):
                raise RuntimeError("Stage 2.6 did not produce a sufficient common gene space and no previous valid files exist.")
            print("WARNING: new Stage 2.6 mapping is insufficient; keeping the last valid common-space files.")

    _require_common_space()

    if args.stage28:
        print("Using existing Stage 2.6 common space. Skipping Stage 2.6 and Stage 2.7.")
        print("\nRunning Stage 2.8 cross-dataset diagnostics...")
        result = stage2_8()
        print(f"Stage 2.8 result: {result}")
        return

    if args.stage29:
        print("Using existing Stage 2.6 common space. Skipping Stage 2.6, Stage 2.7, and Stage 2.8.")
        print("\nRunning Stage 2.9 leakage-free invariant module-state validation...")
        from dynamics.stage29 import stage2_9
        result = stage2_9()
        print(f"Stage 2.9 result: {result}")
        return

    if args.stage291:
        print("Using existing Stage 2.6 common space. Skipping Stage 2.6, Stage 2.7, Stage 2.8, and Stage 2.9.")
        print("\nRunning Stage 2.9.1 robust invariant temporal-program validation...")
        from dynamics.stage291 import stage2_9_1
        result = stage2_9_1()
        print(f"Stage 2.9.1 result: {result}")
        return

    if args.stage292:
        print("Using existing Stage 2.6 common space. Skipping Stage 2.6, Stage 2.7, Stage 2.8, Stage 2.9, and Stage 2.9.1.")
        print("\nRunning Stage 2.9.2 latent dataset-invariant reprogramming progress validation...")
        from dynamics.stage292 import stage2_9_2
        result = stage2_9_2()
        print(f"Stage 2.9.2 result: {result}")
        return

    _recover_legacy_metadata(); _write_dataset_roles()
    print("Running Stage 2.7 independent validation using existing Stage 2.6 common space...")
    summary=stage2_7(); print("\nStage 2.7 summary:"); print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
