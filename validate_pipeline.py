"""Run the independent validation milestones added after Stage 2.6."""
from pathlib import Path
import re
import pandas as pd

from dynamics.stage26_v2 import stage2_6_robust
from dynamics.validation import stage2_7
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
    from Dynamics import time_hours
    rows=[]
    for ds,g in metadata.groupby("dataset",sort=True):
        timed=pd.Series([time_hours(ds,s) for s in g["sample"].astype(str)]).dropna(); unique=sorted(set(float(x) for x in timed))
        role="trajectory" if len(unique)>=2 else "single_timepoint" if len(unique)==1 else "context_only"
        rows.append({"dataset":ds,"n_samples":len(g),"n_timed_samples":len(timed),"n_unique_times":len(unique),"time_hours":",".join(map(str,unique)),"role":role})
    out=pd.DataFrame(rows); out.to_csv(ROOT/"results/Dynamics/stage2_7/00_dataset_roles.csv",index=False)
    print("\nStage 2.7 dataset roles:"); print(out.to_string(index=False)); return out


if __name__ == "__main__":
    print("Running bounded Stage 2.6 refresh (MyGene/NCBI; no giant BioMart requests)...")
    refresh=stage2_6_robust()
    if refresh.get("status")!="sufficient_common_human_gene_space":
        if COMMON_MATRIX.exists() and COMMON_METADATA.exists(): print("WARNING: new Stage 2.6 mapping is insufficient; continuing with the last valid common-space files.")
        else: raise RuntimeError("Stage 2.6 did not produce a sufficient common gene space and no previous valid files exist.")
    _recover_legacy_metadata(); _write_dataset_roles()
    print("Running Stage 2.7 independent validation...")
    summary=stage2_7(); print("\nStage 2.7 summary:"); print(summary.to_string(index=False))
