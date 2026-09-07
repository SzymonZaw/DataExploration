"""Run independent validation milestones."""
from pathlib import Path
import argparse
import re
import pandas as pd
from dynamics.validation import stage2_7
from dynamics.stage28 import stage2_8
from Dynamics import PCA_FILES, GSE28688_ROW_SAMPLE
ROOT=Path(__file__).resolve().parent
STAGE26=ROOT/"results"/"Dynamics"/"stage2_6"
COMMON_MATRIX=STAGE26/"06_common_human_gene_matrix.csv"
COMMON_METADATA=STAGE26/"07_common_gene_sample_metadata.csv"

def _recover_legacy_metadata():
    if not COMMON_MATRIX.exists() or not COMMON_METADATA.exists(): return 0
    matrix=pd.read_csv(COMMON_MATRIX,index_col=0); metadata=pd.read_csv(COMMON_METADATA)
    if not {"dataset","sample"}.issubset(metadata.columns): return 0
    changed=0
    for ds in metadata["dataset"].astype(str).unique():
        actual=[str(c) for c in matrix.columns if str(c).startswith(ds+"__")]
        if not actual or not all(re.fullmatch(rf"{re.escape(ds)}__\d+",c) for c in actual): continue
        samples=[]; pca=PCA_FILES.get(ds)
        if pca is not None and pca.exists():
            try:samples=[str(v) for v in pd.read_csv(pca,index_col=0).index]
            except Exception:samples=[]
        if ds=="GSE28688" and len(actual)==len(GSE28688_ROW_SAMPLE): samples=list(GSE28688_ROW_SAMPLE)
        if len(samples)!=len(actual):continue
        mask=metadata["dataset"].astype(str).eq(ds)
        for idx in metadata.index[mask]:
            sample=str(metadata.at[idx,"sample"]);m=re.search(r"(?:__)?(\d+)$",sample)
            if not m and "matrix_column" in metadata.columns:m=re.search(r"__(\d+)$",str(metadata.at[idx,"matrix_column"]))
            if m and int(m.group(1))<len(samples):
                n=m.group(1);metadata.at[idx,"sample"]=samples[int(n)];metadata.at[idx,"matrix_column"]=f"{ds}__{n}";changed+=1
    if changed:metadata.to_csv(COMMON_METADATA,index=False);print(f"Stage 2.7: recovered {changed} legacy sample-to-column mappings from PCA/GEO order.")
    return changed

def _write_dataset_roles():
    if not COMMON_METADATA.exists():return pd.DataFrame()
    metadata=pd.read_csv(COMMON_METADATA);from dynamics.validation import _time_hours_for_validation,_strip_dataset_prefix
    rows=[]
    for ds,g in metadata.groupby("dataset",sort=True):
        times=[]
        for row_index,sample in enumerate(g["sample"].astype(str)):
            idx=row_index if ds=="GSE28688" else None;t=_time_hours_for_validation(ds,_strip_dataset_prefix(sample),idx)
            if pd.notna(t):times.append(float(t))
        unique=sorted(set(times));role="trajectory" if len(unique)>=2 else "single_timepoint" if len(unique)==1 else "context_only"
        rows.append({"dataset":ds,"n_samples":len(g),"n_timed_samples":len(times),"n_unique_times":len(unique),"time_hours":",".join(map(str,unique)),"role":role})
    out=pd.DataFrame(rows);out.to_csv(ROOT/"results/Dynamics/stage2_7/00_dataset_roles.csv",index=False);print("\nStage 2.7 dataset roles:");print(out.to_string(index=False));return out

def _parse_args():
    p=argparse.ArgumentParser(description="Run DataExploration validation and model experiments.")
    for name in ("refresh","stage28","stage29","stage291","stage292","stage293","stage294","stage295","stage296","stage297","stage298","stage299","stage2910","stage2911","stage2912","stage2913","stage2914","stage2915","stage2916","stage2917","stage2918","stage2919","stage2920","stage2921","stage2922","stage2923","stage2924","stage2925","stage2926","stage2927","stage2928","stage2929","stage2930","stage2931","stage2932","stage210","dynamic-state-model","dynamic-state-model-v2","dynamic-state-memory"):p.add_argument(f"--{name}",action="store_true")
    return p.parse_args()

def _require_common_space():
    if not COMMON_MATRIX.exists() or not COMMON_METADATA.exists():raise RuntimeError("Stage 2.6 common-space files are missing. Run: python validate_pipeline.py --refresh")

def main():
    args=_parse_args()
    if args.refresh:
        from dynamics.stage26_v2 import stage2_6_robust
        print("Running Stage 2.6 refresh (MyGene/NCBI; no giant BioMart requests)...");refresh=stage2_6_robust()
        if refresh.get("status")!="sufficient_common_human_gene_space":
            if not (COMMON_MATRIX.exists() and COMMON_METADATA.exists()):raise RuntimeError("Stage 2.6 did not produce a sufficient common gene space and no previous valid files exist.")
            print("WARNING: new Stage 2.6 mapping is insufficient; keeping the last valid common-space files.")
    _require_common_space()
    if getattr(args,"dynamic_state_model"):
        print("Running central DynamicStateModel: leakage-free leave-one-dataset-out prediction...");from dynamics.run_dynamic_state_lodo import run;print(run());return
    if getattr(args,"dynamic_state_model_v2"):
        print("Running DynamicStateModel v2: delta-t conditioned leakage-free LODO prefix-to-future forecasting...");from dynamics.run_dynamic_state_v2 import run;print(run());return
    if getattr(args,"dynamic_state_memory"):
        print("Running memory-aware DynamicStateModel: Markov vs history-aware LODO prefix-to-future forecasting...");from dynamics.run_dynamic_state_memory import run;print(run());return
    if args.stage210:
        print("Running Stage 2.10 leakage-free conserved biological transition module validation; no ODE/state-space model...");from dynamics.stage210 import run;print(f"Stage 2.10 result: {run()}");return
    if args.stage28:
        print("Running Stage 2.8 cross-dataset diagnostics...");print(f"Stage 2.8 result: {stage2_8()}");return
    if args.stage29:
        from dynamics.stage29 import stage2_9;print(f"Stage 2.9 result: {stage2_9()}");return
    if args.stage291:
        from dynamics.stage291 import stage2_9_1;print(f"Stage 2.9.1 result: {stage2_9_1()}");return
    if args.stage292:
        from dynamics.stage292 import run;print(f"Stage 2.9.2 result: {run()}");return
    if args.stage293:
        from dynamics.stage293_progress import run;print(f"Stage 2.9.3 result: {run()}");return
    if args.stage294:
        from dynamics.stage294 import run;print(f"Stage 2.9.4 result: {run()}");return
    if args.stage295:
        from dynamics.stage295 import run;print(f"Stage 2.9.5 result: {run()}");return
    if args.stage296:
        from dynamics.stage296 import run;print(f"Stage 2.9.6 result: {run()}");return
    if args.stage297:
        from dynamics.stage297 import run;print(f"Stage 2.9.7 result: {run()}");return
    if args.stage298:
        from dynamics.stage298 import run;print(f"Stage 2.9.8 result: {run()}");return
    if args.stage299:
        from dynamics.stage299 import run;print(f"Stage 2.9.9 result: {run()}");return
    if args.stage2910:
        from dynamics.stage2910 import run;print(f"Stage 2.9.10 result: {run()}");return
    if args.stage2911:
        from dynamics.stage2911 import run;print(f"Stage 2.9.11 result: {run()}");return
    if args.stage2912:
        from dynamics.stage2912 import run;print(f"Stage 2.9.12 result: {run()}");return
    if args.stage2913:
        from dynamics.stage2913 import run;print(f"Stage 2.9.13 result: {run()}");return
    if args.stage2914:
        from dynamics.stage2914 import run;print(f"Stage 2.9.14 result: {run()}");return
    if args.stage2915:
        from dynamics.stage2915 import run;print(f"Stage 2.9.15 result: {run()}");return
    if args.stage2916:
        from dynamics.stage2916 import run;print(f"Stage 2.9.16 result: {run()}");return
    if args.stage2917:
        from dynamics.stage2917 import run;print(f"Stage 2.9.17 result: {run()}");return
    if args.stage2918:
        from dynamics.stage2918 import run;print(f"Stage 2.9.18 result: {run()}");return
    if args.stage2919:
        from dynamics.stage2919 import run;print(f"Stage 2.9.19 result: {run()}");return
    if args.stage2920:
        from dynamics.stage2920 import run;print(f"Stage 2.9.20 result: {run()}");return
    if args.stage2921:
        from dynamics.stage2921 import run;print(f"Stage 2.9.21 result: {run()}");return
    if args.stage2922:
        from dynamics.stage2922 import run;print(f"Stage 2.9.22 result: {run()}");return
    if args.stage2923:
        from dynamics.stage2923 import run;print(f"Stage 2.9.23 result: {run()}");return
    if args.stage2924:
        from dynamics.stage2924 import run;print(f"Stage 2.9.24 result: {run()}");return
    if args.stage2925:
        from dynamics.stage2925 import run;print(f"Stage 2.9.25 result: {run()}");return
    if args.stage2926:
        from dynamics.stage2926 import run;print(f"Stage 2.9.26 result: {run()}");return
    if args.stage2927:
        from dynamics.stage2927 import run;print(f"Stage 2.9.27 result: {run()}");return
    if args.stage2928:
        from dynamics.stage2928 import run;print(f"Stage 2.9.28 result: {run()}");return
    if args.stage2929:
        from dynamics.stage2929 import run;print(f"Stage 2.9.29 result: {run()}");return
    if args.stage2930:
        from dynamics.stage2930 import run;print(f"Stage 2.9.30 result: {run()}");return
    if args.stage2931:
        from dynamics.stage2931 import run;print(f"Stage 2.9.31 result: {run()}");return
    if args.stage2932:
        from dynamics.stage2932 import run;print(f"Stage 2.9.32 result: {run()}");return
    _recover_legacy_metadata();_write_dataset_roles();summary=stage2_7();print(summary.to_string(index=False))
if __name__=="__main__":main()
