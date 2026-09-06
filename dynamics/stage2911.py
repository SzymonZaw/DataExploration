"""Stage 2.9.11: diagnose instability of leakage-free biological programs.

This stage does not relax Stage 2.9.10 thresholds and does not fit an ODE.
It separates four possible failure modes: stable-gene selection, identifier
coverage, enrichment significance, and program-selection redundancy/filters.
It also measures recurrence of biological terms across training folds.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "Dynamics" / "stage2_9_11"
OUT.mkdir(parents=True, exist_ok=True)

def log(x):
    print(f"Stage 2.9.11: {x}", flush=True)

def _norm(x):
    return str(x).strip().upper()

def _load_result(path):
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None

def _term_rows(data, query):
    rows = []
    if not data:
        return rows
    for r in data.get("result", []):
        try:
            p = float(r.get("p_value", np.nan))
        except Exception:
            p = np.nan
        n = int(r.get("intersection_size", 0) or 0)
        name = str(r.get("name", "") or "")
        native = str(r.get("native", "") or "")
        intersections = r.get("intersections", []) or []
        genes = [str(query[i]) for i, hit in enumerate(intersections)
                 if hit and i < len(query)]
        rows.append({
            "term_id": native,
            "term_name": name,
            "p_value": p,
            "intersection_size": n,
            "intersection_genes": ",".join(sorted(set(genes))),
            "source": str(r.get("source", "")),
        })
    return rows

def _program_filter_reason(name, p, n):
    low = str(name).lower()
    if not np.isfinite(p) or p >= .05:
        return "fdr_or_p_threshold"
    if n < 5:
        return "too_small_intersection"
    blocked = ("hiv", "viral messenger", "metabolism of rna",
               "gene expression", "protein metabolic process",
               "rna polymerase ii transcription",
               "processing of capped intron-containing pre-mrna")
    if any(k in low for k in blocked):
        return "broad_or_artifactual_term_filter"
    return "candidate"

def run():
    log("starting diagnostic; no ODE/state-space model")
    from dynamics.stage2910 import load_space, load_fold_genes, CACHE
    matrix, meta = load_space()
    background = matrix.index.astype(str).tolist()
    fold_genes = load_fold_genes()
    log(f"loaded {len(background):,} background genes and {len(fold_genes)} training folds")

    fold_diag, term_rows, program_diag = [], [], []
    fold_terms = {}
    for i, held in enumerate(sorted(fold_genes), 1):
        query = [str(x) for x in fold_genes[held]]
        cache = CACHE / f"gprofiler_{held}.json"
        data = _load_result(cache)
        if data is None:
            log(f"FOLD {i}: {held}: cached enrichment missing; rerun Stage 2.9.10 first")
            fold_diag.append({"held_out_dataset": held, "n_training_genes": len(query), "cache_present": False,
                              "n_result_terms": 0, "n_fdr05_terms": 0, "n_go_bp_fdr05": 0,
                              "n_reactome_fdr05": 0, "n_kegg_fdr05": 0, "n_candidate_program_terms": 0,
                              "n_selected_programs": 0})
            continue
        rows = _term_rows(data, query)
        frame = pd.DataFrame(rows)
        if len(frame):
            frame["candidate_reason"] = [_program_filter_reason(n, p, s)
                                          for n, p, s in zip(frame.term_name, frame.p_value, frame.intersection_size)]
            frame["held_out_dataset"] = held
            frame["fdr05"] = frame.p_value < .05
            frame["program_candidate"] = frame.candidate_reason.eq("candidate")
            term_rows.extend(frame.to_dict("records"))
            fold_terms[held] = frame[frame.program_candidate].copy()
            selected = frame[frame.program_candidate].sort_values(["p_value", "intersection_size"], ascending=[True, False])
            # Mirror the redundancy rule used by Stage 2.9.10 to show whether
            # program selection, rather than enrichment, is the bottleneck.
            kept = []
            for _, r in selected.iterrows():
                gs = set(filter(None, str(r.intersection_genes).split(",")))
                if any(len(gs & set(filter(None, str(q.intersection_genes).split(",")))) /
                       max(1, len(gs | set(filter(None, str(q.intersection_genes).split(","))))) >= .75 for q in kept):
                    continue
                kept.append(r)
                if len(kept) >= 8:
                    break
            n_go = int(((frame.fdr05) & frame.source.eq("GO:BP")).sum())
            n_reac = int(((frame.fdr05) & frame.source.eq("REAC")).sum())
            n_kegg = int(((frame.fdr05) & frame.source.eq("KEGG")).sum())
            fold_diag.append({"held_out_dataset": held, "n_training_genes": len(query), "cache_present": True,
                              "n_result_terms": len(frame), "n_fdr05_terms": int(frame.fdr05.sum()),
                              "n_go_bp_fdr05": n_go, "n_reactome_fdr05": n_reac, "n_kegg_fdr05": n_kegg,
                              "n_candidate_program_terms": int(frame.program_candidate.sum()),
                              "n_selected_programs": len(kept)})
            for _, r in frame[frame.program_candidate].sort_values("p_value").head(20).iterrows():
                program_diag.append({"held_out_dataset": held, "term_id": r.term_id, "term_name": r.term_name,
                                     "source": r.source, "p_value": r.p_value,
                                     "intersection_size": r.intersection_size,
                                     "candidate_reason": r.candidate_reason})
        else:
            fold_diag.append({"held_out_dataset": held, "n_training_genes": len(query), "cache_present": True,
                              "n_result_terms": 0, "n_fdr05_terms": 0, "n_go_bp_fdr05": 0,
                              "n_reactome_fdr05": 0, "n_kegg_fdr05": 0,
                              "n_candidate_program_terms": 0, "n_selected_programs": 0})

    terms = pd.DataFrame(term_rows)
    if len(terms):
        key = terms[terms["program_candidate"]].copy()
        recurrence = (key.groupby(["source", "term_id", "term_name"], dropna=False)
                      .agg(n_folds=("held_out_dataset", "nunique"),
                           datasets=("held_out_dataset", lambda x: ",".join(sorted(set(x)))))
                      .reset_index().sort_values(["n_folds", "term_name"], ascending=[False, True]))
    else:
        recurrence = pd.DataFrame(columns=["source", "term_id", "term_name", "n_folds", "datasets"])

    fd = pd.DataFrame(fold_diag)
    pd.DataFrame(term_rows).to_csv(OUT / "01_fold_term_diagnostics.csv", index=False)
    fd.to_csv(OUT / "02_fold_enrichment_diagnostics.csv", index=False)
    pd.DataFrame(program_diag).to_csv(OUT / "03_program_selection_diagnostics.csv", index=False)
    recurrence.to_csv(OUT / "04_cross_fold_term_recurrence.csv", index=False)

    if len(fd):
        summary = pd.DataFrame([{
            "n_folds": len(fd),
            "n_folds_with_cache": int(fd.cache_present.sum()),
            "n_folds_with_any_fdr05_term": int((fd.n_fdr05_terms > 0).sum()),
            "n_folds_with_candidate_program": int((fd.n_candidate_program_terms > 0).sum()),
            "n_folds_with_8_selected_programs": int((fd.n_selected_programs >= 2).sum()),
            "mean_fdr05_terms": fd.n_fdr05_terms.mean(),
            "mean_candidate_program_terms": fd.n_candidate_program_terms.mean(),
            "mean_selected_programs": fd.n_selected_programs.mean(),
            "n_terms_recurrent_across_2plus_folds": int((recurrence.n_folds >= 2).sum()) if len(recurrence) else 0,
            "n_terms_recurrent_across_3plus_folds": int((recurrence.n_folds >= 3).sum()) if len(recurrence) else 0,
        }])
    else:
        summary = pd.DataFrame()
    summary.to_csv(OUT / "05_stage2911_summary.csv", index=False)

    log("complete")
    print("\nStage 2.9.11 fold diagnostics:", flush=True)
    print(fd.to_string(index=False), flush=True)
    print("\nStage 2.9.11 overall:", flush=True)
    print(summary.to_string(index=False), flush=True)
    if len(recurrence):
        print("\nMost recurrent candidate terms:", flush=True)
        print(recurrence.head(15).to_string(index=False), flush=True)
    return summary

if __name__ == "__main__":
    run()
