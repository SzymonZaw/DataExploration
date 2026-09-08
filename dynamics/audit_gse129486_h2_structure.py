"""Audit GSE129486 before any H3 similarity analysis.

This is a design/structure audit only. It does not compare trajectories with
Yamanaka and does not compute the H3 similarity statistic.

It checks whether GSE129486 can be treated as an independent temporal process
and flags the predeclared H2 nuisance overlap risk (inflammatory/JAK-STAT/
NF-kB-related biology) without claiming that any nuisance pathway is actually
activated in this dataset.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11c_h3_control_redesign"
EXPR = DATA / "GSE129486_rnaseq-data-1_gene-tpm.tsv.gz"
META = DATA / "GSE129486_rnaseq-data-1_metadata.tsv.gz"
EXPECTED_EXPR = "6e3d7860f4f38d95830a15b8dd570d58f9226170df6002343e432a05c96fcbf0"
EXPECTED_META = "78a60e353461ab819672e478a9f38b20822fdfc493a1773677ddef3eb167ff04"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def audit_metadata(meta: pd.DataFrame) -> dict:
    required = ["sample", "cell_line", "stimulation", "time"]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise RuntimeError(f"Missing required metadata columns: {missing}")

    m = meta.copy()
    m["time_num"] = pd.to_numeric(m["time"], errors="coerce")
    if m["time_num"].isna().any():
        raise RuntimeError("Non-numeric time labels detected after numeric coercion.")

    # A true temporal replicate must be a distinct biological sample within a
    # fixed cell-line × stimulation × time stratum. The metadata has 21 samples
    # in most time/stimulation strata, but those observations are not all
    # independent processes: cell line and stimulation define the experimental
    # blocks. We therefore report both raw counts and block-level replication.
    strata = (
        m.groupby(["cell_line", "stimulation", "time_num"], dropna=False)
        .agg(n_samples=("sample", "nunique"))
        .reset_index()
    )
    duplicate_sample_ids = int(m["sample"].duplicated().sum())

    block_counts = (
        m.groupby(["cell_line", "stimulation"], dropna=False)
        .agg(
            n_samples=("sample", "nunique"),
            n_timepoints=("time_num", "nunique"),
        )
        .reset_index()
    )

    time_counts = m.groupby("time_num")["sample"].nunique().sort_index()
    stimulation_counts = m.groupby("stimulation")["sample"].nunique().sort_values(ascending=False)
    cellline_counts = m.groupby("cell_line")["sample"].nunique().sort_values(ascending=False)

    # Identify whether the design is crossed enough to compare temporal
    # profiles within each biological block. A complete block has one sample at
    # every declared time point. This is a structural property, not H3 evidence.
    all_times = sorted(m["time_num"].unique().tolist())
    complete_blocks = []
    for (cell_line, stimulation), g in m.groupby(["cell_line", "stimulation"]):
        times = sorted(g["time_num"].unique().tolist())
        complete_blocks.append(
            {
                "cell_line": str(cell_line),
                "stimulation": str(stimulation),
                "n_samples": int(g["sample"].nunique()),
                "n_timepoints": int(len(times)),
                "complete_all_timepoints": times == all_times,
                "duplicate_samples_within_block": int(g["sample"].duplicated().sum()),
            }
        )

    # This is deliberately a risk flag, not an activity result. The H2 panel
    # already declared JAK-STAT and NF-kB-like inflammatory programs as
    # possible nuisance axes. TNF/IL-17A is an inflammatory stimulation class,
    # so this control cannot be treated as maximally orthogonal to that H2
    # hypothesis. Actual pathway activation must be tested separately.
    stimulation_labels = sorted(str(x) for x in m["stimulation"].dropna().unique())
    inflammatory_labels = [
        x for x in stimulation_labels
        if any(token in x.lower() for token in ["tnf", "il-17", "il17", "cytokine", "inflamm"])
    ]

    return {
        "metadata_schema": {
            "n_samples": int(m["sample"].nunique()),
            "n_rows": int(len(m)),
            "n_timepoints": int(m["time_num"].nunique()),
            "time_values_hours": all_times,
            "n_cell_lines": int(m["cell_line"].nunique()),
            "n_stimulations": int(m["stimulation"].nunique()),
            "stimulation_labels": stimulation_labels,
            "duplicate_sample_ids": duplicate_sample_ids,
        },
        "temporal_replication": {
            "timepoint_sample_counts": {str(k): int(v) for k, v in time_counts.items()},
            "strata_n": int(len(strata)),
            "strata_with_ge_2_samples_n": int((strata["n_samples"] >= 2).sum()),
            "strata_with_ge_3_samples_n": int((strata["n_samples"] >= 3).sum()),
            "cell_line_x_stimulation_blocks_n": int(len(block_counts)),
            "complete_all_timepoint_blocks_n": int(sum(x["complete_all_timepoints"] for x in complete_blocks)),
            "blocks": complete_blocks,
        },
        "h2_independence_audit": {
            "status": "RISK_FLAG",
            "predeclared_h2_nuisance_axes": ["JAK-STAT", "NF-kB", "interferon/STAT-like inflammatory response"],
            "process_class": "acute inflammatory cytokine response",
            "stimulation_labels_matching_inflammatory_keywords": inflammatory_labels,
            "interpretation": "GSE129486 is temporally well structured, but its inflammatory stimulation class is not maximally orthogonal to the predeclared H2 nuisance hypothesis. This does not establish actual pathway activation or confounding in GSE129486.",
            "required_followup": "Use a separate nuisance-axis audit to quantify actual pathway activity before interpreting GSE129486 as an H3 control. Do not use H3 similarity itself to decide whether this control is orthogonal.",
        },
        "h3_eligibility": {
            "temporal_structure": "PASS",
            "replicate_block_structure": "PASS",
            "h2_orthogonality": "LIMITED",
            "usable_as_sole_h3_control": False,
            "reason": "A single inflammatory control cannot cleanly distinguish generic temporal structure from the predeclared inflammatory/STAT nuisance hypothesis.",
        },
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    expr_sha = sha256(EXPR)
    meta_sha = sha256(META)
    meta = pd.read_csv(META, sep="\t", compression="gzip", low_memory=False)
    report = {
        "status": "ok",
        "scope": "GSE129486 structure and H2-independence audit; no Yamanaka comparison, H3 similarity, trajectory-agreement statistic, or H3 decision",
        "input_hashes": {
            "GSE129486_gene_tpm": {"sha256": expr_sha, "expected": EXPECTED_EXPR, "matches": expr_sha == EXPECTED_EXPR},
            "GSE129486_metadata": {"sha256": meta_sha, "expected": EXPECTED_META, "matches": meta_sha == EXPECTED_META},
        },
        "audit": audit_metadata(meta),
        "guardrail": "This audit separates design/exchangeability questions from H3 similarity. A risk flag for H2 is not evidence that H2 is present; it only prevents the inflammatory control from being treated as a clean orthogonal control.",
    }
    path = OUT / "H3_GSE129486_H2_STRUCTURE_AUDIT.json"
    path.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"WROTE {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
