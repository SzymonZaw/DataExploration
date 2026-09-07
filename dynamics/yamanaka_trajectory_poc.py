"""Stage 2.11b: time-resolved trajectory proof-of-feasibility.

This module is deliberately conservative. It converts the local GSE297234 RDS
objects to sample-level pseudobulk through R/Seurat, then asks whether the
0/3/7/10-day OSKM trajectory is reproducible across two human fibroblast
samples and whether pathway/TF activities change monotonically or in defined
phases. It is not a cell-lineage or causal-mechanism claim.

The large RDS files remain local and are never copied into Git.
"""
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "Data"
OUT = ROOT / "results" / "Dynamics" / "stage2_11b_yamanaka_trajectory"

RDS_FILES = {
    "GM00731": DATA / "GSE297234_GM00731_SEVOSKM.rds",
    "HFIB_COMBINED": DATA / "GSE297234_HFIB_COMBINED_SEVOSKM.rds",
}


def _find_rscript() -> str:
    for name in ("Rscript", "Rscript.exe"):
        path = shutil.which(name)
        if path:
            return path
    raise RuntimeError("Rscript was not found. Install R and ensure Rscript is on PATH.")


def _write_r_helper(path: Path) -> None:
    path.write_text(r'''args <- commandArgs(trailingOnly=TRUE)
infile <- args[[1]]
outfile <- args[[2]]

# These RDS files contain Seurat objects. Seurat itself is required to
# deserialize the Seurat class; SeuratObject provides the v5 assay API.
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

obj <- readRDS(infile)

if (inherits(obj, "Seurat")) {
  meta <- obj@meta.data
  assay <- obj@active.assay
  if (is.null(assay) || !nzchar(assay)) assay <- names(obj@assays)[1]
  if (!"counts" %in% SeuratObject::Layers(obj, assay=assay)) {
    stop("No counts layer found in assay ", assay,
         ". Available layers: ", paste(SeuratObject::Layers(obj, assay=assay), collapse=","))
  }
  counts <- SeuratObject::LayerData(obj, assay=assay, layer="counts")
} else if (inherits(obj, "SingleCellExperiment")) {
  suppressPackageStartupMessages(library(SummarizedExperiment))
  meta <- as.data.frame(SummarizedExperiment::colData(obj))
  counts <- SummarizedExperiment::assay(obj, "counts")
} else {
  stop("Unsupported RDS class: ", paste(class(obj), collapse=","))
}

# Find a day/time variable conservatively. GEO/Seurat objects do not always
# use a column literally named day/time; sample identifiers such as
# "..._D3" or "day7" are common. Therefore we score every metadata column
# by how consistently its values contain one of the expected days.
extract_day <- function(x) {
  x <- as.character(x)
  m <- regmatches(x, regexpr("(^|[^0-9])([0-9]{1,3})([^0-9]|$)", x, perl=TRUE))
  if (!length(m) || !nzchar(m)) return(NA_real_)
  z <- sub("^[^0-9]*([0-9]{1,3}).*$", "\\1", m)
  d <- suppressWarnings(as.numeric(z))
  if (d %in% c(0,3,7,10)) return(d)
  NA_real_
}

score_time_column <- function(nm) {
  vals <- meta[[nm]]
  if (is.factor(vals)) vals <- as.character(vals)
  vals <- as.character(vals)
  if (!length(vals)) return(list(score=-Inf, day=rep(NA_real_, length(vals)), column=nm))
  parsed <- vapply(vals, extract_day, numeric(1))
  coverage <- mean(!is.na(parsed))
  unique_days <- length(unique(parsed[!is.na(parsed)]))
  if (coverage < 0.50 || unique_days < 2) {
    return(list(score=-Inf, day=parsed, column=nm))
  }
  low <- tolower(nm)
  name_bonus <- if (grepl("day|time|hour|dox|oskm|ident|sample|condition|group", low)) 0.25 else 0
  score <- coverage + min(unique_days, 4) * 0.05 + name_bonus
  list(score=score, day=parsed, column=nm)
}

candidates <- lapply(colnames(meta), score_time_column)
scores <- vapply(candidates, function(x) x$score, numeric(1))
if (all(!is.finite(scores))) {
  # Write a compact metadata diagnostic next to the expected output so a
  # future dataset-specific failure is immediately inspectable.
  diag <- data.frame(
    column=colnames(meta),
    class=vapply(meta, function(x) class(x)[1], character(1)),
    n_unique=vapply(meta, function(x) length(unique(x)), integer(1)),
    example=vapply(meta, function(x) paste(head(unique(as.character(x)), 5), collapse=" || "), character(1)),
    stringsAsFactors=FALSE
  )
  write.csv(diag, paste0(outfile, ".metadata_diagnostic.csv"), row.names=FALSE, quote=TRUE)
  stop("Could not identify a metadata column encoding at least two of day 0/3/7/10")
}

best <- candidates[[which.max(scores)]]
time_col <- best$column
day <- best$day

# Pseudobulk within each unique metadata condition/time group. For this POC we
# preserve donor/cell-source columns where available to keep groups separable.
meta$.__day <- day
candidate_group <- c("sample", "orig.ident", "donor", "patient", "condition", "group", "cell_source")
group_col <- candidate_group[candidate_group %in% colnames(meta)]
if (!length(group_col)) group_col <- time_col
meta$.__group <- do.call(paste, c(meta[group_col], sep="|"))
valid <- !is.na(meta$.__day)
meta <- meta[valid, , drop=FALSE]
counts <- counts[, valid, drop=FALSE]

groups <- unique(meta$.__group)
out <- vector("list", length(groups))
for (i in seq_along(groups)) {
  idx <- which(meta$.__group == groups[i])
  vec <- Matrix::rowSums(counts[, idx, drop=FALSE])
  out[[i]] <- vec
}
mat <- do.call(cbind, out)
colnames(mat) <- groups
rownames(mat) <- rownames(counts)

df <- data.frame(gene=rownames(mat), as.matrix(mat), check.names=FALSE)
write.csv(df, outfile, row.names=FALSE, quote=FALSE)

info <- data.frame(group=groups, day=vapply(groups, function(g) meta$.__day[match(g, meta$.__group)], numeric(1)), stringsAsFactors=FALSE)
write.csv(info, paste0(outfile, ".meta.csv"), row.names=FALSE, quote=FALSE)
''', encoding="utf-8")


def _convert_rds(ds: str, rds: Path, tmp: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not rds.exists():
        raise FileNotFoundError(f"Missing local RDS for {ds}: {rds}")
    rscript = _find_rscript()
    helper = tmp / "extract_gse297234.R"
    _write_r_helper(helper)
    out_csv = tmp / f"{ds}.csv"
    cmd = [rscript, str(helper), str(rds), str(out_csv)]
    # R packages may emit locale-specific text on Windows. Decode explicitly
    # and tolerate one malformed byte so the real R error is still reported.
    proc = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if proc.returncode != 0:
        stderr = (proc.stderr or "").strip()
        stdout = (proc.stdout or "").strip()
        detail = stderr or stdout or f"Rscript exited with code {proc.returncode}"
        raise RuntimeError(f"RDS extraction failed for {ds}: {detail}")
    if not out_csv.exists():
        raise RuntimeError(f"RDS extraction produced no output CSV for {ds}")
    X = pd.read_csv(out_csv, index_col=0)
    meta = pd.read_csv(str(out_csv) + ".meta.csv")
    X.index = X.index.astype(str).str.upper().str.replace(r"\.\d+$", "", regex=True)
    X = X.groupby(level=0).sum()
    return X, meta


def _log_cpm(X: pd.DataFrame) -> pd.DataFrame:
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0).clip(lower=0)
    lib = X.sum(axis=0).replace(0, np.nan)
    return np.log1p(X.div(lib, axis=1) * 1e6)


def _trajectory_metrics(X: pd.DataFrame, meta: pd.DataFrame, ds: str) -> tuple[pd.DataFrame, dict]:
    X = _log_cpm(X)
    meta = meta.copy()
    meta["day"] = pd.to_numeric(meta["day"], errors="coerce")
    meta = meta.dropna(subset=["day"])
    meta["day"] = meta["day"].astype(float)
    meta = meta[meta.day.isin([0, 3, 7, 10])]
    cols = [c for c in meta.group.astype(str) if c in X.columns]
    meta = meta[meta.group.astype(str).isin(cols)].copy()
    meta = meta.sort_values(["day", "group"])
    X = X.loc[:, meta.group.astype(str).tolist()]

    # Training-only feature selection is not meaningful for a single dataset,
    # so select genes by variance across the observed trajectory only and label
    # this explicitly as descriptive rather than predictive.
    var = X.var(axis=1).sort_values(ascending=False)
    keep = var.head(min(2000, len(var))).index
    Z = StandardScaler().fit_transform(X.loc[keep].T)
    pca = PCA(n_components=min(3, Z.shape[0], Z.shape[1]), random_state=0)
    coords = pca.fit_transform(Z)
    traj = meta.reset_index(drop=True).copy()
    for j in range(coords.shape[1]):
        traj[f"PC{j+1}"] = coords[:, j]
    explained = pca.explained_variance_ratio_.tolist()

    # Summaries by day, preserving replicates rather than averaging them away.
    day_summary = traj.groupby("day", as_index=False)[[c for c in traj.columns if c.startswith("PC")]].mean()
    rows = []
    for pc in [c for c in traj.columns if c.startswith("PC")]:
        y = day_summary[pc].to_numpy(float)
        t = day_summary.day.to_numpy(float)
        rho = pd.Series(t).corr(pd.Series(y), method="spearman") if len(y) >= 3 else np.nan
        rows.append({"dataset": ds, "component": pc, "spearman_day": float(rho) if pd.notna(rho) else np.nan})
    metrics = pd.DataFrame(rows)
    traj.to_csv(OUT / f"{ds}_trajectory_coordinates.csv", index=False)
    day_summary.to_csv(OUT / f"{ds}_day_centroids.csv", index=False)
    return metrics, {"n_genes": int(X.shape[0]), "n_groups": int(X.shape[1]), "days": sorted(meta.day.unique().tolist()), "pca_explained": explained}


def _replication(coords: dict[str, pd.DataFrame]) -> dict:
    # Compare ordered pairwise distances between consecutive day centroids when
    # both datasets contain the same days. This avoids claiming aligned absolute
    # coordinates across independently fitted PCAs.
    if len(coords) < 2:
        return {"status": "insufficient_datasets"}
    names = sorted(coords)
    a, b = names[:2]
    A, B = coords[a], coords[b]
    common_days = sorted(set(A.day) & set(B.day))
    if len(common_days) < 3:
        return {"status": "insufficient_common_days", "common_days": common_days}
    # Compare within-dataset temporal distance profiles after normalizing each
    # profile; this is invariant to arbitrary PCA rotation/scale.
    profiles = {}
    for name, df in ((a, A), (b, B)):
        df = df[df.day.isin(common_days)].sort_values("day")
        pcs = [c for c in df.columns if c.startswith("PC")]
        d = np.linalg.norm(np.diff(df[pcs].to_numpy(float), axis=0), axis=1)
        profiles[name] = d / (np.linalg.norm(d) or 1.0)
    corr = float(np.corrcoef(profiles[a], profiles[b])[0,1]) if len(profiles[a]) >= 2 else np.nan
    return {"status": "ok", "dataset_a": a, "dataset_b": b, "common_days": common_days, "temporal_distance_profile_correlation": corr}


def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="gse297234_", dir=str(ROOT)) as td:
        tmp = Path(td)
        audits, metrics, centroids = [], [], {}
        for ds, rds in RDS_FILES.items():
            X, meta = _convert_rds(ds, rds, tmp)
            m, audit = _trajectory_metrics(X, meta, ds)
            metrics.append(m)
            centroids[ds] = pd.read_csv(OUT / f"{ds}_day_centroids.csv")
            audits.append({"dataset": ds, **audit})
        metric_df = pd.concat(metrics, ignore_index=True)
        metric_df.to_csv(OUT / "01_trajectory_metrics.csv", index=False)
        replication = _replication(centroids)
        summary = {
            "status": "ok",
            "datasets": audits,
            "trajectory_metrics": metric_df.to_dict(orient="records"),
            "replication": replication,
            "interpretation_guardrail": "Trajectory structure is descriptive evidence of reproducible state change; it is not evidence of cell lineage, causality, or a discovered mechanism.",
            "next_step": "Add pathway/TF activity and perturbation-conditioned forecasting only after trajectory reproducibility is established.",
        }
        (OUT / "02_poc_summary.json").write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
        return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.parse_args()
    print(json.dumps(run(), indent=2, default=str))
