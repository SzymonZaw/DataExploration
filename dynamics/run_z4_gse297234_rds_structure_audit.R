#!/usr/bin/env Rscript

# Structural audit of local GSE297234 Seurat RDS objects.
# Diagnostic only: no normalization, integration, feature selection, transfer score,
# orthology mapping, or endpoint scoring is performed here.

suppressPackageStartupMessages({
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript dynamics/run_z4_gse297234_rds_structure_audit.R --rds <file1> [file2 ...] --output <json>")
}

get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[[i + 1]]
}

rds_idx <- which(args == "--rds")
out_path <- get_arg("--output")
if (length(rds_idx) != 1) stop("Exactly one --rds argument is required")
out_idx <- which(args == "--output")
rds_paths <- args[(rds_idx + 1):(out_idx - 1)]
if (length(rds_paths) < 1) stop("At least one RDS path is required")

json_escape <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub('"', '\\"', x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x
}

json_array <- function(x) paste0("[", paste(sprintf('"%s"', vapply(x, json_escape, character(1))), collapse = ","), "]")

summarise_object <- function(path) {
  if (!file.exists(path)) stop(paste("Missing RDS:", path))
  x <- readRDS(path)
  meta <- x[[]]
  assays <- Assays(x)
  feature_examples <- list()
  for (a in assays) {
    feature_examples[[a]] <- head(rownames(x[[a]]), 20)
  }
  meta_cols <- colnames(meta)
  selected_meta <- meta[, intersect(c(
    "orig.ident", "age_group", "age_ident", "AGE_UP1", "AGE_DOWN1",
    "cell_state", "ordered_clusters", "PartialReprog1", "NonReprog1",
    "EarlyPluripotency1", "Pluripotency1", "HALLMARK_EMT1", "HALLMARK_TGFB1",
    "batch", "treatment", "day"
  ), meta_cols), drop = FALSE]
  list(
    path = path,
    class = class(x),
    dims = dim(x),
    assays = assays,
    default_assay = DefaultAssay(x),
    n_cells = ncol(x),
    n_features = nrow(x),
    feature_examples = feature_examples,
    metadata_columns = meta_cols,
    selected_metadata = selected_meta
  )
}

results <- lapply(rds_paths, summarise_object)

# JSON is deliberately produced through jsonlite if available; otherwise stop rather than
# silently emitting an invalid machine-readable artifact.
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required. Install it once with install.packages('jsonlite').")
}

out <- list(
  candidate = "GSE297234",
  purpose = "RDS_STRUCTURE_AUDIT_ONLY",
  transfer_score_evaluated = FALSE,
  orthology_mapping_evaluated = FALSE,
  endpoint_evaluated = FALSE,
  objects = results
)

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(out, out_path, pretty = TRUE, auto_unbox = TRUE, dataframe = "rows", na = "null")
cat(jsonlite::toJSON(out, pretty = TRUE, auto_unbox = TRUE, dataframe = "rows", na = "null"))
cat("\n")
