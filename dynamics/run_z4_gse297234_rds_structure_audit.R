#!/usr/bin/env Rscript

# Structural audit of local GSE297234 Seurat RDS objects.
# Diagnostic only: no normalization, integration, feature selection, transfer score,
# orthology mapping, or endpoint scoring is performed here.

suppressPackageStartupMessages({
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)

# Canonical defaults: the normal workflow requires no arguments.
default_rds <- c(
  file.path("Data", "GSE297234_GM00731_SEVOSKM.rds"),
  file.path("Data", "GSE297234_HFIB_COMBINED_SEVOSKM.rds")
)
default_output <- file.path(
  "results", "Dynamics", "z4_gse297234_rds_structure_audit", "structure.json"
)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("Missing", flag))
  args[[i + 1]]
}

# Optional overrides are retained for reproducibility/debugging.
out_path <- get_arg("--output", default_output)
rds_idx <- which(args == "--rds")
if (length(rds_idx) > 1) stop("Only one --rds option may be supplied")

if (length(rds_idx) == 1) {
  out_idx <- which(args == "--output")
  end_idx <- if (length(out_idx) == 1) out_idx - 1 else length(args)
  rds_paths <- args[(rds_idx + 1):end_idx]
  if (length(rds_paths) < 1) stop("At least one RDS path is required")
} else {
  rds_paths <- default_rds
}

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
  selected_names <- intersect(c(
    "orig.ident", "age_group", "age_ident", "AGE_UP1", "AGE_DOWN1",
    "cell_state", "ordered_clusters", "PartialReprog1", "NonReprog1",
    "EarlyPluripotency1", "Pluripotency1", "HALLMARK_EMT1", "HALLMARK_TGFB1",
    "batch", "treatment", "day"
  ), meta_cols)
  selected_meta <- if (length(selected_names)) meta[, selected_names, drop = FALSE] else NULL
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

# Concise console output; the complete machine-readable report is in structure.json.
cat("GSE297234 Z4 RDS STRUCTURE AUDIT\n")
for (obj in results) {
  cat(sprintf(
    "%s | %s | %d features x %d cells | assays: %s | default: %s\n",
    basename(obj$path), paste(obj$class, collapse = ";"), obj$n_features, obj$n_cells,
    paste(obj$assays, collapse = ";"), obj$default_assay
  ))
}
cat(sprintf("output: %s\n", out_path))
cat("transfer score: NOT EVALUATED\n")
cat("orthology mapping: NOT EVALUATED\n")
cat("endpoint: NOT EVALUATED\n")
