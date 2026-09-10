suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Matrix))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i+1]]
}
rds <- get_arg("--rds", "Data/GSE297234_HFIB_COMBINED_SEVOSKM.rds")
out_expr <- get_arg("--output-expression", "results/Dynamics/z4_gse297234_frozen_transfer_audit/01_gse297234_pseudobulk_human.csv")
out_meta <- get_arg("--output-metadata", "results/Dynamics/z4_gse297234_frozen_transfer_audit/02_gse297234_sample_metadata.csv")

dir.create(dirname(out_expr), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(out_meta), recursive=TRUE, showWarnings=FALSE)

obj <- readRDS(rds)
meta <- obj[[]]
assay <- if ("RNA" %in% names(obj@assays)) "RNA" else DefaultAssay(obj)

# Prefer normalized RNA data. If unavailable, fall back to counts and let the
# Python audit operate on the resulting sample-level matrix.
mat <- tryCatch(GetAssayData(obj, assay=assay, layer="data"), error=function(e) NULL)
if (is.null(mat)) mat <- GetAssayData(obj, assay=assay, layer="counts")

sample_col <- if ("orig.ident" %in% colnames(meta)) "orig.ident" else colnames(meta)[1]
samples <- as.character(meta[[sample_col]])
keep <- !is.na(samples) & nzchar(samples)
mat <- mat[, keep, drop=FALSE]
samples <- samples[keep]

sample_levels <- unique(samples)
message("Seurat assay: ", assay)
message("cells: ", ncol(mat), "; genes: ", nrow(mat), "; samples: ", length(sample_levels))

# Sparse sample pseudobulk mean expression.
idx <- split(seq_along(samples), samples)
pb <- matrix(0, nrow=nrow(mat), ncol=length(idx), dimnames=list(rownames(mat), names(idx)))
for (s in names(idx)) pb[, s] <- Matrix::rowMeans(mat[, idx[[s]], drop=FALSE])

write.csv(as.data.frame(pb), out_expr, quote=FALSE)

meta_out <- data.frame(sample=names(idx), n_cells=vapply(idx, length, integer(1)), stringsAsFactors=FALSE)
# Preserve donor/time/state fields when they are constant within sample.
for (col in intersect(c("age_group","age_ident","cell_state","PartialReprog1","NonReprog1","EarlyPluripotency1","Pluripotency1","orig.ident"), colnames(meta))) {
  vals <- tapply(as.character(meta[[col]][keep]), samples, function(x) paste(unique(x), collapse=";"))
  meta_out[[col]] <- unname(vals[meta_out$sample])
}
write.csv(meta_out, out_meta, row.names=FALSE, quote=FALSE)
