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

# Seurat v5 may expose an existing but empty `data` layer. Treat an empty
# layer as unavailable and fall back explicitly to counts. Align cells by
# cell names rather than positional logical indexing because assay layers can
# have a different cell universe/order from the object metadata.
layer_names <- Layers(obj[[assay]])
preferred <- intersect(c("data", "counts"), layer_names)
mat <- NULL
used_layer <- NULL
for (layer in preferred) {
  candidate <- tryCatch(LayerData(obj, assay=assay, layer=layer), error=function(e) NULL)
  if (!is.null(candidate) && nrow(candidate) > 0L && ncol(candidate) > 0L) {
    mat <- candidate
    used_layer <- layer
    break
  }
}
if (is.null(mat)) stop("No non-empty RNA assay layer found; checked: ", paste(preferred, collapse=", "))

common_cells <- intersect(colnames(mat), rownames(meta))
if (length(common_cells) == 0L) {
  stop("No overlapping cell names between assay layer and Seurat metadata")
}
mat <- mat[, common_cells, drop=FALSE]
meta <- meta[common_cells, , drop=FALSE]

sample_col <- if ("orig.ident" %in% colnames(meta)) "orig.ident" else colnames(meta)[1]
samples <- as.character(meta[[sample_col]])
keep <- !is.na(samples) & nzchar(samples)
mat <- mat[, keep, drop=FALSE]
meta <- meta[keep, , drop=FALSE]
samples <- samples[keep]

sample_levels <- unique(samples)
message("Seurat assay: ", assay)
message("layer: ", used_layer)
message("cells: ", ncol(mat), "; genes: ", nrow(mat), "; samples: ", length(sample_levels))

# Sparse sample pseudobulk mean expression.
idx <- split(seq_along(samples), samples)
pb <- matrix(0, nrow=nrow(mat), ncol=length(idx), dimnames=list(rownames(mat), names(idx)))
for (s in names(idx)) pb[, s] <- Matrix::rowMeans(mat[, idx[[s]], drop=FALSE])

write.csv(as.data.frame(pb), out_expr, quote=FALSE)

meta_out <- data.frame(sample=names(idx), n_cells=vapply(idx, length, integer(1)), stringsAsFactors=FALSE)
# Preserve donor/time/state fields when they are constant within sample.
for (col in intersect(c("age_group","age_ident","cell_state","PartialReprog1","NonReprog1","EarlyPluripotency1","Pluripotency1","orig.ident"), colnames(meta))) {
  vals <- tapply(as.character(meta[[col]]), samples, function(x) paste(unique(x), collapse=";"))
  meta_out[[col]] <- unname(vals[meta_out$sample])
}
write.csv(meta_out, out_meta, row.names=FALSE, quote=FALSE)
