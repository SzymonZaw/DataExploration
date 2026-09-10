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
assay_obj <- obj[[assay]]

layer_names <- tryCatch(Layers(assay_obj), error=function(e) character(0))
message("assay class: ", paste(class(assay_obj), collapse=";"))
message("reported layers: ", if (length(layer_names)) paste(layer_names, collapse=",") else "<none>")

# GSE297234 is stored in Seurat v5 as one counts layer per GEO sample,
# e.g. counts.GM00731_D0. Therefore selecting the first non-empty layer would
# silently reduce the audit to one sample. Instead, load every non-empty
# counts/data layer and construct one pseudobulk column per layer/sample.
if (length(layer_names)) {
  data_layers <- layer_names[grepl("^data($|\\.)", layer_names, ignore.case=TRUE)]
  count_layers <- layer_names[grepl("^counts($|\\.)", layer_names, ignore.case=TRUE)]
  candidate_layers <- unique(c(data_layers, count_layers))
} else {
  candidate_layers <- character(0)
}

layer_matrices <- list()
layer_cells <- list()
for (layer in candidate_layers) {
  candidate <- tryCatch(LayerData(obj, assay=assay, layer=layer), error=function(e) NULL)
  if (!is.null(candidate) && nrow(candidate) > 0L && ncol(candidate) > 0L) {
    common_cells <- intersect(colnames(candidate), rownames(meta))
    if (length(common_cells) > 0L) {
      layer_matrices[[layer]] <- candidate[, common_cells, drop=FALSE]
      layer_cells[[layer]] <- common_cells
    }
  }
}

# Legacy Seurat Assay fallback: one matrix containing all cells.
if (length(layer_matrices) == 0L) {
  for (slot_name in c("data", "counts")) {
    candidate <- tryCatch(GetAssayData(obj, assay=assay, slot=slot_name), error=function(e) NULL)
    if (!is.null(candidate) && nrow(candidate) > 0L && ncol(candidate) > 0L) {
      common_cells <- intersect(colnames(candidate), rownames(meta))
      if (length(common_cells) > 0L) {
        layer_matrices[[paste0("legacy_slot:", slot_name)]] <- candidate[, common_cells, drop=FALSE]
        layer_cells[[paste0("legacy_slot:", slot_name)]] <- common_cells
        break
      }
    }
  }
}

if (length(layer_matrices) == 0L) {
  stop("No non-empty RNA expression matrix found; reported layers: ",
       if (length(layer_names)) paste(layer_names, collapse=",") else "<none>")
}

# For split v5 layers, preserve the layer/sample identity explicitly. Do not
# combine cells across layers because each layer already represents one GEO
# sample. For a legacy unsplit assay, use orig.ident as the sample grouping.
pseudobulk <- list()
meta_rows <- list()
for (layer in names(layer_matrices)) {
  mat <- layer_matrices[[layer]]
  cells <- layer_cells[[layer]]
  meta_layer <- meta[cells, , drop=FALSE]

  if (startsWith(layer, "legacy_slot:")) {
    sample_col <- if ("orig.ident" %in% colnames(meta_layer)) "orig.ident" else colnames(meta_layer)[1]
    samples <- as.character(meta_layer[[sample_col]])
    valid <- !is.na(samples) & nzchar(samples)
    mat <- mat[, valid, drop=FALSE]
    meta_layer <- meta_layer[valid, , drop=FALSE]
    samples <- samples[valid]
    idx <- split(seq_along(samples), samples)
    for (s in names(idx)) {
      pseudobulk[[s]] <- Matrix::rowMeans(mat[, idx[[s]], drop=FALSE])
      row <- data.frame(sample=s, n_cells=length(idx[[s]]), stringsAsFactors=FALSE)
      for (col in intersect(c("age_group","age_ident","cell_state","PartialReprog1","NonReprog1","EarlyPluripotency1","Pluripotency1","orig.ident"), colnames(meta_layer))) {
        vals <- unique(as.character(meta_layer[[col]][idx[[s]]]))
        row[[col]] <- paste(vals[!is.na(vals) & nzchar(vals)], collapse=";")
      }
      meta_rows[[s]] <- row
    }
  } else {
    sample_name <- sub("^(counts|data)\\.", "", layer, ignore.case=TRUE)
    pseudobulk[[sample_name]] <- Matrix::rowMeans(mat)
    row <- data.frame(sample=sample_name, n_cells=ncol(mat), stringsAsFactors=FALSE)
    for (col in intersect(c("age_group","age_ident","cell_state","PartialReprog1","NonReprog1","EarlyPluripotency1","Pluripotency1","orig.ident"), colnames(meta_layer))) {
      vals <- unique(as.character(meta_layer[[col]]))
      row[[col]] <- paste(vals[!is.na(vals) & nzchar(vals)], collapse=";")
    }
    meta_rows[[sample_name]] <- row
  }
}

# Build a rectangular expression table over the union of genes. Missing genes
# in a sample-specific layer are represented as zero, preserving sparse counts
# semantics without inventing normalized values.
genes <- sort(unique(unlist(lapply(pseudobulk, names))))
samples <- names(pseudobulk)
pb <- matrix(0, nrow=length(genes), ncol=length(samples), dimnames=list(genes, samples))
for (s in samples) pb[names(pseudobulk[[s]]), s] <- pseudobulk[[s]]

write.csv(as.data.frame(pb), out_expr, quote=FALSE)
meta_out <- do.call(rbind, meta_rows)
rownames(meta_out) <- NULL
write.csv(meta_out, out_meta, row.names=FALSE, quote=FALSE)

message("Seurat assay: ", assay)
message("expression sources: ", length(layer_matrices), " non-empty layer(s)")
message("cells: ", sum(vapply(layer_matrices, ncol, integer(1))), "; genes union: ", length(genes), "; samples: ", length(samples))
message("samples: ", paste(samples, collapse=", "))
