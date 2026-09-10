# Z4 GSE242421 analysis-product audit

## Purpose

The GSE242421 preflight passed. The next step is provenance and format verification of the published analysis products before any orthogonal validation is run.

Do **not** download the 17.3 GB raw GEO tar for this gate. The published analysis-product record provides a 1.6 GB `scATAC.zip` and a 58.6 MB `scATAC_scRNA_integration.zip`, which are sufficient for the intended audit path.

## Frozen inputs

The following remain frozen and must not be refit on GSE242421:

- GSE67462 expression-derived module membership
- frozen human↔mouse orthology mapping
- module direction definitions
- training-time feature selection and scaling rules

## Required published products

From Zenodo record `10.5281/zenodo.8313962`:

### `scATAC.zip`

Required members:

- `cells.tsv`
- `peaks.bed`
- `features.tsv`
- `cell_x_peak.mtx.gz`

The published metadata states that `cells.tsv` contains sample/time-point labels, `peaks.bed` contains the common 500 bp peak set, and `cell_x_peak.mtx.gz` contains the sparse fragment-count matrix with rows matching peaks and columns matching cells.

### `scATAC_scRNA_integration.zip`

Required members:

- `peak_gene_links_fdr1e-4.tsv`
- `harmony.cca.30.feat.tsv`
- `harmony.cca.metadata.tsv`

The peak-gene links are useful for a frozen mapping from genes to regulatory peaks without training a new gene-activity model on the validation experiment.

## What is intentionally not accepted as an independent endpoint

The published scATAC cell clusters, UMAP coordinates, and pseudotime are analysis products derived from the same chromatin dataset. They can be used as contextual diagnostics, but they are not independent validation endpoints.

Likewise, a generic correlation between the frozen modules and the first chromatin principal component is not sufficient for Z4 support.

## Next decision gate

If both archives contain the required members, proceed to a deterministic gene/module-activity derivation:

1. verify the nine expected samples (`D0,D2,D4,D6,D8,D10,D12,D14,iPSC`);
2. verify hg38 peak coordinates;
3. map frozen human module genes to published peak-gene links;
4. aggregate accessibility over linked peaks per frozen module;
5. test temporal concordance against the frozen expression-module trajectory;
6. use a time-label permutation null while keeping the frozen module definitions unchanged.

A gene-activity matrix should therefore be **derived from published chromatin products**, not downloaded as an opaque fitted representation.

## Source

Published analysis products: Zenodo `10.5281/zenodo.8313962`.
