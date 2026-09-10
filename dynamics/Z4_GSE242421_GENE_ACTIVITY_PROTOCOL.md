# Z4 — GSE242421 orthogonal gene-activity derivation

## Purpose

Construct an independent chromatin-accessibility representation for GSE242421 without fitting or redefining any frozen GSE67462 module.

## Frozen inputs

- GSE67462 frozen module membership and direction.
- GSE67462 feature-selection and human↔mouse orthology decisions.
- GSE242421 published peak set and peak→gene links.

## Derivation

The published `scATAC.zip` contains a sparse peak-by-cell matrix, `peaks.bed`, and `cells.tsv`. The published integration product contains `peak_gene_links_fdr1e-4.tsv`.

1. Read the sparse peak×cell matrix without densifying it.
2. Resolve each cell to one of D0,D2,D4,D6,D8,D10,D12,D14,iPSC using the published `cells.tsv` sample field.
3. Aggregate accessibility across cells within each sample. This produces a 9×peak matrix and avoids constructing a large dense cell×gene matrix.
4. Parse the published peak→gene links.
5. Apply the prespecified paper-level correlation filter `abs(correlation) >= 0.45` when a correlation column is available. If the link file has no correlation column, use all FDR 1e-4 links and record that fact explicitly.
6. Construct a binary peak→gene incidence matrix. Gene activity is the sum of accessibility over linked peaks; no expression data and no frozen module information are used in this step.
7. Log-transform the sample-level gene activity as `log1p(CPM)` for comparability across samples.
8. Write the gene×sample matrix, mapping/provenance summary, and coverage diagnostics.

## Scientific boundary

This is an orthogonal molecular representation, not an independent clinical or phenotypic endpoint. A later positive result can strengthen interpretation of the frozen biological-state representation but cannot establish causality by itself.

## Leakage prevention

No GSE242421 feature selection, module redefinition, module threshold tuning, or direction learning is permitted. The only GSE242421-specific operations are deterministic parsing, aggregation, and the prespecified peak→gene mapping.

## Expected output

`results/Dynamics/z4_gse242421_gene_activity/`

- `01_gene_activity_log1p_cpm.tsv.gz`
- `02_peak_sample_cpm.tsv.gz`
- `03_peak_gene_mapping_summary.json`
- `04_derivation_summary.json`

The output is intentionally sample-level because the first Z4 orthogonal test concerns temporal molecular-state concordance, not cell-level prediction.
