# Z6 GSE67462 temporal-module enrichment protocol

## Scientific purpose

Test whether the reproducible temporal modules identified in the GSE67462 multimodal core have statistically supported functional enrichment, without selecting genes manually or assigning pathway names from trajectory shape alone.

This is a **mechanistic interpretation / hypothesis-generation audit**. It does not modify frozen Z6 predictive-support decisions or the multimodal-support threshold.

## Input

- `03_gene_module_assignments.csv` from the temporal-module audit.
- The provenance-validated GSE67462 common-space universe (11,048 genes in the current analysis).
- Explicitly versioned GMT collections for:
  - GO Biological Process,
  - Reactome,
  - KEGG (optional),
  - nuisance/context gene sets (optional).

The GMT files must be supplied by the user and retained with the analysis provenance. The script never silently downloads pathway databases.

## Tested modules

The primary report evaluates modules with at least 20 genes. In the current result this excludes M2 (4 genes) and M5 (10 genes), which are too small for stable pathway inference.

Primary modules:

- M1: 196 genes
- M3: 170 genes
- M4: 258 genes
- M6: 497 genes

Small modules remain in the machine-readable assignment output but are not used for the primary enrichment conclusions.

## Statistical test

For each module and gene set, use a one-sided hypergeometric over-representation test:

`P(X >= k)`,

where the background is the full provenance-validated expression universe, not the whole genome and not the union of selected module genes.

For each database collection:

1. restrict every gene set to the validated background,
2. require a minimum overlap of 5 genes for reporting,
3. calculate the hypergeometric p-value,
4. apply Benjamini-Hochberg FDR across all module × term tests in that collection,
5. retain the raw p-value, global collection FDR, and module-local FDR.

## Interpretation rules

A pathway is a **candidate enriched process** when it has adequate overlap and FDR < 0.05. Biological interpretation must also consider:

- module temporal shape,
- enrichment direction/absence in neighboring modules,
- redundancy among highly overlapping terms,
- nuisance/context enrichment,
- whether the process is expected for reprogramming independently of this dataset.

An enriched term is not evidence that the pathway causes the transition.

## Nuisance/context audit

If nuisance gene sets are supplied, they are tested against the same validated background. Examples include broad proliferation/cell-cycle, generic stress, interferon/innate immune response, apoptosis, hypoxia, and housekeeping/ribosomal programs.

Nuisance enrichment is not automatically a failure. It determines whether a mechanistic interpretation needs to be framed as a generic cellular response rather than a transition-specific program.

## Multiple testing

The primary FDR correction is performed across all module × term hypotheses within each collection. GO, Reactome and KEGG are treated as separate collections because their term universes and redundancy structures differ.

No uncorrected p-value is used as a biological conclusion.

## Negative control

H3K27me3 remains a contrast modality from the upstream multimodal analysis. It is not used as a gene-set background or silently removed from the expression universe.

## Reproducibility

The exact GMT files, their source/version metadata, SHA-256 hashes, background size, module sizes, thresholds and script commit must be recorded in the generated manifest.

The enrichment audit is diagnostic only and does not change the frozen Z6 predictive or multimodal support results.
