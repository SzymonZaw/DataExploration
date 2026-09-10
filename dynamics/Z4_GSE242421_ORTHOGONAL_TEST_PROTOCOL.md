# Z4 — GSE242421 orthogonal frozen-module test

## Purpose

Test whether the six frozen GSE67462 temporal modules have concordant temporal behavior in an independently derived human scATAC gene-activity representation from GSE242421.

This is an orthogonal molecular corroboration test, not a clinical endpoint and not a causal test.

## Frozen inputs

- GSE67462 module membership is read from the existing frozen module-assignment table.
- Human↔mouse mapping remains the frozen one-to-one mapping already audited for Z4.
- GSE242421 gene activity is read from `results/Dynamics/z4_gse242421_gene_activity/01_gene_activity_log1p_cpm.tsv.gz` and is not refit here.
- No target-data feature selection, module redefinition, threshold tuning or sign selection is permitted.

## Analysis

For each frozen module:

1. Map its frozen mouse genes to human genes using the frozen one-to-one mapping.
2. Intersect with GSE242421 gene activity.
3. Compute the unweighted mean activity per timepoint using the frozen gene membership.
4. Evaluate temporal Spearman correlation across D0,D2,D4,D6,D8,D10,D12,D14. iPSC is retained as an endpoint and reported separately because it is qualitatively distinct from the continuous induction series.
5. Compare the observed temporal statistic against a prespecified time-label permutation null.
6. Compare the target direction with the frozen GSE67462 module direction supplied by the reference trajectory file.

## Null

The target scores remain fixed and the observed time labels are permuted. For each module, 5,000 permutations are used by default. The direction-specific empirical p-value is the fraction of null statistics at least as extreme in the frozen reference direction, with +1 pseudocount calibration.

## Decision logic

`Z4_ORTHO_SUPPORTED`: all evaluated modules have adequate coverage, the frozen direction is available, and the aggregate module-level evidence passes the prespecified direction-specific permutation threshold without relying on iPSC alone.

`Z4_ORTHO_PARTIAL`: adequate coverage and valid null, but only a subset of modules passes.

`Z4_ORTHO_UNRESOLVED`: frozen reference direction or sample/time provenance is unavailable, or the null cannot be calibrated.

`Z4_ORTHO_FAILED`: adequate coverage, frozen directions and valid null, but no module provides directionally concordant evidence and the result is systematically discordant.

## Boundary

A positive result strengthens cross-modality biological-state interpretation. It does not establish biological specificity, causality, or predictive transferability to an outcome.
