# Z6 GSE67462 multimodal validation protocol

## Purpose

The GSE67462 expression-only audits now show that the positive within-system Z6 result is predominantly carried by a reproducible common dynamic component. The next question is whether that component corresponds to biologically interpretable molecular remodeling rather than a generic expression trajectory.

GSE67462 is paired with GSE67520 from the same secondary-reprogramming program. GEO describes GSE67520 as genome-wide profiling of Oct4 binding, H3K4me1, H3K27ac, H3K4me3, H3K27me3 and RNAPII across the reprogramming series. GSE67462 contains the corresponding bulk expression time series with two samples per time point. The two accessions therefore provide an appropriate within-program multimodal validation target, but not an independent biological system.

## Frozen scientific question

Does the expression-derived common dynamic component identified in GSE67462 show concordant temporal structure in independent molecular modalities from the same reprogramming program?

## Modalities

Primary:

- GSE67520 total Oct4 ChIP-seq;
- GSE67520 H3K4me1;
- GSE67520 H3K27ac;
- GSE67520 H3K4me3;
- GSE67520 H3K27me3;
- GSE67520 RNAPII.

Optional, if processed files are available locally:

- 3xFlag Oct4 ChIP-seq.

## Time alignment

Use the common time grid supported by both expression and ChIP-seq datasets:

`day0, day1, day3, day5, day7, day11, day15, day18, iPSC`

Do not interpolate missing molecular measurements.

## Validation strategy

1. Recover gene-associated or promoter/enhancer-associated modality scores using the processed GSE67520 files where available.
2. Aggregate the expression-derived common dynamic component into gene sets or ranked gene scores.
3. Test temporal concordance between expression dynamics and each chromatin/Oct4 modality.
4. Separate monotonic concordance from stage-specific concordance; do not require every modality to be monotonic.
5. Compare observed concordance with time-preserving permutation nulls.
6. Report effect sizes and empirical p-values; do not treat p-values alone as validation.

## Required controls

- time-label permutation preserving the observed marginal distribution;
- endpoint-only comparison versus full trajectory comparison;
- random gene-set control matched for gene-set size;
- if modality coverage permits, promoter versus distal/enhancer stratification;
- explicit missingness/coverage report.

## Interpretation hierarchy

- `MULTIMODAL_DYNAMIC_SUPPORT`: concordance with at least two biologically distinct modalities across intermediate timepoints and significant against time-preserving nulls.
- `PARTIAL_MULTIMODAL_SUPPORT`: concordance in one modality or endpoint-dominated concordance.
- `EXPRESSION_ONLY`: no reproducible multimodal concordance.
- `NOT_RUN_DATA_UNAVAILABLE`: required processed modality data are not locally available.

## What this does not establish

A positive result is still not an independent validation of biological specificity because GSE67462 and GSE67520 belong to the same experimental program. It supports mechanistic coherence within the system. Independent-system validation remains necessary for transferability.

## Frozen Z6 rule

This audit must not modify the Phase 0 or Phase 0.1 predictive-support decision and must not be used to tune the predictive models.
