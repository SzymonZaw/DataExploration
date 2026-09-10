# Z6 GSE67462/GSE67520 regulatory-element audit

## Purpose

The first multimodal validation used a deliberately conservative promoter-only mapping (peak midpoint within +/-2 kb of one representative TSS). It yielded only 2-6 genes per modality and therefore is not an adequate biological falsification test.

This audit tests whether the GSE67462 common dynamic expression component has temporal concordance with GSE67520 regulatory signal when distal regulatory elements are allowed. It is diagnostic and does not alter the frozen Z6 predictive-support rule.

## Data

- Expression: GSE67462, the same two-branch common dynamic representation used in the signal-attribution audit.
- Regulatory data: GSE67520 peak calls for total Oct4, H3K4me1, H3K27ac, H3K4me3, H3K27me3 and RNAPII.
- Genome: mouse mm9 / NCBI37.
- Annotation: local `Data/GSE67520/mm9.refGene.gtf.gz`, generated from the UCSC mm9 `refGene.txt.gz` table.

UCSC documents RefSeq/refGene as one of the available gene-model sources and notes that different gene-model sets have different transcript coverage. Therefore the annotation source and assembly are explicitly retained in the audit provenance.

## Temporal alignment

GSE67520 labels are mapped to the GSE67462 hour grid without interpolation:

- d0 -> 0 h
- d1 -> 24 h
- d3 -> 72 h
- d5 -> 120 h
- d7 -> 168 h
- d11 -> 264 h
- d15 -> 360 h
- d18 -> 432 h

## Regulatory mappings

Four representations are compared:

1. `promoter_2kb`: peak midpoint within +/-2 kb of a representative TSS.
2. `nearest_tss_25kb`: nearest TSS, maximum distance 25 kb.
3. `nearest_tss_50kb`: nearest TSS, maximum distance 50 kb.
4. `nearest_tss_100kb`: nearest TSS, maximum distance 100 kb.

A peak is assigned to at most one gene. This is a sensitivity analysis, not a claim that every distal peak regulates the assigned nearest gene.

The key question is whether concordance is stable as the mapping radius expands, rather than whether one arbitrary radius produces significance.

## Statistics

For each modality and mapping radius:

- number of mapped peaks;
- number of genes represented;
- median and mean gene-level temporal Spearman correlation;
- global flattened temporal Spearman correlation;
- 1000 time-label permutations preserving the modality values and expression time grid;
- permutation p-value and null 95th percentile.

A mapping is considered `SUPPORTED` only if:

- at least 50 genes are shared;
- observed global Spearman > permutation null q95;
- permutation p < 0.05.

This threshold is deliberately a data-adequacy gate, not a Z6 support criterion. It prevents a result based on a handful of genes from being promoted to mechanistic evidence.

## Interpretation

- `ROBUST_MULTIMODAL_SUPPORT`: at least two modalities are supported at two or more mapping radii, with support including at least one distal-capable radius (25/50/100 kb).
- `DISTAL_SENSITIVE_SUPPORT`: at least two modalities are supported only after allowing distal mappings.
- `PARTIAL_REGULATORY_SUPPORT`: exactly one modality satisfies the support rule, or support is inconsistent across radii.
- `NO_REGULATORY_SUPPORT`: no modality satisfies the support rule despite adequate gene counts.
- `MAPPING_INADEQUATE`: no modality/radius reaches the 50-gene adequacy gate.

## What this does not establish

Positive concordance is mechanistic coherence within the same OSKM reprogramming program, not independent biological validation. It does not establish causal regulation, enhancer-gene assignment correctness, or cross-system transferability.

A negative result after adequate distal mapping is more informative than the previous promoter-only negative result, but still does not by itself falsify the predictive expression representation.
