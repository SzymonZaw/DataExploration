# Z6 GSE67462 — Mechanistic interpretation audit

## Scientific purpose

This audit asks whether the multimodal dynamic signal identified in GSE67462/GSE67520 can be reduced to a coherent set of genes and temporal regulatory patterns that are plausible candidates for mechanistic interpretation.

It is a **diagnostic interpretation layer**, not a causal test. It does not modify the frozen Z6 predictive benchmark or the multimodal support thresholds.

## Evidence hierarchy

The audit separates four claims:

1. **Multimodal coherence** — the same genes show concordant temporal behavior across several independent molecular layers.
2. **Regulatory ordering** — regulatory measurements at an earlier sampled time are associated with subsequent expression changes more strongly than with contemporaneous expression.
3. **Functional enrichment** — multimodally coherent genes are enriched in an externally supplied pathway/gene-set collection.
4. **Specificity** — candidate genes are not explained equally well by predefined nuisance/context gene sets.

Only (1) is expected from the current GSE67462/GSE67520 evidence. Claims (2)–(4) require additional evidence and are reported separately.

## Core analysis

For every gene present in the validated GPL19972 → symbol → TSS namespace:

- calculate expression-vs-modality Spearman correlation for H3K27ac, H3K4me3, RNAPII, OCT4 and H3K27me3;
- count the number of supported active modalities with positive correlation;
- calculate a robust multimodal score as the median positive correlation across supported active modalities;
- classify genes as `MULTIMODAL_CORE`, `MULTIMODAL_SUPPORT`, `SINGLE_MODAL`, or `UNRESOLVED`;
- calculate endpoint concordance separately from full-trajectory concordance;
- calculate an adjacent-time regulatory-lead diagnostic where the regulatory modality at time `t` is compared with expression at the next sampled timepoint.

The lead analysis is explicitly labelled **hypothesis-generating** because the experiment contains only eight main timepoints and unequal temporal intervals.

## Mechanistic interpretation rules

A gene is a `MULTIMODAL_CORE` candidate when it has positive temporal correlation with at least three of the four supported active modalities and its median correlation exceeds 0.30. This is a ranking rule, not a statistical significance claim.

A gene is a `REGULATORY_LEAD_CANDIDATE` only when an active regulatory modality has a stronger adjacent-time correlation with subsequent expression than with contemporaneous expression by at least 0.10. No causal direction is inferred.

H3K27me3 is retained as a negative-control/contrasting modality and is never counted toward the active multimodal core score.

## Functional enrichment

The executable audit accepts an optional GMT file. If supplied, it performs a simple overlap enrichment against the ranked `MULTIMODAL_CORE` gene set using a hypergeometric test with Benjamini–Hochberg correction. The background is the validated expression universe, not the whole genome.

No external pathway database is silently downloaded. This keeps the analysis reproducible and prevents database-version drift. Recommended future inputs include versioned PROGENy/DoRothEA-derived gene sets or a versioned MSigDB GMT.

## Nuisance/context specificity

The executable accepts an optional second GMT file containing nuisance/context gene sets (for example proliferation, interferon/stress, apoptosis or generic cell-cycle sets). Enrichment is reported as a diagnostic comparison only. Absence of enrichment is not evidence of biological specificity.

## Required non-claims

The audit must not state that:

- OCT4 causes the expression transition;
- enhancer/promoter remodeling proves causal regulation;
- a pathway is activated merely because its genes are enriched;
- the candidate genes are universal markers of cellular state;
- the GSE67462 mechanism transfers to GSE297234 or another system.

## Reproducibility

The audit consumes the same validated local inputs used by the multimodal robustness analysis:

- GSE67462 common-space expression;
- GPL19972 SOFT annotation;
- `mm9.refGene.gtf.gz` derived from the local UCSC `refGene.txt.gz` annotation;
- GSE67520 regulatory peak files.

The output is intended to be machine-readable and manuscript-ready, with a compact Markdown findings report and CSV tables containing the ranked candidate genes.
