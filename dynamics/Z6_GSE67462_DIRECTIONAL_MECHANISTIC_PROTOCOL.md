# Z6 GSE67462 directional mechanistic audit

## Purpose

Test whether regulatory/chromatin trajectories show a temporal lead over expression trajectories within the already validated GSE67462/GSE67520 multimodal framework.

This is a **diagnostic mechanistic audit**. It does not alter frozen Z6 predictive support, temporal module assignments, validation splits, or multimodal-support decisions.

## Biological question

For each temporal module M1, M4 and M6, ask whether a regulatory signal at time `t` is more concordant with expression at the next sampled timepoint than with expression at the same timepoint:

`lead_gain = rho(regulator(t), expression(t+1)) - rho(regulator(t), expression(t))`

The analysis is motivated by the original GSE67462 study, which reports stage-specific Oct4 binding and histone-mark/enhancer signatures concordant with expression changes, including enhancer activation preceding later pluripotency-network activation. citeturn0search0turn0search1

## Modalities

Active modalities:

- H3K27ac
- H3K4me3
- RNAPII
- OCT4

Negative-control modality:

- H3K27me3

The negative control is not assumed to be completely independent of expression. H3K27me3 is treated as a repressive mark, so its expected directional relationship is opposite to active marks.

## Tests

For every validated gene in M1, M4 and M6:

1. calculate concurrent Spearman correlation across sampled timepoints;
2. calculate one-step lead Spearman correlation between regulatory signal at `t` and expression at `t+1`;
3. calculate `lead_gain`;
4. calculate directional agreement between regulatory change at `t` and expression change at `t+1`, using positive expected direction for active modalities and negative expected direction for H3K27me3;
5. compare the observed module-level lead gain with a gene-wise circular time-shift null.

The null preserves each gene's regulatory trajectory but destroys its temporal alignment with expression. Shift 0 is excluded from the null.

## Interpretation

A positive lead gain and strong directional agreement support a **temporal-order hypothesis**, not causality.

Evidence is stronger when:

- the lead effect is positive at module level;
- the lead gain exceeds the circular-shift null;
- the effect is present across multiple active modalities;
- the effect is reproducible across the three modules rather than driven by one module;
- H3K27me3 behaves differently from the active marks or provides a coherent inverse relationship.

A failure to show lead gain does not prove absence of regulation because the experiment has sparse sampling and bulk measurements.

## Multiple-testing boundary

This audit is hypothesis-generating. Per-modality/module permutation p-values are reported descriptively and are not used to modify frozen Z6 support.

## Required outputs

- `01_gene_directional_effects.csv`
- `02_modality_summary.csv`
- `03_module_directional_summary.csv`
- `04_permutation_summary.csv`
- `05_manifest.json`
- `GSE67462_Z6_directional_mechanistic_report.md`

## Scientific boundary

The strongest claim supported by this audit is that a regulatory trajectory is **temporally ordered with subsequent expression change**. It is not evidence that the regulatory mark is necessary or sufficient. Perturbation remains required for causal inference.
