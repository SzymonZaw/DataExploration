# Z6 GSE67462 — directional mechanistic ordering audit

## Purpose

Test whether the multimodal GSE67462/GSE67520 signal contains a reproducible directional temporal ordering in which regulatory/chromatin-associated measurements at time `t` are associated with subsequent expression at `t+1`, rather than merely contemporaneous correlation.

This is a diagnostic mechanistic audit. It does not modify frozen Z6 predictive support, multimodal support, temporal modules, or validation splits.

## Hypotheses

- H-D1: Active regulatory modalities (H3K27ac, H3K4me3, RNAPII, OCT4) show stronger lagged regulatory→expression concordance than contemporaneous concordance.
- H-D2: The lagged signal is concentrated in mechanistic core genes and the principal modules M1/M4/M6 rather than being diffuse across all genes.
- H-D3: The direction of the dominant regulatory trajectory is concordant with the subsequent expression trajectory for a substantial fraction of candidate genes.
- H-D4: The H3K27me3 negative-control modality should not show the same systematic lagged pattern.

## Operational definitions

For each gene and modality, using the validated common gene universe:

- contemporaneous rho: Spearman(`modality(t)`, `expression(t)`)
- lagged rho: Spearman(`modality(t)`, `expression(t+1)`)
- lead gain: `lagged_rho - contemporaneous_rho`
- lead candidate: lead gain >= 0.10 and lagged rho >= 0.30
- directional agreement: sign of modality trajectory slope equals sign of subsequent-expression trajectory slope

The slope sign is descriptive; it is not treated as causal evidence.

## Stability checks

The audit reports:

1. candidate counts by module and modality;
2. fraction of module genes represented by candidates;
3. fraction of candidates with >=3 active modalities;
4. directional agreement;
5. active-modality versus H3K27me3 negative-control comparison;
6. M1/M4/M6 summary.

## Scientific boundary

A positive result supports a temporal ordering hypothesis, not causal precedence. Independent perturbation/intervention data remain required to establish necessity or causality.
