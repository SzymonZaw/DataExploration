# Z4 GSE242421 peak→gene-link robustness audit

## Purpose

Test whether the orthogonal corroboration of frozen GSE67462 modules 4 and 6 in GSE242421 depends on the primary peak→gene-link threshold used to derive gene activity.

GSE242421 is an independent human fibroblast reprogramming scATAC time series with D0–D14 and iPSC samples. Published analysis products provide the scATAC matrix and peak→gene links. The present audit reuses the already derived sample-level peak CPM matrix and changes only the absolute peak→gene correlation threshold. urlNCBI GSE242421https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242421 urlZenodo analysis productshttps://zenodo.org/records/8313962

## Frozen inputs

- `results/Dynamics/z4_gse242421_gene_activity/02_peak_sample_cpm.tsv.gz`
- `Data/GSE242421/scATAC_scRNA_integration.zip`
- `results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv`
- `results/Dynamics/z6_gse67462_temporal_modules/02_temporal_modules.csv`
- `results/Dynamics/z4_frozen_orthology_mapping/01_frozen_human_mouse_mapping.csv`

No feature selection, module fitting, or target-derived direction inference is allowed.

## Thresholds

Primary threshold: `|correlation| >= 0.45`.

Robustness thresholds: `0.50`, `0.60`, `0.70`.

The binary peak→gene aggregation rule is kept identical to the primary derivation. Only links passing the stricter threshold are retained.

## Test

For modules 4 and 6, which have frozen `late_rising` direction in GSE67462:

1. retain genes mapped through the frozen human↔mouse mapping;
2. aggregate linked peak CPM to gene activity;
3. aggregate module gene activity into a module score;
4. correlate D0–D14 module score with ordered time;
5. use 5,000 time-label permutations for a directional null;
6. report iPSC separately and do not let it determine the continuous-time result.

## Decision

`Z4_ORTHO_ROBUST` only if **both modules 4 and 6 pass at every threshold**.

Otherwise: `Z4_ORTHO_THRESHOLD_SENSITIVE`.

This is a robustness audit, not an independent validation endpoint. Even a robust result remains `ORTHOGONAL_MOLECULAR_CORROBORATION_ONLY` because GSE242421 has one sample per time point.
