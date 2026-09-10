# Z4 GSE242421 orthogonal validation — findings

## Status

**Decision: `Z4_ORTHO_PARTIAL`**

The frozen GSE67462 module set contains six modules, but only modules 4 and 6 have a pre-existing unambiguous directional label suitable for a directional temporal test. Modules 1, 2 and 5 are labelled `early_transient`, and module 3 is labelled `mixed`; these labels do not define a single monotonic direction and are therefore not converted into an artificial +/− direction.

This is important: the two modules with an a-priori frozen direction (modules 4 and 6) both pass the independent GSE242421 chromatin-activity directional permutation test.

## Results

| Module | Frozen label | Validated human genes | Coverage | GSE242421 Spearman rho | Frozen direction | Directional p |
|---|---|---:|---:|---:|---:|---:|
| 1 | early_transient | 160 | 81.6% | −0.190 | not directional | — |
| 2 | early_transient | 4 | 100.0% | +0.071 | not directional | — |
| 3 | mixed | 117 | 68.8% | +0.476 | not directional | — |
| 4 | late_rising | 218 | 84.5% | **+0.881** | +1 | **0.0024** |
| 5 | early_transient | 8 | 80.0% | +0.357 | not directional | — |
| 6 | late_rising | 395 | 79.5% | **+1.000** | +1 | **0.0002** |

The continuous temporal test used D0–D14 (8 ordered timepoints). iPSC was retained as a separate terminal observation and was not allowed to determine the continuous-time result.

## Interpretation

The result supports **orthogonal molecular corroboration** of two frozen late-rising modules in an independent human fibroblast reprogramming experiment measured by scATAC-seq-derived gene activity. It does **not** establish that all six modules transfer, and it does not by itself establish biological specificity or a universally transferable state trajectory.

The target experiment is genuinely orthogonal at the molecular measurement level: GSE242421 is an independent human fibroblast OSKM reprogramming series with scATAC-seq at D0, D2, D4, D6, D8, D10, D12, D14 and iPSC. citeturn0search0turn0search2

## Methodological safeguards

- Frozen GSE67462 module membership was not refit on GSE242421.
- Frozen mouse→human 1:1 orthology was used.
- No target-data feature selection was performed.
- Gene activity was derived from the published GSE242421 peak–gene links rather than fitted to the target data.
- The temporal test used a 5,000-permutation time-label null.
- iPSC was reported separately from the continuous D0–D14 test.
- Modules without a pre-existing directional interpretation were not assigned an artificial direction.

## Scientific conclusion

The appropriate thesis-level statement is:

> In an independent human fibroblast reprogramming experiment, the two frozen modules with an a-priori late-rising direction showed concordant temporal dynamics in independently derived chromatin-accessibility gene activity. This provides partial orthogonal corroboration of the frozen biological-state representation, while the remaining transient/mixed modules and the absence of replicate timepoints prevent claiming full external validation or universal transferability.

`Z4_ORTHO_PARTIAL` remains the frozen decision.