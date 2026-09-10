# Z4 GSE242421 orthogonal validation — findings

## Status

**Decision: `Z4_ORTHO_THRESHOLD_SENSITIVE`**

The primary GSE242421 test showed orthogonal corroboration for the two frozen modules with an a-priori unambiguous directional label (modules 4 and 6). However, the prespecified peak→gene link threshold robustness audit shows that this corroboration is **threshold-sensitive** and therefore cannot be treated as a robust external validation result.

Modules 1, 2 and 5 are labelled `early_transient`, and module 3 is labelled `mixed`; these labels do not define a single monotonic direction and were not converted into artificial +/− directions.

## Primary orthogonal test

| Module | Frozen label | Validated human genes | Coverage | GSE242421 Spearman rho | Frozen direction | Directional p |
|---|---|---:|---:|---:|---:|---:|
| 1 | early_transient | 160 | 81.6% | −0.190 | not directional | — |
| 2 | early_transient | 4 | 100.0% | +0.071 | not directional | — |
| 3 | mixed | 117 | 68.8% | +0.476 | not directional | — |
| 4 | late_rising | 218 | 84.5% | **+0.881** | +1 | **0.0024** |
| 5 | early_transient | 8 | 80.0% | +0.357 | not directional | — |
| 6 | late_rising | 395 | 79.5% | **+1.000** | +1 | **0.0002** |

The continuous temporal test used D0–D14 (8 ordered timepoints). iPSC was retained as a separate terminal observation and was not allowed to determine the continuous-time result.

## Peak→gene link robustness

The published GSE242421 peak–gene resource was re-analysed with absolute correlation thresholds of 0.45, 0.50, 0.60 and 0.70. The published analysis itself describes filtering FDR < 1e-4 links to absolute correlation > 0.45. citeturn0search2

| Link threshold | Module | Validated human genes | Coverage | rho D0–D14 | Direction concordant | Directional p | Pass |
|---:|---:|---:|---:|---:|---|---:|---|
| 0.45 | 4 | 218 | 84.5% | −0.119 | **No** | 0.6331 | No |
| 0.45 | 6 | 395 | 79.5% | +0.619 | **Yes** | 0.0576 | No |
| 0.50 | 4 | 204 | 79.1% | −0.071 | **No** | 0.5825 | No |
| 0.50 | 6 | 355 | 71.4% | +0.690 | **Yes** | 0.0368 | **Yes** |
| 0.60 | 4 | 167 | 64.7% | +0.143 | **Yes** | 0.3695 | No |
| 0.60 | 6 | 277 | 55.7% | +0.881 | **Yes** | 0.0042 | **Yes** |
| 0.70 | 4 | 112 | 43.4% | +0.214 | **Yes** | 0.3085 | No |
| 0.70 | 6 | 187 | 37.6% | +0.905 | **Yes** | 0.0022 | **Yes** |

Thus module 4 does not provide directional support at any tested threshold, while module 6 passes at 0.50–0.70 but fails the primary 0.45 threshold because its directional permutation p-value is 0.0576. The aggregate result therefore fails the frozen robustness criterion at every threshold: **0/4 thresholds pass both modules**.

This is not evidence that the underlying biological signal is absent. It shows that the external chromatin corroboration is sensitive to the definition of the peak→gene link set, so the current result is not sufficiently stable to promote to robust Z4 validation.

## Module 6 biological interpretation audit

An interpretation-only enrichment audit of frozen module 6 found concordant biological themes in both the original mouse module and its frozen validated human projection. Significant terms included skin/epithelium/tissue development, epithelial differentiation, keratinization/cornified-envelope processes, apoptosis-related cell-adhesion cleavage, and cholesterol biosynthesis. The human validated module reproduced the core skin/epithelial themes, with 28/28 genes in the reported `skin development` intersection and 12/12 in the reported cornified-envelope and keratinization intersections.

The target temporal core contained 109 genes but did not yield significant enrichment terms under the same audit. This does **not** contradict the module-level interpretation: the core was deliberately defined as the top 25% of frozen validated module genes by signed D0–D14 Spearman in the target and is therefore a selected temporal subset, while the enrichment test is descriptive only.

The enrichment audit does not change the Z4 decision and must not be used to retroactively select, refit, or justify module 6.

## Interpretation

The strongest defensible statement is now:

> In an independent human fibroblast reprogramming experiment, the frozen late-rising module set showed partial orthogonal molecular corroboration in scATAC-derived gene activity, with module 6 showing stable directional behavior under stricter peak–gene link thresholds. However, the result is sensitive to the prespecified link threshold: module 6 does not pass at the primary 0.45 threshold and module 4 does not pass at any tested threshold. Therefore the orthogonal corroboration is supportive but not robust enough to claim full external validation of the frozen representation.

GSE242421 is an independent human fibroblast OSKM reprogramming series with scATAC-seq at D0, D2, D4, D6, D8, D10, D12, D14 and iPSC. citeturn0search0turn0search2

## Methodological safeguards

- Frozen GSE67462 module membership was not refit on GSE242421.
- Frozen mouse→human 1:1 orthology was used.
- No target-data feature selection was performed for the module test.
- Gene activity was derived from the published GSE242421 peak–gene links rather than fitted to the target data.
- The temporal test used a 5,000-permutation time-label null.
- iPSC was reported separately from the continuous D0–D14 test.
- Modules without a pre-existing directional interpretation were not assigned an artificial direction.
- Link-threshold robustness was evaluated without changing frozen module membership or orthology.
- Module 6 enrichment was interpretation-only and did not alter the frozen decision.

## Scientific conclusion

The appropriate thesis-level statement is:

> The independent GSE242421 chromatin-accessibility experiment provides partial orthogonal corroboration for the frozen biological-state representation, particularly for a late-rising module with coherent epithelial/tissue-developmental annotation. Nevertheless, sensitivity to the peak–gene link threshold and the absence of replicate timepoints prevent claiming robust external validation or universal transferability.

`Z4_ORTHO_THRESHOLD_SENSITIVE` remains the frozen decision.