# Z6 Phase 0.1 — Within-system transferability audit

## Purpose

Phase 0 showed no frozen predictive support under leave-one-dataset-out (LODO) validation. This audit separates two hypotheses:

- **H-Z6a:** the current common-state representation contains insufficient information for future-state prediction even within a biological system.
- **H-Z6b:** predictive information exists within a biological system but does not transfer across datasets/platforms/organisms.

The audit therefore keeps the model, preprocessing, free-rollout forecast, baselines and predictive-support rule unchanged, while replacing the LODO unit with an independent branch/replicate within the same dataset.

## Datasets and held-out branches

Only datasets with at least two independently identifiable temporal branches are used:

- **GSE28688:** replicate branches `a` and `b` from the paired HFF1/OSKM samples; 4 timed points per branch (0, 24, 48, 72 h).
- **GSE67462:** replicate branches `1` and `2`; 8 timed points per branch (0, 24, 72, 120, 168, 264, 360, 432 h).
- **GSE297234:** donor branches `aged` (GM00731) and `young` (GM23815); 4 timed points per donor (0, 72, 168, 240 h). These are treated explicitly as **donor holdout**, not technical replication.

GSE297234 donor assignment follows the GEO sample design: GM00731 and GM23815 are the two fibroblast donors sampled at days 0, 3, 7 and 10.

## Validation design

For each dataset separately:

1. Construct one trajectory per independent branch.
2. Train on all branches except the held-out branch.
3. Apply the same training-only feature selection, imputation and scaling as Phase 0.
4. Observe the first 60% of the held-out trajectory.
5. Perform genuine free rollout over all remaining timepoints.
6. Compare against persistence, nearest-time and linear extrapolation baselines.
7. Apply the same permutation test: held-out observations are permuted while the time vector remains fixed.
8. Do not tune thresholds or select the better branch after observing results.

## Frozen predictive-support rule

A model has predictive support only if all of the following hold:

- mean improvement versus persistence > 0;
- q05 improvement versus persistence > 0;
- mean improvement versus nearest-time > 0;
- q05 improvement versus nearest-time > 0;
- mean improvement versus linear > 0;
- q05 improvement versus linear > 0;
- permutation p < 0.05.

The PCA evaluator uses the corrected inverse transformation: PCA inverse followed by the benchmark scaler inverse, so predictions and observations are evaluated in the same space.

## Interpretation matrix

| Within-system result | LODO result | Interpretation |
|---|---|---|
| Positive | Negative | Supports H-Z6b: predictive information may exist but is not cross-dataset transferable. |
| Negative | Negative | Supports H-Z6a: current representation/dynamics lack robust predictive information even within systems. |
| Positive | Positive | Strongest evidence that the representation supports transferable prediction. |
| Mixed | Negative | Representation likely contains context-dependent predictive information; inspect dataset/branch dependence before model expansion. |

A positive result is not treated as mechanistic evidence. It establishes predictive utility under the specified validation design only.

## Non-goals

This audit does not:

- remove difficult datasets from Phase 0;
- change the frozen thresholds;
- optimize architecture/hyperparameters against the held-out branch;
- claim causality;
- replace the LODO benchmark;
- establish a universal biological trajectory.

The output is a diagnostic for interpreting the negative LODO result and deciding whether a representation-ablation study is scientifically justified.
