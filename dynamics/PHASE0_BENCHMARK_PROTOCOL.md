# Phase 0 — Frozen Predictive Benchmark Protocol

## Scientific question

> Does the raw gene-level common space contain a representation in which experimental temporal ordering provides predictive information that generalizes to an unseen dataset?

The benchmark is a **representation-level falsification experiment**, not a claim that any particular neural architecture must succeed.

## Models

All five candidates use the same held-out-dataset splits:

1. PCA + latent transition
2. Static autoencoder + latent Ridge transition
3. Markov DynamicStateModel
4. Delta-t conditioned DynamicStateModel
5. History-aware DynamicStateModel with GRU memory

The latent dimension is currently fixed at 8 for Phase 0. Changing it is a separate representation/dimensionality ablation and must not be mixed into this frozen benchmark.

## Data splitting

Evaluation uses leave-one-dataset-out (LODO). For every held-out dataset:

- feature selection is fit only on training datasets;
- median imputation is fit only on training datasets;
- scaling statistics are fit only on training datasets;
- PCA is fit only on training datasets;
- neural models are trained only on training datasets;
- the held-out trajectory is used only as an observed prefix and future ground truth.

No future held-out observation is fed back into a model rollout.

## Forecasting task

The observed prefix is approximately 60% of the ordered held-out trajectory, with at least two observed points. The model then performs a genuine free rollout:

`observed prefix -> predicted state -> predicted state -> ...`

The true future trajectory is never used to update the model state.

For the memory model, history is updated with predicted observations only.

## Baselines

Every model is evaluated against the same three baselines in the same transformed observation space:

- **Persistence:** future state equals the final observed prefix state.
- **Nearest-time:** mean training-dataset observation at the closest available training time.
- **Linear:** extrapolation from the last two observed prefix points.

This symmetry is required. A model is not considered superior merely because it beats one selected baseline.

## Metrics

For every fold/model combination we record:

- future-point count;
- RMSE of the model;
- RMSE of persistence, nearest-time and linear baselines;
- improvement over each baseline;
- training loss where applicable.

The summary reports mean, median and 5th/95th percentile improvements.

## Temporal permutation null

For each fitted model, the held-out observations are randomly permuted while the held-out time vector is kept fixed. The fitted model is then evaluated using the same free-rollout procedure.

The permutation p-value tests whether the observed temporal ordering produces a larger improvement over persistence than the null distribution.

The null is not a test of whether the model can reconstruct the data. It is a test of whether the predictive advantage depends on temporal correspondence.

## Predictive-support criterion

A model receives `predictive_support=True` only if all of the following hold:

1. mean improvement over persistence > 0;
2. 5th percentile improvement over persistence > 0;
3. mean improvement over nearest-time > 0;
4. 5th percentile improvement over nearest-time > 0;
5. mean improvement over linear > 0;
6. 5th percentile improvement over linear > 0;
7. permutation p < 0.05.

This is deliberately conservative. `p < 0.05` is necessary but not sufficient.

## Interpretation of failure

If all models fail, the conclusion is **not** that cellular dynamics do not exist.

The narrower conclusion is:

> Temporal generalization has not been demonstrated in the raw gene-level common space under the frozen benchmark protocol.

Possible explanations remain open, including inadequate biological representation, platform/context heterogeneity, insufficient temporal coverage, cell-state heterogeneity, or limited identifiability from the available datasets.

## What is frozen

For the official Phase 0 result, do not change the following after inspecting results:

- dataset list and LODO splits;
- preprocessing order;
- baseline definitions;
- free-rollout rule;
- permutation definition;
- support criterion;
- primary interpretation.

Implementation bugs and methodological inconsistencies may still be corrected. Any such correction requires rerunning the complete benchmark and replacing the official results.

## Next phase

Only after Phase 0 is frozen should the project move to **Phase 1 — Ablation of biological state representation**:

`genes -> pathways -> TF activity -> single-cell latent -> multimodal latent`

The same forecasting benchmark should be reused so that representation, rather than benchmark changes, is the experimental variable.

LLM reasoning, mechanistic hypothesis generation, Neural ODE/SDE and Digital Biological Twin work are intentionally outside Phase 0.
