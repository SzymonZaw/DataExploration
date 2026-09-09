# Z6 — Predictive modeling of cellular-state transitions

## Scientific objective

Test whether a learned biological-state representation supports genuine out-of-sample prediction of future cellular state:

\[
z(t) \rightarrow z(t+\Delta t).
\]

The primary claim tested by Z6 is predictive, not geometric: a model must forecast observations that were not used for fitting and must outperform meaningful persistence, nearest-time and linear-extrapolation baselines.

## Current implementation

Z6 reuses the frozen Phase 0 predictive benchmark because it already implements the required leakage-free forecasting design:

- leave-one-dataset-out evaluation;
- training-only feature selection, imputation and scaling;
- approximately 60% observed prefix of the held-out trajectory;
- genuine free rollout over the complete future;
- persistence, nearest-time and linear baselines;
- PCA transition, static autoencoder + Ridge transition, Markov DynamicStateModel, delta-t-conditioned DynamicStateModel and history-aware DynamicStateModel;
- temporal permutation null;
- conservative predictive-support criterion.

The benchmark is therefore the first operational experiment for Z6, not a claim that Z6 is already scientifically completed.

## Primary datasets

The current frozen benchmark uses the three trajectories already supported by the common-state pipeline:

- `GSE67462`
- `GSE28688`
- `GSE297234`

Evaluation is leave-one-dataset-out. The held-out dataset is never used for fitting preprocessing or model parameters.

## Forecasting task

For each held-out trajectory:

```text
observed prefix
      ↓
state inference
      ↓
transition model
      ↓
free multi-step rollout
      ↓
predicted future state
      ↓
comparison with held-out future observations
```

The true future is never fed back into the rollout.

## Baselines

Every model is compared against:

1. persistence — future equals the last observed state;
2. nearest-time — training observations closest to the target time;
3. linear extrapolation — extrapolation from the final two observed states.

## Primary metrics

For every held-out dataset and model report:

- future-point RMSE;
- RMSE improvement over persistence;
- RMSE improvement over nearest-time;
- RMSE improvement over linear extrapolation;
- mean, median and 5th/95th percentile improvements across folds/seeds;
- temporal permutation p-value.

## Predictive-support rule

A model receives `predictive_support=True` only when all of the following are satisfied:

- mean improvement over persistence > 0;
- 5th percentile improvement over persistence > 0;
- mean improvement over nearest-time > 0;
- 5th percentile improvement over nearest-time > 0;
- mean improvement over linear > 0;
- 5th percentile improvement over linear > 0;
- permutation p < 0.05.

This prevents a model from being called predictive because it only wins on average or only beats one weak baseline.

## Interpretation

### Positive Z6 result

A positive result supports the narrower statement:

> The tested state representation contains information that allows future biological observations to be predicted out-of-sample beyond the predefined baselines under the frozen protocol.

It does **not** by itself establish causality or a universal biological trajectory.

### Negative Z6 result

A negative result means temporal prediction has not been demonstrated for the tested representation and datasets. It does not imply that cellular dynamics are absent.

The next scientific variables are representation quality, context dependence, temporal coverage, perturbation information and history dependence.

## Relation to Z4

Z6 is evaluated independently of the unresolved H3 control question. A predictive result must subsequently be audited for nuisance/context dependence. In particular, predictive performance cannot be interpreted as biological specificity until Z4 controls are adequate.

## Execution

Run the frozen benchmark through:

```powershell
python run_z6_predictive_transition.py
```

For a fast smoke test:

```powershell
python run_z6_predictive_transition.py --epochs 25 --permutation-n 50 --seeds 411
```

The official run uses the defaults defined in the wrapper and writes results to:

`results/Dynamics/z6_predictive_transition/`

The smoke test is diagnostic only and must never be used as the official Z6 result.
