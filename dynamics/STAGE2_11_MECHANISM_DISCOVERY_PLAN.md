# Stage 2.11 — Perturbation-conditioned mechanism discovery

## Purpose

This stage converts the trajectory proof-of-feasibility into a falsifiable mechanism-discovery experiment. It is deliberately not a claim of causal discovery. The objective is to test whether reproducible temporal structure can be represented by biologically interpretable pathway and transcription-factor activity states and whether perturbation identity improves prediction of future state.

## Current evidence entering Stage 2.11

The Yamanaka trajectory POC now reproduces the same 0/3/7/10-day structure in two independent human fibroblast datasets:

- GM00731: 26,852 genes, 4 day groups.
- HFIB_COMBINED: 26,382 genes, 8 sample groups.
- Temporal-distance-profile correlation between datasets: **0.921**.

This is descriptive evidence of reproducible state change. It does not establish lineage, causality, or a mechanism.

## Scientific question

Does a low-dimensional, biologically interpretable activity state preserve enough information to predict future reprogramming state across independent experiments, and does explicit perturbation information explain residual differences between trajectories?

## Working model

For dataset `d` and time `t`:

`x_d(t) -> z_d(t)`

`z_d(t+Δt) = F(z_d(t), u_d(t), h_d(t), θ_d)`

where:

- `x` = measured molecular data;
- `z` = latent biological state;
- `u` = perturbation/intervention state (OSK, OSKM, OCT4, SOX2, etc.);
- `h` = optional history/memory state;
- `θ` = context-specific parameters;
- `F` = learned dynamical transition function.

## Experimental sequence

### 1. Representation layer

Construct pathway/TF activity matrices from the existing harmonized gene expression data using the already established PROGENy and DoRothEA pipeline. Keep raw-gene representations as a baseline rather than discarding them.

### 2. Perturbation audit

Create an explicit perturbation table per sample/group. Do not infer perturbation identity from the outcome. Missing or ambiguous perturbations remain `OTHER/UNKNOWN` and are excluded from perturbation-specific claims.

### 3. Forecasting benchmark

Compare:

1. persistence;
2. nearest-time baseline;
3. linear transition model;
4. state-only neural dynamics;
5. state + perturbation dynamics;
6. state + perturbation + short history model.

Use leave-one-dataset-out evaluation as the primary generalization test and leave-one-replicate-out as a secondary stability test.

### 4. Mechanism candidates

For transitions with reproducible predictive signal, rank pathway/TF activities by their contribution to the predicted state change. A candidate mechanism must satisfy all of:

- reproducible direction across independent datasets;
- temporal ordering consistent with the candidate transition;
- predictive contribution above matched nulls;
- stability across random seeds;
- biological interpretability;
- a falsifiable downstream prediction.

### 5. Null tests

The mechanism claim is rejected if performance survives only because of dataset identity, static condition labels, leakage, or arbitrary temporal ordering. Required nulls include temporal permutation and perturbation-label permutation.

## AI / LLM role

LLMs are used only after quantitative candidate generation. For each candidate, the AI layer should produce:

`candidate -> evidence FOR -> evidence AGAINST -> predicted consequence -> discriminating experiment`

The LLM is not treated as ground truth and cannot promote a candidate to mechanism without quantitative evidence.

## Success criteria

Stage 2.11 is successful if at least one of the following is demonstrated on held-out data:

- perturbation-conditioned dynamics outperform state-only dynamics and persistence;
- a pathway/TF activity component consistently predicts a future transition across independent datasets;
- a candidate mechanism generates a quantitative prediction that can be tested experimentally.

A negative result is also scientifically useful: it identifies which representation or context variable is insufficient for cross-experiment dynamical generalization.

## Guardrails

- No cell-lineage inference from pseudobulk trajectories.
- No causal language from observational time-course data alone.
- No mechanism claim based on a single dataset.
- No tuning on the held-out dataset.
- Digital Biological Twin remains a final demonstrator, not the core evidence criterion.
