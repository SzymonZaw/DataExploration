# Z6 — Results and scientific findings

## Status

This document records the current empirical findings for objective Z6: testing whether the learned biological-state representation supports genuine out-of-sample prediction of future cellular state.

The results below are **versioned experimental findings**, not a claim that Z6 is scientifically complete.

## 1. Experimental logic

Z6 is evaluated as a forecasting problem:

```text
observed trajectory prefix
        ↓
state inference
        ↓
transition model
        ↓
free multi-step rollout
        ↓
predicted future state
        ↓
held-out future observation
```

The frozen predictive-support rule requires all of the following:

- positive mean improvement versus persistence;
- positive q05 improvement versus persistence;
- positive mean and q05 improvement versus nearest-time;
- positive mean and q05 improvement versus linear extrapolation;
- temporal permutation p < 0.05.

No dataset was removed from the primary Phase 0 benchmark because of a negative result, and the Phase 0.1 audit does not change the model, preprocessing, baselines or support thresholds.

## 2. Phase 0 — leave-one-dataset-out benchmark

### Result

The frozen LODO benchmark did **not** establish predictive support for any tested model.

| Model | RMSE | Improvement vs persistence | q05 vs persistence | Improvement vs nearest | q05 vs nearest | Improvement vs linear | q05 vs linear | permutation p | Support |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| autoencoder | 0.901629 | 0.159065 | -0.052136 | 0.125055 | 0.033247 | 1.166434 | 0.435963 | 0.992008 | false |
| delta_t | 0.846834 | 0.213860 | -0.046771 | 0.179850 | 0.107191 | 1.221228 | 0.441327 | 0.992008 | false |
| markov | 0.843220 | 0.217474 | -0.058432 | 0.183464 | 0.113749 | 1.224843 | 0.429667 | 0.979021 | false |
| memory | 0.913102 | 0.147592 | -0.174806 | 0.113582 | 0.016180 | 1.154961 | 0.313292 | 0.996004 | false |
| pca | 1.107466 | -0.046772 | -0.316259 | -0.080782 | -0.417815 | 0.960597 | 0.256811 | 1.000000 | false |

### Interpretation

Phase 0 is a valid negative result under strict leave-one-dataset-out validation. It does **not** show that cellular dynamics are absent and does not show that Z6 is impossible. It shows that the current representation and models do not demonstrate cross-dataset future-state prediction under the frozen protocol.

The result motivates separating two hypotheses:

- **H-Z6a:** the representation contains insufficient predictive information even within a biological system.
- **H-Z6b:** predictive information exists within a biological system but does not transfer across independent datasets/platforms/organisms.

## 3. Phase 0.1 — within-system transferability audit

### Purpose

Phase 0.1 changes only the validation unit. Instead of training on all but one dataset, the audit trains on one independent temporal branch and predicts the other branch within each dataset. The model, preprocessing, free-rollout task, baselines, permutation test and support rule remain frozen.

The three eligible systems were:

- **GSE28688:** paired HFF1/OSKM branches `a` and `b`, 4 timepoints each (0, 24, 48, 72 h).
- **GSE67462:** replicate branches `1` and `2`, 8 timepoints each (0, 24, 72, 120, 168, 264, 360, 432 h).
- **GSE297234:** donor branches `aged` and `young`, 4 timepoints each (0, 72, 168, 240 h). These are treated as donor holdouts, not technical replicates.

### Results by system

| System | Autoencoder | Delta-t | Markov | Memory | PCA |
|---|---|---|---|---|---|
| GSE28688 | false | false | false | false | false |
| GSE67462 | **true** | false | false | false | **true** |
| GSE297234 | false | false | false | false | false |

### GSE28688

No model reached predictive support. The models nevertheless showed positive mean improvements over persistence, but none passed the complete frozen criterion, including the required lower-tail and permutation constraints.

### GSE67462

Two representations reached full predictive support:

**Autoencoder**

- mean improvement vs persistence: **0.408949**
- q05 vs persistence: **0.343518**
- mean improvement vs nearest-time: **0.111096**
- q05 vs nearest-time: **0.016143**
- mean improvement vs linear: **2.767563**
- q05 vs linear: **2.551851**
- permutation p: **0.000999**

**PCA**

- mean improvement vs persistence: **0.456518**
- q05 vs persistence: **0.447309**
- mean improvement vs nearest-time: **0.158665**
- q05 vs nearest-time: **0.129416**
- mean improvement vs linear: **2.815133**
- q05 vs linear: **2.665125**
- permutation p: **0.000999**

Both therefore satisfy the frozen `predictive_support=True` rule.

### GSE297234

No model reached predictive support. Improvements over persistence were small or negative, and the lower-tail criteria and permutation tests did not support robust prediction.

## 4. Transferability summary

| Model | Systems with positive mean | Systems with positive q05 | Systems with permutation support |
|---|---:|---:|---:|
| autoencoder | 2/3 | 2/3 | 1/3 |
| delta_t | 3/3 | 2/3 | 1/3 |
| markov | 3/3 | 1/3 | 0/3 |
| memory | 3/3 | 2/3 | 0/3 |
| pca | 2/3 | 2/3 | 1/3 |

The decisive observation is that predictive support is **system-specific rather than consistently reproduced across all three systems**.

## 5. Current scientific conclusion

### H-Z6a is not supported as a general explanation

The within-system positive result in GSE67462 is a counterexample to the strong statement that the current representation contains no predictive information even within a biological system.

Thus, the current evidence does **not** justify concluding that the representation is intrinsically non-predictive.

### H-Z6b is currently the preferred working interpretation

The combination of:

- negative cross-dataset Phase 0 results;
- positive within-system support in GSE67462;
- negative within-system results in GSE28688 and GSE297234;

supports the more specific interpretation that **predictive information can exist in the learned representation within a biological system but is not demonstrably transferable across the tested independent systems**.

This should be treated as a **working empirical conclusion**, not as a universal biological law.

## 6. What the result does not establish

The current Z6 evidence does **not** establish:

- a universal biological trajectory;
- causal state transitions;
- mechanistic validity of the latent dynamics;
- biological specificity of the predictive signal;
- transferability across species, platforms or experimental protocols;
- superiority of the learned representation in every biological context.

In particular, predictive performance alone is insufficient evidence of biological specificity. The unresolved Z4/H3 control problem remains relevant before a positive predictive result can be interpreted as specifically biological.

## 7. Scientific implications for the thesis

The current result shifts the Z6 question from:

> Can one learned state representation predict cellular trajectories universally?

 toward the more falsifiable question:

> Under which biological, experimental and observational conditions does a learned state representation contain information that supports prediction of future cellular state, and which components of that information transfer across independent systems?

This is consistent with the broader thesis objective of separating biological state from context-dependent and technical effects rather than treating every temporally predictive signal as biological state.

## 8. Phase 0.1 technical audit closure

During the initial Phase 0.1 run, `model_benchmark.py` emitted `RuntimeWarning: All-NaN slice encountered` from `training_statistics()` while computing the feature-wise median used for training-set imputation.

The warning source was identified as columns that were entirely non-finite within the training trajectories. The implementation was changed to detect all-NaN columns explicitly, compute the median only for columns with at least one finite value, and retain the existing downstream fallback of `0.0` for all-NaN columns.

The corrected benchmark was rerun with the complete Phase 0.1 protocol. The warning disappeared and **all reported numerical results remained identical** to the pre-fix run, including RMSE values, baseline improvements, q05 statistics, permutation p-values and predictive-support decisions for every system/model combination.

Therefore the warning was a preprocessing diagnostic artifact and did not affect the Phase 0.1 scientific conclusions. The predictive-support rule, model architecture, preprocessing protocol, baselines, validation branches and dataset inclusion were not relaxed or changed.

Phase 0.1 is therefore considered **technically closed and stable**, while its scientific interpretation remains subject to the limitations documented above.

## 9. Next required audit before extending the model

With the NaN warning resolved, the next step should **not** be architecture expansion or threshold relaxation. The scientifically relevant follow-up is to investigate representation/context dependence and the unresolved biological-specificity problem (Z4/H3), especially whether the positive GSE67462 prediction reflects transferable biological state information or system-specific temporal/context structure.

Any follow-up must preserve the frozen Phase 0/0.1 results as a reference point and must not retroactively redefine predictive support.

## 10. Provenance

Primary machine-readable outputs from Phase 0.1 are stored under:

`results/Dynamics/z6_within_system_audit/`

including branch inventory, fold-level metrics, model summaries, support decisions and the protocol snapshot.

The methodological specification is recorded in:

`dynamics/Z6_WITHIN_SYSTEM_TRANSFERABILITY_PROTOCOL.md`

The primary Z6 protocol is recorded in:

`dynamics/Z6_PREDICTIVE_TRANSITION_PROTOCOL.md`

This document is the human-readable scientific interpretation of those versioned experimental artifacts.
