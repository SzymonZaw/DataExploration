# Z6 GSE67462 signal attribution protocol

## Purpose

This audit follows the positive Z6 Phase 0.1 result in GSE67462 and the subsequent context-specificity audit. Its purpose is to decompose the reproducible transition signal into:

1. common biological transition structure;
2. branch/context-specific structure;
3. technical or nuisance structure;
4. endpoint-only structure.

The audit is diagnostic. It does not modify the frozen Z6 predictive-support decision.

## Dataset context

GSE67462 is a mouse secondary-reprogramming time series with bulk expression profiling and two samples per time point. GEO describes a time-series design spanning day 0 through day 18, with OSKM induction and subsequent collection of reprogrammed/iPS cells. The current Z6 common-space mapping uses the eight shared timed points represented in the benchmark branches. The original GEO study also reports complementary Oct4 binding and histone-modification profiling in the same biological program.

## Attribution logic

For each gene and each branch, define the centered trajectory relative to day 0. Estimate four components:

### A. Common transition component

The branch mean trajectory:

`C_g(t) = mean(X_g,branch1(t), X_g,branch2(t))`

A signal is considered common when both branches contribute concordantly and the common trajectory explains a substantial fraction of the branch-specific trajectory variance.

### B. Branch/context component

The branch contrast:

`B_g(t) = X_g,branch1(t) - X_g,branch2(t)`

Persistent or time-structured branch contrasts are treated as context dependence rather than automatically as biological state.

### C. Technical/nuisance component

Technical attribution is evaluated using the available sample-level nuisance proxies and branch structure. No biological interpretation is assigned solely from temporal prediction. If a candidate signal is strongly associated with library/measurement intensity, missingness, or branch identity after accounting for time, it is flagged as nuisance/context-sensitive.

### D. Endpoint-only component

Endpoint-only signals are identified by comparing full-trajectory reproducibility against endpoint effect concordance. A high day-0-to-final correlation without intermediate-time consistency is not treated as evidence for a dynamic state trajectory.

## Required outputs

The executable must produce:

- gene-level common-transition scores;
- branch-contrast scores;
- endpoint-only scores;
- nuisance-association diagnostics where available;
- pathway/regulatory aggregation for PROGENy and DoRothEA where existing common-space outputs permit it;
- permutation nulls preserving time labels;
- a system-level attribution table.

## Frozen rules

- Do not change the Phase 0 or Phase 0.1 predictive-support rule.
- Do not tune model architecture.
- Do not remove GSE67462 because it is positive.
- Do not use the positive prediction result itself as evidence of biological specificity.
- Do not call branch reproducibility causal.
- Preserve all intermediate outputs for provenance.

## Interpretation hierarchy

The preferred interpretation hierarchy is:

1. **COMMON_DYNAMIC_SIGNAL** — reproducible across branches across intermediate timepoints;
2. **COMMON_ENDPOINT_SIGNAL** — reproducible mainly at the endpoint;
3. **BRANCH_CONTEXT_SIGNAL** — predictive/reproducible structure concentrated in branch differences;
4. **NUISANCE_CANDIDATE** — strongly associated with technical/sample-intensity proxies;
5. **UNRESOLVED** — insufficient evidence to attribute the signal.

The audit can support biological interpretation only when a signal survives the dynamic, branch, and nuisance checks. Even then, the result is associative rather than causal.

## Scientific decision

The central question is not whether GSE67462 is predictive; that has already been established for AE and PCA under the frozen Phase 0.1 criterion. The question is whether the predictive information is carried predominantly by a reproducible dynamic molecular component rather than by branch identity, endpoint separation, or technical/context structure.
