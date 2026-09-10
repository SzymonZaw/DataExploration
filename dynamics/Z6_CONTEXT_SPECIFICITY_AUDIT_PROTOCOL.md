# Z6 context-specificity audit protocol

## Purpose

The Phase 0.1 result establishes predictive support only for autoencoder and PCA in GSE67462. Before treating that result as evidence for a biologically meaningful state representation, this audit tests whether the two temporal branches of each system reproduce the same underlying molecular transition.

The audit is deliberately orthogonal to the predictive-support test. It does not change models, thresholds, train/test splits or the Phase 0/0.1 conclusions.

## Primary question

Is the positive GSE67462 forecasting result accompanied by a reproducible molecular transition across its two independent branches, or is it compatible with branch-specific temporal/context structure?

## Data

Use the same Stage 2.7 common gene-space matrix and the same branch definitions as Z6 Phase 0.1:

- GSE28688: branches `a` and `b`;
- GSE67462: replicate branches `1` and `2`;
- GSE297234: donor branches `aged` and `young`.

GSE297234 donors remain biologically independent donor holdouts, not technical replicates.

## Locked measurements

For each system, align the two branch trajectories on their shared native timepoints and calculate:

1. **Feature-wise temporal concordance**: Spearman correlation across time for every common gene, then median absolute correlation.
2. **Directional concordance**: fraction of genes whose endpoint change has the same sign in both branches.
3. **Endpoint-effect concordance**: Pearson correlation of branch-specific endpoint-minus-baseline gene effects.
4. **Per-timepoint molecular concordance**: Spearman correlation across genes between branches at each shared timepoint; report median across timepoints.
5. **Permutation diagnostic**: shuffle gene identities in one branch and recompute endpoint-effect correlation to provide a distributional reference.

No threshold is used to relabel predictive support. These are descriptive/context-dependence diagnostics.

## Interpretation

The audit distinguishes:

- **REPRODUCIBLE_TRANSITION**: high branch concordance across multiple metrics;
- **MIXED_CONTEXT_DEPENDENCE**: predictive support occurs but branch concordance is heterogeneous;
- **BRANCH_SPECIFIC_STRUCTURE**: predictive support occurs alongside weak branch concordance.

These labels are not claims of causality or biological specificity. In particular, high concordance can still represent a generic temporal response, stress, proliferation or other nuisance process.

## Decision rule

The positive GSE67462 prediction remains a valid Z6 predictive result regardless of this audit. The audit only determines how strongly it can be interpreted as a reproducible state trajectory.

If GSE67462 is reproducible while GSE28688/GSE297234 are not, the preferred interpretation is **system-specific reproducible dynamics**, not universal biological state transfer.

If GSE67462 is weakly concordant across its two branches, the positive forecast is treated as **context-dependent predictive structure** and should not be used as evidence of a transferable biological state without further Z4/H3 controls.

## Non-goals

This audit does not:

- change the frozen predictive-support rule;
- tune the forecasting models;
- remove any dataset;
- use predictive performance to choose favorable features;
- call generic temporal concordance biological specificity;
- replace the prospective H3 control program.

## Output

`results/Dynamics/z6_context_specificity_audit/`

Expected files:

- `01_branch_temporal_concordance.csv`
- `02_branch_directional_concordance.csv`
- `03_branch_endpoint_correlation.csv`
- `04_branch_timepoint_correlation.csv`
- `05_permutation_endpoint_null.csv`
- `06_system_interpretation.csv`
- `07_protocol.json`
