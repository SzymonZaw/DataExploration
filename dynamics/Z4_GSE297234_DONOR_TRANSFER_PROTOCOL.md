# Z4 GSE297234 Donor-Stratified Transfer Protocol

## Purpose

Evaluate whether the frozen GSE67462 temporal modules retain temporal direction independently in the two human donors represented in GSE297234:

- GM00731 — aged donor
- GM23815 — young donor

This is an external transfer/replication diagnostic. It does **not** by itself satisfy the frozen Z4 acceptance rule.

## Frozen inputs

- GSE67462 temporal-module assignments
- frozen human↔mouse one-to-one orthology mapping
- GSE297234 sample-level frozen module scores

No module discovery, feature selection, scaling, or parameter tuning is performed on GSE297234 in this audit.

## Analysis

For each frozen module and donor:

1. order the four samples by inferred time (D0, D3, D7, D10),
2. calculate Spearman correlation between time and frozen module score,
3. classify direction as increasing/decreasing/unresolved,
4. compare the direction between aged and young donors.

The four timepoints are sparse, so this is descriptive transfer evidence rather than a causal or predictive test.

## Interpretation

Strong donor agreement can support the bounded statement that a frozen temporal representation transfers reproducibly across the two independent human donors in GSE297234.

It must not be interpreted as:

- proof of biological specificity,
- proof of causal mechanism,
- proof of universal transferability,
- an independent endpoint,
- or formal Z4 support.

Z4 remains unresolved until an independent biological endpoint permits target-vs-nontarget specificity and context-robustness testing under the frozen rule.
