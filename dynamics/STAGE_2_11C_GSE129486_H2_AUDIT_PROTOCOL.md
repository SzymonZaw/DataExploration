# Stage 2.11C — GSE129486 H2 Independence Audit

## Status

Prospective design audit. No Yamanaka-vs-control similarity is computed here.

## Purpose

Determine whether GSE129486 can legitimately serve as an H3 unrelated temporal control, and separately quantify whether its inflammatory design overlaps the predeclared H2 nuisance hypothesis.

The key distinction is:

- **H3:** generic temporal/representation structure;
- **H2:** nuisance biology, including the predeclared inflammatory/JAK-STAT/interferon-STAT-like axes.

A control that is temporally excellent but shares the suspected nuisance axis is not a clean discriminator between H2 and H3.

## Required structural audit

The audit must report:

1. exact input SHA-256 hashes;
2. sample IDs and duplicate IDs;
3. exact time labels and samples per time point;
4. `cell_line × stimulation` biological blocks;
5. number of time points represented in each block;
6. whether blocks are complete across the common time grid;
7. whether replication is biological/sample-level rather than cell-level pseudoreplication;
8. identifier namespace and mapping integrity;
9. treatment labels and process class.

## H2 independence rule

GSE129486 is an acute inflammatory cytokine-response experiment. Because the predeclared H2 nuisance hypothesis includes inflammatory/JAK-STAT/interferon-STAT-like programs, GSE129486 receives an **H2 orthogonality risk flag** by design.

This flag is not evidence that those pathways are activated in GSE129486. Actual activity must be quantified separately from H3 similarity.

## Decision rule before H3 similarity

- `STRUCTURALLY_EXECUTABLE`: temporal and block structure is adequate.
- `H2_ORTHOGONALITY_LIMITED`: the process class overlaps the predeclared nuisance hypothesis; it cannot be the sole H3 control.
- `H3_CONTROL_ELIGIBLE_WITH_SECOND_ORTHOGONAL_CONTROL`: may remain in the prospective control set only if a second, biologically distinct temporal control is also locked and the nuisance-axis audit does not make the interpretation circular.
- `H3_UNRESOLVED`: if structural or independence requirements fail.

## Anti-circularity rule

Do not inspect Yamanaka-vs-GSE129486 trajectory similarity and then decide whether GSE129486 is sufficiently orthogonal. Control eligibility is determined from metadata/design and predeclared nuisance biology before H3 similarity.

Likewise, actual H2 pathway activity is an independent attribution analysis; it must not be used to relabel an H3 result after observing the H3 statistic.

## Current expected interpretation

GSE129486 has strong temporal replication and an explicit `cell_line × stimulation × time` structure, but it should be treated as **limited-orthogonality rather than a clean H3 negative control** until the nuisance-axis audit is completed.
