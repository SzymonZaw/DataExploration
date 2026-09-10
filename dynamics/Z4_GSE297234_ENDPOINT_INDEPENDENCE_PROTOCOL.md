# Z4 GSE297234 endpoint-independence audit

## Purpose

Determine whether variables available in the independent GSE297234 Seurat object can serve as an **independent biological endpoint** for the frozen Z4 external-specificity test.

This audit is deliberately conservative. A variable being biologically plausible, statistically associated with the transferred representation, or present in the published/processed object is **not sufficient** to establish endpoint independence.

## Frozen Z4 requirement

The endpoint must be external to the representation being tested and must support a target-vs-non-target specificity question. The audit therefore distinguishes:

1. **Raw experimental metadata** — donor, sample identity, intervention, time, platform.
2. **Derived annotations** — cell-state scores, pathway scores, classifier labels, ordered clusters, or other variables plausibly computed from the same expression data.
3. **Independent biological endpoint** — a measurement generated independently of the expression representation and with explicit provenance that makes it suitable as an outcome for the Z4 test.

Only category 3 can satisfy endpoint independence.

## GSE297234 candidate fields

The current RDS contains fields including:

- `orig.ident`
- `age_group`, `age_ident`
- `cell_state`
- `PartialReprog1`, `NonReprog1`
- `EarlyPluripotency1`, `Pluripotency1`
- `ordered_clusters`
- `Fibroblast1`
- `PI16_fibtype1`, `LRRC15_fibtype1`, `COL3A1_fibtype1`
- `HALLMARK_EMT1`, `HALLMARK_TGFB1`

The audit does **not** assume that any of these are independent endpoints. In particular, expression-derived scores and state labels are treated as derived until independent provenance is demonstrated.

## Required checks

### E1 — Sample resolution

All eight GSE297234 samples must be represented exactly once in the transfer audit:

`GM00731_D0`, `GM00731_D3`, `GM00731_D7`, `GM00731_D10`, `GM23815_D0`, `GM23815_D3`, `GM23815_D7`, `GM23815_D10`.

### E2 — Candidate provenance

For each candidate endpoint, record whether its provenance is:

- `raw_metadata`
- `derived_from_same_expression_object`
- `unknown`
- `independent_external_measurement`

A candidate with `derived_from_same_expression_object` or `unknown` provenance cannot establish endpoint independence.

### E3 — Temporal structure

Determine whether the candidate endpoint is available at sample level across the four intervention timepoints and both donors.

This is a descriptive check only. A temporal trend does not establish independence.

### E4 — Target/non-target contrast

Determine whether the candidate distinguishes the intended biological transition from a non-target state/process. Merely increasing over OSKM time is insufficient.

### E5 — Leakage exclusion

Check whether the candidate is a direct transformation of expression features, pathway scores, clusters, or labels generated from the same Seurat object. If so, it is not an independent endpoint for Z4.

### E6 — Decision

The frozen decision remains `Z4_UNRESOLVED` unless at least one endpoint passes E2, E3, E4 and E5 with documented independent provenance.

No threshold, module definition, mapping rule, or representation parameter may be changed to obtain a positive result.

## Interpretation policy

A result such as `TRANSFER_ASSOCIATED_BUT_ENDPOINT_NONINDEPENDENT` is a useful external-transfer observation, but it is **not** Z4 support.

A result of `ENDPOINT_INDEPENDENCE_UNESTABLISHED` is not a failed biological hypothesis. It means the current dataset/object does not provide the independent outcome needed to test the hypothesis under the frozen Z4 rule.

## Expected output

The audit writes:

- candidate endpoint table,
- sample-level availability table,
- provenance/leakage assessment,
- transfer-association diagnostics where already computed,
- final endpoint-independence decision.

The frozen Z4 acceptance rule is unchanged.
