# Stage 2.11C — Audit and Prospective Revision Protocol

## Purpose

This document records the methodological audit of the existing Stage 2.11C result. It is intentionally separate from the original protocol so that the historical protocol remains immutable and its provenance can be inspected.

The audit does **not** rerun the full Stage 2.11C pipeline.

## Historical result under audit

The existing local run reported:

- Tier A: 28 candidates;
- Tier B: 1 candidate;
- survival: 1/28 = 3.57%;
- historical decision label: `CONFOUNDED_UNSUPPORTED`.

The audit replaces the interpretation of that label with:

- `SPECIFICITY_UNSUPPORTED`;
- `confound_source: unresolved`.

The historical output is not modified in-place.

## Finding 1 — context control does not identify a unique confounder

GSE304042 is a heterologous biological/context comparison. It differs from the original GSE297233 system in multiple dimensions, including cell type, perturbation structure (OSK versus OSKM), experimental design and delivery/context.

Therefore the comparison can establish lack of cross-system generalization, but cannot identify Sendai, interferon/STAT, Myc, cell type, or any single factor as the causal confounder.

Required language:

> The candidate signal does not generalize across the tested biological contexts; the source of the context dependence remains unresolved.

Forbidden upgrade:

> Sendai is the identified confounder.

## Finding 2 — gene-ID mapping is technically resolved but not independently locked

The GSE297233 diagnostic established that the raw matrix uses Ensembl gene identifiers and that mapping to HGNC symbols creates substantial network overlap.

Recorded values include:

- 62,703 input gene IDs;
- 45,906 mapped IDs;
- 44,554 unique mapped symbols;
- 1,352 mapped-row collisions (`45,906 - 44,554`);
- 1,281 PROGENy targets shared after mapping;
- 9,143 DoRothEA targets shared after mapping.

The current implementation averages duplicate mapped HGNC symbols. This is a defined transformation, but the exact external MyGene mapping is not persisted as a versioned mapping table. Consequently the mapping is considered **partially validated**, not fully reproducibility-locked.

Prospective requirement:

1. persist the exact Ensembl → HGNC mapping table used for an analysis;
2. record mapping-service/annotation version where available;
3. record counts of one-to-one mappings, many-to-one collisions, unmapped IDs and post-aggregation symbols;
4. retain the mapping artifact together with the result manifest.

## Finding 3 — the original protocol and implementation disagree on Null 1

The original protocol specified an empirical permutation p-value `< 0.01`. With four time points and exclusion of the observed ordering, there are only 23 null permutations and the minimum +1-corrected empirical p-value is `1/24 ≈ 0.0417`.

The implementation instead uses `p <= 0.05` together with the exact 23-permutation null.

This mismatch means the historical run cannot be described as strictly compliant with the written Null-1 decision rule.

Prospective correction:

- retain the exact 23-permutation reference for four time points;
- use an attainable threshold consistent with that resolution, explicitly `p <= 0.05` if that is the intended rule;
- record the exact attainable p-value grid in machine-readable output;
- never claim p < 0.01 resolution from this four-time-point permutation space.

The corrected rule must be committed **before** any prospective rerun whose result is intended as confirmatory evidence.

## Finding 4 — multiple testing applies to candidate-level inferential tests, not to the deterministic Tier-B survival count

Tier B is currently a deterministic conjunction of:

- concordant direction;
- non-zero activity in both contexts;
- top-25% rank filter;
- exclusion of the predefined delivery-confound panel.

There is no Tier-B alpha-level p-value. Therefore Bonferroni/FDR should not be mechanically applied to the number `1/28`.

However, Null-1 generates candidate-level empirical p-values. Those p-values form an inferential family and require multiplicity control where the candidate family is large enough to support such inference.

The lightweight audit therefore calculates Benjamini–Hochberg q-values from the **existing saved Null-1 p-values** without recomputing the pipeline.

Important limitation: with only 23 null permutations, p-values are discrete and coarse. BH correction can therefore eliminate all apparent significance even when several raw p-values equal the minimum attainable value. This is a valid reason to report the result as exploratory rather than evidence of mechanism.

## Finding 5 — preregistration provenance

The original Stage 2.11C protocol was committed in Git before the later Ensembl-normalization commits. The protocol commit is:

`f90972f5f312b0e4ba39ec64693ccd63b7807c79`

with commit time `2026-09-07T19:40:25Z`.

The later Ensembl normalization commit is:

`331adb99acc233297195ec721f2f18d89b331da4`

with commit time `2026-09-08T10:35:58Z`.

The audit runner also compares the protocol commit with the local saved-result timestamp when available. This provides provenance evidence that the decision thresholds were recorded before the observed local result, but a filesystem timestamp is not a cryptographic execution log.

If the local output timestamp is unavailable or inconsistent, the result must be classified as lacking independently verifiable execution-time provenance.

## Prospective decision policy

Until the protocol/implementation mismatch and mapping reproducibility issue are resolved:

- the historical Stage 2.11C result is **exploratory**;
- `1/28` is not treated as evidence for a true biological regulator;
- the heterologous context result is interpreted as non-generalization, not identification of a confounder;
- no Tier-B candidate is promoted to a mechanistic claim.

After the prospective protocol is locked, a future validation dataset/control should be used to test the revised decision rule without changing thresholds in response to its result.

## Scientific guardrail

The statement

> Reproducibility of a state representation does not imply reproducibility of a mechanism.

is retained only as a **hypothesis to be tested**, not as a conclusion established by the current 28-to-1 result.

The immediate scientific objective is therefore:

> determine whether candidate temporal signals survive technically valid mapping, correctly specified inferential tests, multiplicity control and an identifiable biological control design.
