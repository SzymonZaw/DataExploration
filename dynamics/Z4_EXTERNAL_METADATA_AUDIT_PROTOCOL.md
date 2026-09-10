# Z4 External Metadata Audit Protocol

## Purpose

Prepare independent GEO candidates for the frozen Z4 commensurability audit **without evaluating the frozen transfer score**.

The audit extracts sample-level provenance from GEO SOFT metadata and checks the pre-specified design gates needed before any human-to-mouse feature mapping or frozen transfer is attempted.

## Candidates

- `GSE297234`: primary homologous external target — human fibroblast partial reprogramming with Sendai OSKM, days 0/3/7/10, two donor lines.
- `GSE28688`: predefined specificity/context challenge — human HFF1 fibroblasts, OSKM transduction at 24/48/72 h plus iPSC/hESC endpoints.

## No-look rule

This script must not consume frozen transfer scores, learned embeddings, validation predictions, or downstream enrichment results. Candidate eligibility is determined from provenance/design metadata only.

## Required outputs

For each candidate:

- `metadata.csv`: one row per GEO sample with accession, title, time, donor, intervention, platform, batch and raw characteristics.
- `audit.json`: design-level commensurability summary.
- `metadata_raw.txt`: cached GEO SOFT source for provenance.

## Gate logic

### C1 — temporal commensurability

At least 3 ordered biological time points must be identifiable for a transition candidate.

### C2 — cell-state commensurability

The source metadata must support a pre-specified relationship to fibroblast-state change/reprogramming. This is a human biological target, not a literal mouse trajectory replication.

### C3 — intervention commensurability

OSKM intervention must be explicitly identifiable. Sendai OSKM is accepted for GSE297234 as a homologous intervention; viral OSKM transduction is accepted for GSE28688 as a context challenge.

### C4 — feature-space commensurability

This metadata audit does not manufacture orthology mappings. A separate deterministic human→mouse mapping table must be supplied before frozen transfer. Mapping must be frozen before looking at transfer results.

### C5 — endpoint independence

The biological endpoint must be independently defined. It must not be constructed from the transferred feature score.

### C6 — nuisance observability

Donor/age/batch/platform/time information must be observable where available. Missing nuisance metadata does not prove confounding, but it blocks a clean specificity claim.

## Important boundary

Passing this audit means only **eligible for frozen transfer preparation**. It does not establish Z4 specificity, and it does not evaluate the predictive representation.
