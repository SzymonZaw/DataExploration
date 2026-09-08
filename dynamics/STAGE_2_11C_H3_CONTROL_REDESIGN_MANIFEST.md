# Stage 2.11C — H3 Control Redesign Manifest

Status: **prospectively locked — final H3 controls selected before similarity analysis**

This manifest supersedes the *primary analytical role* of the historical H3 controls. The historical manifest remains immutable and is retained for provenance.

## 1. Locked purpose

Select at least two biologically unrelated, time-resolved processes to test whether the state-representation pipeline produces target-like trajectory agreement outside Yamanaka reprogramming.

The controls are negative controls for **generic temporal/representation structure**, not controls for shared biology.

## 2. Mandatory selection criteria

A candidate can enter the final H3 control set only if all criteria are satisfied before similarity analysis:

- genuine time-resolved biological process;
- no known direct role as a reprogramming intervention, intermediate, pluripotency induction, or aging-reversal paradigm;
- not a reprogramming experiment;
- not selected using observed similarity to Yamanaka trajectories;
- sufficient temporal points and biological replication for the locked null procedure;
- compatible transcriptomic measurement;
- sufficient metadata for time and biological replicate/block definition;
- processed data can be passed through the same representation/preprocessing contract;
- biological relationship to reprogramming can be documented as absent or sufficiently remote for the H3 purpose.

## 3. Orthogonality requirement

The final set must contain **at least two independent process classes**.

The two controls must not both be replicative senescence/cellular lifespan datasets. This explicitly prevents the previous GSE175533 + GSE179848 pair from being reused as the primary H3 negative-control set.

## 4. Candidate-selection audit

The candidate audit is committed separately in `STAGE_2_11C_H3_CANDIDATE_AUDIT.md`. Selection was based on metadata/design properties only; no H3 similarity, trajectory correlation, classifier performance, or representation score was inspected for candidate selection.

## 5. Prospective statistical defaults

These values are binding for the final H3 run:

- `delta_H3 = 0.10` absolute trajectory-agreement units;
- one-sided empirical null criterion: `p < 0.05`;
- BH correction: `q < 0.05` for predeclared feature-level inferential families;
- biological replicate/dataset block is the exchangeability unit;
- no cell-level permutation when the representation is pseudobulk/sample-level;
- if attainable p-value resolution is too coarse for the declared criterion, return `H3_UNRESOLVED`.

The 0.10 margin is a prospective operational threshold, not a biological law.

## 6. Final locked controls

The following controls are locked **before any H3 similarity analysis**:

| Role | Accession | Process | Expected local file | Locked input SHA-256 |
|---|---|---|---|---|
| H3-Control-1 | GSE3945 | acute serum-response / wound-healing program in human fibroblasts | `Data/GSE3945_series_matrix.txt.gz` | `a79ec8c7b5bb88742a43bba7495856dd6a3d3e82b802588ea58d9773ef0598eb` |
| H3-Control-2a | GSE129486 | acute inflammatory cytokine response in human synovial fibroblasts — gene TPM matrix | `Data/GSE129486_rnaseq-data-1_gene-tpm.tsv.gz` | `6e3d7860f4f38d95830a15b8dd570d58f9226170df6002343e432a05c96fcbf0` |
| H3-Control-2b | GSE129486 | metadata paired with the gene TPM matrix | `Data/GSE129486_rnaseq-data-1_metadata.tsv.gz` | `78a60e353461ab819672e478a9f38b20822fdfc493a1773677ddef3eb167ff04` |

These hashes were recorded from the exact local resources acquired by `dynamics/download_h3_controls.py`, before any H3 similarity analysis.

## 7. Historical controls — retained, not primary

The earlier selection remains documented separately:

- GSE175533 — replicative cellular senescence / cellular lifespan;
- GSE179848 — replicative cellular lifespan / aging in primary human fibroblasts.

These datasets must not be retroactively treated as if they had satisfied this redesigned orthogonality criterion.

## 8. Execution lock

Before analysis, persist SHA-256 hashes of the exact local input files and the exact mapping artifact/version used by the representation pipeline. Also persist a machine-readable technical audit with matrix orientation, time labels, replicate/block structure, identifier namespace, mapping/collision statistics, retained genes, shared representation genes, zero/missing handling, preprocessing parameters, and representation/network versions.

If either final control fails technical validation, do not replace it after inspecting H3 similarity. Return `H3_UNRESOLVED` and create a new prospective manifest if a replacement design is scientifically justified.

## 9. Interpretation lock

- `H3_SUPPORTED`: unrelated temporal processes show target-like agreement under the locked criterion; generic temporal/representation structure remains plausible.
- `H3_NOT_SUPPORTED`: both unrelated controls are materially separated from the target under the locked criterion; H3 is weakened, but process specificity is not thereby proven.
- `H3_UNRESOLVED`: insufficient or ambiguous evidence; no mechanistic or predictive escalation.

## 10. Provenance

This manifest is a new prospective gate. It does not modify the historical Stage 2.11C result, the historical audit, or the historical GSE175533/GSE179848 selection record.
