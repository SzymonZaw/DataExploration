# Stage 2.11C — H3 Control Redesign Protocol

## Status

**Prospective redesign. No H3 results may be analyzed under this protocol until the control manifest and decision thresholds are committed.**

This document does not alter the historical Stage 2.11C result or its audit. It replaces the scientific design of the *future H3 test* because the previously selected controls are not sufficiently orthogonal to the reprogramming hypothesis.

## 1. Scientific question

The unresolved question is:

> **Does the state representation `S(t)` contain information specific to the target process, rather than merely stable temporal structure that appears across unrelated biological processes?**

H3 concerns the latter possibility. It is not a test of whether two biological processes share any biology.

## 2. Why the previous H3 controls are retired from primary use

The previously locked controls GSE175533 and GSE179848 are retained as historical design records, but they are **not valid primary negative controls for the redesigned H3 gate**.

Replicative senescence / cellular lifespan is biologically related to cell-state regulation and shares regulatory axes with reprogramming. Therefore a high trajectory agreement would be ambiguous between:

1. generic representation/aggregation structure; and
2. genuine overlap between two biologically related processes.

Such ambiguity prevents the control from cleanly falsifying H3.

The historical manifest must remain immutable for provenance. No historical result may be relabeled as if these controls had been tested under the redesigned protocol.

## 3. H3 control definition

A primary H3 negative control must satisfy all of the following **before analysis**:

1. a genuine time-resolved biological process;
2. no known direct role as a reprogramming intervention, reprogramming intermediate, pluripotency induction, or aging-reversal paradigm;
3. no selection based on observed similarity to or dissimilarity from the Yamanaka result;
4. sufficient time points and biological replication for the declared null procedure;
5. transcriptomic data compatible with the locked representation pipeline;
6. metadata sufficient to define the temporal unit and replicate/block structure;
7. biological interpretation that does not depend on treating the process as a reprogramming-like state transition.

Examples of candidate classes include acute stimulation-response trajectories, wound/repair responses, or other temporally resolved cellular responses with no mechanistic dependence on reprogramming. Examples are candidate classes, not approved datasets.

## 4. Required number and orthogonality of controls

Use **at least two independent unrelated temporal processes**.

Prefer two biologically distinct process classes rather than two variants of the same process. At minimum, the selected pair should not both belong to cellular aging/senescence.

The controls are not required to be biologically identical to the Yamanaka system. Their purpose is to test whether the representation pipeline generates high temporal agreement outside the target process.

If fewer than two valid controls can be locked before analysis, the H3 gate is **not executable** and the primary status remains `SIGNAL_UNRESOLVED`.

## 5. Selection procedure and anti-cherry-picking rule

Dataset selection must be based only on criteria observable **before computing H3 similarity**.

For every candidate considered, record:

- accession;
- organism and cell system;
- process class;
- number of time points;
- number of biological replicates;
- measurement type;
- processed-resource availability;
- identifier namespace;
- expected preprocessing compatibility;
- explicit biological relationship to reprogramming;
- inclusion/exclusion decision and reason.

Do not inspect trajectory similarity before the final selection.

A candidate rejected because it appears too similar to Yamanaka is evidence of selection bias and invalidates that selection round; create a new prospective manifest instead of silently replacing the candidate.

## 6. Statistical design

The H3 comparison must use the same representation family and preprocessing rules as the target analysis.

### 6.1 Primary quantity

For each unrelated process `u`, compute the cross-dataset trajectory-agreement statistic using the same feature identities and scoring procedure used for the Yamanaka comparison.

Report:

- observed target-process agreement;
- observed agreement for each unrelated process;
- difference from the target-process value;
- empirical null distribution;
- exact/Monte-Carlo attainable p-value resolution where applicable.

The raw correlation value alone is never a decision criterion.

### 6.2 Null model

The null must preserve the internal temporal structure of every dataset while breaking the cross-dataset identity correspondence being tested.

Where samples/replicates are the unit of exchangeability, permutations must respect biological replicate/block structure. Cells must never be treated as independent exchangeable units when the analysis is pseudobulk or sample-level.

If the number of temporal points makes an exact permutation test too coarse to support the declared decision, the protocol must say so and return `SIGNAL_UNRESOLVED`; do not compensate by changing alpha after seeing the result.

### 6.3 Multiple testing

If candidate-level inferential p-values are produced, apply Benjamini–Hochberg correction across the **predeclared inferential family** before feature-level claims are made.

The number of controls or the deterministic comparison count is not itself an alpha-level statistical test.

## 7. Prespecified decision threshold

H3 must not use an ad hoc rule such as “high correlation means artifact.”

Before analysis, define a **relative-control margin** `δ_H3` and a null-tail criterion.

Primary H3-supporting criterion:

> The target-process agreement is **not materially higher** than the unrelated-control agreement distribution under the locked null, i.e. the target-to-control separation is at or below `δ_H3` and the empirical null does not support a target-specific excess.

Primary evidence against H3:

> The target-process agreement exceeds the unrelated-control distribution by **more than `δ_H3`**, with the separation supported by the locked permutation/block-bootstrap procedure.

Recommended prospective default for the new manifest: `δ_H3 = 0.10` absolute trajectory-agreement units, with a one-sided empirical null criterion of `p < 0.05` and BH `q < 0.05` where feature-level inference is performed.

These values are **prospective defaults only**. They become binding only when copied into and committed with the final control manifest before analysis.

## 8. Decision states

Use exactly one primary state:

### `H3_SUPPORTED`

At least one independent unrelated temporal process produces target-like agreement within the prespecified `δ_H3` margin, and the result is not explained by an execution/data-quality failure.

Interpretation: generic temporal/representation structure is a plausible contributor to the observed agreement.

### `H3_NOT_SUPPORTED`

The target process is materially separated from **both** unrelated controls by more than `δ_H3` under the locked null procedure.

Interpretation: H3 is weakened, but process specificity is not automatically proven; matched nuisance/context attribution remains necessary.

### `H3_UNRESOLVED`

The controls are insufficient, statistically underpowered, technically incompatible, or yield ambiguous evidence under the locked criteria.

Interpretation: no claim about process specificity is permitted.

## 9. What H3 can and cannot establish

A successful H3 falsification does **not** prove H1.

Specifically:

- `H3_NOT_SUPPORTED` does not establish causality;
- low similarity to unrelated processes does not identify a molecular mechanism;
- high similarity to a biologically related process is not automatically an aggregation artifact;
- the H3 test cannot identify Sendai, Myc, interferon/STAT, or any other nuisance source without an isolating control.

The scientific logic is:

```text
H3 supported      -> generic temporal structure remains plausible
H3 not supported  -> H3 weakened; specificity still requires H2/process discrimination
H3 unresolved     -> no escalation
```

## 10. Technical data-quality gate

Before any H3 statistic is computed, persist a machine-readable audit containing:

- file SHA-256;
- dataset accession;
- sample count;
- time-point count and exact labels;
- replicate/block count;
- matrix orientation;
- gene identifier namespace;
- mapping rate and collision count where mapping occurs;
- number of genes retained;
- number of shared representation genes;
- missing-value/zero handling;
- preprocessing parameters;
- representation/network versions.

For the legacy GSE3945 microarray resource, missing expression cells are permitted at the technical-validation stage. The locked structural criterion is that at least **90% of non-missing expression cells in the first 200 expression rows are numeric**, while missingness is reported explicitly. This is a data-format gate, not a biological or inferential decision threshold.

A parser repair is not considered sufficient evidence of data validity. The reconstructed matrix structure must be independently checked against the source metadata.

## 11. Historical boundary

The following remain unchanged:

- historical Stage 2.11C result: 28 Tier A → 1 Tier B;
- historical audit conclusion: `SPECIFICITY_UNSUPPORTED`;
- historical GSE175533/GSE179848 selection record;
- historical signal-attribution status: `SIGNAL_UNRESOLVED`.

This redesign is a **new prospective gate**, not a retroactive correction of historical results.

## 12. Stopping rule

Do not keep adding unrelated controls until the result becomes favorable.

The redesigned H3 analysis is limited to the precommitted control set. If the result is `H3_UNRESOLVED`, the correct response is to document the limitation and decide whether a genuinely new prospective protocol is warranted.

Two consecutive prospective designs that fail their locked discrimination criteria trigger the existing methodological stopping rule rather than another post-hoc control swap.
