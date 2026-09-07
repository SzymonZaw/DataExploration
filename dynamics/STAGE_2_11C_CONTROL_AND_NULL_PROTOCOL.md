# Stage 2.11C — Control and Null Protocol

## Status

**Protocol-first document.** This file defines the decision rules before implementation and before inspecting the Stage 2.11C results.

The purpose is to determine whether the temporal pathway/TF signals observed in the Yamanaka proof-of-feasibility are (a) temporally ordered, (b) reproducible, and (c) sufficiently specific to justify progression toward perturbation-conditioned dynamical modeling.

This protocol is not a causal mechanism-discovery protocol. It is a falsification gate between the current proof-of-feasibility and the next modelling stage.

---

## 1. Scientific question

Do the temporal pathway/TF activity patterns observed during OSKM reprogramming contain reproducible biological information that cannot be explained adequately by generic monotonic change, arbitrary temporal ordering, or the experimental delivery/system context?

The central downstream question remains:

> Can a biologically interpretable state representation preserve information needed to predict future cellular state under explicit perturbation?

Stage 2.11C does **not** attempt to answer that question directly. It determines whether the current Yamanaka signal is sufficiently credible to justify that experiment.

---

## 2. Starting evidence

The current Yamanaka trajectory POC uses two independent human fibroblast objects from GSE297234 with days 0/3/7/10 and sample-level pseudobulk aggregation.

The observed temporal-distance-profile correlation is approximately **0.921**.

This is treated as an **observation requiring validation**, not as evidence of statistical significance or mechanism.

Stage 2.11 mechanism discovery also found highly reproducible pathway/TF activity trajectories. These are referred to here as **candidate temporal regulators**, not mechanisms.

---

## 3. Terminology and separation of evidence types

Three evidence classes must not be conflated.

### 3.1 Statistical null

A constructed reference distribution asks whether the observed statistic can arise when the relevant structure is removed while preserving as much of the data structure as possible.

### 3.2 Shape control

A control asks whether the result is stronger than what would be expected from generic monotonic temporal processes. This addresses the possibility that two curves agree simply because both change monotonically with time.

### 3.3 Biological/context control

A heterogeneous biological comparison asks whether candidate activity patterns persist outside the exact original experimental context. GSE304042 is **not** a clean Sendai-only control because cell type, perturbation and experimental context differ. Its result therefore cannot by itself establish a causal Sendai effect.

---

## 4. Control hierarchy

The analysis must report all applicable levels separately:

```text
Observed signal
      |
      +--> Null 1: temporal-order null
      |
      +--> Null 2: monotonic-shape control
      |
      +--> Context/biological control
      |
      +--> Candidate-level specificity
```

No single level may be presented as proof of the others.

---

## 5. Null 1 — temporal-order null

### Question

Does correct temporal ordering contain information beyond the same observations with time labels disrupted?

### Construction

For each dataset, permute the ordering of the four valid day labels while preserving the activity values. The observed statistic is calculated once using the prespecified real ordering.

Because the primary Yamanaka trajectory has exactly four unique days, there are only **24 exact orderings**. The observed ordering is excluded from the null reference, leaving **23 null permutations**. No Monte Carlo approximation is used for this four-day case.

### Decision rule

The observed statistic passes Null 1 only if:

- it lies above the **99th percentile** of the absolute-value temporal-order null distribution, and
- the +1-corrected empirical two-sided permutation p-value is **≤ 0.05**.

With 23 null permutations, the minimum attainable +1-corrected empirical p-value is **1/24 ≈ 0.0417**. Therefore a 0.01 threshold is not statistically resolvable for this design and is not used.

### Interpretation

Passing means that the observed ordering is unusual under arbitrary time assignment. It does **not** mean that the process is specific to reprogramming.

---

## 6. Null 2 — monotonic-shape control

### Question

Is cross-dataset agreement stronger than expected for generic monotonic temporal processes?

A simple time-label permutation is insufficient for this question because it destroys monotonicity and can create an artificially weak null.

### Primary construction

Construct matched monotonic control trajectories using the same four time points (0, 3, 7, 10) and the same feature count as the observed candidate set.

The preferred control is a **biological but non-reprogramming temporal signal** available in the analyzed data, such as a prespecified cell-cycle or generic stress-response activity score. The exact feature set must be fixed before running the analysis and documented in the output.

If no defensible non-reprogramming temporal control exists in the input data, use a synthetic monotonic null generated from the empirical marginal distributions, with monotonic ordering imposed independently in each replicate. The synthetic construction must be fixed in code and reported explicitly.

### Decision rule

For each observed candidate statistic, compare against the matched monotonic-control distribution.

A candidate passes Null 2 if its statistic exceeds the **95th percentile** of the matched monotonic-control distribution.

At the aggregate level, the Stage 2.11C signal passes only if the median candidate statistic is above the 95th percentile of the aggregate control distribution.

### Interpretation

Passing means that the result is not adequately explained by generic monotonic shape alone.

It does **not** establish reprogramming specificity or causality.

---

## 7. Biological/context control — GSE304042

### Purpose

Test whether candidate pathway/TF activities are reproducible under a different biological system and perturbation structure.

GSE304042 contains human RPE samples with GFP, OCT4, SOX2, KLF4 and OSK conditions.

### Critical limitation

GSE304042 differs from GSE297234 in at least:

- cell type;
- experimental design;
- perturbation structure;
- delivery/context.

Therefore it is a **heterologous biological/context control**, not a clean Sendai-vs-non-Sendai experiment.

### Prespecified candidate comparison

Candidates are grouped into:

1. **conserved-context candidates** — activity direction is compatible across systems;
2. **context-dependent candidates** — signal is present in the original system but not reproduced in the heterologous system;
3. **inconclusive candidates** — the comparison cannot distinguish biology from technical/context differences.

### Candidate survival threshold

A candidate enters the conserved-context set only if:

- its direction is concordant in both applicable datasets;
- the magnitude is not effectively zero in either dataset under the prespecified activity normalization;
- it remains within the top **25%** of candidates after applying the same scoring procedure to the heterologous context.

The top-25% rule is a ranking filter, not a statistical significance claim.

### Stage-level decision rule

- **≥60%** of Tier-A candidates survive: proceed with the conserved/context-specific decomposition.
- **30–59%** survive: classify the result as **MIXED** and redesign the representation/model to explicitly separate conserved and context-dependent components.
- **<30%** survive: classify the current candidate interpretation as **CONFOUNDED/UNSUPPORTED** and do not promote these candidates to the next mechanism-discovery stage.

These thresholds are fixed before inspecting Stage 2.11C results.

---

## 8. Candidate tiers

### Tier A — reproducible temporal signal

A candidate must satisfy all of:

1. valid temporal measurements in both primary trajectory datasets;
2. concordant direction;
3. temporal ordering consistent with the observed transition;
4. pass Null 1 at the candidate level;
5. pass Null 2 at the candidate level;
6. stable ranking across the prespecified bootstrap/replicate resampling procedure.

Tier A means **reproducible temporal candidate**, not mechanism.

### Tier B — context-supported candidate

Tier A plus:

1. survives the heterologous biological/context comparison;
2. remains in the prespecified top-25% ranking band;
3. is not classified as an obvious generic delivery/stress/cell-cycle confound;
4. has an explicit evidence-for/evidence-against record.

Tier B means **candidate temporal regulator with cross-context support**.

### Tier C — mechanistic candidate

Tier C is intentionally **not attainable from Stage 2.11C alone**.

It additionally requires:

1. intervention-specific temporal data;
2. quantitative prediction of an intervention-specific state change;
3. evidence that perturbing the candidate changes the predicted transition in the expected direction;
4. preferably independent experimental validation.

Only Tier C candidates may later be described as mechanistic candidates.

---

## 9. Candidate-level statistical multiplicity

PROGENy pathways and DoRothEA regulons are not independent observations. Therefore:

- `14/14 pathways` and `257/294 TFs` must remain descriptive summaries;
- they must not be interpreted as 308 independent confirmations;
- candidate-level inferential tests must control multiplicity using **Benjamini–Hochberg FDR < 0.05** where enough independent units exist for a meaningful test;
- correlation structure among regulons should be reported rather than ignored.

If the effective number of independent features cannot be estimated robustly, the result must remain explicitly descriptive.

---

## 10. Stability criterion

A candidate is considered stable only if its rank remains in the same broad tier in at least **80% of 1,000 prespecified resampling replicates**.

The resampling unit must be the independent biological sample/group, not individual cells created by pseudobulk expansion.

No resampling result may be used to alter the thresholds defined in this protocol.

---

## 11. Decision tree

```text
                    Stage 2.11C
                         |
             +-----------+-----------+
             |                       |
        Null 1 passes?          Null 1 fails
             |                       |
             v                       v
        Null 2 passes?          UNSUPPORTED
             |
       +-----+-----+
       |           |
      yes          no
       |           |
       v           v
 biological      SHAPE-
 context         EXPLAINED
 control
       |
 +-----+----------------+
 |                      |
 >=60%                  30-59%
 |                      |
 v                      v
PROCEED                MIXED
 |                      |
 v                      v
conserved/context      redesign
 decomposition         representation
 |
 +----------+
 |          |
 <30%       |
 |
 v
CONFOUNDED / UNSUPPORTED
```

A result classified as MIXED is not a failure. It means the next model must explicitly represent context dependence rather than assuming invariance.

---

## 12. Stop/redesign rules

### CLEAN / PROCEED

Requirements:

- Null 1 passes;
- Null 2 passes;
- ≥60% of Tier-A candidates survive the heterologous context control;
- no dominant known confound explains the majority of Tier-B candidates.

Action:

Proceed to intervention-conditioned temporal modeling.

### MIXED

Any of the following:

- 30–59% candidate survival;
- strong conserved signal in some pathways/TFs but clear context dependence in others;
- evidence that the signal contains both biological and delivery/context components.

Action:

Do not discard the project. Redesign the representation to separate conserved and context-specific components before mechanism ranking.

### CONFOUNDED / UNSUPPORTED

Any of the following:

- Null 1 fails;
- Null 2 fails for the aggregate signal;
- <30% candidate survival;
- the majority of high-ranked candidates map to a single unsupported confound class.

Action:

Do not promote the current candidates. Revisit the representation and control strategy.

### INCONCLUSIVE

If metadata, sample matching, or biological controls are insufficient to apply the decision rules reliably:

Action:

Do not force a positive or negative conclusion. Improve metadata/control coverage and repeat the prespecified analysis.

---

## 13. What this protocol does NOT establish

Passing Stage 2.11C does **not** establish:

- cell lineage;
- causality;
- a discovered molecular mechanism;
- universal reprogramming dynamics;
- validity of a Digital Biological Twin.

It establishes only that the current temporal candidate signals are sufficiently robust and specific to justify the next experimental question.

---

## 14. Required output of the implementation

The future implementation must write a machine-readable result containing at least:

- protocol version/hash;
- random seeds;
- dataset/sample identifiers;
- exact control feature definitions;
- number of valid permutations/resamples;
- observed statistics;
- null/control quantiles;
- empirical p-values where applicable;
- BH-adjusted q-values where applicable;
- candidate tier;
- candidate survival rate;
- aggregate decision: `PROCEED`, `MIXED`, `CONFOUNDED_UNSUPPORTED`, or `INCONCLUSIVE`.

The implementation must never modify the thresholds in this document based on observed results.

---

## 15. Scientific interpretation template

If the result is positive, use language of the form:

> “The analysis supports reproducible, temporally ordered and partially context-robust pathway/TF activity patterns. These features are treated as candidate temporal regulators and motivate intervention-conditioned forecasting; they do not establish causal mechanisms.”

If mixed:

> “The analysis indicates that the observed temporal program contains both conserved and context-dependent components. This supports an explicit context-aware state representation rather than a universal mechanism assumption.”

If unsupported:

> “The current candidate signals do not survive the prespecified temporal/shape/context controls. We therefore do not promote them to mechanistic candidates and treat this as evidence that the current representation or control strategy is insufficient.”

---

## 16. Relationship to the PhD hypothesis

Stage 2.11C is a **gate**, not the PhD itself.

The central PhD contribution remains the development and validation of dynamic biological state representations capable of separating:

- conserved state-transition structure;
- context-specific structure;
- perturbation effects;
- history-dependent effects;

and using these representations for out-of-context future-state prediction.

The Yamanaka system is one biologically meaningful test case for this methodology.

Digital Biological Twin remains a downstream demonstrator rather than the primary scientific claim.
