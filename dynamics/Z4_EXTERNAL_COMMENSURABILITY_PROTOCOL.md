# Z4 External Commensurability Protocol

## Purpose

This stage precedes any frozen transfer score. It asks whether GSE297234 and GSE28688 are sufficiently commensurate with the frozen GSE67462 discovery object to make a specificity claim scientifically interpretable.

**No transfer score, target similarity, module score or validation p-value may be inspected before these gates are frozen.**

## Frozen discovery object

- Discovery: GSE67462 + GSE67520.
- Biological target: OSKM-associated cellular-state transition, with GSE67462 as the discovery trajectory.
- Frozen multimodal core: H3K27ac, H3K4me3, RNAPII, OCT4.
- Frozen feature universe: the validated GSE67462 expression universe mapped through GPL19972 to mm9 TSS space.
- Frozen temporal modules: M1-M6; primary modules M1, M3, M4, M6.
- Z4 does not refit these objects.

## Candidate-specific target constructs

### GSE297234

Target construct: **human partial OSKM reprogramming / mesenchymal-drift reversal**.

This is a homologous biological transition, not a literal replication of full mouse secondary reprogramming.

Required evidence:

1. samples span at least two ordered states/timepoints;
2. OSKM exposure is documented;
3. fibroblast/somatic starting state is documented;
4. donor identity is available;
5. independent phenotype or source-study endpoint exists for mesenchymal-drift/rejuvenation;
6. the endpoint can be kept independent from the transferred score.

### GSE28688

Role: **predefined specificity/context challenge**, not a positive target replication.

Required evidence:

1. OSKM exposure and early temporal sampling;
2. biological replicate structure;
3. documented viral/transduction context;
4. redox/stress/p53 context is observable or documented;
5. endpoint samples exist but are not treated as equivalent to the discovery transition without an explicit biological rationale.

## Gates

### C1 — temporal commensurability

Pass if the external dataset has an ordered biological transition with >=3 ordered states or timepoints. A single endpoint comparison is insufficient for the primary Z4 transfer.

### C2 — cell-state commensurability

Pass if starting cells are comparable somatic/fibroblast states and the transition is induced rather than merely comparing unrelated cell types.

### C3 — intervention commensurability

Pass if OSKM or a clearly documented homologous reprogramming intervention is present. For GSE297234, Sendai-OSKM is acceptable as homologous but not identical to GSE67462.

### C4 — feature-space commensurability

Pass only after deterministic human→mouse orthology mapping is available. Mapping must be defined before validation results are inspected.

Required mapping properties:

- one frozen mapping release/source;
- case normalization only where deterministic;
- one-to-many orthologs handled by a pre-specified rule;
- no validation-outcome-based mapping selection;
- mapping coverage reported over the frozen GSE67462 feature universe.

### C5 — endpoint independence

Pass if the external biological endpoint is not mathematically derived from the frozen transfer features. Source-study phenotypes, donor-defined states or independently measured markers are preferred.

### C6 — nuisance observability

Pass if donor, time, protocol/intervention and platform/technical metadata are available. For GSE297234, donor age must remain an explicit covariate rather than being silently absorbed into the target definition.

### C7 — specificity-control compatibility

For GSE28688, the dataset must share the broad OSKM/time structure while exposing a documented acute stress/redox/p53 response. This makes it a meaningful specificity challenge.

### C8 — frozen decision

The final commensurability decision is made from provenance/design metadata and mapping coverage only. No transfer metric is allowed to change the decision.

## Decision states

- `ELIGIBLE_FOR_FROZEN_TRANSFER`: all mandatory gates pass.
- `ELIGIBLE_WITH_LIMITATIONS`: target is scientifically useful but one non-critical dimension is weaker; must be labelled in advance.
- `INELIGIBLE`: a mandatory gate fails.
- `UNRESOLVED`: required provenance/mapping information is missing.

## Z4 promotion remains separate

Passing commensurability does **not** support Z4. It only authorizes the next locked analysis.

The subsequent frozen transfer must separately test:

1. positive external target association;
2. target versus predefined non-target/context challenge;
3. nuisance robustness;
4. biological replicate consistency.

## Current expected decision

GSE297234 should be evaluated first as the primary homologous external target.

GSE28688 should be evaluated as the predefined context challenge.

GSE67462 remains discovery-only.
