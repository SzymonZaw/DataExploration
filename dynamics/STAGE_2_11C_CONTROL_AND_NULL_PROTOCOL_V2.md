# Stage 2.11C v2 — Prospective Control and Null Protocol

## Status

**Prospective revision.** This document does not rewrite the historical Stage 2.11C protocol. It defines the corrected rules for any future validation run after the audit.

The historical `28 -> 1` result remains exploratory because the v1 document and implementation disagreed on the attainable Null-1 p-value threshold.

## 1. Scientific question

Determine whether pathway/TF activity signals in the Yamanaka trajectory are reproducible, temporally ordered, statistically supported and sufficiently specific to justify intervention-conditioned dynamical modeling.

This protocol does not establish causality, lineage, a molecular mechanism, or a Digital Biological Twin.

## 2. Null 1 — temporal-order test

With four time points, enumerate all 24 permutations and exclude the observed ordering, leaving 23 null permutations.

For each candidate:

- compute the observed cross-dataset activity correlation;
- compute the exact null distribution over the 23 alternative orderings;
- report the +1-corrected empirical p-value;
- report the attainable p-value grid.

Decision rule:

- `passes_null1_raw` requires the observed absolute correlation to exceed the 99th percentile of the absolute null distribution and `p <= 0.05`;
- `passes_null1_fdr` additionally requires Benjamini–Hochberg `q < 0.05` across the prespecified candidate family.

No claim of `p < 0.01` is permitted with only 23 null permutations.

## 3. Null 2 — monotonic-shape control

Retain the existing synthetic monotonic control construction as an explicitly synthetic null until a biological non-reprogramming control is formally prespecified.

A candidate passes Null 2 when its observed correlation exceeds the 95th percentile of its matched monotonic-control distribution.

The output must state that this is a shape control, not a biological control.

## 4. Gene-ID mapping integrity

GSE297233 is Ensembl-based and must be mapped to HGNC symbols before PROGENy/DoRothEA scoring.

The implementation must report:

- input Ensembl IDs;
- mapped IDs;
- unmapped IDs;
- unique HGNC symbols;
- many-to-one collision count;
- aggregation rule;
- mapping artifact hash or persisted mapping table;
- annotation/service version where available.

Duplicate HGNC symbols are aggregated by the prespecified mean rule. A future confirmatory run must use a persisted mapping artifact rather than relying only on a live annotation query.

## 5. Context control

GSE304042 is a heterologous context control. It is **not** a Sendai-only control.

Because cell type, perturbation structure (OSK versus OSKM), delivery and experimental design differ, the comparison can test cross-system generalization but cannot identify a unique confounder.

Required output field:

`confound_source = unresolved`

The analysis must not label the observed failure as a Sendai effect, interferon effect, Myc effect, or cell-type effect without an experiment that separates those axes.

## 6. Candidate tiers

### Tier A

A candidate must satisfy:

1. valid measurements in both primary trajectory datasets;
2. concordant direction;
3. pass Null 1 raw;
4. pass Null 1 BH/FDR;
5. pass Null 2;
6. stable direction under the prespecified bootstrap procedure.

Tier A means **statistically supported reproducible temporal candidate**, not mechanism.

### Tier B

Tier A plus:

1. concordant activity direction in the heterologous context;
2. worst-context rank fraction <= 0.25;
3. not in the predefined delivery-confound panel;
4. explicit evidence-for/evidence-against record.

Tier B means **cross-context candidate**, not causal mechanism.

### Tier C

Not attainable from Stage 2.11C alone. Requires intervention-specific temporal validation and predictive intervention response.

## 7. Stage-level decision

The survival threshold is retained as an **operational prospective rule**, not a statistically derived effect-size threshold:

- `>= 60%` Tier-B/Tier-A: `PROCEED`;
- `30–59.9%`: `MIXED`;
- `< 30%`: `SPECIFICITY_UNSUPPORTED`;
- insufficient control/metadata: `INCONCLUSIVE`.

The thresholds are locked here before any prospective v2 validation run.

If no Tier-A candidates remain after FDR correction, report `SPECIFICITY_UNSUPPORTED` with an explicit reason that no statistically supported candidate survived the inferential gate; do not infer a particular confounder.

## 8. Multiplicity

PROGENy pathways and DoRothEA regulons are correlated and are not treated as independent biological confirmations.

BH-FDR is applied to the prespecified Null-1 candidate p-value family. The number `1/28` is not itself an inferential p-value and must not be Bonferroni-corrected as though it were.

## 9. Provenance

The original v1 protocol was committed before the Ensembl normalization work. The audit must preserve that historical fact.

This v2 protocol is prospective and must not be used to retroactively reclassify the historical `28 -> 1` run as confirmatory.

## 10. Required machine-readable fields

Every future v2 run must include:

- `protocol_version`;
- `protocol_hash`;
- `decision_rule_status`;
- `confound_source`;
- mapping integrity metrics;
- exact permutation count and attainable p-value grid;
- raw p-values and BH q-values;
- candidate tier;
- Tier-A and Tier-B counts;
- survival rate;
- aggregate decision;
- explicit reason when no candidate survives.

## 11. Interpretation guardrail

Use:

> The revised analysis tests whether candidate temporal signals survive technically validated gene mapping, exact temporal-order testing, multiplicity control, shape controls and an explicitly limited heterologous context comparison.

Do not use:

> The analysis identified the confounder responsible for the temporal program.
