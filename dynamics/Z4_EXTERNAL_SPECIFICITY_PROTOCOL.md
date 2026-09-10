# Z4 External Biological Specificity Protocol

## 1. Purpose

This protocol defines the experiment required to move Z4 from **PARTIALLY SUPPORTED / UNRESOLVED** to **SUPPORTED**.

The question is not whether a temporal signature is reproducible. GSE67462 already demonstrates reproducible temporal and multimodal structure. The Z4 question is stronger:

> **Does the learned state/transition signature correspond to a biological process that remains identifiable in an independent experiment, rather than merely encoding the original experimental context, protocol, batch or generic time structure?**

External validation must therefore be performed on a genuinely independent specimen set with the representation and decision rule locked before inspection of the validation outcome. Independent external validation and strict lock-down are established safeguards against optimistic omics interpretation. citeturn0search0turn0search11

## 2. Current Z4 boundary

Current evidence supports:

- reproducible common dynamic structure in GSE67462;
- multimodal concordance across H3K27ac, H3K4me3, RNAPII and OCT4;
- stable temporal modules;
- predictive information within GSE67462.

Current evidence does **not** establish that this structure is specific to biological reprogramming independently of experimental context.

The Z4 experiment must therefore be designed as an **external specificity test**, not as another within-dataset correlation analysis.

## 3. Frozen discovery object

The discovery system is GSE67462.

The following objects are frozen and must not be refit using the external validation data:

1. expression feature universe and identifier mapping;
2. GPL19972 → gene symbol → mm9 TSS mapping;
3. multimodal core definition;
4. six-module temporal solution and module assignments;
5. selected biological signature(s);
6. normalization/scaling parameters used by the transfer procedure;
7. state-distance / signature-score definition;
8. biological endpoint definition;
9. negative-control definitions;
10. acceptance thresholds.

Any modification after validation data are inspected constitutes a new model version and requires a new validation set.

## 4. Required external validation design

### 4.1 Minimum design

The validation experiment should contain at least:

- one independent biological experiment;
- at least one independent biological replicate per condition/timepoint, preferably multiple biological replicates;
- the same biological transition or a clearly defined homologous transition;
- a biological endpoint independent of the learned signature;
- technical metadata sufficient to test batch/protocol effects;
- a negative-control process that shares generic temporal or stress structure but is not expected to reproduce the target transition.

The strongest design uses a different experimental batch and, where feasible, a different cell preparation/protocol or laboratory while preserving the biological question. Technical variation between experiments can materially alter omics measurements, so counterbalancing and explicit batch accounting are essential. citeturn0search4turn0search2

### 4.2 Preferred design: crossed specificity matrix

The ideal validation matrix is:

| System | Target transition | Non-target process | Purpose |
|---|---:|---:|---|
| Discovery GSE67462 | + | | model construction only |
| Independent reprogramming experiment | + | | positive external validation |
| Independent differentiation/stress/proliferation experiment | | + | process-specificity control |
| Technical/context perturbation | | + | nuisance/context control |

The non-target controls should be selected before analysis and should share plausible nuisance structure such as time dependence, proliferation, stress, or major protocol changes. A useful specificity test asks whether the representation identifies the target while remaining uninformative for these alternative explanations.

## 5. Modalities

### Tier 1 — sufficient for the core Z4 test

A new experiment with **RNA expression plus an independent biological endpoint** can test whether the learned state representation transfers to a biological state.

### Tier 2 — stronger multimodal validation

If available, add one or more regulatory modalities such as:

- H3K27ac;
- H3K4me3;
- RNAPII;
- OCT4 or another transition-specific regulator.

The regulatory modalities should be treated as validation measurements, not as evidence merely because they were used to define the discovery signature.

### Important constraint

A second RNA-only dataset can validate state specificity, but it cannot replicate the full GSE67520 regulatory measurement chain. This distinction must remain explicit in the final interpretation.

## 6. Frozen transfer procedure

For each external sample:

1. apply the frozen identifier mapping;
2. restrict to the frozen validated feature universe;
3. apply frozen preprocessing parameters where mathematically appropriate;
4. compute the frozen state representation/signature score;
5. calculate distance or similarity to the frozen target-state trajectory;
6. evaluate against the independently defined biological endpoint;
7. evaluate target-vs-control discrimination;
8. test residual association with batch, protocol, donor, sequencing run and other nuisance variables.

No clustering, feature selection, module discovery, threshold optimization or latent-space refitting may use the validation samples.

## 7. Primary Z4 endpoint

The primary endpoint is **biological specificity**, not prediction RMSE.

A successful result should demonstrate all three:

### A. Positive external validity

The frozen signature/state score identifies the target biological transition in the independent experiment.

### B. Negative-control discrimination

The same score does not identify the non-target process at a comparable rate or effect size.

### C. Context robustness

The target association remains after accounting for experimental context and technical nuisance variables.

These three conditions are jointly required. A positive A alone is insufficient because a context-specific signature may also transfer to a similar protocol.

## 8. Quantitative acceptance rule

Z4 is promoted to **SUPPORTED** only if all mandatory gates pass.

### Gate Z4-A — external target association

The frozen state score must discriminate target biological state/transition from its external control with:

- pre-registered effect direction;
- permutation or appropriate exact/randomization p < 0.05;
- effect remaining positive under bootstrap uncertainty (95% CI excluding the null, where applicable).

### Gate Z4-B — target versus non-target specificity

The target effect must be significantly stronger than the strongest predefined non-target process:

`effect_target > effect_non_target`

with a pre-registered paired/permutation comparison p < 0.05.

If several non-target controls exist, the primary comparison uses the maximum non-target effect; secondary results are multiplicity-adjusted.

### Gate Z4-C — context robustness

After adjustment/stratification for batch, donor, protocol and other predefined nuisance variables:

- target association remains in the same direction;
- the adjusted effect retains at least 50% of the unadjusted effect;
- no single nuisance variable explains the majority of the target score variance.

The 50% rule is a pragmatic pre-specified robustness criterion, not a universal statistical law.

### Gate Z4-D — replication

If the external experiment has ≥2 biological replicates, the target effect must have the same direction in each replicate and must not be driven by one replicate.

If only one independent biological replicate exists, the result may be labelled **Z4_PROMISING_EXTERNAL_SIGNAL**, but cannot promote Z4 to fully supported.

## 9. What does NOT count as Z4 support

The following are insufficient on their own:

- significant correlation with time;
- similarity of PCA trajectories;
- high correlation of gene-wise trajectories;
- successful prediction in another dataset without a biological endpoint;
- multimodal correlation within the same experiment;
- enrichment of target-related gene sets;
- preservation of temporal module labels;
- improvement after tuning on the validation dataset;
- absence of an obvious batch effect in PCA/UMAP alone.

These can be supporting diagnostics but not the primary specificity claim.

## 10. Failure interpretations

### A — external target fails

Interpretation: the learned representation may be context-specific or insufficiently transferable.

### B — target and non-target both score positively

Interpretation: the signature likely captures a generic temporal/stress/proliferative/context component rather than target-specific biology.

### C — target survives but disappears after nuisance adjustment

Interpretation: the observed specificity was substantially explained by experimental context.

### D — target transfers across experiments but not across biological systems

Interpretation: **contextual biological specificity**, not universal biological specificity.

### E — target transfers and discriminates non-target processes

Interpretation: strongest current evidence for **biological specificity of the representation**.

## 11. Relation to Z5 and Z7

Z4 is not identical to Z5 or Z7.

**Z4:** Does the signal correspond to the target biological process rather than a generic/context signal?

**Z5:** Does the representation transfer across independent biological systems?

**Z7:** Does an intervention establish causal mechanism?

Therefore:

```text
Z4 specificity
    ≠ Z5 universal transferability
    ≠ Z7 causality
```

A representation can be biologically specific within a defined process while remaining non-universal across species, cell types or protocols.

## 12. Recommended experimental priority

The preferred next experiment is **not another GSE67462 re-analysis**.

Priority should be:

1. identify an independent OSKM/reprogramming dataset with a biological endpoint;
2. freeze the GSE67462-derived representation before opening the validation data;
3. add a biologically distinct temporal control with comparable sampling structure;
4. execute the frozen transfer;
5. report positive transfer, negative-control discrimination and nuisance robustness together.

If no suitable public dataset satisfies these constraints, the scientifically cleaner option is to label Z4 unresolved and design a prospective external experiment rather than relax the acceptance rule.

## 13. Decision states

```text
Z4_NOT_TESTED
      │
      ▼
Z4_EXTERNAL_VALIDATION_READY
      │
      ├── target fails ───────────────► Z4_UNSUPPORTED_IN_CURRENT_FORM
      │
      ├── target positive only ───────► Z4_EXTERNAL_SIGNAL
      │
      ├── target + controls pass ─────► Z4_CONTEXT_SPECIFIC_SUPPORTED
      │
      └── target + controls + robust
          external replication ───────► Z4_SUPPORTED
```

## 14. Final thesis criterion

The thesis may state **Z4 = SUPPORTED** only if an independent validation experiment demonstrates that the frozen representation tracks the target biological transition, discriminates predefined non-target processes, and remains robust to measured experimental context.

Otherwise the thesis should retain the more conservative statement:

> **The framework establishes an explicit test for biological specificity, but current evidence does not demonstrate that the observed predictive structure is independent of experimental context.**

This conservative boundary is intentional. In high-dimensional biological representation learning, IID hold-out performance can remain strong while out-of-distribution generalization fails because nuisance variation is entangled with the biological signal. citeturn0search2turn0search6
