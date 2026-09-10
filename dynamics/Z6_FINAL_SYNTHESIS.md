# Z6 — Final scientific synthesis

## Status

This document closes the current Z6 experimental program. It integrates the frozen predictive benchmark, within-system transferability audit, context-specificity analysis, signal attribution, multimodal validation, temporal module analysis, functional enrichment, mechanistic hypothesis audit, aggregate directional falsification and replicate-level directional falsification.

The purpose is not to claim that all questions about dynamic biological-state representation are solved. The purpose is to state precisely what the current evidence supports and where inference must stop.

## 1. Scientific question

Z6 asks whether a learned representation of cellular state contains information that supports genuine out-of-sample prediction of future cellular state and whether such predictive information is transferable and biologically interpretable.

The final evidence shows that these are separable questions. Predictive information can be present within a biological system without being demonstrably transferable across systems, and predictive or multimodal association does not by itself establish causal mechanism.

## 2. Evidence chain

### 2.1 Cross-dataset prediction

The frozen Phase 0 leave-one-dataset-out benchmark did not establish predictive support for any tested model. This is a valid negative result under the predefined protocol and was not rescued by changing thresholds, removing difficult datasets or relaxing baselines.

Therefore the current representation does not support a claim of universal cross-system forecasting.

### 2.2 Within-system prediction

Phase 0.1 changed only the validation unit while keeping the predictive task, preprocessing, models, baselines, permutation test and support rule frozen.

GSE67462 provided formal predictive support for the autoencoder and PCA representations. GSE28688 and GSE297234 did not.

This falsifies the strong interpretation that the representation is intrinsically non-predictive even within biological systems. At the same time, the system-to-system inconsistency prevents a claim of universal transferability.

The preferred empirical interpretation is therefore:

> Predictive information can exist in the learned representation within a biological system, but its transferability across independent biological and experimental systems is not demonstrated and appears to depend on context.

This is a working empirical conclusion, not a universal biological law.

### 2.3 Context specificity and signal attribution

Within GSE67462, predictive information is predominantly associated with a reproducible common dynamic expression component rather than branch-specific structure. The context-specificity audit classifies GSE67462 as a reproducible transition, while the broader Z6 comparison shows mixed context dependence across systems.

This supports the interpretation that the positive GSE67462 prediction is not simply an arbitrary branch artifact. It does not, however, establish that the common dynamic component is biologically specific rather than a reproducible property of this experimental system.

### 2.4 Multimodal validation

The GSE67462 transition shows robust multimodal dynamic concordance for:

- H3K27ac;
- H3K4me3;
- RNAPII;
- OCT4.

The validated expression universe contains 11,048 of 11,899 common-space genes (92.85%) after the explicit GPL19972 → gene-symbol → mm9 TSS mapping chain.

Assignment-robustness analysis supports the main multimodal result. H3K27me3 was retained as a negative control and did not provide positive multimodal support.

Thus the transition has a reproducible multimodal regulatory/transcriptional signature. This is stronger than expression-only trajectory similarity, but it remains an observational association.

### 2.5 Temporal modules and functional structure

The multimodal core resolves into six stable temporal modules, with M1, M3, M4 and M6 sufficiently large for primary functional interpretation.

The main patterns are:

- **M1:** early transient ECM/mesenchymal remodeling, growth-factor and receptor/adhesion signaling; strong EMT-Hallmark similarity but insufficient evidence to assign transition direction.
- **M4:** late-rising ECM, cytoskeletal, trafficking, RAC1/RHO, integrin and RTK-associated remodeling.
- **M6:** late-rising module with strong OCT4 and RNAPII association, together with epithelialization/cornified-envelope and cholesterol-metabolism signals. Pluripotency-related enrichment is present but not sufficiently specific to support a pluripotency claim.
- **M3:** mixed functional structure with weaker category-level support.

These results support a multi-phase description of the GSE67462 transition rather than a single scalar trajectory.

### 2.6 Mechanistic hypothesis audit

Mechanistic falsification supports descriptive hypotheses for M1 and M4, while M6 remains mixed and epithelialization-dominant.

Accordingly:

- M1 may be described as an ECM/mesenchymal-remodeling program;
- M4 may be described as structural/trafficking remodeling;
- M6 may be described as a late epithelialization-dominant program with OCT4-associated transcriptional structure and competing metabolic/pluripotency interpretations.

These are mechanistic hypotheses, not causal mechanisms.

## 3. Directional falsification

### 3.1 Aggregate directional audit

The directional audit tested whether regulatory trajectories at time `t` improve prediction of expression at `t+1` beyond same-time association and circular time-shift null expectations.

No primary module showed a positive lead-gain result that exceeded the null with statistical support.

The absence of a positive lead signal was observed despite strong concurrent multimodal correlations in several modules.

This distinction is important: **concurrent multimodal association was reproducible, but temporal precedence was not demonstrated.**

### 3.2 Replicate-level directional audit

The replicate-level audit independently tested the expression response in the two GSE67462 expression replicates.

For all active modalities in M1, M4 and M6:

- same-sign lead effects were generally reproducible between replicates;
- neither replicate consistently exceeded the circular-shift null;
- no active modality achieved `both_above_null=True` and `both_p<0.05=True`.

The sole formal positive case was H3K27me3 in M4. Because H3K27me3 is the predefined negative-control modality, this result cannot be used as evidence for an activating regulatory mechanism and instead reinforces the need for the negative-control boundary.

The replicate audit therefore strengthens the negative directional conclusion rather than weakening it.

### 3.3 Interpretation boundary

The correct conclusion is:

> No reproducible evidence was obtained for temporal precedence of the tested regulatory signals over later transcriptional changes in GSE67462 under the specified lead-lag methodology.

This does **not** establish absence of regulation. Sparse bulk timepoints, measurement resolution, aggregation, imperfect temporal alignment and the lack of independent regulatory replicates limit the sensitivity of the test.

## 4. Final Z6 conclusion

The complete Z6 evidence supports the following bounded conclusion:

> **The learned biological-state representation contains predictive information that is reproducible within selected biological systems, but this predictive information is not demonstrably transferable across the independent systems tested. In GSE67462, predictive information is associated with a reproducible multimodal temporal structure that resolves into stable functional modules. However, neither aggregate nor replicate-level directional analyses provide reproducible evidence that the measured regulatory signals temporally precede transcriptional changes. The representation therefore supports context-dependent predictive state modeling, while causal or mechanistic interpretation of the inferred transitions remains unsupported by the present observational evidence.**

## 5. What Z6 establishes

Z6 provides evidence that:

1. future-state prediction can succeed within a biological system under strict out-of-sample validation;
2. such predictive support is not uniformly transferable across independent systems;
3. predictive information in GSE67462 is associated with a reproducible common dynamic component;
4. GSE67462 contains a robust multimodal temporal signature across enhancer, chromatin, transcriptional and OCT4-associated measurements;
5. this signature can be decomposed into stable temporal modules with interpretable functional structure;
6. functional and mechanistic analyses can generate bounded hypotheses without being promoted to causal claims;
7. directional and replicate-level falsification provide evidence against claiming a validated regulator→expression temporal mechanism from these data alone.

## 6. What Z6 does not establish

Z6 does not establish:

- a universal biological trajectory;
- a universally transferable latent state;
- causal state transitions;
- a validated regulator→chromatin→expression chain;
- biological specificity of every predictive feature;
- pluripotency as the unique interpretation of M6;
- EMT or MET directionality from enrichment similarity alone;
- transferability across species, platforms or protocols;
- that negative directional evidence proves absence of molecular regulation.

## 7. Consequence for the thesis

The scientific contribution of Z6 is therefore methodological rather than the discovery of a universal trajectory.

The important result is that **predictive representation, biological specificity and mechanistic validity must be treated as separate validation layers**.

A state representation can be:

- internally predictive but context-dependent;
- multimodally coherent without being causally validated;
- temporally structured without supporting a universal ordering of molecular events.

This provides a concrete methodological basis for the thesis objective of separating biological state from context-dependent and technical information.

## 8. Final validation hierarchy

The Z6 workflow should be interpreted in the following hierarchy:

```text
Prediction
  ↓
Transferability
  ↓
Context dependence
  ↓
Multimodal concordance
  ↓
Temporal module stability
  ↓
Functional interpretation
  ↓
Mechanistic hypothesis
  ↓
Directional falsification
  ↓
Causal claim only with independent intervention evidence
```

Passing an earlier layer does not imply passing a later layer.

In particular:

```text
predictive ≠ biological-specific
multimodal ≠ causal
correlated ≠ temporally preceding
stable module ≠ mechanism
```

## 9. Frozen boundaries

The following were not changed by the final mechanistic audits:

- Phase 0 predictive-support thresholds;
- Phase 0.1 validation splits;
- model architectures;
- preprocessing rules;
- baseline definitions;
- permutation framework;
- dataset inclusion;
- multimodal support criterion;
- temporal module assignment;
- earlier positive/negative decisions.

The final synthesis is therefore a post hoc integration of versioned results, not a redefinition of the experiments.

## 10. Reproducibility and provenance

Primary predictive outputs remain under:

`results/Dynamics/z6_predictive_transition_audit/`

Within-system outputs remain under:

`results/Dynamics/z6_within_system_audit/`

GSE67462 context and signal-attribution outputs remain under the corresponding `results/Dynamics/z6_*` directories.

Multimodal, temporal-module, enrichment, specificity, mechanistic and directional outputs are likewise retained as machine-readable CSV/JSON artifacts alongside their protocol documents.

The replicate-level directional audit is stored under:

`results/Dynamics/z6_gse67462_replicate_directional_audit/`

with its protocol in:

`dynamics/Z6_GSE67462_REPLICATE_DIRECTIONAL_PROTOCOL.md`

## 11. Recommended next step

No further lead-lag optimization should be performed on GSE67462 as part of the current Z6 evidence chain.

The next scientific task is **cross-objective synthesis of Z1–Z6**: identify which claims survive all validation layers, define the methodological contribution of the thesis, and formulate the final model of context-aware dynamic biological-state representation.

The Digital Biological Twin remains a downstream demonstrator rather than evidence used to strengthen the biological claims above.
