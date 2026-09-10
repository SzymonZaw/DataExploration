# Z6 — Results and scientific findings

## Status

This document records the versioned empirical findings for objective Z6 and is now **scientifically closed for the current experimental program**. The final synthesis is recorded separately in `dynamics/Z6_FINAL_SYNTHESIS.md`.

Z6 tests whether a learned biological-state representation supports genuine out-of-sample prediction of future cellular state, whether predictive information transfers across systems, and whether the observed predictive/multimodal structure supports a defensible mechanistic interpretation.

## 1. Experimental logic

Z6 is evaluated as a forecasting problem:

```text
observed trajectory prefix
        ↓
state inference
        ↓
transition model
        ↓
free multi-step rollout
        ↓
predicted future state
        ↓
held-out future observation
```

The frozen predictive-support rule requires all of the following:

- positive mean improvement versus persistence;
- positive q05 improvement versus persistence;
- positive mean and q05 improvement versus nearest-time;
- positive mean and q05 improvement versus linear extrapolation;
- temporal permutation p < 0.05.

No dataset was removed from the primary benchmark because of a negative result, and subsequent audits did not change the model, preprocessing, baselines or support thresholds.

## 2. Phase 0 — leave-one-dataset-out benchmark

The frozen LODO benchmark did **not** establish predictive support for any tested model.

| Model | RMSE | Improvement vs persistence | q05 vs persistence | Improvement vs nearest | q05 vs nearest | Improvement vs linear | q05 vs linear | permutation p | Support |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| autoencoder | 0.901629 | 0.159065 | -0.052136 | 0.125055 | 0.033247 | 1.166434 | 0.435963 | 0.992008 | false |
| delta_t | 0.846834 | 0.213860 | -0.046771 | 0.179850 | 0.107191 | 1.221228 | 0.441327 | 0.992008 | false |
| markov | 0.843220 | 0.217474 | -0.058432 | 0.183464 | 0.113749 | 1.224843 | 0.429667 | 0.979021 | false |
| memory | 0.913102 | 0.147592 | -0.174806 | 0.113582 | 0.016180 | 1.154961 | 0.313292 | 0.996004 | false |
| pca | 1.107466 | -0.046772 | -0.316259 | -0.080782 | -0.417815 | 0.960597 | 0.256811 | 1.000000 | false |

Phase 0 is therefore a valid negative result under strict leave-one-dataset-out validation. It does not show that cellular dynamics are absent or that Z6 is impossible; it shows that the current representation and models do not demonstrate cross-dataset future-state prediction under the frozen protocol.

## 3. Phase 0.1 — within-system transferability

Phase 0.1 changed only the validation unit. Training on one independent temporal branch and predicting the other was performed while keeping the predictive task, preprocessing, models, baselines, permutation test and support rule frozen.

| System | Autoencoder | Delta-t | Markov | Memory | PCA |
|---|---|---|---|---|---|
| GSE28688 | false | false | false | false | false |
| GSE67462 | **true** | false | false | false | **true** |
| GSE297234 | false | false | false | false | false |

The decisive result is that predictive support is **system-specific rather than consistently reproduced across all three systems**.

The preferred empirical interpretation is therefore:

> Predictive information can exist in the learned representation within a biological system, but its transferability across independent biological and experimental systems is not demonstrated and appears to depend on context.

This is a working empirical conclusion, not a universal biological law.

## 4. Context specificity and signal attribution

Within GSE67462, predictive information is predominantly associated with a reproducible common dynamic expression component rather than branch-specific structure. The context-specificity audit classified GSE67462 as a reproducible transition, while the broader comparison showed mixed context dependence across systems.

This supports the interpretation that the positive GSE67462 prediction is not simply an arbitrary branch artifact. It does not establish that the common dynamic component is biologically specific rather than a reproducible property of this experimental system.

## 5. Multimodal validation

The GSE67462 transition shows robust multimodal dynamic concordance for:

- H3K27ac;
- H3K4me3;
- RNAPII;
- OCT4.

The validated expression universe contains 11,048 of 11,899 common-space genes (92.85%) after the explicit GPL19972 → gene-symbol → mm9 TSS mapping chain.

Assignment-robustness analysis supports the main multimodal result. H3K27me3 was retained as a negative control and did not provide positive multimodal support.

This establishes a reproducible multimodal temporal signature, but it remains observational association rather than causal evidence.

## 6. Temporal modules and functional structure

The multimodal core resolves into six stable temporal modules, with M1, M3, M4 and M6 sufficiently large for primary functional interpretation.

- **M1:** early transient ECM/mesenchymal remodeling, growth-factor and receptor/adhesion signaling. Strong EMT-Hallmark similarity is present, but transition direction is not established.
- **M4:** late-rising ECM, cytoskeletal, trafficking, RAC1/RHO, integrin and RTK-associated remodeling.
- **M6:** late-rising module with strong OCT4 and RNAPII association, together with epithelialization/cornified-envelope and cholesterol-metabolism signals. Pluripotency-related enrichment is present but not sufficiently specific to support a pluripotency claim.
- **M3:** mixed functional structure with weaker category-level support.

These results support a multi-phase description of the GSE67462 transition rather than a single scalar trajectory.

## 7. Mechanistic hypothesis and falsification audits

Mechanistic falsification supports descriptive hypotheses for M1 and M4, while M6 remains mixed and epithelialization-dominant.

Accordingly:

- M1 may be described as an ECM/mesenchymal-remodeling program;
- M4 may be described as structural/trafficking remodeling;
- M6 may be described as a late epithelialization-dominant program with OCT4-associated transcriptional structure and competing metabolic/pluripotency interpretations.

These are mechanistic hypotheses, not causal mechanisms.

### 7.1 Aggregate directional falsification

The directional audit tested whether regulatory trajectories at time `t` improve prediction of expression at `t+1` beyond same-time association and circular time-shift null expectations.

No primary module showed a positive lead-gain result that exceeded the null with statistical support.

Thus, concurrent multimodal association was reproducible, but temporal precedence was not demonstrated.

### 7.2 Replicate-level directional falsification

The replicate-level audit independently tested the expression response in the two GSE67462 expression replicates. For all active modalities in M1, M4 and M6, no active modality achieved both `above_null=True` in both replicates and `p<0.05` in both replicates.

The sole formal positive case was H3K27me3 in M4. Because H3K27me3 is the predefined negative-control modality, this result cannot be used as evidence for an activating regulatory mechanism.

The replicate audit therefore strengthens the negative directional conclusion.

The exact machine-readable outputs are retained under:

`results/Dynamics/z6_gse67462_replicate_directional_audit/`

with protocol:

`dynamics/Z6_GSE67462_REPLICATE_DIRECTIONAL_PROTOCOL.md`

## 8. Final scientific conclusion

> **The learned biological-state representation contains predictive information that is reproducible within selected biological systems, but this predictive information is not demonstrably transferable across the independent systems tested. In GSE67462, predictive information is associated with a reproducible multimodal temporal structure that resolves into stable functional modules. However, neither aggregate nor replicate-level directional analyses provide reproducible evidence that the measured regulatory signals temporally precede transcriptional changes. The representation therefore supports context-dependent predictive state modeling, while causal or mechanistic interpretation of the inferred transitions remains unsupported by the present observational evidence.**

## 9. What Z6 establishes

Z6 provides evidence that:

1. future-state prediction can succeed within a biological system under strict out-of-sample validation;
2. such predictive support is not uniformly transferable across independent systems;
3. predictive information in GSE67462 is associated with a reproducible common dynamic component;
4. GSE67462 contains a robust multimodal temporal signature across enhancer, chromatin, transcriptional and OCT4-associated measurements;
5. this signature can be decomposed into stable temporal modules with interpretable functional structure;
6. functional and mechanistic analyses can generate bounded hypotheses without being promoted to causal claims;
7. directional and replicate-level falsification provide evidence against claiming a validated regulator→expression temporal mechanism from these data alone.

## 10. What Z6 does not establish

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

## 11. Scientific implication for the thesis

The main methodological result of Z6 is not discovery of a universal trajectory. It is the demonstration that **predictive representation, biological specificity and mechanistic validity are separate validation layers**.

A state representation can be internally predictive but context-dependent; multimodally coherent without being causally validated; and temporally structured without supporting a universal ordering of molecular events.

This directly supports the thesis objective of separating biological state from context-dependent and technical information.

## 12. Frozen boundaries

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

The final synthesis is therefore an integration of versioned results, not a redefinition of the experiments.

## 13. Provenance

Primary predictive outputs remain under the corresponding `results/Dynamics/z6_*` directories. The methodological specifications and final synthesis are versioned under `dynamics/`.

See also:

`dynamics/Z6_FINAL_SYNTHESIS.md`

This document and the final synthesis together constitute the current closure of the Z6 experimental program.
