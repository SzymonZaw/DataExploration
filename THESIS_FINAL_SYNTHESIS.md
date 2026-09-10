# Thesis final synthesis — Z1–Z6

## 1. Purpose

This document consolidates the current empirical status of objectives Z1–Z6 into one methodological argument for the doctoral thesis:

> **Reprezentacja dynamicznego stanu biologicznego na podstawie heterogenicznych danych omicznych powinna być oceniana nie jako pojedynczy model predykcyjny, lecz jako wieloetapowy system walidacji obejmujący reprezentację, integrację, dynamikę, specyficzność biologiczną, transferowalność i falsyfikację mechanistyczną.**

The synthesis deliberately separates three claims that are often conflated:

1. **representational validity** — whether heterogeneous molecular measurements can be placed in a reproducible state space;
2. **predictive validity** — whether that state representation contains information useful for predicting future observations;
3. **biological/mechanistic validity** — whether the predictive structure corresponds to transferable biological state rather than context-specific or technical structure, and whether proposed mechanisms are supported by temporal evidence.

The distinction is central to the thesis. External validation is a stricter test of generalizability than performance within the development dataset, and failures of external validation can themselves provide scientific information. This principle is particularly important for high-dimensional omics data, where platform, batch and experimental-context effects can limit transferability. citeturn0search0turn0search2

## 2. Central thesis claim

The current evidence supports the following bounded methodological claim:

> **A dynamic biological-state representation can be constructed from heterogeneous omics measurements and can recover reproducible temporal structure and context-specific predictive information. However, predictive information is not automatically equivalent to transferable biological state, and temporal or multimodal association is not sufficient to establish causal mechanism. Therefore, robust representation of dynamic biological state requires explicit separation of representation, prediction, specificity, transferability and mechanistic evidence.**

This is stronger and more defensible than claiming that one universal latent trajectory has been discovered.

## 3. Objective-by-objective synthesis

| Objective | Scientific question | Current status | Main evidence | Boundary |
|---|---|---|---|---|
| **Z1** | Can biological state be represented dynamically? | **SUPPORTED AS METHODOLOGICAL FRAMEWORK** | reproducible state-space/trajectory analyses across selected systems | representation is not automatically biologically specific |
| **Z2** | Can heterogeneous omics be integrated into a common representation? | **SUPPORTED** | multimodal concordance and validated identifier/provenance chain in GSE67462 | integration establishes concordance, not causal coupling |
| **Z3** | Can temporal structure be modeled reproducibly? | **SUPPORTED WITH CONTEXT LIMITATION** | stable temporal modules in GSE67462; reproducible trajectories in independent datasets | temporal structure can be system-specific |
| **Z4** | Can biological signal be separated from technical/context effects? | **PARTIALLY SUPPORTED / UNRESOLVED** | negative controls, signal attribution, context audits, explicit confounding tests | universal biological specificity has not been established |
| **Z5** | Does the representation transfer across independent biological systems? | **NOT SUPPORTED AS UNIVERSAL** | strict cross-dataset LODO failure; heterogeneous within-system results | transferability is context-dependent |
| **Z6** | Does the representation support prediction of future cellular state? | **SUPPORTED CONDITIONALLY** | GSE67462 within-system predictive support; strict LODO negative across systems | prediction is system-specific and does not establish mechanism |

## 4. Z1 — dynamic biological-state representation

### Claim supported

The project demonstrates a reproducible framework for representing temporal biological state using molecular measurements rather than treating each time point as an isolated sample.

The representation is based on a state variable whose observed modalities are treated as noisy projections of an underlying dynamic biological state:

```text
x_d(t) = G_d(z(t), u_d(t), h_d(t)) + ε_d(t)
```

where `z(t)` denotes the latent biological state, `u_d(t)` biological inputs, `h_d(t)` context-specific effects and `ε_d(t)` measurement noise.

The practical contribution of Z1 is therefore not the claim that `z(t)` has been uniquely identified. It is the establishment of a reproducible analytical object on which temporal and predictive hypotheses can subsequently be tested.

### Boundary

A low-dimensional representation can be stable while still encoding context, batch or experiment-specific structure. Consequently, representation stability is a prerequisite for biological interpretation, not proof of it.

## 5. Z2 — integration of heterogeneous omics

### Claim supported

GSE67462 provides the strongest current demonstration of multimodal integration. After resolving the GPL19972 feature namespace and linking expression identifiers to mm9 transcript/TSS coordinates, 11,048 of 11,899 common-space expression genes were validated (92.85%).

The regulatory validation identified reproducible dynamic concordance for:

- H3K27ac;
- H3K4me3;
- RNAPII;
- total OCT4.

H3K27me3 was retained as a negative-control modality and did not support the active multimodal interpretation.

The assignment remained robust over multiple TSS-centered mapping radii for the supported modalities.

### Methodological contribution

The important contribution is the **provenance-aware integration chain**:

```text
expression feature
      ↓
GPL19972 platform identifier
      ↓
gene symbol
      ↓
mm9 transcript/TSS
      ↓
regulatory element assignment
      ↓
multimodal temporal concordance
```

This makes identifier provenance part of the scientific evidence rather than an invisible preprocessing step. This is consistent with established recommendations that omics analyses document data provenance, processing, software and sources of variation. citeturn0search1turn0search3

### Boundary

Multimodal concordance does not demonstrate that one modality causes another. It establishes that the modalities carry coherent information about the same temporal transition under the tested experimental context.

## 6. Z3 — temporal structure

### Claim supported

The GSE67462 multimodal core can be partitioned into reproducible temporal modules. The selected six-module solution showed high clustering stability, with the principal biological modules represented by M1, M3, M4 and M6.

The most interpretable modules were:

- **M1:** early transient, dominated by ECM/mesenchymal remodeling, growth-factor and receptor-associated programs;
- **M4:** late-rising, associated with ECM, cytoskeletal, trafficking and RAC1/RHO/integrin/RTK-related programs;
- **M6:** late-rising, large module with strong RNAPII and OCT4 association and concurrent epithelialization/metabolic signatures.

Functional enrichment and category-specific analyses support these as **temporal hypotheses**, not mechanistic labels.

### Boundary

Temporal clustering does not imply a universal ordering of biological events. In particular, the M6 OCT4 association is compatible with reprogramming but does not establish a pluripotency mechanism, because epithelialization, cornified-envelope and cholesterol-metabolism signals provide competing interpretations.

## 7. Z4 — biological specificity versus context

### Claim supported

The project establishes a framework for explicitly testing whether predictive or temporal structure is specific to biological state rather than simply being reproducible temporal structure.

The GSE67462 signal-attribution audit showed that predictive information is predominantly carried by a common dynamic component rather than branch-specific expression structure. This is evidence against a trivial branch-specific explanation, but it is not proof of biological specificity.

The broader H3/control program also demonstrated why positive temporal similarity alone is insufficient: the circadian control had a different temporal geometry, while the target commensurability audit showed that not every candidate control is suitable for formal comparison.

### Current status

Z4 therefore remains **partially resolved rather than closed as proven**.

The correct methodological conclusion is:

> reproducibility of a temporal signal is necessary but not sufficient for biological specificity.

This distinction is consistent with severe-testing approaches in high-dimensional omics, which emphasize explicit hypotheses, falsification and protection against incremental corroboration of weakly specified claims. citeturn0search8

## 8. Z5 — cross-system transferability

### Claim supported

Strict leave-one-dataset-out validation failed to establish universal future-state prediction for any tested representation/model.

The within-system audit then separated two hypotheses:

- H-Z6a: no predictive information even within a biological system;
- H-Z6b: predictive information exists within systems but does not transfer universally.

GSE67462 provided formal within-system support for both autoencoder and PCA, whereas GSE28688 and GSE297234 did not. Therefore H-Z6a is too strong as a general explanation.

### Current interpretation

The preferred interpretation is **context-dependent transferability**:

> predictive information can exist in a learned state representation within a biological system, while remaining non-transferable across the tested independent systems.

This is not an anomalous outcome for omics prediction. Cross-platform and cross-experiment transferability is a known challenge because molecular measurements can vary systematically between datasets, platforms and experimental environments. citeturn0search9

### Boundary

The result does not prove that transferability is impossible. It establishes that universal transferability was not demonstrated under the frozen protocol.

## 9. Z6 — prediction of future state

### Claim supported

Z6 demonstrates **conditional predictive validity**.

The frozen Phase 0 benchmark failed under strict leave-one-dataset-out validation. However, Phase 0.1 demonstrated genuine within-system predictive support in GSE67462 for autoencoder and PCA under the unchanged support rule, while GSE28688 and GSE297234 remained unsupported.

Thus the representation is not intrinsically non-predictive, but its predictive information is context-dependent.

### Mechanistic falsification

The mechanistic extension deliberately imposed a stronger requirement: regulatory information should provide reproducible evidence of temporal precedence relative to later expression.

That requirement was not met:

- aggregate directional analysis showed no active modality/module with reproducible gain above the circular-shift null;
- replicate-level directional analysis reproduced the absence of positive temporal precedence across the two GSE67462 expression replicates;
- the only formal positive directional result occurred for H3K27me3 in M4, which is the negative-control modality and therefore cannot support the proposed activating mechanism.

Consequently:

> **Z6 supports context-dependent prediction, but does not support a causal or mechanistic interpretation of the inferred transitions.**

This separation is methodologically important because predictive association and mechanistic validity are distinct evidentiary claims.

## 10. The resulting methodological architecture

The completed workflow can now be represented as a validation ladder:

```text
                 HETEROGENEOUS OMICS
                         │
                         ▼
              Z1  STATE REPRESENTATION
                         │
                         ▼
              Z2  MULTIMODAL INTEGRATION
                         │
                         ▼
              Z3  TEMPORAL STRUCTURE
                         │
                         ▼
          Z4  BIOLOGICAL SPECIFICITY AUDIT
                         │
                         ▼
            Z5  EXTERNAL TRANSFERABILITY
                         │
                         ▼
              Z6  FUTURE-STATE PREDICTION
                         │
                         ▼
             MECHANISTIC FALSIFICATION
```

The crucial property of this architecture is that **success at one level does not automatically validate the next level**.

For example:

```text
multimodal concordance
        ≠
causal mechanism

predictive performance
        ≠
biological specificity

within-system prediction
        ≠
cross-system transferability

temporal module stability
        ≠
universal biological trajectory
```

This hierarchy is the central methodological result of the project.

## 11. Main scientific contribution

The primary contribution should therefore be framed as a **methodological framework for evidence-controlled representation of dynamic biological state**, rather than as a single predictive model.

The framework contributes five principles:

### P1 — State representation must precede prediction

A prediction model should operate on an explicitly defined state representation rather than on unconstrained high-dimensional temporal measurements.

### P2 — Multimodal integration requires provenance

Cross-omics integration must retain identifier, platform, coordinate and mapping provenance so that apparent biological concordance cannot arise from an undocumented namespace transformation.

### P3 — Predictive validity must be separated from biological validity

A model can predict a transition because it learned context, experiment, protocol or other reproducible structure. Predictive success therefore requires independent-context testing before being interpreted biologically.

### P4 — Transferability is an empirical property, not an assumption

A representation should be tested across independent systems rather than assuming that a successful within-system model captures a universal biological state variable.

### P5 — Mechanistic interpretation requires a stronger test than prediction

Multimodal correlation, enrichment and temporal association generate hypotheses. Causal interpretation requires additional evidence such as perturbation, temporal precedence that survives appropriate null models, intervention-response relationships or independent mechanistic validation.

## 12. What the thesis does not claim

The thesis should explicitly avoid the following claims:

- that one universal latent trajectory of cellular reprogramming has been discovered;
- that GSE67462 establishes a causal OCT4/chromatin mechanism;
- that predictive information is necessarily biological;
- that the learned representation transfers across species, platforms or protocols;
- that temporal modules are themselves causal biological entities;
- that statistical significance of multimodal association establishes mechanism;
- that failure of external prediction means biological dynamics are absent.

These boundaries are scientifically important. In omics predictor development, locked-down validation, independent external evaluation, data provenance and explicit documentation of sources of variation are recognized as essential safeguards against optimistic interpretation. citeturn0search3turn0search4

## 13. Final thesis-level conclusion

The current evidence supports the following final statement:

> **Dynamic biological state can be represented from heterogeneous omics data in a reproducible, provenance-aware state space that captures temporal and multimodal structure. Such representations can contain genuine predictive information within specific biological systems, but predictive information is not universally transferable and cannot be assumed to represent biological state independently of experimental context. In the tested GSE67462 system, multimodal temporal modules are reproducible and biologically interpretable at the level of hypotheses, while mechanistic temporal-precedence tests do not support a causal regulator-to-expression interpretation. The resulting methodology therefore treats biological-state representation as a hierarchy of separately validated claims—representation, integration, dynamics, specificity, transferability and mechanism—rather than as a single model-performance endpoint.**

## 14. Implications for the Digital Biological Twin

The Digital Biological Twin remains a **downstream application**, not the primary scientific claim.

The current evidence suggests that a future twin should not be represented as a single universal trajectory. Instead, it should maintain:

- a state representation;
- context metadata;
- uncertainty estimates;
- provenance for each modality;
- a transferability profile;
- model applicability boundaries;
- explicit distinction between observational prediction and mechanistic knowledge.

This follows directly from the Z1–Z6 evidence: the useful object is not merely a predicted future state, but a **prediction accompanied by knowledge of when, why and under which context the prediction is valid**.

## 15. Recommended thesis structure

A coherent dissertation structure emerging from the current evidence is:

1. **Problem formulation** — dynamic biological state as a latent, context-dependent object.
2. **Heterogeneous omics integration** — common state representation and provenance.
3. **Temporal modeling** — dynamic structure and module discovery.
4. **Specificity and confounding** — biological versus technical/contextual signal.
5. **Transferability** — independent-system validation.
6. **Prediction** — future-state forecasting under frozen validation.
7. **Mechanistic falsification** — stronger tests of temporal precedence and intervention relevance.
8. **Methodological synthesis** — evidence hierarchy for dynamic biological-state modeling.
9. **Digital Biological Twin demonstrator** — downstream application constrained by the validated evidence hierarchy.

## 16. Current project status

**Z1–Z3:** sufficiently mature for thesis-level synthesis.

**Z4:** partially resolved; biological specificity remains bounded and context-dependent.

**Z5:** universal transferability not supported; context-dependent transferability is the preferred interpretation.

**Z6:** conditional predictive validity supported; universal prediction and mechanistic validity not supported.

The next development phase should therefore prioritize **thesis-level consolidation, reproducibility packaging and explicit evidence accounting**, rather than further model complexity or repeated hypothesis-specific tuning of the same GSE67462 directional analysis.
