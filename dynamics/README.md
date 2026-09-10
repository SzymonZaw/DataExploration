# Dynamic state model

This directory now contains the first implementation of the central AI model proposed for the PhD project.

## Scientific role

The goal is to learn, rather than manually specify, a low-dimensional representation of cellular state and its temporal dynamics from harmonized biological observations.

## Z6 final synthesis

The complete Z6 evidence chain is summarized in:

- `Z6_RESULTS_AND_FINDINGS.md` — versioned empirical results and bounded interpretation.
- `Z6_FINAL_SYNTHESIS.md` — final synthesis of predictive, multimodal, temporal and mechanistic evidence.
- `../THESIS_FINAL_SYNTHESIS.md` — thesis-wide synthesis mapping Z1–Z6 into one evidence hierarchy.

The methodological conclusion is that representation, multimodal integration, temporal structure, biological specificity, transferability, prediction and mechanism are **separate validation layers**. Success at one layer does not automatically validate the next.

## Z4 external specificity

The GSE297234 donor-stratified transfer audit is now complete as a diagnostic layer. Four of six frozen modules showed the same temporal direction in the aged and young donors, and three of six were positive in both donors. This supports partial donor-level transferability, but does not establish biological specificity because no independent external endpoint was identified in the audited GSE297234 metadata.

The next frozen validation stage is an orthogonal chromatin-accessibility test using GSE242424/GSE242421. The prespecified protocol is:

- `Z4_GSE242424_ORTHOGONAL_ENDPOINT_PROTOCOL.md`
- `Z4_GSE242421_GENE_ACTIVITY_PROTOCOL.md`

The orthogonal assay is deliberately kept separate from the frozen expression representation. No module redefinition, feature selection or threshold tuning may use the target dataset. Gene activity is derived deterministically from the published peak-by-cell matrix and independent peak-gene links.

## Mechanistic falsification closure

The aggregate and replicate-level directional audits are now considered closed. They tested whether regulatory signals reproducibly precede later expression changes under a circular time-shift null.

The audits did **not** support a reproducible regulator-to-expression temporal-precedence effect for H3K27ac, H3K4me3, RNAPII or OCT4 in M1, M4 or M6. The only formal positive directional result was H3K27me3 in M4, which is the negative-control modality and is therefore not interpreted as an activating mechanism.

The replicate-level audit further showed that the absence of positive lead gain for the active modalities was stable across the two GSE67462 expression replicates. Because GSE67520 does not provide independent regulatory replicate trajectories, this is expression-response reproducibility rather than full replication of the regulatory measurements.

Therefore Z6 is **closed at the current mechanistic evidence level**: predictive validity is supported conditionally, but causal or mechanistic interpretation is not.

## Research constraints

The implementation must be evaluated with strict dataset-level and time-based holdouts. Preprocessing statistics, feature selection and imputation parameters must be fitted on training data only.

A good reconstruction loss is not sufficient evidence of a biological state. The main evaluation target is prediction on observations that were not used to learn the representation.

Predictive performance alone is also not evidence of biological specificity. Z6 results must be audited for context-dependent and technical nuisance signals under the Z4 framework.

## Planned model progression

The architecture is deliberately simple at first. It can later be extended with:

1. probabilistic latent states and uncertainty,
2. continuous-time / Neural ODE dynamics,
3. modality-specific encoders for RNA, scRNA and regulatory data,
4. perturbation embeddings,
5. explicit history/memory models,
6. sparse/Jacobian-based interpretability,
7. symbolic regression over learned latent dynamics.

These are model variants, not separate scientific stages. The scientific question remains constant: **can AI learn a biologically meaningful state representation whose dynamics generalize across independent experiments?**
