# Z6 — Context-specificity audit findings

## Status

This audit is a **diagnostic follow-up to Z6 Phase 0.1**. It does not change the frozen predictive-support rule, model benchmark, dataset inclusion, or previously reported Z6 conclusions.

The purpose is to test whether the positive within-system prediction observed in GSE67462 is accompanied by reproducible molecular transition structure between the two independent branches, rather than being treated automatically as evidence of biological specificity.

## 1. Execution

The audit was executed successfully with:

```powershell
python -m dynamics.run_z6_context_specificity_audit
```

The branch definitions are inherited directly from Z6 Phase 0.1. The run completed without a preprocessing warning or execution failure.

The underlying GSE67462 experiment is a mouse secondary-reprogramming time series with two samples per time point; GEO describes the design as bulk expression profiling during OSKM-mediated reprogramming, with samples collected across the reprogramming trajectory and a corresponding iPSC endpoint. citeturn0search0

## 2. Results

| Dataset | Context interpretation | Phase 0.1 predictive support | Positive models | Endpoint-effect permutation p |
|---|---|---:|---|---:|
| GSE28688 | MIXED_CONTEXT_DEPENDENCE | false | — | 0.000999 |
| GSE67462 | **REPRODUCIBLE_TRANSITION** | true | autoencoder; PCA | 0.000999 |
| GSE297234 | MIXED_CONTEXT_DEPENDENCE | false | — | 0.000999 |

The audit therefore identifies **GSE67462 as the only system in which predictive support coincides with a reproducible transition pattern under the predefined descriptive concordance criteria**.

## 3. Interpretation

The positive GSE67462 result is strengthened relative to a purely model-level observation: the two independent branches show sufficiently concordant molecular transition structure to meet the audit's descriptive `REPRODUCIBLE_TRANSITION` criteria.

However, this is **not equivalent to biological specificity**. A reproducible transition can still contain context-dependent, platform-dependent, experimental, or other nuisance structure. In particular, the audit establishes branch concordance, not causality or invariance across unrelated experiments.

GSE28688 and GSE297234 are classified as `MIXED_CONTEXT_DEPENDENCE`, not as failed biology. The label means that the observed branch agreement does not satisfy the predefined strong reproducibility criteria and does not meet the stricter branch-specificity criteria either. These systems therefore remain informative about context dependence but should not be overinterpreted as negative controls.

## 4. Scientific consequence for H-Z6b

The context audit modifies the wording of the current Z6 working interpretation.

The evidence still supports the statement that predictive information is **not demonstrably transferable across the tested independent systems**. However, for GSE67462 there is now evidence that the positive predictive result is accompanied by reproducible branch-level molecular transition structure.

Therefore the preferred formulation is:

> Predictive information can be reproducible within a biological system, but its transferability across independent systems remains unestablished and may depend on biological and experimental context.

This is stronger and more precise than describing the GSE67462 result simply as system-specific prediction.

## 5. Important limitation

The endpoint permutation p-value is `0.000999` for each of the three systems. This is the minimum attainable p-value for the implemented +1 correction with 1000 permutations. It should therefore be reported as `p = 0.000999` rather than as an exact smaller probability, and it should not be treated as evidence that the three systems have equivalent statistical strength.

The endpoint permutation is a **diagnostic concordance test**, not the Z6 predictive permutation test. It does not replace the frozen Phase 0.1 predictive-support criterion.

## 6. What this audit establishes

It establishes that:

- the audit is executable and reproducible under the frozen branch definitions;
- GSE67462 meets the predefined descriptive criteria for reproducible transition structure;
- the GSE67462 predictive result is accompanied by branch-level molecular concordance;
- GSE28688 and GSE297234 show intermediate/mixed context dependence rather than meeting the strong reproducibility criterion;
- the original Z6 predictive-support decisions remain unchanged.

## 7. What this audit does not establish

It does not establish:

- causal biological mechanisms;
- biological specificity of the predictive signal;
- universality of the transition;
- cross-species or cross-platform transferability;
- that GSE67462 represents a canonical cellular-state trajectory;
- that the learned latent state is invariant to experimental context.

Those claims require independent biological controls and mechanistic/perturbational evidence under Z4/H3.

## 8. Next scientific step

Do **not** expand the model architecture or relax predictive thresholds based on this result.

The next useful analysis is a **GSE67462 specificity decomposition**: determine which components of the reproducible branch transition are shared with known reprogramming biology and which components are associated with experimental context. The analysis should be performed without changing the frozen Z6 benchmark and should explicitly separate:

1. shared temporal biological signal;
2. branch-specific signal;
3. technical/contextual signal;
4. endpoint-only effects;
5. features that remain concordant under perturbation or independent validation.

The current audit therefore moves Z6 from a question of whether the GSE67462 result is merely branch-specific toward the more informative question of **which reproducible components are biologically specific**.

## 9. Provenance

Implementation:

`dynamics/run_z6_context_specificity_audit.py`

Protocol:

`dynamics/Z6_CONTEXT_SPECIFICITY_AUDIT_PROTOCOL.md`

Machine-readable outputs:

`results/Dynamics/z6_context_specificity_audit/`

Human-readable interpretation:

`dynamics/Z6_CONTEXT_SPECIFICITY_FINDINGS.md`
