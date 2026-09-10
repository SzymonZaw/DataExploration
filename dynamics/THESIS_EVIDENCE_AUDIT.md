# Thesis-wide evidence audit — Z1–Z8

## 1. Purpose

This audit checks whether each doctoral objective has a defensible claim, machine-readable evidence, an explicit limitation, and a reproducible route back to the result. It is intentionally stricter than a narrative summary: a claim is considered complete only when the evidence level matches the wording of the claim.

The audit follows the same evidence hierarchy already established in the project: representation → integration → temporal structure → specificity/context → transferability → prediction → mechanism. This separation is consistent with current reproducibility guidance emphasizing transparent methods, provenance, independent validation and explicit reporting of failed reproducibility tests. citeturn0search0turn0search2

## 2. Executive verdict

| Objective | Status | Thesis-level verdict |
|---|---|---|
| Z1 | **SUPPORTED AS METHODOLOGICAL FRAMEWORK** | Ready for thesis synthesis |
| Z2 | **SUPPORTED** | Strongest positive evidence is GSE67462 multimodal integration |
| Z3 | **SUPPORTED WITH CONTEXT LIMITATION** | Stable temporal modules are demonstrated, not a universal trajectory |
| Z4 | **PARTIALLY SUPPORTED / UNRESOLVED** | Biological specificity remains bounded and context-dependent |
| Z5 | **NOT SUPPORTED AS UNIVERSAL** | Cross-system transferability was not demonstrated |
| Z6 | **SUPPORTED CONDITIONALLY** | Within-system prediction exists; universal prediction and mechanism do not |
| Z7 | **SUPPORTED AS HYPOTHESIS GENERATION, NOT AS CAUSAL DISCOVERY** | Mechanistic interpretation remains observational |
| Z8 | **DOWNSTREAM / NOT YET AN EMPIRICAL THESIS RESULT** | Digital Biological Twin should be a demonstrator constrained by Z1–Z7 evidence |

### Overall conclusion

**No additional experiment is required to make the core methodological thesis claim defensible.** The present evidence is sufficient to support a thesis centered on evidence-controlled dynamic-state representation.

There are, however, two genuine scientific gaps if the thesis were to claim more than that:

1. **Z4:** biological specificity has not been established independently of experimental context.
2. **Z7:** causal mechanism has not been demonstrated by intervention or independent mechanistic validation.

These are not reasons to keep tuning the existing GSE67462 analysis. They are explicit boundaries of the current contribution.

## 3. Evidence matrix

| Objective | Claim that survives audit | Primary evidence | Quantitative / decision evidence | Validation level | Limitation / forbidden overclaim | Thesis use |
|---|---|---|---|---|---|---|
| **Z1** | A reproducible dynamic state representation can be constructed as an analytical object. | Stage 2.4–2.9 validation; trajectory/state-space analyses; `THESIS_FINAL_SYNTHESIS.md` | Reproducible temporal structure across selected systems; leakage-controlled validation | Representation / methodological | Stable representation may encode context or technical structure | Core methodological foundation |
| **Z2** | Heterogeneous expression/regulatory observations can be integrated with explicit identifier provenance. | `run_z6_gse67462_identifier_audit.py`; mapped multimodal validation; multimodal robustness audit | 11,048 / 11,899 genes = **92.85%** validated; H3K27ac, H3K4me3, RNAPII, OCT4 supported; H3K27me3 negative control | Multimodal concordance + provenance | Concordance is not causal coupling | Strong empirical methods chapter |
| **Z3** | GSE67462 contains stable multi-phase temporal structure. | `run_z6_gse67462_temporal_modules.py`; temporal enrichment | 6 modules; median pairwise ARI **0.9464**; M1/M3/M4/M6 large enough for primary interpretation | Temporal reproducibility + functional interpretation | Module stability does not establish universal ordering or mechanism | Results + methodological synthesis |
| **Z4** | The project can explicitly test biological specificity versus nuisance/context, but specificity is not fully resolved. | H3 controls; GSE263713; target commensurability audit; context-specificity; signal attribution | GSE67462 common dynamic fraction median **0.9759**; branch fraction median **0.0241**; context results heterogeneous | Negative controls + attribution + context audit | Common dynamic signal may still be system-specific; no orthogonal biological specificity proof | Must be presented as partial/unresolved |
| **Z5** | Universal transferability is not demonstrated under frozen validation. | Phase 0 LODO + Phase 0.1 within-system audit | Phase 0: no model supported; Phase 0.1: GSE67462 AE/PCA supported, GSE28688/GSE297234 unsupported | External / cross-system validation | Failure of universal transfer does not prove impossibility | Important negative result |
| **Z6** | Future-state prediction can succeed conditionally within a biological system. | `Z6_RESULTS_AND_FINDINGS.md`; `Z6_FINAL_SYNTHESIS.md` | GSE67462 AE/PCA support; frozen LODO negative across datasets | Out-of-sample prediction | Prediction is context-dependent and not mechanistic | Core predictive chapter |
| **Z7** | Mechanistic analyses can generate bounded hypotheses and falsify an overstrong regulator→expression interpretation. | mechanistic interpretation, directional audit, replicate directional audit | No active modality/module showed reproducible lead gain above null; H3K27me3 M4 positive result is a negative control | Mechanistic hypothesis / falsification | No intervention, no independent regulatory replicates, no causal identification | Mechanistic limits chapter |
| **Z8** | A Digital Biological Twin is a valid downstream application architecture for the validated state representation. | `THESIS_FINAL_SYNTHESIS.md`; project README | No separate empirical performance claim required at current stage | Application / demonstrator | Must not be used to inflate biological evidence | Final demonstrator, not primary evidence |

## 4. Claim-to-evidence audit

### Z1 — representation

**Survives:** “The project establishes a reproducible methodological framework for dynamic biological-state representation.”

**Does not survive:** “The project uniquely identifies the true latent biological state.”

Reason: representation stability is necessary for subsequent testing but is not itself biological specificity.

### Z2 — integration

**Survives:** “Heterogeneous modalities can be placed into a common, provenance-aware analytical space.”

The strongest concrete chain is:

```text
GSE67462 expression feature
        ↓
GPL19972 platform identifier
        ↓
gene symbol
        ↓
mm9 transcript/TSS
        ↓
regulatory-element assignment
        ↓
multimodal temporal concordance
```

The 92.85% validated mapping coverage is sufficient for the stated integration claim, provided the thesis reports the mapping chain and its provenance.

**Does not survive:** “The modalities causally encode one common state variable.”

### Z3 — temporal dynamics

**Survives:** “GSE67462 contains reproducible multi-phase temporal modules.”

**Does not survive:** “The module sequence is a conserved universal biological trajectory.”

The current evidence supports a context-rich temporal description, not a universal ordering.

### Z4 — specificity

**Survives:** “The framework distinguishes reproducibility from biological specificity and provides explicit nuisance/control audits.”

**Does not survive:** “The predictive signal has been proven biologically specific.”

The common dynamic component is reproducible within GSE67462, but that is not equivalent to proving biological specificity. This is the main unresolved thesis objective.

### Z5 — transferability

**Survives:** “Universal cross-system transferability was not demonstrated under the frozen protocol.”

**Does not survive:** “The representation generalizes across independent biological systems.”

This negative result is scientifically useful rather than a failure of the project. External validation is precisely where context dependence becomes visible. citeturn0search2turn0search6

### Z6 — prediction

**Survives:** “The representation contains predictive information in selected systems.”

**Does not survive:** “The representation is a universally predictive state model.”

The correct wording is **conditional predictive validity**.

### Z7 — mechanism

**Survives:** “Mechanistic analyses generate testable hypotheses and can falsify stronger interpretations.”

**Does not survive:** “GSE67462 demonstrates a causal OCT4/chromatin→expression mechanism.”

The aggregate and replicate directional audits provide negative evidence against claiming temporal precedence under the applied methodology, while not proving absence of regulation.

### Z8 — Digital Biological Twin

**Survives:** “The evidence hierarchy defines constraints for a future Digital Biological Twin.”

**Does not survive:** “The Digital Biological Twin validates the biological model.”

The Twin remains a downstream application. It should expose state, context, uncertainty, provenance and applicability boundaries rather than collapse these into a single prediction.

## 5. Contradiction / stale-wording audit

### Resolved

The current `dynamics/README.md` explicitly states that the aggregate and replicate directional audits are closed and that Z6 is closed at the current mechanistic evidence level. fileciteturn289file0

`Z6_RESULTS_AND_FINDINGS.md` and `Z6_FINAL_SYNTHESIS.md` consistently state that:

- strict LODO did not support universal forecasting;
- GSE67462 provides conditional within-system prediction;
- multimodal concordance is observational;
- temporal precedence was not demonstrated;
- H3K27me3 is a negative control;
- no further lead-lag optimization is warranted. fileciteturn280file0turn285file0

### Remaining documentation issue

The root README still describes the **Stage 2.11 perturbation/mechanism-discovery program** as the current development direction. This is not scientifically wrong, but after Z6 closure it should be read as a historical/programmatic framing rather than evidence that another GSE67462 mechanism-tuning cycle is required. The thesis evidence audit therefore recommends **consolidation before new modeling**.

No result-level contradiction was found between the final Z6 documents and the thesis synthesis.

## 6. Reproducibility completeness

The core evidence chain has the required artifact classes:

- protocol documents;
- executable analysis scripts;
- machine-readable CSV/JSON results;
- bounded scientific findings;
- final synthesis documents;
- explicit negative controls and null models;
- identifier/provenance audit for multimodal mapping.

The project therefore satisfies the core structure expected for a reproducible computational evidence chain. Current reproducibility literature likewise emphasizes planning, methods, data management, analysis transparency and dissemination as minimum components rather than treating reproducibility as a final cosmetic step. citeturn0search0turn0search4

## 7. Scientific gaps that are real — and gaps that are not

### Real gaps

**G1 — independent biological specificity test**

Z4 remains partially unresolved. A future experiment that could genuinely close it would require an orthogonal biological control or perturbation that is expected to preserve/alter the biological transition while decoupling the major nuisance/context response.

**G2 — causal mechanism**

Z7 remains hypothesis-level. A causal claim would require intervention or another independent mechanistic validation, not another correlation/enrichment/lead-lag optimization on the same bulk time series.

### Not a required gap for the core thesis claim

**Not G3 — another neural architecture.**

Changing AE/PCA/Markov/memory architecture without a new biological validation question would not strengthen the central thesis claim.

**Not G4 — more permutations of the same test.**

The current negative directional result is already replicated at aggregate and expression-replicate levels.

**Not G5 — another candidate control without target commensurability.**

The H3 audit demonstrated that control selection must first pass target commensurability; adding datasets blindly would weaken, not strengthen, the logic.

**Not G6 — a Digital Biological Twin performance benchmark.**

The Twin is downstream of the validated methodology and should not be used as surrogate evidence for Z1–Z7.

## 8. Final evidence ladder for the dissertation

```text
Z1  REPRESENTATION
 │
 ▼
Z2  MULTIMODAL INTEGRATION + PROVENANCE
 │
 ▼
Z3  TEMPORAL STRUCTURE
 │
 ▼
Z4  SPECIFICITY / CONTEXT AUDIT
 │
 ▼
Z5  EXTERNAL TRANSFERABILITY
 │
 ▼
Z6  FUTURE-STATE PREDICTION
 │
 ▼
Z7  MECHANISTIC HYPOTHESIS + FALSIFICATION
 │
 ▼
Z8  DIGITAL BIOLOGICAL TWIN DEMONSTRATOR
```

The thesis should make explicit that this is **not a ladder where every rung must be positive**. A scientifically useful result may be a negative or bounded decision that prevents an overstrong downstream claim.

## 9. Recommended thesis wording

Use:

> **The study develops and validates an evidence-controlled framework for representing dynamic biological state from heterogeneous omics data. The framework demonstrates reproducible temporal and multimodal structure and conditional future-state predictive information within selected biological systems, while showing that such predictive information is not universally transferable across independent systems. Mechanistic analyses further demonstrate that multimodal association and temporal structure do not suffice to establish causal regulatory precedence. The principal contribution is therefore a methodology that explicitly separates representation, integration, specificity, transferability, prediction and mechanistic validity.**

Avoid:

> “We discovered the universal latent trajectory of cellular reprogramming.”

> “GSE67462 proves the OCT4/chromatin mechanism.”

> “The predictive representation is biologically specific.”

> “The failure of external prediction shows that biological dynamics are absent.”

## 10. Decision on further experiments

### Decision: **NO CORE EXPERIMENT IS MISSING**

The current evidence is sufficient for the core methodological thesis claim.

If additional scientific work is performed, it should be selected because it closes one of the two explicit boundaries:

1. an orthogonal specificity/perturbation experiment for Z4; or
2. an intervention-based mechanistic validation for Z7.

Everything else should be treated as optional model development or downstream application work.

## 11. Reproducibility note

This audit is a synthesis artifact. It does not modify frozen predictive thresholds, model architectures, preprocessing, dataset inclusion, multimodal support rules or earlier scientific decisions. Its purpose is to align thesis wording with the evidence already generated.

## 12. Current audit status

**THESIS EVIDENCE AUDIT: COMPLETE**

**Core thesis claim: SUPPORTED WITH EXPLICIT BOUNDARIES**

**Primary unresolved objective: Z4 biological specificity**

**Primary unclosed scientific claim: Z7 causal mechanism**

**Z5 universal transferability: NOT SUPPORTED**

**Z6 conditional prediction: SUPPORTED**

**Z8 Digital Biological Twin: DOWNSTREAM DEMONSTRATOR**
