# Z6 GSE67462 — mechanistic falsification protocol

## Purpose

This audit tests whether the strongest temporal modules in GSE67462 are better explained by specific multimodal transition programs than by generic temporal drift, stress/inflammation, proliferation, or metabolism.

It is **diagnostic only**. It must not modify frozen Z6 predictive support, multimodal support, temporal module assignments, or validation thresholds.

## Hypotheses

### H-M1 — early ECM/mesenchymal remodeling
M1 should show coherent early-transition behavior across expression, active chromatin, RNAPII and OCT4, with stronger support for ECM/signaling/adhesion genes than for generic nuisance programs.

### H-M4 — late structural remodeling
M4 should show coherent late-transition behavior concentrated in ECM, integrin, RHO/RAC/RAB and cytoskeletal/trafficking genes. Evidence should not depend on one modality alone.

### H-M6 — epithelialization versus pluripotency
M6 should be tested as a competing-hypothesis problem:

- H-M6a: epithelialization/cornification is the dominant interpretation;
- H-M6b: pluripotency-associated transition is dominant;
- H-M6c: OCT4-associated chromatin/transcriptional activity is present but not sufficient to establish pluripotency.

The audit must not label M6 as a pluripotency module solely because OCT4 correlates with its centroid.

## Falsification tests

1. **Multimodal coherence**
   - require positive temporal correlation across at least 3 active modalities for a candidate core gene;
   - report the negative-control H3K27me3 separately.

2. **Temporal direction**
   - calculate signed trajectory concordance with the module centroid;
   - distinguish early-transient from late-rising programs;
   - do not use correlation magnitude without direction.

3. **Specificity contrast**
   - compare candidate programs against stress/inflammation, proliferation/cell-cycle and metabolism categories;
   - report category overlap rather than filtering nuisance categories out.

4. **Expression-only versus multimodal support**
   - every candidate is classified as multimodal-core, multimodal-support, single-modal or unresolved;
   - a pathway supported only by expression is not treated as mechanistic evidence.

5. **Regulatory lead diagnostic**
   - compare contemporaneous modality→expression correlation with adjacent-time modality(t)→expression(t+1) correlation;
   - report only as hypothesis generation because the dataset has eight principal timepoints and irregular temporal spacing.

6. **Competing M6 interpretation**
   - quantify epithelialization, pluripotency/stem, proliferation and metabolism signals separately;
   - determine whether the leading pluripotency-associated genes are concentrated in a small minority of the module or form a broad coherent component;
   - explicitly report when evidence is insufficient to discriminate H-M6a/H-M6b/H-M6c.

## Decision language

- **SUPPORTED_AS_HYPOTHESIS**: multiple independent evidence channels agree and no tested nuisance explanation dominates.
- **MIXED_EVIDENCE**: biological and nuisance/context programs are both substantial.
- **UNSUPPORTED**: the proposed interpretation is not supported by the multimodal evidence.
- **UNRESOLVED**: evidence is insufficient to discriminate competing interpretations.

No category is declared absent merely because it fails enrichment.

## Required outputs

- candidate-level modality correlations;
- module-centroid direction scores;
- leading genes by competing interpretation;
- modality coverage;
- regulatory-lead diagnostic;
- M6 competing-hypothesis scorecard;
- machine-readable JSON/CSV outputs;
- bounded Markdown report.

## Provenance

The audit must reuse the validated GPL19972 → gene-symbol → mm9.refGene mapping already established for GSE67462/GSE67520. It must record input SHA256 values where available.
