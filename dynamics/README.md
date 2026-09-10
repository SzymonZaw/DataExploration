# Dynamic state model

This directory now contains the first implementation of the central AI model proposed for the PhD project.

## Scientific role

The goal is to learn, rather than manually specify, a low-dimensional representation of cellular state and its temporal dynamics from harmonized biological observations.

The model is designed around:

```text
observation x(t)
      ↓
encoder
      ↓
latent state z(t)
      ↓
transition model
      ↓
predicted z(t+Δt)
      ↓
decoder
      ↓
reconstructed observation
```

Optional context/perturbation `u(t)` and history `h(t)` can be provided to the transition model. This allows future experiments to test whether cell-state dynamics are Markovian or require information about previous states.

## Modules

- `dynamic_state_model.py` — encoder, decoder, temporal transition and optional history/context representation.
- `losses.py` — reconstruction and latent prediction objectives, including missing-value masks.
- `train_dynamic_state.py` — leakage-aware training scaffold. Data splitting, normalization and imputation deliberately remain outside the trainer.
- `model_benchmark.py` — unified leakage-free forecasting benchmark used by Z6.
- `run_model_benchmark.py` — frozen Phase 0 benchmark runner.
- `Z6_PREDICTIVE_TRANSITION_PROTOCOL.md` — doctoral Z6 protocol for out-of-sample prediction of future cellular state.
- `Z6_WITHIN_SYSTEM_TRANSFERABILITY_PROTOCOL.md` — frozen Phase 0.1 within-system validation protocol.
- `Z6_CONTEXT_SPECIFICITY_AUDIT_PROTOCOL.md` — diagnostic protocol for testing reproducibility versus branch/context dependence.
- `run_z6_context_specificity_audit.py` — executable branch-concordance audit following Phase 0.1.
- `run_z6_gse67462_identifier_mapping_audit.py` — diagnostic audit of GPL19972/RefSeq feature provenance and mapping to GSE67520 gene/TSS identifiers.
- `Z6_GSE67462_IDENTIFIER_MAPPING_AUDIT_PROTOCOL.md` — protocol and interpretation rules for the identifier-provenance audit.
- `run_z6_gse67462_mechanistic_interpretation.py` — diagnostic ranking of multimodally coherent genes and optional pathway/nuisance enrichment.
- `Z6_GSE67462_MECHANISTIC_INTERPRETATION_PROTOCOL.md` — protocol and interpretation boundaries for mechanistic hypothesis generation.
- `run_z6_gse67462_temporal_modules.py` — clusters the multimodal mechanistic core into reproducible temporal expression modules.
- `Z6_GSE67462_TEMPORAL_MODULE_PROTOCOL.md` — protocol for temporal-module stability and interpretation.
- `run_z6_gse67462_temporal_enrichment.py` — reproducible GO/Reactome/KEGG/nuisance over-representation audit for temporal modules using explicit GMT files.
- `Z6_GSE67462_TEMPORAL_ENRICHMENT_PROTOCOL.md` — statistical and provenance protocol for temporal-module enrichment.
- `../run_z6_predictive_transition.py` — dedicated Z6 executable wrapper.

## Z6 predictive transition benchmark

Z6 is now operationally defined as a genuine future-state forecasting experiment:

```text
held-out trajectory prefix
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

The benchmark uses leave-one-dataset-out validation, training-only preprocessing, persistence/nearest-time/linear baselines and a temporal permutation null. Candidate models are evaluated under the same forecasting protocol so that the representation/model, rather than the evaluation procedure, is the experimental variable.

Run the official benchmark with:

```powershell
python run_z6_predictive_transition.py
```

For a diagnostic smoke test only:

```powershell
python run_z6_predictive_transition.py --epochs 25 --permutation-n 50 --seeds 411
```

The official Z6 result is not established until the full benchmark has been executed and its predefined predictive-support criteria have been evaluated.

## Z6 context-specificity audit

After Phase 0.1, predictive support was observed only in GSE67462 for autoencoder and PCA. Predictive performance alone is not treated as biological specificity. The context audit therefore compares the two independent branches within each system using gene-level temporal concordance, endpoint direction concordance, endpoint-effect correlation and per-timepoint molecular concordance.

Run it with:

```powershell
python -m dynamics.run_z6_context_specificity_audit
```

The audit is diagnostic only: it does not alter Phase 0/0.1 predictive-support decisions or thresholds. Weak branch concordance in a predictive system is evidence for context-dependent predictive structure, not evidence of biological absence.

## GSE67462/GSE67520 identifier provenance audit

GSE67462 is GPL19972, Brainarray `MoGene10stv1_Mm_REFSEQ version 18`, and its processed GEO feature IDs are RefSeq-like identifiers such as `NM_001001130.2_at`. The Stage 2.6 common-space matrix currently uses gene symbols. The multimodal regulatory analysis therefore requires an explicit feature-provenance audit before interpreting expression/ChIP concordance biologically.

Run:

```powershell
python -m dynamics.run_z6_gse67462_identifier_mapping_audit `
  --gtf Data\GSE67520\mm9.refGene.gtf.gz `
  --platform-soft Data\GPL19972_family.soft.gz
```

This audit is diagnostic only. It does not change frozen Z6 predictive support or any threshold.

## GSE67462 mechanistic interpretation audit

The current multimodal result establishes robust within-system molecular coherence, but not causality. The mechanistic audit ranks genes that are jointly coherent across H3K27ac, H3K4me3, RNAPII and OCT4, while keeping H3K27me3 as a contrasting negative-control modality.

Run:

```powershell
python -m dynamics.run_z6_gse67462_mechanistic_interpretation `
  --gtf Data\GSE67520\mm9.refGene.gtf.gz `
  --platform-soft Data\GPL19972_family.soft.gz
```

The audit produces gene-level multimodal scores, candidate regulatory-lead genes, and optional pathway/nuisance enrichment. Pathway enrichment is performed only when an explicitly versioned GMT file is supplied; no external database is silently downloaded.

The output is **hypothesis-generating**. An OCT4-associated or chromatin-associated trajectory does not establish causal regulation.

## GSE67462 temporal modules

The next interpretation layer groups the multimodal core into temporal modules before assigning biological pathway labels. This prevents individual-gene cherry-picking and separates early, transient and late expression programs using the actual observed trajectory.

Run:

```powershell
python -m dynamics.run_z6_gse67462_temporal_modules `
  --gtf Data\GSE67520\mm9.refGene.gtf.gz `
  --platform-soft Data\GPL19972_family.soft.gz
```

The analysis tests several cluster counts and reports adjusted-Rand agreement between them. Module labels such as `early_declining`, `late_rising` and `transient` are descriptive summaries of the observed centroids, not mechanistic assignments. Pathway interpretation should follow this stability audit and use an explicitly versioned gene-set collection.

## GSE67462 temporal-module enrichment

Functional interpretation is performed only for modules with at least 20 genes. The current six-module solution therefore treats M1, M3, M4 and M6 as primary enrichment targets and keeps M2/M5 as small modules without stable pathway inference.

The enrichment audit uses the full provenance-validated GSE67462 expression universe as the hypergeometric background, a minimum overlap of five genes, and Benjamini-Hochberg correction across all module × term hypotheses within each gene-set collection. GO Biological Process, Reactome, KEGG and nuisance/context collections are supplied explicitly as versioned GMT files.

Example:

```powershell
python -m dynamics.run_z6_gse67462_temporal_enrichment `
  --gtf Data\GSE67520\mm9.refGene.gtf.gz `
  --platform-soft Data\GPL19972_family.soft.gz `
  --go Data\GeneSets\GO_BP_mouse.gmt `
  --reactome Data\GeneSets\Reactome_mouse.gmt `
  --kegg Data\GeneSets\KEGG_mouse.gmt `
  --nuisance Data\GeneSets\nuisance_mouse.gmt
```

The exact GMT files and SHA-256 hashes are recorded in the generated manifest. No pathway database is silently downloaded, and enrichment does not modify frozen Z6 support.

## Research constraints

The implementation must be evaluated with strict dataset-level and time-based holdouts. In particular, preprocessing statistics, feature selection and imputation parameters must be fitted on training data only.

A good reconstruction loss is not sufficient evidence of a biological state. The main evaluation target is prediction on observations that were not used to learn the representation.

Predictive performance alone is also not evidence of biological specificity. Z6 results must subsequently be audited for context-dependent and technical nuisance signals under the Z4 framework.

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
