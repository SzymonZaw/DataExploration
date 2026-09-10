# Z6 GSE67462 identifier mapping audit

## Purpose

This is a diagnostic audit for the GSE67462/GSE67520 multimodal validation path. It does **not** modify frozen Z6 predictive support, thresholds, model selection, or prior Phase 0/0.1 conclusions.

The audit was added because the initial GSE67520 regulatory analysis produced only 2–13 overlapping genes under promoter/TSS mapping. A subsequent namespace check showed that the Stage 2.6 common-space expression matrix uses uppercase gene symbols, whereas GSE67520 TSS annotation uses mixed-case mouse symbols. More importantly, NCBI documents GSE67462 as GPL19972, Brainarray `MoGene10stv1_Mm_REFSEQ version 18`, with processed feature IDs such as `NM_001001130.2_at`. Therefore the correct question is feature provenance, not simple string overlap.

NCBI sources:
- GSE67462: GPL19972, Brainarray MoGene10stv1_Mm_REFSEQ v18, 18 samples.
- GPL19972: platform table exposes RefSeq-like feature IDs.
- GSM1647454 and other samples: processed values use `ID_REF` such as `NM_001001130.2_at`.

## Command

```powershell
python -m dynamics.run_z6_gse67462_identifier_mapping_audit `
  --gtf Data\GSE67520\mm9.refGene.gtf.gz
```

If a local GPL19972 annotation table is available:

```powershell
python -m dynamics.run_z6_gse67462_identifier_mapping_audit `
  --gtf Data\GSE67520\mm9.refGene.gtf.gz `
  --platform-table Data\GSE67462\GPL19972.txt.gz
```

## Outputs

`results/Dynamics/z6_gse67462_identifier_mapping_audit/`

- `01_feature_mapping_audit.csv` — feature-level namespace and mapping diagnostics.
- `02_mapping_report.json` — provenance, counts, SHA-256 values and audit status.

## Interpretation rules

1. Exact overlap between uppercase Stage 2.6 symbols and mixed-case TSS symbols is not treated as biological evidence.
2. Case-insensitive symbol overlap is useful as a namespace diagnostic, but is not sufficient to establish that the expression feature was correctly mapped from GPL19972.
3. RefSeq/probe-to-symbol mapping is considered the preferred provenance path when a GPL19972 annotation table is available.
4. One-to-many mappings are not silently collapsed into a biological claim; they remain visible in the audit table.
5. Missing local GPL19972 annotation is reported as `PARTIAL_PLATFORM_TABLE_MISSING`, not as a negative biological result.
6. The audit must be completed before interpreting the GSE67520 multimodal concordance as a biological validation result.

## Scientific status

The prior multimodal result remains `EXPRESSION_ONLY` / `MAPPING_INADEQUATE` until this provenance question is resolved. No Z6 support criterion is changed by this audit.
