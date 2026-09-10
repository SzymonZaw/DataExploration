# Z4 — GSE242424 orthogonal endpoint validation protocol

## Scientific purpose

The donor-stratified GSE297234 audit established partial transferability of frozen expression modules across two independent human donors, but it did not establish an independent biological endpoint. This protocol therefore tests a different validation layer: whether frozen expression modules from GSE67462 have reproducible counterparts in an **orthogonal chromatin-accessibility assay** from an independent human fibroblast reprogramming experiment.

The target resource is GSE242424 / GSE242421, which contains scATAC-seq across D0, D2, D4, D6, D8, D10, D12, D14 and iPSC, plus a paired scRNA subseries and D1/D2 multiome samples. NCBI documents the scATAC series as fragment BED files and the scRNA series as processed count matrices.

## Frozen hypothesis

H-Z4-ORTHO:

> Frozen expression-state modules identified in GSE67462 should show reproducible temporal association with an orthogonal chromatin-accessibility representation of the same reprogramming process, without fitting the module definition to GSE242424.

This is **not** a claim that chromatin accessibility is an endpoint in the clinical/outcome sense. It is an orthogonal molecular validation layer. A positive result strengthens biological-state interpretation; it does not by itself establish causality.

## Dataset requirements

### Required

- GSE242421 scATAC samples D0–D14 and iPSC.
- hg38 genomic coordinates.
- A reproducible gene-activity or regulatory-element representation derived without using the frozen GSE67462 module definitions.
- Explicit sample/time metadata.

### Preferred

Use processed gene-activity data if available. Otherwise derive gene activity from the published fragment BED files with a frozen hg38 gene annotation. NCBI reports the scATAC supplementary data as fragment BED files and identifies Chromap-based processing on hg38.

## Frozen analysis order

1. Audit file/sample provenance.
2. Construct orthogonal gene-activity representation independently of GSE67462 modules.
3. Map frozen GSE67462 module genes into the GSE242424 gene space.
4. Compute module scores per timepoint using the same frozen gene membership.
5. Compare temporal direction/rank with the GSE67462 frozen trajectory.
6. Use time-shift/permutation nulls; do not treat raw correlation as evidence.
7. Report module-level and aggregate concordance.

## Leakage prevention

The following must remain frozen:

- GSE67462 module membership.
- GSE67462 training-time feature selection.
- GSE67462 human↔mouse orthology mapping.
- Module direction definitions.

No feature selection, module refinement, threshold tuning or module redefinition may use GSE242424.

## Decision logic

### `Z4_ORTHO_SUPPORTED`
Requires:
- sufficient gene-space coverage,
- reproducible orthogonal representation,
- positive evidence for the frozen modules under the prespecified null,
- and no dependence on a single timepoint or single donor/branch.

### `Z4_ORTHO_PARTIAL`
Used when only a subset of modules passes the orthogonal validation while the remaining modules are unresolved or discordant.

### `Z4_ORTHO_UNRESOLVED`
Used when coverage, sample provenance, orthogonal representation or null calibration is inadequate.

### `Z4_ORTHO_FAILED`
Used only when the frozen modules show systematic lack of concordance despite adequate coverage and a valid null test.

## Important boundary

GSE297234 remains `Z4_UNRESOLVED` because it has no independent external endpoint in the currently audited metadata. The donor result (4/6 modules directionally concordant; 3/6 positive in both donors) is retained as **partial donor-level transferability**, not as specificity evidence.

The GSE242424 experiment is attractive because it provides a different molecular layer and a dense 14-day human fibroblast reprogramming time course. It should therefore be treated as the next orthogonal validation candidate, not retroactively used to change the frozen representation.
