# Z4 External Validation Candidate Audit

## Current status

Candidate selection is frozen. The next stage is a **commensurability audit only**; no frozen transfer score may be evaluated before the gates in `Z4_EXTERNAL_COMMENSURABILITY_PROTOCOL.md` pass.

## Frozen candidate roles

- **GSE297234** — primary homologous external target: human fibroblast partial OSKM reprogramming / mesenchymal-drift reversal.
- **GSE28688** — predefined specificity/context challenge: human early OSKM induction with a documented acute viral/redox/p53 stress component.
- **GSE67462** — discovery only.

## Why GSE297234 is the first target

NCBI GEO describes GSE297234 as 10x Genomics scRNA-seq of human fibroblasts from young and aged donors treated with Sendai-virus OSKM, with samples at D0, D3, D7 and D10. The study explicitly frames partial reprogramming as reducing mesenchymal drift before dedifferentiation/pluripotency acquisition. citeturn0search1

This is a biologically homologous transition, but not a literal replication of mouse secondary reprogramming. Species, donor age, Sendai delivery, partial-reprogramming endpoint and single-cell technology are deliberate domain shifts.

## Why GSE28688 is a specificity challenge

NCBI GEO describes GSE28688 as human foreskin fibroblasts measured at baseline and 24/48/72 h after OSKM transduction, plus iPSC and hESC samples. The source study specifically reports that viral transduction perturbs redox homeostasis and activates p53, senescence and apoptosis. citeturn0search0turn0search4

Therefore it is useful as a predeclared challenge against the hypothesis that a representation merely detects generic OSKM exposure, acute stress or elapsed time.

## Mandatory next action

Run `run_z4_external_commensurability_audit.py` separately for GSE297234 and GSE28688 after obtaining:

1. sample metadata tables;
2. frozen human→mouse feature mapping;
3. the frozen GSE67462 feature list;
4. pre-specified endpoint/nuisance metadata.

The script must return `ELIGIBLE_FOR_FROZEN_TRANSFER` before the transfer stage begins.

A failed commensurability gate is **not** a negative biological result. It means the dataset cannot answer Z4 under the frozen design.

## Scientific boundary

Passing commensurability does not support Z4. It only authorizes the locked transfer analysis. Z4 support additionally requires positive external target association, discrimination from the predefined specificity challenge, nuisance robustness and biological-replicate consistency.
