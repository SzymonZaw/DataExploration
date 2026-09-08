# Stage 2.11C — Signal Attribution Dataset Manifest

Status: **prospectively locked before H3 analysis**

## Purpose

This manifest fixes the unrelated temporal negative controls for the Signal Attribution Audit. These datasets must not be replaced after inspection of their H3 results.

## Rejected datasets already present in `Data/`

The following local datasets were audited for suitability:

- **GSE148158** — human Yamanaka reprogramming; not unrelated to the target process.
- **GSE28688** — human OSKM reprogramming time course; useful as independent reprogramming/context evidence, but not an unrelated temporal control.
- **GSE52052** — human reprogramming cocktail comparison at day 11 only; not temporal and not unrelated.
- **GSE67462** — mouse secondary OSKM reprogramming time series; independent reprogramming evidence, not H3.
- **GSE67520** — mouse Oct4 ChIP-seq / chromatin data from the GSE67462 study; not an expression trajectory suitable for the primary H3 test.

## Locked unrelated temporal controls

### H3-Control-1 — GSE175533

**Process:** replicative cellular senescence / cellular lifespan

**Why selected:** independent biological process, human fibroblast system, high-resolution temporal design, and transcriptomic measurements. The study profiles replicative senescence in WI-38 human fetal lung fibroblasts across population-doubling time points. GEO provides a processed transcriptomic table suitable for downstream analysis.

**GEO accession:** GSE175533

**Primary processed resource:** `GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx`

**Expected local path:** `Data/GSE175533_hTERT.RS.RIS.CD.TPM_table.xlsx`

### H3-Control-2 — GSE179848

**Process:** replicative cellular lifespan / aging in primary human fibroblasts

**Why selected:** independent biological process, human primary fibroblasts, repeated longitudinal sampling at approximately 11-day intervals, multiple donors, and directly available processed gene-level expression matrices. This provides a stronger temporal negative control than a two-time-point aging comparison.

**GEO accession:** GSE179848

**Primary processed resource:** `GSE179848_processed_cell_lifespan_RNAseq_data.csv.gz`

**Expected local path:** `Data/GSE179848_processed_cell_lifespan_RNAseq_data.csv.gz`

## Exclusion rule

These controls are intentionally **not** matched to Yamanaka reprogramming in cell identity or perturbation. Their role is specifically to test H3: whether the same state-representation procedure generates similarly high cross-dataset temporal agreement on unrelated biological processes.

A high agreement in either unrelated control is evidence against interpreting generic trajectory agreement as process-specific.

## Analysis lock

No H3 result may be used to change this manifest. If either dataset cannot be processed under the locked representation/preprocessing rules, the primary decision remains `SIGNAL_UNRESOLVED` rather than replacing the control.
