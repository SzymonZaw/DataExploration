# Stage 2.11C — H3 Control Candidate Audit

Status: **prospectively locked candidate selection; no H3 similarity has been computed**

## Selection rule

Candidates were assessed using only metadata/design properties, before any trajectory similarity or representation score was computed. The primary negative-control class must be temporally resolved, biologically unrelated to reprogramming, sufficiently replicated, and technically compatible with the representation contract.

## Candidate table

| Accession | Process class | Organism/system | Temporal design | Replication | Measurement | Decision | Reason |
|---|---|---|---|---|---|---|---|
| GSE3945 | acute serum-response / wound-healing program | human fibroblasts | 18 samples across time points up to 36 h | series-level repeated time-course design | expression microarray, log2 ratios | **FINAL H3-Control-1** | Acute stimulation response rather than reprogramming/aging; temporally resolved; processed GEO data available; sufficiently orthogonal biological process for H3 |
| GSE129486 | acute inflammatory cytokine response | human synovial fibroblasts | 0, 2, 4, 6, 8, 10, 12, 18, 24 h | multiple donors / cell-level time-course; donor/block structure available in metadata | RNA-seq / TPM, SMART-seq2 | **FINAL H3-Control-2** | Cytokine-response process distinct from pluripotency induction and replicative aging; dense temporal sampling; processed transcriptomic data available |
| GSE241132 | human wound healing | human skin wound tissue | intact, day 1, day 7, day 30 | multiple wound stages; donor structure | single-cell + spatial transcriptomics | **EXCLUDED** | Strong biological candidate, but raw sequencing is controlled-access and the historical analysis contract is sample/pseudobulk expression; would require a separate technical ingestion design |
| GSE56 | serum stimulation of fibroblasts | human fibroblasts | 0, 0.5, 2, 4 h | 4 samples | expression microarray | **EXCLUDED** | Temporal resolution is too coarse for the locked H3 null/decision design and would make attainable inference weak |
| GSE14256 | PDGF-BB / b-FGF stimulation | human fibroblasts | 1 h and 24 h | two experiments | expression microarray | **EXCLUDED** | Only two temporal points; insufficient for the declared temporal null procedure |
| GSE175533 | replicative senescence / cellular lifespan | human WI-38 fibroblasts | PDL time course | replicated | TPM | **RETIRED FROM PRIMARY H3** | Biologically related to cell-state regulation and aging/reprogramming; similarity would be ambiguous between generic structure and genuine biological overlap |
| GSE179848 | replicative cellular lifespan / aging | primary human fibroblasts | approximately every 11 days, ~7 time points/donor | 4 healthy + 3 SURF1 donors | RNA-seq | **RETIRED FROM PRIMARY H3** | Same ambiguity class as GSE175533; retained only as historical provenance |

## Final prospective control set

The prospective H3 control set is:

1. **GSE3945 — acute serum-response / wound-healing program**
2. **GSE129486 — acute inflammatory cytokine response in synovial fibroblasts**

These represent two distinct acute-response process classes and are not replicative-aging datasets.

## Technical gate before H3 execution

The final selection is biological/prospective, not a claim that ingestion has already succeeded. Before any H3 statistic is computed, both datasets must pass the technical gate in `STAGE_2_11C_H3_CONTROL_REDESIGN_PROTOCOL.md`:

- exact source-file SHA-256 persisted;
- matrix orientation verified;
- time labels and biological replicate/block structure verified against source metadata;
- gene/probe identifier namespace recorded;
- mapping/collision statistics recorded;
- representation gene overlap recorded;
- missing/zero handling recorded;
- same representation/preprocessing contract as Yamanaka applied;
- no post-selection replacement based on observed H3 similarity.

If either final control fails this gate, the H3 decision is `H3_UNRESOLVED`; do not substitute another dataset after looking at similarity results. A new control design would require a new prospective manifest.

## External metadata provenance

- GSE3945 is the human fibroblast serum-response time course, with 18 samples and measurements through 36 h. NCBI GEO describes it as a fibroblast serum-response time course and relates the program to wound healing.
- GSE129486 contains human synovial-fibroblast time-course transcriptomics after TNF / TNF+IL-17A stimulation at 0, 2, 4, 6, 8, 10, 12, 18 and 24 h, with donor information and processed transcriptomic data.

No H3 similarity result has been inspected or used for this selection.
