# Z4 — GSE242421 orthogonal validation preflight

## Decision

Do **not** download `GSE242421_RAW.tar` (17.3 GB) at this stage.

NCBI reports nine scATAC samples covering D0, D2, D4, D6, D8, D10, D12, D14 and iPSC, processed on hg38 with Chromap, with fragment BED files supplied as processed supplementary data. The GEO record explicitly states that raw data are not provided. citeturn0search0turn0search1

The broader GSE242424 SuperSeries additionally contains paired scRNA and D1/D2 multiome data. citeturn0search3

## Preferred next source

Use the study's published **analysis products** rather than the 17.3 GB fragment archive. The paper states that the analysis products include integrated scATAC/scRNA/multiome count matrices, cell representations, cluster assignments, and peak calls. citeturn0search6

The study's project browser also exposes mapped products and analysis products. citeturn0search7

## Scientific consequence

This makes the candidate substantially more practical than initially expected: the orthogonal assay has a dense human fibroblast reprogramming time course and an independently generated chromatin layer, while avoiding a multi-gigabyte raw-fragment workflow.

The next gate is therefore **analysis-product provenance and format audit**, not data download.
