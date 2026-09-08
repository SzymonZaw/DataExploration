# Stage 2.11C — GSE263713 Representation Gate

## Status

Prospective technical gate committed before any GSE263713/Yamanaka similarity calculation.

## Scope

This gate may inspect only:

- input SHA-256;
- matrix orientation and sample structure;
- gene identifier namespace;
- identifier-to-HGNC mapping;
- mapping collisions;
- overlap with the locked PROGENy and DoRothEA target universes.

It must **not** calculate Yamanaka similarity, trajectory agreement, H3 null statistics, or a biological decision.

## Locked input

- Accession: `GSE263713`
- Expected SHA-256: `8ce1e16a0039d93449e70aa6e827bbc8510a9b15dfbb78dc920cc2a2de5615a6`
- Expected design: 6 individuals × 13 time points (0–48 h, every 4 h)
- Expected RNA-seq sample count: 78

NCBI GEO independently describes GSE263713 as RNA-seq from six healthy primary fibroblast cell lines sampled at 13 time points over 48 hours. urlNCBI GEO GSE263713https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263713

## Identifier gate

The full gene-ID column must be audited, not merely a prefix/sample of rows.

The canonical namespace must be resolved to HGNC symbols before network overlap. For Ensembl IDs, remove version suffixes before mapping. Duplicate canonical HGNC symbols are retained in the mapping audit and must be averaged before network scoring.

Minimum mapping rate: **90% of unique input gene IDs**.

## Representation gate

The candidate must overlap the locked human network universes at least as follows:

- PROGENy top-100: **≥80% of network targets**;
- DoRothEA A/B/C: **≥80% of network targets**.

These are technical compatibility thresholds, not biological effect thresholds.

## Sample/matrix gate

All of the following are required:

- gene × sample orientation;
- exactly 78 RNA-seq sample columns;
- six inferred individuals;
- 13 inferred time points: 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48 h;
- every individual has a complete 13-point trajectory;
- ≥95% numeric content in the first 20 expression rows after excluding gene-annotation columns.

## Decision

Only if **all** gates pass may the candidate enter the separate H3 similarity/null analysis.

A pass means **technical eligibility only**. It does not support H3, H3_NOT_SUPPORTED, process specificity, causality, or mechanism.
