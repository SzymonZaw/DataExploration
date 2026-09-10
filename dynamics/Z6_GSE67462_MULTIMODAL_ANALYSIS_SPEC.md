# Z6 GSE67462/GSE67520 multimodal analysis specification

## Scope

This specification operationalizes the existing multimodal validation protocol without changing frozen Z6 predictive support.

GSE67520 processed peak calls are analyzed in mm9 coordinates. The primary peak-to-gene representation is promoter-associated burden: peak midpoint within +/-2 kb of one representative transcript TSS per gene, with non-negative peak scores summed per gene and time point. Distal peaks are not assigned to nearest genes in the primary analysis.

The expression reference is the day-0-centered common branch trajectory already defined for GSE67462. Primary temporal validation uses day 0, 1, 3, 5, 7, 11, 15 and 18. GSE67520 iPSC files are retained in provenance but excluded from the primary statistic because the current expression trajectory does not provide a validated iPSC point.

## Annotation

The runner defaults to `Data/GSE67520/mm9.refGene.gtf.gz`. If absent, it downloads the UCSC mm9 RefSeq GTF to that path. A local annotation can be supplied with `--gtf`. The annotation source is recorded in the output JSON.

## Statistics

For each of six modalities (Oct4, H3K4me1, H3K27ac, H3K4me3, H3K27me3, RNAPII):

1. parse all processed peak calls;
2. recover sample time from the filename;
3. map peak midpoint to an overlapping promoter;
4. sum peak scores per gene/time point;
5. intersect with the GSE67462 common expression gene space;
6. calculate per-gene temporal Spearman concordance;
7. calculate a global gene-time Spearman effect;
8. create a time-label permutation null with the gene values unchanged;
9. report observed effect, null mean, null 95th percentile and empirical one-sided p-value.

A modality is supported only when its global observed concordance is above the permutation 95th percentile and empirical p < 0.05. The overall label is `MULTIMODAL_DYNAMIC_SUPPORT` only when at least two distinct modalities satisfy both criteria.

## Interpretation constraints

This is mechanistic/multimodal coherence within one experimental reprogramming program. It is not independent-system validation and does not establish biological specificity or cross-system transferability. A positive result cannot alter the frozen Z6 predictive-support decision.
