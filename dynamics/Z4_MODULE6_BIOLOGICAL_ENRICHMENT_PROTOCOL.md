# Z4 Module 6 Biological Enrichment Protocol

## Purpose

Characterize the biological content of frozen Z6 module 6 after the GSE242421 orthogonal test identified module 6 as the most stable candidate across peak-to-gene link thresholds.

This is an **interpretation audit**, not a new validation gate. The Z4 decision remains `Z4_ORTHO_THRESHOLD_SENSITIVE` regardless of enrichment results.

## Frozen inputs

- Z6 GSE67462 module assignments: `results/Dynamics/z6_gse67462_temporal_modules/03_gene_module_assignments.csv`
- frozen mouse-human 1:1 orthology: `results/Dynamics/z4_frozen_orthology_mapping/01_frozen_human_mouse_mapping.csv`
- GSE242421 derived gene activity: `results/Dynamics/z4_gse242421_gene_activity/01_gene_activity_log1p_cpm.tsv.gz`
- target temporal samples: D0, D2, D4, D6, D8, D10, D12, D14

## Analysis

1. Keep module 6 exactly as defined in GSE67462.
2. Map mouse genes to human only through the frozen 1:1 orthology table.
3. Retain only human genes present in the GSE242421 gene-activity matrix.
4. Perform GO Biological Process and Reactome enrichment with g:Profiler.
5. Use the frozen orthology universe as the enrichment background rather than the whole genome.
6. Rank validated module-6 genes by Spearman correlation with D0-D14 target gene activity.
7. Define an exploratory target-side temporal core as the top 25% by signed Spearman correlation. This is a descriptive subset, not a refitted module.
8. Run the same enrichment analysis on the target-side temporal core.

## Interpretation rules

- Enrichment identifies biological processes associated with the gene set; it does not demonstrate pathway activation or causal mechanism.
- Leading/core genes should be treated as hypotheses for mechanism, not independently validated causal genes.
- Because GSE242421 has one sample per time point, target-side enrichment cannot establish biological replication by itself.
- Enrichment must not be used to alter the frozen module definition or the Z4 statistical decision.

GSEA documentation describes the leading-edge subset as the genes contributing most to an enrichment score and recommends examining those genes when interpreting enriched sets. urlGSEA User Guidehttps://docs.gsea-msigdb.org/GSEA/GSEA_User_Guide/

## Outputs

- `01_module6_target_gene_ranking.csv`
- `02_module6_frozen_orthology.csv`
- `03_module6_core_temporal_genes.csv`
- `04_mouse_frozen_module_enrichment.csv`
- `04_human_validated_module_enrichment.csv`
- `04_human_target_core_enrichment.csv`
- `05_summary.json`

## Decision

This analysis does not change `Z4_ORTHO_THRESHOLD_SENSITIVE`. Its purpose is to determine whether the stable module-6 signal has coherent biological interpretation and whether the same broad processes are visible after orthogonal transfer.
