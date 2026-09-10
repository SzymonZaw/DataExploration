# Z4 GSE297234 RDS structure audit protocol

## Purpose

Inspect the two locally available GSE297234 Seurat RDS objects before constructing the frozen human→mouse mapping for external specificity validation.

This is a **diagnostic pre-transfer step**. It must not compute a Z4 transfer score, fit a model, select features using validation results, or define the final biological endpoint.

## Questions

For each RDS:

1. Is the object a Seurat object?
2. Which assays are present and which assay is default?
3. What are the cell and feature dimensions?
4. What are the actual feature identifiers stored in each assay?
5. Which metadata fields define donor/age, time, treatment, cell state and nuisance/context?
6. Are the two RDS objects compatible with a deterministic common feature namespace?

## Required outputs

The audit writes one JSON manifest containing:

- object class and dimensions;
- assay names and default assay;
- representative feature identifiers for every assay;
- all metadata column names;
- selected metadata fields relevant to donor, time, intervention, cell state and nuisance context.

## Scientific boundary

GSE297234 is an independent human partial-reprogramming experiment: NCBI GEO describes young and aged human fibroblasts treated with Sendai-virus OSKM at days 0, 3, 7 and 10. The study uses 10x scRNA-seq and GRCh38. These facts establish external commensurability as plausible, but do not establish specificity by themselves.

The downstream Z4 decision still requires a locked human→mouse mapping, an independently defined biological endpoint, positive target transfer, and a specificity/context challenge.

## Reproducibility

Run from the repository root with R 4.6.x and SeuratObject available:

```powershell
Rscript dynamics/run_z4_gse297234_rds_structure_audit.R `
  --rds Data\GSE297234_GM00731_SEVOSKM.rds Data\GSE297234_HFIB_COMBINED_SEVOSKM.rds `
  --output results\Dynamics\z4_gse297234_rds_structure_audit\structure.json
```

The script loads one RDS at a time and performs no expression-level modelling.
