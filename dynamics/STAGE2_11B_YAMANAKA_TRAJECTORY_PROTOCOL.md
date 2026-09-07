# Stage 2.11b — GSE297234 trajectory proof-of-feasibility

## Scientific question

Can recent human OSKM single-cell data support a reproducible low-dimensional trajectory of cellular state across days 0, 3, 7 and 10?

This is a feasibility test for the PhD hypothesis. It does **not** claim cell lineage, causality, or mechanism discovery.

## Input

- `Data/GSE297234_GM00731_SEVOSKM.rds`
- `Data/GSE297234_HFIB_COMBINED_SEVOSKM.rds`

These files remain local and are not committed to Git.

## Method

1. Read the local RDS object with `Rscript`.
2. Support Seurat and SingleCellExperiment objects.
3. Identify the metadata field containing day/time information and retain only parsable day 0/3/7/10 observations.
4. Aggregate cells to sample/metadata-group pseudobulk; this deliberately avoids claiming cell-level lineage.
5. Apply log-CPM normalization.
6. Build a descriptive PCA trajectory from the 2,000 most variable genes.
7. Save coordinates and day centroids for each input object.
8. Compare **temporal distance profiles** between datasets rather than comparing absolute PCA coordinates, because independent PCA axes are not directly identifiable.
9. Report day-component monotonicity and reproducibility.

## Interpretation

Evidence supporting feasibility should include:

- all four expected days recovered,
- coherent state change across time,
- similar ordering/temporal-distance structure across independent human fibroblast inputs,
- no dependence on a single sample or one arbitrary PCA axis.

Failure means that this particular representation/aggregation does not recover robust trajectory structure. It does not imply that reprogramming lacks dynamics.

## Next milestone

If trajectory reproducibility is positive, add pathway/TF activity and perturbation-conditioned forecasting. Only then introduce explicit candidate mechanisms and experimental discrimination.
