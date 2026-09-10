# Z4 GSE297234 Preflight Protocol

## Purpose

Prepare an independent human scRNA-seq candidate for frozen Z4 external-specificity validation without evaluating the frozen representation score.

## Gates

1. Verify local GSE297234 RDS/sample objects and sample-level expression identifiers.
2. Freeze human→mouse orthology mapping independently of validation results.
3. Quantify mapping coverage against the frozen GSE67462 discovery feature universe.
4. Define biological state/time metadata from GEO annotations, without using expression-derived outcomes.
5. Identify an endpoint that is independent of the transferred representation score.
6. Record donor, batch, platform, treatment, and time as nuisance/context variables.

## Biological commensurability

GSE297234 is treated as an external target only if its OSKM-induced transition can be operationally related to the discovery transition without requiring result-driven remapping. The GEO sample metadata document human fibroblasts treated with Sendai OSKM at D0/D3/D7/D10 and processing against GRCh38. The preflight therefore records, but does not infer, biological state from the transferred score.

## Orthology policy

Use a deterministic one-to-one human→mouse orthology table. Prefer NCBI/Ensembl-supported ortholog relationships and retain only unambiguous mappings. Do not select mappings because they improve transfer performance.

## Endpoint independence

A candidate endpoint must be defined before frozen transfer and must not be constructed from the transferred state score or from genes selected using transfer results. If no independent endpoint can be established, Z4 remains unresolved even if frozen transfer is positive.

## Frozen boundary

This preflight must not train, refit, cluster, select features, optimize thresholds, or compute the Z4 transfer score. A successful preflight only authorizes the next stage.
