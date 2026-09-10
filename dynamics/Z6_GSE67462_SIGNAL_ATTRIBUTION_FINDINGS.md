# Z6 GSE67462 signal attribution findings

## Audit status

The GSE67462 signal-attribution audit was executed on `stage-2-11c-z6-within-system-audit` with:

```powershell
python -m dynamics.run_z6_gse67462_signal_attribution
```

The audit is diagnostic and does not modify the frozen Z6 Phase 0.1 predictive-support decision.

## Dataset context

GSE67462 is a mouse secondary-reprogramming bulk-expression time series with two samples per time point. GEO describes a time-series design spanning day 0 through day 18 and a corresponding iPS-cell endpoint. The experiment used OSKM induction and expression profiling by array. These characteristics make the two branch trajectories useful for a reproducibility/context decomposition, but they do not by themselves establish biological specificity. See NCBI GEO accession GSE67462 for the experimental design and sample structure.

## Observed attribution

The audit recovered 11,899 genes across the two branches and the eight shared time points used by the Z6 representation.

| quantity | value |
|---|---:|
| median common fraction | 0.975871 |
| median branch fraction | 0.024129 |
| fraction common dynamic | 0.750903 |
| fraction common endpoint | 0.142449 |
| fraction branch context | 0.000840 |
| fraction unresolved | 0.105807 |
| endpoint observed Pearson | 0.875416 |
| endpoint permutation p | 0.000999 |

The endpoint permutation p is at the minimum attainable value with the current 1,000-permutation plus-one correction. It should therefore be reported as `p=0.000999` / permutation floor, not as an exact probability with greater numerical precision.

## Interpretation

The decomposition strongly favors a **common dynamic molecular component** over a branch-specific component within GSE67462. In particular:

- the median common fraction is ~97.6%, while the median branch fraction is ~2.4%;
- ~75.1% of genes are classified as carrying common dynamic structure;
- only ~0.084% are classified as branch-context structure;
- ~14.2% are classified as common endpoint structure;
- ~10.6% remain unresolved;
- endpoint effects are strongly concordant between branches (Pearson r=0.8754), with the permutation statistic reaching its floor.

This is materially stronger than the earlier statement that GSE67462 merely contains a reproducible branch trajectory. The positive Z6 result is **not obviously explained by branch identity alone**.

However, the result must not be promoted to biological specificity or causality. The decomposition is expression-only and associative. A low branch fraction does not prove that the common component is biological rather than a shared experimental/contextual effect. In particular, both branches belong to the same GSE67462 experiment, organism, reprogramming system, platform, and treatment regime. Shared technical and biological context therefore remains inseparable without an external perturbation, orthogonal molecular layer, or independent experiment.

## Refined Z6 interpretation

The current evidence supports the following working formulation:

> Within GSE67462, the predictive signal is predominantly carried by a reproducible common dynamic expression component rather than branch-specific structure. Its biological specificity and transferability beyond this experimental context remain unestablished.

This refines H-Z6b rather than rejecting it. Cross-system transfer still failed under the frozen Phase 0 LODO criterion, while within-system transfer was supported for AE and PCA in GSE67462. The present decomposition explains why that within-system result should not be dismissed as simple branch-specific overfitting, but it does not establish that the same state representation is portable across independent biological systems.

## What remains unresolved

The next validation target is now **orthogonal specificity**, not more architecture search. Priority tests are:

1. aggregate the common dynamic component into PROGENy/DoRothEA or other mechanistically interpretable programs;
2. compare the common component against the original GSE67462 Oct4/chromatin measurements where compatible data are available;
3. test whether the same common programs are recovered in an independent OSKM experiment;
4. use an independent perturbation/context dataset to distinguish biological state from shared treatment/time structure;
5. only then consider a stronger Z4/H3 claim.

No model threshold, architecture, dataset inclusion rule, or predictive-support criterion should be changed based on this result.
