# GSE67462 temporal-module specificity findings

## Status

Diagnostic mechanistic-interpretation audit downstream of the frozen Z6 predictive analysis. It does not modify the Z6 predictive-support criterion, multimodal-support decision, temporal module assignments, or clustering procedure.

The audit uses the provenance-validated GSE67462 expression universe (11,048 genes), primary temporal modules M1, M3, M4 and M6 (>=20 genes), and versioned mouse GO-BP, Reactome and Hallmark gene sets.

## Main scientific conclusion

The temporal modules are **not adequately described by a single generic pluripotency trajectory**. The strongest evidence supports a structured transition from early somatic/mesenchymal remodeling toward later epithelialization-compatible and OCT4-associated regulatory activity, while retaining substantial context-dependent stress, proliferation and metabolic signals.

| Module | Shape | Primary interpretation | Qualification |
|---|---|---|---|
| M1 | early transient | ECM/mesenchymal remodeling and receptor/growth-factor signaling | strong EMT/mesenchymal and stress/inflammation context; weak pluripotency signal |
| M3 | mixed | exploratory intermediate program | pathway evidence is weaker and less coherent |
| M4 | late rising | ECM remodeling, trafficking and cytoskeletal/RHO-RAC-RAB signaling | stress/proliferation are secondary context |
| M6 | late rising | epithelialization-compatible, OCT4-associated late transition | keratinization/cornified-envelope and metabolic signals make a specific pluripotency claim premature |

## M1: early remodeling

M1 contains 196 genes. GO/Reactome enrichment is dominated by extracellular-matrix organization, vasculature/development, growth-factor and receptor signaling, cell adhesion and wound/tissue remodeling. Hallmark EMT is exceptionally strong (global FDR ~6.36e-08). Leading multimodally coherent genes include **EMILIN1, WNT5A, IL6, DDR2, MMP2, DCN, EFEMP2, ANGPT2, FGF2 and CAV1**.

This is consistent with substantial remodeling of the starting somatic/mesenchymal state. It should **not** be labelled proof of EMT or MET direction: enrichment is associative and the module was defined from temporal expression geometry rather than a causal state label.

## M4: late structural remodeling

M4 contains 258 genes and is late-rising. Reactome terms include **RAB regulation of trafficking, extracellular-matrix organization, RAC1/RHO GTPase cycles, integrin interactions and receptor tyrosine kinase signaling**. Leading genes include **SPP1, EPHA2, COL6A1, PAK1, CDC42EP1, COL5A3, PLD2 and ITGA5**.

The defensible interpretation is late structural/signaling and trafficking reorganization, not a causal pathway assignment.

## M6: epithelialization-compatible OCT4-associated transition

M6 is the largest module (497 genes; ~44% of the multimodal core) and is late-rising. It has strong OCT4/RNAPII association and enrichment for epithelialization-related Hallmark terms plus Reactome **formation of the cornified envelope** and **keratinization**. A smaller pluripotency/stem-cell signal includes highly coherent candidates such as **NODAL, SPINT1, LEF1, KRT8, SOX15, BMP7 and PDGFB**.

Therefore M6 is best described as **compatible with epithelialization and an OCT4-associated late transition**, rather than as a pure pluripotency module. The keratinization/cornified-envelope and metabolic signals indicate broader cell-state remodeling.

## M3: exploratory

M3 contains 170 genes and has a mixed trajectory. No specificity category reaches the same global-FDR support as the other major modules. It should remain exploratory rather than receive a mechanistic label.

## Specificity audit

The keyword-based category screen returned 68 category-matched significant terms and 409 leading gene-category pairs.

- **M1:** ECM/mesenchymal remodeling dominates; stress/inflammation, proliferation and metabolism are meaningful secondary signals. Epithelialization is present but weaker; pluripotency/stem-cell terms are sparse.
- **M4:** ECM/mesenchymal remodeling remains prominent, with epithelialization and stress/proliferation as secondary context.
- **M6:** epithelialization dominates; a smaller pluripotency/stem-cell component and metabolic context remain.
- **M3:** no category-specific enrichment passes the same global-FDR threshold.

This pattern is more consistent with a **multi-process state transition** than with a one-dimensional pluripotency score.

## Important qualifiers

The Hallmark nuisance screen is scientifically informative and was not filtered away. Generic stress/inflammation, proliferation/cell-cycle and metabolism occur in several modules. Therefore the recovered multimodal signal cannot be described as purely lineage-specific or pluripotency-specific.

The evidence does **not** establish that any enriched pathway is causal, that M1/M4 are definitively EMT/MET states, that M6 is a pure pluripotency state, or that these modules are universal or transferable across independent reprogramming systems.

## Relation to multimodal evidence

This interpretation is downstream of the multimodal validation in which **H3K27ac, H3K4me3, RNAPII and total OCT4** passed the frozen support criterion. **H3K27me3** remained a negative/contrast modality and did not satisfy support. Thus the mechanistic interpretation is not expression-only, although enrichment itself remains associative and does not prove direct enhancer-to-gene regulation.

## Thesis-level wording

> In the GSE67462 secondary OSKM reprogramming system, the provenance-validated multimodal state representation resolves temporally stable modules whose functional annotations suggest a transition from early ECM/mesenchymal and growth-factor-associated remodeling toward later epithelialization-compatible, OCT4-associated regulatory activity. The modules retain substantial context-dependent stress, proliferation and metabolic signals, indicating that the recovered state representation is biologically structured but not reducible to a single universal pluripotency axis.

## Next scientific step

The appropriate next step is **mechanistic falsification rather than another enrichment pass**: test whether the M1 ECM/mesenchymal program, M4 RHO/RAC/RAB-integrin remodeling and M6 epithelialization/OCT4-associated program preserve their ordering in an independent system or perturbation context, and whether these candidates outperform generic proliferation/stress markers in predicting the transition.

## Provenance

- validated background: 11,048 genes
- primary modules: M1, M3, M4, M6
- minimum module size: 20 genes
- minimum term overlap: 5 genes
- test: one-sided hypergeometric ORA
- multiple testing: Benjamini-Hochberg FDR
- gene sets: MSigDB 2026.1.Mm mouse symbols (GO-BP, Reactome, Hallmark)
- frozen Z6 predictive-support decision unchanged

NCBI GEO describes GSE67462 as a mouse secondary OSKM reprogramming time series with two replicates per timed point and integrated expression, Oct4-binding and histone-mark profiling; its summary likewise describes regression of the somatic program followed by gradual acquisition of a pluripotent program. citeturn0search0