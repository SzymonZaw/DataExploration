# Phase 1 — Biological state representation ablation

## Scientific question

> Does replacing the raw gene-level common space with biologically informed representations improve temporal generalization across independent reprogramming datasets?

Phase 1 changes **only the state representation**. The forecasting benchmark remains the Phase 0 benchmark.

## Representations

1. **Genes** — the existing 11,899-gene common human space; Phase 0 reference.
2. **PROGENy** — pathway activity scores inferred from the harmonized gene matrix using the human PROGENy prior.
3. **DoRothEA** — TF activity scores inferred from the harmonized gene matrix using human DoRothEA confidence levels A/B/C.

The prior-knowledge networks are fixed external resources. No representation parameters are fitted on the held-out dataset.

## Forecasting protocol

For every representation:

- same trajectory datasets: GSE67462, GSE28688 and GSE297234;
- leave-one-dataset-out evaluation;
- same five random seeds;
- training-only feature selection, imputation and scaling inside the forecasting benchmark;
- approximately 60% observed prefix;
- genuine free rollout without feeding true future observations back into the model;
- identical persistence, nearest-time and linear baselines;
- identical temporal permutation null.

The representation is therefore the experimental variable.

## Prior knowledge

PROGENy and DoRothEA are retrieved through `decoupler`/OmniPath. The Phase 1 implementation uses the academic-license resources explicitly. The resulting sample-by-representation matrices are saved under:

`results/Dynamics/phase1_representation_ablation/`

## Primary evidence package

For each representation/model combination we retain:

- RMSE;
- improvement versus persistence;
- improvement versus nearest-time;
- improvement versus linear extrapolation;
- 5th and 95th percentile improvements;
- permutation p-value;
- conservative `predictive_support` flag.

A representation is considered promising only if the predictive signal is positive and robust across the predefined criteria. A low p-value alone is not sufficient.

## Interpretation

If pathway or TF representations outperform the gene representation, the result supports the hypothesis that **biological abstraction improves temporal identifiability**.

If they do not, the next candidate explanation is not automatically a more complex neural architecture. The project should then investigate single-cell representations, temporal coverage, experimental heterogeneity and identifiability limits.

This phase deliberately excludes LLM reasoning, Neural ODE/SDE and Digital Biological Twin construction.
