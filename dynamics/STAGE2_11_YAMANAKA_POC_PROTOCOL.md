# Stage 2.11 — Yamanaka perturbation proof-of-concept

## Purpose

This is a **proof-of-feasibility**, not a final mechanistic model. The goal is to test whether recent human Yamanaka-factor datasets contain a reproducible perturbation signal that is strong enough to justify the PhD research programme.

## Scientific question

> Can expression data distinguish the effects of OSK from individual Yamanaka factors, and does the OSK perturbation produce a reproducible direction of change across independent datasets?

## Datasets

1. **GSE297233** — recent human bulk RNA-seq perturbation experiment. Tracked matrix: `Data/GSE297233_raw_counts_matrix.csv.gz`.
2. **GSE304042** — recent human RPE experiment containing control, individual factors and OSK conditions. Tracked matrix: `Data/GSE304042_5_ARPE_single_triple_OSK_30-1011800743.csv.gz`.
3. **GSE297234** — recent human scRNA-seq partial reprogramming time course. The large RDS files remain local and ignored by Git; they are the next temporal-validation input, not a prerequisite for this POC.

## What is tested

- conservative log-CPM normalization;
- automatic discovery of gene and numeric sample columns;
- transparent sample-condition parsing;
- factor-vs-control expression signatures;
- cross-dataset OSK signature direction concordance;
- leave-one-sample-out perturbation classification when sample counts permit it.

## What is deliberately *not* claimed

- classification accuracy is not evidence of mechanism;
- cross-dataset correlation is not causal identification;
- the POC does not infer biology from a handful of marker genes;
- the POC does not use an LLM as ground truth;
- no large common-gene harmonization is required, so the result is complementary to the Stage 2.6/Phase 0 failure.

## Success criteria for the PhD feasibility discussion

The result is useful if it establishes a reproducible perturbation signal and identifies at least one clear next mechanistic question. Stronger evidence consists of:

1. OSK vs control produces a non-trivial signature in both datasets;
2. the OSK signature has positive cross-dataset direction agreement;
3. individual-factor signatures are distinguishable from OSK;
4. results are not explained by a single extreme gene;
5. the same analysis can later be repeated on GSE297234 pseudobulk trajectories.

A negative result is also informative: it would narrow the proposed thesis toward better representation learning and context-aware models rather than claiming that the raw transcriptome directly exposes a universal mechanism.

## Run

```powershell
python validate_pipeline.py --yamanaka-poc
```

Outputs are written to:

`results/Dynamics/stage2_11_yamanaka_poc/`

- `01_input_audit.json`
- `02_perturbation_signatures.csv`
- `03_poc_summary.json`

The large GSE297234 RDS files are intentionally not committed. When temporal single-cell validation is added, it should use pseudobulk/sample-level summaries first, preserving the existing leakage-free validation philosophy.
