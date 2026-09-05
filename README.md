# Mathematical Modeling and Prediction of Cellular Reprogramming Dynamics

## Overview

This repository contains the computational part of a research project focused on reconstructing and predicting the dynamic state of a biological cell from heterogeneous, incomplete and independently generated biological data.

The central scientific problem is not the construction of a Digital Biological Twin itself. The primary objective is to determine whether independent multimodal experiments can be integrated into a biologically meaningful state space, whether cellular trajectories can be reconstructed in that space, and whether their dynamics can be modeled and predicted for data that were not used to build the model.

The Digital Biological Twin is treated as a downstream demonstrator of the resulting methodology.

---

## Proposed PhD topic

> **Mathematical modeling and prediction of cellular reprogramming dynamics based on heterogeneous multimodal biological data.**

An alternative broader formulation is:

> **Reconstruction and prediction of the dynamic cellular state based on heterogeneous multimodal data using machine learning, state-space geometry and dynamical-system modeling.**

The main biological process used as a research case is **OSKM-mediated cellular reprogramming**, including transitions from fibroblast-like states toward pluripotency.

---

## Main research question

> **Can a reliable dynamic trajectory of cellular state be reconstructed from independent, heterogeneous and incomplete multimodal biological datasets, and can this trajectory be used to predict future cellular states?**

A more ambitious methodological question is:

> **Can a common biological state space and an interpretable mathematical model of reprogramming dynamics be identified from heterogeneous experiments and generalize to independent datasets?**

---

## Main objective

The main objective is to develop and validate a computational methodology for:

1. harmonizing heterogeneous biological measurements,
2. constructing a biologically meaningful common cellular state space,
3. reconstructing longitudinal cellular trajectories,
4. modeling the dynamics of state transitions,
5. identifying changes in stability and possible transition points,
6. discovering interpretable mathematical relationships governing the dynamics,
7. validating predictions on observations and entire datasets excluded from model construction,
8. representing the resulting state, trajectory, prediction and uncertainty in a Digital Biological Twin.

The key requirement is **generalization**, not merely obtaining a visually convincing trajectory or a low in-sample reconstruction error.

---

## Research hypothesis

> **Integration of independent multimodal biological measurements in a biologically defined state space can produce stable and interpretable cellular trajectories and enable prediction of future cellular states, provided that heterogeneity between experiments, uncertainty and missing observations are explicitly modeled.**

A stronger hypothesis to be investigated is:

> **Cellular reprogramming dynamics can be approximated by a low-dimensional, mathematically interpretable dynamical model whose predictive ability can be validated on independent experiments.**

These hypotheses are empirical and must be tested against appropriate baselines and held-out data.

---

# Specific research tasks

## Task 1 — Integration and characterization of heterogeneous biological data

Develop a reproducible data-processing framework for datasets representing different biological modalities and experimental designs, including:

- bulk RNA-seq,
- microarray expression data,
- single-cell RNA-seq,
- ChIP-seq / regulatory measurements,
- longitudinal and time-course experiments.

The methodology must explicitly account for:

- different platforms,
- different organisms,
- different measurement resolutions,
- different sampling times,
- biological and technical replicates,
- incomplete observations,
- experimental perturbations and branches.

**Expected result:** a standardized and traceable representation of the input datasets and their metadata.

---

## Task 2 — Construction of a common biological state space

This is the central methodological task.

The current PCA + Procrustes approach is treated as an experimental diagnostic, not as the final biological representation. A proper common state space should be based on biological features that can be meaningfully compared between studies.

Potential components include:

- shared genes,
- human–mouse orthologs,
- gene modules and pathways,
- regulatory features,
- pseudobulk representations of single-cell data,
- platform-aware normalization,
- batch-effect correction,
- latent-factor models and other dimensionality-reduction methods.

The state of a biological system will be represented as a vector such as:

\[
z(t) = [z_1(t), z_2(t), \ldots, z_n(t)]
\]

where the coordinates have a reproducible biological or latent interpretation.

**Expected result:** a common representation that can be tested independently rather than imposed through time alignment alone.

---

## Task 3 — Reconstruction of dynamic cellular trajectories

Reconstruct trajectories of cellular state:

\[
z(t_0) \rightarrow z(t_1) \rightarrow z(t_2) \rightarrow \cdots \rightarrow z(t_n)
\]

The analysis will quantify:

- direction of state change,
- rate of change,
- trajectory length,
- trajectory similarity,
- branching,
- divergence between experimental conditions,
- stability and variability of replicate measurements.

Special attention will be paid to the distinction between genuine biological trajectories and trajectories introduced by preprocessing or alignment methods.

**Expected result:** quantitative representations of reprogramming trajectories together with uncertainty and replicate variability.

---

## Task 4 — Mathematical modeling of cellular-state dynamics

The dynamics may be formulated as a state-space or dynamical-system problem, for example:

\[
\frac{dz}{dt} = F(z,t,u,\theta)
\]

where:

- \(z\) is the cellular state,
- \(t\) is biological time,
- \(u\) represents experimental perturbations such as OSKM,
- \(\theta\) represents model parameters,
- \(F\) describes the state transition dynamics.

Candidate model families may include:

- ordinary differential equations,
- state-space models,
- probabilistic dynamical models,
- Neural ODEs,
- geometric models of cellular state space,
- other appropriate machine-learning dynamical models.

The model family should be selected empirically rather than assumed in advance.

**Expected result:** a quantitative model capable of representing cellular-state transitions and producing testable predictions.

---

## Task 5 — Identification of state transitions and changes in stability

Investigate whether reprogramming contains mathematically detectable changes in dynamical stability before the appearance of classical pluripotency markers.

Potential indicators include:

- increasing variance,
- increasing autocorrelation,
- critical slowing down,
- changes in local velocity,
- changes in the Jacobian,
- local geometric curvature,
- changes in potential-like landscapes,
- bifurcation-like behavior.

A key biological question is:

> **Can an approaching transition toward pluripotency be detected before the cell reaches a classical pluripotent state?**

**Expected result:** quantitative indicators of state transitions and their predictive value.

---

## Task 6 — Mathematical model discovery and symbolic regression

Once a validated common state representation and reliable trajectories have been obtained, symbolic regression / genetic programming can be used to search for interpretable dynamical relationships.

For example:

\[
\frac{dz_1}{dt} = f_1(z_1,z_2,\ldots,z_n)
\]

\[
\frac{dz_2}{dt} = f_2(z_1,z_2,\ldots,z_n)
\]

The goal is not simply to maximize predictive accuracy. The analysis should investigate whether relatively simple equations can explain the observed dynamics and generalize beyond the datasets used for model discovery.

Symbolic models will be treated as hypotheses about the system and will require independent validation. They will not automatically be interpreted as causal mechanisms.

**Expected result:** interpretable candidate equations describing aspects of cellular-state dynamics.

---

## Task 7 — Out-of-sample and cross-dataset validation

Generalization is a central criterion of the project.

Validation strategies will include:

- leave-one-replicate-out validation,
- leave-one-timepoint-out validation,
- leave-one-dataset-out validation,
- prediction on an entire independent experiment,
- comparison against naive and conventional baselines.

The fundamental test is:

\[
D_1 + D_2 + D_3 \rightarrow \text{model}
\]

followed by:

\[
\text{model} \rightarrow D_4
\]

where \(D_4\) was not used during model construction.

A low in-sample alignment or reconstruction error will not be considered sufficient evidence of validity.

**Expected result:** quantitative evidence for, or against, cross-study generalization.

---

## Task 8 — Digital Biological Twin demonstrator

The final methodology will be integrated into the existing Digital Biological Twin prototype (`testHP`).

The Twin is intended to demonstrate how the scientific methodology can represent:

```text
OBSERVED
    ↓
MEASURED
    ↓
DERIVED
    ↓
INFERRED
    ↓
PREDICTED STATE
```

Each inferred or predicted element should be associated with:

- uncertainty,
- provenance,
- temporal context,
- supporting observations,
- confidence/evidence level.

The Twin therefore serves as an application and demonstrator of the research methodology rather than being the primary scientific contribution.

---

# Current computational strategy

The project currently follows this conceptual pipeline:

```text
HETEROGENEOUS BIOLOGICAL DATA
RNA-seq / microarray / scRNA-seq / ChIP-seq / time-course
                         │
                         ▼
              DATA HARMONIZATION
       genes / orthologs / modules / batches
                         │
                         ▼
             COMMON STATE SPACE
                       z(t)
                         │
                         ▼
             TRAJECTORY RECONSTRUCTION
                         │
                         ▼
                GEOMETRY + DYNAMICS
                         │
                         ▼
               DYNAMICAL MODEL
                    dz/dt = F(z)
                         │
               ┌─────────┴─────────┐
               ▼                   ▼
          PREDICTION         MODEL DISCOVERY
                               symbolic regression
               │                   │
               └─────────┬─────────┘
                         ▼
                  OUT-OF-SAMPLE
                    VALIDATION
                         │
                         ▼
              INDEPENDENT DATASETS
                         │
                         ▼
             DIGITAL BIOLOGICAL TWIN
```

---

# Current datasets

The repository currently contains exploratory pipelines for complementary datasets related to cellular reprogramming and cell-state regulation.

| Dataset | Organism | Modality | Main role |
|---|---|---|---|
| GSE148158 | Human | RNA-seq | Reprogramming / OSKM and GFP controls |
| GSE28688 | Human | Microarray | Early reprogramming time course |
| GSE52052 | Human | Agilent microarray | Day-11 reprogramming conditions |
| GSE67462 | Mouse | Affymetrix microarray | Detailed OSKM reprogramming time course |
| GSE67520 | Mouse | ChIP-seq | Regulatory / chromatin-state information |
| GSE297234 | Human | scRNA-seq | Single-cell reprogramming trajectories |

The datasets differ substantially in organism, platform, modality, sampling schedule and experimental design. This heterogeneity is not treated as an inconvenience to hide; it is one of the central methodological challenges of the project.

---

# Current computational components

## Data exploration

The repository contains dedicated scripts for dataset-level exploratory analysis, including preprocessing, quality control, correlations, variance analysis, PCA and visualization.

Current main scripts include:

- `GSE148158.py`
- `GSE28688.py`
- `GSE52052.py`
- `GSE67462.py`
- `GSE67520.py`
- `GSE297234.py`
- `Summary.py`

## Dynamics analysis

`Dynamics.py` contains the current experimental framework for trajectory reconstruction, common-space diagnostics, validation and subsequent dynamic-analysis stages.

The development has proceeded through several validation stages rather than assuming that a common latent space is valid from the beginning.

### Stage 2.1 — Time-anchored common trajectory geometry

An initial common trajectory representation was constructed from time-aware PCA trajectories using interpolation and Procrustes alignment.

This stage demonstrated that trajectories from some datasets can be geometrically aligned, but it does **not** prove that their coordinates represent the same biological variables.

### Stage 2.2 — Common-space diagnostics

The framework evaluates:

- pairwise trajectory correlation,
- aligned RMSE,
- path-length ratios,
- leave-one-dataset-out comparisons,
- cross-dataset dispersion,
- reference sensitivity.

### Stage 2.3 — Within-time residual validation

Stage 2.3 applies transformations to actual sample-level latent points and retains replicate variability.

This exposed an important confounder in GSE148158: GFP and OSKM observations at the same time cannot simply be treated as biological replicates of one trajectory. They are represented as separate branches.

This stage is therefore an example of the project principle that validation should be used to discover failure modes rather than only to confirm a desired model.

### Stage 2.4 — Out-of-sample validation

Stage 2.4 performs leave-one-timepoint-out validation and compares aligned prediction error with a naive unaligned baseline.

The current results show partial generalization, but the behavior is heterogeneous across dataset pairs. Therefore the present PCA + Procrustes representation is **not yet accepted as the final biological common state space**.

### Stage 2.5 — Planned next step

The next methodological step is construction of a feature-level biological common representation based on comparable biological quantities rather than time-anchored PCA alignment alone.

This is expected to include:

1. shared gene space,
2. human–mouse ortholog mapping where appropriate,
3. scRNA-seq pseudobulk or carefully defined sample-level summaries,
4. platform-aware preprocessing,
5. batch-aware integration,
6. common latent-factor modeling,
7. renewed out-of-sample validation.

Only after this representation passes appropriate validation should the project rely on symbolic regression and higher-level dynamic inference.

---

# Important methodological principles

## 1. Alignment is not biological validation

Two trajectories can become highly similar after mathematical alignment even when the underlying biological states are not equivalent.

Therefore:

> **High post-alignment correlation is a diagnostic result, not evidence of biological equivalence.**

## 2. Time must not define the state by itself

A representation constructed directly from the expected time trajectory can artificially make samples at the same time appear similar.

Actual sample-level variation must be retained and independently evaluated.

## 3. Replicates are not interchangeable conditions

Different perturbations, controls and treatment branches must not be treated as technical or biological replicates merely because they were measured at the same time.

## 4. Cross-dataset prediction is more important than visual similarity

A method that reconstructs the training datasets beautifully but fails on an independent experiment is not sufficient for the main research objective.

## 5. Symbolic equations are hypotheses

A symbolic-regression equation is an interpretable model candidate. It does not automatically establish causality or biological mechanism.

## 6. Uncertainty must be explicit

The final methodology should distinguish measured observations from derived quantities, inferred states and predictions, together with uncertainty and provenance.

---

# Scientific contribution sought

The intended scientific contribution is **not** simply a software implementation of a Digital Biological Twin.

The primary contribution should be a validated methodology that answers whether heterogeneous biological observations can be transformed into a common dynamic representation that:

- preserves meaningful biological variation,
- reconstructs cellular trajectories,
- identifies changes in state and stability,
- supports prediction of future states,
- generalizes to independent experiments,
- and can yield interpretable mathematical descriptions of dynamics.

A successful outcome would therefore be a chain of evidence:

```text
heterogeneous data
      ↓
biological harmonization
      ↓
common state representation
      ↓
trajectory reconstruction
      ↓
dynamical model
      ↓
prediction
      ↓
independent validation
      ↓
interpretable model / equations
      ↓
Digital Biological Twin demonstrator
```

An unsuccessful result at any stage is also scientifically informative if it identifies a limitation of the proposed representation or modeling assumptions.

---

# Long-term research direction

A central long-term question is whether cellular reprogramming has a reproducible low-dimensional dynamical structure that can be learned from independent experiments:

\[
\text{independent heterogeneous experiments}
\rightarrow
\text{shared biological state space}
\rightarrow
\text{geometry}
\rightarrow
\text{dynamics}
\rightarrow
\text{symbolic model}
\rightarrow
\text{external prediction}
\]

An especially interesting biological question is whether a mathematically detectable loss of stability precedes the transition toward pluripotency and whether this information can be used for early prediction of cell fate.

These questions remain hypotheses to be tested rather than established conclusions.

---

# Repository outputs

Each analysis writes results to a dataset-specific directory under:

```text
results/
```

The `Dynamics.py` pipeline writes its outputs under:

```text
results/Dynamics/
```

including stage-specific diagnostic and validation outputs.

The repository is intended to keep the analysis reproducible and traceable from source datasets through intermediate representations to final model diagnostics.

---

# Reproducibility

The project is developed and executed in a dedicated Python virtual environment.

Example:

```text
.venv\Scripts\python.exe Dynamics.py
```

Input datasets are expected to be available locally in the project data directories according to the individual analysis scripts.

Large raw archives and supplementary datasets should not be committed to the repository. The repository should contain code, lightweight metadata and derived results necessary for reproducibility and interpretation.

---

# Status

**Current status:** methodological development and validation.

The project has successfully established:

- reproducible exploratory analysis of multiple reprogramming datasets,
- time-aware trajectory extraction,
- branch-aware handling of experimental conditions,
- common-trajectory diagnostics,
- within-time replicate validation,
- leave-one-timepoint-out validation.

The current PCA + Procrustes common representation remains an **experimental baseline**. It has not yet been accepted as a biologically validated universal state space.

The immediate research priority is therefore to develop and validate a feature-level biological common representation before relying on downstream symbolic dynamics.
