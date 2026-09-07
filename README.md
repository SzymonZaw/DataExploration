# AI-driven discovery of dynamical mechanisms of cellular state transitions

## Research project

This repository contains the computational research framework for a PhD project on **discovering and validating dynamical mechanisms of cellular state transitions from heterogeneous biological data**.

Cellular reprogramming is the primary biological case study, but the methodological goal is broader: to understand how a biological system moves between states, which components of that transition are conserved across experiments, which depend on context or perturbation, and whether the learned representation can support prediction and mechanistic hypothesis testing.

The project combines heterogeneous multi-omics data, machine learning and representation learning, biological knowledge and regulatory information, mathematical dynamical systems, perturbation-response modeling, uncertainty-aware prediction, symbolic/model discovery, and AI/LLM-assisted mechanistic reasoning.

The **Digital Biological Twin** is a downstream demonstrator of the methodology, not the primary scientific claim.

---

## PhD topic

### Primary formulation

> **AI-driven discovery of conserved and context-dependent dynamical mechanisms of cellular reprogramming from heterogeneous multi-omics data.**

### Broader methodological formulation

> **AI-driven discovery of dynamical mechanisms of cell-state transitions from heterogeneous biological data.**

The second formulation is intentionally broader. It allows cellular reprogramming to remain the main experimental case while making the scientific contribution applicable to other biological state transitions, perturbation responses and cell-fate decisions.

---

## Main scientific problem

Modern biological experiments observe cellular state through different technologies, organisms, laboratories, sampling schedules and perturbations. These datasets are therefore heterogeneous, incomplete and often incompatible at the level of individual measurements.

A central question is whether these observations can be transformed into a **biologically meaningful dynamical representation** in which we can distinguish:

1. conserved components of cellular-state dynamics,
2. context- and experiment-specific components,
3. effects of perturbations,
4. effects of cellular history or memory,
5. genuine predictive information about future state.

The project therefore does **not** assume that all experiments follow one universal trajectory.

The working dynamical formulation is conceptually:

\[
x_d(t) = G_d\big(z_d(t),u_d(t),h_d(t)\big)+\epsilon_d(t)
\]

with

\[
\frac{dz_d}{dt}=F\big(z_d(t),u_d(t),h_d(t),\theta_d\big).
\]

Here \(x_d(t)\) is the observed data, \(z_d(t)\) the latent biological state, \(u_d(t)\) the experimental perturbation, \(h_d(t)\) relevant history or memory, \(\theta_d\) context-dependent parameters, and \(F\) the unknown dynamical mechanism.

One of the main scientific goals is to determine **which parts of this model are actually supported by independent data**.

---

## Central research question

> **Can artificial intelligence, biological knowledge and mathematical dynamical modeling be combined to identify a minimal, interpretable state representation and discover mechanisms that predict cellular-state transitions across heterogeneous experiments?**

Secondary questions include:

- What information is conserved between independent experiments?
- What information is specific to cell type, protocol, perturbation or environment?
- When is the current cellular state sufficient to predict its future, and when is historical information required?
- Can early molecular signals identify an approaching state transition before the final phenotype is established?
- Which genes, regulatory modules or molecular processes are most informative about the future state?
- Can competing mechanistic hypotheses be generated, quantitatively compared and falsified against independent observations?
- Can AI identify the most informative perturbation or sampling time for distinguishing competing mechanisms?

---

## Main objective

> **To develop and validate an AI-assisted mathematical framework for discovering dynamical mechanisms of cellular-state transitions by integrating heterogeneous biological measurements, biological knowledge and perturbation information, while explicitly modeling uncertainty, context dependence and cellular history.**

The framework should move from observations toward increasingly mechanistic representations:

```text
HETEROGENEOUS BIOLOGICAL DATA
            ↓
   BIOLOGICAL HARMONIZATION
            ↓
   AI / KNOWLEDGE-GUIDED
      REPRESENTATION
            ↓
       LATENT STATE z
            ↓
     STATE TRANSITIONS
            ↓
    DYNAMICAL MODEL F
            ↓
  COMPETING MECHANISTIC
       HYPOTHESES
            ↓
 PREDICTION + UNCERTAINTY
            ↓
       FALSIFICATION
            ↓
 INFORMATIVE PERTURBATION /
       TIME POINT
            ↓
 DIGITAL BIOLOGICAL TWIN
```

The Digital Biological Twin represents the final application layer of this chain.

---

# Scientific objectives

## Objective 1 — Construct a biologically meaningful state representation

Develop methods that transform heterogeneous measurements into a representation of cellular state:

\[
x(t) \rightarrow z(t).
\]

Potential inputs include bulk RNA-seq, microarray expression, scRNA-seq, ATAC-seq, ChIP-seq/CUT&Tag, histone marks, regulatory features, proteomics and perturbation metadata.

Candidate approaches include latent-factor models, variational models, contrastive learning, foundation-model embeddings, knowledge-guided representation learning, multimodal integration and module/pathway representations.

The representation must be evaluated by **out-of-sample biological prediction**, not only by reconstruction quality or visualization.

---

## Objective 2 — Separate conserved and context-dependent dynamics

Instead of assuming one universal trajectory, investigate a decomposition such as:

\[
F_d = F_{shared} + F_{context,d}.
\]

Determine which dynamical features are reproducible across independent experiments, cell types, organisms, protocols, perturbations and measurement modalities.

A central criterion is whether a component learned from some datasets remains informative for a completely held-out dataset.

---

## Objective 3 — Identify conserved biological transition modules

Determine whether small biological modules exhibit reproducible temporal behavior even when global trajectories differ.

Candidate processes include:

- epithelial–mesenchymal transition / mesenchymal–epithelial transition,
- proliferation and MYC-associated programs,
- metabolic remodeling,
- chromatin and enhancer remodeling,
- pluripotency networks,
- stress responses,
- extracellular matrix and adhesion,
- growth-factor signaling.

Quantify activation onset, peak time, rate of change, ordering, interactions, persistence and reproducibility across experiments.

Biological modules are treated as hypotheses or priors, not as predefined proof of mechanism.

---

## Objective 4 — Model perturbation-dependent dynamics

Extend the state representation from \(z(t)\) to \(z(t,u)\), where \(u\) represents a perturbation such as transcription-factor induction, drug treatment, environmental change or genetic intervention.

The goal is to determine whether AI can predict how perturbation changes the future distribution or trajectory of cellular states.

This provides a bridge between observational trajectory reconstruction and mechanistic experimentation.

---

## Objective 5 — Determine when cellular history matters

Test whether cellular dynamics can be modeled as approximately Markovian:

\[
z_{t+1}=F(z_t,u_t),
\]

or whether prediction requires historical information:

\[
z_{t+1}=F(z_t,z_{t-1},\ldots,u_t).
\]

Potential sources of memory include chromatin state, epigenetic memory, persistent regulatory programs, metabolic state and previous perturbations.

This creates a mathematically explicit framework for studying when a current-state representation is sufficient and when a state-plus-history representation is necessary.

---

## Objective 6 — Discover and falsify mechanistic models

AI systems will be used to construct **competing mechanistic hypotheses** rather than only a single predictive model.

A hypothesis should be represented as:

```text
biological evidence
       ↓
mechanistic hypothesis
       ↓
mathematical representation
       ↓
predictions
       ↓
uncertainty / counter-evidence
       ↓
independent validation
```

Candidate model classes include ordinary and stochastic differential equations, state-space models, probabilistic dynamical systems, Neural ODEs and related neural dynamical models, optimal-transport models of population-state transitions, and symbolic regression/sparse equation discovery.

LLMs may assist with literature-grounded hypothesis generation, biological interpretation and evidence synthesis, but **LLM output is never treated as biological ground truth**.

---

## Objective 7 — Active selection of informative experiments

A key long-term objective is to move from passive analysis toward **active scientific discovery**.

For competing hypotheses \(H_1,H_2,\ldots,H_k\), an experiment \(e\) can be evaluated by expected information gain:

\[
e^* = \arg\max_e \operatorname{EIG}(e).
\]

Candidate experiments may differ in gene/module perturbation, perturbation strength, perturbation timing, measurement modality, sampling time and cell population.

The aim is to answer questions such as:

> **Which gene/module and which time point would most strongly discriminate between two competing mechanisms?**

This creates a direct bridge between computational inference and future wet-lab validation.

---

# AI and LLM role

AI is a methodological component of the project, not an end in itself.

### Machine learning

Used for multimodal representation learning, latent state inference, trajectory reconstruction, perturbation-response prediction, uncertainty estimation and nonlinear dynamical modeling.

### Biological foundation models

Potentially used to provide pretrained representations of genes, cells and molecular states, particularly for single-cell and multimodal data.

### LLMs

Used as a **knowledge and reasoning layer** for literature-grounded biological interpretation, candidate mechanism generation, mapping between genes/pathways/regulatory processes, evidence and counter-evidence extraction, comparison of competing explanations and generation of testable predictions.

The system should be self-critical: a proposed mechanism must produce predictions that can be evaluated against data.

### Mathematical AI

AI-generated hypotheses should ultimately be translated into quantitative models that can be simulated, fitted, compared and falsified.

---

# Methodological principles

## 1. Alignment is not biological validation

Two trajectories can be made geometrically similar without representing the same biological process.

> **Post-alignment similarity is evidence about geometry, not proof of biological equivalence.**

Time warping and trajectory alignment are therefore treated as diagnostic or comparative tools unless they improve independent prediction.

## 2. Time must not define the state by itself

A method must not obtain its apparent biological signal simply because samples measured at similar times were forced to be similar.

## 3. Validation must be leakage-free

Any operation that learns from data — including feature selection, imputation, normalization, dimensionality reduction, module selection or model fitting — must be performed inside the appropriate training fold.

## 4. Independent experiments are the strongest test

The preferred validation hierarchy is:

```text
replicate holdout
      ↓
timepoint holdout
      ↓
dataset holdout
      ↓
independent experiment
      ↓
new perturbation / wet-lab validation
```

## 5. Prediction must beat meaningful baselines

A model is not considered predictive merely because it correlates with time. It should be compared against persistence, nearest-time observation, linear extrapolation and appropriate conventional/statistical models.

## 6. Statistical significance is not sufficient

A small permutation p-value can coexist with poor predictive performance. Model support therefore requires both statistical evidence and meaningful predictive improvement over appropriate baselines.

## 7. Symbolic equations are hypotheses

An interpretable equation is a candidate mathematical explanation. It does not establish causality without appropriate intervention and validation.

## 8. Uncertainty and provenance are first-class objects

The framework should distinguish:

```text
OBSERVED
   ↓
MEASURED
   ↓
DERIVED
   ↓
INFERRED
   ↓
PREDICTED
```

Every inferred or predicted quantity should retain uncertainty, provenance and the evidence used to derive it.

---

# Current research program

The repository has deliberately progressed through falsifiable stages. Earlier stages tested whether heterogeneous datasets could be aligned into a common trajectory. Later stages progressively tested whether this apparent structure represented a transferable biological state.

## Stage 2.1–2.2 — trajectory geometry

Initial PCA trajectories and Procrustes alignment were used to investigate whether datasets could be geometrically compared.

These analyses established the usefulness and limitations of trajectory geometry but did not establish biological equivalence.

## Stage 2.4 — leakage-free validation

Leave-one-timepoint-out validation was implemented with PCA and scaling fitted only on training timepoints.

## Stage 2.6 — common human gene space

A cross-dataset gene-level harmonization framework was developed using gene identifiers, orthology and platform annotation. The current common space contains approximately 11,899 genes across five expression datasets.

The resulting representation is treated as a research object to validate, not as an assumed biological truth.

## Stage 2.7–2.9 — increasingly rigorous common-state tests

Multiple analyses tested leave-one-dataset-out prediction, leave-one-replicate-out prediction, permutation nulls, biological program representations, latent progress coordinates, temporal concordance, lagged relationships, dataset-specific effects, shared and context-specific dynamics, predictive state representations, Markov sufficiency, invariant biological coordinates, shared temporal components and trajectory-family structure.

The key conclusion is methodological:

> **The data contain reproducible temporal structure, but current analyses do not support the stronger claim that all experiments share one transferable predictive biological trajectory.**

This motivates the shift from searching for a universal trajectory toward **discovering conserved mechanisms and context-dependent dynamics**.

## Stage 2.10 — conserved transition modules

Stage 2.10 tested whether fixed biological programs exhibited a reproducible ordering of activation across independent datasets.

The latest analysis did **not** support a conserved ordering strongly enough to establish a transferable module sequence.

This motivates the next research direction rather than being treated as a failure of the project.

## Stage 2.11 — perturbation and mechanism discovery

The current development branch is dedicated to moving from trajectory similarity toward **perturbation-dependent mechanism discovery**.

The target workflow is:

```text
state representation
       ↓
perturbation response
       ↓
candidate mechanism
       ↓
mathematical model
       ↓
prediction
       ↓
falsification
```

---

# Biological case study: cellular reprogramming

Cellular reprogramming provides a useful test system because it contains large changes in cellular identity, transient intermediate states, heterogeneous cell responses, transcriptional and epigenetic remodeling, metabolic and proliferative changes, perturbation-dependent trajectories and experimentally accessible time courses.

OSKM-mediated reprogramming is an important historical case in the repository, but the project is intentionally expanding toward other forms of cell-state transition and reprogramming.

Potential mechanisms include:

```text
somatic identity
      ↓
cellular plasticity
      ↓
EMT / MET remodeling
      ↓
chromatin / enhancer remodeling
      ↓
pluripotency network
      ↓
stabilized cell state
```

This is a biological hypothesis, not an assumed causal ordering. The project must determine which relationships are supported by data and which are context-dependent.

---

# Current datasets

| Dataset | Organism | Modality | Role |
|---|---|---|---|
| GSE148158 | Human | RNA-seq | Reprogramming / OSKM and controls |
| GSE28688 | Human | Microarray | Early reprogramming time course |
| GSE52052 | Human | Agilent microarray | Reprogramming conditions |
| GSE67462 | Mouse | Affymetrix microarray | Detailed OSKM time course |
| GSE67520 | Mouse | ChIP-seq | Regulatory / chromatin information |
| GSE297234 | Human | scRNA-seq | Single-cell reprogramming trajectories |

The dataset collection should expand toward newer longitudinal and perturbational multi-omics experiments. New datasets should be added because they answer a specific scientific question — not simply because they are large or recent.

Priority should be given to datasets providing new information about perturbations, single-cell heterogeneity, chromatin regulation, protein-level state, alternative cell fates, in vivo transitions and independent experimental validation.

---

# Repository structure

Important current components include:

```text
Data/
Dynamics.py
GSE148158.py
GSE28688.py
GSE297234.py
GSE52052.py
GSE67462.py
GSE67520.py
Summary.py
validate_pipeline.py

dynamics/
    stage24.py
    stage26_v2.py
    stage28.py
    stage2928.py
    stage2929.py
    stage2930.py
    stage2931.py
    stage2932.py
    stage210.py
    validation.py
```

The exact set of experimental stages may evolve as hypotheses are falsified and replaced.

---

# Reproducibility and validation

Analysis outputs are stored under `results/`, with Dynamics outputs under `results/Dynamics/`.

Validation reports should identify datasets used for training, held-out datasets or timepoints, selected features, preprocessing operations, model parameters, baselines, uncertainty estimates, permutation/bootstrap results and the final scientific interpretation.

The project favors explicit negative results over hidden methodological failures.

---

# Expected scientific contribution

The intended contribution is a **validated methodology for AI-assisted discovery of biological dynamics**, rather than simply a new predictor or a software package.

A successful outcome would demonstrate a chain such as:

```text
heterogeneous observations
          ↓
biological representation
          ↓
latent cellular state
          ↓
conserved + context-specific dynamics
          ↓
perturbation response
          ↓
competing mechanisms
          ↓
mathematical models
          ↓
predictions + uncertainty
          ↓
falsification
          ↓
informative experiment
```

The final Digital Biological Twin would provide an interpretable interface to this chain.

A scientifically successful result does **not** require a single universal trajectory. It may instead show that different biological contexts require different dynamical mechanisms while sharing a smaller set of conserved principles.

---

# Long-term vision

The long-term objective is a system that can take heterogeneous observations of a biological system and reason about its future state in a scientifically testable way:

```text
                 BIOLOGICAL DATA
                       │
                       ▼
              KNOWLEDGE + AI
                       │
                       ▼
                CELLULAR STATE
                       │
              ┌────────┴────────┐
              ▼                 ▼
       SHARED DYNAMICS     CONTEXT / MEMORY
              │                 │
              └────────┬────────┘
                       ▼
              MECHANISTIC MODELS
                       │
                       ▼
             PREDICTION / UNCERTAINTY
                       │
                       ▼
              MODEL FALSIFICATION
                       │
                       ▼
             OPTIMAL PERTURBATION
                       │
                       ▼
             NEW BIOLOGICAL DATA
                       │
                       └──────────↺
```

This closed loop is the conceptual foundation for a future **Digital Biological Twin** and for AI systems that assist, rather than merely imitate, the scientific process.

---

## Status

This repository is an active research project. Conclusions are intentionally updated as new validation experiments are performed.

The current evidence supports reproducible temporal and biological structure in the analyzed datasets, but **does not yet establish a universally transferable predictive state or mechanism**.

That distinction is central to the project: the goal is not to force the data to support a predetermined story, but to build methods capable of discovering what is conserved, what is context-dependent, and what can genuinely be predicted and experimentally tested.