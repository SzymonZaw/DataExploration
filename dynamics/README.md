# Dynamic state model

This directory now contains the first implementation of the central AI model proposed for the PhD project.

## Scientific role

The goal is to learn, rather than manually specify, a low-dimensional representation of cellular state and its temporal dynamics from harmonized biological observations.

The model is designed around:

```text
observation x(t)
      ↓
encoder
      ↓
latent state z(t)
      ↓
transition model
      ↓
predicted z(t+Δt)
      ↓
decoder
      ↓
reconstructed observation
```

Optional context/perturbation `u(t)` and history `h(t)` can be provided to the transition model. This allows future experiments to test whether cell-state dynamics are Markovian or require information about previous states.

## Modules

- `dynamic_state_model.py` — encoder, decoder, temporal transition and optional history/context representation.
- `losses.py` — reconstruction and latent prediction objectives, including missing-value masks.
- `train_dynamic_state.py` — leakage-aware training scaffold. Data splitting, normalization and imputation deliberately remain outside the trainer.
- `model_benchmark.py` — unified leakage-free forecasting benchmark used by Z6.
- `run_model_benchmark.py` — frozen Phase 0 benchmark runner.
- `Z6_PREDICTIVE_TRANSITION_PROTOCOL.md` — doctoral Z6 protocol for out-of-sample prediction of future cellular state.
- `../run_z6_predictive_transition.py` — dedicated Z6 executable wrapper.

## Z6 predictive transition benchmark

Z6 is now operationally defined as a genuine future-state forecasting experiment:

```text
held-out trajectory prefix
          ↓
state inference
          ↓
transition model
          ↓
free multi-step rollout
          ↓
predicted future state
          ↓
held-out future observation
```

The benchmark uses leave-one-dataset-out validation, training-only preprocessing, persistence/nearest-time/linear baselines and a temporal permutation null. Candidate models are evaluated under the same forecasting protocol so that the representation/model, rather than the evaluation procedure, is the experimental variable.

Run the official benchmark with:

```powershell
python run_z6_predictive_transition.py
```

For a diagnostic smoke test only:

```powershell
python run_z6_predictive_transition.py --epochs 25 --permutation-n 50 --seeds 411
```

The official Z6 result is not established until the full benchmark has been executed and its predefined predictive-support criteria have been evaluated.

## Research constraints

The implementation must be evaluated with strict dataset-level and time-based holdouts. In particular, preprocessing statistics, feature selection and imputation parameters must be fitted on training data only.

A good reconstruction loss is not sufficient evidence of a biological state. The main evaluation target is prediction on observations that were not used to learn the representation.

Predictive performance alone is also not evidence of biological specificity. Z6 results must subsequently be audited for context-dependent and technical nuisance signals under the Z4 framework.

## Planned model progression

The architecture is deliberately simple at first. It can later be extended with:

1. probabilistic latent states and uncertainty,
2. continuous-time / Neural ODE dynamics,
3. modality-specific encoders for RNA, scRNA and regulatory data,
4. perturbation embeddings,
5. explicit history/memory models,
6. sparse/Jacobian-based interpretability,
7. symbolic regression over learned latent dynamics.

These are model variants, not separate scientific stages. The scientific question remains constant: **can AI learn a biologically meaningful state representation whose dynamics generalize across independent experiments?**
