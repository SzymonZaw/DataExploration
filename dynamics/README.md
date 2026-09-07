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

## Research constraints

The implementation must be evaluated with strict dataset-level and time-based holdouts. In particular, preprocessing statistics, feature selection and imputation parameters must be fitted on training data only.

A good reconstruction loss is not sufficient evidence of a biological state. The main evaluation target is prediction on observations that were not used to learn the representation.

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
