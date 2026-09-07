"""Training utilities for the central dynamic state model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import torch
from torch import nn

from .dynamic_state_model import DynamicStateModel
from .losses import total_loss


@dataclass
class TrainConfig:
    epochs: int = 200
    learning_rate: float = 1e-3
    reconstruction_weight: float = 1.0
    prediction_weight: float = 0.25
    future_weight: float = 1.0
    grad_clip: Optional[float] = 1.0


def train_epoch(model: DynamicStateModel, batches: Iterable[dict[str, torch.Tensor]],
                optimizer: torch.optim.Optimizer, config: TrainConfig) -> float:
    """Run one epoch over already split batches.

    Normalization, imputation, feature selection and train/test splitting must be
    performed before this function using training data only.
    """
    model.train()
    losses = []
    for batch in batches:
        optimizer.zero_grad(set_to_none=True)
        z_t = model.encode(batch["x_t"])
        z_next = model.encode(batch["x_next"]).detach()
        reconstruction = model.decode(z_t)
        predicted = model.transition(z_t, context=batch.get("context"), history=batch.get("history"))
        predicted_observation = model.decode(predicted)
        loss = total_loss(
            reconstruction, batch["x_t"], predicted, z_next,
            predicted_observation=predicted_observation,
            next_observation=batch["x_next"],
            reconstruction_weight=config.reconstruction_weight,
            prediction_weight=config.prediction_weight,
            future_weight=config.future_weight,
        )
        loss.backward()
        if config.grad_clip is not None:
            nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        optimizer.step()
        losses.append(float(loss.detach().cpu()))
    return sum(losses) / max(len(losses), 1)


def build_model(input_dim: int, state_dim: int = 16, hidden_dim: int = 128,
                context_dim: int = 0, history_dim: int = 0) -> DynamicStateModel:
    return DynamicStateModel(
        input_dim=input_dim, state_dim=state_dim, hidden_dim=hidden_dim,
        context_dim=context_dim, history_dim=history_dim,
    )
