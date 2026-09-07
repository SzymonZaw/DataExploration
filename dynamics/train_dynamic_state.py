"""Training utilities for the central dynamic state model.

This is intentionally a scaffold: data loading and experiment-specific batching stay
outside this module so that leakage-free train/validation/test splits remain explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional

import torch
from torch import nn

from .dynamic_state_model import DynamicStateModel
from .losses import total_loss


@dataclass
class TrainConfig:
    epochs: int = 100
    learning_rate: float = 1e-3
    reconstruction_weight: float = 1.0
    prediction_weight: float = 1.0
    grad_clip: Optional[float] = 1.0


def train_epoch(
    model: DynamicStateModel,
    batches: Iterable[dict[str, torch.Tensor]],
    optimizer: torch.optim.Optimizer,
    config: TrainConfig,
) -> float:
    """Run one training epoch over pre-split batches.

    Expected batch keys are ``x_t`` and ``x_next``. Optional keys are ``context``,
    ``history`` and ``mask``. Splitting, normalization and imputation must happen
    before this function using training data only.
    """
    model.train()
    losses = []
    for batch in batches:
        optimizer.zero_grad(set_to_none=True)
        z_t = model.encode(batch["x_t"])
        z_next = model.encode(batch["x_next"]).detach()
        reconstruction = model.decode(z_t)
        predicted = model.transition(
            z_t,
            context=batch.get("context"),
            history=batch.get("history"),
        )
        loss = total_loss(
            reconstruction,
            batch["x_t"],
            predicted,
            z_next,
            reconstruction_weight=config.reconstruction_weight,
            prediction_weight=config.prediction_weight,
            mask=batch.get("mask"),
        )
        loss.backward()
        if config.grad_clip is not None:
            nn.utils.clip_grad_norm_(model.parameters(), config.grad_clip)
        optimizer.step()
        losses.append(float(loss.detach().cpu()))
    return sum(losses) / max(len(losses), 1)


def build_model(input_dim: int, state_dim: int = 16, hidden_dim: int = 128, context_dim: int = 0, history_dim: int = 0) -> DynamicStateModel:
    """Construct the default research model."""
    return DynamicStateModel(
        input_dim=input_dim,
        state_dim=state_dim,
        hidden_dim=hidden_dim,
        context_dim=context_dim,
        history_dim=history_dim,
    )
