"""Loss functions for DynamicStateModel experiments."""

from __future__ import annotations

from typing import Optional

import torch
import torch.nn.functional as F


def masked_mse(pred: torch.Tensor, target: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
    """Mean squared error while ignoring missing observations when a mask is given."""
    if pred.shape != target.shape:
        raise ValueError(f"shape mismatch: {pred.shape} != {target.shape}")
    loss = (pred - target) ** 2
    if mask is None:
        return loss.mean()
    mask = mask.to(dtype=loss.dtype)
    if mask.shape != loss.shape:
        mask = mask.expand_as(loss)
    denom = mask.sum().clamp_min(1.0)
    return (loss * mask).sum() / denom


def state_prediction_loss(predicted_state: torch.Tensor, target_state: torch.Tensor) -> torch.Tensor:
    """Loss for one-step latent-state prediction."""
    return F.mse_loss(predicted_state, target_state)


def reconstruction_loss(reconstruction: torch.Tensor, observation: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
    """Observation reconstruction loss."""
    return masked_mse(reconstruction, observation, mask)


def total_loss(
    reconstruction: Optional[torch.Tensor],
    observation: torch.Tensor,
    predicted_state: torch.Tensor,
    target_state: torch.Tensor,
    reconstruction_weight: float = 1.0,
    prediction_weight: float = 1.0,
    mask: Optional[torch.Tensor] = None,
) -> torch.Tensor:
    """Combine reconstruction and temporal prediction objectives."""
    if reconstruction_weight < 0 or prediction_weight < 0:
        raise ValueError("loss weights must be non-negative")
    loss = prediction_weight * state_prediction_loss(predicted_state, target_state)
    if reconstruction is not None:
        loss = loss + reconstruction_weight * reconstruction_loss(reconstruction, observation, mask)
    return loss
