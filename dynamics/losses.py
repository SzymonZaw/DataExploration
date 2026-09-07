"""Loss functions for dynamic state learning."""

from __future__ import annotations

from typing import Optional

import torch
import torch.nn.functional as F


def masked_mse(pred: torch.Tensor, target: torch.Tensor, mask: Optional[torch.Tensor] = None) -> torch.Tensor:
    if pred.shape != target.shape:
        raise ValueError(f"shape mismatch: {pred.shape} != {target.shape}")
    loss = (pred - target) ** 2
    if mask is None:
        return loss.mean()
    mask = mask.to(dtype=loss.dtype)
    if mask.shape != loss.shape:
        mask = mask.expand_as(loss)
    return (loss * mask).sum() / mask.sum().clamp_min(1.0)


def state_prediction_loss(predicted_state: torch.Tensor, target_state: torch.Tensor) -> torch.Tensor:
    return F.mse_loss(predicted_state, target_state)


def reconstruction_loss(reconstruction: torch.Tensor, observation: torch.Tensor,
                        mask: Optional[torch.Tensor] = None) -> torch.Tensor:
    return masked_mse(reconstruction, observation, mask)


def future_observation_loss(predicted_observation: torch.Tensor, next_observation: torch.Tensor,
                            mask: Optional[torch.Tensor] = None) -> torch.Tensor:
    """Require predicted latent state to retain information about the future."""
    return masked_mse(predicted_observation, next_observation, mask)


def total_loss(reconstruction: Optional[torch.Tensor], observation: torch.Tensor,
               predicted_state: torch.Tensor, target_state: torch.Tensor,
               predicted_observation: Optional[torch.Tensor] = None,
               next_observation: Optional[torch.Tensor] = None,
               reconstruction_weight: float = 1.0,
               prediction_weight: float = 1.0,
               future_weight: float = 1.0,
               mask: Optional[torch.Tensor] = None) -> torch.Tensor:
    if min(reconstruction_weight, prediction_weight, future_weight) < 0:
        raise ValueError("loss weights must be non-negative")
    loss = prediction_weight * state_prediction_loss(predicted_state, target_state)
    if reconstruction is not None:
        loss = loss + reconstruction_weight * reconstruction_loss(reconstruction, observation, mask)
    if predicted_observation is not None and next_observation is not None:
        loss = loss + future_weight * future_observation_loss(predicted_observation, next_observation, mask)
    return loss
