"""Deep-learning model for learning cellular state and its dynamics.

This module contains a model architecture, not a biological claim. The model learns
a low-dimensional latent state from harmonized observations and predicts future
observations. Interpretation and causal claims require independent validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
from torch import nn


@dataclass
class DynamicStateOutput:
    state: torch.Tensor
    predicted_state: Optional[torch.Tensor] = None
    reconstruction: Optional[torch.Tensor] = None
    predicted_observation: Optional[torch.Tensor] = None


class MLP(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int, dropout: float = 0.0):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden_dim, output_dim)
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class DynamicStateModel(nn.Module):
    """Learn latent cellular state and its temporal transition.

    The model can operate in three modes: current observation -> latent state,
    latent state -> future latent state, and future latent state -> observation.
    Optional context/perturbation and history inputs allow later experiments to
    test context dependence and memory without changing the core API.
    """

    def __init__(self, input_dim: int, state_dim: int = 16, hidden_dim: int = 128,
                 context_dim: int = 0, history_dim: int = 0, dropout: float = 0.1) -> None:
        super().__init__()
        if input_dim <= 0 or state_dim <= 0 or hidden_dim <= 0:
            raise ValueError("input_dim, state_dim and hidden_dim must be positive")
        if context_dim < 0 or history_dim < 0:
            raise ValueError("context_dim and history_dim must be non-negative")
        self.input_dim, self.state_dim = input_dim, state_dim
        self.context_dim, self.history_dim = context_dim, history_dim
        self.encoder = MLP(input_dim, hidden_dim, state_dim, dropout)
        self.decoder = MLP(state_dim, hidden_dim, input_dim, dropout)
        self.history_encoder = nn.GRU(state_dim, history_dim, batch_first=True) if history_dim else None
        self.dynamics = MLP(state_dim + context_dim + history_dim, hidden_dim, state_dim, dropout)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def encode_history(self, history: torch.Tensor) -> torch.Tensor:
        if self.history_encoder is None:
            raise ValueError("history_dim was set to zero")
        _, hidden = self.history_encoder(history)
        return hidden[-1]

    def transition(self, z: torch.Tensor, context: Optional[torch.Tensor] = None,
                   history: Optional[torch.Tensor] = None) -> torch.Tensor:
        parts = [z]
        if self.context_dim:
            if context is None:
                raise ValueError("context is required when context_dim > 0")
            parts.append(context)
        elif context is not None:
            raise ValueError("context was supplied but context_dim is zero")
        if self.history_dim:
            if history is None:
                raise ValueError("history is required when history_dim > 0")
            parts.append(self.encode_history(history))
        elif history is not None:
            raise ValueError("history was supplied but history_dim is zero")
        return self.dynamics(torch.cat(parts, dim=-1))

    def predict(self, x: torch.Tensor, context: Optional[torch.Tensor] = None,
                history: Optional[torch.Tensor] = None) -> torch.Tensor:
        return self.transition(self.encode(x), context=context, history=history)

    def predict_observation(self, x: torch.Tensor, context: Optional[torch.Tensor] = None,
                            history: Optional[torch.Tensor] = None) -> torch.Tensor:
        """Predict the next observation directly through the learned latent state."""
        return self.decode(self.predict(x, context=context, history=history))

    def forward(self, x: torch.Tensor, context: Optional[torch.Tensor] = None,
                history: Optional[torch.Tensor] = None, reconstruct: bool = True) -> DynamicStateOutput:
        z = self.encode(x)
        predicted = self.transition(z, context=context, history=history)
        return DynamicStateOutput(
            state=z,
            predicted_state=predicted,
            reconstruction=self.decode(z) if reconstruct else None,
            predicted_observation=self.decode(predicted),
        )
