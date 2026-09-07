"""Deep-learning model for learning cellular state and its dynamics.

This module intentionally contains a model architecture, not a biological claim.
The model learns a low-dimensional latent state from heterogeneous observations and
predicts future latent states. Interpretation and causal claims must be validated
out-of-sample.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import torch
from torch import nn


@dataclass
class DynamicStateOutput:
    """Outputs returned by :class:`DynamicStateModel`."""

    state: torch.Tensor
    predicted_state: Optional[torch.Tensor] = None
    reconstruction: Optional[torch.Tensor] = None


class MLP(nn.Module):
    """Small configurable multilayer perceptron."""

    def __init__(self, input_dim: int, hidden_dim: int, output_dim: int, dropout: float = 0.0):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, output_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class DynamicStateModel(nn.Module):
    """Learn a latent biological state and its temporal transition.

    Parameters
    ----------
    input_dim:
        Number of input features in the harmonized observation space.
    state_dim:
        Dimension of the learned latent state.
    hidden_dim:
        Hidden width of encoder/decoder/dynamics networks.
    context_dim:
        Optional dimension of a perturbation/context vector ``u``.
    history_dim:
        Optional dimension of a history representation. If non-zero, history is
        represented by a GRU over previous latent states.
    dropout:
        Dropout used in feed-forward networks.

    Notes
    -----
    The first implementation is deliberately modest. It provides a stable API for
    experiments before introducing more expressive architectures such as Neural
    ODEs, transformers or stochastic state-space models.
    """

    def __init__(
        self,
        input_dim: int,
        state_dim: int = 16,
        hidden_dim: int = 128,
        context_dim: int = 0,
        history_dim: int = 0,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if input_dim <= 0 or state_dim <= 0 or hidden_dim <= 0:
            raise ValueError("input_dim, state_dim and hidden_dim must be positive")
        if context_dim < 0 or history_dim < 0:
            raise ValueError("context_dim and history_dim must be non-negative")

        self.input_dim = input_dim
        self.state_dim = state_dim
        self.context_dim = context_dim
        self.history_dim = history_dim

        self.encoder = MLP(input_dim, hidden_dim, state_dim, dropout)
        self.decoder = MLP(state_dim, hidden_dim, input_dim, dropout)

        self.history_encoder = None
        dynamics_input = state_dim + context_dim
        if history_dim:
            self.history_encoder = nn.GRU(state_dim, history_dim, batch_first=True)
            dynamics_input += history_dim

        self.dynamics = MLP(dynamics_input, hidden_dim, state_dim, dropout)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Map observations ``x`` to latent state ``z``."""
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Reconstruct observations from latent state."""
        return self.decoder(z)

    def encode_history(self, history: torch.Tensor) -> torch.Tensor:
        """Encode a sequence of previous latent states into a history vector."""
        if self.history_encoder is None:
            raise ValueError("history_dim was set to zero")
        _, hidden = self.history_encoder(history)
        return hidden[-1]

    def transition(
        self,
        z: torch.Tensor,
        context: Optional[torch.Tensor] = None,
        history: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Predict the next latent state from current state and optional context/history."""
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

    def predict(
        self,
        x: torch.Tensor,
        context: Optional[torch.Tensor] = None,
        history: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        """Encode ``x`` and predict the next latent state."""
        z = self.encode(x)
        return self.transition(z, context=context, history=history)

    def forward(
        self,
        x: torch.Tensor,
        context: Optional[torch.Tensor] = None,
        history: Optional[torch.Tensor] = None,
        reconstruct: bool = True,
    ) -> DynamicStateOutput:
        z = self.encode(x)
        reconstruction = self.decode(z) if reconstruct else None
        predicted = self.transition(z, context=context, history=history)
        return DynamicStateOutput(
            state=z,
            predicted_state=predicted,
            reconstruction=reconstruction,
        )
