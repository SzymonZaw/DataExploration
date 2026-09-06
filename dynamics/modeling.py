"""Dynamical modeling components (Stage 3 and later).

This module is intentionally only a boundary for now. No dynamical model is
introduced by the structural refactor.
"""


def infer_state_dynamics(*args, **kwargs):
    """Future entry point for estimating dz/dt = F(z, t, u, theta)."""
    raise NotImplementedError("Stage 3 dynamics has not been implemented yet")


def fit_symbolic_model(*args, **kwargs):
    """Future entry point for interpretable model discovery."""
    raise NotImplementedError("Symbolic model discovery has not been implemented yet")


__all__ = ["infer_state_dynamics", "fit_symbolic_model"]
