"""Trajectory analysis components (Stage 2.1–2.4).

During the refactor this module provides the stable home for trajectory logic.
The existing implementations remain in Dynamics.py until they can be moved
without changing behavior.
"""

from Dynamics import apply_transform, curve, fit_transform, rotation
from Dynamics import stage2_1, stage2_2, stage2_4

__all__ = [
    "apply_transform",
    "curve",
    "fit_transform",
    "rotation",
    "stage2_1",
    "stage2_2",
    "stage2_4",
]
