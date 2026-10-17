"""Periodic finite-difference operators on uniform 3D grids."""

from __future__ import annotations

import numpy as np


def periodic_derivative(field: np.ndarray, axis: int, spacing: float) -> np.ndarray:
    """Return the centered periodic derivative along one axis."""

    return (np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis)) / (
        2.0 * spacing
    )


def gradient(
    field: np.ndarray, dx: float, dy: float, dz: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return periodic centered gradients along x, y, z."""

    return (
        periodic_derivative(field, axis=0, spacing=dx),
        periodic_derivative(field, axis=1, spacing=dy),
        periodic_derivative(field, axis=2, spacing=dz),
    )


def divergence(
    vx: np.ndarray,
    vy: np.ndarray,
    vz: np.ndarray,
    dx: float,
    dy: float,
    dz: float,
) -> np.ndarray:
    """Return periodic centered divergence of a vector field."""

    return (
        periodic_derivative(vx, axis=0, spacing=dx)
        + periodic_derivative(vy, axis=1, spacing=dy)
        + periodic_derivative(vz, axis=2, spacing=dz)
    )


def magnitude(gx: np.ndarray, gy: np.ndarray, gz: np.ndarray) -> np.ndarray:
    """Return Euclidean magnitude without squaring extreme finite components."""

    return np.hypot(np.hypot(gx, gy), gz)
