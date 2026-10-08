"""Helpers for candidate shock criteria."""

from __future__ import annotations

import numpy as np


def aligned_dot(
    ax: np.ndarray,
    ay: np.ndarray,
    az: np.ndarray,
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
) -> np.ndarray:
    """Return vector dot product."""

    return ax * bx + ay * by + az * bz


def alignment_sign(
    ax: np.ndarray, ay: np.ndarray, az: np.ndarray,
    bx: np.ndarray, by: np.ndarray, bz: np.ndarray,
) -> np.ndarray:
    """Return a bounded dot product with the same sign for finite vectors.

    Candidate predicates need direction agreement, not dot magnitude. Divide
    each vector by its largest absolute component before multiplying: positive
    scaling preserves the sign while avoiding overflow/cancellation of infinities.
    A zero vector stays zero; aligned_dot retains its public raw-dot semantics.
    """
    scale_a = np.maximum.reduce([np.abs(ax), np.abs(ay), np.abs(az)])
    scale_b = np.maximum.reduce([np.abs(bx), np.abs(by), np.abs(bz)])
    scale_a = np.where(scale_a > 0, scale_a, 1.)
    scale_b = np.where(scale_b > 0, scale_b, 1.)
    return aligned_dot(ax/scale_a, ay/scale_a, az/scale_a,
                       bx/scale_b, by/scale_b, bz/scale_b)
