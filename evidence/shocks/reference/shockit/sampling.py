"""Sampling helpers for upstream/downstream states."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import map_coordinates


@dataclass
class SampledJumps:
    """Paired upstream/downstream states and ratios across the configured width."""

    pressure_jump: np.ndarray
    temperature_jump: np.ndarray
    density_jump: np.ndarray
    pressure_up: np.ndarray
    pressure_down: np.ndarray
    temperature_up: np.ndarray
    temperature_down: np.ndarray
    density_up: np.ndarray
    density_down: np.ndarray


def dominant_axis(nx: np.ndarray, ny: np.ndarray, nz: np.ndarray) -> np.ndarray:
    """Return the largest absolute normal component; ties prefer x, then y, then z."""

    components = np.stack([np.abs(nx), np.abs(ny), np.abs(nz)], axis=0)
    return np.argmax(components, axis=0)


def axis_signs(
    nx: np.ndarray, ny: np.ndarray, nz: np.ndarray, axes: np.ndarray
) -> np.ndarray:
    """Return the chosen normal component's sign, treating zero as positive."""

    signs = np.sign(nx)
    signs = np.where(axes == 1, np.sign(ny), signs)
    signs = np.where(axes == 2, np.sign(nz), signs)
    return np.where(signs == 0.0, 1.0, signs)


def sample_axis_offsets(
    field: np.ndarray, axes: np.ndarray, signs: np.ndarray, width: int
) -> tuple[np.ndarray, np.ndarray]:
    """Sample plus/minus states along the nearest-axis normal direction."""

    plus = np.empty_like(field, dtype=float)
    minus = np.empty_like(field, dtype=float)
    for axis in range(3):
        axis_mask = axes == axis
        if not np.any(axis_mask):
            continue
        positive = axis_mask & (signs >= 0.0)
        negative = axis_mask & (signs < 0.0)
        if np.any(positive):
            plus_field = np.roll(field, -width, axis=axis)
            minus_field = np.roll(field, width, axis=axis)
            plus[positive] = plus_field[positive]
            minus[positive] = minus_field[positive]
        if np.any(negative):
            plus_field = np.roll(field, width, axis=axis)
            minus_field = np.roll(field, -width, axis=axis)
            plus[negative] = plus_field[negative]
            minus[negative] = minus_field[negative]
    return plus, minus


def sample_jumps(
    pressure: np.ndarray,
    temperature: np.ndarray,
    rho: np.ndarray,
    nx: np.ndarray,
    ny: np.ndarray,
    nz: np.ndarray,
    width: int,
    method: str = "nearest_axis",
    candidate_mask: np.ndarray | None = None,
    *,
    physical_distance: float | None = None,
    spacing: tuple[float, float, float] | None = None,
) -> SampledJumps:
    """Sample paired states using a width measured in local cell-index units.

    nx, ny, nz are normal components. Trilinear offsets use width times each
    component even for unequal physical spacing; this is the existing index-space
    convention. Only trilinear sampling uses candidate_mask to limit work.
    With physical_distance, width is ignored: endpoints have that physical reach
    along the normal (trilinear) or selected axis (nearest_axis), using spacing.
    Fixed endpoints do not locate plateaus or exclude intervening features.
    """

    if physical_distance is not None:
        if (not np.isscalar(physical_distance) or isinstance(physical_distance, (bool, np.bool_))
            or not np.isfinite(physical_distance) or physical_distance <= 0):
            raise ValueError("physical_distance must be finite and positive")
        steps = np.asarray(spacing, dtype=float)
        if steps.shape != (3,) or not np.all(np.isfinite(steps)) or np.any(steps <= 0):
            raise ValueError("physical sampling requires three finite positive spacing values")
        if not np.all(np.isfinite(physical_distance / steps)):
            raise ValueError("physical sampling distance exceeds representable index offsets")
        directions = (nx, ny, nz)
        if method == "nearest_axis":
            axes = dominant_axis(nx, ny, nz)
            signs = axis_signs(nx, ny, nz, axes)
            directions = tuple(np.where(axes == axis, signs, 0.) for axis in range(3))
        elif method != "trilinear":
            raise ValueError(f"Unsupported sampling method: {method}")
        # Normals are physical unit vectors. Convert each endpoint displacement
        # independently to cell-index units; do NOT renormalize after division.
        # Even nearest-axis uses linear interpolation for a fractional distance.
        # Form distance/spacing first: 1/spacing can overflow even when the
        # final displacement is small (e.g. both in subnormal physical units).
        offsets = tuple(direction * cells for direction, cells in zip(directions, physical_distance / steps))
        return _sample_jumps_trilinear(
            pressure, temperature, rho, *offsets, 1.,
            candidate_mask=candidate_mask if method == "trilinear" else None,
        )
    if method == "nearest_axis":
        return _sample_jumps_nearest_axis(pressure, temperature, rho, nx, ny, nz, width)
    if method == "trilinear":
        return _sample_jumps_trilinear(
            pressure, temperature, rho, nx, ny, nz, width, candidate_mask=candidate_mask
        )
    raise ValueError(f"Unsupported sampling method: {method}")


def _sample_jumps_nearest_axis(
    pressure: np.ndarray,
    temperature: np.ndarray,
    rho: np.ndarray,
    nx: np.ndarray,
    ny: np.ndarray,
    nz: np.ndarray,
    width: int,
) -> SampledJumps:
    """Sample jump ratios by stepping along the dominant grid axis."""

    axes = dominant_axis(nx, ny, nz)
    signs = axis_signs(nx, ny, nz, axes)

    pressure_plus, pressure_minus = sample_axis_offsets(pressure, axes, signs, width)
    temperature_plus, temperature_minus = sample_axis_offsets(
        temperature, axes, signs, width
    )
    rho_plus, rho_minus = sample_axis_offsets(rho, axes, signs, width)

    # Pressure chooses the downstream side; keep its density and temperature
    # together as one state, including the defined plus-side tie preference.
    downstream_is_plus = pressure_plus >= pressure_minus
    pressure_down = np.where(downstream_is_plus, pressure_plus, pressure_minus)
    pressure_up = np.where(downstream_is_plus, pressure_minus, pressure_plus)
    temperature_down = np.where(downstream_is_plus, temperature_plus, temperature_minus)
    temperature_up = np.where(downstream_is_plus, temperature_minus, temperature_plus)
    rho_down = np.where(downstream_is_plus, rho_plus, rho_minus)
    rho_up = np.where(downstream_is_plus, rho_minus, rho_plus)

    return SampledJumps(
        pressure_jump=_safe_ratio(pressure_down, pressure_up),
        temperature_jump=_safe_ratio(temperature_down, temperature_up),
        density_jump=_safe_ratio(rho_down, rho_up),
        pressure_up=pressure_up,
        pressure_down=pressure_down,
        temperature_up=temperature_up,
        temperature_down=temperature_down,
        density_up=rho_up,
        density_down=rho_down,
    )


def _safe_ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """Divide only by positive upstream values; leave invalid ratios as NaN."""
    valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0.0)
    ratio = np.full(numerator.shape, np.nan, dtype=float)
    ratio[valid] = numerator[valid] / denominator[valid]
    return ratio


def _sample_trilinear_points(
    field: np.ndarray, coordinates: tuple[np.ndarray, np.ndarray, np.ndarray]
) -> np.ndarray:
    """Interpolate eight periodic corners at floating xyz index coordinates."""
    # Period is N cells on each axis, including the last-to-first interpolation
    # interval. SciPy wrap instead overlaps endpoints with period N-1.
    stacked = np.vstack(coordinates)
    sampled = map_coordinates(field, stacked, order=1, mode="grid-wrap", output=float)
    return np.asarray(sampled, dtype=float)


def _sample_jumps_trilinear(
    pressure: np.ndarray,
    temperature: np.ndarray,
    rho: np.ndarray,
    nx: np.ndarray,
    ny: np.ndarray,
    nz: np.ndarray,
    width: int,
    candidate_mask: np.ndarray | None,
) -> SampledJumps:
    """Sample jump ratios by trilinear interpolation along the local normal."""

    active_mask = (
        np.ones(pressure.shape, dtype=bool)
        if candidate_mask is None
        else candidate_mask.astype(bool, copy=False)
    )
    active_indices = np.argwhere(active_mask)
    empty = np.full(pressure.shape, np.nan, dtype=float)
    if active_indices.size == 0:
        return SampledJumps(
            pressure_jump=empty.copy(),
            temperature_jump=empty.copy(),
            density_jump=empty.copy(),
            pressure_up=empty.copy(),
            pressure_down=empty.copy(),
            temperature_up=empty.copy(),
            temperature_down=empty.copy(),
            density_up=empty.copy(),
            density_down=empty.copy(),
        )

    base_x = active_indices[:, 0].astype(float)
    base_y = active_indices[:, 1].astype(float)
    base_z = active_indices[:, 2].astype(float)
    offset_x = width * nx[active_mask]
    offset_y = width * ny[active_mask]
    offset_z = width * nz[active_mask]

    plus_coordinates = (base_x + offset_x, base_y + offset_y, base_z + offset_z)
    minus_coordinates = (base_x - offset_x, base_y - offset_y, base_z - offset_z)

    pressure_plus_values = _sample_trilinear_points(pressure, plus_coordinates)
    pressure_minus_values = _sample_trilinear_points(pressure, minus_coordinates)
    temperature_plus_values = _sample_trilinear_points(temperature, plus_coordinates)
    temperature_minus_values = _sample_trilinear_points(temperature, minus_coordinates)
    rho_plus_values = _sample_trilinear_points(rho, plus_coordinates)
    rho_minus_values = _sample_trilinear_points(rho, minus_coordinates)

    pressure_plus = empty.copy()
    pressure_minus = empty.copy()
    temperature_plus = empty.copy()
    temperature_minus = empty.copy()
    rho_plus = empty.copy()
    rho_minus = empty.copy()

    pressure_plus[active_mask] = pressure_plus_values
    pressure_minus[active_mask] = pressure_minus_values
    temperature_plus[active_mask] = temperature_plus_values
    temperature_minus[active_mask] = temperature_minus_values
    rho_plus[active_mask] = rho_plus_values
    rho_minus[active_mask] = rho_minus_values

    # Pressure chooses the downstream side; keep its density and temperature
    # together as one state, including the defined plus-side tie preference.
    downstream_is_plus = pressure_plus >= pressure_minus
    pressure_down = np.where(downstream_is_plus, pressure_plus, pressure_minus)
    pressure_up = np.where(downstream_is_plus, pressure_minus, pressure_plus)
    temperature_down = np.where(downstream_is_plus, temperature_plus, temperature_minus)
    temperature_up = np.where(downstream_is_plus, temperature_minus, temperature_plus)
    rho_down = np.where(downstream_is_plus, rho_plus, rho_minus)
    rho_up = np.where(downstream_is_plus, rho_minus, rho_plus)

    return SampledJumps(
        pressure_jump=_safe_ratio(pressure_down, pressure_up),
        temperature_jump=_safe_ratio(temperature_down, temperature_up),
        density_jump=_safe_ratio(rho_down, rho_up),
        pressure_up=pressure_up,
        pressure_down=pressure_down,
        temperature_up=temperature_up,
        temperature_down=temperature_down,
        density_up=rho_up,
        density_down=rho_down,
    )
