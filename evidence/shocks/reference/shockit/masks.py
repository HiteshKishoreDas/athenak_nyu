"""Mask construction helpers for the shock finder."""

from __future__ import annotations

import numpy as np
from scipy import ndimage

from .config import ShockFinderConfig
from .criteria import alignment_sign
from .derived import DerivedFields, GradientFields, MachFields
from .sampling import SampledJumps


def build_shock_zone_mask(
    derived: DerivedFields,
    gradients: GradientFields,
    config: ShockFinderConfig,
) -> np.ndarray:
    """Build the local shock-zone candidate mask."""

    # Candidate compression and gradient thresholds are strict inequalities.
    mask = derived.div_v < config.div_v_threshold
    mask &= gradients.grad_t_mag > config.grad_T_min
    if config.require_gradT_gradS_alignment:
        mask &= (
            alignment_sign(
                gradients.grad_tx,
                gradients.grad_ty,
                gradients.grad_tz,
                gradients.grad_sx,
                gradients.grad_sy,
                gradients.grad_sz,
            )
            > 0.0
        )
    if config.require_gradT_gradRho_alignment:
        mask &= (
            alignment_sign(
                gradients.grad_tx,
                gradients.grad_ty,
                gradients.grad_tz,
                gradients.grad_rx,
                gradients.grad_ry,
                gradients.grad_rz,
            )
            > 0.0
        )
    return mask


def build_final_shock_mask(
    shock_zone_mask: np.ndarray,
    jumps: SampledJumps,
    mach_fields: MachFields,
    config: ShockFinderConfig,
    gamma: float,
) -> np.ndarray:
    """Build the unreduced shock mask after jump and Mach filtering."""

    jump_mask = np.ones(shock_zone_mask.shape, dtype=bool)
    if config.require_pressure_jump:
        jump_mask &= jumps.pressure_jump > 1.0
    if config.require_temperature_jump:
        jump_mask &= jumps.temperature_jump > 1.0
    if config.require_density_jump:
        jump_mask &= jumps.density_jump > 1.0
    jump_mask &= _physical_jump_mask(jumps, gamma)
    jump_mask &= _upstream_floor_mask(jumps, config)

    # Either Mach estimate may pass; equality at min_mach is accepted.
    mach_mask = np.nan_to_num(mach_fields.mach_pressure, nan=-np.inf) >= config.min_mach
    mach_mask |= (
        np.nan_to_num(mach_fields.mach_temperature, nan=-np.inf) >= config.min_mach
    )
    return shock_zone_mask & jump_mask & mach_mask


def reduce_to_centers(
    shock_mask: np.ndarray,
    div_v: np.ndarray,
    mach_fields: MachFields,
    center_score_mode: str,
) -> tuple[np.ndarray, int]:
    """Choose one representative per nonperiodic six-neighbor component.

    This is component selection, not geometric sheet thinning. Simple-score
    ties follow scipy.ndimage.maximum_position; combined-score ties choose the
    first cell in C-order via argmax. Input arrays are not mutated.
    """

    # SciPy defaults to six-neighbor connectivity with no periodic seam joins.
    labels, count = ndimage.label(shock_mask)
    if count == 0:
        return np.zeros_like(shock_mask, dtype=bool), 0

    reduced = np.zeros_like(shock_mask, dtype=bool)
    simple_score_field = _simple_center_score_field(
        div_v, mach_fields, center_score_mode
    )
    if simple_score_field is not None:
        positions = ndimage.maximum_position(
            simple_score_field, labels=labels, index=np.arange(1, count + 1)
        )
        for position in positions:
            reduced[position] = True
        return reduced, count

    label_slices = ndimage.find_objects(labels)
    for label, component_slice in enumerate(label_slices, start=1):
        if component_slice is None:
            continue
        component = labels[component_slice] == label
        local_mach_fields = MachFields(
            mach_pressure=mach_fields.mach_pressure[component_slice],
            mach_temperature=mach_fields.mach_temperature[component_slice],
        )
        score = center_score(
            component, div_v[component_slice], local_mach_fields, center_score_mode
        )
        best_flat_index = int(np.argmax(score))
        best_local_index = np.unravel_index(best_flat_index, score.shape)
        best_index = tuple(
            component_slice[axis].start + best_local_index[axis]
            for axis in range(len(best_local_index))
        )
        reduced[best_index] = True
    return reduced, count


def center_score(
    component: np.ndarray,
    div_v: np.ndarray,
    mach_fields: MachFields,
    mode: str,
) -> np.ndarray:
    """Return a component-local center score for representative cell selection."""

    masked_out = -np.inf
    compression = np.where(component, -div_v, masked_out)
    finite_pressure = np.where(
        component, np.nan_to_num(mach_fields.mach_pressure, nan=masked_out), masked_out
    )
    finite_temperature = np.where(
        component,
        np.nan_to_num(mach_fields.mach_temperature, nan=masked_out),
        masked_out,
    )

    if mode == "compression":
        return compression
    if mode == "mach_pressure":
        return finite_pressure
    if mode == "mach_temperature":
        return finite_temperature

    combined_mach = np.maximum(finite_pressure, finite_temperature)
    compression_norm = _normalize_component_values(compression, component)
    mach_norm = _normalize_component_values(combined_mach, component)
    return np.where(component, compression_norm + mach_norm, masked_out)


def _normalize_component_values(
    values: np.ndarray, component: np.ndarray
) -> np.ndarray:
    """Normalize finite component values onto [0, 1], safely handling constants."""

    normalized = np.zeros(values.shape, dtype=float)
    finite_component = component & np.isfinite(values)
    if not np.any(finite_component):
        return normalized

    component_values = values[finite_component]
    vmin = float(np.min(component_values))
    vmax = float(np.max(component_values))
    if vmax <= vmin:
        normalized[finite_component] = 1.0
        return normalized

    normalized[finite_component] = (component_values - vmin) / (vmax - vmin)
    return normalized


def _simple_center_score_field(
    div_v: np.ndarray,
    mach_fields: MachFields,
    mode: str,
) -> np.ndarray | None:
    """Return a direct score field for modes that do not need per-component normalization."""

    if mode == "compression":
        return -div_v
    if mode == "mach_pressure":
        return np.nan_to_num(mach_fields.mach_pressure, nan=-np.inf)
    if mode == "mach_temperature":
        return np.nan_to_num(mach_fields.mach_temperature, nan=-np.inf)
    return None


def _physical_jump_mask(jumps: SampledJumps, gamma: float) -> np.ndarray:
    """Reject jump combinations that violate simple ideal-gas shock relations."""

    density_limit = (gamma + 1.0) / (gamma - 1.0)
    valid_density = np.isfinite(jumps.density_jump) & (
        jumps.density_jump <= density_limit
    )
    expected_temperature = np.full(jumps.pressure_jump.shape, np.nan, dtype=float)
    valid_expected = (
        np.isfinite(jumps.pressure_jump)
        & np.isfinite(jumps.density_jump)
        & (jumps.density_jump > 0.0)
    )
    expected_temperature[valid_expected] = (
        jumps.pressure_jump[valid_expected] / jumps.density_jump[valid_expected]
    )
    # Existing acceptance tolerance: T jump must agree with P jump / rho jump
    # within 25% of that expected ratio (zero absolute tolerance). This is
    # an empirical inherited tolerance; its original calibration is unknown.
    consistent_triplet = np.isclose(
        jumps.temperature_jump,
        expected_temperature,
        rtol=0.25,
        atol=0.0,
        equal_nan=False,
    )
    return valid_density & consistent_triplet


def _upstream_floor_mask(jumps: SampledJumps, config: ShockFinderConfig) -> np.ndarray:
    """Optionally reject jumps whose upstream sampled state falls below configured floors."""

    valid = np.ones(jumps.pressure_jump.shape, dtype=bool)
    if config.upstream_pressure_floor is not None:
        valid &= np.isfinite(jumps.pressure_up) & (
            jumps.pressure_up >= config.upstream_pressure_floor
        )
    if config.upstream_temperature_floor is not None:
        valid &= np.isfinite(jumps.temperature_up) & (
            jumps.temperature_up >= config.upstream_temperature_floor
        )
    if config.upstream_density_floor is not None:
        valid &= np.isfinite(jumps.density_up) & (
            jumps.density_up >= config.upstream_density_floor
        )
    return valid
