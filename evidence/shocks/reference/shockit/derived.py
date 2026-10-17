"""Derived fields and gradient helpers used by the shock finder."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fields import FluidCube
from .gradients import divergence, gradient, magnitude
from .sampling import SampledJumps
from .mach import mach_from_temperature_jump


@dataclass
class DerivedFields:
    """Thermodynamic proxy fields used by the shock finder."""

    temperature: np.ndarray
    entropy: np.ndarray
    div_v: np.ndarray
    compression: np.ndarray


@dataclass
class GradientFields:
    """Spatial gradients: t=temperature, s=entropy, p=pressure, r=density.

    Each suffix x/y/z names the derivative axis, not a velocity component.
    """

    grad_tx: np.ndarray
    grad_ty: np.ndarray
    grad_tz: np.ndarray
    grad_sx: np.ndarray
    grad_sy: np.ndarray
    grad_sz: np.ndarray
    grad_px: np.ndarray
    grad_py: np.ndarray
    grad_pz: np.ndarray
    grad_rx: np.ndarray
    grad_ry: np.ndarray
    grad_rz: np.ndarray
    grad_t_mag: np.ndarray
    grad_p_mag: np.ndarray


@dataclass
class NormalFields:
    """Normalized shock-normal field components."""

    normal_x: np.ndarray
    normal_y: np.ndarray
    normal_z: np.ndarray


@dataclass
class MachFields:
    """Mach fields inferred from pressure and temperature jumps."""

    mach_pressure: np.ndarray
    mach_temperature: np.ndarray


def compute_derived_fields(cube: FluidCube) -> DerivedFields:
    """Return float64 xyz proxies without mutating primitives.

    T=P/rho; S=log(P)-gamma*log(rho); div_v has inverse-time units.
    Absolute entropy depends on caller units, but its gradient does not depend
    on a constant reference offset. Reject nonfinite results before masking.
    """

    # Ideal-gas proxies omit constant unit factors: T=P/rho and S=log(P/rho^gamma).
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        temperature = cube.pressure / cube.rho
        entropy = np.log(cube.pressure) - cube.gamma * np.log(cube.rho)
        div_v = divergence(cube.vx, cube.vy, cube.vz, cube.dx, cube.dy, cube.dz)
    _require_finite("derived thermodynamics", temperature=temperature, entropy=entropy)
    _require_finite("derived velocity", div_v=div_v)
    compression = -div_v
    return DerivedFields(
        temperature=temperature,
        entropy=entropy,
        div_v=div_v,
        compression=compression,
    )


def compute_gradients(
    cube: FluidCube, temperature: np.ndarray, entropy: np.ndarray
) -> GradientFields:
    """Compute periodic centered gradients of the relevant fields."""

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        grad_tx, grad_ty, grad_tz = gradient(temperature, cube.dx, cube.dy, cube.dz)
        grad_sx, grad_sy, grad_sz = gradient(entropy, cube.dx, cube.dy, cube.dz)
        grad_px, grad_py, grad_pz = gradient(cube.pressure, cube.dx, cube.dy, cube.dz)
        grad_rx, grad_ry, grad_rz = gradient(cube.rho, cube.dx, cube.dy, cube.dz)
    fields = GradientFields(
        grad_tx=grad_tx,
        grad_ty=grad_ty,
        grad_tz=grad_tz,
        grad_sx=grad_sx,
        grad_sy=grad_sy,
        grad_sz=grad_sz,
        grad_px=grad_px,
        grad_py=grad_py,
        grad_pz=grad_pz,
        grad_rx=grad_rx,
        grad_ry=grad_ry,
        grad_rz=grad_rz,
        grad_t_mag=magnitude(grad_tx, grad_ty, grad_tz),
        grad_p_mag=magnitude(grad_px, grad_py, grad_pz),
    )
    _require_finite("gradient", **vars(fields))
    return fields


def compute_normals(gradients: GradientFields, normal_field: str) -> NormalFields:
    """Compute unit shock normals from the requested gradient field."""

    if normal_field == "temperature":
        normal_x_raw = gradients.grad_tx
        normal_y_raw = gradients.grad_ty
        normal_z_raw = gradients.grad_tz
        normal_mag = gradients.grad_t_mag
    else:
        normal_x_raw = gradients.grad_px
        normal_y_raw = gradients.grad_py
        normal_z_raw = gradients.grad_pz
        normal_mag = gradients.grad_p_mag

    # A zero gradient keeps a finite zero normal instead of dividing by zero.
    safe_mag = np.where(normal_mag > 0.0, normal_mag, 1.0)
    fields = NormalFields(
        normal_x=normal_x_raw / safe_mag,
        normal_y=normal_y_raw / safe_mag,
        normal_z=normal_z_raw / safe_mag,
    )
    _require_finite("normal", **vars(fields))
    return fields


def _require_finite(context: str, **arrays: np.ndarray) -> None:
    """Reject invalid computed fields before masks can hide numeric-range failures."""
    for name, array in arrays.items():
        if not np.all(np.isfinite(array)):
            raise ValueError(
                f"Nonfinite {context} field {name} for shape {array.shape}; "
                "rescale primitive units or grid spacing to a representable range."
            )


def compute_mach_fields(
    jumps: SampledJumps,
    gamma: float,
    mach_max: float,
    candidate_mask: np.ndarray | None = None,
) -> MachFields:
    """Compute Mach estimates from pressure and temperature jumps."""

    shape = jumps.pressure_jump.shape
    mach_pressure = np.full(shape, np.nan, dtype=float)
    valid_pressure = jumps.pressure_jump > 1.0
    if candidate_mask is not None:
        valid_pressure &= candidate_mask
    if np.any(valid_pressure):
        mach_pressure[valid_pressure] = np.sqrt(
            ((jumps.pressure_jump[valid_pressure] * (gamma + 1.0)) + (gamma - 1.0))
            / (2.0 * gamma)
        )

    # mach_max belongs to temperature inversion; pressure inversion is uncapped.
    mach_temperature = np.full(shape, np.nan, dtype=float)
    valid_temperature = jumps.temperature_jump > 1.0
    if candidate_mask is not None:
        valid_temperature &= candidate_mask
    if np.any(valid_temperature):
        values = [
            mach_from_temperature_jump(value, gamma, mach_max=mach_max)
            for value in jumps.temperature_jump[valid_temperature]
        ]
        mach_temperature[valid_temperature] = np.asarray(values, dtype=float)

    return MachFields(mach_pressure=mach_pressure, mach_temperature=mach_temperature)
