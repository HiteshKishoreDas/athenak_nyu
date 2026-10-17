"""Shared local shock physics and execution contracts.

No runner or public dispatcher imports: uniform, chunked and AMR all depend here.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
from typing import Literal

import numpy as np

from .config import ShockFinderConfig
from .fields import FluidCube
from .derived import (
    DerivedFields,
    MachFields,
    NormalFields,
    compute_derived_fields,
    compute_gradients,
    compute_mach_fields,
    compute_normals,
)
from .masks import build_final_shock_mask, build_shock_zone_mask
from .sampling import SampledJumps, sample_jumps

ProgressCallback = Callable[[str], None]
ProgressReporter = ProgressCallback | Literal[True, False] | None


@dataclass
class FinderPass:
    """Intermediate physics arrays before output cropping or center reduction."""

    derived: DerivedFields
    normals: NormalFields
    shock_zone_mask: np.ndarray
    jumps: SampledJumps
    mach_fields: MachFields
    full_shock_mask: np.ndarray


def analyze_cube_pass(
    cube: FluidCube,
    config: ShockFinderConfig,
    *,
    report: ProgressCallback | None = None,
) -> FinderPass:
    """Run the local shock-finder physics on one cube."""

    config.validate()
    total_steps = 9 if config.reduce_to_centers else 8
    if report is not None:
        report(f"[1/{total_steps}] Computing derived fields")
    derived = compute_derived_fields(cube)
    if report is not None:
        report(f"[2/{total_steps}] Computing gradients")
    gradients = compute_gradients(cube, derived.temperature, derived.entropy)
    if report is not None:
        report(f"[3/{total_steps}] Building shock-zone mask")
    shock_zone_mask = build_shock_zone_mask(derived, gradients, config)
    if report is not None:
        report(f"[4/{total_steps}] Computing shock normals")
    normals = compute_normals(gradients, config.normal_field)
    if report is not None:
        report(f"[5/{total_steps}] Sampling upstream and downstream jumps")
    jumps = sample_jumps(
        pressure=cube.pressure,
        temperature=derived.temperature,
        rho=cube.rho,
        nx=normals.normal_x,
        ny=normals.normal_y,
        nz=normals.normal_z,
        width=config.shock_width_cells,
        method=config.sampling_method,
        candidate_mask=shock_zone_mask,
        physical_distance=config.sampling_distance,
        spacing=(cube.dx, cube.dy, cube.dz),
    )
    if report is not None:
        report(f"[6/{total_steps}] Computing Mach fields")
    mach_fields = compute_mach_fields(
        jumps,
        cube.gamma,
        mach_max=config.mach_max,
        candidate_mask=shock_zone_mask,
    )
    if report is not None:
        report(f"[7/{total_steps}] Applying final shock filters")
    full_shock_mask = build_final_shock_mask(
        shock_zone_mask,
        jumps,
        mach_fields,
        config,
        cube.gamma,
    )
    return FinderPass(
        derived=derived,
        normals=normals,
        shock_zone_mask=shock_zone_mask,
        jumps=jumps,
        mach_fields=mach_fields,
        full_shock_mask=full_shock_mask,
    )


def build_reporter(progress: ProgressReporter) -> ProgressCallback:
    """Resolve print, silence or a caller callback; None preserves printed progress."""
    if progress is True or progress is None:
        return print
    if progress is False:
        return _noop_report
    if not callable(progress):
        raise TypeError("progress must be bool, None, or callable.")
    return progress


def _noop_report(_: str) -> None:
    """Consume progress messages when reporting is disabled."""
    return None


def required_halo(config: ShockFinderConfig, spacing=None) -> int:
    """Bound gradients and both endpoint/interpolation supports in cell units.

    A physical unit-vector component is at most one, so ceil(distance/min(dx))
    bounds every endpoint corner. AMR may use the finest spacing for all tiles.
    """
    if config.sampling_distance is None:
        return max(1, config.shock_width_cells)
    if spacing is None or len(spacing) != 3 or any(not math.isfinite(step) or step <= 0 for step in spacing):
        raise ValueError("physical sampling requires three finite positive spacing values")
    cells = config.sampling_distance / min(spacing)
    if not math.isfinite(cells):
        raise ValueError("sampling_distance exceeds representable halo size")
    return max(1, math.ceil(cells))


def pass_fields(
    local: FinderPass, crop=None, *, shock_mask=None
) -> dict[str, np.ndarray]:
    """Map a local pass to the common result schema; preserve final-mask ownership.

    With no final mask supplied, return an independent copy of the full mask.
    Cropping happens before that copy so streamed runners copy only the core.
    """
    values = dict(
        shock_zone_mask=local.shock_zone_mask,
        full_shock_mask=local.full_shock_mask,
        mach_temperature=local.mach_fields.mach_temperature,
        mach_pressure=local.mach_fields.mach_pressure,
        compression=local.derived.compression,
        div_v=local.derived.div_v,
        temperature=local.derived.temperature,
        entropy=local.derived.entropy,
        temperature_jump=local.jumps.temperature_jump,
        pressure_jump=local.jumps.pressure_jump,
        density_jump=local.jumps.density_jump,
        normal_x=local.normals.normal_x,
        normal_y=local.normals.normal_y,
        normal_z=local.normals.normal_z,
    )
    if crop is not None:
        values = {name: array[crop] for name, array in values.items()}
    values["shock_mask"] = (
        values["full_shock_mask"].copy() if shock_mask is None else shock_mask
    )
    return values
