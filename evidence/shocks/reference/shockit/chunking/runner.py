"""Explicit chunked runner for streamed shock finding."""

from __future__ import annotations


from copy import deepcopy

import numpy as np

from ..provenance import build_provenance
from .layout import require_empty_output
from ..config import ShockFinderConfig
from ..fields import FluidCube
from ..local_pass import (
    analyze_cube_pass,
    pass_fields,
    required_halo,
    build_reporter,
    ProgressCallback,
    ProgressReporter,
)
from ..result import RESULT_FIELD_NAMES
from .centers import finalize_chunked_shock_outputs
from .reader import ChunkedInputStore
from .writer import ChunkedOutputStore

OUTPUT_FIELDS = RESULT_FIELD_NAMES


def run_chunked_shock_finder(
    input_store: ChunkedInputStore,
    output_store: ChunkedOutputStore,
    *,
    config: ShockFinderConfig,
    dx: float = 1.0,
    dy: float = 1.0,
    dz: float = 1.0,
    gamma: float = 5.0 / 3.0,
    progress: ProgressReporter = True,
) -> dict[str, object]:
    """Run periodic xyz core/halo analysis, publishing a summary only on success.

    Input cores remain lazy. Exact medians retain detected Mach values. Legacy
    simple-score centers allocate at least 28 bytes/global cell before temporaries;
    this is a lower bound, not a peak RSS estimate. Caller inputs are not modified.
    """

    config.validate()
    config = deepcopy(config)
    report = build_reporter(progress)
    require_empty_output(output_store.root)

    layout = input_store.layout
    halo = required_halo(config, (dx, dy, dz))
    total_steps = 2 if config.reduce_to_centers else 1
    report(
        f"[1/{total_steps}] Processing {len(layout.specs)} chunk files with halo width {halo}"
    )

    primitive_readers = {
        name: input_store.field_reader(name)
        for name in ("rho", "pressure", "vx", "vy", "vz")
    }
    shock_zone_cells = 0
    full_shock_cells = 0
    full_mach_pressure_values: list[np.ndarray] = []
    full_mach_temperature_values: list[np.ndarray] = []

    for chunk_index, spec in enumerate(layout.specs, start=1):
        report(
            f"[1/{total_steps}] Chunk {chunk_index}/{len(layout.specs)} {spec.start}->{spec.stop}"
        )
        local_cube = FluidCube(
            rho=primitive_readers["rho"].read_with_halo(spec, halo),
            pressure=primitive_readers["pressure"].read_with_halo(spec, halo),
            vx=primitive_readers["vx"].read_with_halo(spec, halo),
            vy=primitive_readers["vy"].read_with_halo(spec, halo),
            vz=primitive_readers["vz"].read_with_halo(spec, halo),
            dx=dx,
            dy=dy,
            dz=dz,
            gamma=gamma,
        )
        local_pass = analyze_cube_pass(local_cube, config)
        core_crop = tuple(slice(halo, halo + size) for size in spec.core_shape)
        output_values = pass_fields(local_pass, core_crop)
        for field_name in OUTPUT_FIELDS:
            output_store.write_core(field_name, spec, output_values[field_name])

        full_mask = output_values["full_shock_mask"].astype(bool)
        shock_zone_cells += int(np.count_nonzero(output_values["shock_zone_mask"]))
        full_shock_cells += int(np.count_nonzero(full_mask))
        full_mach_pressure_values.append(output_values["mach_pressure"][full_mask])
        full_mach_temperature_values.append(
            output_values["mach_temperature"][full_mask]
        )

    base_summary = _base_summary(
        shape=layout.shape,
        shock_zone_cells=shock_zone_cells,
        full_shock_cells=full_shock_cells,
        full_mach_pressure_values=full_mach_pressure_values,
        full_mach_temperature_values=full_mach_temperature_values,
        config=config,
        gamma=gamma,
    )
    base_summary["provenance"] = build_provenance(
        config, gamma=gamma, spacing=(dx, dy, dz),
        source_metadata=getattr(input_store, "metadata", None),
    )
    base_summary["center_memory_lower_bound_bytes"] = 28 * int(np.prod(layout.shape))
    base_summary["provenance"]["center_memory"] = {
        "simple_score_lower_bound_bytes": base_summary["center_memory_lower_bound_bytes"],
        "bytes_per_global_cell": 28, "peak_rss": False,
        "exact_medians": "retain detected Mach values",
    }
    if not config.reduce_to_centers:
        output_store.write_summary(base_summary)
        return base_summary
    return finalize_chunked_shock_outputs(
        output_store,
        config,
        base_summary,
        progress=progress,
    )


def _base_summary(
    *,
    shape: tuple[int, int, int],
    shock_zone_cells: int,
    full_shock_cells: int,
    full_mach_pressure_values: list[np.ndarray],
    full_mach_temperature_values: list[np.ndarray],
    config: ShockFinderConfig,
    gamma: float,
) -> dict[str, object]:
    """Summarize unreduced chunk populations; medians require retained Mach samples."""
    total_cells = int(np.prod(shape))
    mach_pressure = _concat(full_mach_pressure_values)
    mach_temperature = _concat(full_mach_temperature_values)
    return {
        "grid_shape": shape,
        "total_cells": total_cells,
        "shock_zone_cells": shock_zone_cells,
        "shock_cells": full_shock_cells,
        "shock_fraction": full_shock_cells / total_cells,
        "full_shock_cells": full_shock_cells,
        "full_shock_fraction": full_shock_cells / total_cells,
        "n_connected_components": None,
        "center_topology": "nonperiodic six-connected components",
        "component_count_topology": "periodic six-connected components when computed",
        "mask_semantics": "shock_mask matches full_shock_mask; center reduction disabled",
        "mach_pressure_min": _finite_stat(mach_pressure, np.min),
        "mach_pressure_median": _finite_stat(mach_pressure, np.median),
        "mach_pressure_max": _finite_stat(mach_pressure, np.max),
        "mach_temperature_min": _finite_stat(mach_temperature, np.min),
        "mach_temperature_median": _finite_stat(mach_temperature, np.median),
        "mach_temperature_max": _finite_stat(mach_temperature, np.max),
        "full_mach_pressure_min": _finite_stat(mach_pressure, np.min),
        "full_mach_pressure_median": _finite_stat(mach_pressure, np.median),
        "full_mach_pressure_max": _finite_stat(mach_pressure, np.max),
        "full_mach_temperature_min": _finite_stat(mach_temperature, np.min),
        "full_mach_temperature_median": _finite_stat(mach_temperature, np.median),
        "full_mach_temperature_max": _finite_stat(mach_temperature, np.max),
        "min_mach": config.min_mach,
        "gamma": gamma,
        "shock_width_cells": config.shock_width_cells,
        "reduce_to_centers": config.reduce_to_centers,
        "center_score": config.center_score,
        "sampling_method": config.sampling_method,
    }


def _concat(values: list[np.ndarray]) -> np.ndarray:
    """Combine stored populations or return an empty float array."""
    non_empty = [array for array in values if array.size > 0]
    if not non_empty:
        return np.array([], dtype=float)
    return np.concatenate(non_empty)


def _finite_stat(values: np.ndarray, reducer) -> float:
    """Apply a statistic to finite values, returning NaN if none remain."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(reducer(finite))
