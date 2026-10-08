"""Adapter for extracted AthenaK metadata; native file parsing stays external."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import math
import numpy as np
from .mesh import AMRBlock, AMRMesh, Region


def amr_mesh_from_athenak(
    metadata: Mapping,
    load_region: Callable[[int, Region], Mapping[str, np.ndarray]],
    *,
    gamma: float = 5 / 3,
    boundaries=(("periodic", "periodic"),) * 3,
) -> AMRMesh:
    """Adapt extracted AthenaK metadata and a selective zyx primitive loader.

    Geometry entries are physical bounds, not centers. Infer levels from root
    spacing so meshes without active root blocks still normalize correctly.
    Pressure must be supplied explicitly; construction reads no field arrays.
    """
    try:
        root_shape = tuple(metadata[f"Nx{i}"] for i in (1, 2, 3))
        block_shape = tuple(metadata[f"nx{i}_mb"] for i in (1, 2, 3))
        domain_bounds = tuple(
            (metadata[f"x{i}min"], metadata[f"x{i}max"]) for i in (1, 2, 3)
        )
        logical_locations = np.asarray(metadata["mb_logical"])
        block_geometry = np.asarray(metadata["mb_geometry"], dtype=float)
        block_count = metadata["n_mbs"]
    except KeyError as exc:
        raise ValueError(f"missing AthenaK metadata: {exc}") from exc
    if (
        logical_locations.shape != (block_count, 4)
        or block_geometry.shape != (block_count, 6)
        or block_count < 1
    ):
        raise ValueError("invalid AthenaK metadata counts or shapes")
    if (
        not np.isfinite(logical_locations).all()
        or not np.equal(logical_locations, np.floor(logical_locations)).all()
    ):
        raise ValueError("AthenaK logical locations and levels must be integers")
    for i, n in enumerate(block_shape, 1):
        if metadata.get(f"nx{i}_out_mb", n) != n:
            raise ValueError("sliced AthenaK snapshots are unsupported")
    if "mb_index" in metadata:
        index = np.asarray(metadata["mb_index"])
        if (
            index.shape != (block_count, 6)
            or not np.array_equal(index[:, ::2], np.zeros((block_count, 3)))
            or not np.all(index[:, 1::2] == np.array(block_shape) - 1)
        ):
            raise ValueError("sliced AthenaK mb_index is unsupported")
    blocks = []
    level_offsets = set()
    for row in range(block_count):
        block_bounds = tuple(
            tuple(block_geometry[row, 2 * axis : 2 * axis + 2]) for axis in range(3)
        )
        if any(n <= 0 for n in block_shape + root_shape) or any(
            hi <= lo for lo, hi in block_bounds + domain_bounds
        ):
            raise ValueError("invalid AthenaK geometry")
        spacing_ratios = [
            ((domain_bounds[axis][1] - domain_bounds[axis][0]) / root_shape[axis])
            / ((block_bounds[axis][1] - block_bounds[axis][0]) / block_shape[axis])
            for axis in range(3)
        ]
        if any(
            not np.isfinite(spacing_ratio) or spacing_ratio <= 0
            for spacing_ratio in spacing_ratios
        ):
            raise ValueError("invalid AthenaK refinement")
        # Physical root spacing defines level zero, even if every root is refined.
        level = round(math.log2(spacing_ratios[0]))
        if level < 0 or not np.allclose(spacing_ratios, 2.0**level, rtol=1e-10, atol=0):
            raise ValueError("unsupported AthenaK refinement ratio")
        level_offsets.add(int(logical_locations[row, 3]) - level)
        blocks.append(
            AMRBlock(
                str(row),
                level,
                tuple(int(x) for x in logical_locations[row, :3]),
                block_shape,
                block_bounds,
            )
        )
    if len(level_offsets) != 1:
        raise ValueError("inconsistent AthenaK refinement levels")

    def normalized(block_id, slices):
        """Reverse spatial axes and slices only; preserve velocity-component names."""
        values = load_region(int(block_id), slices[::-1])
        return {
            name: np.asarray(value).transpose(2, 1, 0) for name, value in values.items()
        }

    return AMRMesh(blocks, domain_bounds, root_shape, normalized, gamma, boundaries)
