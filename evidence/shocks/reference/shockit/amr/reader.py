"""Bounded primitive region access, boundary mapping and AMR ghost filling."""

from __future__ import annotations

from collections import OrderedDict
from collections.abc import Iterator
from itertools import product
import math

import numpy as np

from .mesh import AMRBlock, AMRMesh, AMRExecutionOptions, Region, validate_slices

FIELDS = ("rho", "pressure", "vx", "vy", "vz")


def minmod(left, right):
    """Keep the smaller same-sign slope; return zero across a slope sign change."""
    return np.where(
        np.sign(left) == np.sign(right),
        np.sign(left) * np.minimum(np.abs(left), np.abs(right)),
        0.0,
    )


class AMRReader:
    """Fill xyz tiles from primitive cell averages using bounded source reads.

    Internal sample matrices have shape (cells, 5), with columns in FIELDS
    order. Coordinate indices are global at the requested level unless a
    method explicitly accepts block-local indices.
    """

    def __init__(self, mesh: AMRMesh, options: AMRExecutionOptions | None = None):
        """Set source-tile and cache limits without invoking the caller loader."""
        self.mesh = mesh
        self.options = options or AMRExecutionOptions()
        self.read_shape = np.array(self.options.fitted_shape(), dtype=np.int64)
        self.cache_limit = min(
            self.options.cache_bytes, self.options.memory_budget_bytes // 4
        )
        self._cache = OrderedDict()
        self.cache_nbytes = 0

    def _load(self, block_id: str, slices: Region) -> dict[str, np.ndarray]:
        """Load, validate and own one interior region; cache only entries that fit."""
        cache_key = (block_id, tuple((s.start, s.stop) for s in slices))
        if cache_key in self._cache:
            self._cache.move_to_end(cache_key)
            return self._cache[cache_key]
        try:
            loaded_fields = self.mesh.load_region(block_id, slices)
        except Exception as exc:
            raise OSError(f"block {block_id} region {slices}: {exc}") from exc
        region_shape = tuple(s.stop - s.start for s in slices)
        fields = {}
        for name in FIELDS:
            if name not in loaded_fields:
                raise ValueError(
                    f"block {block_id} region {slices}: missing field {name}"
                )
            try:
                array = np.asarray(loaded_fields[name], dtype=float)
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid field {name} in block {block_id}") from exc
            if array.shape != region_shape or not np.isfinite(array).all():
                raise ValueError(
                    f"invalid shape or non-finite {name} in block {block_id} region {slices}"
                )
            if name in ("rho", "pressure") and np.any(array <= 0):
                raise ValueError(
                    f"{name} must be strictly positive in block {block_id}"
                )
            # Own the cached data: caller buffers may be reused on the next read.
            fields[name] = array.copy()
        entry_bytes = sum(a.nbytes for a in fields.values())
        if entry_bytes <= self.cache_limit:
            while self.cache_nbytes + entry_bytes > self.cache_limit:
                _, evicted_fields = self._cache.popitem(last=False)
                self.cache_nbytes -= sum(a.nbytes for a in evicted_fields.values())
            self._cache[cache_key] = fields
            self.cache_nbytes += entry_bytes
        return fields

    def _fetch(self, block: AMRBlock, indices: np.ndarray) -> np.ndarray:
        """Gather cells in bounded source tiles, never their global bounding box."""
        samples = np.empty((len(indices), 5), float)
        tile_keys, tile_groups = np.unique(
            indices // self.read_shape, axis=0, return_inverse=True
        )
        for tile_number, tile_key in enumerate(tile_keys):
            positions = np.flatnonzero(tile_groups == tile_number)
            cell_indices = indices[positions]
            source_start = cell_indices.min(axis=0)
            source_stop = cell_indices.max(axis=0) + 1
            source_slices = tuple(
                slice(int(a), int(b)) for a, b in zip(source_start, source_stop)
            )
            fields = self._load(block.block_id, source_slices)
            local_indices = tuple((cell_indices - source_start).T)
            for field_index, name in enumerate(FIELDS):
                samples[positions, field_index] = fields[name][local_indices]
        return samples

    def _boundary(
        self, level: int, indices: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        """Map level-global cell indices to the domain and return per-component velocity signs.

        Repeated reflections use a period of twice the axis length. Scalars are
        unchanged; only the velocity normal to each crossed face reverses sign.
        """
        mapped_indices = indices.copy()
        velocity_parity = np.ones((len(indices), 3))
        level_shape = np.array(self.mesh.root_shape, dtype=np.int64) * 2**level
        for axis, axis_cells in enumerate(level_shape):
            low_boundary, high_boundary = self.mesh.boundaries[axis]
            axis_indices = mapped_indices[:, axis]
            if low_boundary == "periodic":
                axis_indices[:] = axis_indices % axis_cells
                continue
            if low_boundary == high_boundary == "reflecting":
                period_position = axis_indices % (2 * axis_cells)
                reflected = period_position >= axis_cells
                velocity_parity[reflected, axis] = -1
                axis_indices[:] = np.where(
                    reflected, 2 * axis_cells - 1 - period_position, period_position
                )
            else:
                # With at least one outflow face, at most one reflection is
                # followed by clamping to the opposite physical boundary.
                if low_boundary == "reflecting":
                    crossed_face = axis_indices < 0
                    axis_indices[crossed_face] = -axis_indices[crossed_face] - 1
                    velocity_parity[crossed_face, axis] *= -1
                if high_boundary == "reflecting":
                    crossed_face = axis_indices >= axis_cells
                    axis_indices[crossed_face] = (
                        2 * axis_cells - 1 - axis_indices[crossed_face]
                    )
                    velocity_parity[crossed_face, axis] *= -1
                axis_indices[:] = np.clip(axis_indices, 0, axis_cells - 1)
        return mapped_indices, velocity_parity

    def _donors(self, level: int, indices: np.ndarray) -> Iterator[tuple]:
        """Group virtual cells by their active same-level or coarser donor.

        A final group with no block contains cells requiring restriction.
        Integer level coordinates avoid a finest-grid ownership allocation.
        """
        pending = np.ones(len(indices), bool)
        block_shape = np.array(self.mesh.block_shape, dtype=np.int64)
        for source_level in range(min(level, self.mesh.max_level) + 1):
            selected = np.flatnonzero(pending)
            if not len(selected):
                break
            donor_indices = indices[selected] // 2 ** (level - source_level)
            block_keys, block_groups = np.unique(
                donor_indices // block_shape, axis=0, return_inverse=True
            )
            for group_index, block_key in enumerate(block_keys):
                block = self.mesh.leaf_at(source_level, tuple(map(int, block_key)))
                if block is None:
                    continue
                positions = selected[block_groups == group_index]
                global_donor = indices[positions] // 2 ** (level - source_level)
                pending[positions] = False
                yield positions, block, global_donor
        unresolved = np.flatnonzero(pending)
        if len(unresolved):
            yield unresolved, None, None

    def _prolong(
        self,
        level: int,
        source_level: int,
        indices: np.ndarray,
        global_donor: np.ndarray,
        values: np.ndarray,
    ) -> np.ndarray:
        """Limited linear primitive reconstruction; stencil donors are averages.

        `indices` are requested-level global cell indices, `global_donor` is in
        source-level global cells. `values` has axes (points, five primitives).
        Slopes below have axes (points, xyz slope directions, five primitives)
        and units of primitive per donor cell, not per physical length.
        Slope-neighbor averages disable reconstruction; see _average termination.
        A common limiter couples all five fields, including thermodynamics to
        velocity extrema. Primitive averages, not momentum/energy, are preserved.
        """
        slopes = []
        lower = values.copy()
        upper = values.copy()
        for axis in range(3):
            offset = np.zeros(3, dtype=np.int64)
            offset[axis] = 1
            minus_values = self._average(source_level, global_donor - offset, False)
            plus_values = self._average(source_level, global_donor + offset, False)
            slopes.append(minmod(values - minus_values, plus_values - values))
            lower = np.minimum(lower, np.minimum(minus_values, plus_values))
            upper = np.maximum(upper, np.maximum(minus_values, plus_values))
        slopes = np.stack(slopes, axis=1)
        # Axial extrema alone clip valid diagonal affine gradients at corners.
        for offset in product((-1, 0, 1), repeat=3):
            if sum(abs(x) for x in offset) <= 1:
                continue
            neighbor = self._average(
                source_level, global_donor + np.array(offset), False
            )
            lower = np.minimum(lower, neighbor)
            upper = np.maximum(upper, neighbor)
        # Bound every reconstructed point within the parent cell using all 26
        # neighbors. One limiter per cell scales all primitive slopes together,
        # preserving positive density/pressure and the parent average.
        # For every parent-local offset |offset_axis| <= 1/2, triangle inequality
        # gives |sum(slope * offset)| <= 1/2 * sum(abs(slope)). Taking the minimum
        # ratio over ALL five primitives bounds each by its 26-neighbor extrema.
        # Positive rho/P neighbor minima imply positive reconstructed rho/P.
        # Symmetric child offsets average to zero, so one common parent limiter
        # preserves primitive child means (not conservative momentum or energy).
        max_excursion = 0.5 * np.abs(slopes).sum(axis=1)
        extrema_room = np.minimum(values - lower, upper - values)
        slope_ratios = np.ones_like(extrema_room)
        np.divide(
            extrema_room, max_excursion, out=slope_ratios, where=max_excursion > 0
        )
        limiter = np.minimum(1.0, np.maximum(0.0, slope_ratios.min(axis=1)))
        child_offset = (indices + 0.5) / 2 ** (level - source_level) - (
            global_donor + 0.5
        )
        return (
            values
            + np.sum(slopes * child_offset[:, :, None], axis=1) * limiter[:, None]
        )

    def _restrict(
        self, level: int, indices: np.ndarray, reconstruct: bool
    ) -> np.ndarray:
        """Stream eight equal-volume child contributions in a fixed order."""
        if level >= self.mesh.max_level:
            raise ValueError("uncovered AMR cells")
        average = np.zeros((len(indices), 5))
        for bits in product((0, 1), repeat=3):
            average += (
                self._average(level + 1, indices * 2 + np.array(bits), reconstruct) / 8
            )
        return average

    def _average(
        self, level: int, indices: np.ndarray, reconstruct: bool = True
    ) -> np.ndarray:
        """Map level-global cells, dispatch donor transfers, apply boundary parity.

        Call/termination table:
        exact active donor -> _fetch -> return;
        coarser donor, reconstruct=True -> _prolong -> slope queries with False;
        coarser donor, reconstruct=False -> _fetch parent average -> return;
        finer donors -> _restrict at level+1 -> active leaves by mesh.max_level.

        Thus restriction strictly approaches leaves and slope queries cannot
        recursively reconstruct other slopes. Boundary mapping changes indices,
        not levels. Donor indices become block-local only immediately before fetch.
        """
        indices, velocity_parity = self._boundary(level, indices)
        samples = np.empty((len(indices), 5), float)
        block_shape = np.array(self.mesh.block_shape, dtype=np.int64)
        for positions, block, global_donor in self._donors(level, indices):
            if block is None:
                samples[positions] = self._restrict(
                    level, indices[positions], reconstruct
                )
                continue
            local_indices = (
                global_donor - np.array(block.logical_location) * block_shape
            )
            values = self._fetch(block, local_indices)
            if block.level < level and reconstruct:
                values = self._prolong(
                    level, block.level, indices[positions], global_donor, values
                )
            samples[positions] = values
        # FIELDS starts with two scalars, followed by vx, vy, vz.
        samples[:, 2:] *= velocity_parity
        return samples

    def read_with_halo(
        self, block_id: str, slices: Region, halo: int
    ) -> dict[str, np.ndarray]:
        """Return primitive arrays for a block-local core plus level-local ghost cells.

        Transfer primitives before the physics pass derives temperature and entropy.
        Boundary wrapping uses the physical domain, never the tile or block extent.
        """
        block = self.mesh.block(block_id)
        slices = validate_slices(slices, block.shape)
        if not isinstance(halo, (int, np.integer)) or halo < 0:
            raise ValueError("halo must be a nonnegative integer")
        output_shape = tuple(s.stop - s.start + 2 * halo for s in slices)
        if math.prod(output_shape) * 2048 > self.options.memory_budget_bytes - min(
            self.cache_limit, self.options.memory_budget_bytes // 4
        ):
            raise ValueError("memory budget exceeded by requested halo region")
        # Direct interiors avoid coordinate work and preserve source values exactly.
        if halo == 0 and all(n <= s for n, s in zip(output_shape, self.read_shape)):
            return {
                name: array.copy()
                for name, array in self._load(block_id, slices).items()
            }
        fields = {name: np.empty(output_shape) for name in FIELDS}
        global_start = np.array([s.start - halo for s in slices]) + np.array(
            block.logical_location
        ) * np.array(block.shape)
        flat_fields = [array.ravel() for array in fields.values()]
        # Bound reconstruction batches independently of the compute tile size.
        batch_size = max(1, min(8192, self.options.memory_budget_bytes // 32768))
        for begin in range(0, math.prod(output_shape), batch_size):
            stop = min(begin + batch_size, math.prod(output_shape))
            coordinates = (
                np.array(
                    np.unravel_index(
                        np.arange(begin, stop, dtype=np.int64), output_shape
                    )
                ).T
                + global_start
            )
            values = self._average(block.level, coordinates)
            for field_index, flat_field in enumerate(flat_fields):
                flat_field[begin:stop] = values[:, field_index]
        return fields
