"""Cartesian leaf-block metadata and bounded, lazy compute tiling."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from itertools import product
import math
from numbers import Integral

import numpy as np

Region = tuple[slice, slice, slice]
RegionLoader = Callable[[str, Region], Mapping[str, np.ndarray]]


def _ints(values, name, positive=False):
    """Validate three axis counts or logical indices, excluding booleans."""
    if len(values) != 3 or any(
        isinstance(x, bool) or not isinstance(x, Integral) or x < (1 if positive else 0)
        for x in values
    ):
        raise ValueError(
            f'{name} must contain three {"positive" if positive else "nonnegative"} integers'
        )
    return tuple(int(x) for x in values)


@dataclass(frozen=True)
class AMRBlock:
    """One active leaf: level-local block indices, xyz cell counts and physical bounds."""

    block_id: str
    level: int
    logical_location: tuple[int, int, int]
    shape: tuple[int, int, int]
    bounds: tuple[tuple[float, float], ...]

    def __post_init__(self):
        """Validate and normalize immutable block metadata; allocate no cell arrays."""
        if not isinstance(self.block_id, str) or not self.block_id:
            raise ValueError("block_id must be a nonempty string")
        if (
            isinstance(self.level, bool)
            or not isinstance(self.level, Integral)
            or not 0 <= self.level <= 30
        ):
            raise ValueError("level must be an integer from 0 to 30")
        object.__setattr__(self, "level", int(self.level))
        object.__setattr__(self, "shape", _ints(self.shape, "shape", True))
        object.__setattr__(
            self, "logical_location", _ints(self.logical_location, "logical_location")
        )
        bounds = _bounds(self.bounds)
        object.__setattr__(self, "bounds", bounds)

    @property
    def spacing(self):
        """Physical cell widths (dx, dy, dz) at this block's refinement level."""
        return tuple((hi - lo) / n for (lo, hi), n in zip(self.bounds, self.shape))


def _bounds(bounds):
    """Normalize finite physical bounds in x, y, z order."""
    if len(bounds) != 3 or any(
        len(pair) != 2 or not all(np.isfinite(pair)) or pair[1] <= pair[0]
        for pair in bounds
    ):
        raise ValueError("bounds must be three finite increasing pairs")
    return tuple(tuple(float(x) for x in pair) for pair in bounds)


class AMRMesh:
    """Complete, non-overlapping Cartesian leaves with factor-two, 2:1-balanced refinement.

    The loader accepts block-local xyz interior slices and returns rho, pressure,
    vx, vy, vz. Construction validates metadata only; fields stay caller-owned.
    """

    def __init__(
        self,
        blocks,
        domain_bounds,
        root_shape,
        load_region: RegionLoader,
        gamma=5 / 3,
        boundaries=(("periodic", "periodic"),) * 3,
    ):
        """Validate geometry, ownership, coverage and all neighbor contacts without loading fields."""
        self.blocks = tuple(blocks)
        if not self.blocks:
            raise ValueError("mesh requires active leaf blocks")
        self.domain_bounds = _bounds(domain_bounds)
        self.root_shape = _ints(root_shape, "root_shape", True)
        if not np.isfinite(gamma) or gamma <= 1:
            raise ValueError("gamma must be finite and greater than 1")
        self.gamma = float(gamma)
        if not callable(load_region):
            raise TypeError("load_region must be callable")
        self.load_region = load_region
        if len(boundaries) != 3 or any(
            len(face_pair) != 2
            or any(
                boundary not in ("periodic", "outflow", "reflecting")
                for boundary in face_pair
            )
            or (("periodic" in face_pair) and face_pair[0] != face_pair[1])
            for face_pair in boundaries
        ):
            raise ValueError(
                "boundary faces must be periodic pairs, outflow, or reflecting"
            )
        self.boundaries = tuple(tuple(face_pair) for face_pair in boundaries)
        self.block_shape = self.blocks[0].shape
        if any(
            axis_count % block_cells
            for axis_count, block_cells in zip(self.root_shape, self.block_shape)
        ):
            raise ValueError("root_shape must be divisible by block shape")
        self.root_blocks = tuple(
            axis_count // block_cells
            for axis_count, block_cells in zip(self.root_shape, self.block_shape)
        )
        self.max_level = max(block.level for block in self.blocks)
        if any(
            axis_count * 2**self.max_level > np.iinfo(np.int64).max // 4
            for axis_count in self.root_shape
        ):
            raise ValueError("mesh logical coordinates exceed supported int64 range")
        self._ids = {}
        self._leaves = {}
        self._parents = set()
        for block in self.blocks:
            if block.shape != self.block_shape:
                raise ValueError("all blocks must have the same shape")
            if block.block_id in self._ids:
                raise ValueError("duplicate block_id")
            key = (block.level, *block.logical_location)
            if key in self._leaves:
                raise ValueError("duplicate block extent")
            if any(
                logical_index >= axis_count * 2**block.level
                for logical_index, axis_count in zip(
                    block.logical_location, self.root_blocks
                )
            ):
                raise ValueError("block outside domain")
            self._validate_physical_alignment(block)
            self._ids[block.block_id] = block
            self._leaves[key] = block
            for level in range(block.level):
                self._parents.add(
                    (
                        level,
                        *(
                            x // 2 ** (block.level - level)
                            for x in block.logical_location
                        ),
                    )
                )
        if self._parents.intersection(self._leaves):
            raise ValueError("overlapping active parent and child blocks")
        # Alignment and ancestor checks establish disjoint leaves. Their exact
        # integer volumes must now fill the domain; no finest-grid array is needed.
        covered_units = sum(
            8 ** (self.max_level - block.level) for block in self.blocks
        )
        if covered_units != math.prod(self.root_blocks) * 8**self.max_level:
            raise ValueError("gap in active leaf coverage")
        self._validate_neighbor_balance()

    def _validate_physical_alignment(self, block):
        """Match logical endpoints at local width/ULP precision, without changing bounds.

        Four endpoint ULPs permit representable arithmetic roundoff. If that
        uncertainty exceeds 1/16 cell, physical coordinates cannot reliably
        encode this mesh: translation/rescaling is required, not a larger tolerance.
        """
        for axis, ((low, high), (domain_low, domain_high), logical_index) in enumerate(
            zip(block.bounds, self.domain_bounds, block.logical_location)
        ):
            width = (domain_high - domain_low) / (self.root_blocks[axis] * 2**block.level)
            spacing = width / block.shape[axis]
            actual_spacing = (high - low) / block.shape[axis]
            expected = (
                domain_low + logical_index * width,
                domain_low + (logical_index + 1) * width,
            )
            positive_scales = (width, spacing, actual_spacing)
            if (
                not all(math.isfinite(v) and v > 0 for v in positive_scales)
                or not all(math.isfinite(v) for v in expected)
                or expected[1] <= expected[0]
            ):
                raise ValueError(
                    "nonfinite or collapsed AMR geometry; translate/rescale coordinates"
                )
            endpoint_ulp = max(abs(float(np.spacing(v))) for v in (low, high, *expected))
            allowance = max(1e-12 * width, 4 * endpoint_ulp)
            if allowance > spacing / 16:
                raise ValueError(
                    "AMR geometry precision is coarser than 1/16 cell; "
                    "translate/rescale coordinates"
                )
            if any(
                abs(actual - target) > allowance
                for actual, target in zip((low, high), expected)
            ):
                raise ValueError(
                    "block geometry inconsistent with level and logical location"
                )

    def _validate_neighbor_balance(self):
        """Reject level jumps exceeding one at face, edge and corner contacts.

        Aligned dyadic leaves cannot partially overlap without an ancestor
        overlap. Every forbidden contact has a fine side; probing one quarter
        fine-block width across its boundary lands in the touching coarse leaf.
        Zero-offset axes use the fine midpoint. Thus all 26 directions suffice
        even when a coarse-side midpoint would miss a small fine neighbor.
        """
        # Probe every face/edge/corner from each fine block. The fine-side probe
        # detects every excessive jump; no finest-level occupancy grid is built.
        for block in self.blocks:
            for offset in product((-1, 0, 1), repeat=3):
                if offset == (0, 0, 0):
                    continue
                neighbor_point = []
                outside = False
                for axis, direction in enumerate(offset):
                    coordinate = block.logical_location[axis] + (
                        0.5 if direction == 0 else (-0.25 if direction < 0 else 1.25)
                    )
                    axis_count = self.root_blocks[axis] * 2**block.level
                    if coordinate < 0 or coordinate >= axis_count:
                        if self.boundaries[axis][0] == "periodic":
                            coordinate %= axis_count
                        else:
                            outside = True
                            break
                    neighbor_point.append(coordinate)
                if outside:
                    continue
                for level in range(max(0, block.level - 1)):
                    key = (
                        level,
                        *(
                            math.floor(coordinate / 2 ** (block.level - level))
                            for coordinate in neighbor_point
                        ),
                    )
                    if key in self._leaves:
                        raise ValueError("neighbor level difference exceeds one")

    def block(self, block_id: str) -> AMRBlock:
        """Resolve a stable block ID; raise KeyError for an unknown ID."""
        return self._ids[block_id]

    def leaf_at(
        self, level: int, logical_location: tuple[int, int, int]
    ) -> AMRBlock | None:
        """Exact active-leaf lookup in level-local block coordinates.

        Return None for refined parents or absent locations. Boundary mapping
        and searching other levels belong to the caller, not this lookup.
        """
        return self._leaves.get((level, *logical_location))

    @property
    def total_cells(self):
        """Count active cells with Python integers, without ghosts or covered parents."""
        return sum(math.prod(b.shape) for b in self.blocks)


@dataclass(frozen=True)
class AMRExecutionOptions:
    """Tile dimensions and byte limits for AMR workspace planning and primitive caching."""

    tile_shape: tuple[int, int, int] = (64, 64, 64)
    cache_bytes: int = 64 * 1024**2
    memory_budget_bytes: int = 1024**3

    def __post_init__(self):
        """Reject invalid tile dimensions and non-integer or negative byte limits."""
        object.__setattr__(
            self, "tile_shape", _ints(self.tile_shape, "tile_shape", True)
        )
        for name in ("cache_bytes", "memory_budget_bytes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
        if self.memory_budget_bytes == 0:
            raise ValueError("memory budget must be positive")

    def fitted_shape(self, halo=0):
        """Shrink the longest tile axis until the halo-inclusive workspace estimate fits.

        The 2048-byte cell allowance is a planning estimate, not an RSS guarantee.
        Caller-owned snapshot arrays and baseline interpreter memory are separate.
        """
        tile_shape = list(self.tile_shape)
        # Generous per-cell allowance for local gradients, jumps, normals,
        # result arrays and vectorized coordinate work. Cache is bounded too.
        available_bytes = self.memory_budget_bytes - min(
            self.cache_bytes, self.memory_budget_bytes // 4
        )
        while (
            math.prod(axis_cells + 2 * halo for axis_cells in tile_shape) * 2048
            > available_bytes
        ):
            longest_axis = int(np.argmax(tile_shape))
            if tile_shape[longest_axis] == 1:
                raise ValueError(
                    "memory budget cannot hold minimum tile and halo workspace"
                )
            tile_shape[longest_axis] = max(1, tile_shape[longest_axis] // 2)
        return tuple(tile_shape)


@dataclass(frozen=True)
class AMRTile:
    """One output core with half-open, block-local xyz cell bounds; no ghosts."""

    block_id: str
    start: tuple[int, int, int]
    stop: tuple[int, int, int]

    @property
    def slices(self):
        """Convert the core bounds to three NumPy slices."""
        return tuple(slice(a, b) for a, b in zip(self.start, self.stop))


def iter_tiles(mesh, options=None, halo=0):
    """Yield every active core once, keeping enumeration memory independent of tile count."""
    tile_shape = (options or AMRExecutionOptions()).fitted_shape(halo)
    for block in mesh.blocks:
        # Nested range loops avoid itertools.product materializing huge ranges.
        for x in range(0, block.shape[0], tile_shape[0]):
            for y in range(0, block.shape[1], tile_shape[1]):
                for z in range(0, block.shape[2], tile_shape[2]):
                    start = (x, y, z)
                    yield AMRTile(
                        block.block_id,
                        start,
                        tuple(
                            min(axis_start + tile_cells, limit)
                            for axis_start, tile_cells, limit in zip(
                                start, tile_shape, block.shape
                            )
                        ),
                    )


def validate_slices(slices, shape):
    """Require explicit nonempty interior xyz slices with unit steps."""
    if len(slices) != 3 or any(
        not isinstance(s, slice)
        or s.step not in (None, 1)
        or isinstance(s.start, (bool, np.bool_))
        or isinstance(s.stop, (bool, np.bool_))
        or not isinstance(s.start, Integral)
        or not isinstance(s.stop, Integral)
        or not 0 <= s.start < s.stop <= n
        for s, n in zip(slices, shape)
    ):
        raise ValueError("region requires three nonempty in-bounds unit-step slices")
    return tuple(slice(int(s.start), int(s.stop)) for s in slices)
