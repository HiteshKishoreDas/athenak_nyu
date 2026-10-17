"""Chunk layout primitives for streamed chunked workflows."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math
from numbers import Integral
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ChunkSpec:
    """One chunk file and its core-grid extent."""

    path: Path
    start: tuple[int, int, int]
    stop: tuple[int, int, int]
    full_shape: tuple[int, int, int]
    grid_index: tuple[int, int, int] = (0, 0, 0)

    @property
    def core_shape(self) -> tuple[int, int, int]:
        """Cell counts in this core, excluding halos."""
        return tuple(self.stop[axis] - self.start[axis] for axis in range(3))


ChunkShape = int | tuple[int, int, int]


@dataclass(frozen=True)
class ChunkLayout:
    """Resolved chunk grid shared by one field store."""

    specs: tuple[ChunkSpec, ...]
    axis_slices: tuple[tuple[slice, ...], tuple[slice, ...], tuple[slice, ...]]

    @property
    def shape(self) -> tuple[int, int, int]:
        """Full xyz cell counts shared by all chunk specifications."""
        return self.specs[0].full_shape

    def spec_for_grid_index(self, grid_index: tuple[int, int, int]) -> ChunkSpec:
        """Resolve a chunk's Cartesian position in the chunk lattice."""
        for spec in self.specs:
            if spec.grid_index == grid_index:
                return spec
        raise KeyError(f"Missing chunk at grid index {grid_index}.")

    def iter_face_neighbor_pairs(self, periodic: bool = True):
        """Yield chunk pairs that share a face in the positive axis direction."""

        for spec in self.specs:
            for axis, axis_chunks in enumerate(self.axis_slices):
                if len(axis_chunks) == 1 and not periodic:
                    continue
                neighbor_index = list(spec.grid_index)
                if spec.grid_index[axis] + 1 < len(axis_chunks):
                    neighbor_index[axis] += 1
                elif periodic:
                    neighbor_index[axis] = 0
                else:
                    continue
                yield axis, spec, self.spec_for_grid_index(tuple(neighbor_index))


def balanced_axis_slices(axis_size: int, target_chunk_size: int) -> list[slice]:
    """Split one axis into nearly equal slices with size close to target."""

    if any(isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value <= 0
           for value in (axis_size, target_chunk_size)):
        raise ValueError("Axis and chunk sizes must be positive integers.")
    chunk_count = max(1, math.ceil(axis_size / target_chunk_size))
    boundaries = np.linspace(0, axis_size, num=chunk_count + 1, dtype=int)
    return [
        slice(int(boundaries[index]), int(boundaries[index + 1]))
        for index in range(chunk_count)
    ]


def iter_balanced_slices(shape: tuple[int, int, int], target_chunk_size: ChunkShape):
    """Yield balanced 3D chunk slices for a full cube shape."""

    axis_chunk_sizes = normalize_chunk_shape(target_chunk_size)
    axis_slices = [
        balanced_axis_slices(axis_size, axis_chunk_sizes[axis])
        for axis, axis_size in enumerate(shape)
    ]
    for x_slice, y_slice, z_slice in product(*axis_slices):
        yield x_slice, y_slice, z_slice


def normalize_chunk_shape(chunk_size: ChunkShape) -> tuple[int, int, int]:
    """Return a per-axis chunk shape."""

    if isinstance(chunk_size, Integral) and not isinstance(chunk_size, (bool, np.bool_)):
        chunk_size = (chunk_size,) * 3
    values = integer_triple(chunk_size, "chunk_size")
    if any(value <= 0 for value in values):
        raise ValueError("Chunk sizes must be positive integers.")
    return values


def save_chunked_field(
    field_name: str,
    values: np.ndarray,
    output_dir: str | Path,
    target_chunk_size: ChunkShape,
) -> None:
    """Write one 3D field as a directory of chunked `.npz` files."""

    output_path = Path(output_dir)
    field_dir = output_path / field_name
    values = np.asarray(values)
    validate_payload(values, values.shape, field_dir)
    normalize_chunk_shape(target_chunk_size)
    require_empty_output(field_dir)
    field_dir.mkdir(parents=True, exist_ok=True)

    for chunk_index, chunk_slice in enumerate(
        iter_balanced_slices(values.shape, target_chunk_size)
    ):
        np.savez(
            field_dir / _chunk_filename(chunk_index, chunk_slice),
            values=values[chunk_slice],
            start=np.array(
                [chunk_slice[0].start, chunk_slice[1].start, chunk_slice[2].start],
                dtype=np.int64,
            ),
            stop=np.array(
                [chunk_slice[0].stop, chunk_slice[1].stop, chunk_slice[2].stop],
                dtype=np.int64,
            ),
            full_shape=np.array(values.shape, dtype=np.int64),
        )


def load_chunk_layout(field_dir: str | Path) -> ChunkLayout:
    """Read chunk metadata from one field directory."""

    field_path = Path(field_dir)
    chunk_files = sorted(field_path.glob("chunk_*.npz"))
    if not chunk_files:
        raise FileNotFoundError(f"No chunk files found in {field_path}.")

    raw_specs: list[
        tuple[Path, tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]
    ] = []
    for chunk_path in chunk_files:
        with np.load(chunk_path) as chunk:
            raw_specs.append(
                (
                    chunk_path,
                    integer_triple(chunk["start"], f"{chunk_path}: start"),
                    integer_triple(chunk["stop"], f"{chunk_path}: stop"),
                    integer_triple(chunk["full_shape"], f"{chunk_path}: full_shape"),
                )
            )

    full_shapes = {full_shape for _, _, _, full_shape in raw_specs}
    if len(full_shapes) != 1:
        raise ValueError(f"Inconsistent full shapes in {field_path}.")
    full_shape = next(iter(full_shapes))
    for chunk_path, start, stop, _ in raw_specs:
        _validate_chunk_extent(chunk_path, start, stop, full_shape)

    axis_bounds = tuple(
        tuple(
            slice(start, stop)
            for start, stop in sorted(
                {(spec[1][axis], spec[2][axis]) for spec in raw_specs},
                key=lambda pair: pair[0],
            )
        )
        for axis in range(3)
    )
    for axis, axis_slices in enumerate(axis_bounds):
        _validate_axis_slices(
            axis_slices, full_shape[axis], field_path=field_path, axis=axis
        )
    axis_lookup = [
        {
            (axis_slice.start, axis_slice.stop): index
            for index, axis_slice in enumerate(axis_slices)
        }
        for axis_slices in axis_bounds
    ]
    specs = [
        ChunkSpec(
            path=path,
            start=start,
            stop=stop,
            full_shape=full_shape,
            grid_index=tuple(
                axis_lookup[axis][(start[axis], stop[axis])] for axis in range(3)
            ),
        )
        for path, start, stop, full_shape in raw_specs
    ]
    specs.sort(key=lambda spec: spec.grid_index)

    if len({spec.grid_index for spec in specs}) != len(specs):
        raise ValueError(f"Chunk coverage in {field_path} has duplicate grid indices.")
    expected_chunk_count = math.prod(len(axis_slices) for axis_slices in axis_bounds)
    if len(specs) != expected_chunk_count:
        raise ValueError(
            f"Chunk coverage in {field_path} does not form a complete regular grid."
        )

    return ChunkLayout(specs=tuple(specs), axis_slices=axis_bounds)


def _validate_chunk_extent(
    chunk_path: Path,
    start: tuple[int, int, int],
    stop: tuple[int, int, int],
    full_shape: tuple[int, int, int],
) -> None:
    """Reject inconsistent dimensions or chunk bounds outside the global array."""
    start = integer_triple(start, f"{chunk_path}: start")
    stop = integer_triple(stop, f"{chunk_path}: stop")
    full_shape = integer_triple(full_shape, f"{chunk_path}: full_shape")
    if any(size <= 0 for size in full_shape):
        raise ValueError(f"{chunk_path}: full shape must have positive axes.")
    for axis in range(3):
        if stop[axis] <= start[axis]:
            raise ValueError(
                f"Invalid chunk extent in {chunk_path}: axis {axis} has start={start[axis]} stop={stop[axis]}."
            )
        if start[axis] < 0 or stop[axis] > full_shape[axis]:
            raise ValueError(
                f"Invalid chunk extent in {chunk_path}: axis {axis} must stay within [0, {full_shape[axis]}]."
            )


def _validate_axis_slices(
    axis_slices: tuple[slice, ...],
    axis_size: int,
    *,
    field_path: Path,
    axis: int,
) -> None:
    """Require contiguous, non-overlapping axis coverage and consistent edge chunks."""
    if not axis_slices:
        raise ValueError(
            f"Chunk coverage in {field_path} does not form a complete regular grid."
        )
    for axis_slice in axis_slices:
        if axis_slice.stop <= axis_slice.start:
            raise ValueError(
                f"Chunk coverage in {field_path} has invalid chunk extent on axis {axis}: "
                f"{axis_slice.start}:{axis_slice.stop}."
            )
    if axis_slices[0].start != 0:
        raise ValueError(
            f"Chunk coverage in {field_path} has a gap before the first chunk on axis {axis}: "
            f"starts at {axis_slices[0].start} instead of 0."
        )
    if axis_slices[-1].stop != axis_size:
        raise ValueError(
            f"Chunk coverage in {field_path} does not reach the full domain on axis {axis}: "
            f"ends at {axis_slices[-1].stop} instead of {axis_size}."
        )
    for previous, current in zip(axis_slices, axis_slices[1:], strict=False):
        if previous.stop < current.start:
            raise ValueError(
                f"Chunk coverage in {field_path} has a gap on axis {axis}: "
                f"{previous.stop} -> {current.start}."
            )
        if previous.stop > current.start:
            raise ValueError(
                f"Chunk coverage in {field_path} has an overlap on axis {axis}: "
                f"{previous.stop} overlaps {current.start}."
            )


def _chunk_filename(chunk_index: int, chunk_slice: tuple[slice, slice, slice]) -> str:
    """Build the conventional filename from xyz chunk-grid indices."""
    return (
        f"chunk_{chunk_index:04d}"
        f"_x{chunk_slice[0].start:04d}-{chunk_slice[0].stop:04d}"
        f"_y{chunk_slice[1].start:04d}-{chunk_slice[1].stop:04d}"
        f"_z{chunk_slice[2].start:04d}-{chunk_slice[2].stop:04d}.npz"
    )


def integer_triple(values, context: str) -> tuple[int, int, int]:
    """Reject coercion of fractional/bool/oversized indices before integer conversion."""
    array = np.asarray(values)
    if array.shape != (3,) or array.dtype.kind not in "iu":
        raise ValueError(f"{context}: expected three non-boolean integers.")
    if any(int(value) > np.iinfo(np.int64).max for value in array):
        raise ValueError(f"{context}: integer exceeds supported int64 range.")
    return tuple(int(value) for value in array)


def validate_payload(values: np.ndarray, shape, context, dtype=None) -> None:
    """Require exact xyz core shape and real numeric/bool storage, without coercion."""
    if values.ndim != 3 or any(size <= 0 for size in values.shape) or values.shape != tuple(shape):
        raise ValueError(f"{context}: expected shape {tuple(shape)}, actual shape {values.shape}; axes must be positive.")
    if values.dtype.kind not in "biuf" or (dtype is not None and values.dtype != dtype):
        raise ValueError(f"{context}: expected dtype {dtype or 'real numeric or bool'}, actual dtype {values.dtype}.")


def load_core(field_dir: Path, spec: ChunkSpec, dtype=None) -> np.ndarray:
    """Read only this archive, verify requested geometry and dtype, then close it.

    Callers infer a field dtype from its first core. Every later read checks that
    dtype; untouched cores remain lazy. No input array is modified or coerced.
    """
    path = field_dir / spec.path.name
    _validate_chunk_extent(path, spec.start, spec.stop, spec.full_shape)
    with np.load(path) as archive:
        for name, expected in (("start", spec.start), ("stop", spec.stop), ("full_shape", spec.full_shape)):
            actual = integer_triple(archive[name], f"{path}: {name}")
            if actual != expected:
                raise ValueError(f"{path}: expected {name} {expected}, actual {actual}.")
        try:
            values = np.asarray(archive["values"])
        except ValueError as exc:
            raise ValueError(f"{path}: unsupported payload dtype.") from exc
        validate_payload(values, spec.core_shape, path, dtype)
        return values


def require_empty_output(path: Path) -> None:
    """Preflight a new export/run; low-level writes deliberately allow ongoing work."""
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"Output target must be empty; nonempty target: {path}")
