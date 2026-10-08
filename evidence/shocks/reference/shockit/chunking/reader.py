"""Chunked field readers."""

from __future__ import annotations

from dataclasses import dataclass
from copy import deepcopy
from itertools import product
import os
from pathlib import Path
from typing import TYPE_CHECKING, Protocol

import numpy as np

from .layout import ChunkLayout, ChunkSpec, load_chunk_layout, load_core

if TYPE_CHECKING:
    from ..fields import FluidCube

DEFAULT_FIELD_MAPPING = {
    "rho": "rho",
    "pressure": "prs",
    "vx": "v1",
    "vy": "v2",
    "vz": "v3",
}


class ChunkedFieldReader(Protocol):
    """Abstract chunked field reader."""

    @property
    def shape(self) -> tuple[int, int, int]:
        """Full xyz cell counts from the stored layout."""
        ...

    @property
    def layout(self) -> ChunkLayout:
        """Validated chunk geometry shared by the stored fields."""
        ...

    def read_core(self, spec: ChunkSpec) -> np.ndarray:
        """Read the specified chunk core without halos."""
        ...

    def read_with_halo(self, spec: ChunkSpec, halo: int) -> np.ndarray:
        """Read a core with the requested periodic ghost width."""
        ...


class ChunkedInputStore(Protocol):
    """Abstract chunked input store."""

    @property
    def shape(self) -> tuple[int, int, int]:
        """Full xyz cell counts from the stored layout."""
        ...

    @property
    def layout(self) -> ChunkLayout:
        """Validated chunk geometry shared by the stored fields."""
        ...

    def field_reader(self, field_name: str) -> ChunkedFieldReader:
        """Return a reader for the named stored field."""
        ...


@dataclass(frozen=True)
class NpzFieldReader:
    """Read one field from a chunked `.npz` directory."""

    field_dir: Path
    layout: ChunkLayout
    dtype: np.dtype

    @property
    def shape(self) -> tuple[int, int, int]:
        """Full xyz cell counts from the stored layout."""
        return self.layout.shape

    def read_core(self, spec: ChunkSpec) -> np.ndarray:
        """Read one stored core and close its NPZ archive."""
        return load_core(self.field_dir, spec, self.dtype)

    def read_with_halo(self, spec: ChunkSpec, halo: int) -> np.ndarray:
        """Fill a core and halos with global periodic wrapping across stored chunks."""
        request_slices = tuple(
            slice(spec.start[axis] - halo, spec.stop[axis] + halo) for axis in range(3)
        )
        output = np.empty(
            tuple(request.stop - request.start for request in request_slices),
            dtype=self.dtype,
        )
        axis_segments = [
            _periodic_axis_segments(self.shape[axis], request.start, request.stop)
            for axis, request in enumerate(request_slices)
        ]
        for x_segment, y_segment, z_segment in product(*axis_segments):
            source_slices = (
                slice(x_segment[0], x_segment[1]),
                slice(y_segment[0], y_segment[1]),
                slice(z_segment[0], z_segment[1]),
            )
            dest_slices = (
                slice(x_segment[2], x_segment[3]),
                slice(y_segment[2], y_segment[3]),
                slice(z_segment[2], z_segment[3]),
            )
            _fill_region_from_specs(
                output=output,
                dest_slices=dest_slices,
                source_slices=source_slices,
                field_dir=self.field_dir,
                specs=self.layout.specs,
                dtype=self.dtype,
            )
        return output


class NpzChunkedInput:
    """Chunked `.npz` input store."""

    def __init__(
        self,
        root: str | Path | None = None,
        field_mapping: dict[str, str] | None = None,
        *,
        field_paths: dict[str, str | Path] | None = None,
    ) -> None:
        """Resolve field directories and verify matching chunk layouts."""
        if (root is None) == (field_paths is None):
            raise ValueError(
                "Provide either root=... or field_paths=..., but not both."
            )
        self.root = Path(root) if root is not None else _common_root(field_paths or {})
        self.field_mapping = dict(
            DEFAULT_FIELD_MAPPING if field_mapping is None else field_mapping
        )
        self.field_paths = _resolve_field_paths(
            self.root, self.field_mapping, field_paths
        )
        self.metadata = {}
        self._layout: ChunkLayout | None = None
        self._readers: dict[str, NpzFieldReader] = {}
        self._build_readers()

    @classmethod
    def from_fluid_cube(cls, cube: FluidCube) -> NpzChunkedInput:
        """Build a chunked input store from a chunked `FluidCube`."""

        if not cube.is_chunked:
            raise ValueError(
                "NpzChunkedInput.from_fluid_cube() requires a chunked FluidCube."
            )
        store = cls(
            field_paths={
                "rho": cube.rho.path,
                "pressure": cube.pressure.path,
                "vx": cube.vx.path,
                "vy": cube.vy.path,
                "vz": cube.vz.path,
            }
        )
        store.metadata = deepcopy(cube.metadata)
        return store

    @property
    def shape(self) -> tuple[int, int, int]:
        """Full xyz cell counts from the stored layout."""
        return self.layout.shape

    @property
    def layout(self) -> ChunkLayout:
        """Validated chunk geometry shared by the stored fields."""
        assert self._layout is not None
        return self._layout

    def field_reader(self, field_name: str) -> NpzFieldReader:
        """Return a reader for the named stored field."""
        return self._readers[field_name]

    def _build_readers(self) -> None:
        """Validate shared layouts and inspect each first chunk for its dtype."""
        base_signature: (
            tuple[
                tuple[int, int, int],
                tuple[tuple[tuple[int, int, int], tuple[int, int, int]], ...],
            ]
            | None
        ) = None
        for logical_name, field_dir in self.field_paths.items():
            layout = load_chunk_layout(field_dir)
            signature = _layout_signature(layout)
            if base_signature is None:
                base_signature = signature
                self._layout = layout
            elif signature != base_signature:
                raise ValueError(f"Chunk layout mismatch in {field_dir}.")
            dtype = load_core(field_dir, layout.specs[0]).dtype
            self._readers[logical_name] = NpzFieldReader(
                field_dir=field_dir, layout=layout, dtype=dtype
            )


def _resolve_field_paths(
    root: Path,
    field_mapping: dict[str, str],
    explicit_field_paths: dict[str, str | Path] | None,
) -> dict[str, Path]:
    """Use explicit directories or resolve logical field aliases under the root."""
    if explicit_field_paths is not None:
        return {name: Path(path) for name, path in explicit_field_paths.items()}
    return {
        logical_name: root / directory_name
        for logical_name, directory_name in field_mapping.items()
    }


def _common_root(field_paths: dict[str, str | Path]) -> Path:
    """Find the common directory prefix of explicit field paths."""
    resolved = [str(Path(path)) for path in field_paths.values()]
    return Path(os.path.commonpath(resolved)) if resolved else Path(".")


def _layout_signature(
    layout: ChunkLayout,
) -> tuple[
    tuple[int, int, int], tuple[tuple[tuple[int, int, int], tuple[int, int, int]], ...]
]:
    """Compare global shape and ordered core bounds across fields."""
    return (
        layout.shape,
        tuple((spec.start, spec.stop) for spec in layout.specs),
    )


def _periodic_axis_segments(
    axis_size: int, start: int, stop: int
) -> list[tuple[int, int, int, int]]:
    """Split an axis request into source and destination half-open bounds.

    Each tuple is (source_start, source_stop, destination_start, destination_stop).
    Multiple segments support requests wider than the periodic domain.
    """
    total = stop - start
    segments: list[tuple[int, int, int, int]] = []
    local_start = 0
    global_start = start
    while local_start < total:
        wrapped_start = global_start % axis_size
        available = axis_size - wrapped_start
        segment_length = min(available, total - local_start)
        segments.append(
            (
                wrapped_start,
                wrapped_start + segment_length,
                local_start,
                local_start + segment_length,
            )
        )
        local_start += segment_length
        global_start += segment_length
    return segments


def _fill_region_from_specs(
    output: np.ndarray,
    dest_slices: tuple[slice, slice, slice],
    source_slices: tuple[slice, slice, slice],
    field_dir: Path,
    specs: tuple[ChunkSpec, ...],
    dtype: np.dtype,
) -> None:
    """Copy intersecting global source slices into the requested output offsets."""
    for spec in specs:
        overlap = _intersection_slices(source_slices, spec)
        if overlap is None:
            continue
        source_overlap, chunk_overlap = overlap
        dest_overlap = tuple(
            slice(
                dest_slices[axis].start
                + (source_overlap[axis].start - source_slices[axis].start),
                dest_slices[axis].start
                + (source_overlap[axis].stop - source_slices[axis].start),
            )
            for axis in range(3)
        )
        output[dest_overlap] = load_core(field_dir, spec, dtype)[chunk_overlap]


def _intersection_slices(
    source_slices: tuple[slice, slice, slice],
    spec: ChunkSpec,
) -> tuple[tuple[slice, slice, slice], tuple[slice, slice, slice]] | None:
    """Return matching global and chunk-local slices, or None for disjoint regions."""
    source_overlap: list[slice] = []
    chunk_overlap: list[slice] = []
    for axis in range(3):
        start = max(source_slices[axis].start, spec.start[axis])
        stop = min(source_slices[axis].stop, spec.stop[axis])
        if start >= stop:
            return None
        source_overlap.append(slice(start, stop))
        chunk_overlap.append(slice(start - spec.start[axis], stop - spec.start[axis]))
    return tuple(source_overlap), tuple(chunk_overlap)
