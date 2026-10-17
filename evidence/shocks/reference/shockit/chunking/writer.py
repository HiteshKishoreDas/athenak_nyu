"""Chunked field writers and reconstruction helpers."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Protocol

import numpy as np

from ..provenance import atomic_write_json
from .layout import ChunkLayout, ChunkSpec, load_chunk_layout, load_core, validate_payload, _validate_chunk_extent
from .reader import ChunkedFieldReader, NpzFieldReader


class ChunkedFieldWriter(Protocol):
    """Abstract chunked field writer."""

    def write_core(self, field_name: str, spec: ChunkSpec, values: np.ndarray) -> None:
        """Persist one named field core and its chunk metadata."""
        ...


class ChunkedOutputStore(Protocol):
    """Abstract chunked output store."""

    @property
    def root(self) -> Path:
        """Directory containing the persisted output fields."""
        ...

    @property
    def layout(self) -> ChunkLayout:
        """Validated chunk geometry shared by the stored fields."""
        ...

    def field_reader(self, field_name: str) -> ChunkedFieldReader:
        """Return a reader for the named stored field."""
        ...

    def write_core(self, field_name: str, spec: ChunkSpec, values: np.ndarray) -> None:
        """Persist one named field core and its chunk metadata."""
        ...

    def write_summary(self, summary: dict[str, Any]) -> None:
        """Persist the result summary."""
        ...

    def read_summary(self) -> dict[str, Any]:
        """Load the persisted result summary."""
        ...


class NpzChunkedOutput:
    """Chunked `.npz` output store."""

    def __init__(self, root: str | Path, layout: ChunkLayout) -> None:
        """Attach the output root and chunk layout; create directories when writing."""
        self.root = Path(root)
        self.layout = layout
        self._dtypes: dict[str, np.dtype] = {}

    def write_core(self, field_name: str, spec: ChunkSpec, values: np.ndarray) -> None:
        """Write one field core together with its chunk geometry."""
        expected = self.layout.spec_for_grid_index(spec.grid_index)
        if (spec.start, spec.stop, spec.full_shape, spec.path.name) != (
            expected.start, expected.stop, expected.full_shape, expected.path.name
        ):
            raise ValueError(f"Field {field_name}: core specification does not match output layout.")
        field_dir = self.root / field_name
        values = np.asarray(values)
        _validate_chunk_extent(field_dir / spec.path.name, spec.start, spec.stop, spec.full_shape)
        expected_dtype = self._dtypes.get(field_name)
        first_path = field_dir / self.layout.specs[0].path.name
        if expected_dtype is None and first_path.exists():
            expected_dtype = load_core(field_dir, self.layout.specs[0]).dtype
        validate_payload(values, spec.core_shape, field_dir / spec.path.name, expected_dtype)
        self._dtypes[field_name] = values.dtype
        field_dir.mkdir(parents=True, exist_ok=True)
        np.savez(
            field_dir / spec.path.name,
            values=values,
            start=np.array(spec.start, dtype=np.int64),
            stop=np.array(spec.stop, dtype=np.int64),
            full_shape=np.array(spec.full_shape, dtype=np.int64),
        )

    def field_reader(self, field_name: str) -> NpzFieldReader:
        """Expose a stored field through the chunk reader interface."""
        field_dir = self.root / field_name
        dtype = load_core(field_dir, self.layout.specs[0]).dtype
        return NpzFieldReader(field_dir=field_dir, layout=self.layout, dtype=dtype)

    def write_summary(self, summary: dict[str, Any]) -> None:
        """Serialize the final summary in the output directory."""
        self.root.mkdir(parents=True, exist_ok=True)
        atomic_write_json(self.root / "summary.json", summary)

    def read_summary(self) -> dict[str, Any]:
        """Read the previously persisted summary JSON."""
        if (self.root / "incomplete.json").exists():
            raise ValueError(f"Output is incomplete: {self.root}; retry standalone reduction.")
        summary_path = self.root / "summary.json"
        if not summary_path.exists():
            return {}
        return json.loads(summary_path.read_text())


def join_chunked_field(field_dir: str | Path) -> np.ndarray:
    """Reconstruct a full array from chunked field files."""

    layout = load_chunk_layout(field_dir)
    field_dir = Path(field_dir)
    dtype = load_core(field_dir, layout.specs[0]).dtype
    full_array = np.empty(layout.shape, dtype=dtype)

    for spec in layout.specs:
        chunk_slice = tuple(slice(spec.start[axis], spec.stop[axis]) for axis in range(3))
        full_array[chunk_slice] = load_core(field_dir, spec, dtype)
    return full_array


def join_all_fields(root_dir: str | Path) -> dict[str, np.ndarray]:
    """Reconstruct all chunked fields in one root directory."""

    root_path = Path(root_dir)
    field_arrays: dict[str, np.ndarray] = {}
    for field_dir in sorted(
        path
        for path in root_path.iterdir()
        if path.is_dir() and not path.name.startswith(".")
    ):
        if not list(field_dir.glob("chunk_*.npz")):
            continue
        field_arrays[field_dir.name] = join_chunked_field(field_dir)
    return field_arrays
