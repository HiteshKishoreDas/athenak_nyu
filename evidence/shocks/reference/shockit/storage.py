"""Storage-mode helpers for array-backed and chunked field inputs."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np

from .chunking.layout import ChunkLayout, ChunkSpec, load_chunk_layout, load_core

ArrayLike = Any
StorageMode = Literal["in_memory", "chunked"]


@dataclass(frozen=True)
class ChunkedFieldReference:
    """Lightweight reference to one chunked field directory."""

    path: Path
    layout: ChunkLayout

    @property
    def shape(self) -> tuple[int, int, int]:
        """Full xyz cell counts from metadata; no field arrays are read."""
        return self.layout.shape

    def load(self) -> np.ndarray:
        """Reconstruct the full field array in memory."""

        dtype = load_core(self.path, self.layout.specs[0]).dtype
        full_array = np.empty(self.layout.shape, dtype=dtype)
        for spec, values in self.iter_chunks():
            chunk_slice = tuple(
                slice(spec.start[axis], spec.stop[axis]) for axis in range(3)
            )
            full_array[chunk_slice] = values
        return full_array

    def iter_chunks(self) -> Iterator[tuple[ChunkSpec, np.ndarray]]:
        """Yield `(ChunkSpec, values)` pairs without reconstructing the full field."""

        dtype = load_core(self.path, self.layout.specs[0]).dtype
        for spec in self.layout.specs:
            yield spec, load_core(self.path, spec, dtype)


FieldSource = ArrayLike | str | Path | ChunkedFieldReference


def detect_storage_mode(field_values: dict[str, FieldSource]) -> StorageMode:
    """Return the shared storage mode for primitive fields."""

    pathlike_count = sum(
        _is_chunked_field_source(value) for value in field_values.values()
    )
    if pathlike_count == 0:
        return "in_memory"
    if pathlike_count == len(field_values):
        return "chunked"
    raise ValueError(
        "FluidCube fields must be either all in-memory arrays or all chunked field paths."
    )


def validate_in_memory_fields(
    field_values: dict[str, FieldSource],
) -> dict[str, np.ndarray]:
    """Convert and validate in-memory primitive arrays."""

    arrays = {
        name: _as_float_array(name, value) for name, value in field_values.items()
    }
    shape = arrays["rho"].shape
    if len(shape) != 3:
        raise ValueError("FluidCube fields must be 3D arrays.")
    if any(size <= 0 for size in shape):
        raise ValueError("FluidCube fields must have positive axes; empty arrays are unsupported.")
    for name, value in arrays.items():
        if value.shape != shape:
            raise ValueError(f"Field {name} has shape {value.shape}, expected {shape}.")
        if not np.all(np.isfinite(value)):
            raise ValueError(f"Field {name} must contain only finite values.")
    if np.any(arrays["rho"] <= 0.0):
        raise ValueError("rho must be strictly positive.")
    if np.any(arrays["pressure"] <= 0.0):
        raise ValueError("pressure must be strictly positive.")
    return arrays


def validate_chunked_fields(
    field_values: dict[str, FieldSource],
) -> dict[str, ChunkedFieldReference]:
    """Validate chunked field references without loading full arrays."""

    refs = {
        name: _as_chunked_field_reference(name, value)
        for name, value in field_values.items()
    }
    base_signature: (
        tuple[
            tuple[int, int, int],
            tuple[tuple[tuple[int, int, int], tuple[int, int, int]], ...],
        ]
        | None
    ) = None
    for name, ref in refs.items():
        signature = _layout_signature(ref.layout)
        if base_signature is None:
            base_signature = signature
            continue
        if signature != base_signature:
            raise ValueError(
                f"Chunked field {name} has a layout incompatible with the other primitive fields."
            )
    return refs


def _as_float_array(name: str, value: ArrayLike) -> np.ndarray:
    """Convert one primitive field to a finite float array with useful errors."""

    if isinstance(value, np.lib.npyio.NpzFile):
        keys = ", ".join(sorted(value.files))
        raise ValueError(
            f"Field {name} received a NumPy .npz archive, not an array. "
            f'Load one dataset from the archive first, for example np.load(...)["arr_0"]. '
            f"Available keys: {keys}."
        )
    try:
        return np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Field {name} could not be converted to a float array."
        ) from exc


def _is_chunked_field_source(value: FieldSource) -> bool:
    """Recognize directory paths and existing chunked references."""
    return isinstance(value, (str, Path, ChunkedFieldReference))


def _as_chunked_field_reference(name: str, value: FieldSource) -> ChunkedFieldReference:
    """Validate a field directory and attach its chunk metadata."""
    if isinstance(value, ChunkedFieldReference):
        return value
    if not isinstance(value, (str, Path)):
        raise ValueError(
            f"Field {name} must be a path-like chunk directory when using chunked FluidCube storage."
        )
    field_dir = Path(value)
    if not field_dir.exists():
        raise FileNotFoundError(
            f"Chunked field directory for {name} does not exist: {field_dir}"
        )
    return ChunkedFieldReference(path=field_dir, layout=load_chunk_layout(field_dir))


def _layout_signature(
    layout: ChunkLayout,
) -> tuple[
    tuple[int, int, int], tuple[tuple[tuple[int, int, int], tuple[int, int, int]], ...]
]:
    """Compare global shape and ordered chunk extents, independent of field paths."""
    return (
        layout.shape,
        tuple((spec.start, spec.stop) for spec in layout.specs),
    )
