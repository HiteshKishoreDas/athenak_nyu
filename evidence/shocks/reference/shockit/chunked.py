"""Compatibility wrapper for the explicit `shockit.chunking` package."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .chunking.layout import (
    ChunkShape,
    ChunkSpec,
    balanced_axis_slices,
    iter_balanced_slices,
    load_chunk_layout,
    save_chunked_field,
    validate_payload,
    _validate_chunk_extent,
    integer_triple,
    load_core,
)
from .chunking.writer import join_all_fields, join_chunked_field


def write_chunk_file(
    field_dir: str | Path, spec: ChunkSpec, values: np.ndarray
) -> None:
    """Write one chunk file using an existing core-grid spec."""

    output_dir = Path(field_dir)
    values = np.asarray(values)
    _validate_chunk_extent(output_dir / spec.path.name, spec.start, spec.stop, spec.full_shape)
    dtype = None
    first_path = next(iter(sorted(output_dir.glob("chunk_*.npz"))), None)
    if first_path is not None:
        with np.load(first_path) as archive:
            first_spec = ChunkSpec(first_path, *(integer_triple(archive[key], f"{first_path}: {key}")
                for key in ("start", "stop", "full_shape")))
        dtype = load_core(output_dir, first_spec).dtype
    validate_payload(values, spec.core_shape, output_dir / spec.path.name, dtype)
    output_dir.mkdir(parents=True, exist_ok=True)
    np.savez(
        output_dir / spec.path.name,
        values=values,
        start=np.array(spec.start, dtype=np.int64),
        stop=np.array(spec.stop, dtype=np.int64),
        full_shape=np.array(spec.full_shape, dtype=np.int64),
    )


def load_chunk_specs(field_dir: str | Path) -> list[ChunkSpec]:
    """Compatibility wrapper returning chunk specs from one field directory."""

    return list(load_chunk_layout(field_dir).specs)
