"""Chunk-aware connected-component reduction for chunked outputs."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import shutil
import tempfile
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
from scipy import ndimage

from ..local_pass import build_reporter
from ..provenance import update_postprocessing
from ..config import ShockFinderConfig
from ..derived import MachFields
from ..masks import center_score
from .layout import ChunkLayout, ChunkSpec
from .reader import ChunkedInputStore
from .writer import ChunkedOutputStore, NpzChunkedOutput

ProgressCallback = Callable[[str], None]
ProgressReporter = ProgressCallback | bool | None


@dataclass(frozen=True)
class _LabelRef:
    """Identify a component by both its chunk grid index and local label."""

    grid_index: tuple[int, int, int]
    local_label: int


@dataclass(frozen=True)
class _LocalLabelRecord:
    """Record one local component's best cell and its global tie-breaking index."""

    ref: _LabelRef
    spec: ChunkSpec
    best_local_index: tuple[int, int, int]
    best_global_flat_index: int
    best_score: float
    best_mach_pressure: float
    best_mach_temperature: float


class _UnionFind:
    """Merge chunk-local component labels using parent links and path compression."""

    def __init__(self) -> None:
        """Start with no component labels."""
        self.parent: dict[_LabelRef, _LabelRef] = {}

    def add(self, item: _LabelRef) -> None:
        """Register a component as its own root if it is new."""
        self.parent.setdefault(item, item)

    def find(self, item: _LabelRef) -> _LabelRef:
        """Return the representative label and shorten traversed parent chains."""
        parent = self.parent[item]
        if parent != item:
            parent = self.find(parent)
            self.parent[item] = parent
        return parent

    def union(self, left: _LabelRef, right: _LabelRef) -> None:
        """Join the two component roots, preserving the left root as representative."""
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left != root_right:
            self.parent[root_right] = root_left


def finalize_chunked_shock_outputs(
    output_store: ChunkedOutputStore,
    config: ShockFinderConfig,
    base_summary: dict[str, Any],
    *,
    progress: ProgressReporter = True,
) -> dict[str, Any]:
    """Finalize connected-component accounting and optional center reduction."""

    config.validate()
    config = deepcopy(config)
    report = build_reporter(progress)
    summary = deepcopy(base_summary)
    summary["reduce_to_centers"] = config.reduce_to_centers
    summary["center_score"] = config.center_score
    summary["provenance"] = update_postprocessing(summary.get("provenance"), config)
    summary["center_topology"] = "nonperiodic six-connected components"
    summary["component_count_topology"] = "periodic six-connected components when computed"
    memory_bound = 28 * int(np.prod(output_store.layout.shape))
    summary["center_memory_lower_bound_bytes"] = memory_bound
    summary["provenance"]["center_memory"] = {
        "simple_score_lower_bound_bytes": memory_bound,
        "bytes_per_global_cell": 28, "peak_rss": False,
        "exact_medians": "retain detected Mach values",
    }
    if not config.reduce_to_centers:
        reader = output_store.field_reader("full_shock_mask")
        for spec in output_store.layout.specs:
            output_store.write_core("shock_mask", spec, reader.read_core(spec).astype(bool))
        summary.update(
            shock_cells=summary["full_shock_cells"],
            shock_fraction=summary["full_shock_fraction"],
            n_connected_components=None,
            mask_semantics="shock_mask matches full_shock_mask; center reduction disabled",
        )
        for estimator in ("pressure", "temperature"):
            for statistic in ("min", "median", "max"):
                key = f"mach_{estimator}_{statistic}"
                summary[key] = summary[f"full_{key}"]
        output_store.write_summary(summary)
        return summary
    report("[2/2] Reducing connected components across chunk boundaries")
    if config.center_score != "combined":
        report(f"Legacy center workspace lower bound: {memory_bound} bytes (28 bytes/global cell), before temporaries; not peak RSS")

    fields = output_store
    labels_root = Path(tempfile.mkdtemp(prefix=".local_labels-", dir=output_store.root))
    label_store = NpzChunkedOutput(labels_root, output_store.layout)
    try:
        label_records, base_union_find = _label_local_components(
            fields, label_store, config.center_score
        )
        # Preserve existing semantics: center selection matches nonperiodic uniform
        # labeling, while chunked summary counts also join periodic seam components.
        center_union_find = _copy_union_find(base_union_find)
        summary_union_find = _copy_union_find(base_union_find)
        _merge_boundary_components(
            output_store.layout, center_union_find, label_store, periodic=False
        )
        n_connected_components = _merge_boundary_components(
            output_store.layout,
            summary_union_find,
            label_store,
            periodic=True,
        )

        summary["n_connected_components"] = n_connected_components

        centers = _select_global_centers(
            layout=output_store.layout,
            fields=fields,
            label_store=label_store,
            label_records=label_records,
            union_find=center_union_find,
            center_score_mode=config.center_score,
        )
        _write_center_mask(output_store, centers)
        mach_pressure_values = np.array([center[1] for center in centers], dtype=float)
        mach_temperature_values = np.array([center[2] for center in centers], dtype=float)
        summary.update(
            {
                "shock_cells": len(centers),
                "shock_fraction": len(centers) / summary["total_cells"],
                "mask_semantics": "shock_mask is center-reduced; full_shock_mask is unreduced",
                "mach_pressure_min": _finite_stat(mach_pressure_values, np.min),
                "mach_pressure_median": _finite_stat(mach_pressure_values, np.median),
                "mach_pressure_max": _finite_stat(mach_pressure_values, np.max),
                "mach_temperature_min": _finite_stat(mach_temperature_values, np.min),
                "mach_temperature_median": _finite_stat(mach_temperature_values, np.median),
                "mach_temperature_max": _finite_stat(mach_temperature_values, np.max),
            }
        )
        output_store.write_summary(summary)
        return summary
    finally:
        shutil.rmtree(labels_root, ignore_errors=True)


def _label_local_components(
    fields: ChunkedInputStore,
    label_store: ChunkedOutputStore,
    center_score_mode: str,
) -> tuple[dict[_LabelRef, _LocalLabelRecord], _UnionFind]:
    """Label six-connected core cells and save local center candidates."""
    records: dict[_LabelRef, _LocalLabelRecord] = {}
    union_find = _UnionFind()
    mask_reader = fields.field_reader("full_shock_mask")
    div_v_reader = fields.field_reader("div_v")
    pressure_reader = fields.field_reader("mach_pressure")
    temperature_reader = fields.field_reader("mach_temperature")

    for spec in fields.layout.specs:
        mask = mask_reader.read_core(spec).astype(bool)
        labels, count = ndimage.label(mask)
        label_store.write_core("labels", spec, labels.astype(np.int32))
        if count == 0:
            continue
        div_v = div_v_reader.read_core(spec)
        mach_fields = MachFields(
            mach_pressure=pressure_reader.read_core(spec),
            mach_temperature=temperature_reader.read_core(spec),
        )
        simple_score = _simple_score_field(div_v, mach_fields, center_score_mode)
        for local_label in range(1, count + 1):
            if simple_score is None:
                component = labels == local_label
                score = center_score(
                    component, div_v, mach_fields, mode=center_score_mode
                )
                best_flat_index = int(np.argmax(score))
                best_local_index = np.unravel_index(best_flat_index, score.shape)
                best_score = float(score[best_local_index])
            else:
                best_local_index = tuple(
                    int(value)
                    for value in ndimage.maximum_position(
                        simple_score, labels=labels, index=[local_label]
                    )[0]
                )
                best_score = float(simple_score[best_local_index])
            best_global_index = tuple(
                spec.start[axis] + best_local_index[axis] for axis in range(3)
            )
            ref = _LabelRef(spec.grid_index, local_label)
            union_find.add(ref)
            records[ref] = _LocalLabelRecord(
                ref=ref,
                spec=spec,
                best_local_index=best_local_index,
                best_global_flat_index=int(
                    np.ravel_multi_index(best_global_index, spec.full_shape)
                ),
                best_score=best_score,
                best_mach_pressure=float(mach_fields.mach_pressure[best_local_index]),
                best_mach_temperature=float(
                    mach_fields.mach_temperature[best_local_index]
                ),
            )
    return records, union_find


def _merge_boundary_components(
    layout: ChunkLayout,
    union_find: _UnionFind,
    label_store: ChunkedOutputStore,
    *,
    periodic: bool,
) -> int:
    """Join matching face labels, optionally across periodic seams; return group count."""
    label_reader = label_store.field_reader("labels")
    for axis, spec_a, spec_b in layout.iter_face_neighbor_pairs(periodic=periodic):
        labels_a = label_reader.read_core(spec_a)
        labels_b = label_reader.read_core(spec_b)
        face_a = np.take(labels_a, indices=-1, axis=axis)
        face_b = np.take(labels_b, indices=0, axis=axis)
        overlaps = (face_a > 0) & (face_b > 0)
        if not np.any(overlaps):
            continue
        pairs = np.unique(
            np.stack([face_a[overlaps], face_b[overlaps]], axis=1), axis=0
        )
        for label_a, label_b in pairs:
            union_find.union(
                _LabelRef(spec_a.grid_index, int(label_a)),
                _LabelRef(spec_b.grid_index, int(label_b)),
            )
    return len({union_find.find(item) for item in union_find.parent})


def _copy_union_find(union_find: _UnionFind) -> _UnionFind:
    """Copy parent links so center and summary connectivity can diverge safely."""
    copied = _UnionFind()
    copied.parent = dict(union_find.parent)
    return copied


def _select_global_centers(
    layout: ChunkLayout,
    fields: ChunkedInputStore,
    label_store: ChunkedOutputStore,
    label_records: dict[_LabelRef, _LocalLabelRecord],
    union_find: _UnionFind,
    center_score_mode: str,
) -> list[tuple[tuple[int, int, int], float, float]]:
    """Select one global center per component; normalize combined scores across chunks."""
    groups: dict[_LabelRef, list[_LabelRef]] = {}
    for ref in label_records:
        groups.setdefault(union_find.find(ref), []).append(ref)

    if center_score_mode != "combined":
        return _select_global_centers_simple_mode(
            layout=layout,
            fields=fields,
            label_store=label_store,
            label_records=label_records,
            groups=groups,
            union_find=union_find,
            center_score_mode=center_score_mode,
        )

    label_reader = label_store.field_reader("labels")
    div_v_reader = fields.field_reader("div_v")
    pressure_reader = fields.field_reader("mach_pressure")
    temperature_reader = fields.field_reader("mach_temperature")
    centers = []
    for refs in groups.values():
        compression_min, compression_max, mach_min, mach_max = _combined_score_ranges(
            refs,
            label_records,
            label_reader,
            div_v_reader,
            pressure_reader,
            temperature_reader,
        )
        best_candidate: tuple[float, int, tuple[int, int, int], float, float] | None = (
            None
        )
        for ref in refs:
            record = label_records[ref]
            labels = label_reader.read_core(record.spec)
            component = labels == ref.local_label
            div_v = div_v_reader.read_core(record.spec)
            mach_pressure = pressure_reader.read_core(record.spec)
            mach_temperature = temperature_reader.read_core(record.spec)
            mach = np.maximum(
                np.nan_to_num(mach_pressure, nan=-np.inf),
                np.nan_to_num(mach_temperature, nan=-np.inf),
            )
            score = _combined_component_score(
                component,
                -div_v,
                mach,
                compression_min,
                compression_max,
                mach_min,
                mach_max,
            )
            for local_index in zip(*np.where(component), strict=False):
                local_score = float(score[local_index])
                global_index = tuple(
                    record.spec.start[axis] + local_index[axis] for axis in range(3)
                )
                global_flat = int(np.ravel_multi_index(global_index, layout.shape))
                # Highest score wins; ties use the first global C-order cell.
                candidate = (
                    -local_score,
                    global_flat,
                    global_index,
                    float(mach_pressure[local_index]),
                    float(mach_temperature[local_index]),
                )
                if best_candidate is None or candidate[:2] < best_candidate[:2]:
                    best_candidate = candidate
        assert best_candidate is not None
        centers.append((best_candidate[2], best_candidate[3], best_candidate[4]))
    return centers


def build_base_summary_from_output_store(
    output_store: ChunkedOutputStore,
    config: ShockFinderConfig,
    *,
    gamma: float,
) -> dict[str, Any]:
    """Rebuild base chunked summary stats from stored output fields."""

    shock_zone_reader = output_store.field_reader("shock_zone_mask")
    full_shock_reader = output_store.field_reader("full_shock_mask")
    pressure_reader = output_store.field_reader("mach_pressure")
    temperature_reader = output_store.field_reader("mach_temperature")

    shock_zone_cells = 0
    full_shock_cells = 0
    mach_pressure_values: list[np.ndarray] = []
    mach_temperature_values: list[np.ndarray] = []
    for spec in output_store.layout.specs:
        shock_zone = shock_zone_reader.read_core(spec).astype(bool)
        full_shock = full_shock_reader.read_core(spec).astype(bool)
        mach_pressure = pressure_reader.read_core(spec)
        mach_temperature = temperature_reader.read_core(spec)
        shock_zone_cells += int(np.count_nonzero(shock_zone))
        full_shock_cells += int(np.count_nonzero(full_shock))
        mach_pressure_values.append(mach_pressure[full_shock])
        mach_temperature_values.append(mach_temperature[full_shock])

    total_cells = int(np.prod(output_store.layout.shape))
    mach_pressure = _concat(mach_pressure_values)
    mach_temperature = _concat(mach_temperature_values)
    return {
        "grid_shape": output_store.layout.shape,
        "total_cells": total_cells,
        "shock_zone_cells": shock_zone_cells,
        "shock_cells": full_shock_cells,
        "shock_fraction": full_shock_cells / total_cells,
        "full_shock_cells": full_shock_cells,
        "full_shock_fraction": full_shock_cells / total_cells,
        "n_connected_components": None,
        "mask_semantics": "shock_mask matches full_shock_mask until center reduction runs",
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


def _combined_score_ranges(
    refs: list[_LabelRef],
    label_records: dict[_LabelRef, _LocalLabelRecord],
    label_reader,
    div_v_reader,
    pressure_reader,
    temperature_reader,
) -> tuple[float, float, float, float]:
    """Find finite compression and Mach ranges over an entire connected group."""
    compression_values: list[np.ndarray] = []
    mach_values: list[np.ndarray] = []
    for ref in refs:
        spec = label_records[ref].spec
        labels = label_reader.read_core(spec)
        component = labels == ref.local_label
        compression_values.append((-div_v_reader.read_core(spec))[component])
        mach_pressure = pressure_reader.read_core(spec)
        mach_temperature = temperature_reader.read_core(spec)
        combined_mach = np.maximum(
            np.nan_to_num(mach_pressure, nan=-np.inf),
            np.nan_to_num(mach_temperature, nan=-np.inf),
        )
        finite_mach = combined_mach[component]
        mach_values.append(finite_mach[np.isfinite(finite_mach)])

    compression_concat = (
        np.concatenate(compression_values)
        if compression_values
        else np.array([], dtype=float)
    )
    finite_compression = compression_concat[np.isfinite(compression_concat)]
    mach_concat = (
        np.concatenate(mach_values) if mach_values else np.array([], dtype=float)
    )
    compression_min = (
        float(np.min(finite_compression)) if finite_compression.size else 0.0
    )
    compression_max = (
        float(np.max(finite_compression)) if finite_compression.size else 0.0
    )
    mach_min = float(np.min(mach_concat)) if mach_concat.size else 0.0
    mach_max = float(np.max(mach_concat)) if mach_concat.size else 0.0
    return compression_min, compression_max, mach_min, mach_max


def _select_global_centers_simple_mode(
    *,
    layout: ChunkLayout,
    fields: ChunkedInputStore,
    label_store: ChunkedOutputStore,
    label_records: dict[_LabelRef, _LocalLabelRecord],
    groups: dict[_LabelRef, list[_LabelRef]],
    union_find: _UnionFind,
    center_score_mode: str,
) -> list[tuple[tuple[int, int, int], float, float]]:
    """Reproduce uniform SciPy maximum-position selection for a single score field.

    This legacy path assembles full-grid labels and scores and caches Mach
    chunks; its memory is not bounded by the chunk size or AMR budget.
    """
    label_reader = label_store.field_reader("labels")
    score_reader = fields.field_reader(_score_field_name(center_score_mode))
    pressure_reader = fields.field_reader("mach_pressure")
    temperature_reader = fields.field_reader("mach_temperature")

    group_ids = {
        root_ref: group_index for group_index, root_ref in enumerate(groups, start=1)
    }
    merged_labels = np.zeros(layout.shape, dtype=np.int32)
    score = np.zeros(layout.shape, dtype=float)
    mach_pressure_cache: dict[tuple[int, int, int], np.ndarray] = {}
    mach_temperature_cache: dict[tuple[int, int, int], np.ndarray] = {}
    specs_by_grid = {spec.grid_index: spec for spec in layout.specs}

    for spec in layout.specs:
        chunk_slice = tuple(
            slice(spec.start[axis], spec.stop[axis]) for axis in range(3)
        )
        labels = label_reader.read_core(spec)
        merged_chunk = np.zeros(spec.core_shape, dtype=np.int32)
        for local_label in np.unique(labels):
            if local_label <= 0:
                continue
            root_ref = union_find.find(_LabelRef(spec.grid_index, int(local_label)))
            merged_chunk[labels == local_label] = group_ids[root_ref]
        merged_labels[chunk_slice] = merged_chunk
        score[chunk_slice] = _simple_score_chunk(
            score_reader.read_core(spec), center_score_mode
        )
        mach_pressure_cache[spec.grid_index] = pressure_reader.read_core(spec)
        mach_temperature_cache[spec.grid_index] = temperature_reader.read_core(spec)

    positions = ndimage.maximum_position(
        score,
        labels=merged_labels,
        index=np.arange(1, len(groups) + 1),
    )
    centers: list[tuple[tuple[int, int, int], float, float]] = []
    for position in positions:
        global_index = tuple(int(value) for value in position)
        spec = _spec_for_global_index(layout, specs_by_grid, global_index)
        local_index = tuple(global_index[axis] - spec.start[axis] for axis in range(3))
        centers.append(
            (
                global_index,
                float(mach_pressure_cache[spec.grid_index][local_index]),
                float(mach_temperature_cache[spec.grid_index][local_index]),
            )
        )
    return centers


def _combined_component_score(
    component: np.ndarray,
    compression: np.ndarray,
    mach: np.ndarray,
    compression_min: float,
    compression_max: float,
    mach_min: float,
    mach_max: float,
) -> np.ndarray:
    """Normalize compression and Mach with shared group ranges, masking other cells."""
    score = np.full(component.shape, -np.inf, dtype=float)
    if not np.any(component):
        return score
    compression_norm = _normalize_values(
        compression[component], compression_min, compression_max
    )
    if mach_max <= mach_min:
        mach_norm = np.zeros(np.count_nonzero(component), dtype=float)
    else:
        mach_component = mach[component]
        mach_norm = np.zeros_like(mach_component, dtype=float)
        finite = np.isfinite(mach_component)
        mach_norm[finite] = (mach_component[finite] - mach_min) / (mach_max - mach_min)
    score[component] = compression_norm + mach_norm
    return score


def _normalize_values(values: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    """Scale to a supplied range; return ones when that range is degenerate."""
    if vmax <= vmin:
        return np.ones(values.shape, dtype=float)
    return (values - vmin) / (vmax - vmin)


def _score_field_name(center_score_mode: str) -> str:
    """Map a center-score mode to its stored field name."""
    if center_score_mode == "compression":
        return "div_v"
    if center_score_mode == "mach_pressure":
        return "mach_pressure"
    if center_score_mode == "mach_temperature":
        return "mach_temperature"
    raise ValueError(f"Unsupported simple center-score mode: {center_score_mode}")


def _simple_score_chunk(values: np.ndarray, center_score_mode: str) -> np.ndarray:
    """Convert divergence to compression or replace non-finite Mach scores."""
    if center_score_mode == "compression":
        return -values
    return np.nan_to_num(values, nan=-np.inf)


def _spec_for_global_index(
    layout: ChunkLayout,
    specs_by_grid: dict[tuple[int, int, int], ChunkSpec],
    global_index: tuple[int, int, int],
) -> ChunkSpec:
    """Locate the chunk owning a global cell index."""
    grid_index = []
    for axis, axis_slices in enumerate(layout.axis_slices):
        for axis_index, axis_slice in enumerate(axis_slices):
            if axis_slice.start <= global_index[axis] < axis_slice.stop:
                grid_index.append(axis_index)
                break
        else:
            raise ValueError(f"Index {global_index} falls outside chunk layout.")
    return specs_by_grid[tuple(grid_index)]


def _write_center_mask(
    output_store: ChunkedOutputStore,
    centers: list[tuple[tuple[int, int, int], float, float]],
) -> None:
    """Persist one boolean core mask per chunk, marking selected global centers."""
    centers_by_grid: dict[tuple[int, int, int], list[tuple[int, int, int]]] = {}
    for center_index, _, _ in centers:
        for spec in output_store.layout.specs:
            if all(
                spec.start[axis] <= center_index[axis] < spec.stop[axis]
                for axis in range(3)
            ):
                centers_by_grid.setdefault(spec.grid_index, []).append(
                    tuple(center_index[axis] - spec.start[axis] for axis in range(3))
                )
                break

    for spec in output_store.layout.specs:
        mask = np.zeros(spec.core_shape, dtype=bool)
        for local_index in centers_by_grid.get(spec.grid_index, []):
            mask[local_index] = True
        output_store.write_core("shock_mask", spec, mask)


def _finite_stat(values: np.ndarray, reducer) -> float:
    """Reduce finite values only; return NaN for an empty population."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return float("nan")
    return float(reducer(finite))


def _concat(values: list[np.ndarray]) -> np.ndarray:
    """Concatenate nonempty populations, preserving an empty float result when needed."""
    non_empty = [array for array in values if array.size > 0]
    if not non_empty:
        return np.array([], dtype=float)
    return np.concatenate(non_empty)


def _simple_score_field(
    div_v: np.ndarray,
    mach_fields: MachFields,
    mode: str,
) -> np.ndarray | None:
    """Return a single-field center score, or None for combined scoring."""
    if mode == "compression":
        return -div_v
    if mode == "mach_pressure":
        return np.nan_to_num(mach_fields.mach_pressure, nan=-np.inf)
    if mode == "mach_temperature":
        return np.nan_to_num(mach_fields.mach_temperature, nan=-np.inf)
    return None
