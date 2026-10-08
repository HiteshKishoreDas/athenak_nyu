"""Portable, tile-indexed NPZ results with bounded region reads.

Version 2 adds one byte span per block. Version 1 remains readable by scanning
its index. Metadata is O(blocks); tile entries are never collected in memory.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import asdict
import json
import math
from pathlib import Path
import sqlite3
import tempfile

import numpy as np

from ..result import RESULT_FIELD_NAMES
from .mesh import AMRBlock, AMRMesh, AMRTile, validate_slices


def mesh_metadata(mesh):
    """Serialize geometry and physics metadata without retaining the primitive loader."""
    return dict(
        blocks=[asdict(b) for b in mesh.blocks],
        domain_bounds=mesh.domain_bounds,
        root_shape=mesh.root_shape,
        gamma=mesh.gamma,
        boundaries=mesh.boundaries,
    )


def _no_input(*args):
    """Reject primitive reads from a reopened output-only mesh."""
    raise ValueError("reopened result metadata has no primitive input loader")


class TileIndexWriter:
    """Append contiguous block spans, retaining only the current block's counters."""

    def __init__(self, index, blocks):
        """Track one block span; index must be a binary stream so offsets count bytes."""
        self.index = index
        self.blocks = blocks
        self.block_id = None
        self.start = 0
        self.count = 0

    def finish_block(self):
        """Append the current block span after its last tile; call once at completion."""
        if self.block_id is not None:
            entry = dict(
                block_id=self.block_id,
                start=self.start,
                stop=self.index.tell(),
                count=self.count,
            )
            self.blocks.write(json.dumps(entry) + "\n")

    def write(self, tile: AMRTile, path: str) -> None:
        """Append one tile row; consecutive tiles must be grouped by block ID."""
        if tile.block_id != self.block_id:
            self.finish_block()
            self.block_id = tile.block_id
            self.start = self.index.tell()
            self.count = 0
        row = dict(block_id=tile.block_id, start=tile.start, stop=tile.stop, path=path)
        self.index.write((json.dumps(row) + "\n").encode("utf-8"))
        self.count += 1


class AMRFieldReference:
    """One lazy output field; region coordinates are block-local x,y,z slices."""

    def __init__(self, result: AMRShockFinderResult, name: str):
        """Attach a named field without opening any tile archive."""
        self.result = result
        self.name = name

    def _read(self, row) -> np.ndarray:
        """Read and validate this field from one indexed tile."""
        return self.result._read_fields(row, (self.name,))[self.name]

    def iter_tiles(self) -> Iterator[tuple[AMRTile, np.ndarray]]:
        """Yield one stored core at a time without assembling blocks."""
        for tile, fields in self.result.iter_tiles((self.name,)):
            yield tile, fields[self.name]

    def read_region(
        self,
        block_id: str,
        slices: tuple[slice, slice, slice],
        *,
        max_bytes: int = 128 * 1024**2,
    ) -> np.ndarray:
        """Assemble a bounded block-local region from intersecting stored tiles.

        max_bytes bounds the returned array plus coverage mask. Each intersecting
        NPZ field is loaded as a whole stored tile, adding transient memory.
        """
        block = self.result.mesh.block(block_id)
        slices = validate_slices(slices, block.shape)
        shape = tuple(s.stop - s.start for s in slices)
        dtype = bool if self.name.endswith("mask") else float
        if math.prod(shape) * (np.dtype(dtype).itemsize + 1) > max_bytes:
            raise ValueError(
                "memory limit exceeded: iterate output tiles or request a smaller region"
            )
        output = np.empty(shape, dtype=dtype)
        filled = np.zeros(shape, bool)
        for row in self.result._rows(block_id):
            start = tuple(max(s.start, a) for s, a in zip(slices, row["start"]))
            stop = tuple(min(s.stop, a) for s, a in zip(slices, row["stop"]))
            if any(a >= b for a, b in zip(start, stop)):
                continue
            source_slices = tuple(
                slice(a - origin, b - origin)
                for a, b, origin in zip(start, stop, row["start"])
            )
            output_slices = tuple(
                slice(a - s.start, b - s.start) for a, b, s in zip(start, stop, slices)
            )
            if filled[output_slices].any():
                raise ValueError("overlapping output tile cores")
            output[output_slices] = self._read(row)[source_slices]
            filled[output_slices] = True
        if not filled.all():
            raise ValueError("missing output tile coverage")
        return output


class AMRShockFinderResult:
    """Lazy AMR results. No primitive loader or field reads needed to reopen.

    Use iter_tiles(fields=...) to load several fields with one archive open.
    Each yielded dictionary holds one tile, never a complete block or snapshot.
    """

    shock_mask: AMRFieldReference
    shock_zone_mask: AMRFieldReference
    full_shock_mask: AMRFieldReference
    mach_temperature: AMRFieldReference
    mach_pressure: AMRFieldReference
    compression: AMRFieldReference
    div_v: AMRFieldReference
    temperature: AMRFieldReference
    entropy: AMRFieldReference
    temperature_jump: AMRFieldReference
    pressure_jump: AMRFieldReference
    density_jump: AMRFieldReference
    normal_x: AMRFieldReference
    normal_y: AMRFieldReference
    normal_z: AMRFieldReference

    def __init__(
        self,
        root,
        mesh: AMRMesh,
        summary: dict,
        configuration: dict,
        block_index: dict | None = None,
    ):
        """Attach metadata and lazy field references; perform no field reads."""
        self.root = Path(root)
        self.output_root = self.root
        self.mesh = mesh
        self.summary = summary
        self.configuration = configuration
        self._block_index = block_index
        for name in RESULT_FIELD_NAMES:
            setattr(self, name, AMRFieldReference(self, name))

    @classmethod
    def open(cls, path: str | Path) -> AMRShockFinderResult:
        """Reopen completed metadata and validate block spans without loading tile fields."""
        root = Path(path)
        with (root / "manifest.json").open() as stream:
            manifest = json.load(stream)
        if not isinstance(manifest, dict):
            raise ValueError("invalid AMR manifest object")
        version = manifest.get("version")
        if type(version) is not int or version not in (1, 2):
            raise ValueError("unsupported AMR manifest version")
        if (
            manifest.get("axis_order") != "xyz"
            or manifest.get("index") != "tiles.jsonl"
        ):
            raise ValueError("invalid AMR manifest")
        data = manifest["mesh"]
        mesh = AMRMesh(
            [AMRBlock(**b) for b in data["blocks"]],
            data["domain_bounds"],
            data["root_shape"],
            _no_input,
            gamma=data["gamma"],
            boundaries=data["boundaries"],
        )
        blocks = None
        if version == 2:
            if manifest.get("block_index") != "blocks.jsonl":
                raise ValueError("invalid block index manifest")
            size = (root / "tiles.jsonl").stat().st_size
            if manifest.get("tile_index_bytes") != size:
                raise ValueError("tile index size mismatch")
            blocks = {}
            end = 0
            with (root / "blocks.jsonl").open() as stream:
                for line in stream:
                    entry = json.loads(line)
                    if not isinstance(entry, dict):
                        raise ValueError("invalid block index object")
                    block_id = entry.get("block_id")
                    try:
                        mesh.block(block_id)
                    except (KeyError, TypeError) as exc:
                        raise ValueError("unknown block in block index") from exc
                    values = [entry.get(k) for k in ("start", "stop", "count")]
                    if any(type(v) is not int for v in values):
                        raise ValueError("invalid block index integers")
                    start, stop, count = values
                    if (
                        block_id in blocks
                        or start != end
                        or not start < stop <= size
                        or count < 1
                    ):
                        raise ValueError("invalid block index span or duplicate block")
                    blocks[block_id] = (start, stop, count)
                    end = stop
            if len(blocks) != len(mesh.blocks) or end != size:
                raise ValueError("incomplete block index")
        return cls(root, mesh, manifest["summary"], manifest["configuration"], blocks)

    def validate(self, *, fields: bool = False) -> None:
        """Verify complete index geometry; optionally verify every tile payload.

        Scan V1 rows or V2 byte spans without retaining a tile list. Each row
        must name an existing file and a unique output-relative path. A temporary
        disk SQLite index rejects intersecting half-open integer cores; containment,
        disjointness and exact summed volume prove coverage of each block.

        `fields=False` reads metadata and file status only. `fields=True` also
        loads one tile's required xyz arrays at a time, including NPZ transients.
        The database cache is bounded to 2 MiB; disk usage grows with tile count.
        The exact integer B-tree overlap query can take quadratic time for
        adversarial tile layouts; it never requires a dense coverage array.
        Files/connections are closed and removed even on failure. Raise ValueError
        for invalid metadata/payload or missing files; OS storage failures propagate.
        Stored summary counts are not independently recounted. This explicit full
        scan does not change the lazy, visited-row guarantees of iter_tiles.
        """
        if type(fields) is not bool:
            raise TypeError("fields must be boolean")
        # Reopen to recheck manifest byte lengths/spans, including edits since open.
        checked = type(self).open(self.root)
        with tempfile.TemporaryDirectory(prefix="shockit-amr-index-") as temporary:
            connection = sqlite3.connect(str(Path(temporary) / "geometry.sqlite"))
            try:
                connection.execute("PRAGMA cache_size=-2048")
                connection.execute("PRAGMA temp_store=FILE")
                connection.execute(
                    "CREATE TABLE cores (block TEXT, path TEXT UNIQUE, "
                    "x0 INTEGER, x1 INTEGER, y0 INTEGER, y1 INTEGER, "
                    "z0 INTEGER, z1 INTEGER)"
                )
                connection.execute("CREATE INDEX core_bounds ON cores(block, x0, x1)")
                volumes = {block.block_id: 0 for block in checked.mesh.blocks}
                for row in checked._rows():
                    path = checked.root / row["path"]
                    if not path.is_file():
                        raise ValueError(f"missing tile file: {row['path']}")
                    x0, y0, z0 = row["start"]
                    x1, y1, z1 = row["stop"]
                    block_id = row["block_id"]
                    overlap = connection.execute(
                        "SELECT 1 FROM cores WHERE block=? AND x0<? AND x1>? "
                        "AND y0<? AND y1>? AND z0<? AND z1>? LIMIT 1",
                        (block_id, x1, x0, y1, y0, z1, z0),
                    ).fetchone()
                    if overlap:
                        raise ValueError(f"overlapping output tile cores in block {block_id}")
                    try:
                        connection.execute(
                            "INSERT INTO cores VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                            (block_id, str(Path(row["path"])), x0, x1, y0, y1, z0, z1),
                        )
                    except sqlite3.IntegrityError as exc:
                        raise ValueError(f"duplicate tile path: {row['path']}") from exc
                    volumes[block_id] += (x1 - x0) * (y1 - y0) * (z1 - z0)
                    if fields:
                        checked._read_fields(row, RESULT_FIELD_NAMES)
                for block in checked.mesh.blocks:
                    if volumes[block.block_id] != math.prod(block.shape):
                        raise ValueError(f"missing output tile coverage in block {block.block_id}")
            finally:
                connection.close()

    def _read_fields(self, row, names: Sequence[str]) -> dict[str, np.ndarray]:
        """Open one archive, validate selected field shapes/dtypes, and close on exit."""
        path = self.root / row["path"]
        expected = tuple(b - a for a, b in zip(row["start"], row["stop"]))
        fields = {}
        with np.load(path, allow_pickle=False) as archive:
            for name in names:
                if name not in archive:
                    raise ValueError(f"missing field {name} in {path.name}")
                values = archive[name]
                if values.shape != expected:
                    raise ValueError(f"invalid tile shape for {name}")
                if name.endswith("mask") and values.dtype != bool:
                    raise ValueError("mask tile must have boolean dtype")
                if not name.endswith("mask") and values.dtype.kind not in "iuf":
                    raise ValueError(
                        f"tile field {name} must have real numeric dtype in {path.name}"
                    )
                fields[name] = values
        return fields

    def iter_tiles(
        self, fields: Sequence[str] = RESULT_FIELD_NAMES, *, block_id: str | None = None
    ) -> Iterator[tuple[AMRTile, dict[str, np.ndarray]]]:
        """Yield requested fields, one NPZ open per tile. Optional block filter.

        Only visited rows and selected payloads are checked; overlaps/gaps in the
        complete index require validate(). Memory scales with selected fields in
        one stored tile. For a smaller returned array use an individual field's bounded read_region method.
        """
        names = tuple(fields)
        if (
            not names
            or len(set(names)) != len(names)
            or any(n not in RESULT_FIELD_NAMES for n in names)
        ):
            raise ValueError("fields must be unique known result field names")
        for row in self._rows(block_id):
            tile = AMRTile(row["block_id"], tuple(row["start"]), tuple(row["stop"]))
            yield tile, self._read_fields(row, names)

    def _validate_row(self, row):
        """Check a tile's block-local bounds and portable output-relative path."""
        if not isinstance(row, dict):
            raise ValueError("invalid tile index object")
        try:
            block = self.mesh.block(row["block_id"])
            if len(row["start"]) != 3 or len(row["stop"]) != 3:
                raise ValueError("tile bounds must contain exactly three integers")
            slices = tuple(slice(a, b) for a, b in zip(row["start"], row["stop"]))
            validate_slices(slices, block.shape)
            path = Path(row["path"])
        except (KeyError, TypeError) as exc:
            raise ValueError("malformed or unknown block in tile index") from exc
        if not path.parts or path.is_absolute() or ".." in path.parts:
            raise ValueError("tile paths must be relative to output directory")
        return row

    def _rows(self, block_id: str | None = None):
        """Stream validated tile rows, seeking directly to a block in version 2."""
        if block_id is not None:
            self.mesh.block(block_id)
        with (self.root / "tiles.jsonl").open("rb") as stream:
            if self._block_index is not None:
                spans = (
                    self._block_index.items()
                    if block_id is None
                    else ((block_id, self._block_index[block_id]),)
                )
                for expected_id, (start, stop, count) in spans:
                    stream.seek(start)
                    for _ in range(count):
                        line = stream.readline()
                        if not line or stream.tell() > stop:
                            raise ValueError("truncated or invalid block index span")
                        row = self._validate_row(json.loads(line))
                        if row["block_id"] != expected_id:
                            raise ValueError("block index points to another block")
                        yield row
                    if stream.tell() != stop:
                        raise ValueError("block index count mismatch")
            else:
                # Legacy results retain a streaming full-index scan.
                for line in stream:
                    row = self._validate_row(json.loads(line))
                    if block_id is None or row["block_id"] == block_id:
                        yield row


class StreamingSummary:
    """Only scalar accumulators: no retained detected-cell arrays or quantiles."""

    def __init__(self):
        """Initialize integer counts, physical volumes and optional finite extrema."""
        self.values = {"total_cells": 0, "total_volume": 0.0}
        for prefix in ("shock_zone", "full_shock", "shock"):
            self.values[prefix + "_cells"] = 0
            self.values[prefix + "_volume"] = 0.0
        for prefix in ("", "full_"):
            for name in ("mach_pressure", "mach_temperature"):
                self.values[prefix + name + "_min"] = None
                self.values[prefix + name + "_max"] = None

    def update(self, fields, cell_volume):
        """Accumulate one core using its level-dependent cell volume and each field's mask."""
        count = math.prod(fields["shock_mask"].shape)
        self.values["total_cells"] += count
        self.values["total_volume"] += count * cell_volume
        for prefix in ("shock_zone", "full_shock", "shock"):
            selected_cells = int(np.count_nonzero(fields[prefix + "_mask"]))
            self.values[prefix + "_cells"] += selected_cells
            self.values[prefix + "_volume"] += selected_cells * cell_volume
        for prefix, mask in (
            ("", fields["shock_mask"]),
            ("full_", fields["full_shock_mask"]),
        ):
            for name in ("mach_pressure", "mach_temperature"):
                values = fields[name][mask & np.isfinite(fields[name])]
                if not values.size:
                    continue
                for suffix, combine_extrema in (("_min", min), ("_max", max)):
                    key = prefix + name + suffix
                    value = float(
                        np.min(values) if suffix == "_min" else np.max(values)
                    )
                    self.values[key] = (
                        value
                        if self.values[key] is None
                        else combine_extrema(self.values[key], value)
                    )

    def finish(self, mesh, config):
        """Return count/volume fractions and configuration; absent extrema remain JSON null."""
        summary = dict(self.values)
        for prefix in ("shock_zone", "full_shock", "shock"):
            summary[prefix + "_fraction"] = (
                summary[prefix + "_cells"] / summary["total_cells"]
            )
            summary[prefix + "_volume_fraction"] = (
                summary[prefix + "_volume"] / summary["total_volume"]
            )
        summary.update(
            gamma=mesh.gamma,
            min_mach=config.min_mach,
            shock_width_cells=config.shock_width_cells,
            reduce_to_centers=False,
            sampling_method=config.sampling_method,
            mask_semantics="shock_mask equals full_shock_mask; AMR center reduction is unsupported",
        )
        return summary
