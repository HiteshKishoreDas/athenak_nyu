"""Sub-block execution using the existing uniform local physics pass."""

from __future__ import annotations

from dataclasses import asdict
from copy import deepcopy
import math
from pathlib import Path

import numpy as np

from ..config import ShockFinderConfig
from ..fields import FluidCube
from ..provenance import build_provenance, atomic_write_json
from ..local_pass import (
    analyze_cube_pass,
    pass_fields,
    required_halo,
    build_reporter,
    ProgressReporter,
)
from .mesh import AMRMesh, AMRExecutionOptions, iter_tiles
from .output import (
    AMRShockFinderResult,
    StreamingSummary,
    mesh_metadata,
    TileIndexWriter,
)
from .reader import AMRReader


def find_amr(
    mesh: AMRMesh,
    config: ShockFinderConfig,
    *,
    output: str | Path,
    progress: ProgressReporter = True,
    execution: AMRExecutionOptions | None = None,
) -> AMRShockFinderResult:
    """Stream leaf cores into a new output directory using the shared local physics pass.

    Publish manifest.json only after all tiles and indices finish successfully.
    Failures may leave partial files but cannot publish a completed result.
    AMR center reduction is unsupported; output must be new or empty.
    Workspace planning and the primitive cache are bounded by execution options;
    caller-owned inputs, interpreter baseline and NPZ serialization transients
    are separate. The 2048-byte halo-cell allowance is not a peak RSS guarantee.
    Configuration is copied before callbacks; loader arrays are never mutated.
    """
    config.validate()
    config = deepcopy(config)  # Caller callbacks cannot alter the running detector.
    report = build_reporter(progress)  # Validate before output or input I/O.
    if config.reduce_to_centers:
        raise ValueError(
            "AMR requires reduce_to_centers=False; center reduction is not implemented"
        )
    if output is None:
        raise ValueError("AMR inputs require output=...")
    options = execution or AMRExecutionOptions()
    # Conservative common halo uses the finest leaf spacing. Coarse blocks
    # may over-read, but no tile can wrap at an artificial scratch boundary.
    halo = max(required_halo(config, block.spacing) for block in mesh.blocks)
    options.fitted_shape(halo)  # Validate budget before creating output.
    root = Path(output)
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"AMR output directory must be empty: {root}")
    root.mkdir(parents=True, exist_ok=True)
    (root / "tiles").mkdir()
    reader = AMRReader(mesh, options)
    summary = StreamingSummary()
    report(f"Processing {len(mesh.blocks)} AMR blocks in bounded tiles")
    previous_block_id = None
    with (
        (root / "tiles.jsonl").open("xb") as tile_index,
        (root / "blocks.jsonl").open("x") as block_index,
    ):
        writer = TileIndexWriter(tile_index, block_index)
        for number, tile in enumerate(iter_tiles(mesh, options, halo)):
            if tile.block_id != previous_block_id:
                report(f"AMR block {tile.block_id}")
                previous_block_id = tile.block_id
            block = mesh.block(tile.block_id)
            fields = reader.read_with_halo(tile.block_id, tile.slices, halo)
            cube = FluidCube(
                **fields,
                dx=block.spacing[0],
                dy=block.spacing[1],
                dz=block.spacing[2],
                gamma=mesh.gamma,
            )
            local_pass = analyze_cube_pass(cube, config)
            # Local operators wrap the scratch array; only halo-protected core
            # cells are retained, so scratch edges cannot become domain boundaries.
            crop = tuple(
                slice(halo, halo + b - a) for a, b in zip(tile.start, tile.stop)
            )
            cores = pass_fields(local_pass, crop)
            tile_path = f"tiles/{number:012d}.npz"
            np.savez(root / tile_path, **cores)
            writer.write(tile, tile_path)
            summary.update(cores, math.prod(block.spacing))
            # Do not retain the previous tile while constructing the next pass.
            del cores, local_pass, cube, fields
        writer.finish_block()
    totals = summary.finish(mesh, config)
    totals["provenance"] = build_provenance(
        config, gamma=mesh.gamma, mesh=mesh_metadata(mesh)
    )
    totals["provenance"]["memory_contract"] = {
        "workspace_budget_bytes": options.memory_budget_bytes,
        "primitive_cache_limit_bytes": reader.cache_limit,
        "planning_bytes_per_halo_cell": 2048,
        "excludes": ["caller-owned primitive arrays", "interpreter baseline", "NPZ transients"],
        "peak_rss_guarantee": False,
    }
    manifest = dict(
        version=2,
        axis_order="xyz",
        index="tiles.jsonl",
        block_index="blocks.jsonl",
        tile_index_bytes=(root / "tiles.jsonl").stat().st_size,
        mesh=mesh_metadata(mesh),
        configuration=asdict(config),
        execution=asdict(options),
        summary=totals,
    )
    atomic_write_json(root / "manifest.json", manifest)
    report("AMR output complete")
    result = AMRShockFinderResult.open(root)
    result.validate()
    result.mesh = mesh  # Live results retain the caller loader; reopened ones cannot.
    return result
