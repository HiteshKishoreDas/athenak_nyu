"""Main shock finding API."""

from __future__ import annotations

from pathlib import Path
from copy import deepcopy
from typing import TYPE_CHECKING

import numpy as np

from .config import ShockFinderConfig
from .fields import FluidCube
from .masks import reduce_to_centers
from .result import ChunkedShockFinderResult, RESULT_FIELD_NAMES, ShockFinderResult
from .storage import ChunkedFieldReference
from .summary import make_summary

if TYPE_CHECKING:
    from .amr import AMRMesh, AMRExecutionOptions, AMRShockFinderResult

# Compatibility imports retained for existing callers.
from .local_pass import (
    FinderPass as _FinderPass,
    analyze_cube_pass,
    build_reporter as _build_reporter,
    _noop_report,
    ProgressCallback,
    ProgressReporter,
    pass_fields,
)


class ShockFinder:
    """Skillman-style shock finder for uniform fluid cubes and AMR meshes."""

    def __init__(self, config: ShockFinderConfig | None = None) -> None:
        """Store the supplied configuration or use package defaults."""
        self.config = config or ShockFinderConfig()

    def find(
        self,
        cube: FluidCube | AMRMesh,
        progress: ProgressReporter = True,
        output: str | Path | None = None,
        *,
        execution: AMRExecutionOptions | None = None,
    ) -> ShockFinderResult | ChunkedShockFinderResult | AMRShockFinderResult:
        """Find shocks using a validated snapshot of the mutable configuration.

        Uniform arrays use xyz cell order and periodic boundaries. Chunk stores
        require a new output directory; AMR uses its declared mesh boundaries.
        Legacy chunk center reduction allocates at least 28 bytes/global cell
        before temporaries; exact medians retain detections. These are not peak
        RSS estimates. AMR does not support center reduction.
        """

        self.config.validate()
        config = deepcopy(self.config)
        report = _build_reporter(progress)
        from .amr.mesh import AMRMesh

        if isinstance(cube, AMRMesh):
            from .amr.runner import find_amr

            return find_amr(
                cube, config, output=output, progress=progress, execution=execution
            )
        if execution is not None:
            raise ValueError("execution= is only supported for AMR inputs")
        if not isinstance(cube, FluidCube):
            raise TypeError(
                "ShockFinder.find() expects a FluidCube or AMRMesh instance."
            )

        if cube.is_chunked:
            return self._find_chunked(cube, output=output, progress=progress, config=config)
        if output is not None:
            raise ValueError(
                "output= is only used for chunked FluidCube inputs; "
                "use save_result_hdf5() for in-memory results."
            )
        return self._find_in_memory(cube, report=report, config=config)

    def _find_in_memory(
        self,
        cube: FluidCube,
        *,
        report: ProgressCallback,
        config: ShockFinderConfig,
    ) -> ShockFinderResult:
        """Run uniform physics, optional center reduction and in-memory result assembly."""
        cube_pass = analyze_cube_pass(cube, config, report=report)

        if config.reduce_to_centers:
            report("[8/9] Reducing to center cells")
            shock_mask, _ = reduce_to_centers(
                cube_pass.full_shock_mask,
                cube_pass.derived.div_v,
                cube_pass.mach_fields,
                center_score_mode=config.center_score,
            )
            report("[9/9] Building summary")
        else:
            shock_mask = cube_pass.full_shock_mask.copy()
            report("[8/8] Building summary")

        summary = make_summary(
            cube=cube,
            shock_zone_mask=cube_pass.shock_zone_mask,
            full_shock_mask=cube_pass.full_shock_mask,
            shock_mask=shock_mask,
            mach_fields=cube_pass.mach_fields,
            config=config,
        )
        return ShockFinderResult(
            **pass_fields(cube_pass, shock_mask=shock_mask),
            summary=summary,
        )

    def _find_chunked(
        self,
        cube: FluidCube,
        *,
        output: str | Path | None,
        progress: ProgressReporter,
        config: ShockFinderConfig,
    ) -> ChunkedShockFinderResult:
        """Run chunked storage through the shared physics and chunk finalization paths."""
        if output is None:
            raise ValueError(
                "Chunked FluidCube requires output=... because result fields are written chunked on disk."
            )

        from .chunking import (
            NpzChunkedInput,
            NpzChunkedOutput,
            run_chunked_shock_finder,
        )

        input_store = NpzChunkedInput.from_fluid_cube(cube)
        output_root = Path(output)
        output_store = NpzChunkedOutput(output_root, input_store.layout)
        summary = run_chunked_shock_finder(
            input_store,
            output_store,
            config=config,
            dx=cube.dx,
            dy=cube.dy,
            dz=cube.dz,
            gamma=cube.gamma,
            progress=progress,
        )
        return _build_chunked_result(output_root, input_store.layout, summary)


def _build_chunked_result(
    output_root: Path,
    layout,
    summary: dict[str, object],
) -> ChunkedShockFinderResult:
    """Expose persisted fields as lazy references alongside their summary."""
    field_refs = {
        field_name: ChunkedFieldReference(path=output_root / field_name, layout=layout)
        for field_name in RESULT_FIELD_NAMES
    }
    return ChunkedShockFinderResult(
        output_root=output_root,
        shock_mask=field_refs["shock_mask"],
        shock_zone_mask=field_refs["shock_zone_mask"],
        mach_temperature=field_refs["mach_temperature"],
        mach_pressure=field_refs["mach_pressure"],
        compression=field_refs["compression"],
        div_v=field_refs["div_v"],
        temperature=field_refs["temperature"],
        entropy=field_refs["entropy"],
        temperature_jump=field_refs["temperature_jump"],
        pressure_jump=field_refs["pressure_jump"],
        density_jump=field_refs["density_jump"],
        normal_x=field_refs["normal_x"],
        normal_y=field_refs["normal_y"],
        normal_z=field_refs["normal_z"],
        summary=summary,
        full_shock_mask=field_refs["full_shock_mask"],
    )
