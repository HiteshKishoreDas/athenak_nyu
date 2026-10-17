"""Compatibility wrapper for the explicit chunked runner."""

from __future__ import annotations

from pathlib import Path

from .chunking import NpzChunkedInput, NpzChunkedOutput, run_chunked_shock_finder
from .config import ShockFinderConfig
from .finder import ProgressReporter


def run_chunked_npz_shock_finder(
    input_root: str | Path,
    output_root: str | Path,
    config: ShockFinderConfig,
    *,
    dx: float = 1.0,
    dy: float = 1.0,
    dz: float = 1.0,
    gamma: float = 5.0 / 3.0,
    progress: ProgressReporter = True,
) -> dict[str, object]:
    """Compatibility wrapper around `shockit.chunking.run_chunked_shock_finder()`."""

    config.validate()
    from .local_pass import build_reporter
    build_reporter(progress)
    input_store = NpzChunkedInput(input_root)
    output_store = NpzChunkedOutput(output_root, input_store.layout)
    return run_chunked_shock_finder(
        input_store,
        output_store,
        config=config,
        dx=dx,
        dy=dy,
        dz=dz,
        gamma=gamma,
        progress=progress,
    )
