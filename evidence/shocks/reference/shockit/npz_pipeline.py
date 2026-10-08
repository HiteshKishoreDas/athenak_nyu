"""Legacy `.npz` pipeline helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .chunking import join_chunked_field
from .config import ShockFinderConfig
from .fields import FluidCube
from .finder import ShockFinder
from .provenance import atomic_write_json
from .local_pass import build_reporter


@dataclass(frozen=True)
class NpzLayout:
    """Resolved input layout for a simulation data directory."""

    kind: str
    input_path: Path


def detect_npz_layout(root: str | Path) -> NpzLayout:
    """Detect whether a directory contains plain or chunked shock-finder inputs."""

    root_path = Path(root)
    chunked_root = root_path / "chunked_npz"
    if _has_chunked_inputs(chunked_root):
        return NpzLayout(kind="chunked", input_path=chunked_root)
    if _has_plain_inputs(root_path):
        return NpzLayout(kind="plain", input_path=root_path)
    if _has_chunked_inputs(root_path):
        return NpzLayout(kind="chunked", input_path=root_path)
    raise FileNotFoundError(
        f"Could not find plain .npz inputs or a chunked_npz layout under {root_path}."
    )


def run_npz_shock_finder(
    root: str | Path,
    config: ShockFinderConfig,
    *,
    dx: float = 1.0,
    dy: float = 1.0,
    dz: float = 1.0,
    gamma: float = 5.0 / 3.0,
    progress: bool = True,
) -> dict[str, Any]:
    """Run the shock finder on either plain or chunked `.npz` inputs."""

    config.validate()
    build_reporter(progress)
    layout = detect_npz_layout(root)
    if layout.kind == "chunked":
        cube = FluidCube(
            rho=layout.input_path / "rho",
            pressure=layout.input_path / "prs",
            vx=layout.input_path / "v1",
            vy=layout.input_path / "v2",
            vz=layout.input_path / "v3",
            dx=dx,
            dy=dy,
            dz=dz,
            gamma=gamma,
        )
        result = ShockFinder(config).find(
            cube,
            output=layout.input_path.parent / "shock_chunked_npz",
            progress=progress,
        )
        return result.summary

    _check_plain_destination(layout.input_path)
    cube = _load_plain_cube(
        layout.input_path,
        dx=dx,
        dy=dy,
        dz=dz,
        gamma=gamma,
    )
    result = ShockFinder(config).find(cube, progress=progress)
    _save_plain_outputs(layout.input_path, result)
    return result.summary


def _has_plain_inputs(root: Path) -> bool:
    """Check for the five conventional single-file primitive inputs."""
    return all(
        (root / filename).exists()
        for filename in ("rho.npz", "prs.npz", "v1.npz", "v2.npz", "v3.npz")
    )


def _has_chunked_inputs(root: Path) -> bool:
    """Check for the five conventional chunked primitive directories."""
    return all((root / field).is_dir() for field in ("rho", "prs", "v1", "v2", "v3"))


def _load_plain_cube(
    root: Path,
    *,
    dx: float,
    dy: float,
    dz: float,
    gamma: float,
) -> FluidCube:
    """Load conventional arr_0 primitive arrays using the caller's grid metadata."""
    fields = {}
    for name, filename in (("rho", "rho"), ("pressure", "prs"), ("vx", "v1"),
                           ("vy", "v2"), ("vz", "v3")):
        with np.load(root / f"{filename}.npz") as archive:
            fields[name] = archive["arr_0"]
    return FluidCube(**fields, dx=dx, dy=dy, dz=dz, gamma=gamma,
                     metadata={"source_directory": str(root)})


def _check_plain_destination(root: Path) -> None:
    """Legacy outputs share the input directory; reserve only their result names."""
    for name in ("mask.npz", "mask_full.npz", "mach_T.npz", "mach_prs.npz", "summary.json", "incomplete.json"):
        if (root / name).exists():
            raise ValueError(f"existing result {root / name}; use fresh inputs/output location")


def _save_plain_outputs(root: Path, result) -> None:
    """Persist the conventional masks, Mach fields and summary filenames."""
    _check_plain_destination(root)
    atomic_write_json(root / "incomplete.json", {"operation": "plain_npz_write"})
    np.savez(root / "mask.npz", result.shock_mask)
    np.savez(root / "mask_full.npz", result.full_shock_mask)
    np.savez(root / "mach_T.npz", result.mach_temperature)
    np.savez(root / "mach_prs.npz", result.mach_pressure)
    atomic_write_json(root / "summary.json", result.summary)
    (root / "incomplete.json").unlink()


__all__ = [
    "NpzLayout",
    "detect_npz_layout",
    "run_npz_shock_finder",
    "join_chunked_field",
]
