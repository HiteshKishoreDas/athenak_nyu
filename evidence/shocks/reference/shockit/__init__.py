"""Public package exports for shockit."""

from importlib.metadata import PackageNotFoundError, version

from .config import ShockFinderConfig
from .fields import FluidCube
from .finder import ShockFinder
from .io_hdf5 import load_fluid_cube_from_hdf5
from .output import save_result_hdf5
from .result import ShockFinderResult
from .amr import (
    AMRBlock,
    AMRMesh,
    AMRExecutionOptions,
    AMRShockFinderResult,
    amr_mesh_from_athenak,
)

try:
    __version__ = version("shockit")
except PackageNotFoundError:  # pragma: no cover - fallback for editable local imports
    __version__ = "0+unknown"

__all__ = [
    "AMRBlock",
    "AMRMesh",
    "AMRExecutionOptions",
    "AMRShockFinderResult",
    "amr_mesh_from_athenak",
    "FluidCube",
    "ShockFinder",
    "ShockFinderConfig",
    "ShockFinderResult",
    "load_fluid_cube_from_hdf5",
    "save_result_hdf5",
    "__version__",
]
