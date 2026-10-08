"""HDF5 adapter for loading FluidCube instances."""

from __future__ import annotations

from collections.abc import Sequence

import h5py
import numpy as np

from .fields import FluidCube

ALIASES = {
    "rho": ("rho", "density", "dens"),
    "pressure": ("pressure", "press", "prs"),
    "vx": ("vx", "vel1", "velocity_x"),
    "vy": ("vy", "vel2", "velocity_y"),
    "vz": ("vz", "vel3", "velocity_z"),
}


class AmbiguousFieldError(KeyError):
    """Raised when recursive basename resolution finds multiple datasets."""


def _dataset_exists(handle: h5py.File, path: str) -> bool:
    """Check that an HDF5 path names a dataset rather than a group."""
    try:
        return isinstance(handle[path], h5py.Dataset)
    except KeyError:
        return False


def _find_datasets_by_basename(handle: h5py.File, basename: str) -> list[str]:
    """Find all dataset paths sharing a basename so ambiguous aliases can be rejected."""
    found: list[str] = []

    def visitor(name: str, obj: h5py.Dataset) -> None:
        """Collect matching dataset paths while visiting the HDF5 tree."""
        if (
            isinstance(obj, h5py.Dataset)
            and name.rsplit("/", maxsplit=1)[-1] == basename
        ):
            found.append(name)

    handle.visititems(visitor)
    return found


def _resolve_field_name(
    handle: h5py.File, explicit_name: str | None, aliases: Sequence[str]
) -> str:
    """Resolve explicit paths or a unique basename; reject missing or ambiguous fields."""
    if explicit_name is not None:
        if _dataset_exists(handle, explicit_name):
            return explicit_name
        raise KeyError(f"explicit dataset {explicit_name!r} does not exist or is not a dataset")
    for alias in aliases:
        if _dataset_exists(handle, alias):
            return alias
    for candidate in aliases:
        matches = _find_datasets_by_basename(handle, candidate)
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            match_list = ", ".join(matches)
            raise AmbiguousFieldError(
                f"Ambiguous dataset basename '{candidate}' matched multiple paths: {match_list}. "
                "Pass an explicit dataset path such as 'group/field'."
            )
    raise KeyError(
        f"Could not find any of these fields: {', '.join(aliases)}"
    )


def load_fluid_cube_from_hdf5(
    filename: str,
    rho_field: str | None = None,
    pressure_field: str | None = None,
    vx_field: str | None = None,
    vy_field: str | None = None,
    vz_field: str | None = None,
    dx: float = 1.0,
    dy: float = 1.0,
    dz: float = 1.0,
    gamma: float = 5.0 / 3.0,
) -> FluidCube:
    """Load uniform HDF5 datasets, including nested paths, into memory.

    Arrays must already use xyz order; this adapter does not transpose native
    snapshot axes or convert internal energy into pressure. Omitted field names use
    aliases; explicitly supplied names must identify exact datasets. File handles
    close before return or error; primitive arrays are caller-owned in memory."""

    with h5py.File(filename, "r") as handle:
        rho_name = _resolve_field_name(handle, rho_field, ALIASES["rho"])
        pressure_name = _resolve_field_name(handle, pressure_field, ALIASES["pressure"])
        vx_name = _resolve_field_name(handle, vx_field, ALIASES["vx"])
        vy_name = _resolve_field_name(handle, vy_field, ALIASES["vy"])
        vz_name = _resolve_field_name(handle, vz_field, ALIASES["vz"])
        field_names = {
            "rho": rho_name,
            "pressure": pressure_name,
            "vx": vx_name,
            "vy": vy_name,
            "vz": vz_name,
        }
        field_chunk_shapes = {
            logical_name: _dataset_chunk_shape(handle[dataset_name])
            for logical_name, dataset_name in field_names.items()
        }
        return FluidCube(
            rho=np.asarray(handle[rho_name]),
            pressure=np.asarray(handle[pressure_name]),
            vx=np.asarray(handle[vx_name]),
            vy=np.asarray(handle[vy_name]),
            vz=np.asarray(handle[vz_name]),
            dx=dx,
            dy=dy,
            dz=dz,
            gamma=gamma,
            metadata={
                "source_file": filename,
                "rho_field": rho_name,
                "pressure_field": pressure_name,
                "vx_field": vx_name,
                "vy_field": vy_name,
                "vz_field": vz_name,
                "source_chunk_shape": _shared_chunk_shape(field_chunk_shapes),
                "field_chunk_shapes": field_chunk_shapes,
            },
        )


def _dataset_chunk_shape(dataset: h5py.Dataset) -> tuple[int, int, int] | None:
    """Read valid three-dimensional HDF5 chunk metadata when available."""
    if dataset.chunks is None:
        return None
    return tuple(int(value) for value in dataset.chunks)


def _shared_chunk_shape(
    field_chunk_shapes: dict[str, tuple[int, int, int] | None],
) -> tuple[int, int, int] | None:
    """Return a common stored chunk shape when every primitive field agrees."""
    chunk_shapes = {
        chunk_shape
        for chunk_shape in field_chunk_shapes.values()
    }
    if len(chunk_shapes) != 1:
        return None
    return next(iter(chunk_shapes))
