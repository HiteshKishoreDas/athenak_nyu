"""Output helpers for shock-finder results."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

import json

import h5py
import numpy as np

from .result import ShockFinderResult
from .provenance import json_safe


def _package_version() -> str:
    """Read installed package metadata, or return 0+unknown when unavailable."""
    try:
        return version("shockit")
    except PackageNotFoundError:  # pragma: no cover - editable fallback
        return "0+unknown"


MASK_DESCRIPTIONS = {
    "shock_zone_mask": "Compression/gradient candidate mask before jump and Mach filtering.",
    "full_shock_mask": "Unreduced shock cells after jump and Mach filtering.",
    "shock_mask": "Final public shock mask; may be center-reduced.",
}

MACH_DESCRIPTIONS = {
    "mach_pressure": "Local Mach estimate from pressure jump.",
    "mach_temperature": "Local Mach estimate from temperature jump.",
}


def _write_mask_dataset(handle: h5py.File, name: str, mask: np.ndarray) -> None:
    """Store a boolean mask with its interpretation attached as metadata."""
    dataset = handle.create_dataset(name, data=mask.astype(bool), dtype=np.bool_)
    dataset.attrs["semantic_type"] = "boolean_mask"
    dataset.attrs["stored_dtype"] = str(dataset.dtype)
    dataset.attrs["description"] = MASK_DESCRIPTIONS[name]


def _write_field_dataset(
    handle: h5py.File, name: str, values: np.ndarray
) -> h5py.Dataset:
    """Create a numerical dataset; the caller may attach field-specific metadata."""
    return handle.create_dataset(name, data=values)


def save_result_hdf5(filename: str, result: ShockFinderResult) -> None:
    """Persist in-memory outputs to HDF5, replacing an existing destination file."""

    with h5py.File(filename, "w") as handle:
        _write_mask_dataset(handle, "shock_mask", result.shock_mask)
        _write_mask_dataset(handle, "shock_zone_mask", result.shock_zone_mask)
        if result.full_shock_mask is not None:
            _write_mask_dataset(handle, "full_shock_mask", result.full_shock_mask)
        mach_temperature = _write_field_dataset(
            handle, "mach_temperature", result.mach_temperature
        )
        mach_pressure = _write_field_dataset(
            handle, "mach_pressure", result.mach_pressure
        )
        mach_temperature.attrs["description"] = MACH_DESCRIPTIONS["mach_temperature"]
        mach_pressure.attrs["description"] = MACH_DESCRIPTIONS["mach_pressure"]
        _write_field_dataset(handle, "compression", result.compression)
        _write_field_dataset(handle, "div_v", result.div_v)
        _write_field_dataset(handle, "temperature", result.temperature)
        _write_field_dataset(handle, "entropy", result.entropy)
        _write_field_dataset(handle, "temperature_jump", result.temperature_jump)
        _write_field_dataset(handle, "pressure_jump", result.pressure_jump)
        _write_field_dataset(handle, "density_jump", result.density_jump)
        _write_field_dataset(handle, "normal_x", result.normal_x)
        _write_field_dataset(handle, "normal_y", result.normal_y)
        _write_field_dataset(handle, "normal_z", result.normal_z)
        handle.attrs["package_name"] = "shockit"
        handle.attrs["package_version"] = _package_version()
        handle.attrs["mach_field_note"] = (
            "Mach fields are local cell estimates from sampled jump conditions."
        )
        handle.attrs["summary_json"] = json.dumps(json_safe(result.summary), allow_nan=False)
        handle.attrs["provenance_json"] = json.dumps(json_safe(result.summary.get("provenance")), allow_nan=False)
        for key, value in result.summary.items():
            if value is None or isinstance(value, dict):
                continue
            if isinstance(value, tuple):
                value = np.asarray(value)
            handle.attrs[key] = value
