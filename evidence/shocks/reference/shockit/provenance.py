"""Small JSON/provenance helpers shared by result formats.

Provenance describes the executed settings, not a certificate of scientific accuracy.
Unknown source information stays null; no filesystem scan or input hashing is implied.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict
from importlib.metadata import PackageNotFoundError, version
import json
import math
from pathlib import Path
import tempfile

import numpy as np


def json_safe(value):
    """Copy metadata into JSON values; nonfinite statistics become null.

    Supports mappings, sequences, NumPy scalars/arrays and paths. Reject unknown
    objects rather than silently persisting an unrepeatable string representation.
    """
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    raise TypeError(f"unsupported provenance value: {type(value).__name__}")


def atomic_write_json(path, value):
    """Publish JSON after a complete same-directory write; remove failed scratch files."""
    path = Path(path)
    payload = json.dumps(json_safe(value), indent=2, sort_keys=True, allow_nan=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent,
                                         prefix=f".{path.name}.", delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(payload + "\n")
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def build_provenance(config, *, gamma, spacing=None, source_metadata=None, mesh=None):
    """Snapshot detector settings and caller-supplied geometry/source metadata.

    ``spacing`` is xyz for a uniform grid; AMR uses ``mesh`` instead. Paths identify
    sources but do not fingerprint their content. Constants omitted from ideal-gas
    proxies and unknown source units are not guessed.
    """
    try:
        package_version = version("shockit")
    except PackageNotFoundError:
        package_version = "0+unknown"
    metadata = source_metadata or {}
    return json_safe({
        "schema_version": 1,
        "detection_config": asdict(config),
        "postprocessing_config": {"reduce_to_centers": config.reduce_to_centers,
                                   "center_score": config.center_score},
        "gamma": gamma, "axis_order": "xyz", "spacing": spacing, "mesh": mesh,
        "boundaries": "periodic" if mesh is None else mesh.get("boundaries"),
        "sampling_convention": "index-space offsets; width is half endpoint separation",
        "mask_conventions": {
            "shock_zone_mask": "compression/gradient candidates",
            "full_shock_mask": "unreduced accepted cells",
            "shock_mask": "accepted cells or one representative per component",
        },
        "package_version": package_version,
        "source_metadata": metadata,
        "source_units": metadata.get("units"),
        "source_revision": metadata.get("source_revision"),
    })


def update_postprocessing(provenance, config):
    """Copy provenance and change only the reduction actually applied.

    Legacy outputs retain unknown detection settings instead of retroactively
    attributing the reducer's entire configuration to the original detection.
    """
    updated = deepcopy(provenance) if provenance else {
        "schema_version": 1, "detection_config": None,
    }
    updated["postprocessing_config"] = {
        "reduce_to_centers": config.reduce_to_centers, "center_score": config.center_score,
    }
    return json_safe(updated)
