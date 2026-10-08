"""Data model definitions for uniform fluid cubes."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .storage import (
    FieldSource,
    detect_storage_mode,
    validate_chunked_fields,
    validate_in_memory_fields,
)


@dataclass
class FluidCube:
    """Uniform xyz arrays or five compatible lazy chunk references.

    Density and pressure are positive, velocities are Cartesian components,
    spacing uses caller length units, and gamma is dimensionless and >1.
    Array inputs become finite float64 arrays (existing float64 buffers may be
    shared); construction does not mutate their values. Chunk construction
    checks metadata while payload validity is checked lazily when read.
    Temperature is the proxy P/rho, with no inferred physical unit conversion.
    """

    rho: FieldSource
    pressure: FieldSource
    vx: FieldSource
    vy: FieldSource
    vz: FieldSource
    dx: float = 1.0
    dy: float = 1.0
    dz: float = 1.0
    gamma: float = 5.0 / 3.0
    metadata: dict[str, Any] = field(default_factory=dict)
    storage_mode: str = field(init=False)

    def __post_init__(self) -> None:
        """Validate primitive fields and grid metadata for a uniform cube."""

        field_values = {
            "rho": self.rho,
            "pressure": self.pressure,
            "vx": self.vx,
            "vy": self.vy,
            "vz": self.vz,
        }
        storage_mode = detect_storage_mode(field_values)

        if storage_mode == "in_memory":
            for name, value in validate_in_memory_fields(field_values).items():
                setattr(self, name, value)
        elif storage_mode == "chunked":
            for name, value in validate_chunked_fields(field_values).items():
                setattr(self, name, value)
        else:  # pragma: no cover - defensive guard for future modes
            raise ValueError(f"Unsupported FluidCube storage mode: {storage_mode}")

        self._validate_grid_metadata()
        self.storage_mode = storage_mode

    @property
    def is_in_memory(self) -> bool:
        """Whether primitive fields are arrays rather than chunk directory references."""
        return self.storage_mode == "in_memory"

    @property
    def is_chunked(self) -> bool:
        """Whether all primitive fields use compatible chunked storage."""
        return self.storage_mode == "chunked"

    @property
    def shape(self) -> tuple[int, int, int]:
        """Cell counts in x, y, z order, independent of storage mode."""
        return self.rho.shape

    def _validate_grid_metadata(self) -> None:
        """Require finite positive spacings and a constant ideal-gas gamma greater than one."""
        if (
            not np.isfinite(self.dx)
            or not np.isfinite(self.dy)
            or not np.isfinite(self.dz)
        ):
            raise ValueError("dx, dy, dz must be finite.")
        if not np.isfinite(self.gamma):
            raise ValueError("gamma must be finite.")
        if self.dx <= 0.0 or self.dy <= 0.0 or self.dz <= 0.0:
            raise ValueError("dx, dy, dz must be strictly positive.")
        if self.gamma <= 1.0:
            raise ValueError("gamma must be greater than 1.")
