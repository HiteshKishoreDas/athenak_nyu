"""Configuration for the shock finder."""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral, Real


@dataclass
class ShockFinderConfig:
    """Shock criteria shared by uniform, chunked, and AMR execution.

    Widths use local cell-index units. Upstream floors use the input units
    (temperature means P/rho). mach_max caps temperature inversion only.
    AMR callers must explicitly disable the default center reduction."""

    div_v_threshold: float = 0.0
    grad_T_min: float = 1e-30
    require_gradT_gradS_alignment: bool = True
    require_gradT_gradRho_alignment: bool = True
    normal_field: str = "temperature"
    shock_width_cells: int = 1
    require_pressure_jump: bool = True
    require_temperature_jump: bool = True
    require_density_jump: bool = True
    min_mach: float = 1.1
    mach_max: float = 1000.0
    reduce_to_centers: bool = True
    center_score: str = "compression"
    sampling_method: str = "nearest_axis"
    upstream_pressure_floor: float | None = None
    upstream_temperature_floor: float | None = None
    upstream_density_floor: float | None = None

    # Endpoint distance in the same physical units as dx/dy/dz. None preserves
    # local index-space widths. Explicit distance overrides shock_width_cells.
    sampling_distance: float | None = None

    def __post_init__(self) -> None:
        """Validate initial values; mutable configurations are checked again at execution."""
        self.validate()

    def validate(self) -> None:
        """Raise ValueError for invalid values without coercing or mutating fields.

        Numeric criteria must be finite real scalars; only None disables a floor.
        Width is an integral number of local cells, with booleans excluded.
        """
        if self.sampling_distance is not None:
            if (isinstance(self.sampling_distance, bool)
                or not isinstance(self.sampling_distance, Real)
                or not math.isfinite(self.sampling_distance)
                or self.sampling_distance <= 0):
                raise ValueError("sampling_distance must be finite and strictly positive.")
        for name in (
            "div_v_threshold", "grad_T_min", "min_mach", "mach_max",
            "upstream_pressure_floor", "upstream_temperature_floor",
            "upstream_density_floor",
        ):
            value = getattr(self, name)
            if value is None and name.startswith("upstream_"):
                continue
            if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
                raise ValueError(f"{name} must be a finite real scalar.")
        if isinstance(self.shock_width_cells, bool) or not isinstance(self.shock_width_cells, Integral):
            raise ValueError("shock_width_cells must be a positive integer, excluding bool.")

        if self.normal_field not in {"temperature", "pressure"}:
            raise ValueError("normal_field must be 'temperature' or 'pressure'.")
        if self.grad_T_min < 0.0:
            raise ValueError("grad_T_min must be >= 0.")
        if self.min_mach < 1.0:
            raise ValueError("min_mach must be >= 1.")
        if self.shock_width_cells < 1:
            raise ValueError("shock_width_cells must be >= 1.")
        if self.mach_max <= 1.0:
            raise ValueError("mach_max must be > 1.")
        if self.center_score not in {
            "compression",
            "mach_pressure",
            "mach_temperature",
            "combined",
        }:
            raise ValueError(
                "center_score must be one of: compression, mach_pressure, mach_temperature, combined."
            )
        if self.sampling_method not in {"nearest_axis", "trilinear"}:
            raise ValueError(
                "sampling_method must be either 'nearest_axis' or 'trilinear'."
            )
        for name, value in (
            ("upstream_pressure_floor", self.upstream_pressure_floor),
            ("upstream_temperature_floor", self.upstream_temperature_floor),
            ("upstream_density_floor", self.upstream_density_floor),
        ):
            if value is not None and value < 0.0:
                raise ValueError(f"{name} must be >= 0 when provided.")
