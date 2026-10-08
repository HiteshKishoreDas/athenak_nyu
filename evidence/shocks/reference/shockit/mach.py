"""Rankine-Hugoniot jump relations and inversion utilities."""

from __future__ import annotations

import math


def pressure_jump_from_mach(mach: float, gamma: float) -> float:
    """Return the ideal-gas pressure jump for a normal shock."""

    mach_squared = mach * mach
    return (2.0 * gamma * mach_squared - (gamma - 1.0)) / (gamma + 1.0)


def density_jump_from_mach(mach: float, gamma: float) -> float:
    """Return the ideal-gas density jump for a normal shock."""

    mach_squared = mach * mach
    return ((gamma + 1.0) * mach_squared) / ((gamma - 1.0) * mach_squared + 2.0)


def temperature_jump_from_mach(mach: float, gamma: float) -> float:
    """Return the ideal-gas temperature jump for a normal shock."""

    return pressure_jump_from_mach(mach, gamma) / density_jump_from_mach(mach, gamma)


def mach_from_pressure_jump(pressure_jump: float, gamma: float) -> float:
    """Invert a pressure jump to upstream Mach number."""

    if pressure_jump <= 1.0:
        return math.nan
    mach_squared = ((pressure_jump * (gamma + 1.0)) + (gamma - 1.0)) / (2.0 * gamma)
    return math.sqrt(mach_squared)


def mach_from_temperature_jump(
    temperature_jump: float,
    gamma: float,
    mach_max: float = 1000.0,
) -> float:
    """Invert a temperature jump to upstream Mach number analytically."""

    if temperature_jump <= 1.0:
        return math.nan

    # For R=T_down/T_up > 1 and gamma>1, x=M^2 satisfies
    # 2*gamma*(gamma-1)*x^2 + [-(gamma-1)^2 + 4*gamma
    # - R*(gamma+1)^2]*x - 2*(gamma-1) = 0.
    # Since a>0 and c<0, exactly one root is positive. For physical R>1,
    # b<0, so -b + sqrt(discriminant) has no subtractive cancellation.
    # Reject subsonic roots and roots beyond the configured temperature cap.
    gamma_plus_one = gamma + 1.0
    gamma_minus_one = gamma - 1.0
    coefficient_a = 2.0 * gamma * gamma_minus_one
    coefficient_b = (
        -(gamma * gamma)
        + (6.0 * gamma)
        - 1.0
        - (temperature_jump * gamma_plus_one * gamma_plus_one)
    )
    coefficient_c = -2.0 * gamma_minus_one
    # sqrt(b^2-4ac) = hypot(b, 2*sqrt(a)*sqrt(-c)); the squared
    # discriminant can overflow even when the root and inferred M are finite.
    root_discriminant = math.hypot(
        coefficient_b, 2.0 * math.sqrt(coefficient_a) * math.sqrt(-coefficient_c)
    )
    mach_squared = (-0.5 * coefficient_b + 0.5 * root_discriminant) / coefficient_a
    if mach_squared <= 1.0:
        return math.nan

    mach = math.sqrt(mach_squared)
    # Compare the observed ratio against the forward cap before allowing a
    # root rounded just above mach_max. Monotonic R(M) makes this an exact
    # input-domain cutoff: the next representable ratio above it is rejected.
    # No general tolerance is added to the physical cap.
    cap_ratio = temperature_jump_from_mach(mach_max, gamma)
    if temperature_jump > cap_ratio:
        return math.nan
    if mach > mach_max:
        if temperature_jump <= cap_ratio:
            return mach_max
        return math.nan
    return mach
