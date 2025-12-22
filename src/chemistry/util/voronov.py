"""Utility to load Voronov (1997) collisional ionization coefficients."""

from __future__ import annotations

import functools
from pathlib import Path
from typing import Dict, Tuple

Coeff = Dict[str, float]


def _parse_ion_field(field: str) -> Tuple[str, int]:
    """Return (element, charge) from the 5-char Ion field."""
    parts = field.strip().split()
    if not parts:
        raise ValueError(f"Cannot parse ion field '{field}'")
    element = parts[0]
    charge = 0
    if len(parts) > 1:
        charge_str = parts[1].strip()
        if charge_str.endswith("+"):
            charge_str = charge_str[:-1]
        charge = int(charge_str) if charge_str else 0
    return element, charge


@functools.lru_cache(maxsize=1)
def load_voronov_ci_table(path: str | Path | None = None) -> Dict[Tuple[str, int], Coeff]:
    """Load Voronov CI coefficients from the provided table (ci.dat).

    Returns a mapping {(element, charge_state): {"A": ..., "P": ..., "X": ..., "K": ..., "E0_eV": ...}}
    where charge_state is the ionic charge of the lower stage (e.g., CIII -> CIV uses charge=2).
    """
    if path is None:
        path = Path(__file__).resolve().parent / "ci.dat"
    else:
        path = Path(path)
    coeffs: Dict[Tuple[str, int], Coeff] = {}
    with path.open("r", encoding="ascii") as f:
        for line in f:
            if not line.strip():
                continue
            ion_field = line[0:5]
            element, charge = _parse_ion_field(ion_field)
            try:
                Z = int(line[6:8])
                N = int(line[9:11])
            except ValueError as exc:
                raise ValueError(f"Failed to parse Z/N in line: {line!r}") from exc
            try:
                E0 = float(line[12:19])
                P_param = float(line[20:21])
                A = float(line[22:31])
                X = float(line[32:38])
                K = float(line[39:43])
            except ValueError as exc:
                raise ValueError(f"Failed to parse coefficients in line: {line!r}") from exc

            # Use (Z, N) to catch duplicates; prefer the first occurrence.
            key = (element, charge)
            if key not in coeffs:
                coeffs[key] = {"A": A, "P": P_param, "X": X, "K": K, "E0_eV": E0, "Z": Z, "N": N}
    return coeffs


def get_voronov_coeff(element: str, charge: int) -> Coeff:
    """Fetch coefficients for a given element symbol and ionic charge."""
    table = load_voronov_ci_table()
    key = (element, charge)
    if key not in table:
        raise KeyError(f"No Voronov CI coefficients for {element} charge {charge}+")
    coeff = table[key]
    return {k: coeff[k] for k in ("A", "P", "X", "K", "E0_eV")}
