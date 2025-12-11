"""CLI helper utilities for chemistry scripts.

Provides ANSI-aware table printing and atom-conservation enforcement used
by `build_ode_system.py` to keep the `__main__` section concise.
"""

from __future__ import annotations
from typing import List, Sequence

try:
    import numpy as np
except Exception:
    np = None

import re


ansi_re = re.compile(r"\x1b\[[0-9;]*m")


def strip_ansi(s: str) -> str:
    return ansi_re.sub("", s)


def pad_visible(s: str, width: int, align: str = "left") -> str:
    vis = strip_ansi(s)
    pad = max(0, width - len(vis))
    if align == "left":
        return s + " " * pad
    else:
        return " " * pad + s


def color(text: str, code: str) -> str:
    return f"\033[{code}m{text}\033[0m"


def print_clamped_table(
    species: Sequence[str],
    final: Sequence[float],
    final_clamped: Sequence[float],
    count: int = 10,
    title: str = "Final concentrations (first {count}):",
):
    """Print an ANSI-aware table showing species and clamped values.

    When a value was clamped to zero, the original value is shown inline as
    `(was ...)` and the clamped value is colored red.
    """
    name_w = 22
    clamp_w = 36

    print("\n" + title.format(count=count))
    sep = "+" + "-" * (name_w + 2) + "+" + "-" * (clamp_w + 2) + "+"
    print(sep)
    hdr = f"| {'Species':<{name_w}} | {'Clamped (cm^-3) (was ...)':^{clamp_w}} |"
    print(hdr)
    print(sep)
    for i in range(min(count, len(species))):
        name = species[i]
        raw = float(final[i])
        clamped = float(final_clamped[i])
        clamp_str = f"{clamped:.3e}"
        if clamped == 0.0 and raw < 0.0:
            clamp_disp = color(clamp_str, "91") + f" (was {raw:.3e})"
        else:
            clamp_disp = clamp_str
        print(
            "| "
            + pad_visible(f"{name}", name_w, "left")
            + " | "
            + pad_visible(clamp_disp, clamp_w, "right")
            + " |"
        )
    print(sep)


def enforce_atom_conservation(
    final, y0_use, species: Sequence[str], network
) -> Sequence[float]:
    """Project solution onto atom-conserving subspace if necessary.

    This mirrors the logic formerly embedded in `build_ode_system.py`'s
    `__main__` block. It returns either the unchanged `final` or a corrected
    copy to restore elemental conservation within tolerances.
    """
    element_order = getattr(network, "element_order", None)
    if element_order is None and isinstance(network, dict):
        element_order = network.get("element_order")
    if element_order is None:
        return final

    n_elements = len(element_order)
    species_elements = []
    for s in species:
        spec = network["species"][s] if isinstance(network, dict) else network.species[s]
        row = [spec["elements"].get(e, 0) for e in element_order]
        species_elements.append(row)

    if np is not None:
        sp_el = np.array(species_elements, dtype=int)
        initial_atoms = sp_el.T.dot(y0_use)
        final_atoms = sp_el.T.dot(final)
        deltas = final_atoms - initial_atoms
        for i, e in enumerate(element_order):
            print(
                f"  {e}: initial={initial_atoms[i]:.4e}, final={final_atoms[i]:.4e}, delta={deltas[i]:.2e}"
            )
        threshold = 1e-12 * np.maximum(np.abs(initial_atoms), 1.0)
        if np.any(np.abs(deltas) > threshold):
            print("[WARNING] Atom conservation violated, projecting solution...")
            final_corr = final.copy()
            for i, e in enumerate(element_order):
                if abs(deltas[i]) > threshold[i]:
                    for j, row in enumerate(sp_el):
                        if row[i] > 0:
                            atoms_per_mol = row[i]
                            corr = -deltas[i] / atoms_per_mol
                            final_corr[j] += corr
                            break
            final_atoms_corr = sp_el.T.dot(final_corr)
            print("After projection:")
            for i, e in enumerate(element_order):
                delta_corr = final_atoms_corr[i] - initial_atoms[i]
                print(
                    f"  {e}: initial={initial_atoms[i]:.4e}, final={final_atoms_corr[i]:.4e}, delta={delta_corr:.2e}"
                )
            return final_corr
        return final
    else:
        initial_atoms = [
            sum(row[i] * y0_use[j] for j, row in enumerate(species_elements))
            for i in range(n_elements)
        ]
        final_atoms = [
            sum(row[i] * final[j] for j, row in enumerate(species_elements))
            for i in range(n_elements)
        ]
        deltas = [final_atoms[i] - initial_atoms[i] for i in range(n_elements)]
        for i, e in enumerate(element_order):
            print(
                f"  {e}: initial={initial_atoms[i]:.4e}, final={final_atoms[i]:.4e}, delta={deltas[i]:.2e}"
            )
        threshold = [1e-12 * max(abs(initial_atoms[i]), 1.0) for i in range(n_elements)]
        if any(abs(deltas[i]) > threshold[i] for i in range(n_elements)):
            print("[WARNING] Atom conservation violated, projecting solution...")
            final_corr = list(final)
            for i, e in enumerate(element_order):
                if abs(deltas[i]) > threshold[i]:
                    for j, row in enumerate(species_elements):
                        if row[i] > 0:
                            atoms_per_mol = row[i]
                            corr = -deltas[i] / atoms_per_mol
                            final_corr[j] += corr
                            break
            final_atoms_corr = [
                sum(row[i] * final_corr[j] for j, row in enumerate(species_elements))
                for i in range(n_elements)
            ]
            print("After projection:")
            for i, e in enumerate(element_order):
                delta_corr = final_atoms_corr[i] - initial_atoms[i]
                print(
                    f"  {e}: initial={initial_atoms[i]:.4e}, final={final_atoms_corr[i]:.4e}, delta={delta_corr:.2e}"
                )
            return final_corr
        return final
