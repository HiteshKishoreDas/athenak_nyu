"""Sanity tests for network loading and ODE system construction (functional API)."""

from __future__ import annotations

from pathlib import Path
from typing import Tuple
import sys

import numpy as np

_THIS_DIR = Path(__file__).resolve().parent
_ROOT = _THIS_DIR.parent  # chemistry/
_SRC_ROOT = _ROOT.parent  # src/
for p in (str(_ROOT), str(_SRC_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

from chemistry.custom_network_loader import load_network
from chemistry.custom_network_loader import visualize_network
from chemistry.build_ode_system import make_ode_system, prepare_initial_concentrations
from chemistry.util.jacobian_utils import analytic_jacobian


def assert_finite(arr, name: str):
    if not np.all(np.isfinite(arr)):
        raise AssertionError(f"{name} contains non-finite entries")


def run_framework_smoke(
    network_dir: Path = Path(__file__).resolve().parent.parent / "custom",
) -> Tuple[bool, str]:
    net = load_network(network_dir)
    if not net.get("species") or not net.get("reactions"):
        return False, f"Network directory missing or unreadable: {network_dir}"

    species, S, rate_func, reaction_data = make_ode_system(net)
    if len(species) != len(net["species"]):
        return False, "Species list length mismatch after ODE system build."
    if S.shape != (len(species), len(net["reactions"])):
        return False, f"Stoichiometric matrix shape mismatch: {S.shape}"

    x0 = prepare_initial_concentrations(species)
    rates, rhs = rate_func(x0, t=0.0, T=100.0)

    if len(rates) != len(net["reactions"]):
        return False, "Rates length mismatch."
    if len(rhs) != len(species):
        return False, "RHS length mismatch."

    assert_finite(rates, "rates")
    assert_finite(rhs, "rhs")

    J = analytic_jacobian(x0, reaction_data, species, S, net, T=100.0)
    if J.shape != (len(species), len(species)):
        return False, f"Jacobian shape mismatch: {J.shape}"
    assert_finite(J, "Jacobian")
    return True, "Framework load/build smoke test passed."


def test_visualize_network_smoke():
    net = load_network(Path(__file__).resolve().parent.parent / "custom")
    dot = visualize_network(net, max_reactions=5, temperature=100.0)
    assert isinstance(dot, str)
    assert "digraph ChemicalNetwork" in dot
    if net["reactions"]:
        first_rxn = net["reactions"][0]
        if first_rxn.get("reactants"):
            assert first_rxn["reactants"][0] in dot


def main():
    ok, msg = run_framework_smoke()
    print(f"{'✅' if ok else '❌'} {msg}")

    # Run visualization smoke in main for quick manual check
    try:
        net = load_network(Path(__file__).resolve().parent.parent / "custom")
        dot = visualize_network(net, max_reactions=5, temperature=100.0)
        has_keyword = "digraph ChemicalNetwork" in dot
        has_species = bool(net["reactions"]) and bool(net["reactions"][0].get("reactants")) and net["reactions"][0]["reactants"][0] in dot
        viz_ok = has_keyword and (not net["reactions"] or has_species)
        print(f"{'✅' if viz_ok else '❌'} Visualization DOT generation")
    except Exception as exc:
        print(f"❌ Visualization DOT generation failed: {exc}")


if __name__ == "__main__":
    main()
