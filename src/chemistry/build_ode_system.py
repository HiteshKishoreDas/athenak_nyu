# -*- coding: utf-8 -*-
"""Build ODE system for a chemical network.

This script builds a stoichiometric matrix and a rate-function for the
network loaded by `CustomChemicalNetwork` in `custom_network_loader.py`.

Usage:
    - Import `build_system` and call `make_ode_system(network)`.
    - Use the returned `(species, S, rate_func)` to evaluate RHS for ODE solvers.

The rate function has signature: `rates = rate_func(concs, t, T, av=10.0, cr_zeta=1.3e-17)`
where `concs` is a 1D array of species concentrations ordered as `species`.

A small demo at the bottom shows usage and computes dx/dt for an example
concentration vector.
"""

from __future__ import annotations
from pathlib import Path
import sys
from typing import Dict

import numpy as np

# Make sure the repository `src` parent is on path so `chemistry` can be imported
_THIS_DIR = Path(__file__).resolve().parent
_ROOT = _THIS_DIR.parent

if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

# Detect Sundials (scikits.odes) availability
_HAS_SUNDIALS = False
try:
    from scikits.odes import ode as _scikits_ode  # type: ignore

    _HAS_SUNDIALS = True
except Exception:
    _HAS_SUNDIALS = False

# Import loader helpers
from chemistry.custom_network_loader import calculate_rate_coefficient

# External drivers that are NOT concentration variables (treated as constant drivers)
_EXTERNAL_DRIVERS = {"CR", "Photon", "CRP"}


def make_ode_system(network: Dict):
    """
    Construct species list, stoichiometric matrix, reaction data, and a rate-evaluation function.

    Returns:
        species: list[str] -- ordered list of species names
        S: np.ndarray, shape (n_species, n_reactions) -- stoichiometric matrix
        rate_func: callable(concs, t, T, av=10.0, cr_zeta=1.3e-17) -> (rates, rhs)
            - rates: 1D array length n_reactions with reaction rates (s^-1 or cm^3 s^-1)
            - rhs: 1D array length n_species with d[species]/dt
        reaction_data: list of dicts containing reaction information for Jacobian computation
    """
    # Build ordered species list
    species = sorted(network["species"].keys())
    n_species = len(species)
    n_reactions = len(network["reactions"])

    # Map species -> index
    idx = {s: i for i, s in enumerate(species)}

    # Stoichiometric matrix (species x reactions), numpy only
    S = np.zeros((n_species, n_reactions), dtype=float)

    # For each reaction, record reactant indices and multiplicities
    reaction_data = (
        []
    )  # list of dicts containing reactant_indices, reactant_counts, prod_indices, prod_counts, orig_reaction

    for j, rxn in enumerate(network["reactions"]):
        reactants = rxn.get("reactants", [])
        products = rxn.get("products", [])

        # Count multiplicities (in case of repeated names)
        reac_counts: Dict[str, int] = {}
        prod_counts: Dict[str, int] = {}

        for r in reactants:
            if r in _EXTERNAL_DRIVERS:
                # external driver does not consume a concentration variable
                continue
            reac_counts[r] = reac_counts.get(r, 0) + 1

        for p in products:
            if p in _EXTERNAL_DRIVERS:
                continue
            prod_counts[p] = prod_counts.get(p, 0) + 1

        # Fill stoichiometric matrix columns
        for r_name, count in reac_counts.items():
            if r_name in idx:
                S[idx[r_name], j] -= count

        for p_name, count in prod_counts.items():
            if p_name in idx:
                S[idx[p_name], j] += count

        # Save compact reaction info
        reac_inds = [idx[r] for r in reactants if r in idx]
        # create multiplicity map for vectorized rate law
        reac_mult = [reac_counts[r] for r in reactants if r in idx]

        reaction_data.append(
            {
                "reactants": reactants,
                "products": products,
                "reac_counts": reac_counts,
                "prod_counts": prod_counts,
                "reac_inds": reac_inds,
                "reac_mult": reac_mult,
                "rxn": rxn,
            }
        )

    def rate_func(
        concs: np.ndarray,
        t: float = 0.0,
        T: float = 300.0,
        av: float = 10.0,
        cr_zeta: float = 1.3e-17,
    ):
        """
        Evaluate reaction rates and the RHS (d[species]/dt) for given concentrations.

        Parameters:
            concs: 1D array (n_species,) of concentrations (cm^-3 or normalized units)
            t: time (not used in rate evaluation by default, present for ODE solver signature)
            T: temperature in K
            av: visual extinction (for photorates)
            cr_zeta: cosmic-ray ionization rate

        Returns:
            rates: 1D array (n_reactions,) of reaction rates
            rhs: 1D array (n_species,) of d[species]/dt
        """
        concs = np.asarray(concs, dtype=float)
        rates = np.zeros(n_reactions, dtype=float)

        # Evaluate each reaction rate
        for j, info in enumerate(reaction_data):
            rxn = info["rxn"]
            k = calculate_rate_coefficient(rxn, temperature=T, av=av, cr_zeta=cr_zeta)

            # Determine concentration factor: product [X]^n for each reactant X that's in concentrations
            c_factor = 1.0

            # If there are no reactant species (all reactants are external), treat as zeroth-order: c_factor=1
            for r_name, count in info["reac_counts"].items():
                if r_name in idx:
                    c = concs[idx[r_name]]
                    # negative or zero concentrations will still propagate through; user should ensure non-negative initial conditions
                    c_factor *= c**count

            rates[j] = k * c_factor

        # Compute RHS
        rhs = S.dot(rates)
        return rates, rhs

    return species, S, rate_func, reaction_data


def integrate_with_sundials(
    rate_func,
    y0,
    t0,
    tf,
    T=300.0,
    av=10.0,
    cr_zeta=1.3e-17,
    rtol=1e-6,
    atol=1e-12,
    t_eval=None,
    reaction_data=None,
    species=None,
    S=None,
    network=None,
    use_analytic_jac=False,
):
    """
    Integrate the ODE system using SUNDIALS (via scikits.odes) if available.

    Parameters:
        rate_func: function(concs, t, T, av, cr_zeta) -> (rates, rhs)
        y0: initial concentrations (array-like)
        t0: initial time
        tf: final time
        T, av, cr_zeta: physical params passed to rate_func
        rtol, atol: tolerances
        t_eval: optional array of times to store
        reaction_data, species, S, network: needed if use_analytic_jac=True
        use_analytic_jac: if True, use analytic Jacobian instead of numeric FD

    Returns:
        The solver result object from scikits.odes (if available).
    """
    if not _HAS_SUNDIALS:
        raise RuntimeError(
            "SUNDIALS bindings (scikits.odes) not available. Install via:\n"
            "  pip install scikits.odes\n"
            "See https://github.com/bmcinnes/scikits.odes for build instructions."
        )

    # build RHS for scikits.odes: signature f(t, y, ydot)
    def _f(t, y, ydot):
        # scikits.odes provides y as numpy array; we rely on rate_func to accept it
        _, rhs = rate_func(y, t, T=T, av=av, cr_zeta=cr_zeta)
        # copy into ydot
        for i, val in enumerate(rhs):
            ydot[i] = val
        return 0

    # build numeric Jacobian (finite-difference) for CVODE if requested
    def _jac(t, y, *args):
        """
        Flexible Jacobian callback that supports multiple scikits.odes signatures.

        Common signatures seen in different versions:
          - _jac(t, y, J)                -> J is writable, return 0
          - _jac(t, y, ydot, J)         -> J is writable, return 0
          - _jac(t, y)                  -> return dense Jacobian array

        We treat the last argument as the writable matrix `J` when present,
        otherwise we return a new 2D array (or list-of-lists) containing the
        finite-difference Jacobian.
        """
        # Determine if a writable J was provided (last arg), and ignore any ydot in middle
        J_provided = len(args) >= 1
        J = args[-1] if J_provided else None

        y_arr = np.asarray(y, dtype=float)
        n = y_arr.size
        f0 = np.asarray(
            rate_func(y_arr, t, T=T, av=av, cr_zeta=cr_zeta)[1], dtype=float
        )
        eps = 1e-8
        if not J_provided:
            Jmat = np.zeros((n, n), dtype=float)
        for j in range(n):
            dy = eps * max(1.0, abs(y_arr[j]))
            y_pert = y_arr.copy()
            y_pert[j] += dy
            f1 = np.asarray(
                rate_func(y_pert, t, T=T, av=av, cr_zeta=cr_zeta)[1], dtype=float
            )
            col = (f1 - f0) / dy
            if J_provided:
                for i in range(n):
                    J[i, j] = col[i]
            else:
                Jmat[:, j] = col
        return 0 if J_provided else Jmat

    # Choose Jacobian: analytic if requested and available, else numeric FD
    jac_fn = _jac
    jac_type = "numeric Jacobian"

    if (
        use_analytic_jac
        and reaction_data is not None
        and species is not None
        and S is not None
        and network is not None
    ):
        from chemistry.util.jacobian_utils import analytic_jacobian

        def _jac_analytic(t, y, *args):
            J_provided = len(args) >= 1
            J = args[-1] if J_provided else None
            Jmat = analytic_jacobian(
                y, reaction_data, species, S, network, T=T, av=av, cr_zeta=cr_zeta
            )
            if J_provided:
                for i in range(Jmat.shape[0]):
                    for j in range(Jmat.shape[1]):
                        J[i, j] = Jmat[i, j]
                return 0
            return Jmat

        jac_fn = _jac_analytic
        jac_type = "analytic Jacobian"

    # create solver instance; try to pass jacobian if supported
    jac_passed = False
    try:
        solver = _scikits_ode("cvode", _f, jacfn=jac_fn, atol=atol, rtol=rtol)
        jac_passed = True
    except Exception:
        solver = _scikits_ode("cvode", _f, atol=atol, rtol=rtol)

    # Inform whether we were able to pass the Jacobian callback
    try:
        from sys import stdout

        if jac_passed:
            print(f"[CVODE] {jac_type} callback passed to solver")
        else:
            print(
                "[CVODE] numeric Jacobian not passed (callback not accepted by this build)"
            )
    except Exception:
        pass

    # prepare times
    if t_eval is None:
        # generate 100 output points
        if np is not None:
            t_eval = np.linspace(t0, tf, 100)
        else:
            # simple python linspace
            n = 100
            dt = (tf - t0) / (n - 1)
            t_eval = [t0 + i * dt for i in range(n)]

    # solve: call solver.solve using an appropriate signature
    import inspect

    sig = None
    try:
        sig = inspect.signature(solver.solve)
        params = list(sig.parameters.values())
        nparams = len(params)
    except Exception:
        nparams = None

    # Try to call with the most common forms based on parameter count
    if nparams == 2:
        # expected: (tspan, y0)
        sol = solver.solve(t_eval, y0)
    elif nparams == 3:
        # expected: (t0, y0, tf)
        sol = solver.solve(t0, y0, tf)
    else:
        # Last-resort try several variants
        try:
            sol = solver.solve(t_eval, y0)
        except Exception:
            try:
                sol = solver.solve((t0, tf), y0)
            except Exception:
                sol = solver.solve(t0, y0, tf)

    return sol


def prepare_initial_concentrations(species):
    """Prepare a default initial concentration vector for `species`.

    Returns an array-like matching earlier demo behavior.
    """
    if np is not None:
        x0 = np.zeros(len(species), dtype=float)
    else:
        x0 = [0.0] * len(species)

    for i, s in enumerate(species):
        if s == "e-":
            val = 1e-7
        else:
            val = 1e-8
        if np is not None:
            x0[i] = val
        else:
            x0[i] = val
    return x0


def run_demo(
    rate_func,
    x0,
    species,
    network,
    t0=0.0,
    tf=1e3,
    T=100.0,
    av=10.0,
    count=10,
    S=None,
    reaction_data=None,
    cr_zeta=1.3e-17,
):
    """Run a short demo: sample evaluation and integration, printing results.

    This encapsulates the previous long __main__ demo block so the module's
    entrypoint remains concise.
    """
    # Evaluate rates and RHS at T
    rates, rhs = rate_func(x0, t=0.0, T=T, av=av)

    print("\nSample evaluation at T=100 K")
    print("First 10 reaction rates:")
    for j, r in enumerate(rates[:10]):
        print(f"  r{j+1:3d} = {r:.3e}")

    print("\nFirst 10 entries of dx/dt:")
    for i, val in enumerate(rhs[:10]):
        print(f"  d({species[i]})/dt = {val:.3e}")

    # Try SUNDIALS first, then fall back to SciPy
    if _HAS_SUNDIALS:
        try:
            print("\nRunning integration with SUNDIALS (CVODE)...")
            y0_use = x0
            if np is not None:
                y0_use = np.asarray(x0, dtype=float)

            sol = integrate_with_sundials(
                rate_func,
                y0_use,
                t0,
                tf,
                T=T,
                reaction_data=reaction_data,
                species=species,
                S=S,
                network=network,
                use_analytic_jac=True,
            )
            print("Integration finished with SUNDIALS")

            final = None
            try:
                final = sol.y[-1]
            except Exception:
                try:
                    final = sol.values.y[-1]
                except Exception:
                    final = None

            if final is not None:
                # enforce atom conservation and clamp
                from chemistry.util.cli_utils import (
                    print_clamped_table,
                    enforce_atom_conservation,
                )

                final = enforce_atom_conservation(final, y0_use, species, network)
                if np is not None:
                    final_clamped = np.maximum(final, 0.0)
                else:
                    final_clamped = [max(x, 0.0) for x in final]
                print_clamped_table(species, final, final_clamped, count=count)
            else:
                print("Solver returned result object; inspect `sol` for details.")
        except Exception as e:
            print(f"SUNDIALS integration failed: {e}")
            print("Falling back to SciPy if available...")
            # fall through to SciPy block

    # SciPy fallback or when SUNDIALS not available
    try:
        from scipy.integrate import solve_ivp
        from chemistry.util.jacobian_utils import analytic_jacobian

        print("Running integration with SciPy (solve_ivp)...")

        def rhs_scipy(t, y):
            return rate_func(y, t, T=T, av=av)[1]

        def jac_scipy(t, y):
            return analytic_jacobian(
                y, reaction_data, species, S, network, T=T, av=av, cr_zeta=cr_zeta
            )

        sol = solve_ivp(
            rhs_scipy, (t0, tf), x0, method="BDF", rtol=1e-6, atol=1e-12, jac=jac_scipy
        )
        final = sol.y[:, -1] if sol.y.ndim == 2 else sol.y[-1]
        from chemistry.util.cli_utils import (
            print_clamped_table,
            enforce_atom_conservation,
        )

        final = enforce_atom_conservation(final, x0, species, network)
        if np is not None:
            final_clamped = np.maximum(final, 0.0)
        else:
            final_clamped = [max(x, 0.0) for x in final]
        print_clamped_table(species, final, final_clamped, count=count)
    except Exception as e2:
        print(f"SciPy integration failed or not available: {e2}")


# ---------------------- Demo / CLI portion ----------------------
if __name__ == "__main__":

    # import helpers from cli_utils to keep __main__ concise
    from chemistry.util.cli_utils import print_clamped_table, enforce_atom_conservation

    # Load network from custom dir relative to this file
    curdir = Path(__file__).resolve().parent
    custom_dir = curdir / "custom"

    from chemistry.custom_network_loader import load_network

    network = load_network(custom_dir)
    if not network.get("species") or not network.get("reactions"):
        raise SystemExit("Failed to load network from custom directory")

    species, S, rate_func, reaction_data = make_ode_system(network)

    print("Built ODE system:")
    print(f"  species count = {len(species)}")
    print(f"  reaction count = {len(network['reactions'])}")

    # Prepare default initial concentrations and run the demo (sample + integrate)
    x0 = prepare_initial_concentrations(species)
    run_demo(
        rate_func,
        x0,
        species,
        network,
        t0=0.0,
        tf=1e3,
        T=100.0,
        av=10.0,
        count=10,
        S=S,
        reaction_data=reaction_data,
    )
