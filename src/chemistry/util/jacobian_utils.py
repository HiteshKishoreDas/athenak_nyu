"""Analytic Jacobian computation for ODE systems."""

from __future__ import annotations
from typing import Sequence

import numpy as np
from chemistry.custom_network_loader import calculate_rate_coefficient


def analytic_jacobian(
    concs, reaction_data, species, S, network, T=300.0, av=10.0, cr_zeta=1.3e-17
):
    """Compute analytic Jacobian J (n_species x n_species) from reaction stoichiometry.

    For mass-action kinetics: d(RHS_i)/dc_j = sum_k S[i,k] * d(rate_k)/dc_j
    where d(rate_k)/dc_j = (rate_k / c_j) * n_j if j is a reactant in reaction k, else 0.

    Parameters:
        concs: 1D array-like of concentrations (length n_species)
        reaction_data: list from make_ode_system with keys 'reac_inds','reac_mult','rxn'
        species: list of species names (index ordering)
        S: stoichiometric matrix (n_species x n_reactions)
        network: CustomChemicalNetwork instance
        T, av, cr_zeta: physical parameters for rate coefficient calculation

    Returns:
        Dense numpy array J (n_species x n_species), or pure Python list-of-lists if numpy unavailable.
    """
    concs = np.asarray(concs, dtype=float)
    n_species = len(species)
    n_reactions = len(reaction_data)
    rates = np.zeros(n_reactions, dtype=float)
    for k, info in enumerate(reaction_data):
        rxn = info["rxn"]
        kcoef = calculate_rate_coefficient(rxn, temperature=T, av=av, cr_zeta=cr_zeta)
        prod = 1.0
        for idx, mult in zip(info.get("reac_inds", []), info.get("reac_mult", [])):
            prod *= concs[idx] ** mult
        rates[k] = kcoef * prod
    J = np.zeros((n_species, n_species), dtype=float)
    for k, info in enumerate(reaction_data):
        reac_inds = info.get("reac_inds", [])
        reac_mult = info.get("reac_mult", [])
        if len(reac_inds) == 0:
            continue
        rxn = info["rxn"]
        kcoef = calculate_rate_coefficient(rxn, temperature=T, av=av, cr_zeta=cr_zeta)
        for idx_j, n_j in zip(reac_inds, reac_mult):
            if concs[idx_j] > 1e-30:
                dr_dcj = rates[k] * (n_j / concs[idx_j])
            else:
                pj = (concs[idx_j] ** (n_j - 1)) if n_j > 0 else 0.0
                prod_others = 1.0
                for idx_m, n_m in zip(reac_inds, reac_mult):
                    if idx_m == idx_j:
                        continue
                    prod_others *= concs[idx_m] ** n_m
                dr_dcj = kcoef * n_j * pj * prod_others
            J[:, idx_j] += S[:, k] * dr_dcj
    return J
