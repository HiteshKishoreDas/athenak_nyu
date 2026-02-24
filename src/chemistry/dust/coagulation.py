import numpy as np

import units as un
from dust_consts import *
from ic import *

def cross_sec() -> np.ndarray:
    """
    Calculate the cross-sectional area for collisions between two grains.

    Returns:
    np.ndarray: Pairwise cross-sectional area matrix in kpc^2.
    """

    A, B = np.meshgrid(dust_bin, dust_bin)
    A_kpc = A * um_to_kpc
    B_kpc = B * um_to_kpc

    return np.pi * (A_kpc + B_kpc) ** 2

def v_coag_fn():

    F_stick = 10  # dimensionless

    A, B = np.meshgrid(dust_bin, dust_bin)
    A_cgs = A * um_cgs
    B_cgs = B * um_cgs

    vc = 2.14  # cgs prefactor
    vc *= F_stick
    vc *= np.sqrt((A_cgs**3 + B_cgs**3) / (A_cgs + B_cgs) ** 3)
    vc *= gamma ** (5 / 6)
    vc /= E ** (1 / 3)
    vc /= (A_cgs * B_cgs / (A_cgs + B_cgs)) ** (5 / 6)
    vc /= rho_gr**0.5

    return vc / un.km_s_cgs  # km/s

def maxwell_head_mean(v: float):
    
    v0 = v_turb * np.sqrt(2 / 3)  # km/s

    vmean = np.sqrt(8 / np.pi) * v0  # km/s
    vmean *= (1 + 0.5 * (v / v0) ** 2)  # km/s
    vmean *= np.exp(-0.5 * (v / v0) ** 2)  # km/s

    vmean = np.sqrt(8 / np.pi) * v0 - vmean  # km/s

    return vmean

def coagulation(dist, dt):

    underflow, overflow = 0.0, 0.0  # grains/cm^3

    dist_mat_coag = np.zeros((N_bins, N_bins), dtype=float)  # unordered pair interactions

    bin_sizes = dust_bin_edges[1:] - dust_bin_edges[:-1]  # um
    num_in_bin = dist * bin_sizes * kpc3_per_cm3  # grains/kpc^3
    alp = cross_sec()

    fA, fB = np.meshgrid(num_in_bin, num_in_bin)
    pair_mat = np.triu(fA * fB)
    pair_mat[np.diag_indices(N_bins)] *= 0.5
    dist_mat_coag += pair_mat

    v_coag_kpc_myr = maxwell_head_mean(v_coag_fn()) * km_s_to_kpc_myr  # kpc/Myr

    dist_mat_coag *= alp * v_coag_kpc_myr * dt

    # Conservative limiter: keep removals below available source grains.
    remove_try = (
        np.sum(dist_mat_coag, axis=1)
        + np.sum(dist_mat_coag, axis=0)
    )
    valid = remove_try > 0
    if np.any(valid):
        frac = np.ones_like(remove_try)
        frac[valid] = num_in_bin[valid] / remove_try[valid]
        global_scale = min(1.0, 0.95 * np.min(frac[valid]))
        dist_mat_coag *= global_scale

    # Multiply num_in_bin by dist_mat_shatt and dist_mat_coag to get
    # the exact number of interaction 
    
    #* ======= Coagulation =======

    # Remove coagulated grains
    num_in_bin -= np.sum(dist_mat_coag, axis=1) + np.sum(dist_mat_coag, axis=0)

    # Add coagulated grains
    for i in range(N_bins):
        for j in range(N_bins):
            a_coag = (dust_bin[i]**3 + dust_bin[j]**3)**(1/3)
            k = np.searchsorted(dust_bin_edges, a_coag) - 1
            if k < N_bins:
                num_in_bin[k] += dist_mat_coag[i, j]
            else:
                num_in_bin[i] += dist_mat_coag[i, j]
                num_in_bin[j] += dist_mat_coag[i, j]


    num_in_bin_cm3 = num_in_bin / kpc3_per_cm3  # grains/cm^3
    return num_in_bin_cm3 / bin_sizes, underflow / kpc3_per_cm3, overflow / kpc3_per_cm3
