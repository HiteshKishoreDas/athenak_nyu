# Evolve the grain size instead of the grain density

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell

grain_surface_porosity = 1.0  # Fraction of grain surface that is porous


def sputtering(a: float, T: float, n_H: float, Z: float) -> float:
    """
    Calculate the change in grain size due to sputtering.
    Draine
    Parameters:
    a (float): Current grain size in micrometers.
    T (float): Gas temperature in Kelvin.
    n_H (float): Hydrogen number density in cm^-3.
    Z (float): Metallicity relative to solar.
    dt (float): Time step in years.

    Returns:
    float: Change in grain size in micrometers.
    """
    # if T < 2e6:
    #     da_dt = -1e-6 * n_H  # micrometers per Myr
    # else:
    da_dt = -1.0  # Myr
    da_dt *= n_H / 1.0  # cm^-3
    da_dt /= 1 + (1e6 / T) ** 3  # K
    da_dt *= Z / 1.0  # Solar metallicity
    da_dt *= grain_surface_porosity ** (-2 / 3)

    return da_dt  # micrometers per Myr


def accretion(a: float, T: float, n_H: float, Z: float) -> float:
    """
    Calculate the change in grain size due to accretion.
    Hirashita+2022

    Parameters:
    a (float): Current grain size in micrometers.
    T (float): Gas temperature in Kelvin.
    n_H (float): Hydrogen number density in cm^-3.
    Z (float): Metallicity relative to solar.
    dt (float): Time step in years.

    Returns:
    float: Change in grain size in micrometers.
    """
    da_dt = 0.1 / 537  # micrometers per Myr
    da_dt *= n_H / 1e3  # cm^-3
    da_dt /= (T / 10) ** 0.5  # K
    da_dt *= Z / 1.0  # Solar metallicity
    da_dt *= grain_surface_porosity ** (-2 / 3)

    return da_dt  # micrometers per Myr


def cross_sec(a1: float, a2: float) -> float:
    """
    Calculate the cross-sectional area for collisions between two grains.

    Parameters:
    a1 (float): Size of the first grain in micrometers.
    a2 (float): Size of the second grain in micrometers.

    Returns:
    float: Cross-sectional area in square micrometers.
    """
    return np.pi * (a1 + a2) ** 2  # micrometers^2


def rel_vel_cdf(v: float, v_turb: float) -> float:
    """
    Calculate the cumulative distribution function of relative velocities for a Maxwellian distribution.

    Parameters:
    v (float): Relative velocity in km/s.
    v_turb (float): turbulent velocity in km/s.

    Returns:
    float: probability of v_rel < v
    """
    # Convert turbulent velocity to velocity dispersion
    # Distribtion of v_rel/np.sqrt(2) is same as v1

    v_rel = maxwell.cdf(v / np.sqrt(2), scale=v_turb / np.sqrt(3))
    print(v_rel)

    return v_rel


if __name__ == "__main__":
    a = 0.1  # micrometers
    T = 1e6  # K
    n_H = 1.0  # cm^-3
    Z = 1.0  # Solar metallicity

    da_dt_sputtering = sputtering(a, T, n_H, Z)
    da_dt_accretion = accretion(a, T, n_H, Z)

    print(
        f"Change in grain size due to sputtering: {da_dt_sputtering:.2e} micrometers/Myr"
    )
    print(
        f"Change in grain size due to accretion: {da_dt_accretion:.2e} micrometers/Myr"
    )

    v_turb = 10.0  # km/s

    rel_vel_cdf(10.0, v_turb / np.sqrt(3))
