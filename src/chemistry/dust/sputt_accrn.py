from dust_consts import *

def sputtering() -> float:
    """
    Calculate the grain-size rate due to sputtering.
    Draine
    Parameters:
    a (float): Current grain size in micrometers.
    T (float): Gas temperature in Kelvin.
    n_H (float): Hydrogen number density in cm^-3.
    Z (float): Metallicity relative to solar.
    dt (float): Time step in Myr.

    Returns:
    float: Grain-size rate in micrometers per Myr.
    """
    da_dt = -1.0  # um/Myr
    da_dt *= n_H / 1.0  # cm^-3
    da_dt /= 1 + (2e6 / T) ** 2.5  # K
    da_dt *= Z / 1.0  # Solar metallicity
    da_dt *= grain_surface_porosity ** (-2 / 3)

    return da_dt  # micrometers per Myr

def accretion() -> float:
    """
    Calculate the grain-size rate due to accretion.
    Hirashita+2011

    Exact expression from Hirashita+ 2022, Eqn 22

    Parameters:
    a (float): Current grain size in micrometers.
    T (float): Gas temperature in Kelvin.
    n_H (float): Hydrogen number density in cm^-3.
    Z (float): Metallicity relative to solar.
    dt (float): Time step in Myr.

    Returns:
    float: Grain-size rate in micrometers per Myr.
    """
    da_dt = 0.1 / 537  # micrometers per Myr
    da_dt *= n_H / 1e3  # cm^-3
    da_dt /= (T / 10) ** 0.5  # K
    da_dt *= Z / 1.0  # Solar metallicity
    da_dt *= grain_surface_porosity ** (-2 / 3)

    return da_dt  # micrometers per Myr
