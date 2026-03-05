import numpy as np

import units as un
from dust_consts import *

import utils as ut

v_turb = 10.0  # km/s

dt = 1e0 # Myr

a = 1e-3  # um
b = 1e-1 # um
N_bins = 10  # dimensionless


T = 1e6  # K
n_H = 1.0  # cm^-3
Z = 1.0  # Z_sun

# Simulation parameters
V_cell_pc = 1 # pc^3 (kept for diagnostics; evolution below uses number density)
D = 0.5 # Dust to metal ratio, dimensionless
Z_solar = 0.0134 # Solar metallicity, dimensionless
Zgas = Z_solar  # Gas metallicity, dimensionless

V_cell_cc = V_cell_pc * (un.pc_cgs**3)

# Initial dust density 
rho_d = D * Zgas * n_H  # amu/cm^3 Gas mass density

dist_type = "linear"
# dist_type = "powerlaw"

dist_index = -3.5

# Size edges
dust_bin_edges = np.logspace(np.log10(a), np.log10(b), num=N_bins + 1)  # um

# Bin edge multiplication factor
bin_mult = dust_bin_edges[1] / dust_bin_edges[0]  # dimensionless

# centers
dust_bin = dust_bin_edges[:-1] * np.sqrt(bin_mult)  # um
mass_bin = (4 * np.pi / 3) * rho_gr_um * dust_bin**3  # g
slope_bin = np.zeros_like(dust_bin) + dist_index # dimensionless

dust_bin_sizes = dust_bin_edges[1:] - dust_bin_edges[:-1]  # um

def num_sum(a1, a2, dist):
    return ut.bin_num_integral(a1, 
                            a2, 
                            dist, 
                            slope_bin, 
                            dust_bin_edges, 
                            dust_bin, 
                            amax=a, 
                            amin=b,
                            )
def mass_sum(a1, a2, dist):
    return ut.bin_mass_integral(a1, 
                            a2, 
                            dist, 
                            slope_bin, 
                            dust_bin_edges, 
                            dust_bin, 
                            amax=a, 
                            amin=b,
                            )

# Number density spectrum at bin centers [grains / (um cm^3)]
a_mid = np.sqrt(a*b)
mass_dust = K_dust * a_mid**4 # g / (grains / um / cm^3)
mass_dust *= ut.integrate_power(slope=dist_index+3, x0=a/a_mid, x1=b/a_mid)
mass_dust /= un.atomic_mass_unit_cgs

# Norm for the underlying number distribution
# that is, norm*(x/a_mid)**dist_index
norm = rho_d/mass_dust

a0 = dust_bin_edges[:-1]
a1 = dust_bin_edges[1:]


# Number distribution at bin centers
num_in_bin = norm * a_mid * np.vectorize(ut.integrate_power)(
                                slope=dist_index, 
                                x0=a0/a_mid, 
                                x1=a1/a_mid, 
                                )
dust_dist_init = num_in_bin/dust_bin_sizes

if dist_type=="powerlaw":
    # Fill the slopes
    log_y = np.log(dust_dist_init)
    log_x = np.log(dust_bin)

    deln = np.roll(log_y, -1) - np.roll(log_y, 1)
    da = np.roll(log_x, -1) - np.roll(log_x, 1)
    slope_bin = deln / da  # grains/um^2/cm^3

    slope_bin[0] = (log_y[1] - log_y[0]) / (log_x[1] - log_x[0])  # grains/um^2/cm^3
    slope_bin[-1] = (log_y[-1] - log_y[-2]) / (log_x[-1] - log_x[-2])  # grains/um^2/cm^3
else: # includes dist_type=="linear"
    # Match per-bin mass exactly for n(a) = n_i + m_i (a - a_i):
    # M_i = K_dust * (n_i * I0_i + m_i * I1_i)
    # I0_i = \int_{aL}^{aR} a^3 da
    # I1_i = \int_{aL}^{aR} (a-a_i) a^3 da
    target_mass_bin = (
        norm
        * K_dust
        * a_mid**4
        * np.vectorize(ut.integrate_power)(
            slope=dist_index + 3,
            x0=a0 / a_mid,
            x1=a1 / a_mid,
        )
    )  # g/cm^3 in each bin

    I0 = 0.25 * (a1**4 - a0**4)
    I1 = 0.2 * (a1**5 - a0**5) - 0.25 * dust_bin * (a1**4 - a0**4)

    slope_bin = (target_mass_bin / K_dust - dust_dist_init * I0) / I1


if __name__ == "__main__":
    # Analytic integrals of n(a) = norm * (a/a_mid)^dist_index
    dust_num_sum = norm * a_mid * ut.integrate_power(
        slope=dist_index, x0=a / a_mid, x1=b / a_mid
    )  # /cc
    dust_mass = (
        norm
        * K_dust
        * a_mid**4
        * ut.integrate_power(slope=dist_index + 3, x0=a / a_mid, x1=b / a_mid)
        / un.atomic_mass_unit_cgs
    )  # amu/cc
    

    n_dist = ut.bin_num_integral(a, 
                                 b, 
                                 dust_dist_init, 
                                 slope_bin, 
                                 dust_bin_edges, 
                                 dust_bin,
                                 amin=a,
                                 amax=b,
                                 type=dist_type,
                                 )
    m_dist = ut.bin_mass_integral(a, 
                                 b, 
                                 dust_dist_init, 
                                 slope_bin, 
                                 dust_bin_edges, 
                                 dust_bin,
                                 amin=a,
                                 amax=b,
                                 type=dist_type,
                                 )/un.atomic_mass_unit_cgs

    info_str = f"Dust density: {m_dist:.2e} amu/cc, "
    info_str += f"# density: {n_dist:.2e}/cc, "
    info_str += f"D = {m_dist/n_H/Z_solar:.2f} Zsol"

    print(f"integrated n_dust = {dust_num_sum:.6e} /cc (Dist: {n_dist:.6e} /cc)")
    print(
        f"integrated rho_dust = {dust_mass:.6e} amu/cc "
        f"(Dist: {m_dist:.6e} amu/cc, target: {rho_d:.6e} amu/cc, "
        f"D={dust_mass/(n_H*Z_solar):.6f})"
    )

    if dist_type=="linear":
        plot_fn = ut.plot_piecewise_linear_distribution
    elif dist_type=="powerlaw":
        plot_fn = ut.plot_piecewise_powerlaw_distribution

    plot_fn(
        dust_bin_edges=dust_bin_edges,
        dust_bin=dust_bin,
        dust_dist=dust_dist_init,
        slope_bin=slope_bin,
        info_str=info_str,
    )
