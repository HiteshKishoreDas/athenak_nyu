import numpy as np

import units as un
from dust_consts import *

v_turb = 10.0  # km/s

dt = 1e0 # Myr

dust_bin_edges = np.logspace(np.log10(a), np.log10(b), num=N_bins + 1)  # um
# centers
bin_mult = dust_bin_edges[1] / dust_bin_edges[0]  # dimensionless
dust_bin = dust_bin_edges[:-1] * np.sqrt(bin_mult)  # um
mass_bin = (4 * np.pi / 3) * rho_gr_um * dust_bin**3  # g

dust_bin_sizes = dust_bin_edges[1:] - dust_bin_edges[:-1]  # um

n_dust_cm3 = rho_d / (4/3 * np.pi * rho_gr * (dust_bin * um_cgs)**3)  # grains/cm^3

# Number distribution normalized to n_dust_cm3.
dust_dist_shape = (dust_bin / dust_bin[0]) / dust_bin_sizes  # 1/um
dust_norm = n_dust_cm3 / np.sum(dust_dist_shape * dust_bin_sizes)  # grains/cm^3
dust_dist_init = dust_dist_shape * dust_norm  # grains/um/cm^3

