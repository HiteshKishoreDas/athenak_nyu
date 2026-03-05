# Evolve the grain size instead of the grain density
# Following Hirashita & Yan 2009
# https://doi.org/10.1111/j.1365-2966.2009.14405.x

import numpy as np
import matplotlib.pyplot as plt

import units as un
from dust_consts import *
from ic import *

from rebin import rebin
from sputt_accrn import sputtering, accretion
from coagulation import coagulation
from shattering import shattering
import utils as ut

#* ================ The Run ==========================

da_dt_sputtering = sputtering()  # um/Myr
da_dt_accretion = accretion()  # um/Myr

print(f"Change in grain size due to sputtering: {da_dt_sputtering:.2e} micrometers/Myr")
print(f"Change in grain size due to accretion: {da_dt_accretion:.2e} micrometers/Myr")

net_da_dt = da_dt_accretion + da_dt_sputtering  # um/Myr

sub_dt = dt * 0.01  # Myr
shifted1 = net_da_dt * sub_dt  # um

print(f"Shift in grain size after {sub_dt:.4f} Myr: {shifted1:.5f} micrometers \n\n")

shatt = dust_dist_init.copy()  # grains/um/cm^3
shatt_slopes = slope_bin.copy()

overflow, underflow = 0.0, 0.0  # grains/cm^3

for i in range(5):
    shatt, shatt_slopes, under_add, over_add = rebin(
        dist=shatt,
        slopes=shatt_slopes,
        centers=dust_bin,
        edges=dust_bin_edges,
        shift=shifted1,
    )
    overflow += over_add
    underflow += under_add

    shatt, under_add, over_add = coagulation(shatt, sub_dt)
    overflow += over_add
    underflow += under_add

    shatt, under_add, over_add = shattering(shatt, sub_dt)
    overflow += over_add
    underflow += under_add

dust_dist2 = dust_dist_init.copy()  # grains/um/cm^3
dust_slopes2 = slope_bin.copy()
overflow2, underflow2 = 0.0, 0.0  # grains/cm^3

for i in range(5):
    dust_dist2, shatt_slopes, under_add, over_add = rebin(
        dist=dust_dist2,
        slopes=dust_slopes2,
        centers=dust_bin,
        edges=dust_bin_edges,
        shift=shifted1,
    )
    overflow2 += over_add
    underflow2 += under_add

fig, axs = plt.subplots(1, 1, figsize=(12, 5), sharex=True, sharey=True)
ut.plot_piecewise_powerlaw_distribution(
    dust_bin_edges=dust_bin_edges,
    dust_bin=dust_bin,
    dust_dist=dust_dist_init,
    slope_bin=slope_bin,
    info_str="Original distribution",
    ax=axs,
    show=False,
)
ut.plot_piecewise_powerlaw_distribution(
    dust_bin_edges=dust_bin_edges,
    dust_bin=dust_bin,
    dust_dist=shatt,
    slope_bin=shatt_slopes,
    info_str="Rebinned distribution",
    ax=axs,
    show=False,
)
fig.suptitle(f"Piecewise power law distribution")
fig.tight_layout()
plt.show()
plt.close()

plt.figure()
plt.plot(dust_bin, dust_dist_init, label="Initial")
plt.plot(dust_bin, dust_dist2, label="Just shifted", linestyle="dashed")
plt.plot(dust_bin, shatt, label=f"All processes")
plt.xlabel("Grain size (micrometers)")
plt.ylabel("Number density")
plt.title(f"{underflow} underflow, {overflow:.2f} overflow")
plt.legend()
plt.xscale("log")
plt.yscale("log")
# plt.ylim(1e-13, None)
plt.grid()
plt.show()

plt.figure()
plt.plot(dust_bin, (dust_dist2-shatt)/shatt, label="Diff", linestyle="dashed")
plt.xlabel("Grain size (micrometers)")
plt.ylabel("Number density")
plt.title(f"{underflow} underflow, {overflow:.2f} overflow")
plt.legend()
plt.xscale("log")
plt.yscale("linear")
# plt.ylim(1e-13, None)
plt.grid()
plt.show()

print("Underflow (grains/cm^3):", underflow)
print("Overflow (grains/cm^3):", overflow)
print("Underflow (grains/cm^3):", underflow2)
print("Overflow (grains/cm^3):", overflow2)

mult = 4 * np.pi * rho_gr_um * (dust_bin**3) / 3 # mass per grain [g]
mult /= un.atomic_mass_unit_cgs  # convert mass of grains in amu
# dust_dist * mult gives grain mass density spectrum [amu / (um cm^3)]

plt.figure()
plt.plot(dust_bin, dust_dist_init*mult, label="Initial")
plt.plot(dust_bin, dust_dist2*mult, label="Shifted", linestyle="dashed")
plt.plot(dust_bin, shatt*mult, label=f"After shatt and coag")
plt.xlabel("Grain size (micrometers)")
plt.ylabel("Density")
plt.title(f"{underflow} underflow, {overflow:.2f} overflow")
plt.legend()
plt.xscale("log")
plt.yscale("log")
# plt.ylim(1e-13, None)
plt.grid()
plt.show()

fig, axs = plt.subplots(1, 1, figsize=(12, 5), sharex=True, sharey=True)
ut.plot_piecewise_powerlaw_distribution(
    dust_bin_edges=dust_bin_edges,
    dust_bin=dust_bin,
    dust_dist=dust_dist_init,
    slope_bin=slope_bin,
    info_str="Original distribution",
    ax=axs,
    show=False,
)
ut.plot_piecewise_powerlaw_distribution(
    dust_bin_edges=dust_bin_edges,
    dust_bin=dust_bin,
    dust_dist=dust_dist2,
    slope_bin=dust_slopes2,
    info_str="Rebinned distribution",
    ax=axs,
    show=False,
)
fig.suptitle(f"Piecewise power law distribution")
fig.tight_layout()
plt.show()
plt.close()
