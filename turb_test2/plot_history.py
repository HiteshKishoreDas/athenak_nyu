import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mt
from matplotlib.colors import ListedColormap
import gc
import os

import sys
from pathlib import Path

repo_root = Path.cwd().resolve().parents[0]  # adjust if running in notebook
sys.path.append(str(repo_root / "vis" / "python"))

import cmasher as cm
import bin_convert_new as bc

import own_package.utils.units as ut
import own_package.athena.history as hst

import own_package

theme = "bright"
# theme = "dark"

package_path = os.path.dirname(own_package.__file__)
print(f"own_package path: {package_path}")

style_lib = f"{package_path}/plot/style_lib/"

if theme == "dark":
    pallette = style_lib + "dark_pallette.mplstyle"
elif theme == "bright":
    pallette = style_lib + "bright_pallette.mplstyle"

plot_style = style_lib + "plot_style.mplstyle"
text_style = style_lib + "text.mplstyle"

plt.style.use([pallette, plot_style, text_style])

line_border_color = mt.rcParams["lines.color"]
fig_face_color = mt.rcParams["figure.facecolor"]

gamma_minus_1 = 5.0 / 3.0 - 1.0

Z_solar = 0.0134

# * =========================

hst_read = hst.hst_data(
    ["Turb.hydro.hst", "Turb.user.hst"],
    ncells=[256] * 3,
    box_size=[1.0] * 3,
)

# * =========================
# * Cold gas fraction
# * =========================

plt.figure()
mcut = 1 / 3
cold_frac = hst_read.dict["Mcold"] / hst_read.mass_tot
plt.plot(hst_read.time, cold_frac)
t_stop = hst_read.time[np.argmin(np.abs(cold_frac - mcut))]
plt.axhline(mcut, ls=":")
plt.axvline(t_stop, ls=":")
# plt.yscale("log")
plt.ylabel(f"$m_\mathrm{{cold}}/m_\mathrm{{total}}$")
plt.xlabel(f"$t$ (Myr)")

# * =========================
# * Average temperature
# * =========================

plt.figure()
plt.plot(hst_read.time, hst_read.T_avg)
plt.axvline(t_stop, ls=":")

plt.yscale("log")
plt.ylabel(f"$T_\mathrm{{avg}}$")
plt.xlabel(f"$t$ (Myr)")


# * =========================
# * Scalars: Dust, metal, tracer
# * =========================

nscalars = 8
Ddata = []
for i in range(nscalars):
    Ddata.append(hst_read.dict["D" + str(i)])

Ddata = np.array(Ddata)
Dplot = Ddata / Ddata[:, 0][:, np.newaxis]
Dplot = np.log10(Dplot)

# Sample a set of line colors from the same CMasher colormap used for images
line_colors = cm.take_cmap_colors(
    cm.rainforest, nscalars, cmap_range=(0.05, 0.95), return_fmt="hex"
)

# * Image plot species vs time
plt.figure()
plt.imshow(
    Dplot,
    aspect="auto",
    origin="lower",
    vmin=-2.0,
    vmax=2.0,
    cmap=cm.redshift_r,
)
plt.ylabel("Dust Scalar Index")
# plt.xlabel("$t$ (Myr)")
plt.colorbar(label="$\log_{10}(D_i/D_{i,0})$")

# * Individual species line plot
fig_lines, ax_lines = plt.subplots()
for k in range(nscalars):
    ax_lines.plot(
        hst_read.time,
        Ddata[k, :] / hst_read.mass_tot / Z_solar,
        label=f"$D_{k}/Z_\odot$",
        color=line_colors[k],
    )
ax_lines.axvline(t_stop, ls=":")

ax_lines.set_yscale("log")
ax_lines.set_ylabel(f"$D_\mathrm{{tot}}/Z_\odot$")
ax_lines.set_xlabel(f"$t$ (Myr)")

# Map scalar index -> line color with a compact colorbar
scalar_cmap = ListedColormap(line_colors)
norm = plt.Normalize(vmin=-0.5, vmax=nscalars - 0.5)
cbar = fig_lines.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=scalar_cmap),
    ax=ax_lines,
    pad=0.02,
)
cbar.set_label("Dust scalar index $i$")
cbar.set_ticks(range(nscalars))
cbar.set_ticklabels([str(i) for i in range(nscalars)])

# * Metal evolution line plot
plt.figure()
plt.plot(hst_read.time, hst_read.dict["Zgas"] / hst_read.mass_tot / Z_solar)
plt.axvline(t_stop, ls=":")

plt.yscale("log")
plt.ylabel(f"$Z_\mathrm{{avg}}/Z_\odot$")
plt.xlabel(f"$t$ (Myr)")

# dust_a = 0.001 # Min dust size
# dust_b = 10 # Max dust size

# d_r = np.logspace(np.log10(dust_a), np.log10(dust_b), num=nscalars)
