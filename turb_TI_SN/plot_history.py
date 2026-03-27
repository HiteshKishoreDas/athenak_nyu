import numpy as np
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mt
from matplotlib.colors import ListedColormap
import cmasher as cm

import own_package
import own_package.athena.history as hst
from plot_cgm_cooling import (
    TABLES_PATH,
    USER_CFG,
    compute_units,
    cooling_heating,
    load_tables,
)

script_dir = Path(__file__).resolve().parent
repo_root = script_dir.parent
sys.path.append(str(repo_root / "vis" / "python"))
sys.path.append(str(repo_root / "src" / "chemistry"))

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

# Z_solar = 0.0134
Z_solar = 0.02


def set_log_scale_if_valid(ax, axis, values, label):
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values) & (values > 0.0)
    if not np.any(valid):
        print(f"Skipping {axis}-log scale for {label}: no positive finite values.")
        return False

    if axis == "x":
        ax.set_xscale("log", nonpositive="mask")
    elif axis == "y":
        ax.set_yscale("log", nonpositive="mask")
    else:
        raise ValueError(f"Unsupported axis '{axis}'")

    return True


def parse_athinput_value(raw_value):
    value = raw_value.split("#", 1)[0].strip()
    lower_value = value.lower()
    if lower_value == "true":
        return True
    if lower_value == "false":
        return False

    try:
        if any(char in lower_value for char in (".", "e")):
            return float(value)
        return int(value)
    except ValueError:
        return value


def load_athinput_cfg(path):
    cfg = USER_CFG.copy()
    current_block = None
    tracked_keys = {
        "mesh": {
            "x1min",
            "x1max",
            "nx1",
            "x2min",
            "x2max",
            "nx2",
            "x3min",
            "x3max",
            "nx3",
        },
        "hydro": {
            "hrate",
            "hscale_norm",
            "hscale_flag",
            "hscale_height",
            "hscale_radius",
            "hscale_alpha",
            "T_max",
            "T_cutoff",
        },
        "units": {"length_cgs", "mass_cgs", "time_cgs", "mu"},
        "problem": {"rho_0"},
    }

    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue

        if line.startswith("<") and line.endswith(">"):
            current_block = line[1:-1].strip().lower()
            continue

        if current_block not in tracked_keys or "=" not in line:
            continue

        key, raw_value = [piece.strip() for piece in line.split("=", 1)]
        if key not in tracked_keys[current_block]:
            continue

        parsed_value = parse_athinput_value(raw_value)
        if current_block == "problem" and key == "rho_0":
            cfg["rho0"] = parsed_value
        else:
            cfg[key] = parsed_value

    return cfg


def compute_cooling_time(history_data, cfg, tables, units):
    box_volume = (
        (cfg["x1max"] - cfg["x1min"])
        * (cfg["x2max"] - cfg["x2min"])
        * (cfg["x3max"] - cfg["x3min"])
    )
    mass_tot = np.asarray(history_data.mass_tot, dtype=float)
    mass_safe = np.clip(mass_tot, np.finfo(float).tiny, None)
    temp_avg = np.clip(
        np.asarray(history_data.T_avg, dtype=float), np.finfo(float).tiny, None
    )
    rho_avg = np.clip(mass_safe / box_volume, np.finfo(float).tiny, None)
    metal_mass = np.asarray(history_data.dict["Zgas"], dtype=float)
    metallicity = np.clip(metal_mass / mass_safe / Z_solar, 0.0, None)

    _, _, cool_vol, heat_vol = cooling_heating(
        temp_avg,
        rho_avg,
        metallicity,
        tables,
        cfg,
        units,
    )

    eint = temp_avg * rho_avg / (units["temperature_cgs"] * gamma_minus_1)
    net_cooling = np.abs(cool_vol - heat_vol)
    t_cool = np.full_like(temp_avg, np.inf)
    valid = net_cooling > 0.0
    t_cool[valid] = eint[valid] / net_cooling[valid]
    return t_cool


# * =========================

hst_read = hst.hst_data(
    ["Turb.hydro.hst", "Turb.user.hst"],
    ncells=[256] * 3,
    box_size=[1.0] * 3,
)

cooling_cfg = load_athinput_cfg(script_dir / "turb_dust.athinput")
cooling_tables = load_tables(TABLES_PATH)
cooling_units = compute_units(cooling_cfg)

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
plt.ylabel("$m_\mathrm{cold}/m_\mathrm{total}$")
plt.xlabel("$t$ (Myr)")

# * =========================
# * Cold gas ratio w/ initial
# * =========================

plt.figure()
cold_frac = hst_read.dict["Mcold"] / hst_read.dict["Mcold"][0]
plt.plot(hst_read.time, cold_frac)
t_stop = hst_read.time[np.argmin(np.abs(cold_frac - mcut))]
plt.axvline(t_stop, ls=":")
set_log_scale_if_valid(plt.gca(), "y", cold_frac, "cold gas ratio")
plt.ylabel("$m_\mathrm{cold}/m_\mathrm{cold, 0}$")
plt.xlabel("$t$ (Myr)")

# * =========================
# * Total mass
# * =========================

plt.figure()
plt.plot(hst_read.time, hst_read.mass_tot - hst_read.mass_tot[0])
plt.axvline(t_stop, ls=":")

plt.yscale("log")
plt.ylabel("$m_\mathrm{tot}$")
plt.xlabel("$t$ (Myr)")

# * =========================
# * Total energy
# * =========================

plt.figure()
plt.plot(hst_read.time, hst_read.E_tot)
plt.axvline(t_stop, ls=":")

# plt.yscale("log")
plt.ylabel("$E_\mathrm{tot}$")
plt.xlabel("$t$ (Myr)")

# * =========================
# * turbulent velocity
# * =========================

plt.figure()
plt.plot(hst_read.time, hst_read.turb_vel)
plt.axvline(t_stop, ls=":")

# plt.yscale("log")
plt.ylabel("$v_\mathrm{turb}$")
plt.xlabel("$t$ (myr)")

# * =========================
# * t_TI vs t_eddy
# * =========================

box_size = 0.1  # in kpc
t_eddy = box_size / hst_read.turb_vel  # in Myr

t_cool = compute_cooling_time(hst_read, cooling_cfg, cooling_tables, cooling_units)
t_eddy_arr = np.atleast_1d(np.asarray(t_eddy, dtype=float))
t_cool_arr = np.atleast_1d(np.asarray(t_cool, dtype=float))

plt.figure()

N_last = None  # -2000

plt.plot(hst_read.time[:N_last], t_eddy[:N_last], label=r"$t_\mathrm{eddy}$")
plt.plot(hst_read.time[:N_last], t_cool[:N_last], label=r"$t_\mathrm{cool}$")
plt.axvline(t_stop, ls=":")

set_log_scale_if_valid(
    plt.gca(),
    "y",
    np.concatenate((t_eddy_arr.ravel(), t_cool_arr[np.isfinite(t_cool_arr)].ravel())),
    "eddy and cooling times",
)
plt.ylabel("$t_\mathrm{eddy}, t_\mathrm{cool}$ (Myr)")
plt.xlabel("$t$ (myr)")

plt.tick_params(
    axis="both",  # applies to both x and y axes
    bottom=True,
    top=True,
    left=True,
    right=True,
)

plt.legend()


# * =========================
# * Average temperature
# * =========================

plt.figure()
plt.plot(hst_read.time, hst_read.T_avg)
plt.axvline(t_stop, ls=":")

set_log_scale_if_valid(plt.gca(), "y", hst_read.T_avg, "average temperature")
plt.ylabel("$T_\mathrm{avg}$")
plt.xlabel("$t$ (Myr)")

# * =========================
# * dt
# * =========================

plt.figure()
plt.plot(hst_read.time, hst_read.dt)
plt.axvline(t_stop, ls=":")

set_log_scale_if_valid(plt.gca(), "y", hst_read.T_avg, "average temperature")
plt.ylabel("$dt$")
plt.xlabel("$t$ (Myr)")


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

size_a = 0.01
size_b = 1.0
a_arr = np.logspace(np.log10(size_a), np.log10(size_b), num=len(Dplot[:, 0]) + 1)
a_mid = 0.5 * (a_arr[:-1] + a_arr[1:])

rho_gr = 3.5  # g/cc
K_dust = 4 / 3 * np.pi * rho_gr

da = a_arr[1:] - a_arr[:-1]

dist_plot = 4 * Ddata / K_dust
dist_plot /= (a_arr[1:] ** 4 - a_arr[:-1] ** 4)[:, np.newaxis]


# Sample a set of line colors from the same CMasher colormap used for images
line_colors = cm.take_cmap_colors(
    cm.rainforest, nscalars, cmap_range=(0.05, 0.95), return_fmt="hex"
)
line_colors_time = cm.take_cmap_colors(
    cm.bubblegum, len(Dplot[0, :]), cmap_range=(0.05, 0.95), return_fmt="hex"
)

# * Image plot species vs time
color_limit = max(np.nanmax(np.abs(Dplot)), np.finfo(float).eps)
plt.figure()
plt.imshow(
    Dplot,
    aspect="auto",
    origin="lower",
    vmin=-color_limit,
    vmax=color_limit,
    cmap=cm.redshift_r,
)
plt.ylabel("Dust Scalar Index")
# plt.xlabel("$t$ (Myr)")
plt.colorbar(label="$\log_{10}(D_i/D_{i,0})$")

# * Dust species distribution
fig_lines, ax_lines = plt.subplots()
for k in range(len(Ddata[0, :])):
    ax_lines.plot(
        a_mid,
        dist_plot[:, k],
        color=line_colors_time[k],
        marker="x",
    )
ax_lines.axvline(t_stop, ls=":")

ax_lines.set_xlim(size_a, size_b)

set_log_scale_if_valid(ax_lines, "y", dist_plot, "dust species distribution")
set_log_scale_if_valid(ax_lines, "x", a_mid, "dust grain size")
ax_lines.set_ylabel("$D_\mathrm{tot}/Z_\odot$")
ax_lines.set_xlabel("$a$ ($\mu$m)")

# Map scalar index -> line color with a compact colorbar
scalar_cmap = ListedColormap(line_colors_time)
norm = plt.Normalize(vmin=-0.5, vmax=nscalars - 0.5)
cbar = fig_lines.colorbar(
    plt.cm.ScalarMappable(norm=norm, cmap=scalar_cmap),
    ax=ax_lines,
    pad=0.02,
)
cbar.set_label("$t$ (Myr)")
cbar.set_ticks(range(nscalars))
cbar.set_ticklabels([str(i) for i in range(nscalars)])


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

# ax_lines.set_ylim(1e-30, None)

set_log_scale_if_valid(
    ax_lines, "y", Ddata / hst_read.mass_tot / Z_solar, "dust species abundance"
)
ax_lines.set_ylabel("$D_\mathrm{tot}/Z_\odot$")
ax_lines.set_xlabel("$t$ (Myr)")

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

set_log_scale_if_valid(
    plt.gca(),
    "y",
    hst_read.dict["Zgas"] / hst_read.mass_tot / Z_solar,
    "gas metallicity",
)
plt.ylabel("$Z_\mathrm{avg}/Z_\odot$")
plt.xlabel("$t$ (Myr)")

# dust_a = 0.001 # Min dust size
# dust_b = 10 # Max dust size

# d_r = np.logspace(np.log10(dust_a), np.log10(dust_b), num=nscalars)
