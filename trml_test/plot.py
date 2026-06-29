import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mt
import gc
import os

import sys
from pathlib import Path

repo_root = Path.cwd().resolve().parents[0]  # adjust if running in notebook
sys.path.append(str(repo_root / "vis" / "python"))

import cmasher as cm
import bin_convert_new as bc

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

output_dir = Path("save/temp")
output_dir.mkdir(parents=True, exist_ok=True)


def compute_temperature(out_dict):
    KE = (
        0.5
        * out_dict["dens"]
        * (out_dict["velx"] ** 2 + out_dict["vely"] ** 2 + out_dict["velz"] ** 2)
    )
    TE = out_dict["eint"] - KE

    T = (TE * gamma_minus_1) / out_dict["dens"]

    return T


metal_tot = []
gas_tot = []
time = []

for i in range(0, 101):
    print(f"Processing frame {i}")

    dir_dict = {"x": 2, "y": 1, "z": 0}

    slice_dir = "y"

    out_dict = bc.read_binary_as_athdf(f"bin/TRML.slice_{slice_dir}.{i:05d}.bin")

    i_s = dir_dict[slice_dir]
    a = (i_s + 1) % 3
    b = (i_s + 2) % 3

    slice_list = [0] * 3
    slice_list[a] = slice(None, None, None)
    slice_list[b] = slice(None, None, None)

    arr_shape = np.shape(out_dict["dens"])

    if a < b:
        L_mul = int(arr_shape[a] / arr_shape[b])

        x_arr = np.linspace(0, L_mul, arr_shape[a] + 1)
        y_arr = np.linspace(0, 1, arr_shape[b] + 1)
    else:
        L_mul = int(arr_shape[b] / arr_shape[a])

        x_arr = np.linspace(0, L_mul, arr_shape[b] + 1)
        y_arr = np.linspace(0, 1, arr_shape[a] + 1)

    X, Y = np.meshgrid(y_arr, x_arr)

    fig, ax = plt.subplots(figsize=(12, 12))

    dens = out_dict["dens"][tuple(slice_list)]
    s_01 = out_dict["s_01"][tuple(slice_list)]
    img_arr = np.divide(s_01, dens, out=np.full_like(s_01, np.nan), where=dens != 0) / Z_solar

    metal_tot.append(np.sum(out_dict["s_01"]))
    gas_tot.append(np.sum(out_dict["dens"]))
    time.append(out_dict["Time"])

    pc = ax.pcolormesh(
        X,
        Y,
        img_arr,
        # compute_temperature(out_dict)[*slice_list],
        # out_dict["s_01"][*slice_list],
        # (out_dict["s_00"] / out_dict["dens"])[0, :, :],
        vmin=0.5,
        vmax=1.5,
        cmap=cm.bubblegum,
        # cmap=cm.redshift,
        shading="auto",
    )
    ax.set_aspect("equal", adjustable="box")
    fig.colorbar(pc, ax=ax, label=r"$Z/Z_\odot$")

    # after plotting
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))

    ax.set_ylabel("y")
    ax.set_xlabel("x")

    fig.tight_layout(pad=0.25)

    fig.savefig(output_dir / f"metal_{i:05d}.png", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    # del out_dict, X, Y, fig
    # gc.collect()

# %%

gas_tot = np.array(gas_tot)
metal_tot = np.array(metal_tot)

plt.figure()
plt.plot(time, metal_tot / metal_tot[0], label="Metal mass")
plt.legend()

plt.ylabel(r"$M_Z/M_{Z,0}$")
plt.xlabel(f"$t$ (code units)")

plt.ylim(0, 1.1)
plt.xlim(0, None)


plt.figure()
plt.plot(time, metal_tot / gas_tot / Z_solar, label="Mean metallicity")
plt.legend()

plt.ylabel(r"$Z/Z_\odot$")
plt.xlabel(f"$t$ (code units)")

plt.xlim(0, None)
# %%
