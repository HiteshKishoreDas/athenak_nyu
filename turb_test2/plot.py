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

import own_package.utils.units as ut

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


def compute_temperature(out_dict):

    T = (out_dict["eint"] * gamma_minus_1) / out_dict["dens"] * ut.KELVIN * ut.mu
    # T = (143 * gamma_minus_1) / 100.0 * ut.KELVIN * ut.mu

    # plt.figure()
    # plt.hist(out_dict["eint"].flatten(), bins=100)
    # plt.show()

    # print("===========================")
    # print(f"{T = }")

    return T


dust_tot_small = []
dust_tot_large = []
time = []

for i in range(0, 40):
    print(f"Processing frame {i}")

    dir_dict = {"x": 2, "y": 1, "z": 0}

    slice_dir = "y"

    quantity = "l"  # "s" or "l" or "z" or "ratio"

    out_dict = bc.read_binary_as_athdf(f"bin/Turb.full_hydro_w.{i:05d}.bin")

    box_shape = np.shape(out_dict["dens"])

    i_s = dir_dict[slice_dir]
    a = (i_s + 1) % 3
    b = (i_s + 2) % 3

    slice_list = [box_shape[i_s] // 2] * 3
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

    # key = "dens"
    key = "T"
    # key = "eint"

    slice_flag = True
    # slice_flag = False

    if not slice_flag:
        slice_list = [slice(None, None, None)] * 3

    if key != "T":
        img_arr = out_dict[key][*slice_list]
    else:
        img_arr = compute_temperature(out_dict)[*slice_list]

    if not slice_flag:
        img_arr = np.sum(img_arr, axis=i_s)

    img_arr = np.log10(img_arr)

    # img_arr = out_dict["s_01"][*slice_list]
    # img_arr = out_dict["s_02"][*slice_list]
    # img_arr = out_dict["s_03"][*slice_list]
    # img_arr = out_dict["s_02"][*slice_list] / out_dict["s_03"][*slice_list]
    # img_arr /= Z_solar

    # img_arr = out_dict["s_02"][*slice_list] / img_arr

    dust_tot_small.append(np.sum(out_dict["dens"] * out_dict["s_02"]))
    dust_tot_large.append(np.sum(out_dict["dens"] * out_dict["s_03"]))
    time.append(out_dict["Time"])

    pc = ax.pcolormesh(
        X,
        Y,
        img_arr,
        # compute_temperature(out_dict)[*slice_list],
        # out_dict["s_01"][*slice_list],
        # (out_dict["s_00"] / out_dict["dens"])[0, :, :],
        # vmin=2.5,
        # vmax=3.5,
        vmin=4.0,
        vmax=6.0,
        cmap=cm.bubblegum,
        # cmap=cm.redshift,
        shading="auto",
    )

    ax.set_aspect("equal", adjustable="box")
    fig.colorbar(pc, ax=ax, label="T")
    # fig.colorbar(pc, ax=ax, label=f"$D_l/Z_\odot$")
    # fig.colorbar(pc, ax=ax, label=f"$T$ (K)")
    # fig.colorbar(pc, ax=ax, label=f"rho (amu/cc)")
    # fig.colorbar(pc, ax=ax, label=f"eint")

    # after plotting
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))

    ax.set_ylabel("y")
    ax.set_xlabel("x")

    fig.tight_layout(pad=0.25)

    fig.savefig(
        # f"save/dust/D_r_{i:05d}.png",
        # f"save/dens/dens_{i:05d}.png",
        # f"save/eint/eint_{i:05d}.png",
        f"save/{key}/{key}_{i:05d}.png",
        dpi=600,
        bbox_inches="tight",
        pad_inches=0.02,
    )
    plt.close(fig)

    # del out_dict, X, Y, fig
    # gc.collect()

# %%

dust_tot_large = np.array(dust_tot_large)
dust_tot_small = np.array(dust_tot_small)

plt.figure()
plt.plot(time, dust_tot_small / dust_tot_small[0], label="Small dust")
plt.plot(time, dust_tot_large / dust_tot_large[0], label="Large dust")
plt.legend()

plt.ylabel(r"$D/D_0$")
plt.xlabel(f"$t$ (code units)")

plt.ylim(0, 1.1)
plt.xlim(0, None)


plt.figure()
plt.plot(time, dust_tot_small / dust_tot_large)
plt.legend()

plt.ylabel(r"$D_s/D_l$")
plt.xlabel(f"$t$ (code units)")

plt.xlim(0, None)
# %%
