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

# theme = "bright"
theme = "dark"

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


def compute_temperature(out_dict):
    KE = (
        0.5
        * out_dict["dens"]
        * (out_dict["velx"] ** 2 + out_dict["vely"] ** 2 + out_dict["velz"] ** 2)
    )
    TE = out_dict["eint"] - KE

    T = (TE * gamma_minus_1) / out_dict["dens"]

    return T


for i in range(100, 200):
    print(f"Processing frame {i}")

    fig, ax = plt.subplots(figsize=(16, 3))
    out_dict = bc.read_binary_as_athdf(f"bin/Conv.hydro_w.{i:05d}.bin")

    arr_shape = out_dict["dens"].shape[1:]
    L_mul = int(arr_shape[0] / arr_shape[1])

    x_arr = np.linspace(0, 1, arr_shape[1] + 1)
    y_arr = np.linspace(0, L_mul, arr_shape[0] + 1)

    X, Y = np.meshgrid(y_arr, x_arr)

    pc = ax.pcolormesh(
        X,
        Y,
        compute_temperature(out_dict)[0, :, :].T,
        # (out_dict["s_00"] / out_dict["dens"])[0, :, :],
        vmin=9.5,
        vmax=10.5,
        cmap=cm.bubblegum,
        # cmap=cm.redshift,
        shading="auto",
    )
    ax.set_aspect("equal", adjustable="box")
    fig.colorbar(pc, ax=ax, label="T")

    # after plotting
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(0.25))

    ax.set_ylabel("y")
    ax.set_xlabel("x")

    fig.tight_layout(pad=0.25)

    fig.savefig(
        f"bin/save/temp/temp_{i:05d}.png", dpi=600, bbox_inches="tight", pad_inches=0.02
    )
    plt.close(fig)

    del out_dict, X, Y, fig
    gc.collect()
