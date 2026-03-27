import os
import sys
from pathlib import Path

import cmasher as cm
import matplotlib as mt
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import MultipleLocator


script_dir = Path(__file__).resolve().parent
repo_root = script_dir.parent
sys.path.append(str(repo_root / "vis" / "python"))

import bin_convert_new as bc
import own_package
import own_package.utils.units as ut


theme = "bright"
# theme = "dark"

run_prefix = "Turb.full_hydro_w"
frame_start = 0
frame_stop = 201
slice_dir = "y"

gamma_minus_1 = 5.0 / 3.0 - 1.0

axis_to_index = {"z": 0, "y": 1, "x": 2}
index_to_axis = {value: key for key, value in axis_to_index.items()}

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

quantity_specs = (
    {
        "key": "temperature",
        "title": "Temperature",
        "label": r"$\log_{10} T$ (K)",
        "log": True,
        "cmap": cm.bubblegum,
    },
    {
        "key": "density",
        "title": "Density",
        "label": r"$\log_{10} \rho$ (amu/cc)",
        "log": True,
        "cmap": cm.bubblegum,
    },
    {
        "key": "pressure",
        "title": "Pressure",
        "label": r"$\log_{10} P$ (code units)",
        "log": True,
        "cmap": cm.bubblegum,
    },
    {
        "key": "velx",
        "title": r"$v_x$",
        "label": r"$v_x$ (code units)",
        "log": False,
        "cmap": cm.redshift_r,
    },
    {
        "key": "vely",
        "title": r"$v_y$",
        "label": r"$v_y$ (code units)",
        "log": False,
        "cmap": cm.redshift_r,
    },
    {
        "key": "velz",
        "title": r"$v_z$",
        "label": r"$v_z$ (code units)",
        "log": False,
        "cmap": cm.redshift_r,
    },
)


def compute_temperature(out_dict):
    return (out_dict["eint"] * gamma_minus_1) / out_dict["dens"] * ut.KELVIN * ut.mu


def compute_pressure(out_dict):
    return gamma_minus_1 * out_dict["eint"]


def get_plane_metadata(box_shape, slice_axis):
    los_axis = axis_to_index[slice_axis]
    plane_axes = [axis for axis in range(3) if axis != los_axis]
    row_axis, col_axis = plane_axes

    x_arr = np.linspace(0.0, 1.0, box_shape[col_axis] + 1)
    y_arr = np.linspace(0.0, 1.0, box_shape[row_axis] + 1)
    X, Y = np.meshgrid(x_arr, y_arr)

    return {
        "los_axis": los_axis,
        "xlabel": index_to_axis[col_axis],
        "ylabel": index_to_axis[row_axis],
        "X": X,
        "Y": Y,
    }


def get_slice_data(data_3d, los_axis):
    slice_selector = [slice(None)] * 3
    slice_selector[los_axis] = data_3d.shape[los_axis] // 2
    return data_3d[tuple(slice_selector)]


def get_projection_data(data_3d, los_axis):
    return np.mean(data_3d, axis=los_axis)


def get_plot_data(data_2d, use_log):
    if use_log:
        return np.log10(np.clip(data_2d, 1e-30, None))

    return data_2d


def get_log_limits(plot_arrays):
    finite_chunks = [
        arr[np.isfinite(arr)] for arr in plot_arrays if np.any(np.isfinite(arr))
    ]
    if not finite_chunks:
        return None, None

    finite_values = np.concatenate(finite_chunks)
    if finite_values.size == 0:
        return None, None

    vmin, vmax = np.nanpercentile(finite_values, [1.0, 99.0])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin == vmax:
        return None, None

    return vmin, vmax


def get_velocity_norm(plot_arrays):
    finite_chunks = [
        arr[np.isfinite(arr)] for arr in plot_arrays if np.any(np.isfinite(arr))
    ]
    if not finite_chunks:
        return None

    finite_values = np.concatenate(finite_chunks)
    if finite_values.size == 0:
        return None

    max_abs = np.nanpercentile(np.abs(finite_values), 99.0)
    if max_abs == 0.0:
        return None

    return colors.TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)


def get_quantity_data(out_dict):
    return {
        "temperature": compute_temperature(out_dict),
        "density": out_dict["dens"],
        "pressure": compute_pressure(out_dict),
        "velx": out_dict["velx"],
        "vely": out_dict["vely"],
        "velz": out_dict["velz"],
    }


def make_figure(out_dict, plane_meta):
    quantity_data = get_quantity_data(out_dict)

    fig = plt.figure(figsize=(26.0, 9.8))
    grid = fig.add_gridspec(
        nrows=3,
        ncols=len(quantity_specs),
        height_ratios=(1.0, 1.0, 0.06),
        hspace=0.28,
        wspace=0.30,
    )
    axes = np.empty((2, len(quantity_specs)), dtype=object)
    colorbar_axes = []

    for col_index in range(len(quantity_specs)):
        axes[0, col_index] = fig.add_subplot(grid[0, col_index])
        axes[1, col_index] = fig.add_subplot(grid[1, col_index])
        colorbar_ax = fig.add_subplot(grid[2, col_index])
        cax_pos = colorbar_ax.get_position()
        colorbar_ax.set_position(
            [cax_pos.x0, cax_pos.y0 - 0.012, cax_pos.width, cax_pos.height * 0.85]
        )
        colorbar_axes.append(colorbar_ax)

    row_names = ("Slice", "Projection")

    for col_index, spec in enumerate(quantity_specs):
        data_3d = quantity_data[spec["key"]]
        slice_arr = get_slice_data(data_3d, plane_meta["los_axis"])
        projection_arr = get_projection_data(data_3d, plane_meta["los_axis"])

        plot_arrays = [
            get_plot_data(slice_arr, spec["log"]),
            get_plot_data(projection_arr, spec["log"]),
        ]

        vmin = vmax = norm = None
        if spec["log"]:
            vmin, vmax = get_log_limits(plot_arrays)
        else:
            norm = get_velocity_norm(plot_arrays)

        mappable = None
        for row_index, plot_arr in enumerate(plot_arrays):
            ax = axes[row_index, col_index]
            pc = ax.pcolormesh(
                plane_meta["X"],
                plane_meta["Y"],
                plot_arr,
                cmap=spec["cmap"],
                vmin=vmin,
                vmax=vmax,
                norm=norm,
                shading="auto",
            )
            mappable = pc

            ax.set_aspect("equal", adjustable="box")
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.yaxis.set_major_locator(MultipleLocator(0.25))
            ax.set_ylabel(plane_meta["ylabel"])

            if row_index == 0:
                ax.set_title(spec["title"])
            else:
                ax.set_xlabel(plane_meta["xlabel"])

        axes[0, col_index].text(
            0.02,
            0.98,
            row_names[0],
            transform=axes[0, col_index].transAxes,
            ha="left",
            va="top",
            bbox={
                "facecolor": fig_face_color,
                "edgecolor": line_border_color,
                "pad": 2.5,
            },
        )
        axes[1, col_index].text(
            0.02,
            0.98,
            row_names[1],
            transform=axes[1, col_index].transAxes,
            ha="left",
            va="top",
            bbox={
                "facecolor": fig_face_color,
                "edgecolor": line_border_color,
                "pad": 2.5,
            },
        )

        colorbar = fig.colorbar(
            mappable,
            cax=colorbar_axes[col_index],
            orientation="horizontal",
        )
        colorbar.set_label(spec["label"], labelpad=4)
        colorbar.ax.xaxis.set_label_position("top")
        colorbar.ax.tick_params(axis="x", pad=1)

    fig.suptitle(
        f"Gas properties, frame {out_dict['NumCycles']} | t = {out_dict['Time']:.4f}",
    )

    return fig


def save_figure(fig, frame_index):
    save_dir = script_dir / "save" / "properties"
    save_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        save_dir / f"properties_{frame_index:05d}.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.02,
    )
    plt.close(fig)


def main():
    for frame_index in range(frame_start, frame_stop):
        bin_path = script_dir / "bin" / f"{run_prefix}.{frame_index:05d}.bin"
        if not bin_path.exists():
            print(f"Stopping at missing file: {bin_path}")
            break

        print(f"Processing frame {frame_index}: {bin_path.name}")
        out_dict = bc.read_binary_as_athdf(str(bin_path))

        plane_meta = get_plane_metadata(np.shape(out_dict["dens"]), slice_dir)
        fig = make_figure(out_dict, plane_meta)
        save_figure(fig, frame_index)

        print("Saved property panel.")


if __name__ == "__main__":
    main()
