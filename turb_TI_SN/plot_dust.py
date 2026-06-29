import math
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


theme = "bright"
# theme = "dark"

run_prefix = "Turb.full_hydro_w"
frame_start = 0
frame_stop = 101
slice_dir = "y"
log_flag = True

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

METALLICITY_KEY = "s_00"


def get_dust_keys(out_dict):
    dust_keys = []
    scalar_index = 2

    while True:
        key = f"s_{scalar_index:02d}"
        if key not in out_dict:
            break

        dust_keys.append(key)
        scalar_index += 1

    if not dust_keys:
        raise KeyError("No dust species found. Expected dust scalars starting at s_02.")

    return dust_keys


def get_plane_metadata(box_shape, slice_axis):
    los_axis = axis_to_index[slice_axis]
    plane_axes = [axis for axis in range(3) if axis != los_axis]
    row_axis, col_axis = plane_axes

    x_arr = np.linspace(0.0, 1.0, box_shape[col_axis] + 1)
    y_arr = np.linspace(0.0, 1.0, box_shape[row_axis] + 1)
    X, Y = np.meshgrid(x_arr, y_arr)

    return {
        "los_axis": los_axis,
        "row_axis": row_axis,
        "col_axis": col_axis,
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


def get_ratio_data(current_data, initial_data):
    return np.divide(
        current_data,
        initial_data,
        out=np.full_like(current_data, np.nan, dtype=float),
        where=initial_data != 0,
    )


def get_plot_data(data_2d):
    if not log_flag:
        return data_2d

    return np.log10(np.clip(data_2d, 1e-30, None))


def get_plot_norm(plot_arr):
    finite_values = plot_arr[np.isfinite(plot_arr)]
    if finite_values.size == 0:
        return None

    max_abs = np.max(np.abs(finite_values))
    if max_abs == 0.0:
        return None

    return colors.TwoSlopeNorm(vmin=-max_abs, vcenter=0.0, vmax=max_abs)


def get_panel_configs(dust_keys):
    panel_configs = [{"kind": "metallicity", "title": "Metallicity"}]

    for key in dust_keys:
        species_index = int(key.split("_")[1]) - 2
        panel_configs.append(
            {
                "kind": "dust",
                "key": key,
                "species_index": species_index,
                "title": rf"Dust species {species_index}",
            }
        )

    return panel_configs


def get_colorbar_label(panel_config, plot_mode):
    if panel_config["kind"] == "metallicity":
        if log_flag:
            return r"$\log_{10}(Z/Z_0)$"

        return r"$Z/Z_0$"

    species_index = panel_config["species_index"]

    if plot_mode == "slice":
        if log_flag:
            return rf"$\log_{{10}}(D_{{{species_index}}}/D_{{{species_index},0}})$"

        return rf"$D_{{{species_index}}}/D_{{{species_index},0}}$"

    if log_flag:
        return (
            rf"$\log_{{10}}(\langle D_{{{species_index}}} / D_{{{species_index},0}}"
            rf"\rangle)$"
        )

    return rf"$\langle D_{{{species_index}}} / D_{{{species_index},0}} \rangle$"


def make_figure(
    out_dict,
    initial_metallicity_data,
    initial_dust_data,
    panel_configs,
    plane_meta,
    plot_mode,
):
    n_panels = len(panel_configs)
    ncols = min(4, n_panels)
    nrows = math.ceil(n_panels / ncols)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(5.0 * ncols, 4.5 * nrows),
        squeeze=False,
    )
    axes = axes.flatten()

    for panel_index, panel_config in enumerate(panel_configs):
        ax = axes[panel_index]

        if panel_config["kind"] == "metallicity":
            if METALLICITY_KEY not in out_dict:
                raise KeyError(
                    f"Missing metallicity scalar {METALLICITY_KEY}. "
                    "Expected metallicity data."
                )

            data_3d = get_ratio_data(out_dict[METALLICITY_KEY], initial_metallicity_data)
        else:
            key = panel_config["key"]
            data_3d = get_ratio_data(out_dict[key], initial_dust_data[key])

        if plot_mode == "slice":
            img_arr = get_slice_data(data_3d, plane_meta["los_axis"])
        else:
            img_arr = get_projection_data(data_3d, plane_meta["los_axis"])

        plot_arr = get_plot_data(img_arr)
        plot_norm = get_plot_norm(plot_arr)
        pc = ax.pcolormesh(
            plane_meta["X"],
            plane_meta["Y"],
            plot_arr,
            cmap=cm.redshift_r,
            norm=plot_norm,
            shading="auto",
        )

        ax.set_title(panel_config["title"])
        ax.set_aspect("equal", adjustable="box")
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(MultipleLocator(0.25))
        ax.set_xlabel(plane_meta["xlabel"])
        ax.set_ylabel(plane_meta["ylabel"])

        fig.colorbar(pc, ax=ax, label=get_colorbar_label(panel_config, plot_mode))

    for ax in axes[n_panels:]:
        ax.set_visible(False)

    fig.suptitle(
        f"{plot_mode.title()} dust and metallicity plots, frame "
        f"{out_dict['NumCycles']} | "
        f"t = {out_dict['Time']:.4f}",
    )
    fig.tight_layout(pad=0.3, rect=(0.0, 0.0, 1.0, 0.97))

    return fig


def save_figure(fig, plot_mode, frame_index):
    save_dir = script_dir / "save" / "dust_species" / plot_mode
    save_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        save_dir / f"dust_{plot_mode}_{frame_index:05d}.png",
        dpi=300,
        bbox_inches="tight",
        pad_inches=0.02,
    )
    plt.close(fig)


def main():
    initial_bin_path = script_dir / "bin" / f"{run_prefix}.00000.bin"
    if not initial_bin_path.exists():
        raise FileNotFoundError(f"Missing initial dust snapshot: {initial_bin_path}")

    initial_out_dict = bc.read_binary_as_athdf(str(initial_bin_path))
    dust_keys = get_dust_keys(initial_out_dict)
    panel_configs = get_panel_configs(dust_keys)
    initial_metallicity_data = np.array(initial_out_dict[METALLICITY_KEY], copy=True)
    initial_dust_data = {
        key: np.array(initial_out_dict[key], copy=True) for key in dust_keys
    }

    print(f"Using initial dust snapshot: {initial_bin_path.name}")
    print(f"Dust species keys: {dust_keys}")

    for frame_index in range(frame_start, frame_stop):
        bin_path = script_dir / "bin" / f"{run_prefix}.{frame_index:05d}.bin"
        if not bin_path.exists():
            print(f"Stopping at missing file: {bin_path}")
            break

        print(f"Processing frame {frame_index}: {bin_path.name}")
        out_dict = bc.read_binary_as_athdf(str(bin_path))

        plane_meta = get_plane_metadata(np.shape(out_dict["dens"]), slice_dir)

        slice_fig = make_figure(
            out_dict,
            initial_metallicity_data,
            initial_dust_data,
            panel_configs,
            plane_meta,
            "slice",
        )
        save_figure(slice_fig, "slice", frame_index)

        projection_fig = make_figure(
            out_dict,
            initial_metallicity_data,
            initial_dust_data,
            panel_configs,
            plane_meta,
            "projection",
        )
        save_figure(projection_fig, "projection", frame_index)

        print("Saved slice and projection plots.")


if __name__ == "__main__":
    main()
