import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = Path(__file__).resolve().parent / ".cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)
MATPLOTLIB_CACHE_DIR = CACHE_DIR / "matplotlib"
MATPLOTLIB_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MATPLOTLIB_CACHE_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(CACHE_DIR))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

try:
    import cmasher as cmr
except ImportError:
    cmr = None

VIS_PYTHON = ROOT / "vis" / "python"
if str(VIS_PYTHON) not in sys.path:
    sys.path.insert(0, str(VIS_PYTHON))

from bin_convert import read_binary
from plot_utils import plot_amr_slice_patchwork


DEFAULT_INPUT = Path(__file__).resolve().parent / "gotham.slice_x2.00240.bin"
DEFAULT_OUTPUT_DIR = Path(__file__).resolve().parent / "plots"
DEFAULT_PRIMITIVE_VARIABLES = ["dens", "velx", "vely", "velz", "eint", "s_00"]
DEFAULT_SLICE_POSITION = 0.0
DEFAULT_AXES_LIMIT = 30.0
DEFAULT_TIMESCALE_LIMITS = (1.0e-3, 1.0e4)
DEFAULT_RATIO_LIMITS = (1.0e-3, 1.0e3)
DEFAULT_RATIO_COLORMAP = "cmr.redshift_r" if cmr is not None else "RdBu_r"
LOG_SCALE_VARS = {"dens", "eint"}
COLORMAPS = {
    "dens": "viridis",
    "velx": "coolwarm",
    "vely": "coolwarm",
    "velz": "coolwarm",
    "eint": "magma",
    "s_00": "cividis",
    "T_K": "inferno",
}


def get_arrays_from_bin_file(input_file, variables=None):
    """
    Read a slice binary and return the snapshot metadata plus the requested arrays.

    The returned arrays are the raw per-meshblock arrays stored in the binary:
    ``arrays[variable][mb_id]`` is the NumPy array for one mesh block.
    """

    input_file = Path(input_file).resolve()
    snapshot = read_binary(str(input_file))

    if variables is None:
        variables = [
            variable
            for variable in DEFAULT_PRIMITIVE_VARIABLES
            if variable in snapshot["var_names"]
        ]

    missing = [
        variable for variable in variables if variable not in snapshot["var_names"]
    ]
    if missing:
        raise ValueError(
            f"Requested variables not found in {input_file.name}: {missing}. "
            f"Available variables: {snapshot['var_names']}"
        )

    arrays = {variable: snapshot["mb_data"][variable] for variable in variables}
    return snapshot, arrays


def infer_slice_axis(snapshot_path):
    name = snapshot_path.name
    if ".slice_x1." in name:
        return "x"
    if ".slice_x2." in name:
        return "y"
    if ".slice_x3." in name:
        return "z"
    return "z"


def build_output_path(snapshot_path, output_dir):
    stem = snapshot_path.name.replace(".bin", "")
    return output_dir / f"{stem}_primitive_variables.png"


def build_single_array_snapshot(snapshot, array, variable_name):
    """Return a snapshot-like dict with one variable mapped to the supplied array list."""

    if len(array) != snapshot["n_mbs"]:
        raise ValueError(
            f"Expected {snapshot['n_mbs']} meshblock arrays for '{variable_name}', "
            f"got {len(array)}."
        )

    snapshot_copy = dict(snapshot)
    snapshot_copy["var_names"] = list(snapshot["var_names"])
    snapshot_copy["mb_data"] = dict(snapshot["mb_data"])
    snapshot_copy["mb_data"][variable_name] = array
    if variable_name not in snapshot_copy["var_names"]:
        snapshot_copy["var_names"].append(variable_name)
    return snapshot_copy


def array_supports_lognorm(array):
    """Return True when the meshblock array contains at least one positive value and no negatives."""

    has_positive = False
    for block in array:
        values = np.asarray(block)
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            continue
        if np.any(finite < 0.0):
            return False
        if np.any(finite > 0.0):
            has_positive = True
    return has_positive


def is_timescale_variable(variable_name):
    return variable_name.startswith("t")


def plot_primitive_variables(
    snapshot,
    array,
    variable_name,
    output_path=None,
    input_file=None,
    slice_position=DEFAULT_SLICE_POSITION,
    axes_limit=DEFAULT_AXES_LIMIT,
    log_scale=None,
    colormap=None,
    vmin=None,
    vmax=None,
):
    """Plot one meshblock-array field from a loaded snapshot."""

    input_file = Path(input_file).resolve() if input_file is not None else DEFAULT_INPUT
    plot_snapshot = build_single_array_snapshot(snapshot, array, variable_name)
    slice_axis = infer_slice_axis(input_file)
    plot_limit = axes_limit / 2.0
    fig, ax = plt.subplots(figsize=(7.5, 6.0))

    if log_scale is None:
        log_scale = variable_name in LOG_SCALE_VARS or array_supports_lognorm(array)
    if colormap is None:
        colormap = COLORMAPS.get(variable_name, "viridis")
    if vmin is None and vmax is None and is_timescale_variable(variable_name):
        vmin, vmax = DEFAULT_TIMESCALE_LIMITS

    plot_amr_slice_patchwork(
        plot_snapshot,
        variable=variable_name,
        slice_axis=slice_axis,
        slice_position=slice_position,
        log_scale=log_scale,
        vmin=vmin,
        vmax=vmax,
        glw=0.2,
        show_grid=False,
        colormap=colormap,
        ax=ax,
        add_colorbar=True,
        title=variable_name,
    )
    ax.set_xlim(-plot_limit, plot_limit)
    ax.set_ylim(-plot_limit, plot_limit)
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Saved plot to {output_path}")
    return fig, ax


def plot_timescale_multiplot(
    snapshot,
    arrays,
    accretion_suffix="",
    accretion_label="acc",
    output_path=None,
    input_file=None,
    slice_position=DEFAULT_SLICE_POSITION,
    axes_limit=DEFAULT_AXES_LIMIT,
):
    """Plot the timescale arrays in a fixed 2x3 layout."""

    panel_order = [
        ("tsputt_small", 0, 0),
        ("tsputt_large", 1, 0),
        (f"tacc_small{accretion_suffix}", 0, 1),
        (f"tacc_large{accretion_suffix}", 1, 1),
        ("tcoag", 0, 2),
        ("tshatt", 1, 2),
    ]
    available_names = [name for name, _, _ in panel_order if name in arrays]
    if not available_names:
        raise ValueError("No timescale arrays were provided for the multiplot.")

    input_file = Path(input_file).resolve() if input_file is not None else DEFAULT_INPUT
    slice_axis = infer_slice_axis(input_file)
    plot_limit = axes_limit / 2.0

    fig, axes = plt.subplots(2, 3, figsize=(20.4, 11.6))

    for variable_name, row, col in panel_order:
        ax = axes[row, col]
        if variable_name not in arrays:
            ax.remove()
            continue

        plot_snapshot = build_single_array_snapshot(
            snapshot, arrays[variable_name], variable_name
        )
        plot_amr_slice_patchwork(
            plot_snapshot,
            variable=variable_name,
            slice_axis=slice_axis,
            slice_position=slice_position,
            log_scale=True,
            vmin=DEFAULT_TIMESCALE_LIMITS[0],
            vmax=DEFAULT_TIMESCALE_LIMITS[1],
            glw=0.2,
            show_grid=False,
            colormap=COLORMAPS.get(variable_name, "viridis"),
            ax=ax,
            add_colorbar=True,
            title=variable_name.replace("tacc", accretion_label),
        )
        ax.set_xlim(-plot_limit, plot_limit)
        ax.set_ylim(-plot_limit, plot_limit)

    fig.suptitle(f"{input_file.name} timescales ({accretion_label})", fontsize=16)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.97])

    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Saved plot to {output_path}")
    return fig, axes


def divide_meshblock_arrays(numerator, denominator):
    """Divide two meshblock-array lists safely, returning NaN where the denominator is zero."""

    ratios = []
    for num_block, den_block in zip(numerator, denominator):
        num_values = np.asarray(num_block, dtype=float)
        den_values = np.asarray(den_block, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_block = np.where(den_values != 0.0, num_values / den_values, np.nan)
        ratios.append(ratio_block)
    return ratios


def build_timescale_ratio_arrays(arrays, accretion_suffix=""):
    """Build the three requested ratio fields for one accretion prescription."""

    return {
        "tsputt_small_over_tacc_small": divide_meshblock_arrays(
            arrays["tsputt_small"], arrays[f"tacc_small{accretion_suffix}"]
        ),
        "tsputt_large_over_tacc_large": divide_meshblock_arrays(
            arrays["tsputt_large"], arrays[f"tacc_large{accretion_suffix}"]
        ),
        "tcoag_over_tshatt": divide_meshblock_arrays(arrays["tcoag"], arrays["tshatt"]),
    }


def positive_and_negative_magnitude(array):
    """Split a signed meshblock-array list into positive and negative-magnitude views."""

    positive = []
    negative = []
    for block in array:
        values = np.asarray(block, dtype=float)
        positive.append(np.where(values > 0.0, values, np.nan))
        negative.append(np.where(values < 0.0, -values, np.nan))
    return positive, negative


def plot_timescale_ratio_multiplot(
    snapshot,
    ratio_arrays,
    accretion_label="acc",
    output_path=None,
    input_file=None,
    slice_position=DEFAULT_SLICE_POSITION,
    axes_limit=DEFAULT_AXES_LIMIT,
):
    """Plot the requested timescale ratios in a 1x3 multiplot."""

    panel_order = [
        (
            "tsputt_small_over_tacc_small",
            rf"$t_{{\rm sput,small}} / t_{{\rm {accretion_label},small}}$",
            0,
        ),
        (
            "tsputt_large_over_tacc_large",
            rf"$t_{{\rm sput,large}} / t_{{\rm {accretion_label},large}}$",
            1,
        ),
        ("tcoag_over_tshatt", r"$t_{\rm coag} / t_{\rm shatt}$", 2),
    ]
    available_names = [name for name, _, _ in panel_order if name in ratio_arrays]
    if not available_names:
        raise ValueError("No timescale ratio arrays were provided for the multiplot.")

    input_file = Path(input_file).resolve() if input_file is not None else DEFAULT_INPUT
    slice_axis = infer_slice_axis(input_file)
    plot_limit = axes_limit / 2.0

    fig, axes = plt.subplots(1, 3, figsize=(20.4, 6.2))
    axes = list(getattr(axes, "flat", [axes]))

    for variable_name, title, col in panel_order:
        ax = axes[col]
        if variable_name not in ratio_arrays:
            ax.remove()
            continue

        plot_snapshot = build_single_array_snapshot(
            snapshot, ratio_arrays[variable_name], variable_name
        )
        plot_amr_slice_patchwork(
            plot_snapshot,
            variable=variable_name,
            slice_axis=slice_axis,
            slice_position=slice_position,
            log_scale=True,
            vmin=DEFAULT_RATIO_LIMITS[0],
            vmax=DEFAULT_RATIO_LIMITS[1],
            glw=0.2,
            show_grid=False,
            colormap=DEFAULT_RATIO_COLORMAP,
            ax=ax,
            add_colorbar=True,
            title=title,
        )
        ax.set_xlim(-plot_limit, plot_limit)
        ax.set_ylim(-plot_limit, plot_limit)

    fig.suptitle(f"{input_file.name} timescale ratios", fontsize=16)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Saved plot to {output_path}")
    return fig, axes


def plot_tequiv_multiplot(
    snapshot,
    arrays,
    suffix="",
    label="tequiv",
    output_path=None,
    input_file=None,
    slice_position=DEFAULT_SLICE_POSITION,
    axes_limit=DEFAULT_AXES_LIMIT,
):
    """Plot the equilibrium timescales with positive and negative branches split by row."""

    panel_order = [(f"tequiv_small{suffix}", 0), (f"tequiv_large{suffix}", 1)]
    available_names = [name for name, _ in panel_order if name in arrays]
    if not available_names:
        raise ValueError("No tequiv arrays were provided for the multiplot.")

    input_file = Path(input_file).resolve() if input_file is not None else DEFAULT_INPUT
    slice_axis = infer_slice_axis(input_file)
    plot_limit = axes_limit / 2.0

    fig, axes = plt.subplots(2, 2, figsize=(13.6, 11.8))

    for variable_name, col in panel_order:
        if variable_name not in arrays:
            axes[0, col].remove()
            axes[1, col].remove()
            continue

        positive, negative = positive_and_negative_magnitude(arrays[variable_name])

        positive_snapshot = build_single_array_snapshot(
            snapshot, positive, variable_name
        )
        plot_amr_slice_patchwork(
            positive_snapshot,
            variable=variable_name,
            slice_axis=slice_axis,
            slice_position=slice_position,
            log_scale=True,
            vmin=DEFAULT_TIMESCALE_LIMITS[0],
            vmax=DEFAULT_TIMESCALE_LIMITS[1],
            glw=0.2,
            show_grid=False,
            colormap=COLORMAPS.get(variable_name, "viridis"),
            ax=axes[0, col],
            add_colorbar=True,
            title=variable_name.replace("tequiv", label),
        )
        axes[0, col].set_xlim(-plot_limit, plot_limit)
        axes[0, col].set_ylim(-plot_limit, plot_limit)

        negative_name = f"neg_{variable_name}"
        negative_snapshot = build_single_array_snapshot(snapshot, negative, negative_name)
        plot_amr_slice_patchwork(
            negative_snapshot,
            variable=negative_name,
            slice_axis=slice_axis,
            slice_position=slice_position,
            log_scale=True,
            vmin=DEFAULT_TIMESCALE_LIMITS[0],
            vmax=DEFAULT_TIMESCALE_LIMITS[1],
            glw=0.2,
            show_grid=False,
            colormap=COLORMAPS.get(variable_name, "viridis"),
            ax=axes[1, col],
            add_colorbar=True,
            title=f"-{variable_name.replace('tequiv', label)}",
        )
        axes[1, col].set_xlim(-plot_limit, plot_limit)
        axes[1, col].set_ylim(-plot_limit, plot_limit)

    fig.suptitle(f"{input_file.name} {label}", fontsize=16)
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])

    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"Saved plot to {output_path}")
    return fig, axes


if __name__ == "__main__":
    import units as un
    import dust_rates as dr

    input_file = DEFAULT_INPUT.resolve()
    output_dir = DEFAULT_OUTPUT_DIR.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    snapshot, arrays = get_arrays_from_bin_file(input_file, DEFAULT_PRIMITIVE_VARIABLES)

    dens = arrays["dens"]
    eint = arrays["eint"]
    Zg = arrays["s_00"]

    T_K = [un.temp_cgs * un.gm1 * einti / densi for einti, densi in zip(eint, dens)]
    arrays["T_K"] = T_K

    tsput = [dr.sputtering(rhog, T, Zg) for rhog, T, Zg in zip(dens, T_K, Zg)]
    arrays["tsputt_small"] = [t[0] for t in tsput]
    arrays["tsputt_large"] = [t[1] for t in tsput]

    tacc = [dr.accretion(rhog, T, Zg) for rhog, T, Zg in zip(dens, T_K, Zg)]
    arrays["tacc_small"] = [t[0] for t in tacc]
    arrays["tacc_large"] = [t[1] for t in tacc]

    tacc_new = [dr.new_accretion(rhog, T, Zg) for rhog, T, Zg in zip(dens, T_K, Zg)]
    arrays["tacc_small_new"] = [t[0] for t in tacc_new]
    arrays["tacc_large_new"] = [t[1] for t in tacc_new]

    arrays["tshatt"] = [
        dr.shattering(rhog, T, Zg) for rhog, T, Zg in zip(dens, T_K, Zg)
    ]

    arrays["tcoag"] = [
        dr.coagulation(rhog, T, Zg) for rhog, T, Zg in zip(dens, T_K, Zg)
    ]

    arrays["tequiv_small"] = [
        1 / (-1 / t_sput + 1 / t_acc)
        for t_sput, t_acc in zip(arrays["tsputt_small"], arrays["tacc_small"])
    ]
    arrays["tequiv_large"] = [
        1 / (-1 / t_sput + 1 / t_acc)
        for t_sput, t_acc in zip(arrays["tsputt_large"], arrays["tacc_large"])
    ]
    arrays["tequiv_small_new"] = [
        1 / (-1 / t_sput + 1 / t_acc)
        for t_sput, t_acc in zip(arrays["tsputt_small"], arrays["tacc_small_new"])
    ]
    arrays["tequiv_large_new"] = [
        1 / (-1 / t_sput + 1 / t_acc)
        for t_sput, t_acc in zip(arrays["tsputt_large"], arrays["tacc_large_new"])
    ]

    ratio_arrays = build_timescale_ratio_arrays(arrays, "")
    ratio_arrays_new = build_timescale_ratio_arrays(arrays, "_new")

    for variable_name, array in arrays.items():
        if is_timescale_variable(variable_name):
            continue
        output_path = output_dir / f"{input_file.stem}_{variable_name}.png"
        plot_primitive_variables(
            snapshot=snapshot,
            array=array,
            variable_name=variable_name,
            output_path=output_path,
            input_file=input_file,
        )

    output_path = output_dir / f"{input_file.stem}_timescales.png"
    plot_timescale_multiplot(
        snapshot=snapshot,
        arrays=arrays,
        accretion_suffix="",
        accretion_label="acc",
        output_path=output_path,
        input_file=input_file,
    )

    output_path = output_dir / f"{input_file.stem}_timescales_new.png"
    plot_timescale_multiplot(
        snapshot=snapshot,
        arrays=arrays,
        accretion_suffix="_new",
        accretion_label="acc_new",
        output_path=output_path,
        input_file=input_file,
    )

    output_path = output_dir / f"{input_file.stem}_tequiv.png"
    plot_tequiv_multiplot(
        snapshot=snapshot,
        arrays=arrays,
        suffix="",
        label="tequiv",
        output_path=output_path,
        input_file=input_file,
    )

    output_path = output_dir / f"{input_file.stem}_tequiv_new.png"
    plot_tequiv_multiplot(
        snapshot=snapshot,
        arrays=arrays,
        suffix="_new",
        label="tequiv_new",
        output_path=output_path,
        input_file=input_file,
    )

    output_path = output_dir / f"{input_file.stem}_timescale_ratios.png"
    plot_timescale_ratio_multiplot(
        snapshot=snapshot,
        ratio_arrays=ratio_arrays,
        accretion_label="acc",
        output_path=output_path,
        input_file=input_file,
    )

    output_path = output_dir / f"{input_file.stem}_timescale_ratios_new.png"
    plot_timescale_ratio_multiplot(
        snapshot=snapshot,
        ratio_arrays=ratio_arrays_new,
        accretion_label="acc,new",
        output_path=output_path,
        input_file=input_file,
    )
