import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


def plot_amr_slice_patchwork(
    my_data,
    variable="dens",
    slice_axis="z",
    slice_position=0.0,
    log_scale=True,
    vmin=None,
    vmax=None,
    figsize=(10, 8),
    glw=1.0,
    show_grid=False,
    gamma=5.0 / 3.0,
    unit_temp=1.0,
    scale_factor=1.0,
    colormap="viridis",
    show_quiver=False,
    quiver_color="white",
    quiver_scale=None,
    quiver_stride=2,
    quiver_alpha=0.8,
    quiver_normalize=False,
    ax=None,
    add_colorbar=True,
    title=None,
    norm=None,
):
    """
    Plot a 2D AMR slice by stitching together the intersected mesh blocks.

    Parameters
    ----------
    my_data : dict
        AthenaK binary data loaded via ``read_binary``.
    variable : str, default "dens"
        Variable name to plot. The special value ``"temp"`` derives temperature
        from ``dens`` and ``eint``.
    slice_axis : {"x", "y", "z"}, default "z"
        Axis normal to the plotted slice.
    slice_position : float, default 0.0
        Coordinate along ``slice_axis`` to sample.
    log_scale : bool, default True
        Plot the sampled data with a logarithmic color normalization.
    vmin, vmax : float or None
        Explicit color limits. If omitted, limits are inferred from the slice.
    figsize : tuple[float, float], default (10, 8)
        Figure size used when a new figure is created.
    glw : float, default 1.0
        Grid-line width scaling factor.
    show_grid : bool, default False
        Shrink block extents slightly to reveal AMR block boundaries.
    gamma : float, default 5/3
        Adiabatic index used for the derived temperature.
    unit_temp : float, default 1.0
        Temperature unit conversion factor.
    scale_factor : float, default 1.0
        Multiplicative scale applied to non-temperature variables.
    colormap : str, default "viridis"
        Matplotlib colormap name.
    show_quiver : bool, default False
        Overlay in-plane velocity arrows.
    quiver_color : str, default "white"
        Quiver arrow color.
    quiver_scale : float or None, default None
        Forwarded to ``Axes.quiver``.
    quiver_stride : int, default 2
        Sample every Nth cell when drawing quivers.
    quiver_alpha : float, default 0.8
        Quiver arrow transparency.
    quiver_normalize : bool, default False
        Normalize the velocity arrows before plotting.
    ax : matplotlib.axes.Axes or None, default None
        Existing axes to draw on. A new figure and axes are created if omitted.
    add_colorbar : bool, default True
        Whether to attach a colorbar to the target axes.
    title : str or None, default None
        Explicit axes title. Falls back to the snapshot time when omitted.
    norm : matplotlib.colors.Normalize or None, default None
        Explicit normalization object. If provided, it overrides the default
        linear/log normalization choice.
    """

    if slice_axis not in ["x", "y", "z"]:
        raise ValueError("slice_axis must be 'x', 'y', or 'z'")

    if not isinstance(slice_position, (int, float)):
        raise ValueError("slice_position must be a numeric value (int or float)")
    if not isinstance(show_grid, bool):
        raise ValueError("show_grid must be a boolean value (True or False)")

    if variable == "temp":
        required_vars = ["dens", "eint"]
        missing_vars = [var for var in required_vars if var not in my_data["var_names"]]
        if missing_vars:
            raise ValueError(
                f"To calculate temperature, need variables {required_vars}. Missing: {missing_vars}"
            )
    elif variable not in my_data["var_names"]:
        raise ValueError(
            f"Variable '{variable}' not found in data. Available variables: {my_data['var_names']}"
        )

    if show_quiver:
        vel_vars = ["velx", "vely", "velz"]
        missing_vel = [v for v in vel_vars if v not in my_data["var_names"]]
        if missing_vel:
            raise ValueError(
                f"Quiver requested but velocity variables missing: {missing_vel}"
            )

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    mb_list = list(range(my_data["n_mbs"]))
    mb_data = my_data["mb_data"]
    mb_geometry = my_data["mb_geometry"]

    mesh_bounds = np.array(
        [
            my_data["x1min"],
            my_data["x1max"],
            my_data["x2min"],
            my_data["x2max"],
            my_data["x3min"],
            my_data["x3max"],
        ]
    ).reshape(3, 2)

    axis_map = {"x": 0, "y": 1, "z": 2}
    slice_idx = axis_map[slice_axis]
    plot_axes = [i for i in range(3) if i != slice_idx]
    axis_names = ["x", "y", "z"]

    vel_names = ["velx", "vely", "velz"]
    vel_h_name = vel_names[plot_axes[0]]
    vel_v_name = vel_names[plot_axes[1]]

    valid_blocks = []
    data_min = float("inf")
    data_max = float("-inf")
    mb_logical = my_data.get("mb_logical")

    for mb_id in mb_list:
        bounds = np.array(mb_geometry[mb_id]).reshape(3, 2)

        if not (bounds[slice_idx, 0] <= slice_position < bounds[slice_idx, 1]):
            continue

        reference = variable if variable != "temp" else "dens"
        slice_count = mb_data[reference][mb_id].shape[slice_idx]
        if slice_count == 1:
            slice_block_idx = 0
        else:
            cell_width = (bounds[slice_idx, 1] - bounds[slice_idx, 0]) / slice_count
            slice_coord = bounds[slice_idx, 0] + (
                np.arange(slice_count, dtype=float) + 0.5
            ) * cell_width
            slice_block_idx = np.argmin(np.abs(slice_coord - slice_position))

        def extract_slice(arr):
            if slice_axis == "z":
                return arr[slice_block_idx, :, :]
            if slice_axis == "y":
                return arr[:, slice_block_idx, :]
            return arr[:, :, slice_block_idx]

        if variable == "temp":
            dens_slice = extract_slice(mb_data["dens"][mb_id])
            eint_slice = extract_slice(mb_data["eint"][mb_id])
            slice_data = unit_temp * eint_slice * (gamma - 1.0) / dens_slice
        else:
            slice_data = scale_factor * extract_slice(mb_data[variable][mb_id])

        if log_scale:
            with np.errstate(divide="ignore", invalid="ignore"):
                plot_data = np.where(slice_data > 0.0, slice_data, np.nan)
        else:
            plot_data = slice_data

        finite_data = plot_data[np.isfinite(plot_data)]
        if finite_data.size == 0:
            continue

        dx = (
            glw
            * show_grid
            * abs(mesh_bounds[plot_axes[0], 1] - mesh_bounds[plot_axes[0], 0])
            / 1000.0
        )
        dy = (
            glw
            * show_grid
            * abs(mesh_bounds[plot_axes[1], 1] - mesh_bounds[plot_axes[1], 0])
            / 1000.0
        )
        extent = [
            bounds[plot_axes[0], 0] + dx,
            bounds[plot_axes[0], 1] - dx,
            bounds[plot_axes[1], 0] + dy,
            bounds[plot_axes[1], 1] - dy,
        ]

        data_min = min(data_min, np.nanmin(plot_data))
        data_max = max(data_max, np.nanmax(plot_data))

        level = int(mb_logical[mb_id, 3]) if mb_logical is not None else 0
        cell_area = (
            (bounds[plot_axes[0], 1] - bounds[plot_axes[0], 0]) / plot_data.shape[1]
        ) * (
            (bounds[plot_axes[1], 1] - bounds[plot_axes[1], 0]) / plot_data.shape[0]
        )

        block_info = {
            "plot_data": plot_data,
            "extent": extent,
            "level": level,
            "cell_area": cell_area,
        }

        if show_quiver:
            vel_h_slice = extract_slice(mb_data[vel_h_name][mb_id])
            vel_v_slice = extract_slice(mb_data[vel_v_name][mb_id])

            ny, nx = vel_h_slice.shape
            xc = np.linspace(bounds[plot_axes[0], 0], bounds[plot_axes[0], 1], nx)
            yc = np.linspace(bounds[plot_axes[1], 0], bounds[plot_axes[1], 1], ny)
            xc = xc[::quiver_stride]
            yc = yc[::quiver_stride]
            Xc, Yc = np.meshgrid(xc, yc)

            block_info["quiver"] = {
                "X": Xc,
                "Y": Yc,
                "U": vel_h_slice[::quiver_stride, ::quiver_stride],
                "V": vel_v_slice[::quiver_stride, ::quiver_stride],
            }

        valid_blocks.append(block_info)

    if not valid_blocks:
        raise ValueError(
            f"No mesh blocks intersect slice {slice_axis}={slice_position} for '{variable}'."
        )

    valid_blocks.sort(key=lambda block: (block["level"], -block["cell_area"]))

    vmin = data_min if vmin is None else vmin
    vmax = data_max if vmax is None else vmax
    if norm is None and log_scale:
        norm = colors.LogNorm(vmin=vmin, vmax=vmax)

    im = None
    for block_info in valid_blocks:
        im = ax.imshow(
            block_info["plot_data"],
            extent=block_info["extent"],
            origin="lower",
            aspect="equal",
            alpha=1.0,
            cmap=colormap,
            interpolation="none",
            norm=norm,
            vmin=None if norm is not None else vmin,
            vmax=None if norm is not None else vmax,
        )

    if show_quiver:
        for block_info in valid_blocks:
            q = block_info["quiver"]
            U, V = q["U"], q["V"]
            if quiver_normalize:
                mag = np.sqrt(U**2 + V**2)
                mag[mag == 0] = 1.0
                U, V = U / mag, V / mag
            ax.quiver(
                q["X"],
                q["Y"],
                U,
                V,
                color=quiver_color,
                scale=quiver_scale,
                alpha=quiver_alpha,
                pivot="mid",
                width=0.002,
            )

    ax.set_xlim(mesh_bounds[plot_axes[0], 0], mesh_bounds[plot_axes[0], 1])
    ax.set_ylim(mesh_bounds[plot_axes[1], 0], mesh_bounds[plot_axes[1], 1])

    if variable == "temp":
        cbar_label = "Temperature [K]"
    else:
        cbar_label = variable

    if add_colorbar:
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(cbar_label)

    ax.set_xlabel(f"{axis_names[plot_axes[0]]} [kpc]")
    ax.set_ylabel(f"{axis_names[plot_axes[1]]} [kpc]")
    ax.set_title(title or f"{variable} at time={my_data['time'] / 10.0:.3f} Gyr")

    return fig, ax
