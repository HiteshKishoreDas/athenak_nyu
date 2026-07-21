#!/usr/bin/env python3

"""Plot T-P, T-rho, entropy-P, and entropy-T phase diagrams from turb_drive hydro output.

The default path targets the binary `hydro_w` dumps produced by
`turb_drive/turb.athinput`. The script reads the Athena binary files directly
so it does not depend on `h5py`.

Edit the top-level variables below to point at a different simulation directory
and snapshot, adjust the binning, or change the output path.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, Iterable, Iterator, Tuple

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib"))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402


RUN_DIR = Path(__file__).resolve().parent
SIM = RUN_DIR
SNAPSHOT = "bin/Turb.hydro_w.00010.bin"
ATHINPUT = SIM / "turb.athinput"
OUTPUT = SIM / "phase_t_pres_t_rho.png"

# Plot controls.
NBINS = 128
CMAP = "magma"
# Axis limits are in log10 space.
X_LIM = None  # example: (-1.0, 1.0)
ENTROPY_LIM = None  # example: (-3.0, 1.0)
TP_Y_LIM = None  # example: (-4.0, 2.0)
TRHO_Y_LIM = None  # example: (-4.0, 1.0)
ISOBARIC_STYLE = {"color": "0.2", "linestyle": "--", "linewidth": 2.6, "zorder": 10}
ISOTHERMAL_STYLE = {"color": "0.35", "linestyle": ":", "linewidth": 2.0, "zorder": 10}
ADIABATIC_STYLE = {"color": "0.15", "linestyle": "-.", "linewidth": 1.8, "zorder": 10}
CONTOUR_STYLE = {"colors": "0.05", "linewidths": 0.8, "alpha": 0.8, "zorder": 9}
CONTOUR_TRAIL_STYLE = {"colors": "0.05", "linewidths": 1.0, "zorder": 8}
AVERAGE_MARKER_STYLE = {
    "color": "g",
    "marker": "+",
    "markersize": 12,
    "markeredgewidth": 2.2,
    "linestyle": "None",
    "zorder": 11,
}
AVERAGE_TRAIL_STYLE = {"color": "g", "linewidth": 1.4, "alpha": 0.7, "zorder": 10}
CONTOUR_TRAIL_LEVEL = 1.0e3
CONTOUR_TRAIL_ALPHA_MIN = 0.08
CONTOUR_TRAIL_ALPHA_MAX = 0.55

# If `press` is not present in the output file, pressure is reconstructed from
# the ideal-gas relation `P = (gamma - 1) * eint`.
GAMMA = 1.667
RHO0 = 1.0
PRS0 = 1.0


def read_athinput_value(path: Path, block: str, key: str, default: str | None = None) -> str | None:
    """Read a scalar value from an Athena-style input file."""
    current_block = None
    with path.open("r", encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            if line.startswith("<") and line.endswith(">"):
                current_block = line[1:-1].strip()
                continue
            if current_block != block or "=" not in line:
                continue
            lhs, rhs = line.split("=", 1)
            if lhs.strip() == key:
                return rhs.strip()
    return default


def read_athinput_gamma(path: Path, fallback: float) -> float:
    value = read_athinput_value(path, "hydro", "gamma")
    return float(value) if value is not None else fallback


def read_athinput_reference_state(path: Path) -> tuple[float, float]:
    rho0 = read_athinput_value(path, "problem", "rho0")
    prs0 = read_athinput_value(path, "problem", "prs0")
    return (
        float(rho0) if rho0 is not None else RHO0,
        float(prs0) if prs0 is not None else PRS0,
    )


def parse_binary_meta(path: Path) -> Dict[str, object]:
    """Read the Athena binary header and return metadata needed for streaming."""
    meta: Dict[str, object] = {}
    with path.open("rb") as fp:
        code_header = fp.readline().split()
        if len(code_header) < 1 or code_header[0] != b"Athena":
            raise TypeError(f"bad file format in {path}")
        version = code_header[-1].split(b"=")[-1]
        if version != b"1.1":
            raise TypeError(f"unsupported file format version {version.decode('utf-8')}")

        pheader_count = int(fp.readline().split(b"=")[-1])
        pheader: Dict[str, str] = {}
        for _ in range(pheader_count - 1):
            key, val = [x.strip() for x in fp.readline().decode("utf-8").split("=", 1)]
            pheader[key] = val

        nvars = int(fp.readline().split(b"=")[-1])
        var_names = [v.decode("utf-8") for v in fp.readline().split()[1:]]
        header_size = int(fp.readline().split(b"=")[-1])
        header = [line.decode("utf-8").split("#")[0].strip() for line in fp.read(header_size).split(b"\n")]
        header = [line for line in header if line]

    def get_from_header(block: str, key: str) -> str:
        want_block = block if block.startswith("<") else f"<{block}>"
        if want_block[-1] != ">":
            want_block += ">"
        current = "<none>"
        for line in header:
            if line.startswith("<"):
                current = line
                continue
            if current == want_block and "=" in line:
                lhs, rhs = line.split("=", 1)
                if lhs.strip() == key:
                    return rhs.strip()
        raise KeyError(f"no parameter called {want_block}/{key}")

    meta["path"] = path
    meta["pheader"] = pheader
    meta["var_names"] = var_names
    meta["nvars"] = nvars
    meta["nghost"] = int(get_from_header("<mesh>", "nghost"))
    meta["locsizebytes"] = int(pheader["size of location"])
    meta["varsizebytes"] = int(pheader["size of variable"])
    meta["data_offset"] = None
    return meta


def iter_binary_blocks(path: Path, meta: Dict[str, object]) -> Iterator[Dict[str, np.ndarray]]:
    """Yield each meshblock as a dictionary of variable arrays."""
    var_names = list(meta["var_names"])
    nvars = int(meta["nvars"])
    nghost = int(meta["nghost"])
    locsizebytes = int(meta["locsizebytes"])
    varsizebytes = int(meta["varsizebytes"])
    data_dtype = np.float64 if varsizebytes == 8 else np.float32
    geom_dtype = np.float64 if locsizebytes == 8 else np.float32

    with path.open("rb") as fp:
        fp.readline()
        pheader_count = int(fp.readline().split(b"=")[-1])
        for _ in range(pheader_count - 1):
            fp.readline()
        fp.readline()
        fp.readline()
        header_size = int(fp.readline().split(b"=")[-1])
        fp.read(header_size)

        while True:
            index_bytes = fp.read(24)
            if not index_bytes:
                break
            if len(index_bytes) != 24:
                raise EOFError(f"truncated meshblock index in {path}")

            mb_index = np.frombuffer(index_bytes, dtype=np.int32).astype(np.int64) - nghost
            nx1_out = int(mb_index[1] - mb_index[0] + 1)
            nx2_out = int(mb_index[3] - mb_index[2] + 1)
            nx3_out = int(mb_index[5] - mb_index[4] + 1)

            logical = np.frombuffer(fp.read(16), dtype=np.int32)
            geometry = np.frombuffer(fp.read(6 * locsizebytes), dtype=geom_dtype)
            raw = np.fromfile(fp, dtype=data_dtype, count=nx1_out * nx2_out * nx3_out * nvars)
            if raw.size != nx1_out * nx2_out * nx3_out * nvars:
                raise EOFError(f"truncated meshblock payload in {path}")

            data = raw.reshape(nvars, nx3_out, nx2_out, nx1_out)
            block = {name: data[i] for i, name in enumerate(var_names)}
            block["mb_index"] = mb_index
            block["mb_logical"] = logical
            block["mb_geometry"] = geometry
            yield block


def update_bounds(bounds: Dict[str, Tuple[float, float]], name: str, values: np.ndarray) -> None:
    """Track positive finite bounds for log-scaled histograms."""
    finite = values[np.isfinite(values) & (values > 0.0)]
    if finite.size == 0:
        return
    current = bounds.get(name)
    lo = float(finite.min())
    hi = float(finite.max())
    if current is None:
        bounds[name] = (lo, hi)
    else:
        bounds[name] = (min(current[0], lo), max(current[1], hi))


def compute_thermo_fields(
    block: Dict[str, np.ndarray],
    var_names: set[str],
    gamma: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return density, pressure, temperature, and entropy arrays for a block."""
    dens = block["dens"].astype(np.float64, copy=False)
    if "press" in var_names:
        press = block["press"].astype(np.float64, copy=False)
        temp = block["temperature"].astype(np.float64, copy=False) if "temperature" in var_names else press / dens
    else:
        eint = block["eint"].astype(np.float64, copy=False)
        press = (gamma - 1.0) * eint
        temp = press / dens
    entropy = press / np.power(dens, gamma)
    return dens, press, temp, entropy


def block_cell_volume(block: Dict[str, np.ndarray]) -> float:
    """Return the volume of a single active cell in a meshblock."""
    geometry = np.asarray(block["mb_geometry"], dtype=np.float64)
    if geometry.size != 6:
        raise ValueError("unexpected meshblock geometry payload")

    mb_index = np.asarray(block["mb_index"], dtype=np.int64)
    nx1 = int(mb_index[1] - mb_index[0] + 1)
    nx2 = int(mb_index[3] - mb_index[2] + 1)
    nx3 = int(mb_index[5] - mb_index[4] + 1)

    lengths = np.abs(geometry[1::2] - geometry[0::2])
    return float(np.prod(lengths) / (nx1 * nx2 * nx3))


def average_log_coordinates(hist: Dict[str, np.ndarray]) -> Dict[str, float]:
    """Return the current volume-weighted average coordinates in log space."""
    avg_log_t = np.log10(float(hist["avg_temperature"]))
    avg_log_p = np.log10(float(hist["avg_pressure"]))
    avg_log_rho = np.log10(float(hist["avg_density"]))
    avg_log_entropy = np.log10(float(hist["avg_entropy"]))
    return {
        "tp_x": avg_log_t,
        "tp_y": avg_log_p,
        "trho_x": avg_log_t,
        "trho_y": avg_log_rho,
        "sp_x": avg_log_p,
        "sp_y": avg_log_entropy,
        "st_x": avg_log_t,
        "st_y": avg_log_entropy,
    }


def contour_alpha(index: int, total: int) -> float:
    """Map a historical contour index to a fading alpha."""
    if total <= 0:
        return CONTOUR_TRAIL_ALPHA_MAX
    frac = (index + 1) / total
    return CONTOUR_TRAIL_ALPHA_MIN + frac * (CONTOUR_TRAIL_ALPHA_MAX - CONTOUR_TRAIL_ALPHA_MIN)


def gather_stats(paths: Iterable[Path], gamma: float) -> Dict[str, Tuple[float, float]]:
    """First pass: find histogram bounds for all requested quantities."""
    bounds: Dict[str, Tuple[float, float]] = {}
    for path in paths:
        meta = parse_binary_meta(path)
        var_names = set(meta["var_names"])
        for block in iter_binary_blocks(path, meta):
            dens, press, temp, entropy = compute_thermo_fields(block, var_names, gamma)

            update_bounds(bounds, "density", dens)
            update_bounds(bounds, "pressure", press)
            update_bounds(bounds, "temperature", temp)
            update_bounds(bounds, "entropy", entropy)
    return bounds


def accumulate_histograms(
    paths: Iterable[Path],
    gamma: float,
    bounds: Dict[str, Tuple[float, float]],
) -> Dict[str, np.ndarray]:
    """Second pass: stream cells into fixed 2D histograms."""
    tmin, tmax = bounds["temperature"]
    rmin, rmax = bounds["density"]
    pmin, pmax = bounds["pressure"]
    smin, smax = bounds["entropy"]

    t_edges = np.linspace(np.log10(tmin), np.log10(tmax), NBINS + 1)
    r_edges = np.linspace(np.log10(rmin), np.log10(rmax), NBINS + 1)
    p_edges = np.linspace(np.log10(pmin), np.log10(pmax), NBINS + 1)
    s_edges = np.linspace(np.log10(smin), np.log10(smax), NBINS + 1)

    t_pres = np.zeros((NBINS, NBINS), dtype=np.float64)
    t_rho = np.zeros((NBINS, NBINS), dtype=np.float64)
    s_pres = np.zeros((NBINS, NBINS), dtype=np.float64)
    s_temp = np.zeros((NBINS, NBINS), dtype=np.float64)
    weighted_volume = 0.0
    weighted_density = 0.0
    weighted_pressure = 0.0
    weighted_temperature = 0.0
    weighted_entropy = 0.0

    for path in paths:
        meta = parse_binary_meta(path)
        var_names = set(meta["var_names"])
        for block in iter_binary_blocks(path, meta):
            dens, press, temp, entropy = compute_thermo_fields(block, var_names, gamma)
            cell_volume = block_cell_volume(block)
            dens = dens.ravel()
            press = press.ravel()
            temp = temp.ravel()
            entropy = entropy.ravel()

            valid = np.isfinite(dens) & np.isfinite(press) & np.isfinite(temp) & np.isfinite(entropy)
            valid &= (dens > 0.0) & (press > 0.0) & (temp > 0.0)
            valid &= entropy > 0.0
            if not np.any(valid):
                continue

            log_temp = np.log10(temp[valid])
            log_pres = np.log10(press[valid])
            log_dens = np.log10(dens[valid])
            log_entropy = np.log10(entropy[valid])

            hist_tp, _, _ = np.histogram2d(log_temp, log_pres, bins=[t_edges, p_edges])
            hist_tr, _, _ = np.histogram2d(log_temp, log_dens, bins=[t_edges, r_edges])
            hist_sp, _, _ = np.histogram2d(log_entropy, log_pres, bins=[s_edges, p_edges])
            hist_st, _, _ = np.histogram2d(log_entropy, log_temp, bins=[s_edges, t_edges])
            t_pres += hist_tp
            t_rho += hist_tr
            s_pres += hist_sp
            s_temp += hist_st

            valid_volume = cell_volume * float(valid.sum())
            weighted_volume += valid_volume
            weighted_density += cell_volume * float(dens[valid].sum())
            weighted_pressure += cell_volume * float(press[valid].sum())
            weighted_temperature += cell_volume * float(temp[valid].sum())
            weighted_entropy += cell_volume * float(entropy[valid].sum())

    if not np.any(t_pres) and not np.any(t_rho):
        raise RuntimeError("no valid cells found for phase plot")
    if weighted_volume <= 0.0:
        raise RuntimeError("no valid cell volume found for phase plot")

    avg_density = weighted_density / weighted_volume
    avg_pressure = weighted_pressure / weighted_volume
    avg_temperature = weighted_temperature / weighted_volume
    avg_entropy = weighted_entropy / weighted_volume

    return {
        "t_edges": t_edges,
        "p_edges": p_edges,
        "r_edges": r_edges,
        "s_edges": s_edges,
        "t_pres": t_pres,
        "t_rho": t_rho,
        "s_pres": s_pres,
        "s_temp": s_temp,
        "avg_density": avg_density,
        "avg_pressure": avg_pressure,
        "avg_temperature": avg_temperature,
        "avg_entropy": avg_entropy,
    }


def resolve_input(sim: Path, snapshot: str | Path) -> Path:
    snap = Path(snapshot)
    path = snap if snap.is_absolute() else sim / snap
    if not path.is_file():
        raise FileNotFoundError(f"snapshot not found: {path}")
    return path


def plot_phase_diagrams(
    hist: Dict[str, np.ndarray],
    output: Path,
    show: bool = True,
    x_lim: tuple[float, float] | None = None,
    entropy_lim: tuple[float, float] | None = None,
    tp_y_lim: tuple[float, float] | None = None,
    trho_y_lim: tuple[float, float] | None = None,
    sp_y_lim: tuple[float, float] | None = None,
    st_y_lim: tuple[float, float] | None = None,
    norm_max: float | None = None,
    average_history: list[Dict[str, float]] | None = None,
    contour_history: list[Dict[str, np.ndarray]] | None = None,
) -> None:
    t_edges = hist["t_edges"]
    p_edges = hist["p_edges"]
    r_edges = hist["r_edges"]
    s_edges = hist["s_edges"]
    t_pres = np.ma.masked_less_equal(hist["t_pres"].T, 0.0)
    t_rho = np.ma.masked_less_equal(hist["t_rho"].T, 0.0)
    s_pres = np.ma.masked_less_equal(hist["s_pres"], 0.0)
    s_temp = np.ma.masked_less_equal(hist["s_temp"], 0.0)
    vmax = max(
        float(hist["t_pres"].max()),
        float(hist["t_rho"].max()),
        float(hist["s_pres"].max()),
        float(hist["s_temp"].max()),
    )
    if norm_max is None:
        norm_max = vmax
    norm = LogNorm(vmin=1.0, vmax=max(norm_max, 1.0))
    t0 = PRS0 / RHO0
    log_t0 = np.log10(t0)
    log_rho0 = np.log10(RHO0)
    log_prs0 = np.log10(PRS0)

    isobaric_label = "isobaric"
    isothermal_label = "isothermal"
    adiabatic_label = "adiabatic"
    t_centers = 0.5 * (t_edges[:-1] + t_edges[1:])
    p_centers = 0.5 * (p_edges[:-1] + p_edges[1:])
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    s_centers = 0.5 * (s_edges[:-1] + s_edges[1:])
    current_average = average_log_coordinates(hist)
    average_history = list(average_history or [])
    contour_history = list(contour_history or [])

    def contour_levels(values: np.ndarray) -> np.ndarray:
        vmax_local = float(np.max(values))
        if vmax_local < 1.0e3:
            return np.array([], dtype=np.float64)
        start = int(np.ceil(np.log10(1.0e3)))
        stop = int(np.floor(np.log10(vmax_local)))
        return np.array([10.0 ** n for n in range(start, stop + 1)], dtype=np.float64)

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5), constrained_layout=True)
    fig.suptitle("Turbulence phase diagrams from `hydro_w` output", fontsize=13)

    m0 = axes[0, 0].pcolormesh(t_edges, p_edges, t_pres, shading="auto", cmap=CMAP, norm=norm)
    levels_tp = contour_levels(hist["t_pres"])
    if levels_tp.size > 0:
        axes[0, 0].contour(t_centers, p_centers, hist["t_pres"].T, levels=levels_tp, **CONTOUR_STYLE)
    for index, past_hist in enumerate(contour_history):
        if float(np.max(past_hist["t_pres"])) >= CONTOUR_TRAIL_LEVEL:
            axes[0, 0].contour(
                t_centers,
                p_centers,
                past_hist["t_pres"].T,
                levels=[CONTOUR_TRAIL_LEVEL],
                alpha=contour_alpha(index, len(contour_history)),
                **CONTOUR_TRAIL_STYLE,
            )
    axes[0, 0].set_xlabel(r"$\log_{10}(T)$")
    axes[0, 0].set_ylabel(r"$\log_{10}(P)$")
    axes[0, 0].set_title(r"$T$-$P$")
    if x_lim is None:
        x_lim = X_LIM
    if entropy_lim is None:
        entropy_lim = ENTROPY_LIM
    if tp_y_lim is None:
        tp_y_lim = TP_Y_LIM
    if trho_y_lim is None:
        trho_y_lim = TRHO_Y_LIM
    if sp_y_lim is None:
        sp_y_lim = TP_Y_LIM
    if st_y_lim is None:
        st_y_lim = X_LIM
    if x_lim is not None:
        axes[0, 0].set_xlim(*x_lim)
        axes[0, 1].set_xlim(*x_lim)
    if tp_y_lim is not None:
        axes[0, 0].set_ylim(*tp_y_lim)
    if entropy_lim is not None:
        axes[1, 0].set_ylim(*entropy_lim)
        axes[1, 1].set_ylim(*entropy_lim)
    x0, x1 = axes[0, 0].get_xlim()
    tp_x = np.linspace(x0, x1, 256)
    axes[0, 0].axhline(log_prs0, label=isobaric_label, **ISOBARIC_STYLE)
    axes[0, 0].plot(tp_x, log_prs0 + tp_x * (GAMMA / (GAMMA - 1.0)), label=adiabatic_label, **ADIABATIC_STYLE)
    if len(average_history) > 1:
        axes[0, 0].plot(
            [point["tp_x"] for point in average_history],
            [point["tp_y"] for point in average_history],
            **AVERAGE_TRAIL_STYLE,
        )
    axes[0, 0].plot(current_average["tp_x"], current_average["tp_y"], **AVERAGE_MARKER_STYLE)
    axes[0, 0].legend(loc="best", frameon=False)

    m1 = axes[0, 1].pcolormesh(t_edges, r_edges, t_rho, shading="auto", cmap=CMAP, norm=norm)
    levels_trho = contour_levels(hist["t_rho"])
    if levels_trho.size > 0:
        axes[0, 1].contour(t_centers, r_centers, hist["t_rho"].T, levels=levels_trho, **CONTOUR_STYLE)
    for index, past_hist in enumerate(contour_history):
        if float(np.max(past_hist["t_rho"])) >= CONTOUR_TRAIL_LEVEL:
            axes[0, 1].contour(
                t_centers,
                r_centers,
                past_hist["t_rho"].T,
                levels=[CONTOUR_TRAIL_LEVEL],
                alpha=contour_alpha(index, len(contour_history)),
                **CONTOUR_TRAIL_STYLE,
            )
    axes[0, 1].set_xlabel(r"$\log_{10}(T)$")
    axes[0, 1].set_ylabel(r"$\log_{10}(\rho)$")
    axes[0, 1].set_title(r"$T$-$\rho$")
    if trho_y_lim is not None:
        axes[0, 1].set_ylim(*trho_y_lim)
    x0, x1 = axes[0, 1].get_xlim()
    tr_x = np.linspace(x0, x1, 256)
    axes[0, 1].plot(tr_x, log_prs0 - tr_x, label=isobaric_label, **ISOBARIC_STYLE)
    axes[0, 1].plot(tr_x, log_rho0 + (1 / (GAMMA - 1)) * tr_x, label=adiabatic_label, **ADIABATIC_STYLE)
    if len(average_history) > 1:
        axes[0, 1].plot(
            [point["trho_x"] for point in average_history],
            [point["trho_y"] for point in average_history],
            **AVERAGE_TRAIL_STYLE,
        )
    axes[0, 1].plot(current_average["trho_x"], current_average["trho_y"], **AVERAGE_MARKER_STYLE)
    axes[0, 1].legend(loc="best", frameon=False)

    m2 = axes[1, 0].pcolormesh(p_edges, s_edges, s_pres, shading="auto", cmap=CMAP, norm=norm)
    levels_sp = contour_levels(hist["s_pres"])
    if levels_sp.size > 0:
        axes[1, 0].contour(p_centers, s_centers, hist["s_pres"], levels=levels_sp, **CONTOUR_STYLE)
    for index, past_hist in enumerate(contour_history):
        if float(np.max(past_hist["s_pres"])) >= CONTOUR_TRAIL_LEVEL:
            axes[1, 0].contour(
                p_centers,
                s_centers,
                past_hist["s_pres"],
                levels=[CONTOUR_TRAIL_LEVEL],
                alpha=contour_alpha(index, len(contour_history)),
                **CONTOUR_TRAIL_STYLE,
            )
    axes[1, 0].set_title(r"Entropy-$P$")
    axes[1, 0].set_xlabel(r"$\log_{10}(P)$")
    axes[1, 0].set_ylabel(r"$\log_{10}(P/\rho^\gamma)$")
    p_x0, p_x1 = p_edges[0], p_edges[-1]
    sp_x = np.linspace(p_x0, p_x1, 256)
    axes[1, 0].plot(sp_x, GAMMA* log_t0 + (1.0 - GAMMA) * sp_x, label=isothermal_label, **ISOTHERMAL_STYLE)
    if len(average_history) > 1:
        axes[1, 0].plot(
            [point["sp_x"] for point in average_history],
            [point["sp_y"] for point in average_history],
            **AVERAGE_TRAIL_STYLE,
        )
    axes[1, 0].plot(current_average["sp_x"], current_average["sp_y"], **AVERAGE_MARKER_STYLE)
    axes[1, 0].legend(loc="best", frameon=False)

    m3 = axes[1, 1].pcolormesh(t_edges, s_edges, s_temp, shading="auto", cmap=CMAP, norm=norm)
    levels_st = contour_levels(hist["s_temp"])
    if levels_st.size > 0:
        axes[1, 1].contour(t_centers, s_centers, hist["s_temp"], levels=levels_st, **CONTOUR_STYLE)
    for index, past_hist in enumerate(contour_history):
        if float(np.max(past_hist["s_temp"])) >= CONTOUR_TRAIL_LEVEL:
            axes[1, 1].contour(
                t_centers,
                s_centers,
                past_hist["s_temp"],
                levels=[CONTOUR_TRAIL_LEVEL],
                alpha=contour_alpha(index, len(contour_history)),
                **CONTOUR_TRAIL_STYLE,
            )
    axes[1, 1].set_xlabel(r"$\log_{10}(T)$")
    axes[1, 1].set_ylabel(r"$\log_{10}(P/\rho^\gamma)$")
    axes[1, 1].set_title(r"Entropy-$T$")
    t_x0, t_x1 = t_edges[0], t_edges[-1]
    st_x = np.linspace(t_x0, t_x1, 256)
    axes[1, 1].plot(st_x, GAMMA* st_x + (1.0 - GAMMA) * log_prs0, label=isobaric_label, **ISOBARIC_STYLE)
    if len(average_history) > 1:
        axes[1, 1].plot(
            [point["st_x"] for point in average_history],
            [point["st_y"] for point in average_history],
            **AVERAGE_TRAIL_STYLE,
        )
    axes[1, 1].plot(current_average["st_x"], current_average["st_y"], **AVERAGE_MARKER_STYLE)
    axes[1, 1].legend(loc="best", frameon=False)

    cbar = fig.colorbar(m1, ax=axes, shrink=0.95)
    cbar.set_label("cell count")

    for ax in axes.flat:
        ax.grid(False)

    fig.savefig(output, dpi=200)
    if show:
        plt.show()
    else:
        plt.close(fig)


def main() -> None:
    global GAMMA
    global RHO0
    global PRS0
    sim = Path(SIM)
    if ATHINPUT.is_file():
        GAMMA = read_athinput_gamma(ATHINPUT, GAMMA)
        RHO0, PRS0 = read_athinput_reference_state(ATHINPUT)

    input_path = resolve_input(sim, SNAPSHOT)
    bounds = gather_stats([input_path], GAMMA)

    hist = accumulate_histograms([input_path], GAMMA, bounds)
    plot_phase_diagrams(hist, OUTPUT)


if __name__ == "__main__":
    main()
