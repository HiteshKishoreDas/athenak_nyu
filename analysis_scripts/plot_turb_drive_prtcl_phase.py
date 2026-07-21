#!/usr/bin/env python3

"""Plot a particle thermodynamic trajectory on the turb_drive phase planes."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib"))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402


RUN_DIR = Path(__file__).resolve().parent
ROOT_DIR = RUN_DIR.parent
SCRIPTS_DIR = ROOT_DIR / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from read_prtcl_thermo_history import read_history  # noqa: E402


SIM = RUN_DIR
ATHINPUT = SIM / "turb.athinput"
DEFAULT_HISTORY = SIM / "prtcl_thermo_history" / "Turb.prtcl_thermo_history.thp"
DEFAULT_OUTPUT = SIM / "prtcl_phase_trajectory.png"

GAMMA = 1.667
RHO0 = 1.0
PRS0 = 1.0

ISOBARIC_STYLE = {"color": "0.2", "linestyle": "--", "linewidth": 2.2, "zorder": 4}
ISOTHERMAL_STYLE = {"color": "0.35", "linestyle": ":", "linewidth": 1.8, "zorder": 4}
ADIABATIC_STYLE = {"color": "0.15", "linestyle": "-.", "linewidth": 1.6, "zorder": 4}
START_STYLE = {"marker": "o", "markersize": 7, "color": "white", "markeredgecolor": "black", "zorder": 6}
END_STYLE = {"marker": "s", "markersize": 7, "color": "white", "markeredgecolor": "black", "zorder": 6}


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


def pick_column(data: dict[str, np.ndarray], *names: str) -> np.ndarray:
    """Return the first matching column from a thermo-history table."""
    for name in names:
        if name in data:
            return np.asarray(data[name])
    raise KeyError(f"none of the requested columns are present: {', '.join(names)}")


def resolve_history(path: Path) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"particle thermo-history file not found: {path}")
    return path


def select_particle(
    data: dict[str, np.ndarray],
    tag: int | None = None,
    seed_id: int | None = None,
) -> dict[str, np.ndarray]:
    """Select one particle track from the full thermo-history table."""
    tags = np.asarray(data["tag"])
    seeds = np.asarray(data["seed_id"])
    mask = np.ones(tags.shape, dtype=bool)
    if tag is not None:
        mask &= tags == tag
    if seed_id is not None:
        mask &= seeds == seed_id

    if not np.any(mask):
        raise ValueError("no particle records matched the requested tag/seed_id")

    selected = {name: np.asarray(values)[mask] for name, values in data.items()}
    order = np.argsort(selected["time"], kind="mergesort")
    for name in selected:
        selected[name] = selected[name][order]
    return selected


def choose_default_particle(data: dict[str, np.ndarray]) -> tuple[int, int]:
    """Pick the longest-lived particle track when no explicit tag is supplied."""
    tags = np.asarray(data["tag"])
    seeds = np.asarray(data["seed_id"])
    pairs = np.stack([tags, seeds], axis=1)
    uniq, counts = np.unique(pairs, axis=0, return_counts=True)
    index = int(np.argmax(counts))
    tag, seed_id = uniq[index]
    return int(tag), int(seed_id)


def compute_trajectory_fields(data: dict[str, np.ndarray], gamma: float) -> dict[str, np.ndarray]:
    """Return log-scaled thermodynamic fields for plotting."""
    rho = pick_column(data, "density", "rho").astype(np.float64, copy=False)
    temp = pick_column(data, "temperature", "T").astype(np.float64, copy=False)

    if "pressure" in data:
        press = np.asarray(data["pressure"], dtype=np.float64)
    elif "press" in data:
        press = np.asarray(data["press"], dtype=np.float64)
    else:
        press = rho * temp

    if "specific_entropy" in data:
        entropy = np.asarray(data["specific_entropy"], dtype=np.float64)
    elif "entropy" in data:
        entropy = np.asarray(data["entropy"], dtype=np.float64)
    else:
        entropy = press / np.power(rho, gamma)

    time = np.asarray(data["time"], dtype=np.float64)
    valid = np.isfinite(time) & np.isfinite(rho) & np.isfinite(temp) & np.isfinite(press) & np.isfinite(entropy)
    valid &= (rho > 0.0) & (temp > 0.0) & (press > 0.0) & (entropy > 0.0)
    if not np.any(valid):
        raise RuntimeError("particle track has no valid thermodynamic samples")

    return {
        "time": time[valid],
        "log_rho": np.log10(rho[valid]),
        "log_temp": np.log10(temp[valid]),
        "log_press": np.log10(press[valid]),
        "log_entropy": np.log10(entropy[valid]),
    }


def trajectory_segments(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Build line segments for a colored trajectory."""
    points = np.column_stack([x, y])
    return np.stack([points[:-1], points[1:]], axis=1)


def add_colored_trajectory(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    time: np.ndarray,
    cmap: str,
    norm: Normalize,
    linewidth: float = 2.1,
) -> LineCollection | None:
    """Draw a trajectory as a time-colored line plus start/end markers."""
    if x.size < 2:
        ax.plot(x, y, color=plt.get_cmap(cmap)(0.75), linewidth=linewidth, zorder=5)
        if x.size == 1:
            ax.plot(x[0], y[0], **START_STYLE)
        return None

    segments = trajectory_segments(x, y)
    lc = LineCollection(
        segments,
        cmap=cmap,
        norm=norm,
        linewidth=linewidth,
        zorder=5,
    )
    lc.set_array(time[:-1])
    ax.add_collection(lc)
    ax.plot(x[0], y[0], **START_STYLE)
    ax.plot(x[-1], y[-1], **END_STYLE)
    return lc


def configure_reference_lines(ax_tp: plt.Axes, ax_trho: plt.Axes, ax_sp: plt.Axes, ax_st: plt.Axes, gamma: float, rho0: float, prs0: float) -> None:
    """Overlay the same thermodynamic guide lines used by the hydro phase plots."""
    t0 = prs0 / rho0
    log_t0 = np.log10(t0)
    log_rho0 = np.log10(rho0)
    log_prs0 = np.log10(prs0)

    tp_x = np.linspace(*ax_tp.get_xlim(), 256)
    tr_x = np.linspace(*ax_trho.get_xlim(), 256)
    sp_x = np.linspace(*ax_sp.get_xlim(), 256)
    st_x = np.linspace(*ax_st.get_xlim(), 256)

    ax_tp.axhline(log_prs0, label="isobaric", **ISOBARIC_STYLE)
    ax_tp.plot(tp_x, log_prs0 + tp_x * (gamma / (gamma - 1.0)), label="adiabatic", **ADIABATIC_STYLE)

    ax_trho.plot(tr_x, log_prs0 - tr_x, label="isobaric", **ISOBARIC_STYLE)
    ax_trho.plot(tr_x, log_rho0 + (1.0 / (gamma - 1.0)) * tr_x, label="adiabatic", **ADIABATIC_STYLE)

    ax_sp.plot(sp_x, gamma * log_t0 + (1.0 - gamma) * sp_x, label="isothermal", **ISOTHERMAL_STYLE)
    ax_st.plot(st_x, gamma * st_x + (1.0 - gamma) * log_prs0, label="isobaric", **ISOBARIC_STYLE)


def plot_particle_phase(
    traj: dict[str, np.ndarray],
    output: Path,
    gamma: float,
    rho0: float,
    prs0: float,
) -> None:
    """Render the four phase projections for a single particle trajectory."""
    time = traj["time"]
    log_rho = traj["log_rho"]
    log_temp = traj["log_temp"]
    log_press = traj["log_press"]
    log_entropy = traj["log_entropy"]

    tmin = float(time.min())
    tmax = float(time.max())
    if tmax > tmin:
        norm = Normalize(vmin=tmin, vmax=tmax)
    else:
        norm = Normalize(vmin=tmin - 0.5, vmax=tmax + 0.5)

    cmap = "viridis"
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.5), constrained_layout=True)
    fig.suptitle("Particle thermodynamic trajectory in phase space", fontsize=13)

    lc0 = add_colored_trajectory(axes[0, 0], log_temp, log_press, time, cmap, norm)
    axes[0, 0].set_xlabel(r"$\log_{10}(T)$")
    axes[0, 0].set_ylabel(r"$\log_{10}(P)$")
    axes[0, 0].set_title(r"$T$-$P$")

    lc1 = add_colored_trajectory(axes[0, 1], log_temp, log_rho, time, cmap, norm)
    axes[0, 1].set_xlabel(r"$\log_{10}(T)$")
    axes[0, 1].set_ylabel(r"$\log_{10}(\rho)$")
    axes[0, 1].set_title(r"$T$-$\rho$")

    lc2 = add_colored_trajectory(axes[1, 0], log_press, log_entropy, time, cmap, norm)
    axes[1, 0].set_xlabel(r"$\log_{10}(P)$")
    axes[1, 0].set_ylabel(r"$\log_{10}(P/\rho^\gamma)$")
    axes[1, 0].set_title(r"Entropy-$P$")

    lc3 = add_colored_trajectory(axes[1, 1], log_temp, log_entropy, time, cmap, norm)
    axes[1, 1].set_xlabel(r"$\log_{10}(T)$")
    axes[1, 1].set_ylabel(r"$\log_{10}(P/\rho^\gamma)$")
    axes[1, 1].set_title(r"Entropy-$T$")

    for ax in axes.flat:
        ax.grid(False)
        ax.set_axisbelow(True)

    configure_reference_lines(axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1], gamma, rho0, prs0)

    for ax in axes.flat:
        ax.legend(loc="best", frameon=False)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, shrink=0.95)
    cbar.set_label("time")

    if lc0 is None and lc1 is None and lc2 is None and lc3 is None:
        raise RuntimeError("failed to create trajectory artists")

    fig.savefig(output, dpi=200)
    plt.close(fig)




def main() -> None:
    global GAMMA
    global RHO0
    global PRS0

    hist_path = DEFAULT_HISTORY
    output =DEFAULT_OUTPUT
    tag = None #args.tag
    seed_id = None #args.seed_id

    if ATHINPUT.is_file():
        GAMMA = read_athinput_gamma(ATHINPUT, GAMMA)
        RHO0, PRS0 = read_athinput_reference_state(ATHINPUT)

    history_path = resolve_history(hist_path)
    data = read_history(history_path)

    if tag is None and seed_id is None:
        tag, seed_id = choose_default_particle(data)

    trajectory = select_particle(data, tag=tag, seed_id=seed_id)
    fields = compute_trajectory_fields(trajectory, GAMMA)

    print(
        f"plotting particle trajectory: tag={int(np.asarray(trajectory['tag'])[0])} "
        f"seed_id={int(np.asarray(trajectory['seed_id'])[0])} samples={len(fields['time'])}",
        flush=True,
    )
    plot_particle_phase(fields, output, GAMMA, RHO0, PRS0)
    print(f"wrote {output}", flush=True)


if __name__ == "__main__":
    main()
