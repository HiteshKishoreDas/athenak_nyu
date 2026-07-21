#!/usr/bin/env python3

"""Plot the evolution of the volume-averaged temperature with time.

This script reuses the binary snapshot reader from
`plot_turb_drive_phase.py`, walks all `hydro_w` dumps in a simulation
directory, computes the mean temperature for each snapshot, and plots the
result as a time series.

Edit the top-level variables below to point at a different simulation
directory, output path, or snapshot glob.
"""

from __future__ import annotations

from pathlib import Path

import os

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", str(Path("/tmp") / "matplotlib"))

import matplotlib.pyplot as plt  # noqa: E402

from plot_turb_drive_phase import (  # noqa: E402
    ATHINPUT,
    GAMMA,
    SIM,
    iter_binary_blocks,
    parse_binary_meta,
    read_athinput_gamma,
)


SNAPSHOT_GLOB = "bin/Turb.hydro_w.*.bin"
OUTPUT = Path(SIM) / "temperature_evolution.png"

# Plot controls.
LINE_COLOR = "0.15"
LINE_WIDTH = 2.0
MARKER = "o"
MARKER_SIZE = 3.5

# If `temperature` is not present in the output file, temperature is
# reconstructed from the ideal-gas relation `T = (gamma - 1) * e_int / rho`
# in code units.


def resolve_snapshots(sim: Path, pattern: str) -> list[Path]:
    snapshots = sorted(sim.glob(pattern))
    if not snapshots:
        raise FileNotFoundError(f"no snapshots matched {sim / pattern}")
    return snapshots


def snapshot_time(path: Path) -> float:
    meta = parse_binary_meta(path)
    pheader = meta["pheader"]
    return float(pheader["time"])


def average_temperature(path: Path) -> float:
    meta = parse_binary_meta(path)
    var_names = set(meta["var_names"])
    gamma = GAMMA
    if ATHINPUT.is_file():
        gamma = read_athinput_gamma(ATHINPUT, gamma)
    temp_sum = 0.0
    cell_count = 0

    for block in iter_binary_blocks(path, meta):
        dens = block["dens"].astype(np.float64, copy=False).ravel()
        if "temperature" in var_names:
            temp = block["temperature"].astype(np.float64, copy=False).ravel()
        elif "press" in var_names:
            press = block["press"].astype(np.float64, copy=False).ravel()
            temp = press / dens
        else:
            eint = block["eint"].astype(np.float64, copy=False).ravel()
            temp = (gamma - 1.0) * eint / dens

        valid = np.isfinite(dens) & np.isfinite(temp) & (dens > 0.0) & (temp > 0.0)
        if not np.any(valid):
            continue

        temp_sum += float(temp[valid].sum())
        cell_count += int(valid.sum())

    if cell_count == 0:
        raise RuntimeError(f"no valid cells found in {path}")

    return temp_sum / cell_count


def main() -> None:
    gamma = GAMMA
    if ATHINPUT.is_file():
        gamma = read_athinput_gamma(ATHINPUT, gamma)

    sim = Path(SIM)
    snapshots = resolve_snapshots(sim, SNAPSHOT_GLOB)
    print(f"Found {len(snapshots)} snapshots in {sim / SNAPSHOT_GLOB}")

    times = []
    temps = []
    for index, snapshot in enumerate(snapshots, start=1):
        print(f"[{index}/{len(snapshots)}] processing {snapshot.name}", flush=True)
        times.append(snapshot_time(snapshot))
        temps.append(average_temperature(snapshot))
        print(
            f"[{index}/{len(snapshots)}] time={times[-1]:.6g} "
            f"<T>={temps[-1]:.6g}",
            flush=True,
        )

    times = np.asarray(times, dtype=np.float64)
    temps = np.asarray(temps, dtype=np.float64)

    order = np.argsort(times)
    times = times[order]
    temps = temps[order]

    fig, ax = plt.subplots(figsize=(7, 4.5), constrained_layout=True)
    ax.plot(times, temps, color=LINE_COLOR, linewidth=LINE_WIDTH, marker=MARKER, markersize=MARKER_SIZE)
    ax.set_xlabel("time")
    ax.set_ylabel(r"$\langle T \rangle$")
    ax.set_title("Volume-averaged temperature evolution")
    ax.grid(True, alpha=0.25)

    fig.savefig(OUTPUT, dpi=200)
    plt.show()


if __name__ == "__main__":
    main()
