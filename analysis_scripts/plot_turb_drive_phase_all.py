#!/usr/bin/env python3

"""Create phase plots for every turb_drive hydro snapshot.

This reuses the streaming accumulation and plotting helpers from
`plot_turb_drive_phase.py`, but loops over all matching snapshots and writes one
image per dump.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from plot_turb_drive_phase import (
    ATHINPUT,
    GAMMA,
    SIM,
    accumulate_histograms,
    gather_stats,
    plot_phase_diagrams,
    read_athinput_gamma,
)


SNAPSHOT_GLOB = "bin/Turb.hydro_w.*.bin"
OUTPUT_DIR = Path(SIM) / "phase_plots"


def resolve_snapshots(sim: Path, pattern: str) -> list[Path]:
    snapshots = sorted(sim.glob(pattern))
    if not snapshots:
        raise FileNotFoundError(f"no snapshots matched {sim / pattern}")
    return snapshots


def main() -> None:
    gamma = GAMMA
    if ATHINPUT.is_file():
        gamma = read_athinput_gamma(ATHINPUT, gamma)

    sim = Path(SIM)
    snapshots = resolve_snapshots(sim, SNAPSHOT_GLOB)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Found {len(snapshots)} snapshots in {sim / SNAPSHOT_GLOB}")

    shared_bounds = gather_stats(snapshots, gamma)
    x_lim = tuple(float(v) for v in np.log10(shared_bounds["temperature"]))
    entropy_lim = tuple(float(v) for v in np.log10(shared_bounds["entropy"]))
    tp_y_lim = tuple(float(v) for v in np.log10(shared_bounds["pressure"]))
    trho_y_lim = tuple(float(v) for v in np.log10(shared_bounds["density"]))

    histograms = []
    average_history = []
    norm_max = 1.0

    for index, snapshot in enumerate(snapshots, start=1):
        print(f"[{index}/{len(snapshots)}] processing {snapshot.name}", flush=True)
        hist = accumulate_histograms([snapshot], gamma, shared_bounds)
        average_history.append(
            {
                "tp_x": float(np.log10(float(hist["avg_temperature"]))),
                "tp_y": float(np.log10(float(hist["avg_pressure"]))),
                "trho_x": float(np.log10(float(hist["avg_temperature"]))),
                "trho_y": float(np.log10(float(hist["avg_density"]))),
                "sp_x": float(np.log10(float(hist["avg_pressure"]))),
                "sp_y": float(np.log10(float(hist["avg_entropy"]))),
                "st_x": float(np.log10(float(hist["avg_temperature"]))),
                "st_y": float(np.log10(float(hist["avg_entropy"]))),
            }
        )
        norm_max = max(norm_max, float(hist["t_pres"].max()), float(hist["t_rho"].max()))
        histograms.append((snapshot, hist))

    for index, (snapshot, hist) in enumerate(histograms, start=1):
        output = OUTPUT_DIR / f"{snapshot.stem}.png"
        plot_phase_diagrams(
            hist,
            output,
            show=False,
            x_lim=x_lim,
            entropy_lim=entropy_lim,
            tp_y_lim=tp_y_lim,
            trho_y_lim=trho_y_lim,
            sp_y_lim=tp_y_lim,
            st_y_lim=x_lim,
            norm_max=1e5,
            average_history=average_history[:index],
            contour_history=[item[1] for item in histograms[: index - 1]],
        )
        print(f"[{index}/{len(snapshots)}] wrote {output.name}", flush=True)


if __name__ == "__main__":
    main()
