import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.environ.setdefault("XDG_CACHE_HOME", f"/tmp/xdg-cache-{os.getuid()}")

import matplotlib.pyplot as plt
import numpy as np


script_dir = Path(__file__).resolve().parent
repo_root = script_dir.parent
sys.path.append(str(repo_root / "vis" / "python"))

import bin_convert_new as bc


RUN_PREFIX = "Turb.full_hydro_w"
GAMMA = 5.0 / 3.0
GAMMA_MINUS_1 = GAMMA - 1.0
MU = 0.6
FRAME_START = 0
FRAME_STOP = 30
FRAME_STEP = 1
OUTPUT_PATH = None

LENGTH_CGS = 3.0856775809623245e21
TIME_CGS = 3.15576e13
ATOMIC_MASS_UNIT_CGS = 1.67262192369e-24
KBOLTZ_CGS = 1.3806488e-16

TEMP_SPLIT_K = 787.0
TEMPERATURE_FLOOR = 1.0e-30
DENSITY_FLOOR = 1.0e-30


def configure_style():
    try:
        import own_package

        style_lib = Path(own_package.__file__).resolve().parent / "plot" / "style_lib"
        plt.style.use(
            [
                str(style_lib / "bright_pallette.mplstyle"),
                str(style_lib / "plot_style.mplstyle"),
                str(style_lib / "text.mplstyle"),
            ]
        )
    except Exception:
        pass


def get_temperature_unit():
    velocity_cgs = LENGTH_CGS / TIME_CGS
    return velocity_cgs * velocity_cgs * MU * ATOMIC_MASS_UNIT_CGS / KBOLTZ_CGS


def compute_temperature(out_dict):
    density = np.clip(np.asarray(out_dict["dens"], dtype=float), DENSITY_FLOOR, None)
    eint = np.asarray(out_dict["eint"], dtype=float)
    return GAMMA_MINUS_1 * eint / density * get_temperature_unit()


def compute_speed_squared(out_dict):
    velx = np.asarray(out_dict["velx"], dtype=float)
    vely = np.asarray(out_dict["vely"], dtype=float)
    velz = np.asarray(out_dict["velz"], dtype=float)
    return velx * velx + vely * vely + velz * velz


def compute_sound_speed_squared(out_dict):
    density = np.clip(np.asarray(out_dict["dens"], dtype=float), DENSITY_FLOOR, None)
    eint = np.asarray(out_dict["eint"], dtype=float)
    cs2 = GAMMA * GAMMA_MINUS_1 * eint / density
    return np.clip(cs2, 0.0, None)


def get_phase_stats(mask, weights, speed_sq, sound_speed_sq, temperature):
    cell_count = int(np.count_nonzero(mask))
    if cell_count == 0:
        return {
            "cell_count": 0,
            "mass": 0.0,
            "temperature_mean": np.nan,
            "velocity_rms": np.nan,
            "sound_speed_rms": np.nan,
            "ratio": np.nan,
        }

    masked_weights = weights[mask]
    mass = float(np.sum(masked_weights))
    if mass <= 0.0:
        return {
            "cell_count": cell_count,
            "mass": mass,
            "temperature_mean": np.nan,
            "velocity_rms": np.nan,
            "sound_speed_rms": np.nan,
            "ratio": np.nan,
        }

    velocity_rms = np.sqrt(np.sum(masked_weights * speed_sq[mask]) / mass)
    sound_speed_rms = np.sqrt(np.sum(masked_weights * sound_speed_sq[mask]) / mass)
    ratio = velocity_rms / sound_speed_rms if sound_speed_rms > 0.0 else np.nan

    return {
        "cell_count": cell_count,
        "mass": mass,
        "temperature_mean": np.sum(masked_weights * temperature[mask]) / mass,
        "velocity_rms": float(velocity_rms),
        "sound_speed_rms": float(sound_speed_rms),
        "ratio": float(ratio),
    }


def compute_phase_summary(out_dict, temperature_threshold):
    temperature = np.clip(compute_temperature(out_dict), TEMPERATURE_FLOOR, None)
    weights = np.clip(np.asarray(out_dict["dens"], dtype=float), 0.0, None)
    speed_sq = compute_speed_squared(out_dict)
    sound_speed_sq = compute_sound_speed_squared(out_dict)

    cold_mask = temperature < temperature_threshold
    hot_mask = ~cold_mask

    cold_stats = get_phase_stats(
        cold_mask, weights, speed_sq, sound_speed_sq, temperature
    )
    hot_stats = get_phase_stats(
        hot_mask, weights, speed_sq, sound_speed_sq, temperature
    )

    total_mass = float(np.sum(weights))
    total_cells = int(temperature.size)

    cold_stats["cell_fraction"] = cold_stats["cell_count"] / total_cells
    hot_stats["cell_fraction"] = hot_stats["cell_count"] / total_cells

    if total_mass > 0.0:
        cold_stats["mass_fraction"] = cold_stats["mass"] / total_mass
        hot_stats["mass_fraction"] = hot_stats["mass"] / total_mass
    else:
        cold_stats["mass_fraction"] = np.nan
        hot_stats["mass_fraction"] = np.nan

    return {
        "below": cold_stats,
        "above": hot_stats,
    }


def print_summary(frame_index, out_dict, temperature_threshold, summary):
    print(
        f"Frame {frame_index:05d} | cycle = {out_dict['NumCycles']} | "
        f"time = {out_dict['Time']:.6f}"
    )
    print(f"Temperature split: {temperature_threshold:.3f} K")

    for label, stats in (
        (f"T < {temperature_threshold:.0f} K", summary["below"]),
        (f"T >= {temperature_threshold:.0f} K", summary["above"]),
    ):
        print(label)
        print(
            f"  cells={stats['cell_count']} "
            f"({stats['cell_fraction']:.3%}), "
            f"mass={stats['mass']:.6e} "
            f"({stats['mass_fraction']:.3%})"
        )
        print(
            f"  T_mean={stats['temperature_mean']:.6e} K, "
            f"v_rms={stats['velocity_rms']:.6e}, "
            f"cs_rms={stats['sound_speed_rms']:.6e}, "
            f"v_rms/cs_rms={stats['ratio']:.6e}"
        )


def make_history_figure(times, summary_history, temperature_threshold):
    phase_styles = {
        "below": {
            "label": f"T < {temperature_threshold:.0f} K",
            "color": "#4c78a8",
        },
        "above": {
            "label": f"T >= {temperature_threshold:.0f} K",
            "color": "#f58518",
        },
    }
    quantity_specs = (
        ("velocity_rms", r"$v_\mathrm{rms}$"),
        ("sound_speed_rms", r"$c_{s,\mathrm{rms}}$"),
        ("ratio", r"$v_\mathrm{rms} / c_{s,\mathrm{rms}}$"),
    )

    fig, axes = plt.subplots(
        nrows=len(quantity_specs),
        ncols=1,
        figsize=(9.5, 10.0),
        sharex=True,
    )

    for ax, (key, ylabel) in zip(axes, quantity_specs):
        for phase_name, style in phase_styles.items():
            values = np.asarray(
                [summary[key] for summary in summary_history[phase_name]],
                dtype=float,
            )
            ax.plot(
                times,
                values,
                lw=2.0,
                color=style["color"],
                label=style["label"],
            )

        ax.set_ylabel(ylabel)
        ax.set_title(ylabel)
        ax.legend(loc="best")
        ax.tick_params(which="both", top=True, right=True)

    axes[-1].set_xlabel("Time")
    fig.suptitle(
        f"Phase RMS history | temperature split = {temperature_threshold:.0f} K"
    )
    fig.tight_layout(rect=(0.0, 0.0, 0.97, 0.96))
    return fig


def get_output_path(frame_start, frame_stop):
    if OUTPUT_PATH is not None:
        return Path(OUTPUT_PATH)

    save_dir = script_dir / "save" / "phase"
    save_dir.mkdir(parents=True, exist_ok=True)
    return save_dir / f"phase_history_{frame_start:05d}_{frame_stop - 1:05d}.png"


def main():
    configure_style()
    times = []
    processed_frames = []
    summary_history = {"below": [], "above": []}

    for frame_index in range(FRAME_START, FRAME_STOP, FRAME_STEP):
        bin_path = script_dir / "bin" / f"{RUN_PREFIX}.{frame_index:05d}.bin"
        if not bin_path.exists():
            print(f"Stopping at missing snapshot: {bin_path}")
            break

        print(f"Reading {bin_path}")
        out_dict = bc.read_binary_as_athdf(str(bin_path))

        summary = compute_phase_summary(out_dict, TEMP_SPLIT_K)
        print_summary(frame_index, out_dict, TEMP_SPLIT_K, summary)

        times.append(float(out_dict["Time"]))
        processed_frames.append(frame_index)
        summary_history["below"].append(summary["below"])
        summary_history["above"].append(summary["above"])

    if not times:
        raise FileNotFoundError(
            "No snapshots found for the configured frame range in turb_TI_SN/bin."
        )

    fig = make_history_figure(
        np.asarray(times, dtype=float), summary_history, TEMP_SPLIT_K
    )
    output_path = get_output_path(processed_frames[0], processed_frames[-1] + 1)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()
