#!/usr/bin/env python3
"""Paper Section 4 statistics from a driven, periodic SGS campaign.

Temporal scatter and turnover-time block means are descriptive statistics,
not confidence intervals: eight integral times need not be independent.
"""

import argparse
import json
from pathlib import Path

import numpy as np


PAPER = {
    1.0: {
        64: (5.38, 1.107, 0.053, 0.0726, 0.0055, 1.236, 0.141),
        128: (5.50, 0.787, 0.030, 0.0407, 0.0025, 1.230, 0.113),
        256: (5.55, 0.578, 0.022, 0.0236, 0.0013, 1.213, 0.098),
        512: (5.52, 0.412, 0.012, 0.0138, 0.0012, 1.219, 0.159),
    },
    0.0: {
        64: (5.29, 1.353, 0.049, 0.0235, 0.0044, 0.253, 0.072),
        128: (5.43, 1.040, 0.041, 0.0148, 0.0018, 0.286, 0.047),
        256: (5.57, 0.767, 0.029, 0.0086, 0.0013, 0.293, 0.061),
        512: (5.86, 0.528, 0.023, 0.0048, 0.0008, 0.292, 0.067),
    },
}


def window(t, y, start, stop):
    """Clip a piecewise-linear time series, including interpolated endpoints."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    if (
        len(t) < 2
        or len(y) != len(t)
        or not np.all(np.isfinite(t))
        or not np.all(np.isfinite(y))
        or np.any(np.diff(t) <= 0)
    ):
        raise ValueError("need finite samples with strictly increasing time")
    tolerance = 1e-10 * max(1.0, abs(stop))
    if stop <= start or start < t[0] - tolerance or stop > t[-1] + tolerance:
        raise ValueError("requested statistical window is not covered")
    inside = (t > start) & (t < stop)
    return (
        np.r_[start, t[inside], stop],
        np.r_[np.interp(start, t, y), y[inside], np.interp(stop, t, y)],
    )


def statistics(t, y, start, stop):
    t, y = window(t, y, start, stop)
    dt = np.diff(t)
    duration = stop - start
    mean = float(np.sum(dt * (y[:-1] + y[1:]) / 2) / duration)
    # Exact second moment of the same linear interpolant, centered for stability.
    a, b = y[:-1] - mean, y[1:] - mean
    variance = float(np.sum(dt * (a * a + a * b + b * b) / 3) / duration)
    mid = 0.5 * (start + stop)
    slope = float(
        np.sum(
            dt
            * (
                (t[:-1] - mid) * (2 * y[:-1] + y[1:])
                + (t[1:] - mid) * (y[:-1] + 2 * y[1:])
            )
            / 6
        )
        / (duration**3 / 12)
    )
    return dict(
        mean=mean,
        temporal_std=float(np.sqrt(max(0.0, variance))),
        start=float(start),
        stop=float(stop),
        slope=slope,
        fractional_linear_drift=slope * duration / mean if mean else None,
    )


def counter_rate(t, counter, start, stop):
    t, counter = window(t, counter, start, stop)
    rates = np.diff(counter) / np.diff(t)
    weights = np.diff(t) / (stop - start)
    mean = float((counter[-1] - counter[0]) / (stop - start))
    return dict(
        mean=mean,
        temporal_std=float(np.sqrt(np.sum(weights * (rates - mean) ** 2))),
        start=float(start),
        stop=float(stop),
    )


def summarize(data):
    p, samples, history = data["protocol"], data["samples"], data["history"]
    T, L, V = p["T"], p["L"], p["V"]
    if not samples or not history:
        raise ValueError("no native samples or histories")
    required = ("represented-E", "sgs-production-integral",
                "sgs-dissipation-integral", "force-work")
    if any(key not in row for row in history for key in required):
        raise ValueError("required applied SGS/forcing history counters missing")
    ts = np.array([s["time"] for s in samples])
    th = np.array([h["time"] for h in history])
    rho0 = float(p.get("rho0", samples[0]["rho_mean"]))
    start, stop = 2 * T, 10 * T
    covered = ts[0] <= start and ts[-1] >= stop - 1e-10 * max(1.0, stop)
    covered = covered and th[0] <= start and th[-1] >= stop - 1e-10 * max(1.0, stop)
    scales = dict(
        M_rms=1.0,
        M_sgs=1.0,
        K=1 / (rho0 * V * V),
        epsilon=L / V**3,
        rho_epsilon=L / (rho0 * V**3),
    )
    result = dict(
        protocol=p,
        complete_paper_window=bool(covered),
        numerical_status=data.get("numerical_status", "unreported"),
        requested_window=[start, stop],
        available_window=[float(ts[0]), float(ts[-1])],
        definitions={
            "K": "volume mean K / (rho0 V^2)",
            "epsilon": "volume mean specific epsilon * L / V^3",
            "rho_epsilon": "volume mean rho epsilon * L / (rho0 V^3)",
            "paper_dissipation": "Table 4 labels epsilon; nearby text says rho epsilon",
            "scatter": "time-weighted temporal standard deviation, not SEM",
            "blocks": "eight 1T block means; independence is not assumed",
            "interpolation": "linear between stored diagnostic samples",
            "heating": "corrected heating unavailable; pressure-work absent",
        },
        paper_table4=PAPER.get(float(p.get("sol_fraction", 1.0)), {}).get(p["n"]),
    )
    work = np.array([h["force-work"] for h in history])
    represented = np.array([h.get("represented-E", h["tot-E"]) for h in history])
    residual = represented - represented[0] - (work - work[0])
    work_scale = max(float(np.max(np.abs(work - work[0]))), np.finfo(float).tiny)
    result["energy_residual_over_applied_work"] = float(
        np.max(np.abs(residual)) / work_scale
    )
    result["energy_residual_over_initial_represented"] = float(
        np.max(np.abs(residual)) / max(abs(represented[0]), np.finfo(float).tiny)
    )
    if not covered:
        result["qualification"] = "Incomplete 2T-10T window; pilot only"
        return result
    result["statistics"] = {}
    result["blocks"] = {}
    result["early_late"] = {}
    for key, scale in scales.items():
        values = [s[key] * scale for s in samples]
        result["statistics"][key] = statistics(ts, values, start, stop)
        result["blocks"][key] = [
            statistics(ts, values, i * T, (i + 1) * T)["mean"] for i in range(2, 10)
        ]
        result["early_late"][key] = [
            statistics(ts, values, start, 6 * T)["mean"],
            statistics(ts, values, 6 * T, stop)["mean"],
        ]
    result["applied_rates"] = {}
    for key in ("sgs-production-integral", "sgs-dissipation-integral", "force-work"):
        vals = [h.get(key, 0.0) * L / (rho0 * V**3) for h in history]
        result["applied_rates"][key] = counter_rate(th, vals, start, stop)
    prod = result["applied_rates"]["sgs-production-integral"]["mean"]
    diss = result["applied_rates"]["sgs-dissipation-integral"]["mean"]
    result["production_dissipation_relative_imbalance"] = (
        abs(prod - diss) / diss if diss else None
    )
    result["qualification"] = (
        "Full sampling window; stationarity must be assessed from "
        "trends and blocks, and is not implied by elapsed time"
    )
    if data.get("numerical_status") != "pass":
        result["qualification"] += "; NUMERICAL CHECKS FAILED: unqualified"
    return result


def write_report(output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    output = Path(output)
    data = json.loads((output / "diagnostics.json").read_text())
    result = summarize(data)
    (output / "statistics.json").write_text(json.dumps(result, indent=2) + "\n")
    p, samples, hist = data["protocol"], data["samples"], data["history"]
    T, L, V = p["T"], p["L"], p["V"]
    rho0 = float(p.get("rho0", samples[0]["rho_mean"]))
    t = np.array([s["time"] / T for s in samples])
    ht = np.array([h["time"] / T for h in hist])

    def sample(key):
        return np.array([s[key] for s in samples])

    def history(key):
        return np.array([h.get(key, 0.0) for h in hist])

    target = result["paper_table4"]
    footer = (
        f"N={p['n']} cubed, gamma={p['gamma']}, T={T:g}, "
        f"a_rms={p['amplitude']:g}, seed={p['seed']}; "
    ) + (
        "2T-10T sampled; trends still require assessment"
        if result["complete_paper_window"]
        else "INCOMPLETE 2T-10T WINDOW: PILOT ONLY"
    )
    footer += f"; numerical checks: {data.get('numerical_status', 'unknown')}"
    pages = []
    with PdfPages(output / "report.pdf") as pdf:

        def finish(fig, axes, title, time_axes=True):
            for ax in np.asarray(axes).ravel():
                ax.grid(alpha=0.22)
                if time_axes:
                    ax.set_xlabel("time / T")
                    ax.axvline(2, color=".5", linestyle=":", linewidth=1)
                    ax.axvline(10, color=".5", linestyle=":", linewidth=1)
                if ax.get_legend_handles_labels()[0]:
                    ax.legend(fontsize=8)
            fig.suptitle(title, fontsize=14)
            fig.text(0.5, 0.012, footer, ha="center", fontsize=8)
            fig.tight_layout(rect=(0, 0.04, 1, 0.95))
            path = output / f"page{len(pages) + 1}.png"
            fig.savefig(path, dpi=150)
            pdf.savefig(fig)
            plt.close(fig)
            pages.append(path.name)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for ax, key, scale, label, band in zip(
            axes.flat,
            ("M_rms", "M_sgs", "K", "epsilon"),
            (1, 1, 1 / (rho0 * V**2), L / V**3),
            (
                "Resolved volume RMS Mach",
                "Volume RMS SGS Mach",
                "Mean K / (rho0 V^2)",
                "Volume mean epsilon * L / V^3",
            ),
            ((0, None), (1, 2), (3, 4), (5, 6)),
        ):
            ax.plot(
                t,
                sample(key) * scale,
                label="volume epsilon" if key == "epsilon" else "AthenaK",
            )
            if key == "epsilon":
                ax.plot(
                    t,
                    sample("rho_epsilon") * L / (rho0 * V**3),
                    label="volume rho epsilon / rho0",
                    color="tab:purple",
                )
            ax.set_ylabel(label)
            if target:
                mean = target[band[0]]
                ax.axhline(
                    mean, color="tab:orange", linestyle="--", label="paper Table 4"
                )
                if band[1] is not None:
                    sd = target[band[1]]
                    ax.axhspan(mean - sd, mean + sd, color="tab:orange", alpha=0.15)
        finish(
            fig,
            axes,
            "Driven 3D pseudo-isothermal turbulence: paper statistics\n"
            "Published targets: forcing differs; paper epsilon convention is ambiguous",
        )

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        rates = {}
        midpoint = 0.5 * (ht[:-1] + ht[1:])
        for key, label in [
            ("sgs-production-integral", "applied production"),
            ("sgs-dissipation-integral", "applied rho epsilon"),
            ("force-work", "applied forcing"),
        ]:
            rates[key] = np.diff(history(key)) / np.diff(ht) / T * L / (rho0 * V**3)
            axes[0, 0].plot(midpoint, rates[key], label=label)
        axes[0, 0].plot(
            t,
            sample("rho_epsilon") * L / (rho0 * V**3),
            ":",
            label="snapshot rho epsilon",
        )
        axes[0, 0].set_ylabel("Energy-density rates * L / (rho0 V^3)")
        for key in ("kinetic", "thermal", "represented"):
            values = sample(key)
            axes[0, 1].plot(t, values - values[0], label=key)
        axes[0, 1].set_ylabel("Inventory change from initial state")
        represented = (
            history("represented-E") if "represented-E" in hist[0] else history("tot-E")
        )
        work = history("force-work")
        residual = represented - represented[0] - (work - work[0])
        axes[1, 0].plot(ht, residual, label="delta(E+K) - applied work")
        axes[1, 0].set_ylabel("Represented-energy residual")
        for key, label in [
            ("sgs-stress-limited-integral", "stress suppressed"),
            ("sgs-diffusion-limited-integral", "diffusion suppressed"),
            ("sgs-fallback-stages", "fallback stages"),
        ]:
            axes[1, 1].plot(ht, history(key), label=label)
        axes[1, 1].set_ylabel("Cumulative guard diagnostics (different units)")
        finish(fig, axes, "Applied transfers, energy accounting, and numerical guards")

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        for ax, key, label in zip(
            axes.flat[:3],
            ("K", "M_sgs", "epsilon"),
            ("K / (rho0 V^2)", "SGS Mach", "epsilon * L / V^3"),
        ):
            if result["complete_paper_window"]:
                ax.plot(np.arange(2.5, 10), result["blocks"][key], "o-")
                ax.axhline(
                    result["statistics"][key]["mean"], linestyle="--", color=".4"
                )
            else:
                ax.text(
                    0.5,
                    0.5,
                    "Full 2T-10T window unavailable",
                    ha="center",
                    transform=ax.transAxes,
                )
            ax.set_xlabel("1T block midpoint / T")
            ax.set_ylabel(label)
        cov = sample("velocity_covariance")
        trace = np.trace(cov, axis1=1, axis2=2)
        for d in range(3):
            axes[1, 1].plot(
                t,
                np.divide(
                    cov[:, d, d], trace, out=np.zeros_like(trace), where=trace > 0
                ),
                label=f"component {d + 1}",
            )
        axes[1, 1].axhline(1 / 3, linestyle="--", color=".4")
        axes[1, 1].set_xlabel("time / T")
        axes[1, 1].set_ylabel("Velocity covariance diagonal / trace")
        finish(fig, axes, "Turnover-time blocks and resolved isotropy", time_axes=False)

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        axes[0, 0].plot(t, sample("cs2_mean") / p["c0"] ** 2)
        axes[0, 0].set_ylabel("Volume mean sound speed squared / c0^2")
        for ax, key, label in [
            (axes[0, 1], "density_percentiles", "Density percentiles"),
            (axes[1, 0], "sgs_mach_percentiles", "SGS Mach percentiles"),
        ]:
            values = sample(key)
            for j in range(values.shape[1]):
                levels = data.get("percentile_levels", [1, 10, 50, 90, 99])
                name = (
                    f"{levels[j]}%" if len(levels) == values.shape[1] else f"index {j}"
                )
                ax.plot(t, values[:, j], label=name)
            ax.set_ylabel(label)
        K0 = samples[0]["K"]
        reconstructed = (
            K0
            + history("sgs-production-integral")
            - history("sgs-production-integral")[0]
        )
        reconstructed -= (
            history("sgs-dissipation-integral") - history("sgs-dissipation-integral")[0]
        )
        axes[1, 1].plot(t, sample("K"), label="native snapshot K")
        axes[1, 1].plot(ht, reconstructed, "--", label="initial K + applied P - D")
        axes[1, 1].set_ylabel("Mean K; periodic box volume = 1")
        finish(
            fig,
            axes,
            "Thermal drift, intermittent tails, and independent SGS inventory",
        )

        final = data.get("final", {})
        if "density" in final and "sgs_mach" in final:
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            density = np.asarray(final["density"])
            axes[0].hist(np.log(density / rho0), bins=64, density=True, histtype="step")
            axes[0].set_xlabel("log(rho / rho0)")
            axes[0].set_ylabel("Volume probability density")
            mach = np.sort(final["sgs_mach"])
            survival = (len(mach) - np.arange(len(mach))) / len(mach)
            axes[1].step(mach, survival, where="pre")
            axes[1].set_xlabel("SGS Mach threshold")
            axes[1].set_ylabel("Volume fraction at or above threshold")
            axes[1].set_yscale("log")
            axes[1].axvline(1.0, color=".5", linestyle=":")
            finish(
                fig,
                axes,
                f"Final snapshot distributions at t/T={t[-1]:g}",
                time_axes=False,
            )
    result["plot_pages"] = pages
    (output / "statistics.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def self_check():
    t = np.array([0.0, 0.1, 0.7, 1.0])
    check = statistics(t, t, 0.0, 1.0)
    assert abs(check["mean"] - 0.5) < 1e-14
    assert abs(check["temporal_std"] - np.sqrt(1 / 12)) < 1e-14
    assert abs(check["slope"] - 1.0) < 1e-14
    check = statistics(t, np.full(4, 3.0), 0.2, 0.9)
    assert abs(check["mean"] - 3.0) < 1e-14 and check["temporal_std"] < 1e-14
    assert abs(counter_rate(t, 2 * t, 0.2, 0.9)["mean"] - 2.0) < 1e-14
    try:
        statistics(t, t, 0.0, 2.0)
    except ValueError:
        pass
    else:
        raise AssertionError("incomplete window accepted")
    print("paper statistics self-check passed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_check()
    if args.output:
        write_report(args.output)
    elif not args.self_test:
        parser.error("provide output directory or --self-test")
