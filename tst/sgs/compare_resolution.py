#!/usr/bin/env python3
"""Develop a 2D turbulent snapshot, then compare coarse SGS-on/off to finer hydro.

The high-resolution run is a numerical reference, not a proven DNS. Coarse data
are aligned finite-volume block averages. Removed kinetic variance initializes K;
SGS-off retains that energy as heat. Both lanes conserve the same initial E+K.
"""

import argparse
import html
import json
from pathlib import Path
import re
import time
import subprocess
import tempfile
import textwrap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from filter_reference import error_metrics
from laptop_2d import grid, spectrum
import verify
from verify_forcing import ForcingVerification

plt.switch_backend("Agg")
GAMMA = 5 / 3


def average(data, n):
    size = data.shape[0]
    if data.shape[1] != size or size % n or n > size:
        raise ValueError("require aligned square grids and an integer refinement ratio")
    r = size // n
    return data.reshape(n, r, n, r, *data.shape[2:]).mean(axis=(1, 3))


def restrict(data, n):
    """Conserved averages plus stable density-weighted unresolved kinetic variance."""
    coarse = average(data, n)
    velocity = data[..., 1:4] / data[..., :1]
    favre = coarse[..., 1:4] / coarse[..., :1]
    r = data.shape[0] // n
    expanded = favre.repeat(r, axis=0).repeat(r, axis=1)
    variance = average(
        0.5 * data[..., 0] * np.sum((velocity - expanded) ** 2, axis=-1), n
    )
    return coarse, variance


def physical(data, K):
    rho = data[..., 0]
    velocity = data[..., 1:4] / rho[..., None]
    thermal = data[..., 4] - 0.5 * rho * np.sum(velocity**2, axis=-1)
    dx = 1 / rho.shape[0]
    omega = (
        (np.roll(velocity[..., 1], -1, axis=1) - np.roll(velocity[..., 1], 1, axis=1))
        - (np.roll(velocity[..., 0], -1, axis=0) - np.roll(velocity[..., 0], 1, axis=0))
    ) / (2 * dx)
    return dict(
        rho=rho,
        velocity=velocity,
        speed=np.linalg.norm(velocity, axis=-1),
        pressure=(GAMMA - 1) * thermal,
        K=K,
        omega=omega,
        thermal=thermal,
        resolved_ke=0.5 * rho * np.sum(velocity**2, axis=-1),
    )


def reference_fields(data, n):
    coarse, K = restrict(data, n)
    # E is the resolved reservoir; represented E+K remains the conservative average.
    coarse[..., 4] -= K
    return physical(coarse, K)


def write_state(path, data, K, sgs):
    packed = np.concatenate((data.copy(), np.zeros(data.shape[:2] + (1,))), axis=-1)
    if sgs:
        packed[..., 4] -= K
        packed[..., -1] = K
    with path.open("w") as f:
        f.write(f"{data.shape[1]} {data.shape[0]} {packed.shape[-1]}\n")
        np.savetxt(f, packed.reshape(-1, packed.shape[-1]), fmt="%.17e")
    return packed


def self_test():
    data = np.zeros((2, 2, 7))
    data[..., 0] = [[1, 3], [1, 3]]
    data[..., 1] = [[1, 9], [1, 9]]
    data[..., 4] = 2 + data[..., 1] ** 2 / (2 * data[..., 0])
    coarse, K = restrict(data, 1)
    assert abs(coarse[0, 0, 0] - 2) < 1e-14
    assert abs(coarse[0, 0, 1] / 2 - 2.5) < 1e-14
    assert abs(K[0, 0] - 0.75) < 1e-14
    ref = reference_fields(data, 1)
    assert abs(ref["thermal"][0, 0] - 2) < 1e-14
    boosted = data.copy()
    boosted[..., 1] += 5 * data[..., 0]
    boosted[..., 4] += 5 * data[..., 1] + 0.5 * 25 * data[..., 0]
    assert np.max(abs(restrict(boosted, 1)[1] - K)) < 1e-14
    assert np.array_equal(average(data, 2), data)
    print("Favre hand case, energy partition, alignment and boost checks pass")


def initial_file_checks(exe):
    config = verify.parameters(
        (verify.SOURCE / "inputs/hydro/sgs_verify.athinput").read_text()
    )
    config.update(
        {
            "mesh/nx1": "4",
            "mesh/nx2": "4",
            "mesh/nx3": "1",
            "meshblock/nx1": "4",
            "meshblock/nx2": "4",
            "meshblock/nx3": "1",
            "time/tlim": "0",
            "time/nlim": "0",
        }
    )
    good = np.zeros((4, 4, 8))
    good[..., 0] = 1
    good[..., 4] = 1.5
    good[..., 5:7] = [0.2, 0.25]
    good[..., -1] = 0.01
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for name in (
            "valid",
            "header",
            "truncated",
            "rho",
            "thermal",
            "K",
            "nonfinite",
            "extra",
            "off_K",
        ):
            run = root / name
            run.mkdir()
            data = good.copy()
            if name in ("rho", "thermal", "K", "nonfinite"):
                index = {"rho": 0, "thermal": 4, "K": 7, "nonfinite": 1}[name]
                data[0, 0, index] = {
                    "rho": 0,
                    "thermal": -1,
                    "K": -0.1,
                    "nonfinite": float("nan"),
                }[name]
            path = run / "state.dat"
            with path.open("w") as f:
                f.write("3 4 8\n" if name == "header" else "4 4 8\n")
                np.savetxt(
                    f,
                    data.reshape(-1, 8)[:-1]
                    if name == "truncated"
                    else data.reshape(-1, 8),
                    fmt="%.17e",
                )
                if name == "extra":
                    f.write("unexpected\n")
            case = dict(
                config,
                **{
                    "job/basename": name,
                    "problem/initial_state": str(path),
                    "sgs/enabled": str(name != "off_K").lower(),
                },
            )
            inp = run / "input.athinput"
            inp.write_text(verify.input_text(case))
            process = subprocess.run(
                [str(exe.resolve()), "-i", str(inp)],
                cwd=run,
                capture_output=True,
                text=True,
                timeout=10,
            )
            assert (process.returncode == 0) == (name == "valid"), (
                process.stdout + process.stderr
            )
            if name == "valid":
                loaded = grid(verify.restart(sorted((run / "rst").glob("*.rst"))[0]))
                assert np.array_equal(loaded, good)
            else:
                assert "invalid 2D initial_state" in process.stdout + process.stderr
    print("Native initial-state round-trip and eight malformed-input checks pass")


def run_cases(args):
    suite = ForcingVerification(args.exe, args.output)
    manifest_path = suite.output / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["resolution_comparison"] = dict(
        high=args.high,
        low=args.low,
        spinup=args.spinup,
        decay=args.decay,
        grid_cap=args.high,
        protocol="Developed flow; aligned conservative/Favre restriction; unforced decay",
        numerical_reference_is_validated_DNS=False,
    )
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    # All forcing acts only during preparation. The comparison itself is unforced.
    suite.run_forcing(
        "spinup",
        False,
        0.5,
        ndim=2,
        tlim=args.spinup,
        overrides={
            "mesh/nx1": str(args.high),
            "mesh/nx2": str(args.high),
            "meshblock/nx1": "16",
            "meshblock/nx2": "16",
            "turb_driving/min_kz": "0",
            "turb_driving/max_kz": "0",
            "turb_driving/dedt": "1.0",
            "time/nlim": "20000",
            "output1/dt": str(args.spinup / 20),
            "output4/dt": str(args.spinup),
        },
    )
    if "spinup" not in suite.states or suite.results[-1]["status"] != "pass":
        return suite.results
    fine = grid(suite.states["spinup"][1])
    initial = {}
    initial_errors = {}
    for name, n, enabled in (
        ("reference", args.high, False),
        ("reference_check", args.high // 2, False),
        ("low_off", args.low, False),
        ("low_on", args.low, True),
    ):
        data, K = restrict(fine, n)
        inp = suite.output / f"{name}-initial.dat"
        expected = write_state(inp, data, K, enabled)
        initial[name] = dict(n=n, enabled=enabled, expected=expected)
        suite.run(
            name,
            mode="box",
            ndim=2,
            block=16,
            grid_cap=args.high,
            overrides={
                "mesh/nx1": str(n),
                "mesh/nx2": str(n),
                "time/tlim": str(args.decay),
                "time/nlim": "20000",
                "time/ndiag": "1000",
                "sgs/enabled": str(enabled).lower(),
                "problem/initial_state": str(inp),
                "output1/dt": str(args.decay / 20),
                "output4/dt": str(args.decay / 10),
            },
        )
        suite.results[-1]["initial_state_sha256"] = verify.sha(inp)
        if name in suite.states:
            loaded = grid(suite.states[name][0])
            target = expected if enabled else expected[..., :-1]
            initial_errors[name] = float(np.max(abs(loaded - target)))
    row = dict(name="matched_initial_state", status="blocked")
    suite.results.append(row)
    if len(initial_errors) == 4:
        energies = {
            name: float(np.mean(case["expected"][..., 4] + case["expected"][..., -1]))
            for name, case in initial.items()
        }
        row.update(loaded_linf=initial_errors, represented_energy=energies)
        suite.check(
            row,
            "initial.file_readback",
            max(initial_errors.values()),
            0,
            "ASCII conserved data must round-trip exactly into native checkpoints",
        )
        suite.check(
            row,
            "initial.energy_match",
            max(energies.values()) - min(energies.values()),
            1e-12,
            "all lanes retain the same global represented energy",
        )
        row["status"] = "pass" if all(c["passed"] for c in row["checks"]) else "fail"
    suite.save()
    return suite.results


def report(output):
    rows = json.loads((output / "results.json").read_text())
    states = {
        name: [verify.restart(p) for p in sorted((output / name / "rst").glob("*.rst"))]
        for name in ("reference", "reference_check", "low_off", "low_on")
    }
    n = int(states["low_on"][0]["config"]["mesh/nx1"])
    high_n = int(states["reference"][0]["config"]["mesh/nx1"])
    times = np.linspace(0, states["reference"][-1]["time"], 11)
    spans = {}

    def at(name, t):
        snaps = states[name]
        ts = np.array([s["time"] for s in snaps])
        index = int(np.searchsorted(ts, t))
        if index == 0 or index == len(ts):
            return grid(snaps[0 if index == 0 else -1])
        if abs(ts[index] - t) < 1e-14:
            return grid(snaps[index])
        a, b = snaps[index - 1], snaps[index]
        weight = (t - a["time"]) / (b["time"] - a["time"])
        spans[name] = max(spans.get(name, 0.0), b["time"] - a["time"])
        return (1 - weight) * grid(a) + weight * grid(b)

    def comparison(t):
        reference = reference_fields(at("reference", t), n)
        check = reference_fields(at("reference_check", t), n)
        lanes = {}
        for name in ("low_off", "low_on"):
            data = at(name, t)
            K = data[..., -1] if name == "low_on" else np.zeros((n, n))
            lanes[name] = physical(data[..., :7], K)
        return reference, check, lanes

    comparisons = [comparison(float(t)) for t in times]
    metrics = []
    for t, (ref, check, lanes) in zip(times, comparisons):
        metrics.append(
            dict(
                time=float(t),
                reference_sensitivity={
                    key: error_metrics(ref[key], check[key])
                    for key in ("rho", "velocity", "pressure", "K")
                },
                errors={
                    name: {
                        key: error_metrics(ref[key], lane[key])
                        for key in ("rho", "velocity", "pressure", "K", "omega")
                    }
                    for name, lane in lanes.items()
                },
            )
        )
    plots = output / "plots"
    plots.mkdir(exist_ok=True)
    gallery = []
    pdf = PdfPages(output / "report.pdf")
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.2})

    def save(name, title, caption, shape, draw):
        fig, axes = plt.subplots(
            *shape, figsize=(5 * shape[1], 3.6 * shape[0]), layout="constrained"
        )
        fig.suptitle(title, fontsize=15)
        fig.supxlabel(textwrap.fill(caption, width=125), fontsize=9)
        flat = np.asarray(axes).ravel()
        draw(flat)
        for ax in flat:
            if ax.lines:
                ax.legend(fontsize=8)
        fig.savefig(plots / f"{name}.png", dpi=140)
        pdf.savefig(fig)
        plt.close(fig)
        gallery.append((name, title, caption, len(flat)))

    def image(ax, data, title, bounds=None, difference=False):
        if difference:
            maximum = max(float(np.max(abs(data))), 1e-16)
            bounds = (-maximum, maximum)
        im = ax.imshow(
            data,
            origin="lower",
            extent=(-0.5, 0.5, -0.5, 0.5),
            cmap="RdBu_r" if difference else "viridis",
            vmin=bounds[0] if bounds else None,
            vmax=bounds[1] if bounds else None,
        )
        ax.set(title=title, xlabel="x", ylabel="y")
        ax.grid(False)
        ax.figure.colorbar(im, ax=ax, shrink=0.85)

    ref, check, lanes = comparisons[-1]
    high = physical(
        grid(states["reference"][-1]), np.zeros(grid(states["reference"][-1]).shape[:2])
    )

    def field_maps(axes):
        all_fields = [high, ref, lanes["low_off"], lanes["low_on"]]
        for j, key in enumerate(("rho", "speed", "omega")):
            bounds = (
                min(f[key].min() for f in all_fields),
                max(f[key].max() for f in all_fields),
            )
            for i, (label, f) in enumerate(
                zip(
                    ("High native", "High filtered", "Low SGS off", "Low SGS on"),
                    all_fields,
                )
            ):
                image(axes[3 * i + j], f[key], label + ": " + key, bounds)

    save(
        "01_fields",
        "High-resolution turbulence versus coarse SGS-on/off",
        "Final matched time. The high reference is block-filtered onto "
        "the low grid; column color scales are shared.",
        (4, 3),
        field_maps,
    )

    def errors_map(axes):
        for i, name in enumerate(("low_off", "low_on")):
            for j, key in enumerate(("rho", "speed", "pressure")):
                image(
                    axes[3 * i + j],
                    lanes[name][key] - ref[key],
                    name + " − filtered high: " + key,
                    difference=True,
                )

    save(
        "02_error_maps",
        "Coarse-field errors against the filtered high-resolution run",
        "Signed final errors. Each panel uses its own symmetric color "
        "range; use the quantitative error plots to compare magnitudes.",
        (2, 3),
        errors_map,
    )

    def errors_time(axes):
        for ax, key in zip(axes, ("rho", "velocity", "pressure", "omega")):
            for name in ("low_off", "low_on"):
                ax.plot(
                    times,
                    [m["errors"][name][key]["rms_error"] for m in metrics],
                    label=name,
                )
            if key != "omega":
                ax.plot(
                    times,
                    [m["reference_sensitivity"][key]["rms_error"] for m in metrics],
                    "k--",
                    label=f"{high_n // 2}² vs {high_n}² reference sensitivity",
                )
            ax.set(title=key + " absolute RMS error", xlabel="Decay time")

    save(
        "03_errors",
        "Does SGS improve the coarse resolved fields?",
        "Intermediate diagnostics interpolate conserved states to common "
        "times. The final comparison time is exact in every lane.",
        (2, 2),
        errors_time,
    )

    def energy(axes):
        for name, label in [
            ("filtered_high", "Filtered high"),
            ("low_off", "Low off"),
            ("low_on", "Low on"),
        ]:
            fs = [c[0] if name == "filtered_high" else c[2][name] for c in comparisons]
            for ax, key in zip(axes, ("resolved_ke", "K", "thermal")):
                ax.plot(times, [float(f[key].mean()) for f in fs], label=label)
                ax.set(title=key + " volume mean", xlabel="Decay time")

    save(
        "04_energy",
        "Resolved kinetic, SGS and thermal energy relative to the reference",
        "Filtered-reference K is kinetic variance below the coarse cell "
        "scale. SGS-off retains that initial variance as heat.",
        (1, 3),
        energy,
    )

    def sgs_maps(axes):
        bounds = (0, max(ref["K"].max(), lanes["low_on"]["K"].max()))
        image(axes[0], ref["K"], "Filtered high: unresolved K", bounds)
        image(axes[1], lanes["low_on"]["K"], "Low SGS on: modeled K", bounds)
        image(
            axes[2],
            lanes["low_on"]["K"] - ref["K"],
            "Modeled − reference K",
            difference=True,
        )
        axes[3].scatter(ref["K"].ravel(), lanes["low_on"]["K"].ravel(), s=4, alpha=0.5)
        axes[3].plot(bounds, bounds, "k--", label="Exact agreement")
        axes[3].set(
            xlabel="Reference K",
            ylabel="Modeled K",
            title="Cellwise SGS energy comparison",
        )

    save(
        "05_sgs",
        "Modeled SGS reservoir versus filtered kinetic variance",
        "SGS-on starts from the exact filtered variance. Agreement at the "
        "initial time is initialization, not a prediction.",
        (2, 2),
        sgs_maps,
    )

    def spectra(axes):
        for ax, initial in zip(axes, (True, False)):
            r, _, ls = comparisons[0 if initial else -1]
            for label, f in [
                ("Filtered high", r),
                ("Low off", ls["low_off"]),
                ("Low on", ls["low_on"]),
            ]:
                power = spectrum(f["rho"], f["velocity"])
                ax.semilogy(
                    np.arange(1, len(power)),
                    np.maximum(power[1:], 1e-30),
                    "o-",
                    ms=3,
                    label=label,
                )
            ax.set(
                title="Initial" if initial else "Final",
                xlabel="Radial mode |k|",
                ylabel="Shell energy",
            )

    save(
        "06_spectra",
        "Spectra compared at one common coarse scale",
        "Shell energy of sqrt(rho) velocity fluctuations about the Favre "
        "mean. All plotted fields use the same low grid.",
        (1, 2),
        spectra,
    )

    def reference_check(axes):
        for ax, key in zip(axes, ("rho", "velocity", "pressure", "K")):
            ax.plot(
                times,
                [m["reference_sensitivity"][key]["rms_error"] for m in metrics],
                label=f"{high_n // 2}² vs {high_n}², both filtered to low grid",
            )
            ax.set(title=key + " reference sensitivity", xlabel="Decay time")

    save(
        "07_reference_sensitivity",
        "How much does the numerical reference change with resolution?",
        "A nonzero reference sensitivity limits claims about SGS "
        "accuracy. Neither resolution is assumed to be a validated DNS.",
        (2, 2),
        reference_check,
    )

    def budgets(axes):
        for name in states:
            snaps = states[name]
            ts = [s["time"] for s in snaps]
            e = [
                s["totals"][4] + (s["totals"][-1] if name == "low_on" else 0)
                for s in snaps
            ]
            axes[0].plot(ts, np.array(e) - e[0], label=name)
            axes[1].plot(ts, [s["min_thermal"] for s in snaps], label=name)
            if name == "low_on":
                h = verify.history(output / name / f"{name}.hydro.hst")
                axes[2].plot(
                    [q["time"] for q in h],
                    [q["sgs-production-integral"] for q in h],
                    label="Applied production",
                )
                axes[2].plot(
                    [q["time"] for q in h],
                    [q["sgs-dissipation-integral"] for q in h],
                    label="Dissipation to heat",
                )
                axes[3].plot(
                    [q["time"] for q in h],
                    [q["sgs-fallback-stages"] for q in h],
                    label="Fallback stages",
                )
        for ax, title in zip(
            axes,
            (
                "Represented-energy drift",
                "Minimum thermal energy",
                "Applied SGS transfers",
                "Fallback use",
            ),
        ):
            ax.set(title=title, xlabel="Decay time")

    save(
        "08_budgets",
        "Conservation, energy transfer and admissibility",
        "Unforced decay keeps represented E+K constant. Applied SGS "
        "transfers are separate from implicit numerical heating.",
        (2, 2),
        budgets,
    )

    costs = {}
    for name in ("spinup", *states):
        log = (output / name / "run.log").read_text()
        costs[name] = float(re.search(r"cpu time used\s*=\s*([\deE+.-]+)", log)[1])

    def cost(axes):
        axes[0].bar(costs.keys(), costs.values())
        axes[0].set(
            ylabel="Elapsed solver seconds", title="Local execution cost", xlabel="Run"
        )
        keys = ("rho", "velocity", "pressure", "K")
        for j, name in enumerate(("low_off", "low_on")):
            axes[1].bar(
                np.arange(4) + j * 0.35,
                [metrics[-1]["errors"][name][k]["relative_rms_error"] for k in keys],
                width=0.35,
                label=name,
            )
        axes[1].set(
            xticks=np.arange(4) + 0.175,
            xticklabels=keys,
            ylabel="Relative RMS error",
            title="Final errors",
        )

    save(
        "09_cost",
        "Numerical accuracy and measured local cost",
        "Elapsed timer, single runs with output. An SGS improvement is a "
        "measured outcome, never a condition for passing this test.",
        (1, 2),
        cost,
    )
    pdf.close()
    summary = dict(
        final_time=float(times[-1]),
        high_grid=high["rho"].shape[0],
        low_grid=n,
        figure_pages=len(gallery),
        panels=sum(p[3] for p in gallery),
        final_errors=metrics[-1],
        solver_elapsed_seconds=costs,
        interpolation_max_bracket=spans,
        passing_records=sum(r["status"] == "pass" for r in rows),
        result_records=len(rows),
        qualification="2D decaying numerical comparison from one developed seed; "
        "not validated DNS or 3D turbulence",
    )
    (output / "comparison.json").write_text(
        json.dumps(dict(summary=summary, time_series=metrics), indent=2) + "\n"
    )
    content = "".join(
        f"<h2>{html.escape(title)}</h2><p>{html.escape(caption)}</p>"
        f'<a href="plots/{name}.png"><img src="plots/{name}.png" '
        f'alt="{html.escape(title)}"></a>'
        for name, title, caption, _ in gallery
    )
    (output / "index.html").write_text(
        '<!doctype html><meta charset="utf-8">'
        "<title>High-resolution versus low-resolution SGS turbulence</title>"
        "<style>body{max-width:1200px;margin:40px auto;font:17px/1.5 "
        "system-ui;padding:0 20px}"
        "img{width:100%}pre{white-space:pre-wrap}a{color:#1265a8}</style>"
        "<h1>High-resolution turbulence versus coarse SGS-on/off</h1>"
        "<p>Developed 2D flow, conservative matched initialization, then unforced decay. "
        "The finer hydro run is a numerical reference, not a proven DNS.</p>"
        "<p>SGS-on stores removed kinetic variance as K; SGS-off stores it as heat. "
        "Both start with the same filtered mass, momentum, scalar masses "
        "and represented energy.</p>"
        '<p><a href="report.pdf">Full PDF</a> · <a '
        'href="comparison.json">Measured errors</a> ·'
        '<a href="results.json">Acceptance checks</a> · <a '
        'href="manifest.json">Provenance</a></p>'
        + content
        + "<h2>Measured summary</h2><pre>"
        + html.escape(json.dumps(summary, indent=2))
        + "</pre>"
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--high", type=int, default=128, choices=(64, 128))
    parser.add_argument("--low", type=int, default=32, choices=(16, 32))
    parser.add_argument("--spinup", type=float, default=2.0)
    parser.add_argument("--decay", type=float, default=0.5)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        if args.exe is not None:
            initial_file_checks(args.exe)
        return 0
    if args.output is None or (args.exe is None and not args.plot_only):
        parser.error("--output and (unless --plot-only) --exe are required")
    if args.high < 4 * args.low or not 0 < args.spinup <= 3 or not 0 < args.decay <= 1:
        parser.error("require high >= 4*low, 0 < spinup <= 3 and 0 < decay <= 1")
    args.output = args.output.resolve()
    started = time.monotonic()
    if not args.plot_only:
        rows = run_cases(args)
        if any(r["status"] != "pass" for r in rows) or len(rows) != 6:
            print("Comparison incomplete; inspect retained results.json and logs")
            return 1
    summary = report(args.output)
    summary["workflow_wall_seconds"] = time.monotonic() - started
    (
        args.output / ("plot-execution.json" if args.plot_only else "execution.json")
    ).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    print("Open", args.output / "index.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
