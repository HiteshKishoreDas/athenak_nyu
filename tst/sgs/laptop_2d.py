#!/usr/bin/env python3
"""Run small 2D SGS checks and create a PNG/PDF/HTML visual report.

Uses the existing deterministic problem generator and verification checks. No
solver changes. NumPy and Matplotlib are the only plotting dependencies.
"""

import argparse
import html
import json
from pathlib import Path
import re
import time
import textwrap

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import verify
from verify_forcing import ForcingVerification

plt.switch_backend("Agg")


def grid(state):
    """Assemble active cells in image order (y, x, variable), without ghost cells."""
    nx = int(state["config"]["mesh/nx1"])
    ny = int(state["config"]["mesh/nx2"])
    data = np.empty((ny, nx, len(state["totals"])))
    assert len(state["field"]) == nx * ny
    for (x, y, z), values in state["field"].items():
        assert z == 0 and 0 <= x < nx and 0 <= y < ny
        data[y, x] = values
    return data


def fields(state):
    data = grid(state)
    rho = data[..., 0]
    velocity = data[..., 1:4] / rho[..., None]
    kinetic = 0.5 * rho * np.sum(velocity**2, axis=-1)
    gamma = float(state["config"]["hydro/gamma"])
    K = data[..., -1] if state["min_sgs"] is not None else np.zeros_like(rho)
    dx = 1 / rho.shape[1]
    dy = 1 / rho.shape[0]
    omega = (
        np.roll(velocity[..., 1], -1, axis=1) - np.roll(velocity[..., 1], 1, axis=1)
    ) / (2 * dx) - (
        np.roll(velocity[..., 0], -1, axis=0) - np.roll(velocity[..., 0], 1, axis=0)
    ) / (2 * dy)
    return dict(
        rho=rho,
        pressure=(gamma - 1) * (data[..., 4] - kinetic),
        K=K,
        k=K / rho,
        u=velocity[..., 0],
        velocity=velocity,
        speed=np.linalg.norm(velocity, axis=-1),
        omega=omega,
        scalar=data[..., 5] / rho,
    )


def spectrum(rho, velocity):
    """Shell-integrated spectrum of sqrt(rho)*(u-Favre mean); sum = fluctuation KE."""
    mean = np.sum(rho[..., None] * velocity, axis=(0, 1)) / np.sum(rho)
    modes = np.fft.fft2(
        np.sqrt(rho[..., None]) * (velocity - mean), axes=(0, 1), norm="forward"
    )
    power = 0.5 * np.sum(np.abs(modes) ** 2, axis=-1)
    ny, nx = rho.shape
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    shells = np.rint(np.hypot(kx[None, :], ky[:, None])).astype(int)
    return np.bincount(shells.ravel(), weights=power.ravel())


def self_test():
    n = 16
    x = (np.arange(n) + 0.5) / n
    xx, yy = np.meshgrid(x, x)
    velocity = np.zeros((n, n, 3))
    velocity[..., 0] = np.sin(2 * np.pi * yy)
    spec = spectrum(np.ones((n, n)), velocity)
    assert abs(spec.sum() - 0.25) < 1e-14
    assert abs(spec[1] - 0.25) < 1e-14
    state = dict(
        config={
            "mesh/nx1": str(n),
            "mesh/nx2": str(n),
            "hydro/gamma": "1.6666666666666667",
        },
        totals=[0] * 8,
        min_sgs=0,
        field={
            (i, j, 0): (
                1,
                velocity[j, i, 0],
                0,
                0,
                1.5 + 0.5 * velocity[j, i, 0] ** 2,
                xx[j, i],
                0.25,
                0.01,
            )
            for j in range(n)
            for i in range(n)
        },
    )
    f = fields(state)
    assert np.array_equal(f["scalar"], xx)
    assert np.max(abs(f["pressure"] - 1)) < 1e-14
    expected = -n * np.sin(2 * np.pi / n) * np.cos(2 * np.pi * yy)
    assert np.max(abs(f["omega"] - expected)) < 1e-13
    print("Grid orientation, pressure, periodic vorticity and FFT Parseval checks pass")


def run_cases(exe, output):
    suite = ForcingVerification(exe, output)

    def run(name, mode, n=32, enabled=True, block=16, rst=None, **changes):
        overrides = {
            "mesh/nx1": str(n),
            "mesh/nx2": str(n),
            "sgs/enabled": str(enabled).lower(),
            "time/tlim": "0.2",
            "time/ndiag": "100",
            "output1/dt": "0.01",
            "output4/dt": "0.01",
        }
        overrides.update(changes)
        return suite.run(
            name, mode=mode, ndim=2, block=block, overrides=overrides, rst=rst
        )

    decay = []
    for cfl in ("0.2", "0.1", "0.05"):
        decay.append(
            run(
                "decay_cfl" + cfl,
                "uniform",
                **{"time/cfl_number": cfl, "time/tlim": "0.1", "problem/k0": "0.04"},
            )
        )
    row = dict(name="decay_2d_temporal_order", status="blocked")
    suite.results.append(row)
    if all(r["status"] == "pass" for r in decay):
        errors = [r["decay_error"] for r in decay]
        orders = [np.log2(errors[i] / errors[i + 1]) for i in (0, 1)]
        row.update(
            errors=errors,
            orders=[float(p) for p in orders],
            status="pass"
            if all(
                verify.LIMITS["decay_order_min"] <= p <= verify.LIMITS["decay_order_max"]
                for p in orders
            )
            else "fail",
        )
    suite.save()
    for n in (16, 32):
        for enabled in (False, True):
            run(
                f"shear_n{n}_sgs{int(enabled)}",
                "shear",
                n,
                enabled,
                **{"problem/velocity": "0.6", "problem/k0": "0.04"},
            )
            run(
                f"box_n{n}_sgs{int(enabled)}",
                "box",
                n,
                enabled,
                **{"time/tlim": "0.3", "problem/velocity": "0.6", "problem/k0": "0.04"},
            )
        run(
            f"diffusion_n{n}",
            "diffusion",
            n,
            **{"time/tlim": "0.12", "problem/k0": "0.05", "problem/amplitude": "0.5"},
        )
        run(f"contact_n{n}", "contact", n, **{"sgs/ceps": "0", "sgs/ckappa": "0"})
    run("zero_n32", "uniform", **{"problem/k0": "0", "problem/velocity": "0"})
    run(
        "box_n32_single",
        "box",
        block=32,
        **{"time/tlim": "0.3", "problem/velocity": "0.6", "problem/k0": "0.04"},
    )
    suite.compare(
        "box_n32_sgs1",
        "box_n32_single",
        "I04",
        verify.LIMITS["decomposition_scaled_linf"],
    )
    if "box_n32_sgs1" in suite.states:
        checkpoints = sorted((output / "box_n32_sgs1/rst").glob("*.rst"))
        chosen = next(p for p in checkpoints if 0.1 <= verify.restart(p)["time"] < 0.2)
        run(
            "box_n32_restart",
            "box",
            rst=chosen,
            **{"time/tlim": "0.3", "problem/velocity": "0.6", "problem/k0": "0.04"},
        )
        suite.compare(
            "box_n32_sgs1",
            "box_n32_restart",
            "I13",
            verify.LIMITS["checkpoint_restart_scaled_linf"],
        )
    for label, fraction in (("sol", 1.0), ("comp", 0.0)):
        for enabled in (False, True):
            suite.run_forcing(
                f"forced_{label}_sgs{int(enabled)}",
                enabled,
                fraction,
                ndim=2,
                tlim=1,
                overrides={
                    "mesh/nx1": "32",
                    "mesh/nx2": "32",
                    "meshblock/nx1": "16",
                    "meshblock/nx2": "16",
                    "turb_driving/min_kz": "0",
                    "turb_driving/max_kz": "0",
                    "output1/dt": "0.05",
                    "output4/dt": "0.05",
                },
            )
    return suite.results


def make_plots(output):
    rows = json.loads((output / "results.json").read_text())
    runs = {r["name"]: r for r in rows if "command" in r and r["status"] == "pass"}
    snapshots = {
        name: [verify.restart(p) for p in sorted((output / name / "rst").glob("*.rst"))]
        for name in runs
    }
    histories = {
        name: verify.history(output / name / f"{name}.hydro.hst") for name in runs
    }
    figures = output / "plots"
    figures.mkdir(exist_ok=True)
    gallery = []
    panel_count = 0
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.2})
    pdf = PdfPages(output / "report.pdf")

    def page(name, title, caption, shape=(2, 2), size=(12, 8)):
        fig, axes = plt.subplots(*shape, figsize=size, layout="constrained")
        fig.suptitle(title, fontsize=15)
        fig.supxlabel(textwrap.fill(caption, width=130), fontsize=9)
        return fig, np.asarray(axes).ravel(), name, caption

    def save(item):
        nonlocal panel_count
        fig, axes, name, caption = item
        for ax in axes:
            if ax.lines:
                ax.legend(fontsize=8)
            if ax.has_data() and not ax.images:
                ax.set_xlabel(ax.get_xlabel() or "Time")
        fig.savefig(figures / f"{name}.png", dpi=145)
        pdf.savefig(fig)
        plt.close(fig)
        gallery.append((name, caption))
        panel_count += len(axes)

    def curve(ax, name, label, function):
        h = histories[name]
        ax.plot([q["time"] for q in h], [function(q) for q in h], label=label)

    def image(ax, values, title, limits=None, difference=False):
        if difference:
            bound = max(float(np.max(abs(values))), 1e-16)
            limits = (-bound, bound)
        im = ax.imshow(
            values,
            origin="lower",
            extent=(-0.5, 0.5, -0.5, 0.5),
            cmap="RdBu_r" if difference else "viridis",
            vmin=limits[0] if limits else None,
            vmax=limits[1] if limits else None,
        )
        ax.set(title=title, xlabel="x", ylabel="y")
        ax.grid(False)
        ax.figure.colorbar(im, ax=ax, shrink=0.85)

    def kinetic(q):
        return sum(q[f"{v}-KE"] for v in (1, 2, 3))

    def represented(q):
        return q["tot-E"] + q.get("sgs-E", 0)

    def thermal(q):
        return q["tot-E"] - kinetic(q)

    item = page(
        "01_decay",
        "2D uniform SGS decay: exact solution and RK2 error",
        "Uniform K decays into heat. Dashed curves are analytic, not a fitted reference.",
    )
    ax = item[1]
    errors = []
    for cfl in ("0.2", "0.1", "0.05"):
        name = "decay_cfl" + cfl
        states = snapshots[name]
        times = np.array([s["time"] for s in states])
        exact = 0.04 / (1 + 1.5 * 32 * np.sqrt(0.04) * times / 2) ** 2
        observed = np.array([s["totals"][-1] for s in states])
        ax[0].plot(times, observed, ".", label="CFL " + cfl)
        ax[1].plot(
            times,
            [s["thermal"] - states[0]["thermal"] for s in states],
            ".",
            label="CFL " + cfl,
        )
        ax[2].semilogy(times[1:], abs(observed - exact)[1:], label="CFL " + cfl)
        errors.append(abs(observed[-1] - exact[-1]))
    ax[0].plot(times, exact, "k--", label="Exact K(t)")
    ax[1].plot(times, 0.04 - exact, "k--", label="Exact heat gain")
    ax[3].loglog([0.2, 0.1, 0.05], errors, "o-", label="Final absolute error")
    ax[3].loglog(
        [0.2, 0.05], [errors[0], errors[0] / 16], "k--", label="Second-order guide"
    )
    for a, title in zip(
        ax,
        (
            "SGS energy density",
            "Thermal energy gain",
            "Absolute K error",
            "Time-step refinement",
        ),
    ):
        a.set_title(title)
    ax[3].set_xlabel("CFL (fixed 32² spatial grid)")
    save(item)

    item = page(
        "02_shear",
        "2D shear: SGS-on/off and resolution sensitivity",
        "LLF also damps velocity. Grid refinement changes the SGS filter width; "
        "this comparison is not convergence to one fixed SGS PDE.",
        (2, 3),
        (14, 8),
    )
    ax = item[1]
    for n in (16, 32):
        y = -0.5 + (np.arange(n) + 0.5) / n
        for enabled in (0, 1):
            name = f"shear_n{n}_sgs{enabled}"
            f = fields(snapshots[name][-1])
            label = f"{n}² SGS {'on' if enabled else 'off'}"
            ax[0].plot(y, f["u"].mean(axis=1), label=label)
            curve(ax[3], name, label, kinetic)
            curve(ax[4], name, label, lambda q: thermal(q) - thermal(histories[name][0]))
            if enabled:
                ax[1].plot(y, f["K"].mean(axis=1), label=label)
                ax[2].plot(y, f["pressure"].mean(axis=1), label=label)
                curve(
                    ax[5],
                    name,
                    f"{n}² production",
                    lambda q: q["sgs-production-integral"],
                )
                curve(
                    ax[5],
                    name,
                    f"{n}² dissipation",
                    lambda q: q["sgs-dissipation-integral"],
                )
    ax[0].plot(y, 0.6 * np.sin(2 * np.pi * (y + 0.5)), "k--", label="Initial velocity")
    for a, title in zip(
        ax,
        (
            "Final uₓ(y)",
            "Final K(y)",
            "Final thermal pressure",
            "Resolved kinetic energy",
            "Thermal gain",
            "Applied SGS transfers",
        ),
    ):
        a.set_title(title)
    for a in ax[:3]:
        a.set_xlabel("y")
    save(item)

    item = page(
        "03_diffusion",
        "2D pressure-balanced SGS profile: transport and decay",
        "The full mixed model evolves here; "
        "this is not an exact isolated-diffusion oracle.",
        (2, 3),
        (14, 8),
    )
    ax = item[1]
    for n, offset in ((16, 0), (32, 3)):
        states = snapshots[f"diffusion_n{n}"]
        x = -0.5 + (np.arange(n) + 0.5) / n
        for s in (states[0], states[len(states) // 2], states[-1]):
            f = fields(s)
            label = f"t={s['time']:.3f}"
            ax[offset].plot(x, f["K"].mean(axis=0), label=label)
            ax[offset + 1].plot(x, f["pressure"].mean(axis=0), label=label)
            ax[offset + 2].plot(
                x, (f["pressure"] + 2 * f["K"] / 3).mean(axis=0), label=label
            )
        for a, title in zip(
            ax[offset : offset + 3], ("K", "Thermal pressure", "P + 2K/3")
        ):
            a.set(title=f"{n}² {title}", xlabel="x")
    save(item)

    item = page(
        "04_contact",
        "2D moving contact: exact periodic advection",
        "Dissipation and SGS diffusion are disabled in this isolated check; "
        "constant velocity and K leave an exact advected density profile.",
    )
    ax = item[1]
    contact_errors = []
    for n in (16, 32):
        state = snapshots[f"contact_n{n}"][-1]
        f = fields(state)
        x = -0.5 + (np.arange(n) + 0.5) / n
        exact = 1 + 0.2 * np.cos(2 * np.pi * (x + 0.5 - 0.3 * state["time"]))
        difference = f["rho"] - exact[None, :]
        contact_errors.append(float(np.mean(abs(difference))))
        ax[0].plot(x, f["rho"].mean(axis=0), "o-", ms=3, label=f"{n}² numerical")
        ax[1].plot(x, difference.mean(axis=0), label=f"{n}² numerical − exact")
    ax[0].plot(x, exact, "k--", label="Exact at t=0.2")
    image(ax[2], difference, "32² density error", difference=True)
    ax[3].loglog([16, 32], contact_errors, "o-", label="Mean |density error|")
    ax[3].set(xlabel="Cells per axis", title="Two-grid error decrease")
    ax[0].set(title="Final density", xlabel="x")
    ax[1].set(title="Density error profile", xlabel="x")
    save(item)

    item = page(
        "05_box_fields",
        "2D decaying box: initial and final SGS-on fields",
        "Top: initial; bottom: t=0.3. Matching color scales within each field.",
        (2, 4),
        (16, 8),
    )
    ax = item[1]
    initial = fields(snapshots["box_n32_sgs1"][0])
    final = fields(snapshots["box_n32_sgs1"][-1])
    for j, (key, label) in enumerate(
        (
            ("rho", "Density"),
            ("K", "SGS K"),
            ("speed", "Speed"),
            ("omega", "Vorticity ωz"),
        )
    ):
        bounds = (
            min(initial[key].min(), final[key].min()),
            max(initial[key].max(), final[key].max()),
        )
        image(ax[j], initial[key], "Initial " + label, bounds)
        image(ax[j + 4], final[key], "Final " + label, bounds)
    save(item)

    item = page(
        "06_box_on_off",
        "2D decaying box: matched SGS-on/off fields at t=0.3",
        "Rows: SGS off, SGS on, on minus off. "
        "Differences have separate symmetric scales.",
        (3, 3),
        (14, 12),
    )
    ax = item[1]
    off = fields(snapshots["box_n32_sgs0"][-1])
    for j, (key, label) in enumerate(
        (("rho", "Density"), ("speed", "Speed"), ("omega", "Vorticity ωz"))
    ):
        bounds = (
            min(off[key].min(), final[key].min()),
            max(off[key].max(), final[key].max()),
        )
        image(ax[j], off[key], "Off: " + label, bounds)
        image(ax[j + 3], final[key], "On: " + label, bounds)
        image(ax[j + 6], final[key] - off[key], "Difference: " + label, difference=True)
    save(item)

    item = page(
        "07_box_evolution",
        "2D box: evolution of density, SGS energy and vorticity",
        "Four saved times on a 32² grid. Each row uses one color scale.",
        (3, 4),
        (16, 11),
    )
    ax = item[1]
    states = snapshots["box_n32_sgs1"]
    choices = [states[i] for i in np.linspace(0, len(states) - 1, 4).astype(int)]
    maps = [fields(s) for s in choices]
    for row, key in enumerate(("rho", "K", "omega")):
        bounds = (min(f[key].min() for f in maps), max(f[key].max() for f in maps))
        for col, (s, f) in enumerate(zip(choices, maps)):
            image(ax[4 * row + col], f[key], f"{key}, t={s['time']:.2f}", bounds)
    save(item)

    item = page(
        "08_box_budgets",
        "2D box: energy reservoirs and applied SGS transfers",
        "On/off starts have identical resolved fields; SGS-on adds its stated initial K. "
        "Unforced represented energy E+K should stay constant.",
        (2, 3),
        (14, 8),
    )
    ax = item[1]
    for n in (16, 32):
        for enabled in (0, 1):
            name = f"box_n{n}_sgs{enabled}"
            label = f"{n}² SGS {'on' if enabled else 'off'}"
            h0 = histories[name][0]
            curve(ax[0], name, label, kinetic)
            curve(ax[1], name, label, lambda q: thermal(q) - thermal(h0))
            curve(ax[5], name, label, lambda q: represented(q) - represented(h0))
            if enabled:
                curve(ax[2], name, label, lambda q: q["sgs-E"])
                curve(ax[3], name, label, lambda q: q["sgs-production-integral"])
                curve(ax[4], name, label, lambda q: q["sgs-dissipation-integral"])
    for a, title in zip(
        ax,
        (
            "Resolved kinetic",
            "Thermal gain",
            "SGS K",
            "Applied production",
            "Applied dissipation to heat",
            "E+K drift",
        ),
    ):
        a.set_title(title)
    save(item)

    item = page(
        "09_conservation",
        "Independent native-state inventories: numerical residuals",
        "Dots use checkpoint active cells, with no ghost cells. Transfer residual "
        "compares K change to the recorded applied production minus dissipation.",
    )
    ax = item[1]
    maxima = dict(energy=0.0, mass=0.0, transfer=0.0, scalar=0.0)
    for name in (
        "decay_cfl0.1",
        "shear_n32_sgs1",
        "diffusion_n32",
        "contact_n32",
        "box_n32_sgs1",
    ):
        states = snapshots[name]
        first = states[0]
        ts = [s["time"] for s in states]
        energy = [
            s["totals"][4] + s["totals"][-1] - first["totals"][4] - first["totals"][-1]
            for s in states
        ]
        mass = [s["totals"][0] - first["totals"][0] for s in states]
        transfer = [
            s["totals"][-1]
            - first["totals"][-1]
            - float(s["config"].get("sgs/integrated_production", "0"))
            + float(s["config"].get("sgs/integrated_dissipation", "0"))
            for s in states
        ]
        scalar = [s["totals"][5] - first["totals"][5] for s in states]
        for a, values, key in zip(ax, (energy, mass, transfer, scalar), maxima):
            a.plot(ts, values, ".-", label=name)
            maxima[key] = max(maxima[key], float(np.max(abs(np.array(values)))))
    for a, title in zip(
        ax,
        (
            "E+K inventory drift",
            "Mass inventory drift",
            "K − K₀ − production + dissipation",
            "Passive scalar mass drift",
        ),
    ):
        a.set_title(title)
        a.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    save(item)

    item = page(
        "10_restart",
        "2D restart and MeshBlock decomposition checks",
        "Restart at an intermediate normal checkpoint; compare at matched final time. "
        "Maps use a minimum color range of ±1e−16 to display exact zeros.",
    )
    ax = item[1]
    original = fields(snapshots["box_n32_sgs1"][-1])
    resumed = fields(snapshots["box_n32_restart"][-1])
    image(
        ax[0],
        resumed["rho"] - original["rho"],
        "Restart − uninterrupted: density",
        difference=True,
    )
    image(
        ax[1], resumed["K"] - original["K"], "Restart − uninterrupted: K", difference=True
    )
    a = snapshots["box_n32_sgs1"][-1]
    differences = [
        verify.state_diff(a, snapshots[name][-1])
        for name in ("box_n32_restart", "box_n32_single")
    ]
    ax[2].bar(["Restart", "Single block"], differences)
    ax[2].axhline(
        verify.LIMITS["checkpoint_restart_scaled_linf"],
        color="r",
        ls="--",
        label="Acceptance bound 1e−12",
    )
    ax[2].set(title="Scaled maximum field error", xlabel="Comparison", ylim=(0, 1.2e-12))
    zero = snapshots["zero_n32"]
    ax[3].plot(
        [s["time"] for s in zero],
        [s["totals"][-1] for s in zero],
        "o-",
        label="Exactly zero K",
    )
    ax[3].set(title="Zero-K state stays zero")
    save(item)

    item = page(
        "11_forced_budgets",
        "2D driven boxes: actual injection and energy budgets",
        "Seed 42, forcing modes 1–3, kz=0, tcorr=0.1, nominal dedt=0.1, t=1. "
        "Three velocity components are retained. These are short pilots, "
        "not stationary turbulence.",
        (2, 4),
        (17, 8),
    )
    ax = item[1]
    for row, label in enumerate(("sol", "comp")):
        for enabled in (0, 1):
            name = f"forced_{label}_sgs{enabled}"
            h0 = histories[name][0]
            text = "SGS on" if enabled else "SGS off"
            curve(ax[4 * row], name, text + " kinetic", kinetic)
            curve(
                ax[4 * row + 1],
                name,
                text + " thermal gain",
                lambda q: thermal(q) - thermal(h0),
            )
            curve(
                ax[4 * row + 2], name, text + " applied work", lambda q: q["force-work"]
            )
            curve(
                ax[4 * row + 3],
                name,
                text,
                lambda q: represented(q) - represented(h0) - q["force-work"],
            )
            if enabled:
                curve(ax[4 * row], name, "SGS reservoir K", lambda q: q["sgs-E"])
                curve(
                    ax[4 * row + 1],
                    name,
                    "SGS dissipation",
                    lambda q: q["sgs-dissipation-integral"],
                )
        for a, title in zip(
            ax[4 * row : 4 * row + 4],
            (
                "Kinetic + SGS",
                "Heating + dissipation",
                "Measured work",
                "E+K change − applied work",
            ),
        ):
            a.set_title(label + ": " + title)
    save(item)

    item = page(
        "12_forced_fields",
        "2D driven SGS-on fields at t=1",
        "Rows: solenoidal and compressive forcing. Each map has its own color scale.",
        (2, 4),
        (16, 8),
    )
    ax = item[1]
    for row, label in enumerate(("sol", "comp")):
        f = fields(snapshots[f"forced_{label}_sgs1"][-1])
        for j, key in enumerate(("rho", "K", "speed", "omega")):
            image(ax[4 * row + j], f[key], label + ": " + key)
    save(item)

    item = page(
        "13_forced_on_off",
        "2D driven fields: SGS-on minus SGS-off",
        "Identical forcing parameters and seed. Each lane normalizes its applied "
        "acceleration using its evolving state, "
        "so final acceleration/work need not match.",
        (2, 3),
        (14, 8),
    )
    ax = item[1]
    for row, label in enumerate(("sol", "comp")):
        on = fields(snapshots[f"forced_{label}_sgs1"][-1])
        off = fields(snapshots[f"forced_{label}_sgs0"][-1])
        for j, key in enumerate(("rho", "speed", "omega")):
            image(
                ax[3 * row + j], on[key] - off[key], label + ": Δ" + key, difference=True
            )
    save(item)

    item = page(
        "14_spectra",
        "2D shell spectra: resolved density-weighted velocity fluctuations",
        "Sum over all shells equals fluctuation kinetic energy (Parseval). "
        "No inertial-range fit or physical closure-accuracy claim from 32² pilots.",
        (1, 3),
        (15, 5),
    )
    ax = item[1]
    for a, prefix, title in zip(
        ax,
        ("box_n32_", "forced_sol_", "forced_comp_"),
        ("Decaying box", "Solenoidal driving", "Compressive driving"),
    ):
        for enabled in (0, 1):
            f = fields(snapshots[prefix + f"sgs{enabled}"][-1])
            power = spectrum(f["rho"], f["velocity"])
            a.semilogy(
                np.arange(1, len(power)),
                np.maximum(power[1:], 1e-30),
                "o-",
                ms=3,
                label="SGS " + ("on" if enabled else "off"),
            )
        a.set(
            title=title, xlabel="Radial mode number |k|", ylabel="Shell-integrated energy"
        )
    save(item)

    item = page(
        "15_cost_and_health",
        "Laptop execution cost and solver health",
        "One serial Release measurement per case. Solver elapsed time excludes plotting "
        "and Python inventory parsing; output and timing variability affect these short runs.",
        (2, 2),
        (14, 9),
    )
    ax = item[1]
    names = [
        "shear_n32_sgs0",
        "shear_n32_sgs1",
        "box_n32_sgs0",
        "box_n32_sgs1",
        "forced_sol_sgs0",
        "forced_sol_sgs1",
        "forced_comp_sgs0",
        "forced_comp_sgs1",
    ]
    cpu = []
    throughput = []
    for name in names:
        log = (output / name / "run.log").read_text()
        cpu.append(float(re.search(r"cpu time used\s*=\s*([\deE+.-]+)", log)[1]))
        throughput.append(
            float(re.search(r"zone-cycles/cpu_second\s*=\s*([\deE+.-]+)", log)[1])
        )
    ax[0].barh(names, cpu)
    ax[0].set(xlabel="Solver elapsed seconds", title="Reported solver elapsed time")
    ax[1].barh(names, throughput)
    ax[1].set(xlabel="Zone-cycles / elapsed second", title="Measured throughput")
    for name in ("box_n32_sgs1", "forced_sol_sgs1", "forced_comp_sgs1"):
        states = snapshots[name]
        ax[2].plot(
            [s["time"] for s in states], [s["min_sgs"] for s in states], label=name
        )
        curve(ax[3], name, name + " fallback", lambda q: q["sgs-fallback-stages"])
    ax[2].set(title="Minimum stored K (must be ≥0)")
    ax[3].set(title="Fallback stage count")
    save(item)
    pdf.close()

    summary = dict(
        simulation_runs=len(runs),
        result_records=len(rows),
        passed=sum(r["status"] == "pass" for r in rows),
        figures=len(gallery),
        panels=panel_count,
        contact_l1_errors=contact_errors,
        conservation_max_absolute=maxima,
        restart_and_decomposition_scaled_linf=differences,
        solver_elapsed_seconds=dict(zip(names, cpu)),
        qualification="Small 2D serial CPU checks; "
        "not stationary or reference-validated turbulence",
    )
    (output / "plot-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    links = "".join(
        f'<li><a href="#{name}">{html.escape(caption)}</a></li>'
        for name, caption in gallery
    )
    images = "".join(
        f'<section id="{name}"><h2>{name}</h2><p>{html.escape(caption)}</p>'
        f'<a href="plots/{name}.png"><img src="plots/{name}.png" '
        f'alt="{html.escape(caption)}"></a></section>'
        for name, caption in gallery
    )
    table = "".join(
        "<tr><td>"
        + (
            f'<a href="{html.escape(r["name"])}/input.athinput">'
            f"{html.escape(r['name'])}</a>"
            if "command" in r
            else html.escape(r["name"])
        )
        + "</td>"
        f"<td>{r['status']}</td><td>{r.get('cycles', '')}</td>"
        f"<td>{r.get('seconds', 0):.3f}</td></tr>"
        for r in rows
    )
    (output / "index.html").write_text(
        '<!doctype html><meta charset="utf-8">'
        "<title>AthenaK SGS laptop 2D report</title><style>"
        "body{max-width:1200px;margin:40px auto;padding:0 20px;"
        "font:17px/1.5 system-ui;color:#172b40}"
        "img{width:100%}section{margin:50px 0}td,th{padding:5px 14px;text-align:left}"
        "a{color:#1265a8}pre{white-space:pre-wrap}</style>"
        "<h1>AthenaK SGS: local 2D checks</h1>"
        f"<p>{summary['passed']}/{len(rows)} result records pass; "
        f"{len(runs)} simulations; "
        f"{len(gallery)} figure pages. Grids: 16² and 32². No MPI or GPU.</p>"
        "<p>Short deterministic checks and driven pilots. "
        "Resolution changes the SGS filter width. "
        "These plots do not establish stationary turbulence "
        "or closure accuracy against a physical reference.</p>"
        '<p><a href="report.pdf">Download the full PDF</a> · '
        '<a href="results.json">Numerical acceptance results</a> · '
        '<a href="manifest.json">Source and executable provenance</a></p>'
        "<h2>Plots</h2><ol>"
        + links
        + "</ol>"
        + images
        + "<h2>Executed runs and checks</h2><p>Seconds include "
        "process execution and Python state parsing.</p>"
        "<table><tr><th>Case</th><th>Result</th><th>Cycles</th><th>Seconds</th></tr>"
        + table
        + "</table>"
        "<h2>Measured summary</h2><pre>"
        + html.escape(json.dumps(summary, indent=2))
        + "</pre>"
    )
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return 0
    if args.output is None or (not args.plot_only and args.exe is None):
        parser.error("--output and (unless --plot-only) --exe are required")
    output = args.output.resolve()
    started = time.monotonic()
    if not args.plot_only:
        rows = run_cases(args.exe, output)
        if any(r["status"] != "pass" for r in rows):
            print(
                "One or more checks failed; inspect results.json and retained run logs."
            )
            return 1
    summary = make_plots(output)
    summary["workflow_wall_seconds"] = time.monotonic() - started
    (output / ("plot-execution.json" if args.plot_only else "execution.json")).write_text(
        json.dumps(
            dict(
                workflow_wall_seconds=summary["workflow_wall_seconds"],
                plot_only=args.plot_only,
            ),
            indent=2,
        )
        + "\n"
    )
    print(json.dumps(summary, indent=2))
    print("Open", output / "index.html")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
