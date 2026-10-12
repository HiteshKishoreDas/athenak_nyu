#!/usr/bin/env python3
"""Unforced 3D Mach-9 decay with matched Favre-filtered SGS on/off runs.

Local low-mode random initial conditions, not the imported Kritsuk dataset.
gamma=1.001 approximates isothermality. References are numerical, not proven DNS.
Reuse the existing native initializer, block filter, binary reader and ledgers.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import verify
from compare_resolution import restrict, spectrum, write_state
from heating_plots import budget

sys.path.insert(0, str(verify.SOURCE / "vis/python"))
from bin_convert_new import read_binary

GAMMA = 1.001
plt.switch_backend("Agg")


def initial(n, seed):
    """Real periodic modes 1..4, transverse/longitudinal energy ratio 2:1."""
    rng = np.random.default_rng(seed)
    base = np.fft.fftn(rng.normal(size=(16, 16, 16, 3)), axes=(0, 1, 2))
    index = (np.fft.fftfreq(n) * n).astype(int) % 16
    modes = base[np.ix_(index, index, index)].copy()
    kz, ky, kx = np.meshgrid(*(np.fft.fftfreq(n) * n,) * 3, indexing="ij")
    wave = np.stack((kx, ky, kz), axis=-1)
    norm = np.sum(wave**2, axis=-1)
    dot = np.sum(modes * wave, axis=-1)
    longitudinal = (
        wave * np.divide(dot, norm, out=np.zeros_like(dot), where=norm > 0)[..., None]
    )
    transverse = modes - longitudinal
    keep = (norm >= 1) & (norm <= 16)
    phase = np.exp(1j * np.pi * (kx + ky + kz) / n)
    components = [
        np.fft.ifftn(x * keep[..., None] * phase[..., None], axes=(0, 1, 2)).real
        for x in (transverse, longitudinal)
    ]
    u = (
        sum(
            x * np.sqrt(weight / np.mean(np.sum(x * x, axis=-1)))
            for x, weight in zip(components, (2 / 3, 1 / 3))
        )
        * 9
    )
    data = np.zeros((n, n, n, 7))
    data[..., 0] = 1
    data[..., 1:4] = u
    data[..., 4] = 1 / (GAMMA * (GAMMA - 1)) + 0.5 * np.sum(u * u, axis=-1)
    data[..., 5:7] = (0.2, 0.25)
    return data


def dense(path):
    """Assemble unrefined serial double binary output; reject missing/duplicate cells."""
    raw = read_binary(str(path))
    config = verify.parameters("\n".join(raw["header"]))
    n = raw["Nx1"]
    names = ["dens", "mom1", "mom2", "mom3", "ener", "sgs_energy"]
    names = [v for v in names if v in raw["var_names"]]
    if (
        names[:5] != ["dens", "mom1", "mom2", "mom3", "ener"]
        or raw["Nx2"] != n
        or raw["Nx3"] != n
        or config.get("mesh/refinement", "none") != "none"
        or raw["mb_data"][names[0]][0].dtype != np.float64
    ):
        raise ValueError("require unrefined cubic hydro_u output")
    data = np.empty((n, n, n, len(names)))
    occupied = np.zeros((n, n, n), dtype=bool)
    for m, loc in enumerate(raw["mb_logical"]):
        if loc[3] != raw["mb_logical"][0, 3]:
            raise ValueError("mixed refinement levels")
        shape = raw["mb_data"][names[0]][m].shape
        start = loc[2::-1] * np.array(shape)
        sl = tuple(slice(int(a), int(a + b)) for a, b in zip(start, shape))
        if any(a < 0 or a + b > n for a, b in zip(start, shape)) or occupied[sl].any():
            raise ValueError("duplicate/out-of-domain blocks")
        occupied[sl] = True
        for v, name in enumerate(names):
            data[sl + (v,)] = raw["mb_data"][name][m]
    if not occupied.all() or not np.isfinite(data).all():
        raise ValueError("incomplete or nonfinite volume")
    return raw["time"], data


def inventory(data):
    rho = data[..., 0]
    kinetic = np.sum(data[..., 1:4] ** 2, axis=-1) / (2 * rho)
    thermal = data[..., 4] - kinetic
    K = data[..., 5] if data.shape[-1] == 6 else np.zeros_like(rho)
    cs2 = GAMMA * (GAMMA - 1) * thermal / rho
    if rho.min() <= 0 or thermal.min() <= 0 or K.min() < 0:
        raise ValueError("inadmissible reservoir")
    return dict(
        mass=float(rho.mean()),
        momentum=data[..., 1:4].mean((0, 1, 2)).tolist(),
        kinetic=float(kinetic.mean()),
        thermal=float(thermal.mean()),
        K=float(K.mean()),
        represented=float(np.mean(data[..., 4] + K)),
        rho_min=float(rho.min()),
        thermal_min=float(thermal.min()),
        K_min=float(K.min()),
        cs_rms=float(np.sqrt(cs2.mean())),
        mass_weighted_Mach=float(np.sqrt(np.mean(2 * kinetic / cs2) / rho.mean())),
    )


def run(args):
    suite = verify.Verification(args.exe, args.output)

    def check(row, key, value, limit, description):
        return suite.check(row, key, float(value), limit, description)

    fine = initial(args.high, args.seed)[..., : 5 + args.nscalars].copy()
    T = 1 / 18
    protocol = dict(
        high=args.high,
        reference_check=args.high // 2,
        low=args.low,
        seed=args.seed,
        gamma=GAMMA,
        cs0=1,
        Mach0=9,
        T=T,
        end_T=args.end_T,
        snapshot_dt_T=0.5,
        history_dt_T=0.1,
        cfl=args.cfl,
        block_max=args.block,
        nscalars=args.nscalars,
        checkpoint_dt_T=1 if args.high >= 256 else None,
        threads=int(os.environ.get("OMP_NUM_THREADS", "1")),
        reference_convergence_screen=0.1,
        interpretation_window_T=[1, args.end_T],
        initial_conditions="rho=1; Gaussian random modes1..4, 2:1 solenoidal energy",
        initialization="conservative restriction; exact variance to SGS or heat",
        forcing="absent throughout; no turb_driving block",
        filter="aligned finite-volume/Favre block averages",
        paper="https://arxiv.org/pdf/2410.23339v3",
        departures=[
            "local undeveloped initial field, not Kritsuk imported data",
            "near-isothermal ideal EOS instead of exact isothermal",
            "smaller reference: convergence must be measured",
            "AthenaK RK2/LLF/PPM4; original SGS coefficients unchanged",
        ],
        timeout_seconds=args.timeout,
        maximum_cycles=50000,
    )
    (suite.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    manifest = json.loads((suite.output / "manifest.json").read_text())
    manifest["decay_harness_sha256"] = verify.sha(__file__)
    (suite.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    reused = {}
    if args.reference_output:
        source = args.reference_output.resolve()
        old_protocol = json.loads((source / "protocol.json").read_text())
        for key in protocol.keys() - {"low", "departures", "timeout_seconds"}:
            if old_protocol.get(key) != protocol[key]:
                raise ValueError(f"reference protocol mismatch: {key}")
        old_manifest = json.loads((source / "manifest.json").read_text())
        for key in ("executable_sha256", "source_sha256"):
            if old_manifest[key] != manifest[key]:
                raise ValueError(f"reference provenance mismatch: {key}")
        old_rows = {
            r["name"]: r for r in json.loads((source / "results.json").read_text())
        }
        aliases = {"reference": "reference", "reference_check": "reference_check"}
        if args.high // 2 in args.low:
            aliases[f"n{args.high // 2}_off"] = "reference_check"
        template = verify.parameters(
            (verify.SOURCE / "inputs/hydro/sgs_verify.athinput").read_text()
        )
        for name, original in aliases.items():
            row = old_rows[original]
            wd = source / original
            if row["status"] != "pass" or not all(c["passed"] for c in row["checks"]):
                raise ValueError(f"reference not passed: {original}")
            for filename, key in (
                ("input.athinput", "input_sha256"),
                ("initial.dat", "initial_sha256"),
            ):
                if verify.sha(wd / filename) != row[key]:
                    raise ValueError(f"reference file changed: {wd / filename}")
            config = verify.parameters((wd / "input.athinput").read_text())
            if any(
                config.get(k) != v
                for k, v in template.items()
                if k.startswith("sgs/") and k != "sgs/enabled"
            ):
                raise ValueError(
                    "reference SGS coefficients/settings differ from current template"
                )
            reused[name] = row
        protocol["reference_output"] = str(source)
        (suite.output / "protocol.json").write_text(
            json.dumps(protocol, indent=2) + "\n"
        )
    cases = [
        ("reference", args.high, False),
        ("reference_check", args.high // 2, False),
    ]
    cases += [
        (f"n{n}_{lane}", n, lane == "on") for n in args.low for lane in ("off", "on")
    ]
    for name, n, enabled in cases:
        original = reused.get(name)
        wd = (
            args.reference_output.resolve() / original["name"]
            if original
            else suite.output / name
        )
        if not original:
            wd.mkdir()
        data, K = restrict(fine, n)
        expected = (
            None if original else write_state(wd / "initial.dat", data, K, enabled)
        )
        config = verify.parameters(
            (verify.SOURCE / "inputs/hydro/sgs_verify.athinput").read_text()
        )
        config.update(
            {
                "job/basename": name,
                "mesh/nghost": "3",
                **{f"mesh/nx{d}": str(n) for d in (1, 2, 3)},
                **{f"meshblock/nx{d}": str(min(n, args.block)) for d in (1, 2, 3)},
                "hydro/gamma": str(GAMMA),
                "hydro/nscalars": str(args.nscalars),
                "hydro/reconstruct": "ppm4",
                "time/cfl_number": str(args.cfl),
                "time/tlim": str(args.end_T * T),
                "time/nlim": "50000",
                "time/ndiag": "1000",
                "sgs/enabled": str(enabled).lower(),
                "problem/initial_state": str(wd / "initial.dat"),
                "output1/dt": str(0.1 * T),
                "output1/data_format": "%24.17e",
                "output2/dt": str(0.5 * T),
                "output2/data_precision": "real",
                "output3/dt": "-1",
                # Checkpoint the long256³ campaign; smaller pilots keep their IO low.
                "output4/dt": str(T) if args.high >= 256 else "-1",
            }
        )
        inp = wd / "input.athinput"
        if not original:
            inp.write_text(verify.input_text(config))
        command = (
            original["command"]
            if original
            else [str(args.exe.resolve()), "-i", str(inp)]
        )
        row = dict(
            name=name,
            n=n,
            enabled=enabled,
            command=command,
            status="running",
            input_sha256=verify.sha(inp),
            initial_sha256=verify.sha(wd / "initial.dat"),
        )
        if original:
            row.update(
                reused_from=str(wd), original_wall_seconds=original["wall_seconds"]
            )
        suite.results.append(row)
        suite.save()
        started = time.monotonic()
        try:
            if not original:
                with (wd / "run.log").open("w") as log:
                    subprocess.run(
                        command,
                        cwd=wd,
                        stdout=log,
                        stderr=subprocess.STDOUT,
                        check=True,
                        timeout=args.timeout,
                    )
            hist = verify.history(
                wd / f"{original['name'] if original else name}.hydro.hst"
            )
            if original and hist != original["history"]:
                raise ValueError("reference history changed")
            samples = []
            paths = sorted((wd / "bin").glob("*.bin"))
            if original and [str(p) for p in paths] != [
                s["path"] for s in original["samples"]
            ]:
                raise ValueError("reference snapshot list changed")
            for index, path in enumerate(paths):
                digest = verify.sha(path)
                if original and digest != original["samples"][index]["sha256"]:
                    raise ValueError(f"reference snapshot changed: {path}")
                t, loaded = dense(path)
                values = inventory(loaded)
                values.update(
                    time=t, binary_header_time=t, path=str(path), sha256=digest
                )
                if index == len(paths) - 1:
                    check(
                        row,
                        "final_binary_time",
                        abs(t - hist[-1]["time"]),
                        1e-6,
                        "binary ASCII timestamp is rounded to six significant digits",
                    )
                    values["time"] = hist[-1]["time"]
                if index == 0:
                    target = (
                        np.concatenate((expected[..., :5], expected[..., -1:]), -1)
                        if enabled
                        else data[..., :5]
                    )
                    check(
                        row,
                        "initial_roundtrip",
                        np.max(abs(loaded - target)),
                        0,
                        "native double file exactly equals prescribed initial state",
                    )
                    origin = values
                for key, limit in [("mass", 1e-11), ("represented", 1e-10)]:
                    check(
                        row,
                        f"{index}.{key}",
                        verify.scaled_diff(values[key], origin[key]),
                        limit,
                        "periodic unforced conservation",
                    )
                check(
                    row,
                    f"{index}.momentum",
                    np.max(abs(np.array(values["momentum"]) - origin["momentum"])),
                    1e-10,
                    "periodic unforced momentum",
                )
                if enabled and index == len(paths) - 1:
                    for label, value in [
                        ("kinetic-E", values["kinetic"]),
                        ("thermal-E", values["thermal"]),
                        ("sgs-E", values["K"]),
                    ]:
                        hvalue = (
                            sum(hist[-1][f"{v}-KE"] for v in (1, 2, 3))
                            if label == "kinetic-E"
                            else hist[-1][label]
                        )
                        check(
                            row,
                            "history." + label,
                            verify.scaled_diff(value, hvalue),
                            1e-10,
                            "independent binary inventory vs final history",
                        )
                if name in ("reference", "reference_check"):
                    for low in args.low:
                        coarse, variance = restrict(loaded, low)
                        values[f"K_ref_n{low}"] = float(variance.mean())
                        values[f"resolved_n{low}"] = inventory(coarse)["kinetic"]
                        check(
                            row,
                            f"{index}.kinetic_decomposition_n{low}",
                            verify.scaled_diff(
                                values["kinetic"],
                                values[f"resolved_n{low}"] + variance.mean(),
                            ),
                            1e-11,
                            "Favre resolved plus variance equals fine kinetic energy",
                        )
                samples.append(values)
            check(
                row,
                "completed_time",
                abs(hist[-1]["time"] - args.end_T * T),
                1e-12,
                "reached end time; cycle cap is not completion",
            )
            if enabled:
                final = hist[-1]
                applied = (
                    final["sgs-production-integral"] - final["sgs-dissipation-integral"]
                )
                check(
                    row,
                    "paired_transfer",
                    verify.scaled_diff(samples[-1]["K"] - samples[0]["K"], applied),
                    1e-8,
                    "K change equals applied production minus dissipation",
                )
            row.update(
                status="pass" if all(c["passed"] for c in row["checks"]) else "fail",
                samples=samples,
                history=hist,
                wall_seconds=time.monotonic() - started,
            )
        except Exception as exc:
            row.update(
                status="fail", error=str(exc), wall_seconds=time.monotonic() - started
            )
        suite.save()
        print(name, row["status"], row.get("error", ""), flush=True)
        if row["status"] != "pass":
            return False
    return True


def sample_times(stop):
    return np.unique(np.r_[np.arange(0, stop, 0.5), 1.0, stop])


def report(output):
    rows = json.loads((output / "results.json").read_text())
    runs = {r["name"]: r for r in rows}
    protocol = json.loads((output / "protocol.json").read_text())
    T = protocol["T"]
    times = sample_times(protocol["end_T"])

    def series(name, key):
        s = runs[name]["samples"]
        return np.interp(times, [x["time"] / T for x in s], [x[key] for x in s])

    pdf = PdfPages(output / "report.pdf")
    scores = {}

    def save(name, fig):
        fig.tight_layout()
        fig.savefig(output / f"{name}.png", dpi=160)
        pdf.savefig(fig)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for i, n in enumerate(protocol["low"]):
        ref = series("reference", f"K_ref_n{n}")
        check = series("reference_check", f"K_ref_n{n}")
        model = series(f"n{n}_on", "K")
        mask = times >= 1
        finer_check = protocol["reference_check"] > n
        sensitivity = abs(ref - check) / np.maximum(ref, check) if finer_check else None
        scores[str(n)] = dict(
            integrated_K_ratio_1_to_end=float(
                np.trapezoid(model[mask], times[mask])
                / np.trapezoid(ref[mask], times[mask])
            ),
            final_K_ratio=float(model[-1] / ref[-1]),
            maximum_reference_sensitivity_1_to_end=float(sensitivity[mask].max())
            if finer_check
            else None,
            reference_screen_pass=bool(sensitivity[mask].max() <= 0.1)
            if finer_check
            else None,
            reference_check_finer_than_coarse=finer_check,
        )
        for curve, label in [
            (ref, f"filtered {protocol['high']}³ reference"),
            *(
                [(check, f"filtered {protocol['reference_check']}³ reference")]
                if finer_check
                else []
            ),
            (model, "evolved SGS"),
        ]:
            axes[0, i].plot(times, curve, label=label)
        axes[0, i].set(title=f"{n}³ SGS energy density", ylabel="Volume mean K")
        axes[1, i].plot(times, model / ref, label="SGS / high reference")
        if finer_check:
            axes[1, i].plot(
                times,
                check / ref,
                label=f"filtered {protocol['reference_check']}³ / {protocol['high']}³",
            )
        else:
            axes[1, i].text(
                0.04,
                0.08,
                "No independent finer reference pair",
                transform=axes[1, i].transAxes,
            )
        axes[1, i].axhline(1, color="black", linestyle="--")
        axes[1, i].set(ylabel="Energy ratio", title="Initial equality is imposed")
    for ax in axes.flat:
        ax.set(xlabel="t / initial crossing time T")
        ax.axvspan(0, 1, color="grey", alpha=0.15)
        ax.legend(fontsize=8)
    fig.suptitle(
        f"Unforced Mach-9 decay — {protocol['high']}³ / {protocol['reference_check']}³ references",
        fontsize=14,
    )
    save("energy-comparison", fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for name in runs:
        kinetic = series(name, "kinetic") + series(name, "K")
        axes[0, 0].plot(times, kinetic, label=name)
        axes[0, 1].plot(times, series(name, "mass_weighted_Mach"), label=name)
        represented = series(name, "represented")
        axes[1, 0].plot(
            times, (represented - represented[0]) / represented[0], label=name
        )
        axes[1, 1].plot(times, series(name, "cs_rms"), label=name)
    for ax, title in zip(
        axes.flat,
        [
            "Resolved + SGS kinetic energy",
            "Resolved sonic Mach",
            "Represented-energy relative drift",
            "Sound speed: quasi-isothermal check",
        ],
    ):
        ax.set(title=title, xlabel="t / T")
        ax.legend(fontsize=7)
    save("decay-and-conservation", fig)
    fig, axes = plt.subplots(3, 2, figsize=(12, 11))
    for i, n in enumerate(protocol["low"]):
        h = runs[f"n{n}_on"]["history"]
        ht = np.array([v["time"] / T for v in h])
        for key in [
            "sgs-production-integral",
            "sgs-dissipation-integral",
            "sgs-pressure-work-integral",
            "sgs-stress-limited-integral",
        ]:
            axes[0, i].plot(ht, [v[key] for v in h], label=key.replace("sgs-", ""))
        axes[0, i].set(title=f"{n}³ applied cumulative transfers")
        source = np.array(
            [v["sgs-production-integral"] - v["sgs-dissipation-integral"] for v in h]
        )
        reservoir = np.array([v["sgs-E"] - h[0]["sgs-E"] for v in h])
        axes[1, i].plot(ht, reservoir - source)
        axes[1, i].set(
            title="K change minus production–dissipation", ylabel="Ledger residual"
        )
        augmented = [dict(v, **{"force-work": 0.0}) for v in h]
        bt, transfers, summary = budget(augmented, T, protocol["end_T"] * T)
        scores[str(n)]["heating_window_1_to_end"] = summary
        axes[2, i].plot(bt / T, transfers["sgs"], label="SGS dissipation to heat")
        axes[2, i].plot(bt / T, transfers["numerical"], label="Numerical/RK residual")
        axes[2, i].set(
            title="Heating after thermal pressure work", ylabel="Energy since 1 T"
        )
    for ax in axes.flat:
        ax.set(xlabel="t / T")
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=8)
    save("source-budgets", fig)
    # Local statistics compare distributions, not chaotic pointwise trajectories.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    bins = np.logspace(-3, 3, 30)
    for i, n in enumerate(protocol["low"]):
        for name, label in [("reference", "filtered reference"), (f"n{n}_on", "SGS")]:
            _, data = dense(Path(runs[name]["samples"][-1]["path"]))
            if name == "reference":
                data, K = restrict(data, n)
            else:
                K = data[..., -1]
            rho = data[..., 0]
            sigma = np.sqrt(2 * K / rho)
            median = [
                float(np.median(sigma[(rho >= a) & (rho < b)]))
                if ((rho >= a) & (rho < b)).any()
                else np.nan
                for a, b in zip(bins[:-1], bins[1:])
            ]
            axes[0, i].loglog(np.sqrt(bins[:-1] * bins[1:]), median, label=label)
            u = data[..., 1:4] / rho[..., None]
            p = spectrum(rho, u)
            axes[1, i].loglog(np.arange(1, len(p)), p[1:], label=label)
        _, off = dense(Path(runs[f"n{n}_off"]["samples"][-1]["path"]))
        off_power = spectrum(off[..., 0], off[..., 1:4] / off[..., :1])
        axes[1, i].loglog(np.arange(1, len(off_power)), off_power[1:], label="SGS off")
        axes[0, i].set(
            title=f"{n}³ final density–SGS velocity",
            xlabel="rho / mean rho",
            ylabel="Median sigma_sgs (3D)",
        )
        axes[1, i].set(
            title="Final resolved kinetic spectrum",
            xlabel="Shell k L / (2 pi)",
            ylabel="Shell energy of sqrt(rho) u",
        )
    for ax in axes.flat:
        ax.legend()
    save("local-statistics", fig)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    n = max(protocol["low"])
    maps = []
    for j, (name, label) in enumerate(
        [
            ("reference", "Filtered reference"),
            (f"n{n}_on", "SGS on"),
            (f"n{n}_off", "SGS off"),
        ]
    ):
        _, data = dense(Path(runs[name]["samples"][-1]["path"]))
        if name == "reference":
            data, K = restrict(data, n)
        else:
            K = data[..., -1] if name.endswith("on") else np.zeros(data.shape[:3])
        maps.append((label, data[..., 0], K))
    kmax = max(float(K.max()) for _, _, K in maps)
    for j, (label, rho, K) in enumerate(maps):
        for row, field, title in [
            (0, np.log10(rho), "log10 density"),
            (1, K, "SGS energy density"),
        ]:
            image = axes[row, j].imshow(
                field[n // 2],
                origin="lower",
                extent=(0, 1, 0, 1),
                cmap="viridis",
                vmin=-2 if row == 0 else 0,
                vmax=2 if row == 0 else kmax,
            )
            axes[row, j].set(title=label + ": " + title, xlabel="x / L", ylabel="y / L")
            fig.colorbar(image, ax=axes[row, j], shrink=0.7)
    fig.suptitle(
        f"{n}³ final central-z slices — compare statistics, not pointwise trajectories"
    )
    save("final-slices", fig)
    pdf.close()
    summary = dict(
        protocol=protocol,
        scores=scores,
        analysis_sha256=verify.sha(__file__),
        numerical_checks=sum(len(r.get("checks", [])) for r in rows),
        failed_checks=sum(not c["passed"] for r in rows for c in r.get("checks", [])),
        runs={
            name: dict(
                status=r["status"],
                wall_seconds=r["wall_seconds"],
                reused_from=r.get("reused_from"),
                original_wall_seconds=r.get("original_wall_seconds"),
                final=r["samples"][-1],
                final_counters={
                    key: value
                    for key, value in r["history"][-1].items()
                    if key.startswith("sgs-")
                },
            )
            for name, r in runs.items()
        },
    )
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(scores, indent=2))
    return summary


def self_test():
    for stop in [1.1, 1.25, 10]:
        times = sample_times(stop)
        assert 1.0 in times and times[-1] == stop and sum(times >= 1) >= 2
    fine = initial(16, 42)
    assert abs(inventory(fine[..., :5])["kinetic"] - 40.5) < 1e-12
    assert abs(np.mean(fine[..., 1:4])) < 1e-14
    coarse, K = restrict(fine, 8)
    assert K.min() >= 0
    assert abs(inventory(coarse[..., :5])["kinetic"] + K.mean() - 40.5) < 1e-12
    print("Mach-9 normalization, zero momentum and Favre energy checks pass")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--exe", type=Path)
    p.add_argument(
        "--reference-output",
        type=Path,
        help="reuse passed matching references; verify hashes and recompute filters",
    )
    p.add_argument("--output", type=Path)
    p.add_argument("--high", type=int, choices=[32, 64, 128, 256], default=128)
    p.add_argument("--block", type=int, choices=[16, 32, 64], default=16)
    p.add_argument("--nscalars", type=int, choices=[0, 2], default=2)
    p.add_argument("--low", type=int, nargs=2, default=[8, 16])
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--end-T", type=float, default=10)
    p.add_argument("--cfl", type=float, default=0.2)
    p.add_argument("--timeout", type=float, default=3600)
    p.add_argument("--plot-only", action="store_true")
    p.add_argument("--self-test", action="store_true")
    args = p.parse_args()
    if args.self_test:
        self_test()
        return
    if args.output is None or (not args.plot_only and args.exe is None):
        p.error("require --output and --exe unless --plot-only")
    if (
        len(set(args.low)) != 2
        or any(n < 4 or args.high // 2 % n for n in args.low)
        or not 1 < args.end_T <= 10
        or not 0 < args.cfl <= 0.2
        or args.timeout <= 0
    ):
        p.error(
            "require two distinct coarse divisors >=4 of high/2, 1<end-T<=10, 0<cfl<=.2"
        )
    args.output = args.output.resolve()
    if not args.plot_only and not run(args):
        raise SystemExit("Run failed; retained results and logs identify the failure")
    report(args.output)


if __name__ == "__main__":
    main()
