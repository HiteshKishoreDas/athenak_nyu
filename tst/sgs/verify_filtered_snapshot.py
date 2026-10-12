#!/usr/bin/env python3
"""Paper Section 3 comparison on supplied SGS-off native periodic 3D snapshots.

No forcing parameter enters this comparison. Gaussian width uses the explicit
second-moment convention Delta_G=sqrt(12)*sigma, not a certified paper filter.
Exact covariance K is supplied to the actual C++ stress kernel. This tests the
production closure, not evolved SGS energy or the solver's face discretization.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from filter_reference import (error_metrics, favre_reference,
                              gaussian_filter, gaussian_kernel, mixed_stress,
                              production_metrics)
import verify


def native_stress(kernel, rho, K, delta, gradient):
    if sys.byteorder != "little":
        raise ValueError("native double interchange qualified on little-endian CPU")
    records = np.column_stack((rho.ravel(), K.ravel(),
                               np.full(rho.size, delta), gradient.reshape(-1, 9)))
    with tempfile.TemporaryDirectory(prefix="sgs-filtered-") as directory:
        inp, out = (Path(directory) / name for name in ("in.bin", "out.bin"))
        records.astype(np.float64).tofile(inp)
        subprocess.run([str(kernel), str(inp), str(out)], check=True, timeout=60)
        result = np.fromfile(out, dtype=np.float64)
    if result.size != 10 * rho.size or not np.isfinite(result).all():
        raise ValueError("invalid native stress result")
    result = result.reshape(rho.shape + (10,))
    return result[..., :9].reshape(gradient.shape), result[..., 9]


def direct_covariance(rho, velocity, width, kind, point):
    n = rho.shape[0]
    if kind == "gaussian":
        w = gaussian_kernel(n, width)
        indices = [(p - np.arange(n)) % n for p in point]
    else:
        width = int(width)
        offsets = np.arange(-(width // 2), width - width // 2)
        indices = [(p + offsets) % n for p in point]
        w = np.full(width, 1 / width)
    weights = w[:, None, None] * w[None, :, None] * w[None, None, :]
    density = rho[np.ix_(*indices)]
    u = velocity[np.ix_(*indices)]
    mass = weights * density
    mean = np.sum(mass[..., None] * u, axis=(0, 1, 2)) / mass.sum()
    fluctuation = u - mean
    return -np.einsum("xyz,xyzi,xyzj->ij", mass, fluctuation, fluctuation)


def spectral_gradient(velocity, dx):
    transformed = np.fft.fftn(velocity, axes=(0, 1, 2))
    n = velocity.shape[0]
    wavenumber = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    # The real derivative of the even-grid Nyquist mode is zero at sample points.
    if n % 2 == 0:
        wavenumber[n // 2] = 0
    gradient = np.empty(velocity.shape[:-1] + (3, 3))
    for axis in range(3):
        shape = [1, 1, 1, 1]
        shape[axis] = n
        gradient[..., :, axis] = np.fft.ifftn(
            transformed * (1j * wavenumber.reshape(shape)), axes=(0, 1, 2)).real
    return gradient


def assess(rho, velocity, delta, kind, kernel, detailed=False):
    n = rho.shape[0]
    width = delta * n
    if kind == "box":
        if abs(width - round(width)) > 1e-12:
            raise ValueError("box width must contain an integer number of fine cells")
        width = int(round(width))
    ref = favre_reference(rho, velocity, width, 1 / n, filter_kind=kind)
    if np.any(ref["K"] < 0):
        raise ValueError("negative covariance energy; reference is not clipped")
    tau, production = native_stress(kernel, ref["rho"], ref["K"], delta,
                                     ref["gradient"])
    python_tau = mixed_stress(ref["rho"], ref["K"], delta, ref["gradient"])
    checks = {}

    def check(name, actual, expected, tolerance=2e-11):
        scale = max(float(np.max(np.abs(expected))), 1)
        error = float(np.max(np.abs(np.asarray(actual) - expected))) / scale
        checks[name] = dict(scaled_error=error, tolerance=tolerance, passed=error <= tolerance)
        if error > tolerance:
            raise AssertionError((name, error, tolerance))

    check("native_vs_independent_Python_stress", tau, python_tau, 5e-13)
    check("native_production_contraction", production,
          np.einsum("...ij,...ij->...", python_tau, ref["gradient"]), 5e-13)
    check("mass_preservation", ref["rho"].mean(), rho.mean())
    check("momentum_preservation", np.mean(ref["rho"][..., None] * ref["velocity"],
                                           axis=(0, 1, 2)),
          np.mean(rho[..., None] * velocity, axis=(0, 1, 2)))
    kinetic = 0.5 * np.mean(rho * np.sum(velocity**2, axis=-1))
    resolved = 0.5 * np.mean(ref["rho"] * np.sum(ref["velocity"]**2, axis=-1))
    check("kinetic_energy_decomposition", resolved + ref["K"].mean(), kinetic)
    check("native_trace", np.trace(tau, axis1=-2, axis2=-1), -2 * ref["K"])
    check("exact_stress_symmetry", ref["tau"], ref["tau"].swapaxes(-1, -2))
    minimum_eigenvalue = float(np.linalg.eigvalsh(-ref["tau"]).min())
    check("exact_covariance_positive_semidefinite", min(minimum_eigenvalue, 0), 0)
    for point in [(0, 0, 0), (n // 3, n // 5, n // 7)]:
        check(f"direct_centered_covariance_{point}", ref["tau"][point],
              direct_covariance(rho, velocity, width, kind, point))
    identity = np.eye(3)
    exact_deviatoric = ref["tau"] + (2 / 3) * ref["K"][..., None, None] * identity
    model_deviatoric = tau + (2 / 3) * ref["K"][..., None, None] * identity
    metrics = dict(n=n, delta=delta, filter=kind, width_fine_cells=width,
                   sigma_fine_cells=float(width / np.sqrt(12)) if kind == "gaussian" else None,
                   undersampled_gaussian=kind == "gaussian" and width < 4,
                   K_reference_mean=float(ref["K"].mean()),
                   K_reference_min=float(ref["K"].min()),
                   resolved_kinetic_mean=float(resolved),
                   total_kinetic_mean=float(kinetic),
                   stress=error_metrics(ref["tau"], tau),
                   deviatoric_stress=error_metrics(exact_deviatoric, model_deviatoric),
                   production=production_metrics(ref["production"], production, ref["rho"]),
                   deviatoric_production=production_metrics(
                       np.einsum("...ij,...ij->...", exact_deviatoric, ref["gradient"]),
                       np.einsum("...ij,...ij->...", model_deviatoric, ref["gradient"]), ref["rho"]),
                   checks=checks)
    p = metrics["production"]
    p["mean_ratio"] = p["model_mean"] / p["exact_mean"] if p["exact_mean"] else None
    if detailed:
        boosted = favre_reference(rho, velocity + [1, -2, 0.5], width, 1 / n,
                                  filter_kind=kind)
        check("Galilean_stress", boosted["tau"], ref["tau"])
        check("Galilean_K", boosted["K"], ref["K"])
        check("Galilean_production", boosted["production"], ref["production"])
        spectral = spectral_gradient(ref["velocity"], 1 / n)
        spectral_tau, spectral_model = native_stress(kernel, ref["rho"], ref["K"],
                                                      delta, spectral)
        spectral_exact = np.einsum("...ij,...ij->...", ref["tau"], spectral)
        metrics["spectral_derivative_sensitivity"] = production_metrics(
            spectral_exact, spectral_model, ref["rho"])
    return metrics, ref, tau, production


def self_test(kernel):
    checks = {}
    n = 64
    phase = 2 * np.pi * np.arange(n) / n
    field = np.broadcast_to(np.cos(phase)[:, None, None], (n, n, n))
    for width in [4, 8, 16]:
        w = gaussian_kernel(n, width)
        assert np.all(w >= 0) and abs(w.sum() - 1) < 1e-14
        expected = np.exp(-0.5 * (2 * np.pi * width / (np.sqrt(12) * n))**2)
        error = float(np.max(abs(gaussian_filter(field, width) - expected * field)))
        assert error < 2e-11
        checks[f"Gaussian_mode1_width{width}"] = error
    gradient = spectral_gradient(np.stack((field, field * 0, field * 0), axis=-1), 1 / n)
    assert np.max(abs(gradient[..., 0, 0] + 2 * np.pi * np.sin(phase)[:, None, None])) < 1e-12
    rho = np.ones((1, 1, 1)) * 3
    K = np.ones_like(rho) * 0.04
    a = np.diag([0.2, -0.1, -0.1]).reshape((1, 1, 1, 3, 3))
    tau, production = native_stress(kernel, rho, K, 0.1, a)
    expected = mixed_stress(rho, K, 0.1, a)
    assert np.max(abs(tau - expected)) < 1e-14
    assert abs(float(production.item()) - np.sum(expected * a)) < 1e-14
    with tempfile.TemporaryDirectory(prefix="sgs-invalid-record-") as directory:
        inp, out = (Path(directory) / name for name in ("in.bin", "out.bin"))
        np.ones(11).tofile(inp)
        assert subprocess.run([str(kernel), str(inp), str(out)], timeout=60).returncode != 0
    return dict(status="pass", truncated_native_record="rejected", Fourier_checks=checks,
                native_hand_case="pass", spectral_gradient_axis="pass")


def plots(out, rows, fields):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm
    ref, tau, model, metrics = fields
    exact = ref["production"]
    scale = float(np.sqrt(np.mean(exact**2)))
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    axes[0, 0].hexbin(exact.ravel() / scale, model.ravel() / scale,
                      gridsize=100, bins="log", mincnt=1)
    limit = float(np.quantile(abs(np.concatenate((exact.ravel(), model.ravel()))) / scale, 0.995))
    axes[0, 0].plot([-limit, limit], [-limit, limit], "r--")
    axes[0, 0].set(xlim=(-limit, limit), ylim=(-limit, limit),
                   xlabel="Measured transfer / measured RMS", ylabel="Native model / measured RMS",
                   title="Joint distribution (axes show central 99.5%)")
    edges = np.linspace(-limit, limit, 150)
    for values, label in [(exact, "Measured"), (model, "Native closure")]:
        axes[0, 1].hist(values.ravel() / scale, edges, density=True, histtype="step", label=label)
    axes[0, 1].set(yscale="log", xlabel="Transfer / measured RMS", ylabel="PDF within plotted range",
                   title="Forward transfer and backscatter")
    axes[0, 1].legend()
    exact_dev = ref["tau"] + 2 / 3 * ref["K"][..., None, None] * np.eye(3)
    model_dev = tau + 2 / 3 * ref["K"][..., None, None] * np.eye(3)
    rms = np.sqrt(np.mean(exact_dev**2))
    axes[0, 2].hexbin(exact_dev.ravel() / rms, model_dev.ravel() / rms,
                      gridsize=100, bins="log", mincnt=1)
    axes[0, 2].plot([-5, 5], [-5, 5], "r--")
    axes[0, 2].set(xlim=(-5, 5), ylim=(-5, 5), xlabel="Measured trace-free stress / RMS",
                   ylabel="Native trace-free stress / RMS", title="All tensor components")
    z = exact.shape[2] // 2
    norm = SymLogNorm(linthresh=0.1, vmin=-limit, vmax=limit)
    for ax, values, title in zip(axes[1], [ref["K"], exact / scale, model / scale],
                                 ["Direct filtered unresolved energy", "Measured transfer", "Native modeled transfer"]):
        image = ax.imshow(values[:, :, z].T, origin="lower", extent=[0, 1, 0, 1],
                          cmap="viridis" if title.endswith("energy") else "RdBu_r",
                          norm=None if title.endswith("energy") else norm)
        ax.set(xlabel="x / box length", ylabel="y / box length", title=title)
        fig.colorbar(image, ax=ax, shrink=0.8)
    p = metrics["production"]
    fig.suptitle(f"{metrics['n']}³ SGS-off snapshot: Gaussian width {metrics['width_fine_cells']:g} cells, t={metrics['time_T']:.3f}T\n"
                 f"Native closure with exact K: correlation {p['signed_correlation']:.3f}; "
                 f"mean ratio {p['mean_ratio']:.3f}; RMS error {p['relative_rms_error']:.1%}")
    fig.savefig(out / "snapshot-comparison.png", dpi=160)
    fig.savefig(out / "snapshot-comparison.pdf")
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for index, n in enumerate(sorted(set(row["n"] for row in rows))):
        color = f"C{index}"
        for kind, style in [("gaussian", "o-"), ("box", "s--")]:
            selected = [r for r in rows if r["n"] == n and r["filter"] == kind]
            deltas = sorted(set(r["delta"] for r in selected))
            means = []
            for delta in deltas:
                group = [r for r in selected if r["delta"] == delta]
                exact_mean = np.mean([r["production"]["exact_mean"] for r in group])
                model_mean = np.mean([r["production"]["model_mean"] for r in group])
                means.append([model_mean / exact_mean,
                              np.mean([r["production"]["signed_correlation"] for r in group]),
                              np.mean([r["production"]["relative_rms_error"] for r in group]),
                              np.mean([r["K_reference_mean"] for r in group])])
            for index, ax in enumerate(axes.flat):
                ax.plot(np.asarray(deltas), np.asarray(means)[:, index], style,
                        color=color, label=f"{n}³ {kind}")
    for ax, title, ylabel in zip(axes.flat,
                                ["Mean transfer amplitude", "Signed correlation", "Local transfer error", "Direct unresolved energy"],
                                ["Model mean / measured mean", "Pearson correlation", "RMS error / measured RMS", "Mean K reference"]):
        ax.set(xlabel="Physical filter width / box length", title=title, ylabel=ylabel,
               xscale="log", xticks=[1/16, 1/8, 1/4], xticklabels=["1/16", "1/8", "1/4"])
        ax.minorticks_off()
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    axes[0, 0].axhline(1, color="black", ls=":")
    count = len(set(round(row["time_T"]) for row in rows))
    start, stop = min(row["time_T"] for row in rows), max(row["time_T"] for row in rows)
    fig.suptitle("Filter and reference-resolution sensitivity: fixed coefficients\n"
                 f"Equal-weight {count} nominal-time snapshots, {start:.1f}T-{stop:.1f}T; "
                 "limited scale separation; 32³ primary Gaussian is undersampled")
    fig.savefig(out / "filter-sensitivity.png", dpi=160)
    fig.savefig(out / "filter-sensitivity.pdf")
    plt.close(fig)


    selected = sorted([row for row in rows if row["n"] == metrics["n"]
                       and row["filter"] == "gaussian" and row["delta"] == 1/16],
                      key=lambda row: row["time_T"])
    times = [row["time_T"] for row in selected]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for key, label in [("exact_mean", "Measured"), ("model_mean", "Native closure")]:
        axes[0, 0].plot(times, [row["production"][key] for row in selected], "o-", label=label)
    for key, label in [("production", "Centered derivatives"),
                       ("spectral_derivative_sensitivity", "Spectral derivatives")]:
        axes[0, 1].plot(times, [row[key]["model_mean"] / row[key]["exact_mean"]
                               for row in selected], "o-", label=label)
    axes[0, 1].axhline(1, color="black", ls=":")
    for key, label in [("production", "Full production"),
                       ("deviatoric_production", "Trace-free production")]:
        axes[1, 0].plot(times, [row[key]["signed_correlation"] for row in selected],
                       "o-", label=label)
    for key, label in [("exact_backscatter_fraction", "Measured"),
                       ("model_backscatter_fraction", "Native closure")]:
        axes[1, 1].plot(times, [row["production"][key] for row in selected], "o-", label=label)
    for ax, title, ylabel in zip(axes.flat,
                                ["Mean transfer", "Amplitude and derivative sensitivity",
                                 "Correlation without the shared pressure term", "Backscatter volume fraction"],
                                ["Energy density / time", "Model mean / measured mean",
                                 "Pearson correlation", "Fraction with negative transfer"]):
        ax.set(xlabel="Time / T", title=title, ylabel=ylabel)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    fig.suptitle(f"{metrics['n']}³ Gaussian comparison across saved snapshots: Delta_G=box/16\n"
                 "Measured covariance K supplied; fixed coefficients; no forcing adjustment")
    fig.savefig(out / "time-diagnostics.png", dpi=160)
    fig.savefig(out / "time-diagnostics.pdf")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--snapshots", type=Path, nargs="+")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    kernel = args.kernel.resolve(strict=True)
    if args.self_test:
        print(json.dumps(self_test(kernel), indent=2))
        return
    if not args.snapshots or args.output is None:
        parser.error("--snapshots and --output are required for a snapshot comparison")
    args.output.mkdir(parents=True, exist_ok=False)
    protocol = dict(coefficients=dict(c1=0.02, c2=0.7),
                    filter_widths_box_units=[1/16, 1/8, 1/4],
                    gaussian_width="Delta_G=sqrt(12)*sigma; sampled positive periodic kernel",
                    primary="largest-grid latest snapshot, Gaussian Delta_G=1/16",
                    spatial_sampling="filtered fine grid, physical x/y/z; no downsampling",
                    supplied_K="exact filtered covariance, not evolved SGS reservoir",
                    gradient="centered fine-grid differences; primary spectral sensitivity",
                    kernel_path=str(kernel), kernel_sha256=verify.sha(kernel),
                    source_hashes={name: verify.sha(verify.SOURCE / name) for name in
                                   ["tst/sgs/filter_reference.py", "tst/sgs/verify_filtered_snapshot.py",
                                    "tst/sgs/filtered_stress.cpp", "src/hydro/sgs_closure.hpp"]},
                    arithmetic_thresholds=dict(scaled_default=2e-11, native_Python=5e-13),
                    physical_acceptance="Report bias and scatter; no tuned acceptance threshold",
                    validated_DNS=False, paper_reference_reproduced=False)
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2)+'\n')
    checks = self_test(kernel)
    rows, fields = [], None
    for path in args.snapshots:
        # ponytail: reuse the cell-dictionary reader; large grids need streamed input.
        state = verify.restart(path)
        n = int(state["config"]["mesh/nx1"])
        if (state["config"].get("sgs/enabled", "false") not in ["false", "0"]
                or any(int(state["config"][f"mesh/nx{axis}"]) != n for axis in [2, 3])
                or not state["finite"] or state["min_density"] <= 0 or state["min_thermal"] <= 0):
            raise ValueError("require finite positive-density cubic SGS-off reference")
        if any(state["config"].get(f"coord/{name}", "false") in ["true", "1"]
               for name in ["special_rel", "general_rel"]):
            raise ValueError("Newtonian reference required")
        for axis in [1, 2, 3]:
            for side in ["i", "o"]:
                if state["config"][f"mesh/{side}x{axis}_bc"] != "periodic":
                    raise ValueError("periodic boundaries required")
            if (float(state["config"][f"mesh/x{axis}max"])
                    - float(state["config"][f"mesh/x{axis}min"])) != 1:
                raise ValueError("unit periodic cube required")
        data = np.array([state["field"][key] for key in sorted(state["field"])]).reshape(n, n, n, -1)
        rho, velocity = data[..., 0], data[..., 1:4] / data[..., :1]
        for kind in ["gaussian", "box"]:
            for delta in protocol["filter_widths_box_units"]:
                detailed = delta == 1/16 and kind == "gaussian"
                metrics, ref, tau, production = assess(rho, velocity, delta, kind, kernel, detailed)
                metrics.update(path=str(path.resolve()), sha256=verify.sha(path), time=state["time"],
                               time_T=state["time"] / 0.5,
                               initial_density=state["config"]["problem/density0"],
                               gamma=state["config"]["hydro/gamma"])
                rows.append(metrics)
                if detailed and (fields is None or (n, state["time"]) >=
                                 (fields[-1]["n"], fields[-1]["time"])):
                    fields = ref, tau, production, metrics
                print(n, kind, delta, f"t={metrics['time_T']:.3f}T",
                      {k: metrics["production"][k] for k in ["signed_correlation", "mean_ratio", "relative_rms_error"]}, flush=True)
        del state, data
        (args.output / "results.json").write_text(json.dumps(
            dict(status="running", self_checks=checks, protocol=protocol, rows=rows),
            indent=2, allow_nan=False)+'\n')
    if fields is not None:
        plots(args.output, rows, fields)
    (args.output / "results.json").write_text(json.dumps(
        dict(status="complete", self_checks=checks, protocol=protocol, rows=rows),
        indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
