"""Verification tests with analytic solutions for accuracy benchmarking.

This version keeps the logic minimal: analytic scenarios are plain dictionaries
with functional callbacks for solutions, RHS, and Jacobians. No inheritance or
heavy abstractions remain.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False
    plt = None

try:
    from scipy.integrate import solve_ivp

    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False
    solve_ivp = None

try:
    from scikits.odes import ode as _scikits_ode

    HAS_SUNDIALS = True
except Exception:
    HAS_SUNDIALS = False

_THIS_DIR = Path(__file__).resolve().parent
_ROOT = _THIS_DIR.parent  # chemistry/
_SRC_ROOT = _ROOT.parent  # src/
for p in (str(_ROOT), str(_SRC_ROOT)):
    if p not in sys.path:
        sys.path.insert(0, p)

from chemistry.custom_network_loader import load_network
from chemistry.tests.test_network_framework import run_framework_smoke


# --------------------------------------------------------------------------- #
# Scenario definitions
# --------------------------------------------------------------------------- #


def _exp(val):
    return math.exp(val)


def _scenario_base(name, description, species, params, tol_abs, tol_rel, chem, rate):
    return {
        "name": name,
        "description": description,
        "species": species,
        "params": params,
        "tolerances": {"abs": tol_abs, "rel": tol_rel},
        "equations": {"chem": chem, "rate": rate},
    }


def make_scenarios(overrides: Dict[str, Dict[str, float]]) -> List[Dict]:
    scenarios: List[Dict] = []

    scenarios.append(
        _scenario_base(
            "Exponential Decay",
            "A → B",
            ["A", "B"],
            {"k": overrides.get("ExponentialDecay", {}).get("k", 1.0)},
            1e-7,
            1e-6,
            "A → B",
            "d[A]/dt=-k[A]; d[B]/dt=k[A]",
        )
    )

    scenarios.append(
        _scenario_base(
            "Equilibrium Approach",
            "A ⇌ B",
            ["A", "B"],
            {
                "k_f": overrides.get("EquilibriumApproach", {}).get("k_f", 1.0),
                "k_b": overrides.get("EquilibriumApproach", {}).get("k_b", 0.5),
            },
            1e-7,
            1e-6,
            "A ⇌ B",
            "d[A]/dt=-k_f[A]+k_b[B]; d[B]/dt=k_f[A]-k_b[B]",
        )
    )

    scenarios.append(
        _scenario_base(
            "Near-Equilibrium Stiff",
            "A ⇌ B with k_f ≈ k_b ≫ 1",
            ["A", "B"],
            {
                "k_f": overrides.get("NearEquilibriumStiff", {}).get("k_f", 2500.0),
                "k_b": overrides.get("NearEquilibriumStiff", {}).get("k_b", 2499.0),
            },
            1e-6,
            5e-5,
            "A ⇌ B (fast, near balance)",
            "d[A]/dt=-k_f[A]+k_b[B]; d[B]/dt=k_f[A]-k_b[B]",
        )
    )

    scenarios.append(
        _scenario_base(
            "Sequential Decay",
            "A → B → C",
            ["A", "B", "C"],
            {
                "k1": overrides.get("SequentialDecay", {}).get("k1", 1.0),
                "k2": overrides.get("SequentialDecay", {}).get("k2", 2.0),
            },
            1e-7,
            1e-6,
            "A → B → C",
            "d[A]/dt=-k1[A]; d[B]/dt=k1[A]-k2[B]; d[C]/dt=k2[B]",
        )
    )

    scenarios.append(
        _scenario_base(
            "Stiff System",
            "A →(fast) B →(slow) C",
            ["A", "B", "C"],
            {
                "k_fast": overrides.get("StiffSystem", {}).get("k_fast", 1000.0),
                "k_slow": overrides.get("StiffSystem", {}).get("k_slow", 0.1),
            },
            5e-5,
            1e-4,
            "A →(k_fast) B →(k_slow) C",
            "d[A]/dt=-k_fast[A]; d[B]/dt=k_fast[A]-k_slow[B]; d[C]/dt=k_slow[B]",
        )
    )

    scenarios.append(
        _scenario_base(
            "Quadratic Sink",
            "2A → B (non-linear)",
            ["A", "B"],
            {"k": overrides.get("QuadraticSink", {}).get("k", 0.8)},
            1e-6,
            1e-5,
            "2A → B",
            "d[A]/dt=-2k[A]^2; d[B]/dt=k[A]^2",
        )
    )

    return scenarios


# --------------------------------------------------------------------------- #
# Mock KIDA loader (for reproducible params)
# --------------------------------------------------------------------------- #


def _format_reaction_line(
    r1: str,
    r2: str,
    r3: str,
    products: List[str],
    alpha: float,
    beta: float = 0.0,
    gamma: float = 0.0,
    f: float = 0.0,
    g: float = 0.0,
    itype: int = 0,
    tmin: int = -9999,
    tmax: int = 9999,
) -> str:
    products = (products + [""] * 5)[:5]
    return (
        f"{r1:<11}{r2:<11}{r3:<11} "
        f"{products[0]:<11}{products[1]:<11}{products[2]:<11}{products[3]:<11}{products[4]:<11}"
        f"{alpha:>10.3e} {beta:>10.3e} {gamma:>10.3e} {f:>8.2f} {g:>8.2f} "
        f"{itype:2d} {tmin:6d} {tmax:6d}  0   0  0\n"
    )


def ensure_mock_kida_network(mock_dir: Path) -> None:
    mock_dir.mkdir(parents=True, exist_ok=True)
    spec_readme = mock_dir / "kida_spec_readme_mock.txt"
    if not spec_readme.exists():
        spec_readme.write_text(
            "Mock KIDA network\n"
            "Number of species: 3\n"
            "Number of reactions: 6\n"
            "Elements are the following: H, C, O\n"
        )

    spec_file = mock_dir / "kida_spec_mock.dat"
    if not spec_file.exists():
        spec_file.write_text(
            "! name charge H C O\n"
            "A 0 1 0 0\n"
            "B 0 0 1 0\n"
            "C 0 0 0 1\n"
        )

    reac_file = mock_dir / "kida_reac_mock.dat"
    if not reac_file.exists():
        lines = []
        lines.append(_format_reaction_line("A", "", "", ["B"], 1.0))
        lines.append(_format_reaction_line("B", "", "", ["A"], 0.5))
        lines.append(_format_reaction_line("A", "", "", ["B"], 1000.0))
        lines.append(_format_reaction_line("B", "", "", ["C"], 0.1))
        lines.append(_format_reaction_line("B", "", "", ["C"], 2.0))
        lines.append(_format_reaction_line("A", "A", "", ["B"], 0.8))
        reac_file.write_text("".join(lines))


def load_mock_rates_from_kida(mock_dir: Path) -> Dict[str, Dict[str, float]]:
    ensure_mock_kida_network(mock_dir)
    net = load_network(mock_dir)
    if not net.get("species") or not net.get("reactions"):
        return {}

    def match(reactants, products, alpha):
        for rxn in net.get("reactions", []):
            if rxn["reactants"] == reactants and rxn["products"] == products:
                if abs(rxn["alpha"] - alpha) < 1e-12:
                    return rxn["alpha"]
        return None

    return {
        "ExponentialDecay": {"k": match(["A"], ["B"], 1.0) or 1.0},
        "EquilibriumApproach": {
            "k_f": match(["A"], ["B"], 1.0) or 1.0,
            "k_b": match(["B"], ["A"], 0.5) or 0.5,
        },
        "NearEquilibriumStiff": {"k_f": 2500.0, "k_b": 2499.0},
        "SequentialDecay": {
            "k1": match(["A"], ["B"], 1.0) or 1.0,
            "k2": match(["B"], ["C"], 2.0) or 2.0,
        },
        "StiffSystem": {
            "k_fast": match(["A"], ["B"], 1000.0) or 1000.0,
            "k_slow": match(["B"], ["C"], 0.1) or 0.1,
        },
        "QuadraticSink": {"k": match(["A", "A"], ["B"], 0.8) or 0.8},
    }


# --------------------------------------------------------------------------- #
# Scenario mechanics
# --------------------------------------------------------------------------- #


def scenario_initial_conditions(scenario: Dict) -> Dict[str, float]:
    species = scenario["species"]
    if scenario["name"] == "Exponential Decay":
        return {"A": 1.0, "B": 0.0}
    if scenario["name"] == "Equilibrium Approach":
        return {"A": 1.0, "B": 0.0}
    if scenario["name"] == "Near-Equilibrium Stiff":
        return {"A": 1.0, "B": 0.0}
    if scenario["name"] in ("Sequential Decay", "Stiff System"):
        return {"A": 1.0, "B": 0.0, "C": 0.0}
    if scenario["name"] == "Quadratic Sink":
        return {"A": 1.0, "B": 0.0}
    return {s: 0.0 for s in species}


def scenario_solution(scenario: Dict, t: float) -> Dict[str, float]:
    p = scenario["params"]
    name = scenario["name"]
    exp = _exp
    if name == "Exponential Decay":
        k = p["k"]
        A0, B0 = 1.0, 0.0
        A_t = A0 * exp(-k * t)
        B_t = A0 * (1 - exp(-k * t)) + B0
        return {"A": A_t, "B": B_t}
    if name == "Equilibrium Approach":
        kf, kb = p["k_f"], p["k_b"]
        A0, B0 = 1.0, 0.0
        total = A0 + B0
        A_eq = total * kb / (kf + kb)
        B_eq = total * kf / (kf + kb)
        decay = exp(-(kf + kb) * t)
        A_t = A_eq + (A0 - A_eq) * decay
        B_t = B_eq + (B0 - B_eq) * decay
        return {"A": A_t, "B": B_t}
    if name == "Near-Equilibrium Stiff":
        kf, kb = p["k_f"], p["k_b"]
        A0, B0 = 1.0, 0.0
        total = A0 + B0
        A_eq = total * kb / (kf + kb)
        B_eq = total * kf / (kf + kb)
        decay = exp(-(kf + kb) * t)
        return {"A": A_eq + (A0 - A_eq) * decay, "B": B_eq + (B0 - B_eq) * decay}
    if name == "Sequential Decay":
        k1, k2 = p["k1"], p["k2"]
        A0, B0, C0 = 1.0, 0.0, 0.0
        e1 = exp(-k1 * t)
        e2 = exp(-k2 * t)
        A_t = A0 * e1
        if abs(k2 - k1) > 1e-10:
            B_t = A0 * k1 / (k2 - k1) * (e1 - e2) + B0 * e2
            C_t = (
                A0 * (1 + (k1 * e2 - k2 * e1) / (k2 - k1))
                + B0 * (1 - e2)
                + C0
            )
        else:
            B_t = A0 * k1 * t * e1 + B0 * e2
            C_t = A0 * (1 - (1 + k1 * t) * e1) + B0 * (1 - e2) + C0
        return {"A": A_t, "B": B_t, "C": C_t}
    if name == "Stiff System":
        kf, ks = p["k_fast"], p["k_slow"]
        A0, B0, C0 = 1.0, 0.0, 0.0
        ef = exp(-kf * t)
        es = exp(-ks * t)
        A_t = A0 * ef
        if abs(ks - kf) > 1e-10:
            B_t = A0 * kf / (ks - kf) * (ef - es) + B0 * es
            C_t = (
                A0 * (1 + (kf * es - ks * ef) / (ks - kf))
                + B0 * (1 - es)
                + C0
            )
        else:
            B_t = A0 * kf * t * ef + B0 * es
            C_t = A0 * (1 - (1 + kf * t) * ef) + B0 * (1 - es) + C0
        return {"A": A_t, "B": B_t, "C": C_t}
    if name == "Quadratic Sink":
        k = p["k"]
        A0, B0 = 1.0, 0.0
        denom = 1 + 2 * k * A0 * t
        A_t = A0 / denom
        B_t = B0 + 0.5 * A0 * (1 - 1 / denom)
        return {"A": A_t, "B": B_t}
    return {}


def scenario_rhs(scenario: Dict, concs: List[float]) -> List[float]:
    p = scenario["params"]
    sp = scenario["species"]
    idx = {s: i for i, s in enumerate(sp)}
    y = concs
    rhs = [0.0] * len(sp)
    n = scenario["name"]
    if n == "Exponential Decay":
        rhs[idx["A"]] = -p["k"] * y[idx["A"]]
        rhs[idx["B"]] = p["k"] * y[idx["A"]]
    elif n in ("Equilibrium Approach", "Near-Equilibrium Stiff"):
        rhs[idx["A"]] = -p["k_f"] * y[idx["A"]] + p["k_b"] * y[idx["B"]]
        rhs[idx["B"]] = p["k_f"] * y[idx["A"]] - p["k_b"] * y[idx["B"]]
    elif n == "Sequential Decay":
        rhs[idx["A"]] = -p["k1"] * y[idx["A"]]
        rhs[idx["B"]] = p["k1"] * y[idx["A"]] - p["k2"] * y[idx["B"]]
        rhs[idx["C"]] = p["k2"] * y[idx["B"]]
    elif n == "Stiff System":
        rhs[idx["A"]] = -p["k_fast"] * y[idx["A"]]
        rhs[idx["B"]] = p["k_fast"] * y[idx["A"]] - p["k_slow"] * y[idx["B"]]
        rhs[idx["C"]] = p["k_slow"] * y[idx["B"]]
    elif n == "Quadratic Sink":
        rate = p["k"] * y[idx["A"]] * y[idx["A"]]
        rhs[idx["A"]] = -2 * rate
        rhs[idx["B"]] = rate
    return rhs


def scenario_jacobian(scenario: Dict, concs: List[float]) -> np.ndarray:
    p = scenario["params"]
    sp = scenario["species"]
    idx = {s: i for i, s in enumerate(sp)}
    y = concs
    n = len(sp)
    J = np.zeros((n, n), dtype=float)
    name = scenario["name"]
    if name == "Exponential Decay":
        k = p["k"]
        J[idx["A"], idx["A"]] = -k
        J[idx["B"], idx["A"]] = k
    elif name in ("Equilibrium Approach", "Near-Equilibrium Stiff"):
        kf, kb = p["k_f"], p["k_b"]
        J[idx["A"], idx["A"]] = -kf
        J[idx["A"], idx["B"]] = kb
        J[idx["B"], idx["A"]] = kf
        J[idx["B"], idx["B"]] = -kb
    elif name == "Sequential Decay":
        k1, k2 = p["k1"], p["k2"]
        J[idx["A"], idx["A"]] = -k1
        J[idx["B"], idx["A"]] = k1
        J[idx["B"], idx["B"]] = -k2
        J[idx["C"], idx["B"]] = k2
    elif name == "Stiff System":
        kf, ks = p["k_fast"], p["k_slow"]
        J[idx["A"], idx["A"]] = -kf
        J[idx["B"], idx["A"]] = kf
        J[idx["B"], idx["B"]] = -ks
        J[idx["C"], idx["B"]] = ks
    elif name == "Quadratic Sink":
        k = p["k"]
        d_rate_dA = 2 * k * y[idx["A"]]
        J[idx["A"], idx["A"]] = -2 * d_rate_dA
        J[idx["B"], idx["A"]] = d_rate_dA
    return J


# --------------------------------------------------------------------------- #
# Metrics and test runners
# --------------------------------------------------------------------------- #


def compute_error_metrics(
    times: List[float],
    numerical: List[Dict[str, float]],
    analytic: List[Dict[str, float]],
    species: List[str],
) -> Dict:
    errors = {s: [] for s in species}
    for i in range(len(times)):
        for s in species:
            num_val = numerical[i].get(s, 0.0)
            ana_val = analytic[i].get(s, 0.0)
            abs_err = abs(num_val - ana_val)
            rel_err = abs_err / abs(ana_val) if abs(ana_val) > 1e-30 else abs_err
            errors[s].append({"abs": abs_err, "rel": rel_err})

    stats = {}
    for s in species:
        abs_errs = [e["abs"] for e in errors[s]]
        rel_errs = [e["rel"] for e in errors[s]]
        stats[s] = {
            "abs_max": max(abs_errs) if abs_errs else 0.0,
            "abs_mean": sum(abs_errs) / len(abs_errs) if abs_errs else 0.0,
            "rel_max": max(rel_errs) if rel_errs else 0.0,
            "rel_mean": sum(rel_errs) / len(rel_errs) if rel_errs else 0.0,
            "abs_errors": abs_errs,
            "rel_errors": rel_errs,
        }
    return stats


def validate_errors(error_stats: Dict, abs_tol: float, rel_tol: float) -> Dict[str, Dict]:
    failures = {}
    for species, stats in error_stats.items():
        if stats["abs_max"] > abs_tol or stats["rel_max"] > rel_tol:
            failures[species] = {"abs_max": stats["abs_max"], "rel_max": stats["rel_max"]}
    return failures


def run_scipy(
    scenario: Dict, t_start=0.0, t_end=5.0, n_points=50, rtol=1e-6, atol=1e-9, verbose=True
) -> Tuple[bool, Dict, Dict]:
    if not HAS_SCIPY:
        return False, {"status": "skipped", "reason": "SciPy not available"}, {}

    species = scenario["species"]
    x0_dict = scenario_initial_conditions(scenario)
    x0 = [x0_dict.get(s, 0.0) for s in species]

    def rhs(t, y):
        return scenario_rhs(scenario, y)

    sol = solve_ivp(
        rhs,
        (t_start, t_end),
        x0,
        method="BDF",
        rtol=rtol,
        atol=atol,
        t_eval=np.linspace(t_start, t_end, n_points),
        dense_output=False,
    )
    if not sol.success:
        return False, {"status": "failed", "error": f"Integration status {sol.status}"}, {}

    times = list(sol.t)
    numerical = []
    for i, t in enumerate(times):
        y = sol.y[:, i]
        numerical.append({species[j]: float(y[j]) for j in range(len(species))})

    analytic_vals = [scenario_solution(scenario, t) for t in times]
    error_stats = compute_error_metrics(times, numerical, analytic_vals, species)
    thresholds = scenario["tolerances"]
    failures = validate_errors(error_stats, thresholds["abs"], thresholds["rel"])
    result = {
        "status": "passed" if not failures else "failed",
        "solver": "SciPy BDF",
        "rtol": rtol,
        "atol": atol,
        "time_range": (t_start, t_end),
        "n_points": len(times),
        "error_stats": error_stats,
        "error_thresholds": thresholds,
        "failures": failures,
    }
    data = {"times": times, "numerical": numerical, "analytic": analytic_vals}
    return len(failures) == 0, result, data


def run_sundials(
    scenario: Dict, t_start=0.0, t_end=5.0, n_points=50, rtol=1e-6, atol=1e-9, verbose=True
) -> Tuple[bool, Dict, Dict]:
    if not HAS_SUNDIALS:
        return True, {"status": "skipped", "reason": "SUNDIALS not available"}, {}

    species = scenario["species"]
    x0_dict = scenario_initial_conditions(scenario)
    x0 = np.asarray([x0_dict.get(s, 0.0) for s in species], dtype=float)

    def rhs(t, y, ydot):
        ydot[:] = np.asarray(scenario_rhs(scenario, y), dtype=float)
        return 0

    def jac(t, y, fy, J):
        J[:, :] = scenario_jacobian(scenario, y)
        return 0

    solver = _scikits_ode("cvode", rhs, jacfn=jac, rtol=rtol, atol=atol)
    t_eval = np.linspace(t_start, t_end, n_points)
    sol = solver.solve(t_eval, x0)
    try:
        y_vals = sol.y
    except Exception:
        try:
            y_vals = sol.values.y
        except Exception as exc:
            return False, {"status": "failed", "error": str(exc)}, {}

    times = list(t_eval)
    numerical = []
    for i in range(len(times)):
        y = y_vals[i]
        numerical.append({species[j]: float(y[j]) for j in range(len(species))})

    analytic_vals = [scenario_solution(scenario, t) for t in times]
    error_stats = compute_error_metrics(times, numerical, analytic_vals, species)
    thresholds = scenario["tolerances"]
    failures = validate_errors(error_stats, thresholds["abs"], thresholds["rel"])
    result = {
        "status": "passed" if not failures else "failed",
        "solver": "SUNDIALS CVODE",
        "rtol": rtol,
        "atol": atol,
        "time_range": (t_start, t_end),
        "n_points": len(times),
        "error_stats": error_stats,
        "error_thresholds": thresholds,
        "failures": failures,
    }
    data = {"times": times, "numerical": numerical, "analytic": analytic_vals}
    return len(failures) == 0, result, data


# --------------------------------------------------------------------------- #
# Plotting helpers
# --------------------------------------------------------------------------- #


def plot_comparison(scenario, solutions_by_solver, output_dir):
    """Overlay SciPy and SUNDIALS on one figure per scenario."""
    if not HAS_MATPLOTLIB or not solutions_by_solver:
        return
    species = scenario["species"]
    fig, axes = plt.subplots(len(species), 1, figsize=(10, 3 * len(species)), sharex=False)
    if len(species) == 1:
        axes = [axes]
    chem_eq = scenario["equations"].get("chem", "")
    rate_eq = scenario["equations"].get("rate", "")

    colors = {"scipy": "r", "sundials": "g"}
    linestyles = {"scipy": "--", "sundials": "-."}

    for i, sp in enumerate(species):
        ax = axes[i]
        # Analytic curve at union of all times
        all_times = sorted(
            set(
                t
                for data in solutions_by_solver.values()
                for t in data["times"]
            )
        )
        ana_vals = [scenario_solution(scenario, t).get(sp, 0.0) for t in all_times]
        ax.plot(
            all_times,
            ana_vals,
            "b-",
            linewidth=2.0,
            label="Analytic",
            alpha=0.8,
        )
        # Each solver
        for solver_name, data in solutions_by_solver.items():
            times = data["times"]
            num_vals = [n.get(sp, 0.0) for n in data["numerical"]]
            ax.plot(
                times,
                num_vals,
                color=colors.get(solver_name, "k"),
                linestyle=linestyles.get(solver_name, ":"),
                linewidth=1.5,
                label=solver_name.upper(),
                alpha=0.8,
                marker="o",
                markersize=3,
            )
        ax.set_ylabel(f"[{sp}] (mol/cm³)", fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.legend(fontsize=10, loc="best")
        ax.set_title(
            f"{sp} — {scenario['name']}" + (f" | {chem_eq}" if chem_eq else ""),
            fontsize=12,
            fontweight="bold",
        )
        if chem_eq or rate_eq:
            ax.text(
                0.02,
                0.70,
                f"Chem: {chem_eq}\nRate: {rate_eq}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
            )
    axes[-1].set_xlabel("Time (s)", fontsize=11, fontweight="bold")
    plt.tight_layout()
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    filename = outdir / f"{scenario['name'].replace(' ', '_')}_comparison.png"
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    📊 Plot saved: {filename}")


def plot_combined(scenario, solutions_by_solver, output_dir):
    if not HAS_MATPLOTLIB or not solutions_by_solver:
        return
    species = scenario["species"]
    solvers = list(solutions_by_solver.keys())
    # Combined plot removed (redundant); keep helper for compatibility but no-op.


def plot_error_metrics(scenario, results, solution_data, output_dir):
    if not HAS_MATPLOTLIB:
        return
    species = scenario["species"]
    solvers = [
        s
        for s, r in results.items()
        if r and r.get("error_stats") and s in solution_data
    ]
    if not solvers:
        return
    fig, axes = plt.subplots(len(species), 2, figsize=(14, 4 * len(species)))
    if len(species) == 1:
        axes = axes.reshape(1, 2)
    for sp_idx, sp in enumerate(species):
        ax_abs = axes[sp_idx, 0]
        ax_rel = axes[sp_idx, 1]
        for solver in solvers:
            err_stats = results[solver]["error_stats"][sp]
            times = solution_data[solver]["times"]
            ax_abs.semilogy(
                times,
                err_stats["abs_errors"],
                marker="o",
                label=solver.upper(),
                linewidth=2,
                alpha=0.7,
            )
            ax_rel.semilogy(
                times,
                err_stats["rel_errors"],
                marker="s",
                label=solver.upper(),
                linewidth=2,
                alpha=0.7,
            )
        ax_abs.set_ylabel("Absolute Error", fontsize=10, fontweight="bold")
        ax_rel.set_ylabel("Relative Error", fontsize=10, fontweight="bold")
        ax_abs.set_title(f"{sp} - Absolute Error", fontsize=11, fontweight="bold")
        ax_rel.set_title(f"{sp} - Relative Error", fontsize=11, fontweight="bold")
        ax_abs.grid(True, alpha=0.3, which="both", linestyle="--")
        ax_rel.grid(True, alpha=0.3, which="both", linestyle="--")
        ax_abs.legend(fontsize=9)
        ax_rel.legend(fontsize=9)
    axes[-1, 0].set_xlabel("Time (s)", fontsize=10, fontweight="bold")
    axes[-1, 1].set_xlabel("Time (s)", fontsize=10, fontweight="bold")
    plt.tight_layout()
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)
    filename = outdir / f"{scenario['name'].replace(' ', '_')}_error_metrics.png"
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    📊 Error plot saved: {filename}")


# --------------------------------------------------------------------------- #
# Results markdown
# --------------------------------------------------------------------------- #


def write_results_markdown(test_runs, output_path: Path, framework_result=None):
    def status_cell(result: Optional[Dict]) -> str:
        if not result:
            return "❌ Failed"
        status = result.get("status", "failed")
        if status == "passed":
            return "✅ Passed"
        if status == "skipped":
            reason = result.get("reason", "skipped")
            return f"⏭️ Skipped ({reason})"
        reason = result.get("error")
        failures = result.get("failures") or {}
        if failures and not reason:
            reason = f"tolerance breach: {', '.join(sorted(failures.keys()))}"
        return f"❌ Failed ({reason})" if reason else "❌ Failed"

    def max_errors(result: Optional[Dict]) -> Tuple[str, str]:
        if not result or "error_stats" not in result:
            return ("n/a", "n/a")
        error_stats = result["error_stats"]
        max_abs = max((stats["abs_max"] for stats in error_stats.values()), default=0.0)
        max_rel = max((stats["rel_max"] for stats in error_stats.values()), default=0.0)
        return (f"{max_abs:.3e}", f"{max_rel:.3e}")

    lines = []
    lines.append("# Analytic Solution Test Results")
    lines.append("")
    lines.append("_Auto-generated by test_analytic_solutions.py; re-run the test to refresh._")
    lines.append("")
    lines.append("## Status Legend")
    lines.append("- ✅ Passed — within abs/rel tolerances")
    lines.append("- ❌ Failed — tolerance breach or solver error")
    lines.append("- ⏭️ Skipped — solver unavailable")
    lines.append("")
    if framework_result:
        lines.append("## Framework Smoke Test")
        fr_status = framework_result.get("status", "failed")
        fr_msg = framework_result.get("message", "")
        fr_line = "✅ Passed" if fr_status == "passed" else f"❌ Failed ({fr_msg})"
        lines.append(f"- Network load + ODE build: {fr_line}")
        lines.append("")

    lines.append("## Summary")

    for scenario, results in test_runs:
        scipy_res = results.get("scipy")
        sundials_res = results.get("sundials")
        thresholds = scenario["tolerances"]
        notes = []
        if scipy_res and scipy_res.get("status") == "failed":
            fails = scipy_res.get("failures")
            if fails:
                notes.append(f"SciPy tolerance breach: {', '.join(sorted(fails.keys()))}")
            elif scipy_res.get("error"):
                notes.append(f"SciPy error: {scipy_res['error']}")
        if sundials_res and sundials_res.get("status") == "failed":
            fails = sundials_res.get("failures")
            if fails:
                notes.append(f"SUNDIALS tolerance breach: {', '.join(sorted(fails.keys()))}")
            elif sundials_res.get("error"):
                notes.append(f"SUNDIALS error: {sundials_res['error']}")
        if sundials_res and sundials_res.get("status") == "skipped":
            reason = sundials_res.get("reason", "skipped")
            notes.append(f"SUNDIALS skipped: {reason}")
        scipy_abs, scipy_rel = max_errors(scipy_res)
        sundials_abs, sundials_rel = max_errors(sundials_res)

        lines.append(f"- **{scenario['name']}**")
        lines.append(
            f"  - Tolerances: abs <= {thresholds['abs']:.1e}, rel <= {thresholds['rel']:.1e}"
        )
        lines.append(f"  - SciPy: {status_cell(scipy_res)}")
        lines.append(f"  - SciPy max error: abs={scipy_abs}, rel={scipy_rel}")
        lines.append(f"  - SUNDIALS: {status_cell(sundials_res)}")
        lines.append(f"  - SUNDIALS max error: abs={sundials_abs}, rel={sundials_rel}")
        if notes:
            lines.append(f"  - Notes: {'; '.join(notes)}")
        lines.append("")

    lines.append("## Per-Test Tolerances")
    lines.append("- Exponential, reversible, cascade: `abs <= 1e-7`, `rel <= 1e-6`")
    lines.append("- Near-equilibrium reversible: `abs <= 1e-6`, `rel <= 5e-5`")
    lines.append("- Stiff cascade: `abs <= 5e-5`, `rel <= 1e-4`")
    lines.append("- Quadratic sink: `abs <= 1e-6`, `rel <= 1e-5`")
    lines.append("")
    lines.append("## How to Recreate")
    lines.append("```bash")
    lines.append("python src/chemistry/test_analytic_solutions.py")
    lines.append("```")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines))
    print(f"\n📝 Wrote markdown summary to {output_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #


def run_scenario_tests(scenario: Dict, plotter_output="src/chemistry/test_outputs"):
    results = {}
    solution_data = {}
    scipy_ok, scipy_res, scipy_data = run_scipy(scenario)
    results["scipy"] = scipy_res
    if scipy_data:
        solution_data["scipy"] = scipy_data
    sundials_ok, sundials_res, sundials_data = run_sundials(scenario)
    results["sundials"] = sundials_res
    if sundials_data:
        solution_data["sundials"] = sundials_data

    # Plots
    if solution_data:
        plot_comparison(scenario, solution_data, plotter_output)
        plot_error_metrics(scenario, results, solution_data, plotter_output)

    return scipy_ok and sundials_ok, results, solution_data


def main():
    print("\n" + "=" * 70)
    print("ANALYTIC SOLUTION ACCURACY VERIFICATION TESTS")
    print("=" * 70)

    # Smoke test
    framework_ok, framework_msg = run_framework_smoke()
    print(f"{'✅' if framework_ok else '❌'} Framework smoke test: {framework_msg}")

    # Load mock KIDA rates and build scenarios
    mock_dir = _ROOT / "custom" / "mock_analytic_kida"
    mock_rates = load_mock_rates_from_kida(mock_dir)
    if mock_rates:
        print("ℹ️  Applying mock analytic rate overrides from mock KIDA network")
    scenarios = make_scenarios(mock_rates)

    all_passed = framework_ok
    test_runs = []
    plot_dir = "src/chemistry/test_outputs"

    for scenario in scenarios:
        print(f"\n{'='*70}")
        print(f"ACCURACY VERIFICATION TEST: {scenario['name']}")
        print(f"{'='*70}\n")
        print(f"System: {scenario['description']}")
        print(f"Initial conditions: {scenario_initial_conditions(scenario)}")
        ok, results, _ = run_scenario_tests(scenario, plotter_output=plot_dir)
        test_runs.append((scenario, results))
        scenario_passed = True
        for r in results.values():
            if r.get("status") not in ("passed", "skipped"):
                scenario_passed = False
        if not scenario_passed:
            all_passed = False
        print(f"\n{'='*70}")
        print(f"SUMMARY: {scenario['name']}")
        print(f"{'='*70}")
        for name, res in results.items():
            status = res.get("status")
            label = "✅ PASSED" if status == "passed" else ("⏭️ SKIPPED" if status == "skipped" else "❌ FAILED")
            print(f"  {name:<20} {label}")

    results_md = _ROOT / "tests" / "ANALYTIC_SOLUTIONS_RESULTS.md"
    framework_result = {"status": "passed" if framework_ok else "failed", "message": framework_msg}
    write_results_markdown(test_runs, results_md, framework_result=framework_result)

    print("\n" + "=" * 70)
    if all_passed:
        print("✅ All accuracy tests passed!")
        if HAS_MATPLOTLIB:
            print("📊 Plots saved to 'outputs/' directory")
        sys.exit(0)
    else:
        print("❌ Some accuracy tests failed")
        sys.exit(1)


if __name__ == "__main__":
    main()
