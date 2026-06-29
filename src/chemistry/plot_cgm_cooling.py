#!/usr/bin/env python3
"""
Plot CGM cooling and heating rates using the same tables and formulas
implemented in SourceTerms::CGMCooling (src/srcterms/srcterms.cpp).

Reads the C++ cooling tables directly from src/srcterms/cooling_tables.hpp.
All physical and numerical parameters are set in-code below (no athinput
needed). The script reproduces the interpolation logic,
including CIE/PIE blending, low-temperature fit, shielding factor, and the
high-temperature heating cutoff. Outputs a PNG and shows the plot.
"""

import math
import pathlib
import re
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np


def find_tables_path() -> pathlib.Path:
    """Locate the cooling tables header from this file or the current directory."""
    candidates = []

    if "__file__" in globals():
        this_file = pathlib.Path(__file__).resolve()
        candidates.extend(this_file.parents)

    candidates.extend(pathlib.Path.cwd().resolve().parents)
    candidates.append(pathlib.Path.cwd().resolve())

    seen = set()
    for base in candidates:
        if base in seen:
            continue
        seen.add(base)

        direct = base / "src" / "srcterms" / "cooling_tables.hpp"
        if direct.is_file():
            return direct

        nested = base / "srcterms" / "cooling_tables.hpp"
        if nested.is_file():
            return nested

    raise FileNotFoundError(
        "Could not locate src/srcterms/cooling_tables.hpp from the script or cwd"
    )


TABLES_PATH = find_tables_path()
AMU = 1.67262192369e-24  # g
KB = 1.3806488e-16  # erg/K
X_H = 0.75
Z_SOL = 0.02
T_HOT = 2000  # 1e6

# Temperature grid and environment used for the standalone plot
TMIN = 1e1  # K
TMAX = 1e8  # K
NSAMPLES = 400
Z_REL = 0.2  # Z/Z_sun

# --- User-adjustable parameters (no athinput needed) ---
# Tune these to match a specific run or explore parameter sensitivity.
USER_CFG: Dict[str, float] = {
    # Units (matching turb_dust.athinput defaults)
    "length_cgs": 3.0856775809623245e21,  # cm
    "mass_cgs": 3.036951775493658e40,  # g
    "time_cgs": 3.15576e13,  # s
    "mu": 0.6,
    "rho0": 5.0,
    # Grid extents (for shielding length estimate)
    "x1min": -0.5,
    "x1max": 0.5,
    "nx1": 256,
    # Heating options
    "hrate": 2e-26,  # erg s^-1 per H atom
    "hscale_norm": 1.0,  # N_H * T_HOT / 1e4,
    "hscale_flag": 0.0,
    "hscale_height": 0.0,
    "hscale_radius": 0.0,
    "hscale_alpha": 0.0,
    # Cooling ceiling
    "T_max": 5e8,
    "T_cutoff": 5e5,
}


def load_tables(path: pathlib.Path) -> Dict[str, np.ndarray]:
    """Read the tabulated cooling/heating arrays from the C++ header."""
    text = path.read_text()

    def grab_array(name: str) -> np.ndarray:
        pattern = rf"constexpr\s+Real\s+{name}\s*\[[^\]]+\]\s*=\s*\{{([^}}]+)\}};"
        m = re.search(pattern, text, re.S)
        if not m:
            raise RuntimeError(f"Missing array {name}")
        return np.fromstring(m.group(1), sep=",")

    Tbins = grab_array("Tbins_ARR")
    nHbins = grab_array("nHbins_ARR")

    metal_flat = grab_array("Metal_Cooling_ARR")
    prim_flat = grab_array("H_He_Cooling_ARR")
    metal_cie = grab_array("Metal_Cooling_CIE_ARR")
    prim_cie = grab_array("H_He_Cooling_CIE_ARR")

    nh = nHbins.size
    nt = Tbins.size
    metal = metal_flat.reshape(nt, nh)
    prim = prim_flat.reshape(nt, nh)

    return {
        "Tbins": Tbins,
        "nHbins": nHbins,
        "metal": metal,
        "prim": prim,
        "metal_cie": metal_cie,
        "prim_cie": prim_cie,
    }


def compute_units(cfg: Dict[str, float]) -> Dict[str, float]:
    """Derive CGS conversion factors from the chosen base units."""
    length_cgs = cfg.get("length_cgs", 3.0856775809623245e21)
    mass_cgs = cfg.get("mass_cgs", 3.036951775493658e40)
    time_cgs = cfg.get("time_cgs", 3.15576e13)
    mu = cfg.get("mu", 0.6)

    density_cgs = mass_cgs / length_cgs**3
    velocity_cgs = length_cgs / time_cgs
    temperature_cgs = velocity_cgs**2 * mu * AMU / KB

    nH_unit = density_cgs / AMU
    pressure_cgs = density_cgs * velocity_cgs**2

    return {
        "density_cgs": density_cgs,
        "velocity_cgs": velocity_cgs,
        "temperature_cgs": temperature_cgs,
        "nH_unit": nH_unit,
        "pressure_cgs": pressure_cgs,
        "time_cgs": time_cgs,
        "length_cgs": length_cgs,
    }


def bilinear_interp(
    logT: np.ndarray, logn: np.ndarray, tables: Dict[str, np.ndarray]
) -> Tuple[np.ndarray, np.ndarray]:
    """Bilinear interpolation across the PIE grids in logT–logn space."""
    Tbins = tables["Tbins"]
    nHbins = tables["nHbins"]
    prim = tables["prim"]
    metal = tables["metal"]

    iT = np.clip(np.searchsorted(Tbins, logT, side="right") - 1, 0, Tbins.size - 2)
    jN = np.clip(np.searchsorted(nHbins, logn, side="right") - 1, 0, nHbins.size - 2)

    log_T0 = Tbins[iT]
    log_T1 = Tbins[iT + 1]
    t = (logT - log_T0) / (log_T1 - log_T0)

    log_n0 = nHbins[jN]
    log_n1 = nHbins[jN + 1]
    u = (logn - log_n0) / (log_n1 - log_n0)

    omt, omu = 1.0 - t, 1.0 - u

    prim_PIE = (
        omt * omu * prim[iT, jN]
        + t * omu * prim[iT + 1, jN]
        + omt * u * prim[iT, jN + 1]
        + t * u * prim[iT + 1, jN + 1]
    )
    metal_PIE = (
        omt * omu * metal[iT, jN]
        + t * omu * metal[iT + 1, jN]
        + omt * u * metal[iT, jN + 1]
        + t * u * metal[iT + 1, jN + 1]
    )
    return prim_PIE, metal_PIE


def cooling_heating(
    T: np.ndarray,
    rho: np.ndarray,
    Z: float,
    tables: Dict[str, np.ndarray],
    cfg: Dict[str, float],
    units: Dict[str, float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    nH_unit = units["nH_unit"]
    length_unit = units["length_cgs"]
    pressure_cgs = units["pressure_cgs"]
    time_cgs = units["time_cgs"]

    h_rate = cfg.get("hrate", 0.0)
    h_flag = bool(cfg.get("hscale_flag", 0.0))
    h_norm = cfg.get("hscale_norm", 1.0)
    h_height = cfg.get("hscale_height", 0.0)
    h_radius = cfg.get("hscale_radius", 0.0)
    h_alpha = cfg.get("hscale_alpha", 0.0)
    T_cutoff = cfg.get("T_cutoff", 1e10)

    cooling_unit = pressure_cgs / time_cgs / nH_unit / nH_unit
    heating_unit = pressure_cgs / time_cgs / nH_unit

    Tfloor = 10 ** tables["Tbins"][0]
    Tceil = 10 ** tables["Tbins"][-1]
    nHfloor = 10 ** tables["nHbins"][0]
    nHceil = 10 ** tables["nHbins"][-1]

    nH = X_H * nH_unit * rho
    logT = np.log10(T)
    logn = np.log10(nH)

    m_lowT = (T < Tfloor).astype(float)
    m_Tin = ((T >= Tfloor) & (T <= Tceil)).astype(float)
    m_Nin = ((nH >= nHfloor) & (nH <= nHceil)).astype(float)

    prim_PIE, metal_PIE = bilinear_interp(logT, logn, tables)
    lambda_PIE = m_Tin * m_Nin * (prim_PIE + Z * metal_PIE)

    Tbins = tables["Tbins"]
    idx = np.clip(np.searchsorted(Tbins, logT, side="right") - 1, 0, Tbins.size - 2)
    t = (logT - Tbins[idx]) / (Tbins[idx + 1] - Tbins[idx])
    prim_cie = tables["prim_cie"]
    metal_cie = tables["metal_cie"]
    prim_CIE = prim_cie[idx] + t * (prim_cie[idx + 1] - prim_cie[idx])
    metal_CIE = metal_cie[idx] + t * (metal_cie[idx + 1] - metal_cie[idx])
    lambda_CIE_tab = prim_CIE + Z * metal_CIE

    e1 = np.exp(-1.184e5 / (T + 1.0e3))
    e2 = np.exp(-92.0 / T)
    lambda_lowT = Z * (2.0e-19 * e1 + 2.8e-28 * np.sqrt(T) * e2)

    lambda_CIE = m_Tin * lambda_CIE_tab + (1.0 - m_Tin) * (m_lowT * lambda_lowT)

    if h_flag:
        x1min, x1max = cfg.get("x1min", -0.5), cfg.get("x1max", 0.5)
        x2min, x2max = cfg.get("x2min", -0.5), cfg.get("x2max", 0.5)
        x3min, x3max = cfg.get("x3min", -0.5), cfg.get("x3max", 0.5)
        R = math.sqrt((0.5 * (x1min + x1max)) ** 2 + (0.5 * (x2min + x2max)) ** 2)
        x3v = 0.5 * (x3min + x3max)
        horz_falloff = math.exp(-R / max(h_radius, 1e-12)) if h_radius > 0 else 1.0
        vert_scale2 = h_height**2 + h_alpha * R**2
        vert_falloff = math.exp(-(x3v * x3v) / vert_scale2) if vert_scale2 > 0 else 1.0
        gamma_heating = h_rate * h_norm * X_H * nH_unit * horz_falloff * vert_falloff
    else:
        gamma_heating = h_rate * h_norm * X_H * nH_unit

    m_hot = (T > 1.0e4).astype(float)
    inv_ratio = 1.0e4 / T
    damp_factor = m_hot * inv_ratio**8 + (1.0 - m_hot)
    gamma_heating *= damp_factor

    dx_cgs = (
        (cfg.get("x1max", 1.0) - cfg.get("x1min", 0.0)) / cfg.get("nx1", 1)
    ) * length_unit
    # dx_cgs = 1.0 * length_unit

    neutral_frac = 1.0 - 0.5 * (1.0 + np.tanh((T - 8e3) / 1.5e3))
    tau = neutral_frac * nH * 1.0e-17 * dx_cgs
    frac = np.exp(-tau)

    lambda_cooling = (1.0 - frac) * lambda_CIE + frac * lambda_PIE
    gamma_heating *= 1.0 - frac

    m_cut = (T < T_cutoff).astype(float)
    volumetric_cooling = m_cut * ((X_H * rho) ** 2) * lambda_cooling / cooling_unit
    volumetric_heating = m_cut * gamma_heating * X_H * rho / heating_unit

    return frac, lambda_cooling, volumetric_cooling, volumetric_heating


def main() -> None:
    cfg = USER_CFG.copy()
    tables = load_tables(TABLES_PATH)
    units = compute_units(cfg)

    T = np.logspace(math.log10(TMIN), math.log10(TMAX), NSAMPLES)
    rho = cfg["rho0"] * T_HOT / T

    frac, lambda_cool, cool_vol, heat_vol = cooling_heating(
        T, rho, Z_REL, tables, cfg, units
    )
    net_vol = cool_vol - heat_vol

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    ax = axes.ravel()

    ax[0].plot(T, lambda_cool, label=r"$\Lambda_\mathrm{tot}$", color="C0")
    ax[0].plot(T, -lambda_cool, ls="--", color="C0")

    ax[0].set_yscale("log")
    ax[0].set_xscale("log")

    ax[0].set_ylim(1e-32, 1e-20)

    ax[0].axhline(USER_CFG["hrate"], ls="--")

    ax[0].set_xlabel("Temperature [K]")
    ax[0].set_ylabel(r"Cooling coefficient $\Lambda$ [erg cm$^3$ s$^{-1}$]")
    ax[0].grid(True, which="both", alpha=0.3)
    ax[0].legend()

    ax[1].plot(T, cool_vol, label="Cooling (nH^2 Λ)")
    ax[1].plot(T, heat_vol, label="Heating")

    ax[1].set_ylim(1e-15, None)

    ax[1].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("Temperature [K]")
    ax[1].set_ylabel("Volumetric rate [erg cm$^{-3}$ s$^{-1}$]")
    ax[1].grid(True, which="both", alpha=0.3)
    ax[1].legend()

    nonzero_net = np.abs(net_vol[np.nonzero(net_vol)])
    linthresh = 1e-30 if nonzero_net.size == 0 else np.min(nonzero_net)
    ax[2].plot(T, net_vol, label="Cooling - Heating", color="C3")
    ax[2].axhline(0.0, color="0.4", ls="--", lw=1.0)

    crossing_temps = []
    net_sign = np.sign(net_vol)
    for i in range(1, net_vol.size):
        if net_sign[i] == 0.0:
            crossing_temps.append(T[i])
        elif net_sign[i - 1] == 0.0 or net_sign[i] == net_sign[i - 1]:
            continue
        else:
            log_t0 = np.log10(T[i - 1])
            log_t1 = np.log10(T[i])
            zero_log_t = log_t0 - net_vol[i - 1] * (log_t1 - log_t0) / (
                net_vol[i] - net_vol[i - 1]
            )
            crossing_temps.append(10.0**zero_log_t)

    if crossing_temps:
        ax[2].scatter(
            crossing_temps,
            np.zeros(len(crossing_temps)),
            color="C3",
            s=28,
            zorder=3,
            label="Zero crossings",
        )
        for idx, crossing_temp in enumerate(crossing_temps):
            y_offset = 12 if idx % 2 == 0 else -18
            va = "bottom" if y_offset > 0 else "top"
            ax[2].annotate(
                f"{crossing_temp:.2e} K",
                xy=(crossing_temp, 0.0),
                xytext=(0, y_offset),
                textcoords="offset points",
                ha="center",
                va=va,
                fontsize=8,
            )

    ax[2].set_xscale("log")
    ax[2].set_yscale("symlog", linthresh=linthresh)
    ax[2].set_xlabel("Temperature [K]")
    ax[2].set_ylabel("Net volumetric rate [erg cm$^{-3}$ s$^{-1}$]")
    ax[2].grid(True, which="both", alpha=0.3)
    ax[2].legend()

    ax[3].plot(T, frac, label=r"Shielding factor")

    ax[3].set_xscale("log")
    ax[3].set_xlabel("Temperature [K]")
    ax[3].set_ylabel(r"Shielding factor $f_\mathrm{shield}$")
    ax[3].grid(True, which="both", alpha=0.3)
    ax[3].legend()

    fig.suptitle(f"CGM cooling/heating | rho0={cfg['rho0']:g} cm^-3, Z={Z_REL:g} Z_sun")
    out_path = pathlib.Path(__file__).with_name("cgm_cooling_rates.png")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_path, dpi=200)
    print(f"Saved {out_path}")
    plt.show()


if __name__ == "__main__":
    main()
