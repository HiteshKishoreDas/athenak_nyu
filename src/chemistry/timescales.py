# -*- coding: utf-8 -*-
"""Timescale estimator and contour plotter.

Generates contour plots of equilibrium timescales (in Myr) on a log T vs. log nH
grid for several default photo backgrounds. Outputs are PNGs saved under
src/chemistry/test_outputs/.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, replace
from typing import Callable, Dict
import cmasher as cmr

try:
    import numpy as np
except ImportError as exc:  # pragma: no cover - runtime guard
    raise SystemExit(
        "numpy is required for timescale contours. Install with `pip install numpy`."
    ) from exc

try:
    import matplotlib.pyplot as plt
except ImportError as exc:  # pragma: no cover - runtime guard
    raise SystemExit(
        "matplotlib is required for timescale contours. Install with `pip install matplotlib`."
    ) from exc
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path

if __package__ in (None, ""):  # allow running as a script from a repo checkout
    _src_dir = Path(__file__).resolve().parents[1]
    if str(_src_dir) not in sys.path:
        sys.path.insert(0, str(_src_dir))

from chemistry.util.voronov import get_voronov_coeff

# Solar composition constants (Asplund+ 2009-like) for a simple metal electron term
SOLAR_X = 0.7381  # hydrogen mass fraction
SOLAR_Z = 0.0134  # metals mass fraction
METAL_A = 16.0  # assume oxygen-like metals
METAL_ZNUM = 8.0
K_B_EV = 8.617333262e-5  # Boltzmann constant in eV/K

try:
    from chemistry.cooling.cooling_fn import tcool_calc
    from chemistry.cooling import units as cool_units
except ImportError:  # pragma: no cover - optional dependency
    tcool_calc = None
    cool_units = None


def alpha_B(T: float) -> float:
    """Case B recombination coefficient for H (cm^3 s^-1), Osterbrock fit."""
    return 2.59e-13 * (T / 1.0e4) ** -0.7


def alpha_A(T: float) -> float:
    """Case A recombination coefficient for H (cm^3 s^-1), rough fit.

    For the purposes of timescale comparisons, a simple approximation is fine.
    This keeps alpha_A > alpha_B at fixed T.
    """
    return 4.2e-13 * (T / 1.0e4) ** -0.7


def collisional_ionization_coeff(T: float) -> float:
    """Approximate collisional ionization coefficient for H (cm^3 s^-1).

    Voronov 1997-like fit: C ~ 5.85e-11 * sqrt(T) * exp(-E0/kT) / (1 + sqrt(T/1e5))
    with E0 = 13.6 eV.
    """
    k_B_eV = 8.617333262e-5  # eV/K
    E0 = 13.6  # eV
    sqrt_T = math.sqrt(T)
    return 5.85e-11 * sqrt_T * math.exp(-E0 / (k_B_eV * T)) / (1.0 + sqrt_T / 316.23)


def hydrogenic_alpha_B_scaled(T: float, z_eff: float) -> float:
    """Scale the hydrogenic Case B alpha by effective charge."""
    if z_eff <= 0:
        return alpha_B(T)
    T_eff = T / (z_eff * z_eff)
    return (z_eff * z_eff) * alpha_B(T_eff)


def hydrogenic_alpha_A_scaled(T: float, z_eff: float) -> float:
    """Scale the hydrogenic Case A alpha by effective charge."""
    if z_eff <= 0:
        return alpha_A(T)
    T_eff = T / (z_eff * z_eff)
    return (z_eff * z_eff) * alpha_A(T_eff)


def collisional_ionization_coeff_with_E0(
    T: float, ionization_potential_eV: float
) -> float:
    """Voronov-like collisional ionization coeff with custom ionization energy."""
    k_B_eV = 8.617333262e-5  # eV/K
    sqrt_T = math.sqrt(T)
    return (
        5.85e-11
        * sqrt_T
        * math.exp(-ionization_potential_eV / (k_B_eV * T))
        / (1.0 + sqrt_T / 316.23)
    )


def collisional_ionization_voronov(
    T: float, A: float, P: float, X: float, K: float, E0_eV: float
) -> float:
    """Exact Voronov (1997) collisional ionization rate (cm^3 s^-1).

    C(T) = A * (1 + P*sqrt(U)) / (X + U) * U^K * exp(-U), U = E0 / (k_B T)
    """
    if T <= 0.0:
        return 0.0
    U = E0_eV / (K_B_EV * T)
    if U <= 0.0:
        return 0.0
    denom = X + U
    if denom <= 0.0:
        return 0.0
    return A * math.exp(-U) * (1.0 + P * math.sqrt(U)) / denom * (U**K)


def _verner_rr(
    T: float, A: float, B: float, T0: float, T1: float, T2: float = 0.0
) -> float:
    """Verner & Ferland (1996)-style radiative recombination fit."""
    t = math.sqrt(T / T0)
    tp = math.sqrt(T / T1)
    denom = t * (1.0 + t) ** (1.0 - B) * (1.0 + tp) ** (1.0 + B)
    if T2 > 0.0:
        denom *= 1.0 + math.sqrt(T / T2)
    return A / denom


def alpha_rr_civ(T: float) -> float:
    """Toy radiative recombination for CIV -> CIII (cm^3 s^-1).

    Rough magnitude (few x 1e-12 over 1e4-1e6 K); hand-tuned, not a literature
    fit. Replace with CHIANTI/Badnell coefficients for production use.
    """
    return _verner_rr(T, A=8.0e-12, B=0.50, T0=3.0e5, T1=1.0e7)


def alpha_dr_civ(T: float) -> float:
    """Toy two-term dielectronic recombination for CIV -> CIII (cm^3 s^-1).

    Coefficients are hand-tuned for the right ballpark; replace with tabulated
    DR fits (e.g., CHIANTI/Badnell) for fidelity.
    """
    if T <= 0:
        return 0.0
    # Coefficients chosen to give a mild DR bump around 1e5-1e6 K.
    term1 = 3.0e-12 * math.exp(-2.5e5 / T)
    term2 = 1.0e-12 * math.exp(-7.0e5 / T)
    return term1 + term2


def alpha_total_civ(T: float) -> float:
    """Total recombination coefficient (RR + DR) for CIV -> CIII."""
    return alpha_rr_civ(T) + alpha_dr_civ(T)


def collisional_ionization_ciii(T: float) -> float:
    """CIII -> CIV (Voronov 1997)."""
    c = get_voronov_coeff("C", 2)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def collisional_ionization_civ(T: float) -> float:
    """CIV -> CV (Voronov 1997)."""
    c = get_voronov_coeff("C", 3)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


@dataclass(frozen=True)
class IonLadder:
    """A 1D ionization ladder for one element."""

    element: str
    nstages: int
    # collisional ionization i -> i+1 (cm^3/s)
    C_up: list[Callable[[float], float]]
    # recombination i -> i-1 (cm^3/s); alpha_down[0] unused
    alpha_down: list[Callable[[float], float]]
    # optional photoionization rate i -> i+1 (s^-1)
    Gamma_up: list[Callable[[float, float], float]] | None = None


def ladder_rate_matrix(
    ladder: IonLadder, T: float, n_e: float, nH: float
) -> np.ndarray:
    """Construct the tri-diagonal rate matrix for df/dt = A f."""
    N = ladder.nstages
    A = np.zeros((N, N), dtype=float)

    def gamma_i(i: int) -> float:
        if ladder.Gamma_up is None:
            return 0.0
        return float(ladder.Gamma_up[i](T, nH))

    for i in range(N):
        # loss out of i to i+1
        Rup = 0.0
        if i < N - 1:
            Rup = n_e * float(ladder.C_up[i](T)) + gamma_i(i)

        # loss out of i to i-1
        Rdown = 0.0
        if i > 0:
            Rdown = n_e * float(ladder.alpha_down[i](T))

        A[i, i] = -(Rup + Rdown)

        if i > 0:
            R_from_below = n_e * float(ladder.C_up[i - 1](T)) + gamma_i(i - 1)
            A[i, i - 1] += R_from_below

        if i < N - 1:
            R_from_above = n_e * float(ladder.alpha_down[i + 1](T))
            A[i, i + 1] += R_from_above

    return A


def ladder_equilibration_time(
    ladder: IonLadder, T: float, n_e: float, nH: float
) -> float:
    """Return the slowest non-zero equilibration time (s) for the ladder."""
    A = ladder_rate_matrix(ladder, T=T, n_e=n_e, nH=nH)
    eig = np.linalg.eigvals(A)

    # For this tri-diagonal birth-death matrix the spectrum is real; discard
    # tiny imaginary parts from numerical noise and use a small floor for zero.
    eig_real = np.real(eig)
    eps_zero = 1e-20
    valid = eig_real[np.abs(eig_real) > eps_zero]
    if valid.size == 0:
        return math.inf

    # Physical modes have non-positive real parts; take the one closest to zero
    # from the negative side. If numerical noise pushes all to positive, bail.
    neg = valid[valid < -eps_zero]
    if neg.size == 0:
        return math.inf
    lam_slow = np.max(neg)  # closest to zero (least negative)
    return 1.0 / abs(lam_slow)


def _ladder_with_photo(ladder: IonLadder, gamma_eff: float | None) -> IonLadder:
    """Attach a uniform photo-up rate if ladder lacks stage-specific Gamma."""
    if gamma_eff is None or ladder.Gamma_up is not None:
        return ladder
    gamma_list = [(lambda T, nH, g=gamma_eff: g) for _ in range(ladder.nstages - 1)] + [
        lambda T, nH: 0.0
    ]
    return replace(ladder, Gamma_up=gamma_list)


def make_local_3stage_ladder(
    element: str,
    C_below_to_mid: Callable[[float], float],
    C_mid_to_above: Callable[[float], float],
    alpha_mid_to_below: Callable[[float], float],
    alpha_above_to_mid: Callable[[float], float],
    Gamma_below_to_mid: Callable[[float, float], float] | None = None,
    Gamma_mid_to_above: Callable[[float, float], float] | None = None,
) -> IonLadder:
    """Convenience helper for a 3-stage ladder [below, mid, above]."""
    C_up = [C_below_to_mid, C_mid_to_above, lambda T: 0.0]
    alpha_down = [lambda T: 0.0, alpha_mid_to_below, alpha_above_to_mid]
    Gamma_up = None
    if Gamma_below_to_mid is not None or Gamma_mid_to_above is not None:
        Gamma_up = [
            Gamma_below_to_mid or (lambda T, nH: 0.0),
            Gamma_mid_to_above or (lambda T, nH: 0.0),
            lambda T, nH: 0.0,
        ]
    return IonLadder(
        element=element, nstages=3, C_up=C_up, alpha_down=alpha_down, Gamma_up=Gamma_up
    )


@dataclass(frozen=True)
class RateModel:
    """Container for ion/transition-specific rate functions and cooling params."""

    name: str
    alpha_caseA: Callable[[float], float]
    alpha_caseB: Callable[[float], float]
    coll_ion_fn: Callable[[float], float]
    alpha_total_fn: Callable[[float], float] | None = None
    cooling_Zsol: float = 1.0
    cooling_lambda_fac: float = 1.0
    atomic_number: float = 1.0
    # solar_abundance_by_number: n_species / n_H at solar metallicity
    solar_abundance_by_number: float | None = None
    z_eff: float = 1.0
    ionization_potential_eV: float = 13.6
    supports_case_AB: bool = True
    # Optional per-species photoionization rate (s^-1) or callable(T, nH)
    photo_rate: float | None = None
    photo_rate_fn: Callable[[float, float], float] | None = None
    ladder: IonLadder | None = None

    def alpha(self, T: float, case: str) -> float:
        if not self.supports_case_AB and self.alpha_total_fn is not None:
            return self.alpha_total_fn(T)
        if case == "A":
            return self.alpha_caseA(T)
        return self.alpha_caseB(T)

    def with_overrides(self, **kwargs) -> "RateModel":
        return replace(self, **kwargs)


HYDROGEN_RATES = RateModel(
    name="H",
    alpha_caseA=alpha_A,
    alpha_caseB=alpha_B,
    coll_ion_fn=collisional_ionization_coeff,
    atomic_number=1.0,
    solar_abundance_by_number=1.0,
    z_eff=1.0,
    ionization_potential_eV=13.6,
)


def alpha_total_cv(T: float) -> float:
    """Toy recombination coeff for CV -> CIV."""
    return hydrogenic_alpha_A_scaled(T, 5.0)


def collisional_ionization_sii_from_siI(T: float) -> float:
    """SiI -> SiII (Voronov 1997)."""
    c = get_voronov_coeff("Si", 0)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def collisional_ionization_siii_from_sii(T: float) -> float:
    """SiII -> SiIII (Voronov 1997)."""
    c = get_voronov_coeff("Si", 1)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def alpha_total_sii(T: float) -> float:
    """Toy recombination coeff for SiII -> SiI."""
    return hydrogenic_alpha_B_scaled(T, 1.0)


def alpha_total_siii(T: float) -> float:
    """Toy recombination coeff for SiIII -> SiII."""
    return hydrogenic_alpha_B_scaled(T, 2.0)


def collisional_ionization_oii_to_oiii(T: float) -> float:
    """OII -> OIII (Voronov 1997)."""
    c = get_voronov_coeff("O", 1)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def collisional_ionization_oiii_to_oiv(T: float) -> float:
    """OIII -> OIV (Voronov 1997)."""
    c = get_voronov_coeff("O", 2)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def collisional_ionization_ci_to_cii(T: float) -> float:
    """CI -> CII (Voronov 1997)."""
    c = get_voronov_coeff("C", 0)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def collisional_ionization_cii_to_ciii(T: float) -> float:
    """CII -> CIII (Voronov 1997)."""
    c = get_voronov_coeff("C", 1)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def collisional_ionization_ov_to_ovi(T: float) -> float:
    """OV -> OVI (Voronov 1997)."""
    c = get_voronov_coeff("O", 4)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def collisional_ionization_ovi_to_ovii(T: float) -> float:
    """OVI -> OVII (Voronov 1997)."""
    c = get_voronov_coeff("O", 5)
    return collisional_ionization_voronov(T, c["A"], c["P"], c["X"], c["K"], c["E0_eV"])


def alpha_total_oiii(T: float) -> float:
    """Toy recombination coeff for OIII -> OII."""
    return hydrogenic_alpha_B_scaled(T, 2.0)


def alpha_total_oiv(T: float) -> float:
    """Toy recombination coeff for OIV -> OIII."""
    return hydrogenic_alpha_B_scaled(T, 3.0)


def alpha_total_cii(T: float) -> float:
    """Toy recombination coeff for CII -> CI."""
    return hydrogenic_alpha_B_scaled(T, 1.0)


def alpha_total_ciii_to_cii(T: float) -> float:
    """Toy recombination coeff for CIII -> CII."""
    return hydrogenic_alpha_B_scaled(T, 2.0)


def alpha_total_ovi(T: float) -> float:
    """Toy recombination coeff for OVI -> OV."""
    return hydrogenic_alpha_B_scaled(T, 5.0)


def alpha_total_ovii_to_ovi(T: float) -> float:
    """Toy recombination coeff for OVII -> OVI."""
    return hydrogenic_alpha_B_scaled(T, 6.0)


CIV_LADDER = make_local_3stage_ladder(
    element="C",
    C_below_to_mid=collisional_ionization_ciii,
    C_mid_to_above=collisional_ionization_civ,
    alpha_mid_to_below=alpha_total_civ,
    alpha_above_to_mid=alpha_total_cv,
)

SIII_LADDER = make_local_3stage_ladder(
    element="Si",
    C_below_to_mid=collisional_ionization_sii_from_siI,
    C_mid_to_above=collisional_ionization_siii_from_sii,
    alpha_mid_to_below=alpha_total_sii,
    alpha_above_to_mid=alpha_total_siii,
)

OIII_LADDER = make_local_3stage_ladder(
    element="O",
    C_below_to_mid=collisional_ionization_oii_to_oiii,
    C_mid_to_above=collisional_ionization_oiii_to_oiv,
    alpha_mid_to_below=alpha_total_oiii,
    alpha_above_to_mid=alpha_total_oiv,
)

CII_LADDER = make_local_3stage_ladder(
    element="C",
    C_below_to_mid=collisional_ionization_ci_to_cii,
    C_mid_to_above=collisional_ionization_cii_to_ciii,
    alpha_mid_to_below=alpha_total_cii,
    alpha_above_to_mid=alpha_total_ciii_to_cii,
)

OVI_LADDER = make_local_3stage_ladder(
    element="O",
    C_below_to_mid=collisional_ionization_ov_to_ovi,
    C_mid_to_above=collisional_ionization_ovi_to_ovii,
    alpha_mid_to_below=alpha_total_ovi,
    alpha_above_to_mid=alpha_total_ovii_to_ovi,
)
# Simple placeholders for common ions; they reuse the hydrogenic fits here so the
# plotting pipeline works end-to-end. Swap in ion-specific rate functions as you
# acquire them.
CIV_RATES = HYDROGEN_RATES.with_overrides(
    name="CIV",
    atomic_number=6.0,
    solar_abundance_by_number=2.69e-4,
    z_eff=4.0,
    ionization_potential_eV=64.5,
    alpha_total_fn=alpha_total_civ,
    alpha_caseA=alpha_total_civ,
    alpha_caseB=alpha_total_civ,
    coll_ion_fn=collisional_ionization_civ,
    supports_case_AB=False,
    ladder=CIV_LADDER,
)
CII_RATES = HYDROGEN_RATES.with_overrides(
    name="CII",
    atomic_number=6.0,
    solar_abundance_by_number=2.69e-4,
    z_eff=2.0,
    ionization_potential_eV=24.4,
    alpha_caseA=alpha_total_cii,
    alpha_caseB=alpha_total_cii,
    coll_ion_fn=collisional_ionization_cii_to_ciii,
    supports_case_AB=False,
    ladder=CII_LADDER,
)
SIII_RATES = HYDROGEN_RATES.with_overrides(
    name="SiII",
    atomic_number=14.0,
    solar_abundance_by_number=3.24e-5,
    z_eff=2.0,
    ionization_potential_eV=16.3,
    alpha_caseA=lambda T: hydrogenic_alpha_A_scaled(T, 2.0),
    alpha_caseB=lambda T: hydrogenic_alpha_B_scaled(T, 2.0),
    coll_ion_fn=collisional_ionization_siii_from_sii,
    supports_case_AB=False,
    ladder=SIII_LADDER,
)
OIII_RATES = HYDROGEN_RATES.with_overrides(
    name="OIII",
    atomic_number=8.0,
    solar_abundance_by_number=4.90e-4,
    z_eff=3.0,
    ionization_potential_eV=54.9,
    alpha_caseA=lambda T: hydrogenic_alpha_A_scaled(T, 3.0),
    alpha_caseB=lambda T: hydrogenic_alpha_B_scaled(T, 3.0),
    coll_ion_fn=collisional_ionization_oiii_to_oiv,
    supports_case_AB=False,
    ladder=OIII_LADDER,
)
OVI_RATES = HYDROGEN_RATES.with_overrides(
    name="OVI",
    atomic_number=8.0,
    solar_abundance_by_number=4.90e-4,
    z_eff=6.0,
    ionization_potential_eV=138.1,
    alpha_caseA=alpha_total_ovi,
    alpha_caseB=alpha_total_ovi,
    coll_ion_fn=collisional_ionization_ovi_to_ovii,
    supports_case_AB=False,
    ladder=OVI_LADDER,
)


def _choose_recomb_case(
    nH: float,
    x: float | None,
    length_cm: float,
    sigma_13p6: float = 6.3e-18,
    f_hi_floor: float = 1e-6,
) -> tuple[str, float]:
    """Pick Case A/B from a simple LyC optical depth estimate.

    tau ~ n_HI * sigma_13p6 * L
    If tau >> 1, recombination photons are trapped -> Case B (on-the-spot).
    If tau << 1, photons escape -> Case A.
    """
    if length_cm <= 0 or nH <= 0:
        return "A", 0.0
    if x is None:
        f_hi = f_hi_floor
    else:
        f_hi = max(1.0 - float(x), f_hi_floor)
    n_hi = nH * f_hi
    tau = n_hi * sigma_13p6 * length_cm
    return ("B" if tau >= 1.0 else "A"), float(tau)


def estimate_timescales(
    T: float,
    nH: float,
    Gamma: float | None,
    x: float | None = None,
    recomb_case: str = "B",
    length_cm: float | None = None,
    metallicity: float = 1.0,
    rate_model: RateModel | None = None,
) -> Dict[str, float]:
    """Return characteristic timescales (seconds) for a single ion/transition.

    metallicity is in units of Z/Z_sun; we approximate metals as oxygen-like.
    """
    model = rate_model or HYDROGEN_RATES
    # Crude electron density estimate: hydrogen contribution only. Metals are trace
    # and their free-electron contribution is ignored here for simplicity.
    n_e = nH * (x if x is not None else 1.0)

    case = recomb_case.upper()
    tau_lyc = 0.0
    if not model.supports_case_AB:
        case = "A"  # metals: use optically thin total recombination
    elif case == "AUTO":
        if length_cm is None:
            raise ValueError("recomb_case='auto' requires length_cm.")
        case, tau_lyc = _choose_recomb_case(nH=nH, x=x, length_cm=length_cm)

    if not model.supports_case_AB and model.alpha_total_fn is not None:
        a = model.alpha_total_fn(T)
    elif case == "A":
        a = model.alpha(T, "A")
    elif case == "B":
        a = model.alpha(T, "B")
    else:
        raise ValueError("recomb_case must be 'A', 'B', or 'auto'.")

    C = model.coll_ion_fn(T)

    gamma_eff: float | None = Gamma
    if gamma_eff is None:
        if model.photo_rate_fn is not None:
            gamma_eff = float(model.photo_rate_fn(T, nH))
        elif model.photo_rate is not None:
            gamma_eff = float(model.photo_rate)
    if gamma_eff is None:
        raise ValueError("Gamma must be provided either explicitly or via rate_model.")

    # Guard against pathological zero/negative values to avoid division errors
    t_photo = math.inf if gamma_eff <= 0 else 1.0 / gamma_eff
    t_recomb = math.inf if n_e <= 0 or a <= 0 else 1.0 / (a * n_e)
    t_coll_ion = math.inf if n_e <= 0 or C <= 0 else 1.0 / (C * n_e)
    denom = gamma_eff + (a + C) * n_e
    t_equil_local = math.inf if denom <= 0 else 1.0 / denom
    t_equil_overall = t_equil_local
    if model.ladder is not None:
        ladder = _ladder_with_photo(model.ladder, gamma_eff)
        t_equil_overall = ladder_equilibration_time(ladder, T=T, n_e=n_e, nH=nH)

    return {
        "recomb_case": case,
        "tau_lyc": tau_lyc,
        "alpha": a,
        "coll_ion_coeff": C,
        "n_e": n_e,
        "t_photo": t_photo,
        "t_recomb": t_recomb,
        "t_coll_ion": t_coll_ion,
        "t_equil": t_equil_overall,
        "t_equil_local": t_equil_local,
        "t_equil_overall": t_equil_overall,
    }


def _fmt_seconds(t: float) -> str:
    if not math.isfinite(t):
        return "inf"
    # Use scientific for wide range
    return f"{t:.3e} s"


def _make_grid(logT_range, logn_range, nT=80, nn=80):
    T_vals = np.logspace(logT_range[0], logT_range[1], nT)
    n_vals = np.logspace(logn_range[0], logn_range[1], nn)
    return np.meshgrid(T_vals, n_vals, indexing="ij")


def _timescale_grid(
    T_grid,
    n_grid,
    Gamma,
    x=None,
    recomb_case="B",
    length_cm=None,
    metallicity: float = 1.0,
    rate_model: RateModel | None = None,
):
    """Compute all timescale grids (seconds) for given Gamma."""
    model = rate_model or HYDROGEN_RATES
    t_photo = np.full_like(T_grid, math.inf, dtype=float)
    t_recomb = np.full_like(T_grid, math.inf, dtype=float)
    t_coll = np.full_like(T_grid, math.inf, dtype=float)
    t_equil_local = np.full_like(T_grid, math.inf, dtype=float)
    t_equil_overall = np.full_like(T_grid, math.inf, dtype=float)
    t_cool = np.full_like(T_grid, math.inf, dtype=float)

    for i in range(T_grid.shape[0]):
        for j in range(T_grid.shape[1]):
            T_ij = float(T_grid[i, j])
            nH_ij = float(n_grid[i, j])
            res = estimate_timescales(
                T=T_ij,
                nH=nH_ij,
                Gamma=Gamma,
                x=x,
                recomb_case=recomb_case,
                length_cm=length_cm,
                metallicity=metallicity,
                rate_model=model,
            )
            t_photo[i, j] = res["t_photo"]
            t_recomb[i, j] = res["t_recomb"]
            t_coll[i, j] = res["t_coll_ion"]
            t_equil_local[i, j] = res["t_equil_local"]
            t_equil_overall[i, j] = res["t_equil_overall"]
            if tcool_calc is not None and cool_units is not None:
                # tcool_calc expects rho in code units; unit_density is 1 mp/cm^3,
                # so rho_code = nH * muH (since nH = rho_code / muH).
                rho_code = nH_ij * float(cool_units.muH)
                tc_code = float(
                    tcool_calc(
                        rho_code,
                        T_ij,
                        Zsol=metallicity,
                        Lambda_fac=model.cooling_lambda_fac,
                    )
                )
                t_cool[i, j] = tc_code * float(cool_units.unit_time)
    return {
        "photo": t_photo,
        "recomb": t_recomb,
        "collisional": t_coll,
        "cooling": t_cool,
        "equilibrium": t_equil_overall,
        "equilibrium_local": t_equil_local,
        "equilibrium_overall": t_equil_overall,
    }


def _plot_timescale_grid(T_grid, n_grid, tau_dict_myr, title, out_path):
    # Typical pressures (P/k) in K cm^-3 for reference isobars.
    isobaric_P_over_k = [
        (100.0, "CGM P/k≈1e2"),
        (5.0e3, "ISM P/k≈5e3"),
    ]

    T_min = float(np.nanmin(T_grid))
    T_max = float(np.nanmax(T_grid))
    T_line = np.logspace(math.log10(T_min), math.log10(T_max), 200)

    panels = [
        {"key": "photo", "title": "Photoionization", "cbar": "log10(t / Myr)"},
        {"key": "recomb", "title": "Recombination", "cbar": "log10(t / Myr)"},
        {
            "key": "collisional",
            "title": "Collisional ionization",
            "cbar": "log10(t / Myr)",
        },
        {"key": "cooling", "title": "Cooling", "cbar": "log10(t / Myr)"},
        {
            "key": "equilibrium",
            "title": "Equilibrium (ladder)",
            "cbar": "log10(t / Myr)",
        },
        {
            "key": "tcool_over_teq",
            "title": "t_cool / t_eq",
            "cbar": "log10(t_cool / t_eq)",
        },
    ]
    if "equilibrium_local" in tau_dict_myr:
        panels.insert(
            5,
            {
                "key": "equilibrium_local",
                "title": "Equilibrium (local stage)",
                "cbar": "log10(t / Myr)",
            },
        )
    nplots = len(panels)
    ncols = 3 if nplots > 4 else 2
    nrows = int(math.ceil(nplots / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharex=True, sharey=True
    )
    axes_flat = np.atleast_1d(axes).ravel()

    for ax, panel in zip(axes_flat, panels):
        key = panel["key"]
        lbl = panel["title"]
        data = tau_dict_myr[key]
        masked = np.ma.masked_invalid(np.where(data <= 0, np.nan, data))
        data_log10 = np.ma.log10(masked)

        finite = data_log10.compressed()
        if finite.size:
            # Use percentiles to avoid a few extreme cells dominating the range.
            vmin = float(np.percentile(finite, 5))
            vmax = float(np.percentile(finite, 95))
            if vmax <= vmin:
                vmax = vmin + 1.0
        else:
            vmin, vmax = -6.0, 6.0

        if key == "tcool_over_teq":
            cmap = plt.get_cmap(cmr.redshift)
            # clim = np.max([np.abs(vmin), np.abs(vmax)])
            # vmin, vmax = -clim, clim
            vmin, vmax = -1, 1

        else:
            cmap = plt.get_cmap(cmr.bubblegum)

        levels = np.linspace(vmin, vmax, 21)

        cs = ax.contourf(
            n_grid,
            T_grid,
            data_log10,
            levels=levels,
            cmap=cmap,
            extend="both",
        )

        # Overlay rough isobaric lines (P/k = const).
        for idx_iso, (P_over_k, label) in enumerate(isobaric_P_over_k):
            n_line = P_over_k / T_line
            ax.plot(
                n_line,
                T_line,
                ls="--",
                lw=1.2,
                color="k" if idx_iso == 0 else "0.45",
                alpha=0.9,
                label=label if ax is axes_flat[0] else "_nolabel_",
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("n_H [cm^-3]")
        ax.set_ylabel("T [K]")
        ax.set_title(lbl)
        ax.grid(True, which="both", ls=":", alpha=0.4)

        cbar = fig.colorbar(cs, ax=ax)
        cbar.set_label(panel["cbar"])

    # Single legend for isobars.
    axes_flat[0].legend(loc="lower right", fontsize=8)

    for ax in axes_flat[len(panels) :]:
        ax.set_visible(False)

    fig.suptitle(title, fontsize=12)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def _plot_tcool_over_teq_summary(
    T_grid,
    n_grid,
    ratio_by_scenario,
    gamma_by_scenario,
    ion_name: str,
    out_path,
    pdf_pages: PdfPages | None = None,
):
    """Plot tcool/teq for all scenarios for a single ion."""
    scenarios = list(ratio_by_scenario.keys())
    nplots = len(scenarios)
    ncols = 2 if nplots > 3 else (1 if nplots == 1 else 2)
    nrows = int(math.ceil(nplots / ncols))
    fig, axes = plt.subplots(
        nrows, ncols, figsize=(5 * ncols, 4 * nrows), sharex=True, sharey=True
    )
    axes_flat = np.atleast_1d(axes).ravel()

    # Isobaric guides (P/k in K cm^-3).
    isobaric_P_over_k = [
        (100.0, "CGM P/k≈1e2"),
        (5.0e3, "ISM P/k≈5e3"),
    ]
    T_min = float(np.nanmin(T_grid))
    T_max = float(np.nanmax(T_grid))
    T_line = np.logspace(math.log10(T_min), math.log10(T_max), 200)

    # Determine global vmin/vmax across scenarios using percentiles for robustness.
    all_vals = []
    for data in ratio_by_scenario.values():
        masked = np.ma.masked_invalid(np.where(data <= 0, np.nan, data))
        log_vals = np.ma.log10(masked).compressed()
        if log_vals.size:
            all_vals.append(log_vals)
    # Fix the colorbar limits to [-1, 1] as requested.
    vmin, vmax = -1.0, 1.0

    for ax, name in zip(axes_flat, scenarios):
        data = ratio_by_scenario[name]
        masked = np.ma.masked_invalid(np.where(data <= 0, np.nan, data))
        data_log10 = np.ma.log10(masked)
        levels = np.linspace(vmin, vmax, 21)
        cs = ax.contourf(
            n_grid,
            T_grid,
            data_log10,
            levels=levels,
            cmap=plt.get_cmap(cmr.redshift),
            extend="both",
        )
        # Overlay rough isobaric lines (P/k = const).
        for idx_iso, (P_over_k, label) in enumerate(isobaric_P_over_k):
            n_line = P_over_k / T_line
            ax.plot(
                n_line,
                T_line,
                ls="--",
                lw=1.0,
                color="k" if idx_iso == 0 else "0.45",
                alpha=0.8,
                label=label if ax is axes_flat[0] else "_nolabel_",
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("n_H [cm^-3]")
        ax.set_ylabel("T [K]")
        gamma_val = gamma_by_scenario.get(name, None)
        if gamma_val is not None:
            ax.set_title(f"{name} (Gamma={gamma_val:.1e} s$^{{-1}}$)")
        else:
            ax.set_title(f"{name}")
        ax.grid(True, which="both", ls=":", alpha=0.4)
        fig.colorbar(cs, ax=ax, label="log10(t_cool / t_eq)")

    for ax in axes_flat[len(scenarios) :]:
        ax.set_visible(False)

    # Single legend for isobars.
    axes_flat[0].legend(loc="lower right", fontsize=8)

    fig.suptitle(f"{ion_name}: Cooling vs. Equilibrium across scenarios", fontsize=12)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_path, dpi=200)
    if pdf_pages is not None:
        pdf_pages.savefig(fig)
    plt.close(fig)


def main():
    # Configuration (no CLI args)
    LOGT_RANGE = (2, 7)  # 10^1 - 10^6 K
    LOGN_RANGE = (-4, 4)  # 10^-4 - 10^4 cm^-3
    X_ION = None  # if not None, n_e = x * nH

    pc_cm = 3.085677581e18
    kpc_cm = 1e3 * pc_cm

    # Each scenario uses a fixed Gamma and a characteristic length scale to choose
    # Case A vs B via a simple LyC optical depth estimate.
    scenarios = {
        # diffuse, optically thin-ish (photons more likely to escape)
        "cgm": {
            "Gamma": 1e-14,  # HI photoionization rate
            "length_cm": 10.0 * kpc_cm,
            # Toy species-specific photo rates; metals are lower than HI.
            "Gamma_by_species": {
                "CII": 3e-15,
                "CIV": 3e-15,
                "SiII": 5e-15,
                "OIII": 8e-15,
                "OVI": 8e-15,
            },
        },
        # diffuse ISM; can be optically thick depending on neutral fraction and scale
        "ism": {
            "Gamma": 1e-12,
            "length_cm": 100.0 * pc_cm,
            "Gamma_by_species": {
                "CII": 3e-13,
                "CIV": 3e-13,
                "SiII": 5e-13,
                "OIII": 8e-13,
                "OVI": 8e-13,
            },
        },
        # compact regions where local trapping is more plausible
        "starforming": {
            "Gamma": 1e-10,
            "length_cm": 10.0 * pc_cm,
            "Gamma_by_species": {
                "CII": 3e-11,
                "CIV": 3e-11,
                "SiII": 5e-11,
                "OIII": 8e-11,
                "OVI": 8e-11,
            },
        },
        "near_star": {
            "Gamma": 1e-9,
            "length_cm": 1.0 * pc_cm,
            "Gamma_by_species": {
                "CII": 3e-10,
                "CIV": 3e-10,
                "SiII": 5e-10,
                "OIII": 8e-10,
                "OVI": 8e-10,
            },
        },
    }

    # Metallicity relative to solar (Z/Z_sun); affects electron density (assume metals are O).
    METALLICITY = 1.0

    # Define ions/transitions to evaluate. Add entries here with the appropriate
    # recombination/collisional rate functions and cooling metallicity scalings.
    species_models = [
        HYDROGEN_RATES,
        CII_RATES,
        CIV_RATES,
        SIII_RATES,
        OIII_RATES,
        OVI_RATES,
    ]

    T_grid, n_grid = _make_grid(LOGT_RANGE, LOGN_RANGE)
    myr_factor = 1.0 / (1.0e6 * 365.25 * 24.0 * 3600.0)

    base_out_dir = Path(__file__).resolve().parent / "test_outputs" / "timescales"
    summary_pdf_path = base_out_dir / "tcool_over_teq_summaries.pdf"
    summary_pdf = PdfPages(summary_pdf_path)

    for species in species_models:
        species_dir = base_out_dir / species.name.lower()
        ratio_by_scenario = {}
        gamma_used = {}
        for name, cfg in scenarios.items():
            gamma = float(cfg.get("Gamma", 0.0))
            gamma_map = cfg.get("Gamma_by_species", {})
            if species.name in gamma_map:
                gamma = float(gamma_map[species.name])
            if gamma <= 0.0:
                raise ValueError(
                    f"No Gamma provided for species {species.name} in {name}"
                )
            length_cm = float(cfg["length_cm"])
            t_dict = _timescale_grid(
                T_grid,
                n_grid,
                Gamma=gamma,
                x=X_ION,
                recomb_case="auto",
                length_cm=length_cm,
                metallicity=METALLICITY,
                rate_model=species,
            )
            tau_myr = {k: v * myr_factor for k, v in t_dict.items()}
            tau_myr["tcool_over_teq"] = np.divide(
                t_dict["cooling"],
                t_dict["equilibrium"],
                out=np.full_like(t_dict["equilibrium"], np.nan, dtype=float),
                where=np.isfinite(t_dict["cooling"])
                & np.isfinite(t_dict["equilibrium"])
                & (t_dict["equilibrium"] > 0),
            )
            ratio_by_scenario[name] = tau_myr["tcool_over_teq"]
            gamma_used[name] = gamma
            recomb_label = (
                "recomb auto"
                if species.supports_case_AB
                else "recomb total (optically thin)"
            )
            title = (
                f"{species.name} Timescales "
                f"(Gamma={gamma:.1e} s$^{{-1}}$, L={length_cm/pc_cm:.0f} pc, "
                f"Z={METALLICITY:g} Z$_{{\\odot}}$; {recomb_label})"
            )
            out_path = (
                species_dir / f"{species.name}_Z{METALLICITY:g}_{name}_timescales.png"
            )
            _plot_timescale_grid(T_grid, n_grid, tau_myr, title, out_path)
            print(f"Wrote {out_path}")
        # Summary tcool/teq plot across scenarios for this ion.
        summary_path = (
            species_dir / f"{species.name}_Z{METALLICITY:g}_tcool_over_teq_summary.png"
        )
        _plot_tcool_over_teq_summary(
            T_grid,
            n_grid,
            ratio_by_scenario,
            gamma_used,
            species.name,
            summary_path,
            pdf_pages=summary_pdf,
        )
        print(f"Wrote {summary_path}")
    summary_pdf.close()
    print(f"Wrote {summary_pdf_path}")


if __name__ == "__main__":
    main()
