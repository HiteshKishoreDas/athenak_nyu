import numpy as np
import units as un

a_small = 0.003
a_large = 0.1

Zsol = 0.02
D_small = 0.5 * Zsol  # small-grain dust fraction
D_large = 0.5 * Zsol  # large-grain dust fraction
D_total = D_small + D_large
rho_grain_cgs = 3.0

D_floor = 1.0e-6  # minimum dust fraction for rates to operate

# Dust *mass densities* in code units (what the conserved scalars store).
# const Real rho_Ds = D_small * rho_gas;
# const Real rho_Dl = D_large * rho_gas;

# Accumulators for net conserved-density changes (code units).
# delta_gas_metal = 0.0   # Δ(ρ * Z_gas)
# delta_small     = 0.0   # Δ(ρ * D_small)
# delta_large     = 0.0   # Δ(ρ * D_large)

# Sound speed in km/s (used by shattering timescale)


def sputtering(rhog, T_K, Zg, Dtot=D_total):
    # =================================================================
    #  1.  THERMAL SPUTTERING  (dust → gas-phase metals)
    #      Tsai & Mathews (1995), McKinnon et al. (2017)
    #
    #      Note : τ_sp = t_sp / a where [a] = μm
    #      τ_sp = 1657 Myr × (10^-3 / nH) × [1 + (T / 2×10^6)^{-5/2}]
    #
    #      Δρ_D = -dt × (3 / a) × ρ_D / τ_sp
    #
    #      The factor 3/a arises because grain mass ∝ a³ and the
    #      absolute erosion rate da/dt is size-independent, so
    #      d(ln M)/dt = 3 (da/dt)/a = 3 / (a τ_sp).
    # =================================================================

    nH = un.rho_to_nH * rhog  # cm^-3

    tau_sput_Myr = 1657.0 * (1.0e-3 / nH)
    T_ratio = T_K / 2.0e6
    tau_sput_Myr *= 1.0 + 1.0 / (T_ratio * T_ratio * np.sqrt(T_ratio))

    tsput_small = (tau_sput_Myr * a_small) / 3.0
    tsput_large = (tau_sput_Myr * a_large) / 3.0

    return tsput_small, tsput_large


def new_accretion(rhog, T_K, Zg, Dtot=D_total):
    # =================================================================
    #  2.  GAS-PHASE ACCRETION  (gas-phase metals → dust)
    #      Popping (2017); Asano et al. (2013)
    #
    #      τ_acc = 150 Myr × (100 / n_H) × √(50 K / T)
    #                      × (Z_sol / Z_gas)
    #                      / [1 - Z_gas / (Z_gas + D_total)]
    #
    #      Δρ_D = +dt × ρ_D / (τ_acc × a)
    # =================================================================

    nH = un.rho_to_nH * rhog  # cm^-3
    fref = 0.35

    alp_LB = 1.0 / (1.0 + 1e-4 * (T_K**1.5))

    tau_acc_Myr = 150.0
    tau_acc_Myr *= 100.0 / nH
    tau_acc_Myr *= np.sqrt(50.0 / T_K)
    tau_acc_Myr *= Zsol / Zg / fref
    tau_acc_Myr /= alp_LB
    tau_acc_Myr /= 1.0 - (D_total / (Zg * fref + D_total))

    taccr_small = 0.0
    taccr_large = 0.0

    cond = Zg > D_floor
    # cond *= T_K < 1.0e4  # accretion only in cold gas
    # cond *= nH > 0.1
    cond *= nH < 1e3  # accretion only in dense gas

    taccr_small = tau_acc_Myr * a_small * cond
    taccr_large = tau_acc_Myr * a_large * cond

    return taccr_small, taccr_large


def accretion(rhog, T_K, Zg, Dtot=D_total):
    # =================================================================
    #  2.  GAS-PHASE ACCRETION  (gas-phase metals → dust)
    #      Popping (2017); Asano et al. (2013)
    #
    #      τ_acc = 150 Myr × (100 / n_H) × √(50 K / T)
    #                      × (Z_sol / Z_gas)
    #                      / [1 - Z_gas / (Z_gas + D_total)]
    #
    #      Δρ_D = +dt × ρ_D / (τ_acc × a)
    # =================================================================

    nH = un.rho_to_nH * rhog  # cm^-3

    tau_acc_Myr = 150.0
    tau_acc_Myr *= 100.0 / nH
    tau_acc_Myr *= np.sqrt(50.0 / T_K)
    tau_acc_Myr *= Zsol / Zg
    tau_acc_Myr /= 1.0 - (D_total / (Zg + D_total))

    taccr_small = 0.0
    taccr_large = 0.0

    cond = Zg > D_floor

    taccr_small = tau_acc_Myr * a_small * cond
    taccr_large = tau_acc_Myr * a_large * cond

    return taccr_small, taccr_large


def shattering(rhog, T_K, Zg, Dtot=D_total):
    # =================================================================
    #  3.  SHATTERING  (large grains → small grains)
    #      Dubois et al. (2024)
    #
    #      τ_shatt = 540 Myr × (1 / n_H) × (ρ_gr / 3 g cm^-3)
    #                        × (0.01 / D_large) × (10 km/s / c_s)
    #
    #      Δρ_Dl = -dt × ρ_Dl / (τ_shatt × a_large)
    # =================================================================

    nH = un.rho_to_nH * rhog  # cm^-3

    cs_kms = un.vel_to_kms * np.sqrt(un.gamma * T_K / un.temp_cgs)

    tau_shatt_Myr = 540.0
    tau_shatt_Myr *= 1.0 / nH
    tau_shatt_Myr *= rho_grain_cgs / 3.0
    tau_shatt_Myr *= 0.01 / D_large
    tau_shatt_Myr *= 10.0 / cs_kms

    return tau_shatt_Myr * a_large


def coagulation(rhog, T_K, Zg, Dtot=D_total):
    # =================================================================
    #  4.  COAGULATION  (small grains → large grains)
    #      Dubois et al. (2024)
    #
    #      τ_coag = 2.7 Myr × (ρ_gr / 3 g cm^-3) × (10^3 / n_H)
    #                        × (0.01 / D_small) × (0.1 km/s / v_coag)
    #                        × f_coag
    #
    #      Δρ_Ds = -dt × ρ_Ds / (τ_coag × a_small / 0.05)
    # =================================================================
    f_coag = 0.5
    v_coag_ref_kms = 0.1

    nH = un.rho_to_nH * rhog  # cm^-3

    tau_coag_Myr = 2.7
    tau_coag_Myr *= rho_grain_cgs / 3.0
    tau_coag_Myr *= 1.0e3 / nH
    tau_coag_Myr *= 0.01 / D_small
    tau_coag_Myr *= 0.1 / v_coag_ref_kms
    tau_coag_Myr /= f_coag

    return tau_coag_Myr * (a_small / 0.005)
