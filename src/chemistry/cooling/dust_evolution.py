import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import cmasher as cm

a_min = 1e-7  # cm
a_max = 2e-4  # cm
N_dustbins = 2
dustbins = np.logspace(np.log10(a_min), np.log10(a_max), N_dustbins + 1)

dust_size = 0.5 * (dustbins[:-1] + dustbins[1:])

Z_solar = 0.0134
D_sat = Z_solar

dt = 1
tlim = 200
N_time = int(tlim / dt)
N_T = 100

P = 1e5
# P = 5e3
# P = 1e2
T = np.logspace(1, 8, N_T)  # Temperature array from 10^1 to 10^8 K
rho_g = P / T

init_dust = 0.5
rho_d = init_dust * D_sat * rho_g


def accretion(rho_d, rho_g, dust_size, T, Z):

    t_acc = 200  # Myr
    t_acc *= dust_size / 1e-7  # cm
    t_acc *= (rho_g / 20) ** -1  # cm^-3
    t_acc *= (T / 50) ** -0.5  # K
    t_acc *= Z  # Z_s

    D = rho_d / rho_g
    t_grow = t_acc / (1 - D / D_sat)

    return rho_d / t_grow


def sputtering(rho_d, rho_g, dust_size, T, Z):

    t_sp = 70  # Myr
    t_sp *= dust_size / 1e-7  # cm
    t_sp *= (rho_g / 1e-3) ** -1  # cm^-3
    t_sp *= 1 + (T / 2e6) ** -2.5

    return -rho_d / t_sp


def shattering(rho_d, rho_g, dust_size, T, Z):
    # From Dubois 2024

    D = rho_d / rho_g
    s_i = 3  # grain material density in g cm^-3

    sig_DL = 10  # km s^-1

    t_sh = 54  # Myr
    t_sh *= dust_size / 1e-7  # cm
    t_sh *= s_i / 3  # cm^-3
    t_sh *= (D / 0.01) ** -1
    t_sh *= (rho_g / 1) ** -1  # cm^-3
    t_sh *= (sig_DL / 10) ** -1  # km s^-1

    return rho_d / t_sh


def coagulation(rho_d, rho_g, dust_size, T, Z):
    # From Dubois 2024

    D = rho_d / rho_g
    s_i = 3  # grain material density in g cm^-3
    F = 0.5  # Fudge factor from

    sig_DS = 10  # km s^-1

    t_cog = 0.27  # Myr
    t_cog *= dust_size / 5e-9  # cm
    t_cog *= s_i / 3  # cm^-3
    t_cog /= F
    t_cog *= (D / 0.01) ** -1
    t_cog *= (rho_g / 1e3) ** -1  # cm^-3
    t_cog *= (sig_DS / 0.1) ** -1  # km s^-1

    return rho_d / t_cog


def sne_agb_inj(rho_d, rho_g, dust_size, T, Z):
    # Aoyama 2017

    # averaged dust condensation efficiency in stella rejecta
    f_in = 0.3

    E_Z = Z_solar  # per solar mass of stars
    # Metal injection rate from SNe and AGB stars
    # Can we calculated cell-by-cell

    m_ej_SN = 8.72  # solar mass

    d_inj_per_SN = f_in * E_Z * m_ej_SN  # solar mass

    return d_inj_per_SN


# This should already taken care of by sputtering
def sne_destruction(rho_d, rho_g, dust_size, T, Z):
    # Aoyama 2017

    pass


def plot_rates():

    dust = np.ones((N_dustbins, np.shape(rho_d)[0]), dtype=float)
    dust *= rho_d

    accrn_rhs = np.zeros_like(dust)
    sputt_rhs = np.zeros_like(dust)
    shatt_rhs = np.zeros_like(dust)
    coag_rhs = np.zeros_like(dust)

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    for i in range(N_dustbins):
        accrn_rhs = accretion(dust[i], rho_g, dust_size[i], T, Z_solar)
        sputt_rhs = -sputtering(dust[i], rho_g, dust_size[i], T, Z_solar)
        shatt_rhs = shattering(dust[i], rho_g, dust_size[i], T, Z_solar)
        coag_rhs = coagulation(dust[i], rho_g, dust_size[i], T, Z_solar)

        ax[0, 0].loglog(T, accrn_rhs, label=f"Size {dust_size[i]*1e4:.2f} $\mu$m")
        ax[0, 1].loglog(T, sputt_rhs, label=f"Size {dust_size[i]*1e4:.2f} $\mu$m")
        ax[1, 0].loglog(T, shatt_rhs, label=f"Size {dust_size[i]*1e4:.2f} $\mu$m")
        ax[1, 1].loglog(T, coag_rhs, label=f"Size {dust_size[i]*1e4:.2f} $\mu$m")

    ax[0, 0].set_xlabel("Temperature (K)")
    ax[0, 0].set_ylabel(r"Dust accretion rate (d0 Myr$^{-1}$)")
    ax[0, 0].grid(True, which="major", ls="--")

    ax[0, 1].set_xlabel("Temperature (K)")
    ax[0, 1].set_ylabel(r"Dust sputtering rate (d0 Myr$^{-1}$)")
    ax[0, 1].grid(True, which="major", ls="--")
    ax[0, 1].legend()

    ax[1, 0].set_xlabel("Temperature (K)")
    ax[1, 0].set_ylabel(r"Dust shattering rate (d0$ Myr$^{-1}$)")
    ax[1, 0].grid(True, which="major", ls="--")

    ax[1, 1].set_xlabel("Temperature (K)")
    ax[1, 1].set_ylabel(r"Dust coagulation rate (d0 Myr$^{-1}$)")
    ax[1, 1].grid(True, which="major", ls="--")


def dust_evolution(sne_rate=0.01, dt=0.1, tlim=10, N_T=100):

    dust = np.ones((N_dustbins, np.shape(rho_d)[0]), dtype=float)
    dust *= rho_d

    def rhs(t, dust_flat):
        # * RHS calculation
        dust_local = dust_flat.reshape(N_dustbins, rho_d.shape[0])

        rhs = np.zeros_like(dust_local)
        shatt_rhs = np.zeros(rho_d.shape[0], dtype=float)
        coag_rhs = np.zeros(rho_d.shape[0], dtype=float)

        D_tot = np.sum(dust_local / rho_g, axis=0)

        for i in range(N_dustbins):
            rhs[i] += accretion(dust_local[i], rho_g, dust_size[i], T, Z_solar) * (
                D_tot < D_sat
            )
            rhs[i] += sputtering(dust_local[i], rho_g, dust_size[i], T, Z_solar)

        shatt_rhs += shattering(dust_local[1], rho_g, dust_size[1], T, Z_solar)
        coag_rhs += coagulation(dust_local[0], rho_g, dust_size[0], T, Z_solar)

        rhs[0] += shatt_rhs - coag_rhs
        rhs[1] += coag_rhs - shatt_rhs

        return rhs.reshape(-1)

    t_eval = np.arange(0, tlim, dt)
    sol = solve_ivp(
        rhs,
        [0, tlim],
        dust.reshape(-1),
        t_eval=None,
        vectorized=True,
    )

    return sol


# ================================
# ================================

dust_evolution_sol = dust_evolution(dt=dt, tlim=tlim, N_T=N_T)

dust = dust_evolution_sol.y

dust /= np.hstack((rho_g, rho_g))[:, None]
dust /= D_sat * init_dust

dust1 = dust[:N_T, :]
dust2 = dust[N_T:, :]

# dust1 = dust[0] / rho_g / D_sat
# dust2 = dust[1] / rho_g / D_sat

vmin = -0.1
vmax = 0.1

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(14, 8))

im1 = ax[0, 0].imshow(
    np.log10(dust1),
    # np.log10(rho_stack),
    aspect="auto",
    extent=[0, tlim, np.log10(1e1), np.log10(1e8)],
    origin="lower",
    vmin=vmin,
    vmax=vmax,
    cmap=cm.redshift,
)
ax[0, 0].set_xlabel("Time (Myr)")
ax[0, 0].set_ylabel("log10 T(K)")
fig.colorbar(im1, ax=ax[0, 0], label="log10 Dust mass (d0 Myr)")

ax[0, 1].plot(
    dust_evolution_sol.t,
    init_dust * dust1[10, :],
    label=f"Size {dust_size[0]*1e4:.2f} $\\mu$m, T idx 10",
)
ax[0, 1].plot(
    dust_evolution_sol.t,
    init_dust * (dust1[90, :]),
    label=f"Size {dust_size[0]*1e4:.2f} $\\mu$m, T idx 90",
)
ax[0, 1].set_xlabel("Time (Myr)")
ax[0, 1].set_ylabel("Log Dust mass (d0 Myr)")
ax[0, 1].legend()

im2 = ax[1, 0].imshow(
    np.log10(dust2),
    # np.log10(dust_evolution_sol.y),
    aspect="auto",
    extent=[0, tlim, np.log10(1e1), np.log10(1e8)],
    origin="lower",
    vmin=vmin,
    vmax=vmax,
    cmap=cm.redshift,
)
fig.colorbar(im2, ax=ax[1, 0], label="log10 Dust mass (d0 Myr)")
ax[1, 0].set_xlabel("Time (Myr)")
ax[1, 0].set_ylabel("log10 T(K)")

line_i = np.argmin(np.abs(np.logspace(1, 8, N_T) - 1e4))

ax[1, 1].plot(
    dust_evolution_sol.t,
    init_dust * (dust2[line_i, :]),
    label=f"Size {dust_size[1]*1e4:.2f} $\\mu$m, $T=10^4$K",
)

line_i = np.argmin(np.abs(np.logspace(1, 8, N_T) - 1e6))
ax[1, 1].plot(
    dust_evolution_sol.t,
    init_dust * (dust2[line_i, :]),
    label=f"Size {dust_size[1]*1e4:.2f} $\\mu$m, T idx 90",
)
ax[1, 1].set_xlabel("Time (Myr)")
ax[1, 1].set_ylabel("Log Dust mass (d0 Myr)")

ax[1, 1].legend()

time_step = dust_evolution_sol.t[1:] - dust_evolution_sol.t[:-1]
plt.figure()
plt.plot(time_step)
plt.xlabel("step#")
plt.ylabel("time step (Myr)")


plot_rates()
