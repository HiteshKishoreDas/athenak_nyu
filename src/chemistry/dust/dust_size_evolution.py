# Evolve the grain size instead of the grain density
# Following Hirashita & Yan 2009
# https://doi.org/10.1111/j.1365-2966.2009.14405.x

import numpy as np
import matplotlib.pyplot as plt

import units as un

um_cgs = 1e-4 * un.cm_cgs  # cm/um
um_to_kpc = um_cgs / un.kpc_cgs  # kpc/um
km_s_to_kpc_myr = (un.km_s_cgs * un.myr_cgs) / un.kpc_cgs  # (kpc/Myr)/(km/s)
kpc3_per_cm3 = (un.kpc_cgs / un.cm_cgs) ** 3  # kpc^3/cm^3

grain_surface_porosity = 1.0  # dimensionless

a = 1e-3  # um
b = 1e-1 # um
N_bins = 50  # dimensionless

T = 1e4  # K
n_H = 100.0  # cm^-3
Z = 1.0  # Z_sun

# Dust grain material density
rho_gr_Si = 3.3
rho_gr_C = 2.2
rho_gr = (rho_gr_Si + rho_gr_C) / 2  # g/cm^3
rho_gr_um = rho_gr * (um_cgs / un.cm_cgs) ** 3  # g/um^3

# Shattering threshold velocity
v_shatt_Si = 2.7  # km/s
v_shatt_C = 1.2  # km/s
v_shatt = 0.5 * (v_shatt_C + v_shatt_Si)  # km/s

gamma_Si = 25
gamma_C = 12
gamma = (gamma_Si + gamma_C) / 2  # dimensionless

# Young's modulus
E_Si = 5.4e11 # dyn/cm^2
E_C = 3.4e10 # dyn/cm^2
E = (E_Si + E_C) / 2  # dyn/cm^2

# Material properties
s_Si = 1.2
s_C = 1.9
s = (s_Si + s_C) / 2  # dimensionless

# Grain sound speed
# Tielens+ 1994, Table 1
# https://ui.adsabs.harvard.edu/abs/1994ApJ...431..321T/abstract
c0_Si = 5.0  # km/s
c0_C = 1.8  # km/s
c0 = (c0_Si + c0_C) / 2  # km/s

# Critical Pressure
P1_Si = 3e11  # dyn/cm^2
P1_C = 4e10  # dyn/cm^2
P1 = (P1_Si + P1_C) / 2  # dyn/cm^2

# Critical Pressure of vapourisation
Pv_Si = 5.4e12 # dyn/cm^2
Pv_C = 5.8e12 # dyn/cm^2
Pv = (Pv_Si + Pv_C) / 2  # dyn/cm^2

# Dimensionless number that has to do with jump condition
# at the interface of collision 
R = 1 # dimensionless

# Const in reln for excavation flow
z_const = 3.4  # dimensionless

V_cell = 1e-5 # kpc^3 (kept for diagnostics; evolution below uses number density)
D = 0.5 # Dust to metal ratio, dimensionless
Z_solar = 0.0134 # Solar metallicity, dimensionless
Zgas = Z_solar  # Gas metallicity, dimensionless
rho_d = D * Zgas * n_H  # amu/cm^3 Gas mass density
rho_d *= un.atomic_mass_unit_cgs  # g/cm^3

# Shattered frag distribution
shatt_frag_expo = -3.3  # dimensionless

v_turb = 10.0  # km/s

dt = 1e0 # Myr

dust_bin_edges = np.logspace(np.log10(a), np.log10(b), num=N_bins + 1)  # um
# centers
bin_mult = dust_bin_edges[1] / dust_bin_edges[0]  # dimensionless
dust_bin = dust_bin_edges[:-1] * np.sqrt(bin_mult)  # um
mass_bin = (4 * np.pi / 3) * rho_gr_um * dust_bin**3  # g

dust_bin_sizes = dust_bin_edges[1:] - dust_bin_edges[:-1]  # um

n_dust_cm3 = rho_d / (4/3 * np.pi * rho_gr * (dust_bin * um_cgs)**3)  # grains/cm^3

# Number distribution normalized to n_dust_cm3.
dust_dist_shape = (dust_bin / dust_bin[0]) / dust_bin_sizes  # 1/um
dust_norm = n_dust_cm3 / np.sum(dust_dist_shape * dust_bin_sizes)  # grains/cm^3
dust_dist_init = dust_dist_shape * dust_norm  # grains/um/cm^3


# !!!!!! CHECK UNITS !!!!!!!!

def v_coag_fn():

    F_stick = 10  # dimensionless

    A, B = np.meshgrid(dust_bin, dust_bin)
    A_cgs = A * um_cgs
    B_cgs = B * um_cgs

    vc = 2.14  # cgs prefactor
    vc *= F_stick
    vc *= np.sqrt((A_cgs**3 + B_cgs**3) / (A_cgs + B_cgs) ** 3)
    vc *= gamma ** (5 / 6)
    vc /= E ** (1 / 3)
    vc /= (A_cgs * B_cgs / (A_cgs + B_cgs)) ** (5 / 6)
    vc /= rho_gr**0.5

    return vc / un.km_s_cgs  # km/s

#!!!! Make it conserve mass !!!!
def rebin(dist, bin, edges, shift):

    under, over = 0.0, 0.0  # grains/cm^3
    # shift_back = -shift


    # * Add ghost bins for over and underflow
    dist_pad = np.zeros(np.shape(dist)[0] + 2, dtype=float)  # 1/um
    shifted_bin_pad = np.ones(np.shape(bin)[0] + 2, dtype=float)
    edges_pad = np.ones(np.shape(bin)[0] + 2, dtype=float)
    shifted_edges_pad = np.zeros(np.shape(edges)[0] + 2, dtype=float)

    dist_pad[1:-1] = dist
    shifted_edges_pad[1:-1] = edges + shift  # um
    shifted_bin_pad[1:-1] = bin + shift  # um

    shifted_edges_pad[0] = a * 1e-4 + shift
    shifted_edges_pad[-1] = 1e6 * b + shift
    shifted_bin_pad[0] = a * 1e-2 + shift
    shifted_bin_pad[-1] = b * 1e3 + shift

    edges_pad = shifted_edges_pad - shift  # um
    bin_sizes = edges_pad[1:] - edges_pad[:-1]  # um

    def exact_integral(x1, x2, mass=False):
        # Integrate piecewise-constant dist_pad on shifted bins.
        if x2 <= x1:
            return 0.0

        i1 = np.searchsorted(shifted_edges_pad, x1, side="right") - 1
        i2 = np.searchsorted(shifted_edges_pad, x2, side="left") - 1

        i1 = np.clip(i1, 0, len(dist_pad) - 1)
        i2 = np.clip(i2, 0, len(dist_pad) - 1)

        if i1 == i2:
            return dist_pad[i1] * (x2 - x1)

        next_left_edge = shifted_edges_pad[i1 + 1]
        prev_right_edge = shifted_edges_pad[i2]

        mult=np.ones_like(shifted_bin_pad)  # dimensionless
        if mass:
            mult = 4 * np.pi * rho_gr_um * (shifted_bin_pad**3) / 3  # g/um^3

        integral = 0.0
        integral += dist_pad[i1] * (next_left_edge - x1) * mult[i1]
        integral += dist_pad[i2] * (x2 - prev_right_edge) * mult[i2]
        if i2 > i1 + 1:
            integral += np.sum(dist_pad[i1 + 1 : i2] * bin_sizes[i1 + 1 : i2] * mult[i1 + 1 : i2])
        return integral

    shift_dist = np.zeros_like(dist_pad)  # 1/um
    for i in range(1, len(edges_pad) - 1):
        left_edge = edges_pad[i]
        right_edge = edges_pad[i + 1]
        shift_dist[i] = exact_integral(left_edge, right_edge) / bin_sizes[i]

    under = exact_integral(shifted_edges_pad[0], edges_pad[1], mass=True)
    over = exact_integral(edges_pad[-2], shifted_edges_pad[-1], mass=True)

    return shift_dist[1:-1], under, over

def sputtering() -> float:
    """
    Calculate the grain-size rate due to sputtering.
    Draine
    Parameters:
    a (float): Current grain size in micrometers.
    T (float): Gas temperature in Kelvin.
    n_H (float): Hydrogen number density in cm^-3.
    Z (float): Metallicity relative to solar.
    dt (float): Time step in Myr.

    Returns:
    float: Grain-size rate in micrometers per Myr.
    """
    da_dt = -1.0  # um/Myr
    da_dt *= n_H / 1.0  # cm^-3
    da_dt /= 1 + (2e6 / T) ** 2.5  # K
    da_dt *= Z / 1.0  # Solar metallicity
    da_dt *= grain_surface_porosity ** (-2 / 3)

    return da_dt  # micrometers per Myr

def accretion() -> float:
    """
    Calculate the grain-size rate due to accretion.
    Hirashita+2011

    Exact expression from Hirashita+ 2022, Eqn 22

    Parameters:
    a (float): Current grain size in micrometers.
    T (float): Gas temperature in Kelvin.
    n_H (float): Hydrogen number density in cm^-3.
    Z (float): Metallicity relative to solar.
    dt (float): Time step in Myr.

    Returns:
    float: Grain-size rate in micrometers per Myr.
    """
    da_dt = 0.1 / 537  # micrometers per Myr
    da_dt *= n_H / 1e3  # cm^-3
    da_dt /= (T / 10) ** 0.5  # K
    da_dt *= Z / 1.0  # Solar metallicity
    da_dt *= grain_surface_porosity ** (-2 / 3)

    return da_dt  # micrometers per Myr

def cross_sec() -> np.ndarray:
    """
    Calculate the cross-sectional area for collisions between two grains.

    Returns:
    np.ndarray: Pairwise cross-sectional area matrix in kpc^2.
    """

    A, B = np.meshgrid(dust_bin, dust_bin)
    A_kpc = A * um_to_kpc
    B_kpc = B * um_to_kpc

    return np.pi * (A_kpc + B_kpc) ** 2

def maxwell_tail_mean(v: float):
    
    v0 = v_turb * np.sqrt(2 / 3)  # km/s

    vmean = np.sqrt(8 / np.pi) * v0  # km/s
    vmean *= (1 + 0.5 * (v / v0) ** 2)  # km/s
    vmean *= np.exp(-0.5 * (v / v0) ** 2)  # km/s

    return vmean

def maxwell_head_mean(v: float):
    
    v0 = v_turb * np.sqrt(2 / 3)  # km/s

    vmean = np.sqrt(8 / np.pi) * v0  # km/s
    vmean *= (1 + 0.5 * (v / v0) ** 2)  # km/s
    vmean *= np.exp(-0.5 * (v / v0) ** 2)  # km/s

    vmean = np.sqrt(8 / np.pi) * v0 - vmean  # km/s

    return vmean

def alpha():
    alp = cross_sec()  # kpc^2
    return alp

def sigma_shatt(M):

    M_inv = 1 / M  # dimensionless
    sig = 0.3
    sig *= (s + M_inv - 0.11)**0.13
    sig /= s + M_inv - 1
    
    return sig
    
def M_r (v):
    return v / c0  # dimensionless

def M_1 ():
    c0_cgs = c0 * un.km_s_cgs  # cm/s
    phi_1 = P1 / (rho_gr * c0_cgs**2)  # dimensionless
    
    M1 = 2* phi_1
    M1 /= 1 + np.sqrt(1+4*s*phi_1)

    return M1

def M_shocked (vrel, m_p):
    # Shocked mass

    sigma_r = sigma_shatt(M_r(vrel)/(1+R))
    sigma_1 = sigma_shatt(M_1())
    
    M_sh = (1+2*R)/(1+R)**(9/16)
    M_sh /= sigma_r**(1/9)
    M_sh *= (M_r(vrel)/(sigma_1*M_1()))**(8/9)

    return M_sh * m_p  # g

def a_bound_proj(M_proj, M_targ, v_rel):
    # Minimum and maximum grain size in shattered distrubition
    # in the projectile fragments

    sigma_r = sigma_shatt(M_r(v_rel)/(1+R))
    sigma_1 = sigma_shatt(M_1())

    # Eqn 16 of Hirashita & Yan 2009, with some rearrangement
    v_cat = c0  # km/s
    v_cat *= sigma_1**0.5
    v_cat *= sigma_r**(1/16)
    v_cat *= (1+R)*M_1()
    v_cat *= (M_targ/(1+2*R)/M_proj)**(9/16)

    a_p = (3*M_proj/(4*np.pi*rho_gr_um))**(1/3)
    a_max = 0.22 * a_p
    a_max *= (v_cat/v_rel)

    return 0.03 * a_max, a_max  # um, um

def a_bound_ejec(M_eject):
    # Minimum and maximum grain size in shattered distrubition
    # in the ejected fragments from the target

    a_max = 3/4/np.pi/rho_gr_um
    a_max *= z_const + 1
    a_max /= 4*z_const**3
    a_max /= z_const-2
    a_max *= M_eject

    a_max = a_max**(1/3)

    a_min = a_max * (P1/Pv)**1.47

    return a_min, a_max  # um, um

def shatt_frag_dist_norm(m_frac, a_min, a_max):
    # Normalization constant for the shattered fragment distribution
    # m_frac is the fraction of mass shattered

    expo = shatt_frag_expo+4
    
    Anorm = 3 * expo
    Anorm *= m_frac
    Anorm /= 4 * np.pi * rho_gr_um
    Anorm /= a_max**expo - a_min**expo

    return Anorm

def shatt_frag_integral(a1, a2, a_min, a_max,  Anorm, mass=False):
    # Integrate the shattered fragment distribution from a1 to a2

    ai = max(a1, a_min)
    af = min(a2, a_max)

    expo_extra = 0
    mult = 1.0
    if mass:
        expo_extra = 3
        mult = 4 * np.pi * rho_gr_um / 3

    expo = (shatt_frag_expo+1) + expo_extra
    num_add = Anorm / (expo+1)
    num_add *= af**expo - ai**expo
    num_add *= mult

    
    return num_add * (af > ai)

def shatt_coag(dist, dt):

    underflow, overflow = 0.0, 0.0  # grains/cm^3

    dist_mat_shatt = np.zeros((N_bins, N_bins), dtype=float)  # unordered pair interactions
    dist_mat_coag = np.zeros((N_bins, N_bins), dtype=float)  # unordered pair interactions

    bin_sizes = dust_bin_edges[1:] - dust_bin_edges[:-1]  # um
    num_in_bin = dist * bin_sizes * kpc3_per_cm3  # grains/kpc^3
    alp = alpha()

    fA, fB = np.meshgrid(num_in_bin, num_in_bin)
    pair_mat = np.triu(fA * fB)
    pair_mat[np.diag_indices(N_bins)] *= 0.5
    dist_mat_shatt += pair_mat
    dist_mat_coag += pair_mat

    # n_safe = np.clip(num_in_bin, 0.0, np.inf)
    # pair_counts = np.outer(n_safe, n_safe)
    # pair_counts = np.clip(pair_counts, 0.0, np.finfo(float).max)
    # dist_mat_shatt *= pair_counts
    # dist_mat_coag *= pair_counts

    v_shatt_kpc_myr = maxwell_tail_mean(v_shatt) * km_s_to_kpc_myr  # kpc/Myr
    v_coag_kpc_myr = maxwell_head_mean(v_coag_fn()) * km_s_to_kpc_myr  # kpc/Myr

    dist_mat_shatt *= alp * v_shatt_kpc_myr * dt
    dist_mat_coag *= alp * v_coag_kpc_myr * dt
    # dist_mat_shatt = np.nan_to_num(dist_mat_shatt, nan=0.0, posinf=np.finfo(float).max, neginf=0.0)
    # dist_mat_coag = np.nan_to_num(dist_mat_coag, nan=0.0, posinf=np.finfo(float).max, neginf=0.0)

    # Conservative limiter: keep removals below available source grains.
    remove_try = (
        np.sum(dist_mat_shatt + dist_mat_coag, axis=1)
        + np.sum(dist_mat_shatt + dist_mat_coag, axis=0)
    )
    valid = remove_try > 0
    if np.any(valid):
        frac = np.ones_like(remove_try)
        frac[valid] = num_in_bin[valid] / remove_try[valid]
        global_scale = min(1.0, 0.95 * np.min(frac[valid]))
        dist_mat_shatt *= global_scale
        dist_mat_coag *= global_scale

    # Multiply num_in_bin by dist_mat_shatt and dist_mat_coag to get
    # the exact number of interaction 
    
    #* ======= Coagulation =======

    # Remove coagulated grains
    num_in_bin -= np.sum(dist_mat_coag, axis=1) + np.sum(dist_mat_coag, axis=0)

    # Add coagulated grains
    for i in range(N_bins):
        for j in range(N_bins):
            a_coag = (dust_bin[i]**3 + dust_bin[j]**3)**(1/3)
            k = np.searchsorted(dust_bin_edges, a_coag) - 1
            if k < N_bins:
                num_in_bin[k] += dist_mat_coag[i, j]
            # else:
                # overflow+= dist_mat_coag[i, j]

    #* ======= Shattering =======

    # Remove shattered grains
    num_in_bin -= np.sum(dist_mat_shatt, axis=1) + np.sum(dist_mat_shatt, axis=0)

    mA, mB = np.meshgrid(mass_bin, mass_bin)  # g, g

    # Projectile mass 
    Mproj = np.minimum(mA, mB)  # g

    # Target mass
    Mtarg = np.maximum(mA, mB)  # g

    M_sh = M_shocked(maxwell_tail_mean(v_shatt), Mproj)  # g

    fully_shattered = M_sh > 0.5 * Mtarg  # bool

    M_frac = Mtarg * (fully_shattered)  # g
    M_frac += 0.4 * M_sh * (~fully_shattered)  # g

    # Leftover mass in the target grain after shattering
    M_left = Mtarg - M_frac  # g
    a_left = (3 * M_left / (4 * np.pi * rho_gr_um)) ** (1 / 3)  # um

    a_left_bin = np.searchsorted(dust_bin_edges, a_left) 

    amin_proj, amax_proj = a_bound_proj(Mproj, Mtarg, maxwell_tail_mean(v_shatt))
    amin_ejec, amax_ejec = a_bound_ejec(M_frac)

    Anorm_proj = shatt_frag_dist_norm(Mproj, amin_proj, amax_proj)
    Anorm_ejec = shatt_frag_dist_norm(M_frac, amin_ejec, amax_ejec)

    for i in range(N_bins):
        for j in range(N_bins):

            if dist_mat_shatt[i, j] > 0:

                # Add the leftover grain in the appropriate bin
                if a_left_bin[i, j] == N_bins:
                    overflow += dist_mat_shatt[i, j]* M_left[i, j]
                elif a_left_bin[i, j] == 0:
                    underflow += dist_mat_shatt[i, j] * M_left[i, j] 
                else:    
                    num_in_bin[a_left_bin[i, j]] += dist_mat_shatt[i, j]
                
                # Add the shattered fragments in the appropriate bins
                for k in range(N_bins):
                    a1 = dust_bin_edges[k]
                    a2 = dust_bin_edges[k+1]

                    num_add_proj = shatt_frag_integral(a1, a2, amin_proj[i, j], amax_proj[i, j], Anorm_proj[i, j])
                    num_add_ejec = shatt_frag_integral(a1, a2, amin_ejec[i, j], amax_ejec[i, j], Anorm_ejec[i, j])

                    #! Check!
                    num_in_bin[k] += dist_mat_shatt[i, j]*(num_add_proj + num_add_ejec)


                    # add_factor = np.nan_to_num(num_add_proj + num_add_ejec, nan=0.0, posinf=0.0, neginf=0.0)
                    # num_in_bin[k] += dist_mat_shatt[i, j] * add_factor

                ## Check if mass is conserved for this interaction
                # Shattered mass added
                mass_added = shatt_frag_integral(dust_bin_edges[0], dust_bin_edges[-1], amin_proj[i, j], amax_proj[i, j], Anorm_proj[i, j], mass=True)
                mass_added += shatt_frag_integral(dust_bin_edges[0], dust_bin_edges[-1], amin_ejec[i, j], amax_ejec[i, j], Anorm_ejec[i, j], mass=True)

                mass_underflow = Mtarg[i,j] + Mproj[i,j] - M_left[i,j] - mass_added
                if np.isnan(mass_underflow) or mass_underflow< 0:
                    print(f"Warning: too much mass added, {mass_added}, for interaction ({i}, {j})")
                    print(f"Mass before: {Mtarg[i,j] + Mproj[i,j]}, mass after: {M_left[i,j] + mass_added}, \nleftover mass: {M_left[i,j]}, shattered mass added: {mass_added}")

                    mass_added = 0.0

                underflow += mass_underflow * dist_mat_shatt[i, j]

                

    num_in_bin_cm3 = num_in_bin / kpc3_per_cm3  # grains/cm^3
    return num_in_bin_cm3 / bin_sizes, underflow / kpc3_per_cm3, overflow / kpc3_per_cm3


#* ================ The Run ==========================

da_dt_sputtering = sputtering()  # um/Myr
da_dt_accretion = accretion()  # um/Myr

print(f"Change in grain size due to sputtering: {da_dt_sputtering:.2e} micrometers/Myr")
print(f"Change in grain size due to accretion: {da_dt_accretion:.2e} micrometers/Myr")

net_da_dt = da_dt_accretion + da_dt_sputtering  # um/Myr

sub_dt = dt * 0.01  # Myr
shifted1 = net_da_dt * sub_dt  # um

print(f"Shift in grain size after {sub_dt:.4f} Myr: {shifted1:.5f} micrometers \n\n")

shatt = dust_dist_init.copy()  # grains/um/cm^3

overflow, underflow = 0.0, 0.0  # grains/cm^3

for i in range(50):
    shatt , under_add, over_add = rebin(
        dist=shatt,
        bin=dust_bin,
        edges=dust_bin_edges,
        shift=shifted1,
    )
    overflow += over_add
    underflow += under_add

    shatt, under_add, over_add = shatt_coag(shatt, sub_dt)
    overflow += over_add
    underflow += under_add

dust_dist2 = dust_dist_init.copy()  # grains/um/cm^3
overflow2, underflow2 = 0.0, 0.0  # grains/cm^3

for i in range(50):
    dust_dist2, under_add, over_add = rebin(
        dist=dust_dist2,
        bin=dust_bin,
        edges=dust_bin_edges,
        shift=shifted1,
    )
    overflow2 += over_add
    underflow2 += under_add

    # dust_dist2, under_add, over_add = shatt_coag(dust_dist2, sub_dt)
    # overflow2 += over_add
    # underflow2 += under_add

plt.figure()
plt.plot(dust_bin, dust_dist_init, label="Initial")
plt.plot(dust_bin, dust_dist2, label="Just shifted", linestyle="dashed")
plt.plot(dust_bin, shatt, label=f"All processes")
plt.xlabel("Grain size (micrometers)")
plt.ylabel("Number density")
plt.title(f"{underflow} underflow, {overflow:.2f} overflow")
plt.legend()
plt.xscale("log")
plt.yscale("log")
# plt.ylim(1e-13, None)
plt.grid()
plt.show()

plt.figure()
plt.plot(dust_bin, (dust_dist2-shatt)/shatt, label="Diff", linestyle="dashed")
plt.xlabel("Grain size (micrometers)")
plt.ylabel("Number density")
plt.title(f"{underflow} underflow, {overflow:.2f} overflow")
plt.legend()
plt.xscale("log")
plt.yscale("linear")
# plt.ylim(1e-13, None)
plt.grid()
plt.show()

print("Underflow (grains/cm^3):", underflow)
print("Overflow (grains/cm^3):", overflow)
print("Underflow (grains/cm^3):", underflow2)
print("Overflow (grains/cm^3):", overflow2)

mult = 4 * np.pi * rho_gr_um * (dust_bin**3) / 3 # mass per grain [g]
mult /= un.atomic_mass_unit_cgs  # convert mass of grains in amu
# dust_dist * mult gives grain mass density spectrum [amu / (um cm^3)]

plt.figure()
plt.plot(dust_bin, dust_dist_init*mult, label="Initial")
plt.plot(dust_bin, dust_dist2*mult, label="Shifted", linestyle="dashed")
plt.plot(dust_bin, shatt*mult, label=f"After shatt and coag")
plt.xlabel("Grain size (micrometers)")
plt.ylabel("Density")
plt.title(f"{underflow} underflow, {overflow:.2f} overflow")
plt.legend()
plt.xscale("log")
plt.yscale("log")
# plt.ylim(1e-13, None)
plt.grid()
plt.show()
