import numpy as np

import units as un
from dust_consts import *
from ic import *

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

def shattering(dist, dt):

    underflow, overflow = 0.0, 0.0  # grains/cm^3

    dist_mat_shatt = np.zeros((N_bins, N_bins), dtype=float)  # unordered pair interactions

    bin_sizes = dust_bin_edges[1:] - dust_bin_edges[:-1]  # um
    num_in_bin = dist * bin_sizes * kpc3_per_cm3  # grains/kpc^3
    alp = cross_sec()

    fA, fB = np.meshgrid(num_in_bin, num_in_bin)
    pair_mat = np.triu(fA * fB)
    pair_mat[np.diag_indices(N_bins)] *= 0.5
    dist_mat_shatt += pair_mat

    v_shatt_kpc_myr = maxwell_tail_mean(v_shatt) * km_s_to_kpc_myr  # kpc/Myr

    dist_mat_shatt *= alp * v_shatt_kpc_myr * dt

    # Conservative limiter: keep removals below available source grains.
    remove_try = (
        np.sum(dist_mat_shatt, axis=1)
        + np.sum(dist_mat_shatt, axis=0)
    )
    valid = remove_try > 0
    if np.any(valid):
        frac = np.ones_like(remove_try)
        frac[valid] = num_in_bin[valid] / remove_try[valid]
        global_scale = min(1.0, 0.95 * np.min(frac[valid]))
        dist_mat_shatt *= global_scale

    # Multiply num_in_bin by dist_mat_shatt and dist_mat_coag to get
    # the exact number of interaction 
    
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
