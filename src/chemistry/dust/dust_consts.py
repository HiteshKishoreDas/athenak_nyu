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