# data
# CGS unit per X

cm_cgs = 1.0  # cm
pc_cgs = 3.0856775809623245e18  # cm
kpc_cgs = 3.0856775809623245e21  # cm
g_cgs = 1.0  # g
msun_cgs = 1.98841586e33  # g
atomic_mass_unit_cgs = 1.67262192369e-24  # g 1.660538921
s_cgs = 1.0  # s
yr_cgs = 3.15576e7  # s
myr_cgs = 3.15576e13  # s
cm_s_cgs = 1.0  # cm/s
km_s_cgs = 1.0e5  # cm/s
g_cm3_cgs = 1.0  # g/cm^3
erg_cgs = 1.0  # erg
dyne_cm2_cgs = 1.0  # dyne/cm^2
kelvin_cgs = 1.0  # k

# PHYSICAL CONSTANTS
k_boltzmann_cgs = 1.3806488e-16  # erg/k
grav_constant_cgs = 6.67408e-8  # cm^3/(g*s^2)
speed_of_light_cgs = 2.99792458e10  # cm/s
rad_constant_cgs = 7.56573325e-15  # erg/(cm^3*K^4)
electron_rest_mass_energy_cgs = 5.93e9  # k

length_cgs = 3.0856775809623245e21  # 1 kpc
mass_cgs = 3.036951775493658e40  # Gives n=1 cm^-3 for 1 kpc box
time_cgs = 3.15576e13  # 1 Myr

velocity_cgs = length_cgs / time_cgs

mu = 0.6
gamma = 5 / 3
gm1 = gamma - 1

KELVIN = velocity_cgs * velocity_cgs * mu * atomic_mass_unit_cgs / k_boltzmann_cgs

Eint = 1e6 * 1.0 / KELVIN / gm1
tfloor = 10000.0 / KELVIN

print(f"{Eint = }")
print(f"{tfloor = }")
