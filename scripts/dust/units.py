gamma = 5.0 / 3.0
gm1 = gamma - 1.0

length_cgs = 3.0856775809623245e21  # kpc
mass_cgs = 3.036951775493658e39  # n = 0.1 cm^-3
time_cgs = 3.15576e15  # 100 Myr
mu = 0.62  # mean molecular weight
velocity_cgs = length_cgs / time_cgs

km_s_cgs = 1.0e5  # cm/s
vel_to_kms = velocity_cgs / km_s_cgs

k_boltzmann_cgs = 1.3806488e-16  # erg/k
atomic_mass_unit_cgs = 1.67262192369e-24  # g 1.660538921
temp_cgs = velocity_cgs * velocity_cgs * mu * atomic_mass_unit_cgs / k_boltzmann_cgs

rho_to_nH = 0.75
