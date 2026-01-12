import numpy as np
import matplotlib.pyplot as plt

P = 5e3
G0 = 1.7

T = np.logspace(1, 6, 100)  # Temperature array from 10^1 to 10^8 K
n = P / T

beta = 0.74 / (T**0.068)

cool_pe = 4.65e-30  # erg s^-1
cool_pe *= T**0.94
cool_pe *= (G0 * np.sqrt(T) / n) ** beta
cool_pe *= n

heat_pe = 1e-24 * G0
_den = G0 * np.sqrt(T) / n
eps = 4.9e-2 / (1 + (_den / 1925) ** 0.73)
eps += (3.7e-2 * (T / 1e4) ** 0.7) / (1 + (_den / 5000))

heat_pe *= eps

heat_ISRF = 2e-26 * n

heat_ki02 = 2e-26  # Heating in KI2002
cool_ki02 = 1e7 * np.exp(-114800 / (T + 1000))  # Cooling in KI2002
cool_ki02 = +14 * np.sqrt(T) * np.exp(-92 / T)

cool_ki02 *= heat_ki02

plt.loglog(T, cool_pe, label="Dust cooling", color="C0")
plt.loglog(T, heat_pe / n, label="Dust heating", color="C1")

plt.loglog(T, cool_ki02, label="KI2002 cooling", linestyle="--", color="C0")
plt.loglog(T, heat_ki02 / n, label="KI2002 heating", linestyle="--", color="C1")

plt.loglog(T, heat_ISRF / n, label="ISRF Gotham", linestyle="-.", color="C2")

plt.xlabel("Temperature (K)")
plt.ylabel(r"$\Lambda, \Gamma/n$ (erg cm$^3$ s$^{-1}$)")
plt.grid(True, which="major", ls="--")
plt.legend()

plt.ylim(1e-28, 1e-22)

# Save the plot
plt.savefig("dust_photoelectric_cooling_rate.png")
