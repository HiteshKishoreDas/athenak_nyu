# Calculate the distribution of relative velocity magnitude for a gaussian velocity distribution

import numpy as np
import matplotlib.pyplot as plt
import gc

from scipy.stats import maxwell


vx1 = np.random.normal(0, 1, 100000)
vy1 = np.random.normal(0, 1, 100000)
vz1 = np.random.normal(0, 1, 100000)

vx2 = np.random.normal(0, 1, 100000)
vy2 = np.random.normal(0, 1, 100000)
vz2 = np.random.normal(0, 1, 100000)

veq = np.random.normal(0, np.sqrt(2), 100000)

v1 = np.sqrt(vx1**2 + vy1**2 + vz1**2)
v2 = np.sqrt(vx2**2 + vy2**2 + vz2**2)

v_rel = np.sqrt((vx1 - vx2) ** 2 + (vy1 - vy2) ** 2 + (vz1 - vz2) ** 2)

plt.figure()
plt.hist(vx1 * np.sqrt(2), bins=100, density=True, alpha=0.2, label="$\sqrt{2}v_{x1}$")
plt.hist(vx2 * np.sqrt(2), bins=100, density=True, alpha=0.2, label="$\sqrt{2}v_{x2}$")
plt.hist(veq, bins=100, density=True, alpha=0.8, label="$v_{eq}$", histtype="step")
plt.hist(vx1 - vx2, bins=100, density=True, alpha=0.2, label="$v_{x1}-v_{x2}$")

plt.xlabel("Velocity Component")
plt.ylabel("Probability Density")
plt.title("Distribution of Velocity Components")
plt.legend()
plt.show()

del vy1, vz1, vx2, vy2, vz2
gc.collect()


plt.figure()
plt.hist(
    v_rel / np.sqrt(2), bins=100, density=True, alpha=0.2, label="$|v_1-v_2|/\sqrt{2}$"
)
plt.hist(v1, bins=100, density=True, alpha=0.2, label="$|v_1|$")
plt.hist(v2, bins=100, density=True, alpha=0.2, label="$|v_2|$")
plt.hist(
    v_rel / np.sqrt(2),
    bins=100,
    density=True,
    cumulative=True,
    histtype="step",
    alpha=0.8,
    label="$|v_1 - v_2|/\sqrt{2}$ (CDF)",
)
plt.hist(
    v1,
    bins=100,
    density=True,
    cumulative=True,
    histtype="step",
    alpha=0.8,
    label="$|v_1|$ (CDF)",
)
plt.axhline(0.5, ls=":", color="k", label="CDF=0.5")
plt.axvline(np.sqrt(2), ls=":", color="k")


plt.xlabel("Relative Velocity")
plt.ylabel("Probability Density")
plt.title("Distribution of Relative Velocity Magnitude")
plt.legend(loc="lower right")

plt.show()
