# Resolved spectra: what can be inferred about a Kolmogorov range

Nine integer-T snapshots from 2T through10T for32³/64³ PPM4, with the same
physical input parameters. Spectra are equal-weight snapshot means, not complete
time-weighted2T–10T spectra. Native timestamps and checkpoint hashes are saved.

The analysis subtracts instantaneous Favre mean velocity, transforms either u or
rho^(1/3)u, removes each transformed field's volume mean and takes the full 3D FFT.
Fourier coefficients are normalized by N³. Shell energy is one-half the sum of
all component squared Fourier amplitudes, summed over nearest-integer radial
mode norm. Shell sums preserve Parseval; maximum relative error is 4.34e−16.
No division by shell occupancy is applied, since that would change the spectrum
and its expected slope. The density-weighted variable is a cascade diagnostic;
its spectrum integral is not physical kinetic energy. No unresolved-energy
reference is constructed from the resolved FFT.

Fit shells4–10 were frozen before field analysis. They span less than half a
decade, directly above forcing support1–3 and near numerical damping at these
resolutions. Fits are descriptive, with correlated bins/snapshots and no
independent uncertainty estimate.

| Grid | u spectrum slope | rho^(1/3)u spectrum slope |
|---|---:|---:|
| 32³ | −3.6642 | −3.6413 |
| 64³ | −2.7460 | −2.1797 |

The compensated curves do not show a broad flat K41 range. Spectra steepen toward
the cutoff, particularly32³. The plots are consistent with limited dynamic range
and spatial damping, but do not isolate LLF, reconstruction, SGS stresses, shock
compressibility or forcing as the cause. Lower CFL does not remove a solver's
spatial flux dissipation; the matched CFL experiment separately probes timestep
sensitivity. A near 2/3 two-grid mean-energy exponent alone does not certify a
Kolmogorov cascade or correct SGS amplitude.

[Four spectrum panels](resolved-spectra.pdf), [data and Parseval checks](results.json),
[frozen analysis protocol](protocol.json), [local script](analysis.py).
Slope−5/3 and−2 curves are illustrative guides normalized in the fit range.
Only shells below spherical Nyquist are plotted; the Parseval check includes
all cubical Fourier modes. No fitted tail is reported as SGS ground truth.

The [Schmidt–Federrath paper §4.2](https://arxiv.org/html/1010.4492#S4.SS2) discusses
supersonic spectra and density weighting, and warns that SGS energy is not directly
related to spectra or structure functions. A sufficiently resolved same-forcing
reference, filter convention and resolution sensitivity remain necessary for an
independent SGS-energy validation.
