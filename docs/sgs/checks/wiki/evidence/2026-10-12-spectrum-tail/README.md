# Direct measured-spectrum SGS estimate — 2026-10-12

The requested estimate fits the amplitude of a Kolmogorov spectrum to measured
resolved velocity spectra, then integrates its extrapolated unresolved tail.
This is distinct from the earlier dissipation-amplitude estimate. No new native
simulations, solver changes, coefficient changes or rebuilds were needed.

For Fourier shell mode n, E(n) = A n^(-5/3), and

    A = exp(mean(log(E(n)) + (5/3) log(n)))
    k_unresolved = integral[nc,infinity] A n^(-5/3) dn
                 = (3/2) A nc^(-2/3).

The primary fitting interval is shells 4–10, retained from the earlier spectral
diagnostic. Other intervals show fit sensitivity. nc=N/2 is the ideal spherical
Nyquist cutoff, with box length one. Mode-space integration already includes
the change of variable from physical wavenumber: there is no extra 2pi factor.
Shell sums have width one and already contain the kinetic-energy factor 1/2.
No additional 1/2 is applied. Amplitude comes only from the measured velocity
spectrum, never SGS energy or dissipation. The rho^(1/3)u spectrum is not used:
its shell sum is not kinetic energy.

## Results

Both PPM4 runs have rho0=V=1. Velocity specific energy multiplied by rho0
approximates energy density; all values below are in units rho0 V².
Resolved spectra average nine integer-T snapshots from 2T through 10T.
For a matching average, SGS energy is independently reread from those same
18 native checkpoints, checking their recorded hashes, times, mean mass and
domain volume. The original time-weighted SGS means are also retained.

| Grid | Tail, shells 4–10 | SGS, same nine snapshots | SGS, time-weighted 2T–10T |
| --- | ---: | ---: | ---: |
| 32³ | 0.0299808 | 0.0592459 | 0.0590680 |
| 64³ | 0.0353572 | 0.0374237 | 0.0367689 |

| Fitting shells | 32³ tail, nc=16 | 64³ tail, nc=32 |
| --- | ---: | ---: |
| 4–6 | 0.0555159 | 0.0491582 |
| 4–8 | 0.0393984 | 0.0407708 |
| 4–10 | 0.0299808 | 0.0353572 |
| 6–10 | 0.0215889 | 0.0294707 |

The primary 64³ estimate is 5.52% below the matched SGS average. Under this
cutoff convention, fit-window estimates span 0.02947–0.04916 and encompass
the evolved SGS energy. They do not establish a factor-of-two SGS deficit.
The 32³ primary estimate is roughly half the model energy, illustrating the
limitations of fitting a strongly damped coarse spectrum.

## Assumptions and limits

- There is no broad flat compensated range. Free slopes over 4–10 are -3.664
  at 32³ and -2.746 at 64³. Fixing slope -5/3 across these curved spectra is
  an extrapolation assumption, not evidence that Kolmogorov scaling holds.
  Fitting-window spread is sensitivity, not a statistical confidence interval.
- Halving the assumed cutoff to N/4 increases every tail by 2^(2/3)=1.5874.
  At 64³, the primary estimate becomes 0.0561261; the 4–6 estimate becomes
  0.0780339. The actual finite-volume filter has no verified sharp cutoff.
  This sensitivity is shown explicitly; a cutoff was not calibrated to the
  paper value 0.0726.
- The FFT grid is a cube and contains resolved corner modes outside the
  ideal Nyquist sphere. The spherical tail is an approximation; corner-shell
  energy is recorded in results.json and shown faintly in the figure.
- A velocity spectrum gives volume-weighted specific kinetic energy. Exact
  compressible Favre SGS energy requires density–velocity correlations and
  a defined filter. Multiplication by rho0 is not that exact construction.
- The integral assumes an indefinitely extended power-law tail. A finite
  dissipation cutoff nd replaces the tail with (3/2)A[nc^(-2/3)-nd^(-2/3)].
- The nine spectra and full-window SGS average use different temporal weights;
  the matched nine-snapshot comparison avoids that mismatch. Neither is an
  independent-seed convergence measurement.

This directly answers the spectral-tail question, without settling the published
benchmark discrepancy or replacing a same-forcing resolved Favre reference.

## Checks and reproduction

Synthetic exact-power-law amplitude recovery and independent quadrature after
the change of variable y=(nc/n)^(1/3) pass at relative tolerance 1e-14. Original
spectral Parseval checks pass at maximum relative error 4.34e-16.
An initial script attempt assumed rho0 was present in both summary protocols;
the 64³ summary omitted it. The final script checks density0, mean mass and unit
volume directly from every native checkpoint. No scientific inputs changed.

On the Mac, using the existing plotting environment:

```bash
MPLCONFIGDIR=/private/tmp/sgs-mpl-cache /private/tmp/sgs-plot-venv/bin/python -B /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-spectrum-tail/analysis.py
```

The raw checkpoint paths and hashes are recorded in results.json. The script
uses the existing native restart reader in the multiphase SGS worktree. Source
copies in the wiki packet are archival; execute the original artifact script
at the path above because its input directory layout is fixed.

![Measured spectra, extrapolated tails and sensitivities](spectrum-tail.png)

[Figure PDF](spectrum-tail.pdf) · [Numerical results](results.json)
