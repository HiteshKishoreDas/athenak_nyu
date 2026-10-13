# Forcing-only normalization and covariance oracle

Executed 2026-10-11. **14 checks pass; exact paper forcing equivalence remains
unestablished.** This independent Python calculation reproduces the audited native
forcing algebra. It does not execute AthenaK, evolve hydrodynamics, validate SGS
production, or change any coefficients or native files.

Run locally with the existing NumPy/Matplotlib environment:

```sh
/private/tmp/sgs-plot-venv/bin/python forcing_oracle.py
```

`protocol.json` records the analytical expectations, seed, sample counts,
tolerances and source hashes **before random samples are drawn**. `results.json`
records every error and threshold. Deterministic scaled tolerances are 1e-12;
independent ensemble second-moment and covariance tolerances are 0.02; energy
variance tolerance is 0.03. These thresholds were not changed after sampling.
All statistical samples use a fixed seed and existing dependencies.

The oracle checks transverse projectors, exact signed-support ensemble isotropy,
direct 8³ Fourier rendering against Parseval, constant instantaneous RMS,
density-weighted mean removal, independent Gaussian second and fourth moments,
and the exact discrete OU recurrence. Native independent +k and -k coefficients
must be folded together before applying Parseval. This is explicitly tested.

For the same native radial spectrum and an ensemble acceleration target of 3,
65,536 independent Gaussian fields give `sqrt(mean instantaneous RMS²)=2.99992`.
The fixed modal scale leaves instantaneous RMS scatter **0.18029**. Native
instantaneous normalization instead fixes RMS to exactly 3 and has zero scatter.
Normalized covariance differs from the prescribed diagonal 1/3 by at most 0.000243
for modal scaling and 0.000230 for instantaneous normalization.

![Four forcing-only diagnostic panels](forcing-audit.png)

- Top left: orange is the fluctuating RMS of fixed modal-amplitude OU fields;
  blue marks the native fixed RMS. Their common second moment does not make the
  stochastic laws identical. The comparison deliberately uses the **same native
  radial spectrum** to isolate normalization; orange is not an exact Enzo run.
- Top right: the first 10T of a cheap, one-seed 200T modal recursion. Blue has fixed
  norm; orange fluctuates. This is a forcing sequence, not a turbulence simulation.
- Bottom left: temporal autocorrelations from the 200T sequence compared with
  exp(−lag/T). These curves are descriptive and correlated; their closeness does
  not establish equivalence or an SGS error bound.
- Bottom right: peak-normalized per-mode radial amplitudes. Native uses
  sqrt(1−(n−2)²)×2/n, whose final amplitude maximum is **n=1.5**, even though the
  parabola parameter is npeak=2. Federrath 2010 describes a parabolic amplitude.
  Current public Enzo uses its square; the displayed bandwidth 0.5 is a chosen
  shape comparison giving support 1<n<3, not an established Table 4 input.

The native nonuniform-density fixture uses rho=1+0.4cos(4πx), overlapping an active
n=2 mode. Subtracting the mass-weighted mean acceleration cancels net momentum
injection but can introduce a nonzero volume mean. Instantaneous renormalization
and this density-dependent subtraction mean the applied field is not simply the
original Gaussian modal OU process. The ensemble comparison panels use uniform
density; they do not quantify this additional feedback in a turbulent flow.

## Convention map and provenance

[Schmidt–Federrath §4](https://arxiv.org/html/1010.4492) defines integral L as half
the box and correlation T=L/V. Our actual L=T=0.5 fix V=1. V=0.661 cannot
reinterpret the same forcing; it would require T=0.756. The pure-solenoidal
Table 4 does not specify V/c0 or complete absolute forcing-amplitude parameters.
Table 5's mixed-forcing Mach scale must not be transferred to Table 4. Its
dissipation header says mean epsilon whereas §4.2 identifies mean rho*epsilon;
the density-weighted interpretation is supported by prose, not a verified
correction of the published header.

[Schmidt 2009 §2, equations 7–9](https://arxiv.org/pdf/0809.1321) explicitly
states asymptotic acceleration RMS=3V²/L. This would equal 6 for our L,V.
It is a conditional convention reference, **not a confirmed Table 4 target**.

[Federrath 2010 §§2.1–2.2](https://arxiv.org/pdf/0905.1060) instead states
V=cs*Mrms, a half-box crossing correlation time, a parabolic Fourier amplitude
on 1<n<3 with peak 2, and amplitude adjusted to Mach near 5.5. It does not supply
the complete absolute amplitude for the SGS paper's Enzo realization.

Pinned public Enzo source tree: `18812cdde8851e07592733bda26f7575554ba52c`:

- [Initialization](https://github.com/enzo-project/enzo-dev/blob/18812cdde8851e07592733bda26f7575554ba52c/src/enzo/StochasticForcing_Init.C):
  L=box/alpha, T=L/V; fixed modal normalization
  3V²/[L sqrt(1−2zeta+3zeta²)], with normalized spectral amplitudes.
- [Evolution](https://github.com/enzo-project/enzo-dev/blob/18812cdde8851e07592733bda26f7575554ba52c/src/enzo/StochasticForcing_Evolve.C):
  exact exp(−dt/T) drift and sqrt(1−exp(−2dt/T)) innovations.
- [Renderer](https://github.com/enzo-project/enzo-dev/blob/18812cdde8851e07592733bda26f7575554ba52c/src/enzo/Grid_FTStochasticForcing.C):
  independent half-space modes, sqrt(2) Fourier prefactor, no instantaneous
  normalization. Together these give analytical ensemble RMS=3V²/L.

This implementation's headers record an Enzo 2.4 adaptation in 2014; public
[file introduction](https://github.com/enzo-project/enzo-dev/commit/d368b31235441402748ec7bf9ac819b6a4af7040)
is October 2014. It is useful primary implementation evidence but does not certify
the paper's Enzo 1.5 source/input/seed. The contemporary stochastic-forcing example
is a Grete 2014/2017 MUSCL/MHD test, not the paper parameter file.

Native audit: `src/srcterms/turb_driver.cpp` at source commit
`7bf794adeb45cda1e3640c52d5ae52e5f00da12f`; file SHA256 is recorded in protocol.
Lines 748–750 define the radial spectrum, 828–844 the modal OU recurrence,
1039–1045 the mass-weighted subtraction, and 1114–1119 instantaneous RMS scaling.
There is no ensemble-modal normalization option in the audited native driver.

## Minimal justified follow-up

Preserve the existing 10T continuation and its frozen forcing. Use this oracle to
document the normalization and spectral departures. A future explicitly separate
modal-OU case could reuse the native modal state and exact OU recurrence, applying
one fixed scale from the analytical spectral covariance rather than normalizing
each realization. That would test a specified public convention, not prove exact
Table 4 reproduction. Restarts must encode any such changed semantics; coefficient
tuning, relabeling acceleration RMS 6 as the paper target, or changing an ongoing
campaign would invalidate the controlled comparison.

Ponytail review: no forcing framework, new dependency, native rewrite or duplicated
simulation runner was added. The standalone oracle is an inspectable artifact;
its random temporal curves are descriptive, and its pass status covers only its
prespecified algebra/ensemble checks.

Execution note: the first run passed numerical checks but Matplotlib fell back
from an unwritable default cache. Cache locations were then set inside this
artifact directory. Formatting and the density-overlap fixture were corrected;
the final script passes Ruff E/F and all 14 checks. No numerical thresholds or
random seeds were adjusted. The final figure was visually inspected.
