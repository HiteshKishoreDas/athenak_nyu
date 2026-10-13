# Forcing convention audit, 2026-10-09

Source commit: 7cdff3618726c3ac386a8bb1efa7ec32f0b7d6a2.
Independent read-only specialist review plus coordinator interpretation.

Schmidt–Federrath 2011 section 4 defines T=L/V. Our actual OU correlation
T=0.5 and integral L=0.5 fix V=1. The previous algebraic V=0.661 rescaling
cannot reinterpret this run while preserving its physical forcing correlation.
It would require T=0.756. Retain the published-unit discrepancy; do not correct
old plots or tune coefficients to manufacture agreement.

The pure-solenoidal Table 4 does not give V/c0 or a complete absolute OU-amplitude
prescription. Table 5's Ma=6.67 applies to mixed forcing and is not transferable.

Primary sources:
- https://arxiv.org/html/1010.4492 section 4 and Tables 4/5.
- https://arxiv.org/pdf/astro-ph/0407616 section 2, equations 7–9: external
  amplitude F0=V^2/L and OU time T=L/V; printed force-RMS relation depends on
  forcing mixture. Do not silently replace this with the later convention.
- https://arxiv.org/pdf/0809.1321 section 2, equations 7–9: mixture normalization
  g_zeta=3/sqrt(1-2zeta+3zeta^2), asymptotic f_rms=3V^2/L. For our scales,
  that is 6. The executed native driver imposes instantaneous volume RMS 3.
  The SGS paper does not explicitly establish this later amplitude convention
  for its Table 4 Enzo realization, so 6 is not a confirmed target.
- https://arxiv.org/pdf/0905.1060 sections 2.1–2.2: forcing amplitude calibrated
  to measured Mach near 5.5, V=cs*Mrms, correlation half-box crossing time;
  parabolic Fourier amplitude on 1<|k|<3 peaked at 2. Absolute amplitude needed
  to reproduce the SGS Enzo realization is not given.

The source SGS paper's equation 10 and section 4.2 identify density-weighted
mean rho*epsilon and say its time averages are in Table 4. Both HTML and PDF
Table 4 headers instead print mean epsilon. Density-weighted interpretation
is better supported by prose, but a printed notation error remains an inference.
A journal endpoint returned 403 and figure screenshots failed in the independent
review; no corrected published table or visual confirmation is claimed.

Conclusion: current run is Mach-matched, with incomplete forcing equivalence.
Hold its forcing and closure coefficients fixed for the reconstruction comparison.
Compare raw K, Mach ratio, production/dissipation and thermal evolution. PPM4
versus PLM retains LLF and isolates reconstruction sensitivity only. The smooth
unforced stationary shear separately measures numerical kinetic loss. Thermal
energy increase minus SGS dissipation in turbulence includes pressure work and
cannot be called numerical heating without an additional pressure-work ledger.
