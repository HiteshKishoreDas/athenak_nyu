# SGS energy and Kolmogorov reasoning, 2026-10-12

This is an analysis of completed 32³/64³ PPM4 runs, not a new reference simulation.
The energy shortfall relative to the published table remains unresolved.

## Scaling expectation

For an isotropic, homogeneous, incompressible inertial range, use the radial
three-dimensional spectrum normalized by integral E(kappa) dkappa = specific
kinetic energy = one-half mean velocity squared:

```
E(kappa) = C_K epsilon^(2/3) kappa^(-5/3)
k_sgs = integral[kappa_c,infinity] E(kappa) dkappa
      = (3/2) C_K epsilon^(2/3) kappa_c^(-2/3)
```

A sharp spherical cutoff kappa_c=pi/Delta gives k_sgs proportional to
(epsilon Delta)^(2/3). For density rho0, the same dimensionless estimate is

```
K/(rho0 V²) = (3/2) C_K pi^(-2/3) [Q Delta/L]^(2/3)
Q = epsilon L/V³
```

C_K=1.6 is an illustrative 3D radial constant, with 1.5–1.7 recorded as a chosen
sensitivity range. The one-dimensional longitudinal constant near .5 is not this
3D constant: under local isotropy the conversion is 55/18. This distinction and
the empirical limitations are documented in the primary
[Sreenivasan (1995) paper](https://users.ictp.it/~krs/pdf/1995_003.pdf).
The integral antiderivative and spectral normalization were independently checked.
No additional factor .5 belongs in the SGS production rate.

## Limited measured scale check

Both cases use gamma 1.001, c0=1/5.5, seed 42, instantaneous acceleration RMS 3,
L=T=.5, V=1, RK2/LLF/PPM4, CFL .1, and the published coefficients.
All 19 checked physical/integrator/forcing parameters agree in actual native inputs.
Their domain is a unit cube, hence Delta/L=2/N.

| Grid | Delta/L | Mean K/(rho0 V²),2T–10T |
|---|---:|---:|
| 32³ | 1/16 | .0590679900 |
| 64³ | 1/32 | .0367689001 |

The two-point exponent ln(K32/K64)/ln2 is **.6838905**, near 2/3.
Anchoring a constant-flux Kolmogorov scaling to the32³ measurement predicts
K64=.0372105020, versus measured .0367689001 (about1.19% lower).
This checks relative scale dependence only. It does not establish absolute energy
accuracy or a converged inertial range:32³ is below the paper grid floor, only two
resolutions/one seed exist, and flux/density statistics also change with resolution.
The [SGS paper §4.2](https://arxiv.org/html/1010.4492#S4.SS2) reports a solenoidal
energy exponent .799±.009 across 64³–512³ and explains why density-weighted SGS
energy does not map straightforwardly onto velocity spectra/structure functions.
Its measured law need not equal the incompressible K41 exponent.

## Model-based flux/energy consistency is not independent validation

The SGS dissipation closure gives, with r=Delta/L and mean conserved rho=rho0,

```
Q_SGS = (Cepsilon/r) <(rho/rho0)(k/V²)^(3/2)>
mean K/(rho0 V²) = mass-weighted mean k/V²
mean K/(rho0 V²) <= [Q_SGS r/Cepsilon]^(2/3)
```

The inequality is Jensen convexity of the mass-and-time-weighted moment, not an
assumption that the turbulent field is homogeneous. Cepsilon=1.5, r=1/32,
Q_SGS=.456623089 give homogeneous upper bound .0448963554; the actual mean is
.0367689001. Their implied moment factor
F=<rho k^(3/2)>/[rho0 <K/rho0>^(3/2)] is 1.34926.

Taking the paper's quoted1.236 as density-weighted dissipation yields analogous
F=1.31635. This interpretation is conditional: Table4 labels epsilon while nearby
prose identifies rho epsilon. Its K=.0726 rescaled by our/paper flux ratio^(2/3)
is .0373791612, close to our .0367689001. The measured energy ratio is .50646;
the flux^(2/3) ratio is .51486. Similar moment factors explain their closeness.

**This is not an independent explanation or proof that K is correct.** Our SGS
dissipation is calculated from K via this same closure. The calculation shows
that the low energy and low dissipation are mutually consistent; it does not
identify why the resolved cascade supplies less SGS energy or settle the paper
forcing/normalization discrepancy.

## Conditional absolute spectral estimates

An independent input proxy is the native compression-corrected total heating
(Delta U−C)/duration, including SGS and numerical heating. Its dimensionless
2T–10T value is .654992689. It is measured without using D to define total heat,
but it is not a directly measured inertial-range energy flux. Shock-associated
heating and finite-grid/RK effects remain in it.

| Assumptions/input | Illustrative K41 normalized energy |
|---|---:|
| Sharp kappa_c=pi/Delta, SGS flux .456623, C_K=1.6 | .06582385 |
| Same sharp cutoff, independent total-heat proxy .654993 | .08372124 |
| Gaussian, sigma=Delta/sqrt12, same heat proxy | .10621804 |
| Actual evolved64³ PPM4 SGS energy | .03676890 |

These idealized predictions are larger than the actual energy. They are **not
acceptance targets** for this Mach≈5.5 compressible flow: they assume constant
density, a long inertial range and a specified spectral filter. The implicit
finite-volume filter is not a certified sharp spherical or Gaussian filter.
The Gaussian width is an explicit second-moment convention, not established as
the paper's exact convention. Even the published K=.0726 is below .12784694
obtained from its conditional flux1.236 under the same homogeneous sharp-filter
assumptions. A disagreement with this ideal estimate does not uniquely establish
a bug in either implementation.

For completeness, with Ghat=exp(−sigma² kappa²/2), unresolved mean energy weights
the spectrum by1−Ghat². Integrating the K41 spectrum gives

```
k_sgs = (3/2) C_K Gamma(2/3) (epsilon sigma)^(2/3)
```

Independent logarithmic numerical quadrature, including analytic endpoint tails,
agrees with the Gaussian integral to relative error 5.66e−9. Exact identities and
Jensen-bound checks pass. No model coefficient was changed.

[Four-panel comparison](kolmogorov-comparison.pdf), [numeric results](results.json),
[explicit assumptions](protocol.json), [reproduction script](analysis.py).
Resolved spectra are examined separately; a proper same-forcing high-resolution
Favre-filtered reference remains unexecuted.
