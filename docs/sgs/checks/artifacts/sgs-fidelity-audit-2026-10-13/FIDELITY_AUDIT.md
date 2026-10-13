# Schmidt Federrath SGS implementation fidelity audit

The supported uniform-grid hydrodynamic implementation matches the audited model equations and baseline coefficients. No factor, sign, density convention or RK2 source-weight defect was found. The solver remains unchanged. Numerical regularization and positivity safeguards differ from the original implementation and can affect turbulent results; these differences are identified below.

This audit concerns equation fidelity and deterministic numerical verification. It does not certify agreement with the published forced-turbulence statistics or resolve the measured decay-energy discrepancy.

## Model reference and scope

Reference: Schmidt and Federrath, A&A 528, A106 (2011), [DOI](https://doi.org/10.1051/0004-6361/201015630). Equation numbers and formulas below refer to the [authors' available manuscript](https://arxiv.org/html/1010.4492v1), Eqs. 5–14, Section 3.2 and Section 3.4. The journal and NASA ADS final-PDF endpoints returned HTTP 403; possible changes between the manuscript and final typeset version have not been independently checked. The audited formulas are the same reference used by the implementation and its independent oracles.

Audit date: 2026-10-13. Source branch: `codex/multiphase-schmidt-federrath-hd`; starting commit: `ab54c37d`. Execution: native CPU double precision, serial MPI configuration. The executable and all 324 native source hashes match the previously verified build. No native solver file or default closure coefficient changed during this audit. Term-isolation controls explicitly set selected input coefficients to zero.

The supported scope is Newtonian ideal-gas HD on uniform cubic Cartesian cells with periodic active directions, RK2 and LLF, with DC, PLM or monotone PPM4 reconstruction. The model neglects subgrid gravity and pressure dilatation, as in the reference baseline. Cooling, gravity hooks, AMR/SMR, MHD, relativistic physics and GPU execution are outside this qualification. One- and two-dimensional reductions retain three velocity components and the three-dimensional stress trace.

## Equation to code correspondence

Here K is SGS kinetic energy per volume, k = K/rho is specific SGS energy, A_ij = derivative_j u_i, Q = A_ij A_ij, d = tr(A), and S* = (A+A^T)/2 - d I/3. R denotes the positive dissipation rate per volume.

| Reference term | Implementation | Audit finding |
|---|---|---|
| Eq. 5, energy density and factor 1/2 | `hydro.hpp:82`, `eos/ideal_hyd.cpp:107` and `:165` | Conservative state is K, primitive state is K/rho. K is separate from true passive scalars. No extra factor of two enters conversion. |
| Eq. 8, mixed stress | `hydro/sgs_closure.hpp:28` | tau = 2 C1 Delta sqrt(2 rho K) S* - 2 C2 K A A^T/Q - (2/3)(1-C2)K I. The paper's gradient norm squared is 2Q, so its quadratic numerator factor 2 cancels. |
| Eqs. 5–6, turbulent pressure | `sgs_closure.hpp:67`, `hydro_sgs.cpp:224` | Raw trace is -2K, hence turbulent pressure is 2K/3, independent of gas gamma. Pressure enters through the full traction; it is not added a second time to the thermal Euler flux. |
| Eq. 7, advection and production | `hydro_fluxes.cpp:147`, `hydro_sgs.cpp:262` | K advection uses the existing mass flux times reconstructed k. Production uses the same applied face traction as momentum and energy work. |
| Eq. 10, dissipation | `sgs_closure.hpp:90` | R = Ceps K^(3/2)/(sqrt(rho) Delta). This is an energy-density rate, not a specific-energy rate. |
| Eq. 11, transport | `sgs_closure.hpp:97`, `hydro_sgs.cpp:230` | J = -Ckappa Delta sqrt(rho K) grad(K/rho); K receives -div(J). Diffusion acts on specific energy. It remains in the summed E+K flux; it is not cancelled by an opposite E source. |
| Eq. 13, momentum | `hydro_sgs.cpp:225` | Momentum receives flux -tau.n, equivalent to +div(tau), including the isotropic turbulent-pressure contribution. |
| Eq. 14, resolved energy | `hydro_sgs.cpp:228`, `:340` | Stress-work flux is -u_face dot tau.n. Sources are -Sigma+R in E and +Sigma-R in K. Resolved E excludes K. |
| Eq. 7, compression | `hydro_sgs.cpp:262`, mixed production oracle | The contraction of the full stress already contains -2K div(u)/3. Adding a separate compression term would double count it. |
| Neglected Gamma and rho lambda | Supported configuration and source path | No subgrid gravitational production or pressure-dilatation correction is inserted. Thermal pressure work is a diagnostic, not an additional SGS source. |
| Baseline coefficients | `sgs_closure.hpp:18`, `hydro_sgs.cpp:143` | C1=.02, C2=.7, Ceps=1.5, Ckappa=.65; these match the adopted baseline. |
| Filter length | `hydro_sgs.cpp:193`, `:275` | Delta is the physical cubic cell width dx, not a MeshBlock width, box length or twice dx. |

## Compatible energy exchange and source integration

For an applied face traction t, the discrete production is

```
Sigma_cell = div(u_face dot t) - u_cell dot div(t).
```

The arithmetic face velocity makes this equal to the face-traction times adjacent velocity differences used by `SGSRates`. Thus resolved stress work, momentum work and SGS production use one compatible traction. In the semidiscrete limit, SGS production removes resolved kinetic energy and supplies K without an additional thermal source; R removes K and heats the gas. Finite-step kinetic bookkeeping has a convergent RK error, and the ordinary Euler scheme has numerical heating.

The combined model equation is

```
d(E+K)/dt + div[u(E+K+P) - u dot tau + J] = rho u dot f
```

for the supported periodic, non-self-gravitating, uncooled case. Both Sigma and R cancel in E+K. The sign of diffusion in this summed equation is retained.

`Hydro::Fluxes` computes fluxes and SGS rates from the old stage primitives; `RKUpdate` blends conserved states and applies flux divergence; `HydroSrcTerms` applies the already computed rates before primitive conversion. Each stage therefore evaluates a common RHS. RK2 uses source-update coefficients 1 and 1/2, and the first increment is blended by 1/2 in the second stage. Its final quadrature weights are consequently 1/2 and 1/2. `Driver::SourceTermHistoryWeight` uses these final weights for the source ledger. This is not a source-splitting factor-of-two omission.

Stage-consistent forcing supplies rho*f to momentum and rho*u*f to E using old-stage primitives, with the same final quadrature weights for the forcing ledger. It does not directly supply K. Constant-acceleration quadratic kinetic work emerges from Heun's two-stage update instead of being added as an extra kick per stage.

## Numerical choices that differ from the reference implementation

The reference used Enzo 1.5. This implementation uses AthenaK's RK2/LLF finite-volume machinery. Equation fidelity does not imply identical numerical dissipation, shock thickness, forcing realization or turbulent statistics.

- **LLF and timestep regularization:** the signal estimate uses thermal sound speed squared plus 6k, with related frozen-coefficient, neighbor-density, diffusion and decay restrictions. The 6k term is numerical regularization, not the physical turbulent pressure or a new coefficient in the SGS evolution equation. It affects numerical viscosity.
- **Gradients and face states:** normal face derivatives use adjacent differences; transverse derivatives average centered derivatives. Face rho and K are arithmetic averages. Compatible production uses face work rather than a separate cell-centered stress contraction.
- **Zero-gradient convention:** Eq. 8 has a direction-dependent nonlinear limit at A=0. The implementation explicitly chooses isotropic stress -2K I/3 at exactly zero gradient and zeros closure-only differences below an eight-machine-epsilon conditioning threshold. Nonzero gradients use scaled normalization without a denominator epsilon.
- **Negative-production limiter:** the donor cap scales the entire face traction and its work together. When active, the applied stress need not retain the raw -2K trace. Rejected transfer is recorded. This safeguards positive K but changes the unmodified model locally.
- **Diffusion limiter and fallback:** donor caps can reduce J. An inadmissible full-dt SSP predictor triggers a globally coordinated donor-cell rebuild; a second failure stops. These are numerical safeguards. Kernel and amplitude checks require them inactive; stress/fallback controls exercise them separately.
- **State and initialization policy:** K is appended as a separate active reservoir; the thermal EOS excludes it and cannot silently floor it. Uniform startup seed adds K once, while the verification initializer supplies explicit states. Restarts retain K and source counters. These policies must be matched when comparing benchmarks.

No safeguard was removed or tuned to force agreement with the decay reference. A limited calculation should be described as the published model with the recorded numerical safeguards, rather than an identical Enzo reproduction.

## Executed verification

All six verification suites and the analytic forcing control passed: 134 result records, 2,119 native checks and 174 specialized amplitude/structure checks. Some result records are comparisons or expected configuration rejections rather than new simulations. The original artifact copy of the new diffusion test was also run before the installed repository copy; the count below includes the repository copy only.

| Check | Native checks | Specialized checks | Result |
|---|---:|---:|---|
| Core preservation, conservation, decay, temporal refinement, layouts, restart and safeguard controls | 746 | Included in result records | Pass |
| C1 production and dissipation amplitudes | 144 | 24 | Pass |
| Mixed C2, compression/expansion, traceless strain, signed transfer and energy exchange | 720 | 116 | Pass |
| Thermal pressure work, numerical heating and temporal controls | 205 | 22 | Pass |
| Stage-consistent forcing, analytic constant acceleration and restart contracts | 160 | Included in native checks | Pass |
| New density-varying specific-energy diffusion check | 144 | 12 | Pass |

The independent compiled closure oracle passed 22,575 assertions over 1,025 states. Maximum scaled error was 1.45213e-15 against a frozen 5e-13 tolerance. It checks stress components, symmetry, trace, production, density/energy/length scaling, dissipation, diffusivity, rotations and extreme nonzero gradients without deriving expected values from the production helper.

Uniform decay temporal orders were 1.95053 and 2.00308. Isolated-stress field orders were 2.01919 and 2.01173; each kinetic, thermal and SGS reservoir also approached second order. These establish the tested RK2 behavior rather than accuracy in turbulent closure statistics.

The new `verify_specific_diffusion.py` checks two grids, n=16 and 32, with varying density. Constant k gives exactly zero diffusion-on/off difference. With varying k, a directly assembled published face flux matches the native initial rate; maximum absolute rate errors are 2.34e-11 and 2.55e-11, within frozen O(dt) limits 5.34e-9 and 2.72e-9. Both lanes retain the same isotropic SGS pressure and Euler machinery so their difference isolates diffusion up to finite-step O(dt) effects. All eight lanes complete one unmodified RK2 step with no limiting or fallback.

The [four verification plots](fidelity-checks.png) show diffusion amplitudes and time convergence; they are deterministic equation checks, not turbulent benchmark comparisons. The native constant-acceleration oracle also passed, including the quadratic kinetic work and the 1/2,1/2 final RK2 forcing ledger.

## Reproduction and remaining claim

Build with `PROBLEM=sgs_verify`, double precision, CPU and MPI disabled, then use fresh output directories:

```sh
ctest --test-dir /tmp/sgs-build --output-on-failure
python3 -B tst/sgs/verify.py --exe /tmp/sgs-build/src/athena --output /tmp/sgs-core-new
python3 -B tst/sgs/verify_production.py --exe /tmp/sgs-build/src/athena --output /tmp/sgs-c1-new
python3 -B tst/sgs/verify_mixed_production.py --exe /tmp/sgs-build/src/athena --output /tmp/sgs-mixed-new
python3 -B tst/sgs/verify_heating.py --exe /tmp/sgs-build/src/athena --output /tmp/sgs-heating-new
python3 -B tst/sgs/verify_forcing.py --exe /tmp/sgs-build/src/athena --output /tmp/sgs-forcing-new
python3 -B tst/sgs/verify_forcing.py --oracle-only --exe /tmp/sgs-build/src/athena --output /tmp/sgs-force-oracle-new
python3 -B tst/sgs/verify_specific_diffusion.py --exe /tmp/sgs-build/src/athena --output /tmp/sgs-specific-diffusion-new
```

The mixed-production and specific-diffusion runners need NumPy. Plotting is not required for these checks. Detailed inputs, outputs, fingerprints, logs and machine-readable results are retained in `sgs-fidelity-audit-2026-10-13` under the current task's artifact directory.

The supported implementation passes this equation-fidelity audit against the cited authors' manuscript. No native correction is warranted by the evidence found here. Reproducing the paper's forced-box statistics remains a separate, unfinished claim; recent unforced decay disagreement is retained as additional validation evidence.
