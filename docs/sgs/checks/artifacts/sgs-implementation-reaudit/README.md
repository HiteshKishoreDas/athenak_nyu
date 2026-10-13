# SGS implementation re-audit

Audited production source: `7cdff3618726c3ac386a8bb1efa7ec32f0b7d6a2` on
`codex/multiphase-schmidt-federrath-hd`. Test-only local commit: `3ca288aa`.
Three independent read-only reviews covered equations, solver integration and
paper diagnostics. The coordinator added a bounded native production oracle;
the equation reviewer independently checked its derivation, implementation and
saved results. No production equation or coefficient was changed. No push/merge.

## Findings

No confirmed algebra, sign, missing-half, variable-slot or stage-integration defect
was found in the supported CPU/double/uniform-periodic HD path. This audit does
not establish accuracy of the turbulent SGS closure or paper reproduction.

- `src/hydro/sgs_closure.hpp`: K is energy density, k=K/rho; the nonlinear stress
  coefficient includes the gradient-norm factor correctly. Stress trace is -2K;
  pressure is 2K/3; dissipation is Ceps K^(3/2)/(sqrt(rho) Delta); diffusion acts
  on K/rho. These match [the published equations 7–14](https://arxiv.org/html/1010.4492).
- `src/hydro/hydro_sgs.cpp:215,260`: each adjacent half-face work contribution
  gives the correct total production for a constant gradient. Paired sources
  add P-D to K and subtract it from resolved E, which excludes K.
- `src/driver/driver.cpp:114`, `src/hydro/hydro_update.cpp:79`: RK2 and the source
  stage blend produce dt/2 times each stage rate. Integrated diagnostics use
  the final-solution weights. These weights were verified by code inspection;
  the new almost-constant one-step oracle does not independently resolve them.
- `src/eos/ideal_hyd.cpp:107`, `src/hydro/hydro_fluxes.cpp:134`: primitive k is
  K/rho and advects with mass flux. Its slot follows actual passive scalars;
  ordinary scalar flooring excludes it. Allocation, communication and restart
  carry it through nvars. No aliasing or restart loss was found.
- `tst/sgs/paper_test.py:30`: native checkpoint K, cellwise k, Mach, specific
  dissipation and density-weighted dissipation are correctly distinguished.
  Independent actual-checkpoint and clipped-time-average arithmetic reproduce
  saved diagnostics. Density correlations matter; matching SGS Mach does not
  require matching volume-mean K.

## Native amplitude test added and executed

Earlier global K-change versus its own production ledger checks verified
bookkeeping but did not independently establish native production amplitude.
`verify_production.py` closes a bounded part of that gap without changing solver
machinery. The eight native 16^3/32^3 cases use PLM/PPM4, three ghost zones,
(rho0,k0)=(1,.01) or (3,.08), all default mixed coefficients and exactly one
RK2 step dt=1e-7. Timeout is 30 seconds/case. Threshold 2e-5 was written before
execution and never changed.

Uniform rho and k, and velocity (.3 sin(2pi y), .3 sin(2pi z), .3 sin(2pi x)),
have divergence zero and diagonal A A^T. The closed initial volume rates are

```
dx = 1/N
s = 2*.3*sin(pi/N)/dx
P0 = 3*C1*dx*sqrt(2*rho0^2*k0)*s^2/2
D0 = Ceps*rho0*k0^(3/2)/dx
```

Expectations use these independent Fourier/decay expressions, not solver
stress or gradient helpers. Actual rates are changes in applied counters divided
by actual elapsed time. The small O(dt) difference is the evolved RK2 predictor.

| Executed check | Result |
|---|---|
| Native amplitude/structure checks | 24/24 PASS across eight cases |
| Maximum production relative error | 7.3938233811e-6 (0.0007394%) |
| Maximum dissipation relative error | 1.0169882767e-6 (0.0001017%) |
| Limiting and fallback | Exactly zero in every case |
| Independent C++ closure oracle | 1025 cases, 22575 assertions PASS |
| Maximum kernel scaled error | 1.4521283902e-15; threshold 5e-13 |
| Python lint/parse and staged whitespace | PASS |

The amplitude oracle covers initial C1 and dissipation, rho/K conversion and
aggregate counter normalization. It does not independently verify individual
RK stage weights, general C2 production, transverse-gradient accuracy, variable
coefficient interpolation or diffusion. Cyclic symmetry can hide some axis
permutations. Existing isolated longitudinal C2 wave evidence complements it;
no new general mixed manufactured PDE or filtered turbulent accuracy claim follows.

## Numerical and benchmark limitations

LLF (`src/hydro/rsolvers/llf_hyd_singlestate.hpp:60`) uses sqrt(cs^2+6k) for
numerical damping. Resolved kinetic energy lost numerically can become thermal
energy without feeding K; only explicit stress production feeds K. Increasing K
increases that LLF speed. This is a verified mechanism, not a measurement of its
share of the driven-box deficit. Existing fixed-forcing PPM4 versus PLM results
(K +15.43%, density-weighted dissipation +24.74% at 64^3) support method
sensitivity but do not isolate LLF or resolve the full discrepancy.

Normal gradients use adjacent differences; transverse gradients use averaged
centered differences (`hydro_sgs.cpp:200`). These attenuate near-cutoff modes.
The independent normal Fourier squared-gradient ratio is
[sin(q/2)/(q/2)]^2: .810569 at mode16/64, .509204 at mode28/64.
This is wavelength-dependent, not a uniform factor-of-two error. Global DC
fallback can add damping when active, but fallback and diffusion limiting were
absent in the completed comparison lanes. The negative-production limiter raises
applied global production at fixed state, so it cannot directly explain low K.

The paper comparison remains approximate. The runner assumes V/c0=5.5;
Table4 does not explicitly provide the absolute forcing amplitude or that ratio.
The actual OU correlation T=.5 and L=.5 fix V=L/T=1. Statistics cover 2T–4T,
not the paper's 2T–10T. Mach matching does not prove forcing equivalence.
[The prior primary-source audit](/Users/hitesh/Documents/das-wiki/science/athenak-subgrid-turbulence/evidence/2026-10-09-scheme-comparison/sgs-forcing-convention-audit.md)
remains valid. No existing results were renormalized or coefficients retuned.
The table's epsilon header versus rho-epsilon prose ambiguity remains recorded.

The next useful diagnostic is the measured resolved numerical-loss contribution
in the driven box, with pressure work and SGS transfer distinguished. A matched
forcing/time convention must precede a claim of paper-level quantitative agreement.

## Evidence and reproduction

Run [local-reproduce.sh](local-reproduce.sh) on this Mac using the existing build,
or build the SGS verification problem following the source README. Choose a new
output directory; the harness intentionally refuses existing output directories.
Only standard-library Python is needed for this test.

[Production protocol](native-production/production-protocol.json), [amplitude results](native-production/production-results.json),
[existing budget checks](native-production/results.json), [run manifest](native-production/manifest.json),
[provenance](provenance.json), and [kernel stdout](closure-oracle.log) retain
actual values and hashes. Raw native checkpoints, inputs, histories and logs stay
outside wiki Git at the provenance raw_evidence path. The run manifest records
7cdff361 plus the then-uncommitted test/docs patch; it is not relabelled as having
run after its future commit. Production source and harness hashes match 3ca288aa.
The previous failed runs, recovered campaign and frozen archive were untouched.
