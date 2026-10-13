# Independent SGS forcing contract review

2026-10-08. Read-only review of existing `src/srcterms/turb_driver.cpp`,
task registration, restart writer/reader and proposed SGS source integration.
No forcing code edits, builds or simulation runs in this packet.

## Stage contract decision

The proposed stage RHS is correct for Newtonian resolved forcing:
momentum source=rho_old a; resolved-energy source=rho_old u_old dot a;
SGS source=0. Use unchanged stage primitives, include these same sources in
the full-dt SSP predictor, and apply with beta dt after SGS paired sources and
before ghost/primitive refresh. Skip the existing independent AddForcing task
only on the SGS path to prevent double application. Keep forcing preparation
before the cycle and the acceleration field fixed across its two RK stages.

Do not add a separate quadratic kick-energy term to either RK stage. For constant
rho and a, Heun gives v1=v0+a dt and vnext=v0+a dt, while its energy quadrature
is dt rho[v0 dot a+(v0+a dt) dot a]/2
=dt rho v0 dot a+rho |a|^2 dt^2/2. Thus E=e_int0+rho |v|^2/2,
with unchanged internal energy and K. The RK cross-stage contribution supplies
the required quadratic work naturally. This is the independent constant-force
oracle and a test that detects retaining the legacy quadratic stage term.

Integrate actual work using Heun weights dt/2 at each stage from the same
rho_old u_old dot a used by the RHS. Preserve accumulated work at full binary64
round-trip precision in restart metadata and include it in represented-energy
residuals. Input dedt is not a substitute for measured applied work.

The full-dt FE predictor has a temporary thermal deficit
-rho |a|^2 dt^2/2 for uniform forcing. It can correctly reject a strong/cold
forcing step even when the final analytic Heun state is admissible. Derive or
characterize a source restriction dt<=sqrt[2 e_margin/(rho |a|^2)] or clearly
restrict that regime; fallback reconstruction cannot repair source-only failure.

## Material inherited driver issues found

### Negative-correlation power normalization

`UpdateForcing` defines m0=(dt/2) integral rho|a_raw|^2 and
m1=integral momentum dot a_raw. Its constant-edot normalization should solve
m0 s^2+m1 s-dedt=0. For m1<0, current code uses
s=m1/(2m0)+sqrt[m1^2/(4m0^2)+dedt/m0], whereas the positive root is
s=(-m1+sqrt[m1^2+4m0 dedt])/(2m0).

Independent example: m0=0.5, m1=-1, dedt=0.1 gives current s approximately
0.0954451 and actual work rate approximately -0.0908902, not +0.1.
This is a pre-existing normalization bug, separate from SGS exchange. Fix or
restrict the affected SGS driver mode before claiming constant injected power;
retaining it with an accurate signed-work ledger still cannot imply dedt is met.
Set dedt==0 explicitly to zero forcing: merely choosing the mathematical positive
root can yield nonzero acceleration with zero net kick work when m1<0.
Validate dedt>=0, finite coefficients, positive dt_turb_update and valid mode
counts before numerical operations. The non-constant-edot mode has different
normalization weighting and should not be qualified by constant-edot evidence.

### OU restart does not restore actual process state

Current restart code serializes RNG state and the scaled `pturb->force` field.
It does not serialize `force_tmp1`, the raw OU memory used in
force_tmp1_next=fcorr force_tmp1_old+gcorr innovation. The constructor zeros
that memory. It infers n_turb_updates_yet=floor[(time-start)/dt_update], but
the completed loop stores n_turb_updates_reqd=floor[...] +1.

Consequently the first resumed InitializeModes call can consume an extra RNG
innovation and evolve from zero OU memory; UpdateForcing then overwrites the
restored scaled field. RNG restoration plus integrated work persistence alone
does not reproduce the forced trajectory. Preserve raw OU memory and the exact
update counter in a compatible marked restart layout, or reject SGS-forced
restart until that state lifecycle is implemented and tested. A pre-existing
baseline flaw must not be promoted into a qualified SGS restart claim.

## Driving assumptions and zero-force gates

- Acceleration is held fixed across both RK stages. Preparation is once per cycle;
  the OU update cadence is a separate physical-time discretization. Resolve that
  cadence and test crossings of tdriv_start, duration and update epochs. Do not
  claim deterministic temporal order of a stochastic path without controlling the
  driving history across timestep refinements.
- UpdateForcing currently suppresses forcing on cycle zero. Preserve/document or
  deliberately change this startup policy; tests must establish the actual start.
- Mean acceleration subtraction uses rho at cycle start. Stage-two density can
  change, so integral rho_stage2 a need not vanish. Measure actual accumulated
  momentum injection rather than infer zero mean momentum from initial weighting.
- Zero dedt, pre-start, post-duration for the finite-duration flag, and empty/zero
  mode-amplitude fields must yield finite zero forcing and zero accumulated work.
  No division by zero or artificial amplitude should appear in those cases.
- Spatially weighted/tiled, anisotropic, non-constant-edot and rank-changing
  forcing modes require their own coverage. A minimal periodic 3D untiled,
  unweighted baseline can be qualified first with explicit rejection elsewhere.

Required independent tests: constant deterministic force/no heat; zero force;
positive and negative normalization correlations; actual periodic represented-
energy minus applied-work residual; stage-one and stage-two cold-source trial;
one/many blocks; uninterrupted versus restart across an OU update; work and
momentum histories; resolved forcing cadence sensitivity. Passing an unforced
43-case suite does not execute these new paths or close their gates.

Conclusion: stage source integration is scientifically correct and ready to
implement incrementally. Existing normalization and OU state persistence need
explicit handling before forcing/persistence qualification is claimed.

## Implementation readback after bounded forcing packet

Reviewed new stage integration while implementation was landing. Current
`StageNormalization` uses a scaled hypot expression and the correct positive
quadratic root, with target==0 returning zero before the zero-norm rejection.
`StageForcing` uses unchanged w0, applies only linear RK RHS momentum/work,
leaves K unchanged and integrates work/momentum with half-dt weights. The
SGS full-dt trial includes those same forcing terms before admissibility checks.
Legacy forcing task registration is skipped when stage_consistent is enabled.

The marked restart now substitutes raw `force_tmp1` for scaled acceleration in
the existing three-component payload, preserving binary dimensions, and stores
exact `ou_updates`, `ou_state_version=1`, work and momentum counters. Read paths
restore into the raw memory conditionally on the same mode; main rejects layout
mode changes and old unmarked stage-forcing restart files. The constructor now
applies saved counter after the legacy time-derived inference, so it survives.
Raw OU field plus RNG plus exact counter are sufficient for the next innovation
and acceleration reconstruction at unchanged forcing configuration; retaining
Fourier innovation coefficients is unnecessary when no update is due, because
the already reconstructed raw OU memory is the carried process state.

**Bypass found and resolved:** an explicit stage_consistent=false option initially
could route legacy AddForcing on SGS, omitting the trial and work ledger. Final
source readback confirms TurbulenceDriver now rejects !stage_consistent whenever
phydro exists and isgs>=0, before task registration. The constructor emits
`SGS forcing requires turb_driving/stage_consistent=true`. Runtime rejection
evidence remains the verifier's responsibility; this is no longer an open source
finding.

The timestep cap initially occurred after power normalization. Lifecycle agent
agreed to cap cadence before evaluating m0 proportional to dt. A subsequent
thermal-admissibility dt cap can still change the actual target injection rate;
coordinator explicitly accepts this as a nominal dedt normalization and measures
actual work. Do not describe capped stages as exactly meeting requested dedt.
Fixed-force/no-heat oracle, signed-root cases, zero-power, marked OU restart and
uninterrupted-work tests are still required; no new simulation was independently
run in this five-minute read-only packet.

Final readback also confirms the cadence cap precedes m0 normalization, and
main snapshots every turb_driving parameter string before overrides, then rejects
changed values or additional parameters on marked stage-forcing restarts. This
strict policy includes OU cadence, spectrum/geometry, raw-state counter and work
metadata, avoiding silent reinterpretation of saved process state. Numerically
equivalent reformatted parameter values are intentionally rejected under this
exact-string policy; documentation should state that restriction.
