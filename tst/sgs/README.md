# Schmidt Federrath SGS hydro verification on multiphase_turb

This is a restricted, experimentally qualified uniform-grid HD implementation.
The full M0–M9 scientific release gates are not complete. The wiki implementation
record lists executed instances, failed attempts and remaining coverage.

## Build and run

Use a unique build directory and output directory. CMake includes the independent
closure oracle when `BUILD_TESTING=ON` and double precision is enabled.

```sh
cmake -S . -B /tmp/sgs-build -DPROBLEM=sgs_verify \
  -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=OFF
cmake --build /tmp/sgs-build -j3
ctest --test-dir /tmp/sgs-build --output-on-failure
python3 -B tst/sgs/verify.py --exe /tmp/sgs-build/src/athena \
  --output /tmp/sgs-new-verification
python3 -B tst/sgs/test_configuration.py
python3 -B tst/sgs/configuration.py --exe /tmp/sgs-build/src/athena \
  --output /tmp/sgs-new-configuration
python3 -B tst/sgs/verify_wave.py --exe /tmp/sgs-build/src/athena \
  --output /tmp/sgs-new-wave
python3 -B tst/sgs/verify_forcing.py --exe /tmp/sgs-build/src/athena \
  --output /tmp/sgs-new-forcing
mkdir /tmp/sgs-example
cd /tmp/sgs-example
/tmp/sgs-build/src/athena -i /absolute/checkout/inputs/hydro/sgs_box.athinput
```

The configuration runner saves each attempted result before execution and after
completion. Timeouts and unreadable outputs are recorded as failures; later checks
continue. `test_configuration.py` exercises these failure paths without a build.

The box input is a short energy-evolving example, not a stationary turbulence
validation campaign. `verify_forcing.py` provides capped pilots and measured-work checks.
`filter_reference.py` and `test_filter_reference.py` require NumPy; their synthetic
tests validate the analysis machinery, not the physical closure.

## Model and state

The baseline is Schmidt & Federrath (2011), Eqs. (7)–(14):
<https://arxiv.org/html/1010.4492>. Write `A_ij = derivative_j u_i`,
`Q = sum(A_ij^2)`, and `S* = sym(A) - trace(A) I/3`.

```
tau = 2 C1 Delta sqrt(2 rho K) S*
      - 2 C2 K A A^T / Q - (2/3)(1-C2) K I
Sigma = tau:S
R = Ceps K sqrt(K/rho) / Delta
J = -Ckappa Delta sqrt(rho K) grad(K/rho)
```

Defaults are C1=0.02, C2=0.7, Ceps=1.5, Ckappa=0.65. K is an
energy density; primitive k=K/rho is a specific energy. Existing resolved E
continues to exclude K; its thermal pressure is unchanged. K is appended after
all true passive scalars, whose count, indices and output labels remain intact.
The represented energy is E+K. SGS diffusion transports represented energy.

`<sgs>` enables the model by default; `enabled=false` or omission preserves the
multiphase_turb layout and numerical path. `initialization=uniform` adds K=rho*seed once
on a fresh start (default seed=1e-6), without removing thermal energy. Exactly
zero uniform K stays zero; no repeated seed or floor activates it. Verification
generators can use `initialization=pgen` to supply a specified K profile.
Output selectors `hydro_u_sgs` and `hydro_w_sgs` give energy density and specific
energy. Existing `hydro_sgs` filtering moments retain their previous meaning.

## Numerical contract

The supported subset is double-precision Newtonian ideal-gas HD, uniform cubic
Cartesian cells, periodic active directions, RK2, LLF, DC or PLM, at least two
ghost zones, and CFL in (0,0.2]. One- and two-dimensional reductions retain all
three velocity components and the three-dimensional stress trace.

Face gradients use adjacent differences in the normal direction and averaged
centered transverse differences. One shared traction contributes momentum flux
`-tau.n` and work flux `-u_face.tau.n`, with arithmetic face velocity. Applied
production is `div(u_face.t)-u_cell.div(t)` from the final face traction.
RK2 pairs `-Sigma+R` in E with `Sigma-R` in K. The semidiscrete SGS thermal
source is R; the finite-step kinetic partition has a convergent RK error.

At exact Q=0, raw stress is `-2K I/3`; at K=0 it vanishes. Nonzero Q uses a
scaled norm, with no denominator epsilon. Its limit at zero gradient depends on
direction. Stencil differences below eight machine-epsilon times the magnitude
of the two velocities are treated as zero for closure gradients only. Actual
velocity differences remain in compatible work/production. Large boosts affect
conditioning; finite tested boosts do not prove exact frame invariance.

Negative face production and SGS diffusion have full-dt donor caps. The stress
cap scales the entire traction and work, so the *applied* trace need not equal
the raw -2K identity. History records suppressed stress transfer and diffusion.
Any inadmissible full-dt SSP predictor triggers a globally coordinated donor-cell
flux/stress/source rebuild. A second failure stops with an actionable error.
EOS floors cannot repair SGS states; SGS energy bypasses passive-scalar floors.
True passive scalars retain their existing flooring behavior.

Timestep estimates include thermal-plus-SGS regularization `cs^2=gamma P/rho+6k`,
neighbor SGS energy with receiving density, eddy/diffusion row sums, and explicit
decay. This frozen-shape regularization plus stage checks is not a general
nonlinear stability proof. An attempted k/|A| restriction stalled ordinary
contacts and was removed; the failed evidence remains recorded.

## Forcing, diagnostics and restart

SGS forcing requires `turb_driving/stage_consistent=true`. It uses old-stage
rho*a and rho*u.a RHS values; there is no quadratic energy kick at each RK stage.
Heun supplies the constant-acceleration quadratic kinetic work. Actual integrated
work and momentum sources are separate from production and dissipation. The
source thermal/cadence timestep caps can make injected work differ from nominal
`dedt`; use the measured history, not the requested injection setting.

Ordinary SGS restarts preserve K and all integrated counters. Mode, coefficients
and scalar layout changes are rejected before array loading. Forcing restarts retain the branch-native Fourier-mode OU amplitudes, RNG and
exact update count in `TurbulenceRestartMetadata` version 1. Every forcing parameter string and saved forcing-state
parameter is immutable on resume; even a numerically equivalent reformatted
override is rejected. SGS-off forcing retains the branch-native modal restart format.

History exposes K, E+K, thermal E, applied production/dissipation, rejected stress
and diffusion transfer, fallback count, and `force-work`, `force-mom1..3` when
stage forcing is active. Counters use Heun half-dt quadrature and global MPI
sums; history prevents multiplying global counters by the rank count.

AMR/SMR, MHD, relativistic and multiphysics modes, cooling, viscosity, conduction,
shearing, particles, frame tracking, user source/work hooks and legacy FOFC are rejected with
SGS enabled. Single precision and CUDA/HIP/SYCL SGS execution are unqualified and
rejected. Multi-rank MPI and GPU execution require separate qualification. Prior SGS
qualification on a different base does not qualify this port automatically.

The independent reader supports only single-file, serial, unrefined native-double
HD checkpoints. It checks the modal header and exact payload byte counts. Binary
field output on this branch can use `data_precision=real`; the budget and restart
checks use native checkpoints to include source counters, modal OU and RNG state.

## Evidence limits

The independent kernel oracle detects deliberate sign/factor mutations. Coupled
tests measure E+K, each energy reservoir, scalars, positivity, decomposition,
normal checkpoint restart, actual fallback, stress limiting and temporal decay
order. A phase-shifted 1D pressure oracle isolates the longitudinal C2 branch:
`b=2/3+4C2/3`, `K proportional to rho^(1+b)` and
`cs^2=gamma P/rho+b(1+b)k`; at C2=0.7 its SGS factor is 4.16.
The port runner checks sampled 16/32-cell acoustic refinement. Other pressure-wave cases
remain smoke tests. Exact full mixed manufactured
solutions, compression/expansion PDE oracles, active diffusion limiting, arbitrary
hard-flow stability, multi-rank/GPU performance, converged turbulent reference
accuracy and stationary campaign statistics remain open. Passing budgets is not
evidence of closure accuracy or scientific benefit.
