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

## Laptop 2D visual report

`laptop_2d.py` runs 22 small serial 2D simulations on 16² and 32² grids and
checks 25 result records. It reuses the problem generator, native checkpoint
reader and acceptance thresholds above. It leaves every input, checkpoint,
history, log and result on disk, plus 15 PNG figures (92 panels), a PDF and
an HTML gallery. The C++ solver is unchanged.

From this SGS branch's checkout, build and run locally on a Mac:

```sh
cmake -S . -B build-sgs-laptop -DPROBLEM=sgs_verify \
  -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=OFF
cmake --build build-sgs-laptop -j2
python3 -m venv /tmp/sgs-plot-env
/tmp/sgs-plot-env/bin/python -m pip install numpy matplotlib
/tmp/sgs-plot-env/bin/python -B tst/sgs/laptop_2d.py --self-test
/tmp/sgs-plot-env/bin/python -B tst/sgs/laptop_2d.py \
  --exe build-sgs-laptop/src/athena --output /tmp/sgs-2d-report
open /tmp/sgs-2d-report/index.html
```

Use a new output directory for each simulation run. To regenerate plots from
saved results without rerunning the simulation:

```sh
/tmp/sgs-plot-env/bin/python -B tst/sgs/laptop_2d.py \
  --plot-only --output /tmp/sgs-2d-report
```

Inspect `results.json` for thresholds and measured residuals, `manifest.json`
for source/executable fingerprints, and each case's `input.athinput` to run or
change one case yourself. If running the executable directly, use a fresh working
directory so its output folders do not collide with saved results.

The plots cover analytic uniform SGS decay and heat gain; RK2 time refinement;
2D shear, diffusion and moving-contact profiles; box evolution and SGS-on/off
maps; independent mass/energy/scalar budgets; zero-K, restart and block-layout
checks; solenoidal/compressive driving; density-weighted shell spectra; solver
elapsed time, throughput and admissibility. The standalone check validates image
orientation, pressure conversion, periodic vorticity and FFT energy normalization.

The contact oracle disables dissipation and SGS diffusion to isolate advection.
Grid refinement changes the SGS filter width in the full-model cases, so those
resolution comparisons are sensitivity checks. Driven boxes are t=1 pilots with
three velocity components and kz=0, not stationary or reference-validated
2D turbulence. Runtime bars use AthenaK's elapsed Kokkos timer (labelled `cpu time used`
in the original log), with output included. They are single short local measurements;
they are not a statistically controlled performance benchmark.

## High-resolution turbulence versus low-resolution SGS

`compare_resolution.py` first drives a 128² hydro flow (SGS off) for t=2 with
seed 42, mixed solenoidal/compressive forcing in modes 1–3 and dedt=1. It then
compares unforced decay for another t=0.5 in four lanes: 128² reference, 64²
reference sensitivity check, 32² SGS off and 32² SGS on. This is a developed
2D flow with three velocity components; it is not a validated DNS or a 3D
stationary-turbulence qualification.

All lanes start from the same high-resolution snapshot. Coarse density,
momentum, total energy and scalar masses are aligned finite-volume averages.
The density-weighted velocity variance removed by filtering supplies initial K.
SGS-on subtracts that K from resolved E and stores it separately, preserving
thermal energy and represented E+K. SGS-off retains the conservative average E,
so the removed kinetic variance initially becomes heat. Both lanes have exactly
the same initial filtered density, momentum and represented energy; their initial
thermal pressures differ by the stated energy partition.

Reuse the build and Python environment above:

```sh
/tmp/sgs-plot-env/bin/python -B tst/sgs/compare_resolution.py --self-test \
  --exe build-sgs-laptop/src/athena
/tmp/sgs-plot-env/bin/python -B tst/sgs/compare_resolution.py \
  --exe build-sgs-laptop/src/athena --output /tmp/sgs-resolution-comparison
open /tmp/sgs-resolution-comparison/index.html
```

The self-check covers a hand-calculated Favre variance, energy partition,
alignment, boost invariance, native file round-trip and eight malformed inputs.
The full run keeps every IC, checkpoint, history and log, plus nine plotted PDF
pages, PNGs, HTML and `comparison.json`. To regenerate plots, use `--plot-only`
with the same output directory. The matched-initial-state checks compare file
readback exactly and bound global represented-energy differences at 1e-12.

For a smaller smoke run, use `--high 64 --low 16 --spinup 0.1 --decay 0.01`.
The comparison explicitly raises its grid cap; the original qualification runner
retains its default 32-cell cap and unchanged numerical thresholds.

Final errors compare all lanes at the exact same stop time. Intermediate error
curves linearly interpolate conserved checkpoints to common times, with maximum
brackets recorded. Every high-resolution reference field is block-filtered onto
the low grid before comparing density, Favre velocity, thermal pressure, kinetic
variance and spectra. The 64² versus 128² sensitivity is plotted alongside the
coarse errors. Improvement from SGS is a measured outcome, never a pass criterion.
One seed and two reference resolutions do not prove converged closure accuracy.

The test generator's optional `problem/initial_state` reads a text file containing
`nx ny nvar`, followed by rows in y/x order: rho, three momenta, resolved E, true
scalar densities, then K. E excludes K. SGS-off files must give K=0. Shapes,
finite values, positive density/thermal energy and nonnegative K are checked.
This input path is confined to `sgs_verify.cpp`; solver/restart machinery is unchanged.

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
