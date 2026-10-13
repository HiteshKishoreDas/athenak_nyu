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

## Equation fidelity audit

The branch also includes [all archived SGS check reports](../../docs/sgs/checks/README.md),
including the decay, spatial, forcing, CFL and spectral comparisons.


The [2026-10-13 audit](FIDELITY_AUDIT.md) maps the model equations to code,
records the tested RK2 energy exchange, and identifies numerical differences
from the original implementation. No native equation or coefficient changes
were needed. Agreement with the published forced-turbulence results remains
unestablished.

The additional density-varying diffusion check verifies that transport acts on
specific energy K/rho, including zero diffusion when K/rho is constant:

```sh
python3 -B tst/sgs/verify_specific_diffusion.py \
  --exe /tmp/sgs-build/src/athena --output /tmp/sgs-specific-diffusion-new
```

It requires NumPy and a fresh output directory. It compares diffusion-on/off
initial rates using the same Euler and isotropic SGS-pressure terms.

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
Cartesian cells, periodic active directions, RK2, LLF, DC, PLM or monotone PPM4.
DC/PLM need at least two ghost zones; PPM4 needs three. CFL in (0,0.2]. One- and two-dimensional reductions retain all
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

### 3D comparison

The same runner accepts `--ndim 3`. Laptop defaults become 64³ versus 16³,
with a 32³ reference sensitivity check. Preparation still drives to t=2,
followed by t=.5 unforced decay. All three directions have forcing modes,
and all diagnostics use complete volumes. Field maps show central z-cell
slices: native and coarse cell centers differ; the filtered reference and
both coarse lanes share the same slice. Vorticity errors compare vectors.

Rebuild the executable after updating the test initial-state loader, then run:

```sh
cmake --build build-sgs-laptop -j2
/tmp/sgs-plot-env/bin/python -B tst/sgs/compare_resolution.py --self-test \
  --exe build-sgs-laptop/src/athena
/tmp/sgs-plot-env/bin/python -B tst/sgs/compare_resolution.py --ndim 3 \
  --timeout 900 --exe build-sgs-laptop/src/athena \
  --output /tmp/sgs-3d-resolution-comparison
open /tmp/sgs-3d-resolution-comparison/index.html
```

Pilot with `--spinup .1 --decay .01` in a separate directory. 3D is capped at
64³; the explicit `--timeout` bounds each simulation (at most 3600 seconds).
Existing verification runners retain their original default caps and numerical
thresholds. Timeouts and incomplete runs stay in the result ledger.

The report adds initial/final stress and production checks with exact filtered K
supplied, using the fixed published coefficients. `resolution_closure.py` checks
aligned density-weighted covariance and physical derivative axes. Its coarse,
cell-centered tau:gradient diagnostic tests the mathematical closure; it does
not reproduce the native face-centered production discretization. All statistical
metrics and spectra use full 3D volumes. Ten plot pages contain 47 panels.

The test generator accepts `nx ny nz nvar` initial-state headers with z/y/x rows
for 3D; the 2D `nx ny nvar` format remains supported. No solver changes.
A 64³ reference with a four-to-one filter ratio is not established DNS. These
single-seed decay runs cannot attribute every difference to dimensionality or
establish statistical convergence. To compare the same cell counts in 2D,
run `--ndim 2 --high 64 --low 16` with the same preparation/decay times.

### Larger 2D comparison

The existing comparison also accepts 256²/512² reference grids and 64²/128²
coarse grids. Require at least four reference cells per coarse cell per direction.
The 512²/64² case uses an eight-to-one filter ratio; its 256² reference check
uses four-to-one. Coefficients, solver, CFL, seed, forcing and preparation/decay
intervals remain the same as the smaller comparisons.

```sh
/tmp/sgs-plot-env/bin/python -B tst/sgs/compare_resolution.py --ndim 2 \
  --high 512 --low 64 --timeout 3600 --exe build-sgs-laptop/src/athena \
  --output /tmp/sgs-2d-large
open /tmp/sgs-2d-large/index.html
```

Pilot in a separate folder with `--spinup .1 --decay .01` first. The per-case
cycle cap scales from 20000 to 40000 at 256² and 80000 at 512², and is recorded
in the manifest and saved inputs. The optional timeout is capped at 3600 seconds;
defaults and all numerical acceptance thresholds remain unchanged. No rebuild
is required for this test-runner-only extension. These are resolution sensitivity
checks, not converged DNS or an isolated test of the source of numerical losses.

### Driven 3D paper test

`paper_test.py` attempts the solenoidal Section 4/Table 4 experiment of
[Schmidt & Federrath (2011)](https://arxiv.org/html/1010.4492): uniform periodic
64³ box, gamma=1.001, continuous OU driving, fixed published SGS coefficients,
box side 2L=1, V=1 and T=L/V=0.5. The native forcing uses signed modes in all
three directions, n=1–3 with peak 2. Only forcing amplitude is calibrated;
closure coefficients stay fixed. A tiny positive k0=1e-8 avoids an absorbing
zero-SGS state. The assumed V/c0=5.5 is explicit because Table 4 does not state
this ratio for its pure-solenoidal case.

Pilot before the 10T run. A 32³ pilot estimates cost/stability, but lies below
the paper's lowest meaningful resolution. Use fresh output directories:

```sh
python3 -B tst/sgs/verify_signed_forcing.py --exe /tmp/sgs-build/src/athena \
  --output /tmp/sgs-signed-check
/tmp/sgs-plot-env/bin/python -B tst/sgs/paper_test.py --self-test
/tmp/sgs-plot-env/bin/python -B tst/sgs/paper_plots.py --self-test
/tmp/sgs-plot-env/bin/python -B tst/sgs/paper_test.py \
  --exe /tmp/sgs-build/src/athena --output /tmp/sgs-paper-pilot \
  --n 32 --end-T 4 --accel-rms 2 --timeout 600
/tmp/sgs-plot-env/bin/python -B tst/sgs/paper_test.py \
  --exe /tmp/sgs-build/src/athena --output /tmp/sgs-paper-64 \
  --n 64 --end-T 10 --accel-rms 3 --timeout 7200
open /tmp/sgs-paper-64/report.pdf
```

The full run has a 100000-cycle cap and a 7200-second simulation cap. Set
`OMP_NUM_THREADS=8 OMP_PROC_BIND=false` with an OpenMP build; a serial run
may cost substantially more. The independent native reader remains rank 1,
double precision, unrefined HD. NumPy and Matplotlib are needed for plots.
Checkpoint validation streams through all snapshots; only first/final states
are retained, keeping the 81-window-snapshot case bounded in memory.

The report compares volume RMS resolved/SGS Mach, normalized mean K, and both
volume-mean specific dissipation and density-weighted dissipation. Rates from
applied counters remain separate from checkpoint closure estimates. Statistics
use actual timestamps and time-weighted averages over 2T–10T, temporal scatter,
eight 1T block means, early/late means and trends. Scatter is not an uncertainty
of the mean. Completing the window does not establish stationarity.

Every input, raw checkpoint, history, log, result, protocol and source/executable
fingerprint is retained. `diagnostics.json` stores compact cell-derived statistics;
`statistics.json` records the definitions, sampling window and numeric results.
The five-page 18-panel PDF also has `page1.png`…`page5.png` previews. Regenerate:

```sh
/tmp/sgs-plot-env/bin/python -B tst/sgs/paper_test.py --plot-only \
  --output /tmp/sgs-paper-64
```

This is an approximate reproduction: AthenaK RK2/LLF/PLM differs from Enzo 1.5,
instantaneous acceleration RMS normalization differs from modal OU amplitude,
and equivalence of the native shell-compensated spectral envelope is unproven.
These affect comparisons in V units even after Mach calibration. Earlier frozen paper reports lack an applied resolved pressure-work ledger, so
their compression-corrected heating fraction is unavailable. New builds provide
the ledger and separate heating report described below. No high-resolution filtered CDF or resolution-scaling claim
comes from this single 64³ case. Table 4 is a reference, not a tuned pass threshold.

## Controlled reconstruction comparison

`paper_test.py --reconstruct ppm4` uses the native monotone PPM reconstruction
with LLF and the same SGS flux/source machinery. It requires three ghost zones;
PLM defaults and coefficients are unchanged. PPMX, WENOZ and HLLC remain rejected.
Run the bounded deterministic check before turbulent comparisons:

```sh
python3 -B tst/sgs/verify_ppm.py --exe build-sgs-laptop/src/athena \
  --output /tmp/sgs-ppm-check
python3 -B tst/sgs/paper_test.py --exe build-sgs-laptop/src/athena \
  --output /tmp/sgs-paper-ppm64 --n 64 --end-T 4 --accel-rms 3 \
  --reconstruct ppm4 --timeout 1800
```

The 4T run is a cost/admissibility pilot with a short 2T–4T averaging window.
A 10T run is needed for the full paper sampling interval. A single seed does
not measure ensemble uncertainty. Reconstruction comparisons hold forcing and
closure coefficients fixed; they do not calibrate the paper's forcing units.

## Native production amplitude check

```sh
python3 -B tst/sgs/verify_production.py --exe build-sgs-laptop/src/athena \
  --output /tmp/sgs-production-check
```

Eight one-step 16³/32³ cyclic-shear cases compare the actual applied production
and dissipation counters with independent closed Fourier expressions. Both PLM
and PPM4 retain the default mixed closure coefficients. Density and specific SGS
energy vary between cases; three ghost zones are used for both reconstructions.
The rate tolerance is frozen at 2e-5 before execution, allowing O(dt) stage
changes over dt=1e-7. Exactly one RK2 step and zero limiting/fallback are required.
This closes the gap between kernel formula checks and native source amplitudes.
It does not test general C2 production, transverse-gradient accuracy or diffusion:
initial A A^T is diagonal, divergence is zero and K/rho is uniform. The separate
longitudinal mixed wave check covers an isolated C2 acoustic branch.

## Compression-corrected heating budget

The SGS history now records `sgs-pressure-work-integral`, the cumulative signed
thermal pressure work from the actual reconstructed LLF momentum flux. The SGS
stress, source exchanges and fluid update are unchanged. The history/reduction
capacity grows from 32 to 33 to preserve the forced twelve-scalar case.

```sh
OMP_NUM_THREADS=4 OMP_PROC_BIND=false python3 -B tst/sgs/verify_heating.py \
  --exe build-sgs-laptop/src/athena --output /tmp/sgs-heating-controls
/tmp/sgs-plot-env/bin/python -B tst/sgs/heating_plots.py --self-test
/tmp/sgs-plot-env/bin/python -B tst/sgs/heating_plots.py \
  --runs /tmp/sgs-paper-plm /tmp/sgs-paper-ppm \
  --output /tmp/sgs-heating-report --start-T 2 --end-T 4
```

Build the existing `sgs_verify` problem. `verify_heating.py` tests independent
compression/expansion and DC-shear residual amplitudes, temporal convergence and
restart. Add `--baseline-exe` to check exact field passivity against the previous
binary. Fractions are ratios of integrated increments over the selected window:
SGS share = D/(delta U-C), numerical share = (delta U-C-D)/(delta U-C).
C is resolved thermal-pressure work and D is applied SGS dissipation. The residual
includes spatial and RK effects and shock-associated heating; it is not LLF-only.
Signed residuals are retained. Nonpositive corrected-heating intervals have no
reported fraction. The report uses actual native checkpoint configuration and
linear endpoint clipping of the history, and does not certify stationarity.

`paper_test.py --cfl .05` enables a timestep sensitivity lane. Test continuations
may explicitly override `time/cfl_number`; other physics restart guards remain.
Older checkpoints start the new pressure counter at zero at restart time. Use
fresh output directories, matching counter/inventory increments, and an averaging
window covered by the new ledger. Do not append the new schema to an old history.

Restart tests save the exact checkpoint inventory and counter origin in
`history-origin.json`, covering the gap before the first scheduled history output.
The heating report includes that point when the requested window needs it.

## Full 2T–10T statistics and restart segments

`paper_test.py` accepts `--seed`, pure `--sol-fraction 0` (compressive) or `1`
(solenoidal), grids through 256³, and an explicit `--nlim` cycle cap. The default
physics and existing numerical acceptance limits remain unchanged. Each run
records its protocol before execution; a timeout or cycle-cap stop is a failure.
Larger grids require substantially more time, memory and checkpoint storage.
The local 128³ startup pilot took 528 seconds for only 0.2T on two threads;
startup extrapolation already exceeds seven hours for 10T and understates the
cost of developed turbulence. A 256³ full checkpoint is about 2.79 GB.

For a fresh 64³ run, use the existing plotting environment and serial-double
`sgs_verify` build. These arguments preserve the audited forcing rather than
claiming its equivalence to the paper:

```sh
OMP_NUM_THREADS=4 OMP_PROC_BIND=false /tmp/sgs-plot-env/bin/python -B \
  tst/sgs/paper_test.py --exe /tmp/sgs-build/src/athena \
  --output /tmp/sgs-paper-64-seed43 --n 64 --seed 43 \
  --reconstruct ppm4 --nghost 3 --accel-rms 3 --end-T 10 --timeout 9000
```

`stitch_report.py` is confined to joining one verified diagnostic prefix and one
restart continuation. It requires the prefix's successful numerical record,
its exact field/forcing passivity record, and a successful continuation. It checks
the native checkpoint inventories, cumulative counters, diagnostic boundary,
and stable physics parameters. Counters are never reset or rebased. Original
files are untouched. The output includes combined diagnostics, five PDF/PNG
pages, time-weighted statistics, eight 1T blocks, early/late comparisons,
stationarity screening against a supplied frozen JSON protocol, and hashes.

```sh
/tmp/sgs-plot-env/bin/python -B tst/sgs/test_stitch_report.py
/tmp/sgs-plot-env/bin/python -B tst/sgs/stitch_report.py \
  --prefix-diagnostics /path/to/verified-prefix/diagnostics.json \
  --prefix-run /path/to/pressure-ledger-prefix \
  --continuation /path/to/completed-continuation \
  --checkpoint /path/to/pinned-4T.rst \
  --screen /path/to/stationarity-protocol.json --output /tmp/sgs-full-report
/tmp/sgs-plot-env/bin/python -B tst/sgs/heating_plots.py \
  --runs /tmp/sgs-full-report --start-T 2 --end-T 10 \
  --output /tmp/sgs-full-heating
```

Full window coverage does not imply stationarity or paper reproduction. Temporal
standard deviations and adjacent 1T blocks are descriptive; they are not errors
on the mean. Independent seeds, 64³/128³/256³ resolution statistics and an
adequately resolved Gaussian-filtered reference remain separate requirements.
Native instantaneous forcing normalization and radial weighting differ from
published modal forcing descriptions. Original 2011 input provenance and the
Gaussian width convention have not been recovered; do not tune amplitudes or
SGS coefficients to disguise those gaps.

## Native nonlinear/mixed production amplitude

Reuse the `sgs_verify` binary and installed NumPy environment:

```sh
OMP_NUM_THREADS=2 python3 -B tst/sgs/verify_mixed_production.py --self-test
OMP_NUM_THREADS=2 python3 -B tst/sgs/verify_mixed_production.py \
  --exe build-sgs-laptop/src/athena --output /tmp/sgs-mixed-production
```

The 16³/32³ PLM/PPM4 checks use the existing 3D initial-state loader. Four full
velocity-gradient tensors exercise compression, expansion, trace-free strain
and both signs of C2 transfer. Affine periodic seams are excluded from the local
oracle; only the central core is compared with the analytic equation. Identical
initial states with C2=0/.7 isolate its deviatoric production. Full local K rates
also include conservative advection, C1 and dissipation at the published defaults.

A separate smooth periodic rank-one gradient with modulated K has an exact
shared-face global production sum, including all three velocity-work components.
It checks native counters, independently read K inventories, paired resolved
energy and kinetic stress work. This term-isolation fixture sets C1, dissipation
and diffusion to zero; it does not change their defaults or the turbulent runs.
Every case must complete one RK2 step dt=1e-7 without limiting or fallback.
The rate tolerance 2e-5 is frozen before execution, matching the earlier C1 test.
Results, protocols, all inputs and failed attempts are retained. Passing checks
native amplitudes and exchange, not filtered turbulent closure accuracy.


### Paper Section 3 filtered-snapshot test

`verify_filtered_snapshot.py` applies Gaussian and box Favre filters to supplied
periodic unit-cube SGS-off native 3D checkpoints. It computes exact filtered
stress, covariance energy and signed production, then supplies that exact energy
to the actual C++ `sgs::Stress` kernel through `sgs_filtered_stress`. Coefficients
remain C1=.02/C2=.7. No forcing parameter enters the snapshot comparison.
This tests the production closure, not SGS reservoir evolution or face fluxes.

Build only the test utility with an existing double-precision CPU build:

```sh
cmake -S . -B /tmp/sgs-build -DBUILD_TESTING=ON
cmake --build /tmp/sgs-build --target sgs_filtered_stress -j2
python3 -B tst/sgs/test_filter_reference.py
python3 -B tst/sgs/verify_filtered_snapshot.py --self-test \
  --kernel /tmp/sgs-build/src/sgs_filtered_stress
python3 -B tst/sgs/verify_filtered_snapshot.py \
  --kernel /tmp/sgs-build/src/sgs_filtered_stress \
  --snapshots /absolute/reference.00008.rst /absolute/reference.00040.rst \
  --output /tmp/sgs-new-snapshot-comparison
```

Use the NumPy/Matplotlib environment above instead of `python3` if necessary.
The output directory must be new. The report contains JSON metrics, source/input
hashes and fourteen plot panels (PNG/PDF). All metrics cover full volumes; displayed
field maps are central-z slices. Gaussian width uses Delta_G=sqrt(12)*sigma,
an explicit second-moment convention, not a certified paper-identical filter.
Physical widths are box/16, box/8 and box/4; sampling remains on the fine grid.
A 64³ input thus has Gaussian widths 4/8/16 fine cells. The smallest width at 32³
has only 2 fine cells and is explicitly flagged as undersampled. Do not interpret
its exceptionally high correlation as evidence of reference convergence.

Checks cover independent direct centered covariance, mass/momentum preservation,
kinetic-energy decomposition, covariance positivity, native/Python stress,
trace/symmetry, Galilean invariance and spectral-derivative sensitivity. Reported
physical metrics include full/trace-free stress and production errors, signed
correlation, mean amplitude bias, sign errors, forward/backscatter and density
subsets. Arithmetic checks have frozen tolerances; physical bias is reported
without a tuned pass threshold. Neither limited-grid inputs nor correlation
alone establish a resolved reference or evolved SGS-energy accuracy.


### Unforced supersonic decay: evolved SGS energy

`decay_turbulence.py` runs six matched 3D lanes without a forcing block: 128³ and
64³ SGS-off numerical references, plus 8³/16³ SGS-on/off comparisons. Initial
rho=1 and cs=1; fixed-seed Fourier modes1–4 contain 2/3 solenoidal kinetic energy,
with initial RMS Mach9. Conservative block averages initialize all coarse lanes.
Removed kinetic variance initializes SGS K on, or becomes heat off; represented
energy is identical. Exact initial K agreement is imposed, not a prediction.

This follows the filtered-reference decay method of [Semenov2025](https://arxiv.org/pdf/2410.23339v3),
with declared departures: local undeveloped initial flow instead of imported
Kritsuk data, gamma1.001 ideal-gas approximation instead of exact isothermality,
smaller references, and AthenaK RK2/LLF/PPM4. SGS coefficients remain unchanged.
T=box/(2*initial velocity RMS)=1/18. Assess1–10T separately from the startup.

Use the existing NumPy/Matplotlib environment and native `sgs_verify` binary:

```sh
python3 -B tst/sgs/decay_turbulence.py --self-test
OMP_NUM_THREADS=2 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build-sgs-laptop/src/athena --output /tmp/sgs-decay-new \
  --high 128 --low 8 16 --end-T 10 --timeout 3600
python3 -B tst/sgs/decay_turbulence.py --output /tmp/sgs-decay-new --plot-only
```

The output directory must be new. For a short execution pilot use `--high 32`, `--low 4 8 --end-T 2 --timeout 120`. Largest supported
reference256³; use `--block 32` to reduce ghost-cell memory and pilot cost/memory
before increasing it. `--block` defaults to16 and also accepts32 or64; each
coarse lane uses at most its own grid size. `--nscalars 0` omits the two
ordinary passive tracers to reduce memory and runtime; the default remains2.
SGS is still a separate evolved variable. A pilot should establish physical-field
equivalence before using execution optimizations in a new campaign. Runs keep inputs, logs,
native-double binary fields, native checkpoints when written, histories,
source/input hashes, checks and failure records. Five PDF pages and PNG figures
show energy/ratios, total kinetic decay, sound-speed drift, conservation, applied
source/heating budgets, density-conditioned SGS velocities, spectra and slices.
Native binary time headers have six significant digits; final completion uses
full-precision history time, with a separate binary-time consistency check.

`summary.json` reports time-integrated K ratios and reference sensitivity at the
same physical filter. The predeclared10% maximum reference-sensitivity screen
covers1T to end; it is a resolution screen, not a model-accuracy acceptance test.
A failed screen prevents a firm quantitative SGS-energy verdict. No retuning is
performed. The8³ scale has better reference separation but poorly resolves the
energy-containing modes; a single seed and a passing two-grid screen do not
establish universal closure accuracy or a proven DNS. Heating residuals include
numerical/spatial/RK effects after subtracting thermal pressure work and SGS heat.

For the measured24GB laptop setup, the256³ campaign can use64-cell blocks and
omit unused passive tracers. It writes checkpoints every1T for recovery. The
scientific comparison retains the same IC seed, EOS, CFL, solver and SGS coefficients:

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build/src/athena --output /tmp/sgs-decay256-new \
  --high 256 --low 8 16 --block 64 --nscalars 0 --end-T 10 --timeout 7200
```

A64³ pilot confirms bitwise identical physical fields with/without the ordinary
passives and across16/32/64-cell blocks. This is an execution optimization;
reference convergence remains a measured scientific outcome.


To compare that completed256³ campaign with64³/128³ SGS, reuse its native
references instead of rerunning256³. Use matching seed, EOS, CFL, block, thread,
passive-scalar and time parameters. Reuse checks protocol, native source/executable,
inputs, initial-state and field hashes, and history equality; it recomputes Favre
filters and numerical checks. The matching128³ SGS-off reference also serves as
the128³ control. Original output directories must remain intact.

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build/src/athena --output /tmp/sgs-decay64-128-new \
  --reference-output /tmp/sgs-decay256-new --high 256 --low 64 128 \
  --block 64 --nscalars 0 --end-T 10 --timeout 7200
```

The256³ reference has only4/2 cells per64³/128³ coarse cell. At128³ the
128³ SGS-off control cannot independently qualify unresolved K at its own grid
scale: the second-reference convergence screen is unavailable, not passed or
failed. A finer reference pair is needed for accuracy qualification.

A small native check exercises reuse and that same-grid reporting case:

```sh
OMP_NUM_THREADS=2 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build/src/athena --output /tmp/sgs-reuse-check-reference \
  --high 32 --low 4 8 --block 16 --nscalars 0 --end-T 1.1 --timeout 120
OMP_NUM_THREADS=2 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build/src/athena --output /tmp/sgs-reuse-check-coarse \
  --reference-output /tmp/sgs-reuse-check-reference --high 32 --low 8 16 \
  --block 16 --nscalars 0 --end-T 1.1 --timeout 120
python3 - <<'PY_CHECK'
import json
from pathlib import Path
s = json.loads(Path('/tmp/sgs-reuse-check-coarse/summary.json').read_text())
assert s['failed_checks'] == 0
assert sum(r['reused_from'] is not None for r in s['runs'].values()) == 3
assert s['scores']['16']['reference_screen_pass'] is None
PY_CHECK
```
