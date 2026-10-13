# Quick decaying supersonic turbulence pilot

Source099109b6; core fingerprint1ac1558c43b9ef5c4b338afe57f4c91995f9aa239d66d0a1a4ca876a89a237d5.
Scientific review remains pending. This is a short, low-resolution demonstration,
not a calibrated or statistically stationary turbulent ensemble.

A periodic unit cube has32^3 cells and4096 initially mass-seeded Monte Carlo
tracers. Uniform rho=1 and P=1/1.4 give initial sound speed1. The existing one-time
Fourier initializer adds zero-mean velocity with vector rms3, modes1..3,
spectral_slope2, projection parameter f_solenoidal=.5, seed24680. That parameter
mixes amplitudes; it is not a claim of equal solenoidal/compressive energies.
Hydro uses ideal gamma1.4, RK2/PLM/HLLC; no cooling, forcing, AMR or other physics.
Run to code time.2, about.6 initial box crossing times. The volume rms local Mach
sqrt(<v^2/(gamma P/rho)>) falls from3to1.65 and remains supersonic.

| Case | Completed | Open/suspended at end | Particles with completed episodes | Mean completed per particle |
|---|---:|---:|---:|---:|
| CFL.30 |7123|540|3661/4096|1.7390|
| CFL.15 |7623|530|3672/4096|1.8611|

Both have zero emitted partial events and zero invalid observations. Suspended
views are separate right-censored episodes; never count them as completed.
331baseline particles have no recorded episode including suspended ones. Every
particle survives, and canonical event identities are unique and contiguous.
The initial mass/pressure/Mach oracle passes at1e-12; mass and total energy errors
are at most2.22e-16. Tracking-on/off saved double fluid fields and native particle
observations are exactly equal. No tolerances were tuned to obtain these results.

Halving CFL gives7.02percent more completed episodes. That jointly changes fluid
integration, MC draws and membership sampling; it does not isolate detector cadence
or demonstrate convergence. All episodes have no_plateau and cadence_limited.
Baseline all-recorded flags include946merged flickers and351ambiguous directions;
flags overlap.12087potential downstream reentries were suppressed (counter, not
additional physical shocks). First-inside pressure-Mach median1.540 and10/90percent
quantiles1.113/2.988 are uncalibrated spatial estimates. Fixed endpoints are not
certified plateaus. No truth shock catalogue or exact particle crossing time is
available; do not infer recall or a physical waiting-time law from these counts.

## Inspect the results

Open index.html: five full-resolution plots have adjacent explanations. report.json
has all counts, flag totals, counters, quantiles and field evolution. episodes.csv
has every emitted/suspended identity, status, flags, observed timing bounds and
first-inside strength estimates. status0means completed and1partial; a suspended
view stays unfinished even if its post observation is valid. Blank pre/post times
mean no valid evidence. pre/first/last/post are samples, not exact crossing times.

selected-raw-evidence.zip contains all strict-reader-compatible event/suspended
chunks and inputs/logs for both active cases; four baseline volume snapshots and
half-CFL initial/final volumes. Full local raw runs remain in
/private/tmp/shock-turb-fixed-20261013 and are inventoried in raw-sha256.json.
The archive omits full thermo streams and other volumes; rerun to regenerate them.
The preserved old-source failure and fresh checks/reviews are under regression/.
The earlier28-plot PDF retains its original3de18d97source identity.

## Reproduce locally

From /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DAthena_SHOCK_TESTS=ON
cmake --build build-release -j8
MPLCONFIGDIR=/tmp/shock-turb-mpl \
/Users/hitesh/hitesh/git/shockit/.venv/bin/python tst/shocks-mc/turbulence_pilot.py \
  --executable build-release/src/athena \
  --workdir /tmp/shock-turb-rerun-new \
  --output /tmp/shock-turb-rerun-plots-new
open /tmp/shock-turb-rerun-plots-new/index.html
```

Work/output directories must be new. The runner executes three bounded processes
(180s cap each), checks exact non-interference/IDs/history/conservation, and exports
plots and evidence. Existing NumPy/Matplotlib/h5py and native readers are reused.
No new dependency, solver or detector-threshold change. LocalSerial CPU/double only;
RUSTY/Frontier, driven turbulence, MPI/GPU and scientific acceptance stay unqualified.
