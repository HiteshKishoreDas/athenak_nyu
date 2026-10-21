# Continuous Lagrangian shock carriers — living implementation record

## Accepted requirement, 2026-10-20

The user requires shock encounters to run only on continuous velocity-following
particles, never on Monte Carlo mass-flux tracers. This supersedes historical MC
carrier selections. Curve nodes already follow primitive velocity with SSPRK2.
New standalone `particle_type=lagrangian_tracer`, `pusher=lagrangian_tracer` must
solve the same velocity-following ODE. MC transport remains available with shock
tracking disabled. No implicit conversion of MC input or checkpoint is allowed.

## Scope and decisions

- New branch `codex/shock-lagrangian` starts at curve-history commit `55208570`.
  Durable isolated checkout: `.worktrees/shock-lagrangian` under the main repository.
  Main `feature_dump` and `codex/shock-multiphase` remain untouched.
- Reuse native deterministic seed schedules, tag allocator, real/int particle
  payload, boundary communication, thermodynamic output and paired shock restart.
  Keep `IsLagrangianMC()` strict; `IsFluidTracer()` selects shared infrastructure.
- Share multilinear cell-centered primitive-velocity interpolation with curves.
  Run standalone particles before each fluid RK2 stage and migrate after each push.
  Store predictor velocity in the existing interleaved velocity slots; native
  migration/restart carries every slot. The stage-two increment form avoids
  averaging wrapped coordinates across the periodic seam.
- Preserve continuous coordinates through AMR leaf ownership changes. Only MC
  uses native snapping. Finish regrid/ghost refresh before shock observation and
  checkpoint output. Use the same ordering with tracking on/off.
- Native particle restart header transport discriminator: 1=legacy MC, 2=continuous
  velocity tracer. Reject changed transport even when array dimensions agree.
  Shock contract explicitly identifies velocity interpolation and SSPRK2; source
  hash includes shared interpolation. Old MC shock checkpoints cannot resume.
- Seeding still selects cell centers deterministically, with mass/volume weights
  affecting initial placement only. Continuous motion has no stochastic hops and
  does not promise mass-following statistics in compressible flow. Scheduled births
  start at the actual insertion endpoint, not an unobserved scheduled time.
- Detector/encounter numerical rules remain unchanged. Detector width in leaf
  cells is an experimental AMR default, not scientific calibration.

## Pending validation gates

- [x] LT01: uniform velocity and time-dependent acceleration independent position
  oracles; periodic seam and leaf crossing; explicitly distinguish MC diffusion.
- [x] LT02/CSH02–04: evolved RH shock; one passage per particle, entry/exit timing
  brackets, strength evidence; independent native density-front comparison.
- [x] LT03/CSH05–06: tracking passivity, restart equality, corruption/overflow guards,
  scheduled insertions and transport mismatch rejection.
- [x] LT04/CSH07: dynamic refine/coarsen with standalone particles plus curves;
  position preservation, fresh newborn histories, exact continuation.
- [x] LT05/CSH09: native forced turbulence, continuous-particle/curve passivity.
- [x] Immutable test evidence and explained plots; commands and exact source binding.
- [ ] Scientific human sign-off (always pending until supplied).

## Limits to keep visible

Local Serial/double is the initial supported engineering backend. No RUSTY,
Frontier, GPU, threaded, or genuine multi-rank runtime claim follows from it.
Combined curve+particle multi-rank restart still needs compatible native output
sharding qualification. An RK predictor leaving an outflow domain aborts explicitly before deletion or
exit emission. Even ordinary outward drift can eventually hit this guard.
A corrected stage-two endpoint exit can still be right-censored. Full outflow
trajectory support remains unqualified; no exact crossing time is invented. Periodic tests qualify full two-stage trajectories. Detector shock
membership remains sampled at completed timesteps and may miss unresolved passages.
The inherited native forcing thermal-oracle issue remains scientifically unresolved.

## Commands and evidence

Configuration started with local Kokkos `08ceff92bcf3a828844480bc1e6137eb74028517`:

```sh
cmake -S . -B build-lagrangian -DPROBLEM=built_in_pgens -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian -j8
```

Detailed native run commands, input hashes, executable hashes, outcomes and failures
will be retained by the acceptance runners; no prior MC result qualifies new motion.


## Motion and output contract

Select both matching settings:

```ini
<particles>
particle_type = lagrangian_tracer
pusher = lagrangian_tracer
<shock_tracking>
enabled = true
lineage = choose-a-new-stable-lineage
```

The continuous particle solves \(d\mathbf X/dt=\mathbf u(\mathbf X,t)\).
For a step beginning at \(t_n\), sample the filled primitive field before each
fluid stage. The production update is
\[
\mathbf X^*=\mathbf X^n+\Delta t\,\mathbf u^n(\mathbf X^n),\qquad
\mathbf X^{n+1}=\mathbf X^*+\tfrac12\Delta t
[\mathbf u^*(\mathbf X^*)-\mathbf u^n(\mathbf X^n)].
\]
Migration follows each particle stage. The second expression uses increments,
so periodic wrapping never averages unrelated coordinates on opposite faces.
For a discontinuity this samples the finite-volume/interpolated velocity; it
is not an exact mathematical trajectory through an infinitesimal shock.

The velocity slots are RK scratch/sampled stage velocities, not an independently
integrated particle momentum. Thermodynamic fields and shock membership are
sampled from the containing cell at the synchronized completed endpoint. These
cell-based quantities need not equal interpolated velocity at arbitrary positions.
The native `.thp` files preserve full double-precision positions and unique tags;
VTK position fields remain compatible. Initialization/insertion may select a cell
center, but subsequent motion is continuous and AMR never snaps it to a center.

Particle events retain the existing 124-column strict core; optional traces retain
28 columns. The event identity remains `(lineage, particle_id, sequence)`. Curves
retain their independent node namespace, birth metadata, 133-column events and
33-column traces. A newly inserted curve node never copies parental state.
Inside births begin incomplete, with entry unknown. Every generation has an
immutable chunk, including empty chunks. Right-censored terminal views remain
separate from committed completed episodes. Output capacity exhaustion fails
explicitly rather than dropping events.

Paired restart requires the native fluid/particle checkpoint, `.shock.rankR`,
all event chunks through its generation, and, when curves are present, native
curve geometry and `.curves.shock.rankR`. Same source/configuration, ABI and
partition only. Do not change transport to resume an MC checkpoint. Do not rename
or delete a committed history prefix. A new basename permits deterministic replay
without overwriting an existing immutable checkpoint pair.

## Review fixes and rationale

Three independent read-only reviews covered timing, persistence and minimal
intrusion, using ponytail-review. Applied all concrete recommendations:

1. Make boundary handling stage-aware. A predictor is not a completed endpoint;
   reject unqualified predictor outflow before fabricating physical-exit evidence.
2. Fine-to-coarse diagonal crossings can intentionally lack a diagonal entry in
   the native neighbor table. Reuse the existing global containing-leaf locator
   only for unresolved entries; append their MPI sends before compaction. Ordinary
   neighbor lookup remains the fast path. Rare fallbacks log rank and count.
3. Bind the remap caller and shared interpolation in the source fingerprint.
   Bind resolved AMR controls/criteria and timed callback parameters in restart
   contracts. Freeze again after native refinement-criterion defaults are added;
   this fixed a real active-AMR restart mismatch.
4. Delete unreachable MC shock transfer instrumentation. Preserve native MC draws,
   hops and tracking-off behavior. No parallel tracer/output framework was added.
5. Strengthen tests: all independently swept particles must produce one event;
   exercise nonzero 3D motion and z wrapping; hash actual restart input checkpoints;
   require proof the missing-diagonal fallback executes.

## Final engineering results

Final local source is rerun after the final formatting changes. Exact paths and
hashes live in `evidence/shock-lagrangian/2026-10-20/` and the durable bundle below.

- LT01: 2D and 3D uniform controls, 128 particles each, all cross a leaf boundary;
  3D includes 128 z wraps. Maximum position error \(2.22\times10^{-16}\).
- LT02: time-dependent acceleration and affine shear, every saved position compared
  against independent solutions. Maximum error \(8.88\times10^{-16}\), unchanged
  tolerance \(2\times10^{-12}\).
- LT04: actual refine/coarsen transactions, including mixed levels and simultaneous
  diagonal motion; 128 particles exercise containing-leaf fallback. No position
  snapping. Same unchanged roundoff budget.
- LT05/PHY08/INT06: native evolved RH shock. All 128 seeded particles cross the
  native density midpoint exactly once and record exactly one completed encounter.
  Every event evidence envelope covers the independent sampled crossing bracket.
  Mach estimates are saved as strength evidence, not certified plateau truth.
- LT06/IO01/IO06: outside/active/pending replay is exact. Active particle encounter
  checkpoint at generation 68 on a refined mesh continues through coarsening with
  exact native samples, event identities/evidence and terminal state.
- LT07/EVT02: 128 scheduled particles born inside at actual insertion time
  `0.1006776070884923`, not requested schedule time `0.1`; fresh sequence zero and
  left-censor flags. No prebirth samples.
- LT08/CFG01: MC histories, mismatched pushers and transport changes reject.
  Native MC tracking-off remains runnable and passes its existing moment oracle.
- LT09/CFG06: zero particles yields no invented observations/events.
- LT10/IO08: event overflow and corrupted committed history fail explicitly.
- LT11: predictor-outside/corrector-inside diagnostic fails closed before a false
  physical-exit record can be emitted.
- 20/20 existing curve gates pass on continuous-carrier inputs. Combined live AMR
  now passes: 39 positive-time newborn curve nodes, 17 inside and 22 outside.
- CSH09: accepted native forced \(32^3\) box, 4096 continuous particles and 16 curve
  nodes. Tracking/trace toggles leave native fluid arrays, particle samples and
  geometry exactly unchanged. Curve histories contain 28 complete events.
- 53 detector parity fixtures pass against the pinned exact ShockIt source, with
  exact masks and unchanged absolute/relative tolerance \(10^{-12}\).
- Encounter core EVT01–08/11–14 and strict thermo/curve reader checks pass.

These are engineering results. Scientific sign-off, converged turbulence,
performance budgets, RUSTY and Frontier remain pending. GPU/threaded/single
precision and arbitrary AMR criterion/physics combinations are not qualified.

## Retained failing attempts and corrections

No tolerance was weakened and no failed attempt was relabeled:

- `core-first`: started before the prescribed build existed; its kinematic/AMR
  checks fail with missing executable. Rebuilt first, then reran.
- `lagrangian-first`: birth test exposed a reader omission of the saved flags.
  Reader now exposes the existing field; malformed schema handling remains strict.
- `lagrangian-reviewed` / `-release` / `-qualified` / `-accepted`: stronger mixed-level
  fixture initially did not exercise the diagonal fallback. Inspection showed the
  fixture input repeated `<problem>` and the input editor appended `refine_time=0`
  before a later unchanged `refine_time=.1`, so AthenaK correctly used `.1`.
  Fixed the shared test editor for repeated sections and added an explicit regression.
  Require the native fallback log, not a favorable trajectory result alone.
- Active-AMR replay first failed because native criterion defaults were absent
  from the original constructor-time contract but present after restart. Freeze
  after defaults resolve. Refined active histories now replay exactly.
- First parity invocation lacked SciPy in the visualization environment. Reused
  the existing ShockIt environment and ran all 53 fixtures; no dependency added.
- Multi-rank launcher fails before AthenaK: PMIx runtime `0x60100` versus PRRTE
  build `0x50009`. MPI compilation/singleton are separate evidence, not two-rank
  validation. Combined MPI checkpoint sharding remains unsupported.

## Reproduce locally

From the isolated branch checkout, initialize its exact Kokkos submodule first.
Use an existing Python environment containing NumPy/Matplotlib for plotting and
forcing analysis; parity additionally needs the existing ShockIt/SciPy environment.
The standalone runner otherwise uses the standard library.

```sh
cmake -S . -B build-lagrangian -DPROBLEM=built_in_pgens -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian -j8
cmake -S . -B build-lagrangian-prescribed -DPROBLEM=material_curve_test
cmake --build build-lagrangian-prescribed -j8
python tst/shocks-lagrangian/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/local-particles-UNIQUE
python tst/shocks-curves/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/local-curves-UNIQUE
python tst/shocks-curves/acceptance.py --suite forced \
  --executable build-lagrangian/src/athena \
  --workdir build-lagrangian/evidence/local-forced-UNIQUE
python tst/shocks-lagrangian/visuals.py \
  --lagrangian build-lagrangian/evidence/local-particles-UNIQUE \
  --core build-lagrangian/evidence/local-curves-UNIQUE \
  --forced build-lagrangian/evidence/local-forced-UNIQUE \
  --output build-lagrangian/evidence/local-plots-UNIQUE
```

Each runner requires a new destination and retains inputs, full logs, checkpoints,
raw histories and executable/input/configuration fingerprints. The viewer places
explanations next to all seven figures, including retained development failures.
Current durable evidence/viewer:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/lagrangian-shock-carriers-20261020-final/`.

MPI compile and singleton succeeded on the final source. Singleton native samples,
particle event/terminal records and curve events equal the Serial control exactly.
The first MPI invocation used an incorrect relative executable path and failed
before launch; its log is retained. The corrected singleton command passes.

## Next gates

1. Human scientific review of the exact-code bundle and reported physical oracles.
2. Agree physical-versus-cell-scaled detector width near AMR interfaces; current
   fixtures retain width four leaf cells as an experimental default.
3. Repair/qualify the MPI launcher and combined checkpoint sharding; then RUSTY
   CPU/MPI scale and performance budgets. Frontier GPU requires separate support.
4. Full stage-aware outflow trajectory treatment if requested. Current predictor
   guard deliberately prevents silent false encounters or particle loss.

## Sealed exact-code bundle

Implementation commit: `1c325e5a59d64f743ffdaa4f8678a4805b12608f`.
The subsequent metadata commit changes documentation/evidence only.
Bundle manifest SHA256:
`1e1a497fdb640d64393f79d448d0eb9e325fe0967f7cd6dcc8a7f77279bcc090`.
Verification passes: 332 archived native source files match, 243054 raw archive
members are readable/safely named, and all 130 final archived input/log/report
bindings match their reports. All final native reports match the retained exact
executables. Bundle contains native raw runs, failed attempts, source/dependency
archives, three runtime builds, compiler configurations and seven explained plots.
No scientific sign-off supplied. No new branch push performed in this task.

## Combined visual showcase — 2026-10-21

[Setup, output-only replay, exact detector checks, caveats and reproduction](curve-shock-showcase-20261021.md).
The movie overlays native material-curve geometry and actual node shock flags on
fluid snapshots, then shows retained per-node passages. Suppressed raw flags are
explicitly visible in purple. Analysis only; original scientific implementation
and sealed bundle unchanged. Compact execution/render/replay/media evidence:
`evidence/shock-lagrangian/2026-10-21/`.
