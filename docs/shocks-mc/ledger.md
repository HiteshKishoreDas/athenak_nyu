# Evidence ledger — multiphase MC integration

2026-10-08. Accepted target: user explicitly requested a branch from
`multiphase_turb` and implementation there. New branch `codex/shock-multiphase`,
base `526a187f938e38634e89a88192ebe9731c4d36bb`, Kokkos
`08ceff92bcf3a828844480bc1e6137eb74028517` (4.4.0).
Worktree `/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu`.
Original repository and prior `codex/shock-encounters` worktree retained.
The earlier branch's results do not qualify this new MC implementation.

The actual existing implementation is `particle_type=lagrangian_mc`,
`pusher=lagrangian_mc`; enum `lagrangian_tracer` remains unimplemented.
No velocity tracer was imported. Native MC flux transport, stateless hash RNG,
particle thermo output, seeding schedules and native restart are reused.
Historical assessment/earlier branch are evidence, not new instructions.

## Source/environment baseline (M0/G0)

Baseline Debug: `cmake -S . -B build-baseline -DCMAKE_BUILD_TYPE=Debug` then
`cmake --build build-baseline -j8`; exit0. Apple clang17.0.0, arm64 macOS,
CMake4.1.2, Serial double. Frozen reference includes manifest/source hashes;
local_pass.py SHA2566f5ec8bc89c90ecb1ad2e0135fb8abe7a79f2b43badfc73356fd4ae3f1cc3857.
Complete baseline command/input/executable identities are in
`evidence/shocks-mc/baseline/report.json` and per-run manifests.
Eight fixed seeds,4096tracers each,117steps: independent Bernoulli mean and variance
within predeclared6SE; no deterministic particle path assertion.
Fourblocks:1024IDs survive, all1024changeGID; restart compares102400observations
exactly; zero particles succeed. This is bounded statistical evidence, not a proof
of independent RNG streams. Baseline RNG has additive key aliases (review packet).

Initial strict baseline reader rejected existing duplicate terminal thermo rows;
only identical replay rows are now accepted explicitly for this legacy output.
Conflicting duplicates still reject; shock event identities remain strict.
Original failure and all commands/results retained.

## Numerical/contract checks (M1–M2)

`build-shock/shock_encounter_test` passes EVT01–08/11–14 typed-state checks.
`python tst/shocks/parity.py build-shock/shock_detector_test --reference
 evidence/shocks/reference --report <report>` passes53fixtures, exact masks,
rtol=atol1e-12. Pressure/temperature/entropy, gradients, nearest-axis endpoints,
physical filters and both Mach estimators are compared, with independent RH
and smooth expectations. The first attempt lacked the copied pressure-gradient
 overflow input; restored frozen fixture, no changed numerical tolerance.

## Integration/review discoveries and corrections

First production run died before timestepping (SIGSEGV). This branch allocates
Hydro::psrc only when needed; integration incorrectly assumed nonnull. Guard fixed.
The unsupported old problem-generator hook name also caused an initial build
failure; fixed to the actual branch API. Initial build/run failures are retained.
Independent review found per-rank sidecar resolution incorrectly using rank0's
payload, missing shear/orbital rejection, duplicate seeding IDs that corrupt native
restart schedule matching, self-incompatible oversized contracts, unchecked native
MC writes/close, per-node restart incompatibility, unsafe count products.
Fixes explicitly resolve rank payloads, reject unsupported configurations/IDs,
bind metadata bounds, check native particle/metadata write counts and close status,
and fail before committing a shock sidecar. Fluid payload already checks counts.
Inherited per-step MPI datatype leakage is closed after metadata exchange;
metadata count sums are range checked. Disabled diagnostic error view allocation
removed. Suspended output is also bounded by event_capacity.

Pending results are appended as final runs/review complete. No GPU or remote cluster
was used. RUSTY/Frontier are requested deployment targets, not tested hardware.
