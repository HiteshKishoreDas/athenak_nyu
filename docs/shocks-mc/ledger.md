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

## Final supported engineering evidence (78e5 source)

Production commit `38c55570`, based on `e58c482d` M0 baseline record. Core source
fingerprint `78e5cb1d425daaa8ab2e06deaa40db6c8249de79cd1480c7d0b069fd1587f7c3`.
Independent review verifies26source files and closes identified source blockers;
it does not supply runtime/scientific approval. Exact runtime commands/configs,
executable and build-cache hashes appear in the linked JSON run manifests.

| Milestone/test scope | Actually executed | Result/remaining gate |
| --- | --- | --- |
| M0/G0 | Baseline Debug;8seeds4096MCtracers; zero/fourblock/restart | Statistical bounds and102400exact continuation records pass |
| M1/G1 | Frozen numerical/MC/support definitions and independent RH fixtures | User branch/tracer accepted; physical release tolerances remain experimental |
| M2/G2 | Debug/Release/MPI-compiled native state tests;53case parity suites | Exact masks and1e-12 field tolerances; no changed oracle |
| M3/M4 |20finalRelease integration runs; synchronized .thp comparison; TRC11source/probability/draw/time-bin evidence | Exact on/off fluid/MC records and valid timing pass |
| M5 pilot |8M2seeds128tracers each; independent conditional analytic-front bins |128events/128unique parent-front change bins each; no exact material crossing time claimed |
| M6 local |18actual-boundary fixtures; OUTSIDE/ACTIVE/PENDING restart;12reject and4replay processes;13active domain-exit records | Complete state and local replay/censoring pass; real cross-rank and storage failure remain untested |
| M7 local |3build modes;6existing linear-wave runs each on baseline/final; noninterference; style and encounter sanitizers | Unchanged wave thresholds/error tables pass;0introduced style findings; limited sanitizers pass |
| M7 MPI |Successful MPI build/directsingleton; actual2rank launcher attempt |2rank exit213 before AthenaK: PMIx runtime6.1.0 versus PRRTE build5.0.9; no multirank/GPU pass |
| M8 measured |21heldout physical runs,4frontplots, local small-job timings | Physical acceptance/budget pending; strong/boosted association failures retained |
| M9 engineering |Source review, docs/schema, scoped evidence ledger/manifest | Human exact-code sign-off pending; no release/push |

Final detector command:
`python tst/shocks/parity.py build-release/shock_detector_test --reference
 evidence/shocks/reference --report <report>`; field tolerances unchanged.
Final integration/qualification commands are in `evidence/shocks-mc/README.md`.
`python tst/shocks/system_mc.py --executable build-release/src/athena --baseline
 build-baseline/src/athena --mpi build-mpi/src/athena --workdir <newdir> --report
 <report>` runs the bounded existing hydro subset, small local timings and launcher
checks. Existing hydro regression still names removed parameter `vflow`; the
wrapper maps it to the documented equivalent `vx0` without changing thresholds.
The initial stale-runner failure was retained. Whole-repo style is not clean:
20warnings/19unique inherited findings remain; baseline22unique, introduced0.

Final full-precision replay failure corrected in production: GetOrAddReal encoded
FLT_MIN as1.17549e-38, while active EOS used1.1754943508222875e-38. SetRealPrecise
now writes bound native floors/seed defaults with max_digits10; exact comparisons
remain unchanged. All authentic all-phase continuation cases pass again. Earlier
failed crash/reader-key/floor-serialization fixtures are retained compactly.

HeldoutM5:134events/128unique, six IDs have2episodes, although actual MC parents
show128analytic-front change bins. BoostedM2:132events/127unique plus1suspended;
five IDs have2episodes. Its MC parents show209positive-to-nonpositive front-distance
bins, illustrating stochastic return ambiguity;200bins lie in retained event
 envelopes. These bins describe observed cell parents, not continuous material
paths. M1.2:128events, pressure-Mach error range[-.103117,-.0156752]; M1.5:
127completed+1suspended at cutoff. No detector/exit thresholds were retuned to
remove these outcomes. Strength comes from first-inside fixed endpoints and is
not a certified plateau measurement. `parent_front_bins.json` records the
independent observable-bin definition and every match.

Small quiet-run whole-process medians (3trials,1024particles,64×4/256×4cells)
full/off ratios 3.229/2.599. Includes startup, detector/tracking and
durable generation chunks; empty-particle controls are reported separately.
No production throughput/scaling or acceptable overhead budget is inferred.
RUSTY/Frontier were not accessed; node/backend/rank/cell/particle scale and allowed
runtime/memory/I/O overhead are outstanding deployment decisions.


## Clean archived reproduction and final packet (M9)

Archived `42bcf12c` using `git archive`, extracted to
`/private/tmp/shock-mc-clean-42bcf12c`, and copied the clean pinned Kokkos source
without `.git`. Configured Release Serial with `-DAthena_SHOCK_TESTS=ON` and
built all targets with `-j8`; exit0. Compiled source fingerprint remains78e5.
The clean executable passes the native encounter state suite, all53parity cases
with unchanged tolerances, and all20serial integration processes including
all-phase restart and feature-off noninterference. Exact commands, logs, input
copies, retained fixtures, executable/cache hashes and archive provenance are
in `evidence/shocks-mc/clean-reproduction/`. This is a clean full build and core
suite reproduction, not a repeat of every heldout/cluster gate.

Final packet includes explicit four-mode CMake caches/generated configuration,
actual two-rank pre-launch failure log, source hashes verified against the
independent review, and a noncircular file manifest (manifest itself excluded).
No production source changed after38c55570. Human review, broader physical
association, multi-rank/device execution and deployment budgets remain pending.


Packet collection notes: the first final-manifest invocation used system Python
without NumPy and aborted before writing a manifest; rerun with the recorded
ShockIt virtual environment. Raw copied CMake caches/generated headers/configure
logs retain their exact trailing whitespace/terminal blank lines; whole-packet
`git diff --cached --check` therefore reports these generated-artifact findings.
Source/tests/documentation diff checks remain clean. Artifacts were not normalized,
so their saved hashes continue to match the build inputs actually used.


## Ponytail audit and minimal constructor cleanup — 2026-10-08

User asked whether the feature minimally changes AthenaK machinery. Audited the
full production diff and traced configuration, sampling, transfer, compaction,
outflow and paired restart callers. Solver/flux/mesh/task graph files remain
unchanged against526a187f. Most feature code is new standalone detector/state/
observer code. Driver adds14lines; main9; restart23additions/2deletions. The
largest native integration addition is295lines of guarded particle boundaries.
Its separate enabled path handles physical outflow and typed-state compaction;
merging it into native AMR/general transport would broaden the behavior changed.
No framework or dependency was added, and native MC transport/RNG remain reused.

Cleanup `f0053f2e`: moved the exact existing configuration body out of particles.cpp
into private ConfigureShockTracking(pin), defined in existing shock_particles.cpp;
removed an uncalled InitializeShockState method and its declaration. Native
particles.cpp now adds only one call plus a blank line against the selected base.
No numerical/state/migration/restart algorithms or interfaces changed. This is
less intrusion in core files, not a claim of globally minimal code or better
runtime cost. The pre-cleanup statement about unchanged source since38c55570 is
historical; current fingerprint is80da5390b07a322b83f3ce829f97975c4326b3216a764d77bf530a8f355265c4.

Independent read-only reviewer reconstructs all three old source files exactly
by reversing extraction/include relocation/dead-method deletion. Configuration
body SHA2564ebfaebf08f7e0e4e4f32c3a064d97bec86bc79bd61103b084e5f890bf0cd129
is unchanged, called at the same constructor point. Typed record layout is unchanged.
Builds Debug/Release/MPI exit0; native state tests pass on all3. Release53parity,
20integration,12rejection and4replay checks pass. All21heldout measurements,
including full recorded event evidence, exactly equal previous reports.
Scientific failures and sign-off stay pending; thresholds are unchanged.
Native style wrapper's existing --filter=-build/include_subdir reports0findings
on the three changed files; unfiltered first invocation flags native include
conventions. Diff check on source/docs passes. Commands/input/log/build config/
review hashes are under evidence/shocks-mc/ponytail; no full raw grids committed.

No new cross-rank, boundary-fixture, exit, wave, performance, sanitizer or clean
archive run is claimed for80da: those prior reports retain78e5 and the affected
logic is unchanged by exact source extraction review. Exact source-bound restart
contracts intentionally require the same code fingerprint; older checkpoints
must use their matching executable. This cleanup does not add cross-version
checkpoint compatibility. RUSTY/Frontier and scientific release remain unqualified.


## Three requested Ponytail reviews and implementation — 2026-10-08

User requested exactly three subagents with different prompts and implementation
of suggestions. Numerical/state, native integration, and persistence/output agents
read ponytail-review and relevant flow, produced read-only findings, then independently
reviewed root's sequential integration. Findings were consolidated: duplicate state
representative derives from first; duplicated native restart structs share one same-file
type; remove writer tautology but retain reader schedule comparison; remove only exit-index
sort (final event sort stays); reuse existing ConfigValid and ShockError; share PutEvent;
append nonempty contiguous sidecar state once. No framework/solver/RNG/task changes.

Source commit3de18d97; final fingerprint
3feb3bcdf77d21296516544e67b5ad75f6ac6ca7bc17c4b48be886ddc0fb7cbc.
Forty net production lines removed. Fixture now initializes four stored observations;
all remaining state fields still checked. Compiler ABI probe: native header48bytes,
Observation184, old/new EncounterState1000/816, Event944. All surviving field offsets
and native header offsets match. State bytes/source contract explicitly reject older
sidecars; actual old active checkpoint rejection passed. Public event schema stays124.

Three builds and native state checks pass. Release53parity cases (unchanged1e-12),
20integration,12reject,4replay processes pass. Both actual-method boundary builds and
all18cases pass; actual13domain-exit records exactly match earlier evidence. All21
physical measurement dictionaries match exactly;23748event/suspended/trace/transfer CSVs
are byte-identical after replacing only the implementation fingerprint. Native particle
restart tails are identical in OUTSIDE/ACTIVE/PENDING. Scientific extra episodes remain.

No new style findings: source-text/category comparison has10inherited native MC warnings.
Initial reference comparison failed because temp reference lacked CPPLINT.cfg; rerun
with exact native configuration (same90-column policy/filters) passed. A double blank
line introduced with the shared header was removed before final runtime suites; its
exact reverse reconstruction matches the independent review's hash. No numeric test
threshold, state transition or error-substring was changed. Configuration errors now
use shared fatal prefix/MPI abort; no rejection guard was removed except a tautology.

review3/ contains reviewer findings, accepted changes/hashes, exact commands/configs,
inputs/logs, native ABI/old-layout probe and equivalence reports. Old packets retained
with their own identities; raw grids/binaries are not committed. No new cross-rank,
device, cluster, sustained performance or human scientific sign-off is claimed.

## Requested local2D execution and28plot gallery — 2026-10-08

Fresh commands used existing Release executable and ShockIt Python environment:
`integration.py --executable build-release/src/athena --workdir /private/tmp/shock-mc-laptop-20261008/integration --report /private/tmp/shock-mc-laptop-20261008/integration-report.json`; matching qualification.py arguments use qualification paths.
20engineering passed13.8254s;21physical runs measured28.0548s (run success, not physical approval).
`MPLCONFIGDIR=/private/tmp/shock-mc-mpl .venv/bin/python tst/shocks-mc/laptop_visuals.py --executable build-release/src/athena --workdir /private/tmp/shock-mc-laptop-20261008 --output evidence/shocks-mc/laptop-visuals --reuse`
completed28plots and18sequential benchmark runs. `.venv/bin/python` resolves to
`/Users/hitesh/hitesh/git/shockit/.venv/bin/python`. An initial native snapshot
preflight and failed plot attempt remain in the raw root and/private/tmp. First
Matplotlib failure/repair are recorded in laptop-tests.md and verification.json.
Exact restart equality and all numerical thresholds were retained. Native state
unit passed;210initial packet hashes verified. Benchmark modes rotated order;
medians and every actual command/input/log/digest are in report.json/runs.
Raw1.2GB; compact9MB.12representative figures visually inspected.
Known Mach5/boost/M1.5 outcomes match review3; no production source edits.
Ponytail reused test/readers/installed libraries; no new framework.
Scientific acceptance, resident-memory measurement, sustained performance and
RUSTY/Frontier/MPI/GPU qualification remain pending.

Fresh end-to-end CLI reproduction without `--reuse` passed: 20 engineering runs,
21 physical runs, 28 plots and 18 additional sequential benchmarks. All 21 physical
measurement dictionaries and 27 PNGs match exactly; runtime plot has fresh measurements.
An evidence-packaging assertion initially expected 26 matching figures; actual 27
includes output cost. Corrected count, no numerical/test tolerance changed.
Reproduction reports are in laptop-visuals/clean-reproduction; verification.json
records exact commands and equality scope. Compact native 2D snapshots are explicitly
tracked despite the repository binary-output ignore rule, so the packet is portable.
