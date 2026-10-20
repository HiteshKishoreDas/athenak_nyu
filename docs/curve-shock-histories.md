# Material-curve shock histories — requirements

## 2026-10-17: initial questions (historical)

User requested implementation, with all questions asked before coding. No feature
implementation has started. Current feature_dump base is f27e40a0. Preserve
codex/shock-multiphase unchanged.

Accepted birth rule: newly inserted nodes start their own observation/encounter
history at insertion time. Never copy parents' encounters or fabricate pre-birth
observations. Stable node identity, ancestry and birth time remain distinct.
The treatment of a node already inside a shock at birth is pending below.

Recommended approach: reuse the existing ShockIt-contract detector and encounter
engine, preserving bracketed times and validity/censoring evidence. Existing curve
advection and topology stay intact. Store complete history state through output
and restart. Proposed settings and physics support below are not accepted yet.

Questions for the user:

1. Physics: must the first version include forced supersonic turbulence? Recommend
   uniform ideal-gas Newtonian hydro with forcing included; qualify AMR, cooling
   and MHD separately unless required now.
2. Coexistence: must MC particle histories and curve histories work in the same
   run? Recommend yes, with distinct identities and passive behavior checks.
3. Product: independent per-node histories only, or whole-curve diagnostics too
   (for example instantaneous arc-length fraction in shocks)? Recommend per-node
   history first, curve summaries derived during analysis.
4. Born inside: retain a clearly incomplete/left-censored encounter starting at
   birth, or wait for exit before recording later encounters? Recommend retain
   the incomplete encounter without inventing its entry time.
5. Output: save each node's every-step shock flag as well as event histories by
   default, or make detailed traces optional? Recommend always save events and
   current restart state; every-step traces opt-in. Detection itself remains
   synchronized every complete step and at node birth.
6. Detection: reuse existing MC detector thresholds, strength evidence and
   encounter/rearming rules initially, or request changes? Recommend reuse as
   explicitly experimental defaults, with independent curve-path tests.
7. Existing restarts: must users be able to enable tracking when restarting an
   older curve-only checkpoint? Recommend fresh tracking runs first; support for
   starting new history at restart requires an explicit provenance policy.
8. Deployment: laptop CPU first, RUSTY CPU/MPI, or Frontier GPU required for the
   first supported version? Approximate curves/total nodes and acceptable runtime,
   memory and I/O costs are needed for deployment-scale qualification. Recommend
   laptop engineering gates followed by RUSTY MPI; Frontier remains untested
   until actual hardware evidence exists.

Retain the existing native RK2 forcing thermal-oracle failure as an explicit
physical limitation. Engineering implementation and tests do not supply human
scientific approval. No new physics values, floors or acceptance tolerances
were chosen in this requirements discussion.


## Accepted answers — 2026-10-17

- Include AMR initially, alongside the forced uniform-hydro scope discussed.
- MC particle and curve-node histories must coexist in one simulation.
- Per-node histories first; whole-curve summaries are derived in analysis.
- Initial/inserted nodes already inside shocks retain incomplete left-censored
  encounters beginning at initial observation/insertion. Never copy parent events.
- Always retain events and complete restart state; per-step traces are optional.
  Observe at every completed timestep and at actual insertion.
- Reuse existing detector thresholds, strength evidence and encounter/rearming
  rules initially; independently test curve crossings and AMR effects.
- Fresh tracking runs first; do not silently activate tracking from old curve-only
  checkpoints or invent earlier histories.
- Laptop checks first, followed by RUSTY CPU/MPI qualification. Node counts and
  resource budgets are not yet supplied. Frontier GPU remains outside first gate.

One AMR numerical choice is pending: width in local leaf cells (recommended,
matching existing MC definition) versus fixed physical distance. Continue storage,
birth and uniform-path work independently; do not silently choose AMR scale.

Implementation gates: freeze consumers/state interfaces; retain failing fixtures;
implement node observation/birth tracking; preserve complete curve/MC state through
regrid/migration/checkpoints; qualify uniform and actual refine/coarsen cases,
coexistence, passivity and replay with independent physical references. Existing
MC history rejects AMR and must be extended rather than bypassed. Human scientific
approval and actual RUSTY/hardware performance remain separate.


## Implementation ledger — gate CSH00

Worktree `/private/tmp/athenak-curve-shock-histories-20261017`, branch
`codex/curve-shock-histories`, from requirements commit dd159953. Original
`codex/shock-multiphase` remains 2450e1c2. Reused exact local Kokkos commit
08ceff92bcf3a828844480bc1e6137eb74028517 through a separate Kokkos worktree;
no network/dependency install. No implementation has been declared validated yet.

Frozen ownership: parent owns curve observation/configuration, sampler hooks,
driver/restart integration and headers/build lists. Migration task owns MC AMR
state/boundary/remapping files. Persistence task owns new shared durability helpers,
curve event/trace/checkpoint IO and strict reader. Independent oracle task owns
new fixtures/acceptance runner only. Shared interfaces are declared in
material_curves.hpp; changes integrated sequentially, no competing file edits.

Red baseline /private/tmp/curve-shock-red-v3-20261017 proves actual old binary
rejects curve/MC coexistence and curve-only history requests. Independent tracking-off
native RH control exits0 after107steps:16anchors each have exactly one density-front
passage bracket, without consuming event output. These are moving fluid density
front observations, not invented exact Lagrangian crossing times. Detection scale
at AMR interfaces remains a pending user choice; do not promote a default silently.

In-progress code reuses detector/encounter engine, optional node-ID-owned state
(no geometry/parent copying), native distributed curve sampling route and separate
curve history output namespace. Birth observation is immediate at initial/insertion
state; born-inside is left-censored. Completed-step observations follow completed
AMR transaction and refreshed ghosts. MC state must migrate/regrid intact before
AMR exclusion is lifted. Existing test-only pgen gains a native ShockTube initializer
mode with its already-present timed AMR callback; no shock physics is reimplemented.


## Implemented slice and why — 2026-10-17

This is an engineering increment, **not completion of the accepted AMR goal**.
Production lives in the isolated `codex/curve-shock-histories` worktree. The source
branch `codex/shock-multiphase` remains exactly 2450e1c2; `feature_dump` remains
requirements commit dd159953 until this increment is integrated deliberately.

- Reuse the unchanged C++/Kokkos detector and encounter state machine. One small
  observation converter serves both carriers; no second shock algorithm.
- Store optional state by immutable node ID, outside Node geometry. Adaptation
  cannot copy parent events. Initial and inserted nodes receive their first
  observation at actual birth; inside births have unknown entry and left censoring.
- Extend the existing distributed leaf-owner sampling route rather than create
  a second mesh communication framework. Reuse native synchronized primitive
  ghosts. Save current leaf level and cell spacing as detector-scale evidence.
- Sample complete steps after native regridding, primitive/ghost refresh and curve
  adaptation. Stage-one insertion observes only newborn IDs at the current clock.
  MC observation uses the same advanced clock/generation. Uniform MC regressions
  verify that this ordering change preserves existing observations and trajectories.
- Share the existing immutable output helpers and unchanged MC event core. Curve
  events have a separate directory/schema and typed node IDs, birth and mesh metadata.
  Python identities explicitly include `curve_node` to avoid cross-carrier numeric
  ID collisions. Every-step trace output is optional; events/state are mandatory.
- Pair complete node history state with native fluid and geometry checkpoints.
  Bind source/configuration, ABI, exact IDs/births, rank partition, both native
  files and committed history prefix. Reject corrupt/missing/incompatible state,
  disabled consumers and silent activation from older untracked checkpoints.
- Preserve complete MC encounter state through the existing remap ordering and
  typed MPI payloads; resize detector storage after block-count changes. This
  infrastructure does **not** repair the inherited MC mass-allocation failure.

Accepted requirements above are distinct from unchanged experimental detector
parameters in fixtures: width 4 leaf cells, minimum Mach 1.1 and existing rearming
rules. No physical threshold, floor, numerical tolerance or exact crossing time
was invented to obtain a pass. AMR width policy still awaits the UserQuestion.

## Acceptance ledger — locked production source

Built-in source fingerprint:
`e0e18be4dddc954da31130e884456be8585bf5615e6f7ac789438ad2efe2ddcf`.
Test pgen fingerprint:
`580f268d9337bc1a8826d4d0dc6705ca95a67e58689282e63def837cdafd6d77`.
Both were recomputed from the CMake source list and found in each executable.
Build caches, compiler commands, source hashes, logs and exact inputs are retained
in the evidence bundle. Detector/encounter headers are byte-identical to 2450e1c2.

| Gate | Actual result | Limit / next gate |
|---|---|---|
| CSH00 interfaces/source/build | Serial/double Release built-in and test-pgen builds pass | No new dependency; exact local Kokkos 08ceff92 |
| CSH01 zero/translation/shear | Zero events; maximum affine position error below existing \(2\times10^{-12}\) budget | Prescribed-field engineering control |
| CSH02–03 planar, with/without MC | 16 nodes, 16 complete passages, no extra episodes; all 16 independent density-front brackets covered | Strength errors reported, no owner acceptance tolerance invented |
| CSH04 initial inside | 16 partial, left-censored encounters | Unknown physical entry remains unknown |
| CSH04/07 curve-only live AMR | Actual leaves \(2\to8\to2\); 39 positive-time births, 17 inside and 22 outside; every earliest sample at birth/sequence zero | Local-cell detector experiment, not approved AMR mask parity |
| CSH05 passivity | Native fluid, MC trajectories and curve geometry exactly equal on/off; events unchanged with traces off | Same partition/seed/configuration |
| CSH06 replay | Actual outside/active/pending checkpoints; exact events, full traces, terminal views, MC continuation and node state; relative-prefix `-d` pass | Same ABI/partition; fresh tracking runs first |
| CSH06/07 newborn replay | Checkpoint at generation 13 with active inserted nodes; exact full/birth traces, events, terminal views and final nodes | Curve-only actual refine/coarsen |
| CSH06/08 rejection | Missing/corrupt/future sidecars, prefix corruption, oversized future chunk, changed native source/refinement callback, disabled curve consumer, old untracked activation, nonfinite curves-only forcing, event overflow reject | No silent loss/relabeling |
| CSH07 smooth AMR | Actual \(4\to16\to4\) leaves, 14 positive-time births, zero events; exact leaf geometry audit | Native ghost contract differs from independent ShockIt AMR reader |
| CSH09 native forcing + both carriers | Three runs pass exact fluid/MC/curve passivity; 28 complete curve events | Does not repair inherited RK2 forcing heat-oracle failure |
| Required coexisting live AMR | **Fails explicit MC startup guard** | See inherited occupancy failure and pending physical decision below |

Final serial acceptance reports 19/20 checks passed and returns nonzero overall.
That nonzero status is intentional evidence that the full requested feature is
unfinished. Forced suite separately passes. Existing MC integration: 20 native runs
pass; immutable checkpoint/replay: 4 runs pass. Existing curve suites: 31 core,
7 regression and 14 remaining unsupported-configuration cases pass. Forcing/MC
were removed from the old rejection list because new positive suites exercise them;
no remaining guard assertion was weakened. Strict Python IO, 10 native geometry
reader checks and isolated production C++ geometry persistence checks pass.
Unchanged detector parity: 53 cases, exact masks and original \(10^{-12}\) tolerances;
encounter state checks pass. These are engineering results, not scientific sign-off.

## Independent review and retained failures

Three read-only reviews covered timing/numerics, persistence and minimality.
Applied: shared finite forcing-field check for curves-only and MC, bound resolved
native hydro settings, bounded existing-file replay reads, reused curve-enabled
lookup, removed an unused include and wrapped new lines to repository limits.
New detector/encounter math was not modified. New-file cpplint and diff checks pass.
Reviewers withdrew two stale/missed-code findings (post reset and consumer-disable
bypass); neither justified a production change.

Retain initial failed runs: old baseline rejections, exact wrapped/lifted roundoff
comparison mistake, suppressed checkpoint fixture, birth counting that omitted
completed-step full traces, incorrect manifest suffix and future-chunk generation,
and an invocation with colliding report/raw paths. Test corrections repair those
fixtures; physical settings and frozen tolerances remain unchanged. Final reports
are separate from every failed attempt. The earlier candidate executables were
rebuilt during development; historical logs retain their recorded hashes, but only
locked-source binaries in the final packet are supplied as exact-code evidence.

## Remaining required decisions and qualification

Two UserQuestions remain pending: leaf-cell versus fixed-physical detector width;
and whether to repair inherited MC regrid allocation using actual child fluid
masses, or retain native transport with MC AMR explicitly unqualified. In the
uniform-child oracle, native remap sends all \(4096\) IDs to one child, counts
\([0,0,0,4096]\), versus four equal masses. Preserving encounter state is not a
mass-following tracer qualification. Startup MC AMR rejection remains active.
See `curve-shocks-mc-amr.md` for the reproducer and repair constraints.

AMR detector qualification also needs an independent interface/physical oracle:
ShockIt AMR uses its own conservative halo construction, while this native sampler
uses AthenaK synchronized primitive ghosts. Uniform kernel parity does not prove
whole-AMR-mask parity. Do not claim otherwise.

RUSTY CPU/MPI execution, migration across ranks, node/particle counts, resource
budgets and performance remain pending. The current local MPI launcher fails
before launch with PMIx runtime 0x60100 versus PRRTE build 0x50009. MPI builds and
singleton collectives do not count as multi-rank execution. Combined multi-rank
MC+curve checkpoints currently reject before opening a native output: native MC
requires per-rank files while curves require shared fluid files. A compatible
layout must be implemented and qualified for the RUSTY coexistence restart gate.
Frontier/GPU/threaded/single-precision remain unqualified and guarded. Human
scientific review of the exact-code packet remains pending.

## Laptop reproduction and output

From this checkout, with its exact Kokkos submodule initialized:

```sh
cmake -S . -B build-curve-shock -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_ENABLE_SERIAL=ON -DAthena_ENABLE_MPI=OFF -DAthena_SHOCK_TESTS=ON
cmake --build build-curve-shock -j8
cmake -S . -B build-curve-shock-test -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_ENABLE_SERIAL=ON -DAthena_ENABLE_MPI=OFF -DPROBLEM=material_curve_test
cmake --build build-curve-shock-test -j8
python3 tst/shocks-curves/io_reader.py
python3 tst/shocks-curves/acceptance.py \
  --executable "$PWD/build-curve-shock/src/athena" \
  --prescribed-executable "$PWD/build-curve-shock-test/src/athena" \
  --workdir /private/tmp/curve-shock-check-UNIQUE
python3 tst/shocks-curves/acceptance.py --suite forced \
  --executable "$PWD/build-curve-shock/src/athena" \
  --workdir /private/tmp/curve-shock-forced-UNIQUE
```

The forced suite uses the already installed NumPy dependency of the native binary
reader; no new package was installed. Here Python was the existing ShockIt venv
`/Users/hitesh/hitesh/git/shockit/.venv/bin/python`. The core command intentionally
returns nonzero until the required coexisting AMR gate is implemented/qualified.
Choose new output directories; immutable histories/checkpoints reject collisions.
Input fixtures live in `inputs/tests/shocks-curves/`.

Enable `<material_curves>` normally and `<shock_tracking> enabled=true` with a safe
unique lineage. One global switch tracks both present consumers. MC files remain
under the configured history directory; curves use its `curves/` child. Events are
133 CSV fields, optional full/birth traces 33; strict reader:
`vis/python/curve_shock_history.py`. Complete records use observed entry/exit
brackets; unknown entry and truncated end remain censored. Never merge suspended
views into completed-event counts. Whole-curve statistics are an analysis step.

Restart with the native fluid `.rst`, its `.curves` geometry and
`.curves.shock.rankR` state plus the complete committed history prefix. Use the same
source/configuration/ABI/partition; a new basename avoids immutable checkpoint
collisions. Relative history paths work with native `-d`. Do not enable tracking
from an old untracked checkpoint, disable a recorded consumer, omit a prefix,
copy parent history into inserted nodes, or interpret checksums as authentication.
Detailed wire/persistence contract: `tst/shocks-curves/io-format.md`.


Final MPI build: full `material_curve_test` executable compiled with Homebrew
mpicxx and Serial Kokkos. Singleton MPI native curve core passes 31 cases;
singleton curve-history run repeats 19/20 gates, with the same explicit MC AMR
guard failure. This is **not** an inter-rank test. The two-rank launcher probe
exits 213 before AthenaK starts; exact PMIx/PRRTE error is retained.

Final compact evidence is committed under `evidence/curve-shocks/2026-10-17/`.
Raw test directories/binaries/configuration and historical failure artifacts are
archived outside Git in the exact-code bundle. The compact copied remap log trims
only its final blank line; the raw archive preserves original bytes. No scientific
sign-off, merge into feature_dump, push or original-branch modification is implied.


## Visual validation requirement and fresh rerun — 2026-10-19

User explicitly requires plots with test results, always. Include failures as well
as passes, and put explanations beside each plot. Preserve this presentation
requirement for subsequent work on this feature.

Re-ran the unchanged native core and forced fixtures using the same exact locked
production binaries: 36 core runs give 19/20 checks; the combined MC/curve AMR
startup guard remains the sole failed required gate. Three forced runs pass exact
on/off fluid/MC/node passivity and optional-trace equality. A freshly compiled
actual-method remap harness reproduces [0,0,0,4096] child occupancy. Its process
success reports reproduction, not physical acceptance.

Six reviewed figures show planar passage/strength evidence, independent newborn
histories, a native 32-cubed forced-box projection with both carrier populations,
per-node flags/episodes, exact passivity/restart comparisons, and the inherited
AMR failure alongside all 20 gates. Forced run: 28 complete curve records for
16 nodes; 7063 complete MC records for 4096 particles, plus 556 separately censored
terminal MC views. Different populations sample different spatial regions;
matching their counts is not an acceptance requirement. These are short-run
engineering checks, not converged turbulence statistics, hardware performance
benchmarks, or scientific approval.

Plot reader assertions check actual files, input/log/binary hashes, exact replay
evidence, no pre-birth records and all-zero displayed error curves. Every locked
production source hash and native binary hash was rechecked against the original
exact-code packet; production source is unchanged. Fixed a plot-provenance check
that initially skipped log validation when the runner omitted a directory field;
final plotting validates every explicit run-name directory and input/log hash.

Reviewed output:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-shock-interaction-plots-20261019-reviewed/index.html`.
Full PNGs include nearby explanations; original plot iteration is retained
separately. Raw reruns and exact commands are archived with the reviewed plots;
compact fresh reports and metrics are in `evidence/curve-shocks/2026-10-19/`.

Reproduce the figures, with the already-installed scientific Python environment:

```sh
/Users/hitesh/hitesh/git/shockit/.venv/bin/python tst/shocks-curves/visuals.py \
  --core /private/tmp/curve-shock-plots-core-20261019 \
  --forced /private/tmp/curve-shock-plots-forced-20261019 \
  --remap-report /private/tmp/curve-shock-plots-remap-20261019/report.json \
  --output /private/tmp/curve-shock-figures-UNIQUE
```

The plot script reuses acceptance.py and production readers. It adds no detector,
trajectory reconstruction, dependency, physical parameter, floor or tolerance.
The curve-history implementation remains on codex/curve-shock-histories; it is not
yet merged into the local feature_dump branch. Original shock branch is preserved.


## Why live MC + curve AMR remains blocked — 2026-10-20 clarification

Rechecked current branch505c0304: ValidateShockPhysics still rejects multilevel or
adaptive meshes when MC shock tracking is enabled. Curve-only AMR works through
its tested observation/birth/restart path. The blocking evidence is inherited MC
mass allocation at a regrid, not a loss of the stored EncounterState.

A mass-following MC population colocated at one coarse cell center should split
among refined child cells according to their fluid masses. The focused uniform
2D fixture has4096 distinct identities and four equal-mass children: expected
counts are1024 per child, with sampling fluctuations for a stochastic split.
Actual native remap produces[0,0,0,4096] because it locates the new containing cell,
snaps to its center, then resets PLASTMOVE/PLASTLEVEL. That bypasses redistribution
of this regridded population. Complete IDs/history survival alone cannot establish
correct gas-mass sampling; shock encounter statistics could inherit this bias.

The proposed repair is child-mass-weighted random selection, with stable IDs,
complete histories and reproducible restart/partition behavior. It must preserve
tracking on/off passivity and pass physical refinement/coarsening checks before
the guard is lifted. Physical allocation and AMR detector-width UserQuestions
remain pending. The explanatory chart reuses the actual2026-10-19 rerun; no new
simulation or numerical change was made for this clarification.


## Original MC behavior compared directly — 2026-10-20

User asked how vanilla MC handles this. Here vanilla means the actual local
pre-shock multiphase_turb branch526a187f938e38634e89a88192ebe9731c4d36bb;
no claim is made about another upstream repository/version.

There are three operations in original particles_lagrangian_mc.cpp:

1. PushLagrangianMC(line500): use positive outgoing saved mass flux divided by
   old-cell mass, draw a uniform pseudo-random number keyed by ID/cycle/seed,
   choose one face or remain. Particle positions encode cell-center membership.
   Hydro::SaveFlux stores stage-summed density flux scaled by timestep/cell width.
2. AdjustMeshRefinement(line592): after actual particle boundary transfer, compare
   prior/current levels. Coarse-to-fine crossing chooses a transverse fine entry
   subface using its positive saved flux fraction(two in2D,four in3D). Fine-to-coarse
   crossing applies deterministic coordinate corrections from stored parity/move.
3. RemapAfterMeshRefinement(line734): when the mesh itself changes, locate the new
   leaf from the stored position, snap to its cell center, set PGID, reset last move
   tozero and last level to new level, then migrate real/integer records to owning
   MPI ranks. No child-fluid-mass random redistribution occurs in this routine.

These are separate paths: existing stochastic fine-face entry selection does not
redistribute a coarse-cell population when the mesh is refined underneath it.
This distinction is the reason the integration gate remains blocked.

Direct check: extracted the entire original MC translation unit byte-for-byte
using git show526a187f, compiled it against current headers and existing native
method harness constructor/count stand-ins, and ran occupancy only. Compile0,
process0; child counts[0,0,0,4096], identical to the current history branch. This
proves reproduction of the original positional remap, not scientific acceptance
or a complete original-branch fluid/MPI executable qualification. No original
branch, production source, transport RNG or thresholds were edited.

Source snapshot, full compile command/output, executable and explained comparison
plot: /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/vanilla-mc-amr-20261020.
Compact probe report is retained in evidence/curve-shocks/2026-10-20/.


## Curve nodes versus MC particles — 2026-10-20 clarification

Material-curve geometry uses continuous velocity-following nodes satisfying
\(d\mathbf{X}(a,t)/dt=\mathbf{u}(\mathbf{X}(a,t),t)\), with persistent material
labels. Current MaterialCurves::Sample interpolates native primitive velocity
(bilinear in2D, trilinear in3D); Advance applies the existing two-stage RK2 update.
Curve nodes are not driven by MC cell hops. They do not need to be objects of the
existing particle class to implement continuous Lagrangian advection.

MC particles are an independent optional population following mass flux
statistically. Sharing the detector/encounter code and fluid state does not share
or replace the two motion algorithms. The combined MC/curve AMR qualification
blocker concerns the MC population's regrid allocation, not MC-based curve motion.
This clarification is based on current source inspection; no new simulation,
production change, or scientific qualification was performed.
