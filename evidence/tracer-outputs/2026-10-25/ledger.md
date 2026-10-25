# Existing-output integration

2026-10-25. User resumed output work after signed curve IDs. Starting branch
codex/shock-lagrangian, clean at a051de51. Existing worktree retained; original
branches preserved. Ponytail full and running-sim workflows apply. Existing
accepted physical fixtures/parameters are reused; no new scientific choices.

Accepted output targets: curve thermo and sampled particle/curve shock properties
through THP; mesh shock variables through the native mesh BIN writer. No separate
codec, container or new format selector. Full encounter journals and geometry
remain separate lossless products.

## Frozen first interface

U1: curve thermo uses literal THP v2, including its four signed integer slots
(cycle, node_id, curve_id, gid), standard time/x columns and named real fields.
Birth, immutable material label and sampled refinement level are additional
named fields. IDs never pass through Real. Existing particle THP bytes stay
unchanged. Retain current cycle/rank chunking and CommitFile immutable/replay
semantics; this work does not claim fewer files or append-performance gains.
A curve filename declares cycle, rank/rank-count and curve-bit split, including
empty chunks, using basename.output_id.cycleC.rankRofN.bitsB.thp. The existing
reader normalizes these records to its current curve analysis interface and
retains historical CSV reading. No user conversion step is required.

U2: optional named shock fields use the same native detector evidence and output
clock. Preserve unreduced region membership; distinguish pressure/temperature
jump estimates from fluid Mach. Missing estimates remain explicitly unavailable.
No episode identity/sequence is converted to floating point. Full every-step
events and restart state remain authoritative.

U3: variables=shock uses the existing mesh binary writer and derived-variable
registry. Reuse the synchronized particle detector array when available; the
curve-only periodic case uses the same DetectCell contract at output time.
No new time-integration/regrid ordering, detector criteria or physics support.

## Gates

- [x] U1 native THP curve round trips, exact controls and immutable replay (42 groups).
- [x] U2 sampled fields equal raw traces, synchronized time/AMR/birth tests.
- [x] U3 native BIN mesh fields agree with independent detector/particle samples.
- [x] Historical-reader, empty populations, signed IDs, cadence and restart gates.
- [x] Independent persistence/timing/simplicity review; apply concrete findings.
- [x] Actual explained plots, exact source/configuration and commands retained.
- [ ] Coherent commits and sealed evidence; no push requested.
- [ ] MPI runtime/GPU/cluster performance and human scientific review pending.

## U1 completed

Serial/double builtin and prescribed builds pass. TO15 + CT18 + VB9 =42 groups
pass at build-lagrangian/evidence/thp-U1-{TO,CT,VB}-20261025. Complete physical
observations and events match controls; existing particle THP bytes are exact.
Signed negative IDs, zero-bit live AMR, active restart, empty populations, unchanged
cadence, corrupt-header/truncation rejection and immutable conflicting replay are
covered. Tests that targeted current CSV bytes now target native THP bytes; old
CSV decoding is retained and used by the independent stored controls. No tolerance
or physical parameter changed. Reuse of standard headers adds no full-data copy
to regular particle writes. Binary chunks retain the existing per-cycle/per-rank
file count; no throughput or file-count improvement is claimed.

## U2/U3 implementation and review

Six optional canonical shock fields are shared by THP sampling and variables=shock
mesh BIN output. These are raw detector membership/validity and pressure/temperature
jump estimates plus validity flags. Fluid Mach and all detector criteria remain
unchanged. Unknown estimates stay NaN; no event sequence/identity is converted to
Real. Requesting shock properties requires tracking for each selected population.
Seeding target filters reject shock diagnostics explicitly; the shared fluid
field evaluator rejects accidental use without synchronized evidence.

Mesh output reuses the native writer/precision selection. It reuses the particle
cache when present, preserving physical-boundary validity masking; qualified
periodic curve-only runs evaluate the same DetectCell. Curve thermo requests
native observations without advancing encounter state. Full event journals remain
CSV and geometry remains ASCII; this task adds no independent codec.

Independent timing review found that a zero-step restart had a restored generation
but no restored detector cache. The unchanged detector kernel was extracted into
RefreshShockSignals, called after restart primitive initialization without touching
encounters/counters/journals. A new native regression compares all6 mesh fields and
particle/curve THP values at cycle13 with zero timesteps taken. The original
advancing restart test did not cover this case. The pre-fix defect was identified
from source review; no native pre-fix zero-step failure run was recorded.

Persistence review reproduced empty-block clock validation gaps and inconsistent
partial validity fields. Shared THP reading now validates finite/nonnegative block
time and nonnegative cycle before records, rejects conflicting same-cycle clocks,
and checks every selected validity prerequisite. Generic THP identical NaN replay
is retained. One immutable curve chunk requires exactly one THP block. Rank
completeness checks use bounded observed sets, never allocate a set spanning an
untrusted declared rank count. The 21-case standard-library reader regression
passes and was independently rerun by the reviewer.

The first new mesh parity harness wrongly tuple-indexed the native reader's list
of arrays. Corrected the harness to index its actual structure; failed report and
inputs remain at thp-shock-dev-20261025. No physical tolerance changed. The reviewer
also corrected mesh=False controls to dcycle0/dt-1, so diagnostics are fully off.
Reviewed-code SHOUT01–07/THP01 all8 groups pass: planar/negative-ID live-AMR exact
THP/BIN/raw-trace parity, unchanged fluid/geometry/full histories, curve-only output,
active and zero-step restart, configuration guards, native float32 exact casts,
empty populations and literal THP v2 integer-slot reading. Full regression gates
are being run before the final commit/evidence seal.

## Final engineering gate

Final source-bound builds: builtin Serial/double, prescribed Serial/double and
prescribed MPI compile all pass. The existing implementation fingerprint now
also includes derived_variables.cpp and binary.cpp. Native suites were rerun
on those final binaries. Kokkos and the primary feature_dump checkout remain
clean; the original codex/shock-multiphase branch remains preserved.

119 check groups pass: TO15, CT18, CSH20, LT14, VB9, startup9, reader11,
zero-bit boundaries8, eight-bit boundaries6, and new native output checks9.
Additionally, 21 THP clock/validity/corruption cases, 12 geometry-reader checks,
strict curve journal/trace IO, isolated Serial C++ IO, detector and encounter
state checks, Python compilation and diff checks pass. Independent timing and
persistence/simplicity reviews are recorded in review.md; all concrete findings
were fixed. Exact commands/reproduction are in reproduce.md and raw run reports.
No physical tolerance or initial-condition value was changed.

SHOUT08 repeats the accepted short 32³ forced box with the new fields. It compares
actual THP/BIN properties to raw traces on every saved clock and checks unchanged
fluid arrays, curve geometry and complete histories. The passivity control omits
new fields, mesh output and optional traces. The short fixture qualifies engineering
coexistence; it does not establish converged turbulent statistics.

Eighteen fresh explained figures cover native mesh/particle/curve samples, exact
property parity, zero-step cache restoration, measured file sizes, forced-box
occupancy, standard physical oracles/restart equality, ID limits and curve-reader
modes. Matched thermo sizes: planar CSV 348342 bytes versus THP 237492, 108 files
in each; AMR CSV 2014134 versus THP 1258894, 162 files in each. Records match exactly.
These measurements are not a throughput/cluster-memory benchmark.

Full encounter journals and geometry remain separate and lossless. The new
sampled shock properties follow the configured thermo cadence; no sequence/ID
is converted to Real or inferred from sparse output. MPI runtime, GPU/RUSTY/
Frontier performance budgets and human scientific sign-off remain pending.
No push requested or performed.

## Sealed evidence

Implementation source commit: `89bee32c2e2b65004c056ad5b3d0cd8ab60b5b02`.
The native-thp-bin-outputs-20261025 gallery contains 18 explained figures,
source/Kokkos archives, three binaries/configurations, reports and retained
malformed reader fixtures. All 43,628 raw-file hashes and 190 native-run
executable digests were verified. The 119 local check groups and 21 reader
cases are recorded separately. Manifest SHA256:
`a2af29a5a9b9028369cc92f131dde8d644c808fef67b130a0f13f1fe1ddcfe32`.
The manifest digest was rechecked after packaging; artifact.json records the
absolute gallery path. This is an engineering evidence seal, with the deployment
and scientific gates above still pending.
