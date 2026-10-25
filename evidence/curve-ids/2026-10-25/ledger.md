# Signed32 curve identities

Date: 2026-10-25. User request: use signed IDs for curves too, consistent with
existing THP integer fields. Starting source: `ecd36eff` on `codex/shock-lagrangian`
in `.worktrees/shock-lagrangian`. Existing branch/dirty state inspected: clean.
Original `codex/shock-multiphase` and main `feature_dump` checkout are preserved.

## Decision and implementation

Curve numbers, node IDs, neighbour links and insertion-parent IDs now use
`int32_t`. The requested configurable split still uses all 32 bits. With packed
bits $u$, the signed ID is

\[
s(u)=\begin{cases}u,&u<2^{31},\\u-2^{32},&u\ge2^{31}.\end{cases}
\]

Thus negative IDs are valid, including $-2^{31}$ and $-1$; zero alone is the
missing-link sentinel. Logical curve numbers remain positive and within the
configured namespace. The default remains eight curve bits. Zero curve bits
still permits one logical curve ID1 and $2^{32}-1$ node creations. Existing
max_nodes, MPI batch and IO limits still apply. Capacity does not imply billions
of resident nodes have been tested.

Unsigned bit arithmetic is confined to packing/decoding and creation-order
checks. Signed shifts, signed overflow and implementation-defined unsigned to
signed overflow casts are avoided. The allocator remains uint64 so the exhausted
$2^{32}$ sentinel persists. Parenthood compares the underlying creation order,
not signed numeric order. MPI sampling identity payloads use MPI_INT32_T.
The shared encounter algorithm/state ABI remains unchanged: its uint64 ID slot
holds the canonical uint32 bit pattern. Curve-facing events/traces use the signed
ID, and checkpoint validation maps the bits explicitly back to geometry.

Files remain in their existing families. Geometry/checkpoint ASCII v5 declares
I32 identities/links; thermo CSV v4 and the curve-node-history-v4 contract declare
signed identity semantics with the bit split. The readers retain old geometry
v2/v3/v4, thermo v1/v2/v3 and event versions with their original unsigned identity
values. Current reader id_layout is (32,B,True); prior unsigned32 is (32,B).
Mixed interpretations reject, including empty streams. Default plot-node selection
uses creation order across the sign transition. Standalone signed traces require
curve_id_bits and signed=True unless their manifest is present.

Compatible geometry-only v3 restarts are checked after constructing an explicit
unsigned historical sidecar. This does not promise cross-code tracked restart
migration: exact source/contract/configuration guards remain. New signed live-AMR
checkpoint continuation is exact. The bit split is still immutable on restart.
Regular particle THP v2 and its reader interpretation are unchanged. This change
prepares signed curve identities; it does not add curve THP writers, sampled
particle shock fields or mesh shock BIN variables.

## Acceptance and retained failures

`report.json` records the nine final suite paths, SHA256 digests, check IDs and
results: **110 passing local Serial/double check groups**. Every reported native
executable digest was checked against the final built executable. Raw run folders
retain exact commands, input files, run logs, rejected cases, mutated checkpoints,
compiled source fingerprints and build-cache hashes.

- Output integration 15; curve thermo18; curve shocks20; continuous tracers14;
  startup9; reader/plots11; zero-bit boundaries8; eight-bit boundaries6;
  configurable layouts9.
- VB01 covers actual signed IDs for B=0,1,4,8,16,30, including empty journals.
- Existing physical/oracle and restart tests retain their unchanged tolerances.
  No timing or physical crossing time is manufactured.
- VB03/VB04 compare zero-bit planar/live-AMR observations and complete histories
  exactly against default-eight-bit runs after the explicit ID map.
- VB09 repeats live AMR as curve255 with all negative node IDs: 9644 curve
  samples, 63 emitted events and active cycle13 restart continue exactly after
  mapping identity. Regular particle histories match without an ID map.
- ID08 also constructs a same-source historical unsigned v4/B0 geometry sidecar
  containing UINT32_MAX and next_sequence=2^32, then resumes it natively. The
  signed -1 identity and exhausted counter are preserved. This optional review
  suggestion is now covered.
- ID01 allocates the final bit pattern as -1 and persists counter4294967296 at
  B=0. ID07 creates -2147483648 with positive parents. Production geometry and
  ancestry/link audits accept both. Exhaustion and atomic batch refusal remain.
  These synthetic near-limit counter fixtures do not contain billions of nodes.
- Twelve standard-library geometry-reader checks cover historical formats,
  signed sign-transition ancestry and width rejection. Strict encounter/trace
  reader checks exercise negative identities, unsigned compatibility, exact
  integer range/zero rejection and unchanged particle-reader semantics.
- Isolated Serial C++ IO save/load and malformed/pair rejection pass. That harness
  stubs physics/audits; native suites exercise the production code paths.
- Detector self-check: `8 device numerical/support self-checks PASS`.
  Encounter test: `EVT01-08, EVT11-14 state checks PASS; seed72941; large typed IDs`.
  Python compilation and diff checks pass. Changed C++ source files have no lines
  over100 columns; full-repository style compliance is not claimed.

Independent persistence/simplicity review found that SaveShockRestart temporarily
stored signed map keys in vector<uint64_t>, requiring an implementation-defined
narrowing conversion on lookup. Replaced the temporary vector and loop with
NodeID. This changes no state ABI. Final builds and all native suites were rerun
after the fix; the initial passes are retained as pre-review evidence.

Retained initial failures: CT10 and the isolated C++ harness targeted old version
header strings, making their fault injections ineffective; CR08 expected obsolete
uint64 wording. Fixtures/asserted wording were corrected, preserving equivalent
corruptions and all tolerances. Fresh final suites passed. Earlier development
reports and logs remain; none is relabeled as final evidence. No physical run
failed and no physics parameter was changed to obtain a pass.

## Reproduction

Use the worktree above and fresh output destinations when repeating. Local
builds retain their configured builtin/material_curve_test problem generators,
Serial/double backend and compiler/cache metadata. Kokkos is
`08ceff92bcf3a828844480bc1e6137eb74028517`. MPI build compiles; runtime was not
retried because prior PMIx runtime0x60100 versus PRRTE build0x50009 remains an
unresolved environment conflict. No MPI runtime/GPU/RUSTY/Frontier qualification
or production memory/I/O budget is claimed. Scientific sign-off remains pending.

```sh
cmake --build build-lagrangian -j8
cmake --build build-lagrangian-prescribed -j8
cmake --build build-lagrangian-mpi -j8
python3 tst/tracer-output/check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs build-lagrangian/evidence/native-output-qualified-final-20261023 \
  --workdir build-lagrangian/evidence/signed-id-reviewed-output-20261025
python3 tst/curve-thermo/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/signed-id-reviewed-CT-20261025
python3 tst/shocks-curves/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/signed-id-reviewed-CSH-20261025
python3 tst/shocks-lagrangian/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/signed-id-reviewed-LT-20261025
python3 tst/tracer-output/id_check.py \
  --executable build-lagrangian-prescribed/src/athena --curve-bits 0 \
  --historical-run build-lagrangian/evidence/existing-output-final-20261024/planar \
  --workdir build-lagrangian/evidence/signed-id-sealed-boundary0-20261025
python3 tst/tracer-output/id_check.py \
  --executable build-lagrangian-prescribed/src/athena --curve-bits 8 \
  --historical-run build-lagrangian/evidence/existing-output-final-20261024/planar \
  --workdir build-lagrangian/evidence/signed-id-sealed-boundary8-20261025
python3 tst/tracer-output/variable_id_check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs build-lagrangian/evidence/signed-id-reviewed-output-20261025 \
  --workdir build-lagrangian/evidence/signed-id-reviewed-cases-20261025
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-setup/check.py \
  --executable build-lagrangian-prescribed/src/athena \
  --baseline /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-thermo-synergy-20261022/final-reviewed/runtime/build-lagrangian-prescribed/athena \
  --workdir build-lagrangian/evidence/signed-id-reviewed-setup-20261025
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-output/curve_reader_check.py \
  --runs build-lagrangian/evidence/signed-id-reviewed-output-20261025 \
  --three-d-run build-lagrangian/evidence/signed-id-reviewed-CT-20261025/forced_both \
  --output build-lagrangian/evidence/signed-id-reviewed-reader-20261025
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-output/visuals.py \
  --run build-lagrangian/evidence/signed-id-reviewed-output-20261025 \
  --id-boundaries build-lagrangian/evidence/signed-id-sealed-boundary0-20261025 \
  --variable-layout-run build-lagrangian/evidence/signed-id-reviewed-cases-20261025 \
  --output build-lagrangian/evidence/signed-id-reviewed-plots-20261025
python3 tst/scripts/material_curves_io/test_reader.py
python3 tst/shocks-curves/io_reader.py
python3 tst/scripts/material_curves_io/test_cpp.py \
  --build build-lagrangian --kokkos-library build-lagrangian/kokkos/core/src/libkokkoscore.a
./build-lagrangian/shock_detector_test --self-check
./build-lagrangian/shock_encounter_test
```

The gallery adds CLI pressure plots for actual B=0 node33 and signed curve255
node-16777183 from the final cases. Both begin at their saved birth time
0.009990907573667196, have149 samples and one emitted partial encounter; unknown
entry is retained. Descriptions accompany all figures. No new physical initial
conditions or parameter choices were introduced.

## Next gate

Independent signed persistence/simplicity review is recorded in review.md.
Seal exact source, binaries/cache/configuration, raw hashes and plots after
reviewed-code validation completes. No push requested. Separate next scope: reuse THP for curve
thermo, followed by sampled shock fields and mesh BIN diagnostics.
