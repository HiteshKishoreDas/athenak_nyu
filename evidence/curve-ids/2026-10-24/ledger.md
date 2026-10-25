# Compact curve IDs: high 8 curve bits, low 24 creation bits

2026-10-24. Explicit user instruction: "just use the first 8bits for curve id and
rest 24bits for the sequencve". Worktree `.worktrees/shock-lagrangian`, branch
`codex/shock-lagrangian`, initially clean at `8ae4b2c9`. Original branches preserved.

## Accepted requirement and change

The current packed node identity is now

\[
\mathrm{ID}=2^{24}\,\mathrm{curveID}+\mathrm{creationSequence}.
\]

Curve IDs occupy the high 8 bits and remain one-based: 1–255. Creation sequences
occupy the low 24 bits and remain one-based: 1 through $2^{24}-1$. Zero remains a
missing-link sentinel. IDs never wrap or get reused; the allocator can hold the
exhausted sentinel $2^{24}$. Node IDs, prev/next links, parent links, curve identity
and allocator storage now use `uint32_t`. Counts, simulation generations, hashes,
periodic images and shared particle/curve encounter counters keep their necessary
existing widths. Shared encounter-state ID slots remain uint64; current curve
values are checked before any narrowing lookup.

New node insertion preserves existing identities and updates reciprocal links.
Birth, material label and history rules are unchanged. Exhausted IDs refuse an
insertion and retain the existing quality reason 2; a split batch exceeding the
remaining namespace is rejected atomically. Initial out-of-range curve requests
fail before conversion, including values that would otherwise alias after wrap.
ID-only MPI payloads now use `MPI_UINT32_T`; counts remain `MPI_UINT64_T`.

## Persistence and compatibility

- Current ASCII geometry/checkpoint schema: v3 with U32 identity/link columns.
  Historical v2 geometry remains Python-readable with its original IDs.
- Current curve thermo CSV schema: v2, same fields and filenames, explicit compact
  ID meaning. Python retains historical v1 reading and rejects mixed encodings.
- Curve shock contract begins `curve-node-history-v2`; older v1 contracts retain
  historical decoding. Shared event payload and encounter transition rules remain.
- Old curve checkpoints are explicitly rejected for restart rather than silently
  reassigned. Existing strict shock source/config guards remain in force.
- Particle THP v2 bytes and transport are unchanged. Compact IDs for curves 128–255
  exceed the signed THP tag range; future THP integration must handle unsigned
  identities explicitly. This slice does not add THP curve output or mesh fields.

These are versioned changes to existing file families, not another independent
codec or extension. Historical evidence is preserved. No restart conversion is
claimed. Human scientific sign-off is still pending.

## Reviewed source and tests

`review.md` records the independent persistence/simplicity review and applied fix.
`report.json` lists every final group, raw report path and SHA256. Final native
reports bind executable hashes and compiled source fingerprints. Their binary
hashes were verified against the completed builds. Source/build hashes and raw
inputs/logs remain retained in each run directory.

| Final suite | Passing groups |
| --- | ---: |
| Output/identity-remapped historical controls | 15 |
| Curve thermo | 18 |
| Curve shocks (includes repaired-checksum alias rejection) | 20 |
| Continuous tracers | 14 |
| Startup summaries | 9 |
| Curve reader/plots | 11 |
| Allocation/width/checkpoint boundary gates ID01–ID06 | 6 |
| Total | 93 |

Additional checks passed: 11 Python geometry-reader tests retaining historical
64-bit fixtures, standalone curve IO/adversarial tests, Serial C++ serialization/
parser roundtrip and malformed/pair rejection, detector self-check, encounter
state checks, Python compilation and diff checks. C++ roundtrip has stubbed physics
and audits; production native suites exercise actual sampling and structural audit.

ID01 forces a valid native allocator near exhaustion without creating millions
of nodes. The final valid ID is exactly $2^{32}-1$; links and checkpoint allocator
sentinel are preserved. ID02/03 cover insufficient split capacity and exhausted
counters without wrap/reuse. ID04 rejects curve IDs 0, 256 and 4294967297 and more
than 255 curves. ID05 rejects an out-of-width persistent link before assignment.
ID06 verifies explicit old-checkpoint rejection and unchanged old-analysis IDs.
Each mutated checkpoint used by a final boundary run is retained under that run's
`checkpoint_used/`, with its hash in the report. These are structural fault and
boundary fixtures, not physical initial-condition models.

TO01/02 compare physical records against saved historical controls after a checked
bijective ID encoding remap: no filtering, tolerance widening, field omission or
clock change. Particle THP bytes and planar conserved fluid tables remain exact.
TO03 has zero manufactured thermo error under the unchanged $2\times10^{-12}$
budget. TO21 still has 17 inside-born and 22 outside-born positive-time nodes;
cycle-13 restart exactly preserves all records and final geometry. Birth samples
in this fixture reside in per-step traces; separate birthtrace files are zero and
not newly qualified. No human scientific approval is inferred.

## Commands (from the worktree)

```sh
cmake -S . -B build-lagrangian -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-lagrangian -j8
cmake --build build-lagrangian-prescribed -j8
cmake -S . -B build-lagrangian-mpi -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build-lagrangian-mpi -j8
python3 tst/scripts/material_curves_io/test_reader.py
python3 tst/shocks-curves/io_reader.py
python3 tst/scripts/material_curves_io/test_cpp.py \
  --build build-lagrangian \
  --kokkos-library build-lagrangian/kokkos/core/src/libkokkoscore.a
python3 tst/tracer-output/check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs build-lagrangian/evidence/native-output-qualified-final-20261023 \
  --workdir build-lagrangian/evidence/curve-id32-reviewed-output-20261024
python3 tst/curve-thermo/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/curve-id32-reviewed-CT-20261024
python3 tst/shocks-curves/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/curve-id32-reviewed-CSH-20261024
python3 tst/shocks-lagrangian/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/curve-id32-reviewed-LT-20261024
python3 tst/tracer-output/id_check.py \
  --executable build-lagrangian-prescribed/src/athena \
  --historical-run build-lagrangian/evidence/existing-output-final-20261024/planar \
  --workdir build-lagrangian/evidence/curve-id32-reviewed-boundaries-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-setup/check.py \
  --executable build-lagrangian-prescribed/src/athena \
  --baseline /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-thermo-synergy-20261022/final-reviewed/runtime/build-lagrangian-prescribed/athena \
  --workdir build-lagrangian/evidence/curve-id32-reviewed-setup-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/curve_reader_check.py \
  --runs build-lagrangian/evidence/curve-id32-reviewed-output-20261024 \
  --three-d-run build-lagrangian/evidence/curve-id32-reviewed-CT-20261024/forced_both \
  --output build-lagrangian/evidence/curve-id32-reviewed-reader-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/visuals.py \
  --run build-lagrangian/evidence/curve-id32-reviewed-output-20261024 \
  --id-boundaries build-lagrangian/evidence/curve-id32-reviewed-boundaries-20261024 \
  --output build-lagrangian/evidence/curve-id32-final-reviewed-plots-20261024
./build-lagrangian/shock_detector_test --self-check
./build-lagrangian/shock_encounter_test
```

Use fresh destinations to reproduce. Reuse the recorded source and Kokkos revision
in `report.json`; no physical input values or tolerances were tuned.

## Retained development failures and limits

First CT10 fault injection replaced the historical header version and left the
new file unchanged. Corrected the fixture to mutate v2 and assert the mutation is
effective; all 18 groups then passed, and passed again on reviewed code.
Initial C++ isolated IO invocation lacked exported compiler commands; enabled
that existing CMake option. Its older harness then lacked stubs for the already
existing thermo sampling symbols; added explicit fail-if-called stubs, consistent
with its documented serialization-only purpose. Final C++ check passed.
Independent review found and reproduced the ambiguous substring contract marker;
replaced with an explicit leading version, retained the reproducer in IO checks.
The first bit-plot title incorrectly said neighbour links were unchanged; corrected
to consistent links. Existing IDs stay unchanged; connectivity changes on insertion.
All development run directories and logs remain retained.

Both MPI build attempts pass. The isolated MPI IO launch and a two-rank native
AMR/coexistence launch both fail before AthenaK starts with exit 213:
`PMIx Runtime: 0x60100; Build: 0x50009` (runtime/built PRRTE ABI mismatch).
This is a local toolchain/runtime limitation, not a passed MPI test. No MPI runtime,
GPU, RUSTY, Frontier, production throughput or memory-budget claim is made.
Fix the local PMIx/PRRTE installation or qualify on the intended cluster before
claiming cross-rank support. Fresh Serial/double engineering gates are complete;
human exact-code scientific review remains pending. No push performed.
