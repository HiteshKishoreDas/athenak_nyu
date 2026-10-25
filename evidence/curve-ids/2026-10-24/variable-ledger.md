# Configurable allocation of the 32 curve-node ID bits

2026-10-24. Explicit user request: make curve-number bits an input variable so a
single curve can use the entire 32-bit ID. Worktree `.worktrees/shock-lagrangian`,
branch `codex/shock-lagrangian`, initially clean at `d81ee29e`. Other branches
and prior evidence remain preserved.

## Interface and implementation

```ini
<material_curves>
curve_id_bits = 0
ncurves = 1

<curve1>
id = 1
```

Retain the existing curve geometry/settings in that block. Omitting curve_id_bits
keeps the prior 8/24 default. The option is a checked unsigned integer in 0–30.
For $B>0$, the identity is

\[
\mathrm{ID}=2^{32-B}\,\mathrm{curveID}+\mathrm{creationSequence}.
\]

Curve IDs occupy the high B bits and range 1 through $2^B-1$. Sequences range 1
through $2^{32-B}-1$. With B=0 there is only logical curve ID 1, retained in
metadata; the entire node ID is the sequence, from 1 through $2^{32}-1$. Zero
remains the missing-link sentinel. Values 31/32 cannot represent a valid curve's
minimum two nodes and are rejected. Existing max_nodes (default 1000000), bounded
MPI batch counts and the 100000-curve IO cap still apply. ID capacity is not a
qualification of billions of resident nodes or cluster throughput.

Node IDs and neighbour/parent links stay uint32. Only the per-curve allocator
counter widens to uint64 so it can represent the exhausted $2^{32}$ sentinel.
Encoding/decoding special-case zero bits and compute masks in a wide type, avoiding
undefined shifts by the width of uint32. Allocation never wraps, reuses IDs or
renumbers existing nodes. Namespace, parent/birth/history, AMR and event rules
are unchanged. The startup summary prints the resolved curve and sequence bits.

## Files, readers and restart

Current geometry/checkpoint ASCII v4 adds an `id_bits B` row; next_sequence is
explicitly u64. IDs/links remain u32. The C++ geometry reader also accepts v3 as
fixed 8/24, subject to normal compatibility checks. Historical Python geometry
v2/v3 retain their original identities. Current curve thermo CSV v3 adds B to its
header, retaining v1/v2 analysis. The current curve shock contract starts with
`curve-node-history-v3 B`, retaining historical v1/v2 analysis semantics.

read_run carries the declared layout, including empty streams, and rejects
contradictions across geometry/thermo/journals. Within-stream thermo layout changes
also reject. The plot CLI decodes from this metadata. An isolated 32-bit trace
requires its manifest or an explicit `curve_id_bits=B` argument; no numeric
inference is used to guess a new layout. Missing/invalid metadata fail explicitly.

The split is immutable on restart. B=0 active-AMR checkpoint continuation is exact,
and changing it to B=8 rejects before applying curve state. Geometry v3/default8
compatibility is demonstrated with a checksummed same-source synthetic sidecar;
this does not bypass tracked restart exact-source/config/ABI guards for arbitrary
old checkpoints. No automatic ID conversion or checkpoint migration is claimed.
Particle THP bytes and output transport are unchanged. This increment adds no THP
curve/event writer, mesh shock fields, extra codec or dependency.

## Evidence and acceptance

`variable-report.json` lists every final group, raw report path/hash, additional
checks and limitations. Native reports retain exact commands, inputs, logs,
executable/build-cache hashes and compiled source fingerprints. Every final native
executable digest was verified against the completed builds. Failed/rejected
fixtures and checkpoint mutations remain in their run directories. Independent
persistence/simplicity review and applied findings are in `variable-review.md`.

| Final suite | Passing groups |
| --- | ---: |
| Existing output/control comparisons | 15 |
| Shared curve thermo | 18 |
| Curve shock histories | 20 |
| Continuous tracers | 14 |
| Startup summaries | 9 |
| Curve reader/plots | 11 |
| Fixed-eight-bit ID boundaries | 6 |
| Zero-bit/full-sequence boundaries | 6 |
| Variable-layout cases VB01–VB08 | 8 |
| Total | 107 |

Additional passes: 11 historical/current Python geometry-reader tests; standalone
curve IO including missing/invalid layout, explicit arguments, misleading historical
contracts and empty journals; Serial C++ IO save/load, malformed/pair rejection;
detector self-check; encounter state checks; Python compilation and diff checks.
Line-length filtering leaves only two inherited findings (driver output selector
and sampler allocation); full-repo style compliance is not claimed.

VB01 exercises B=0,1,4,8,16,30 with actual saved geometry/thermo and empty event
journals. VB02 rejects negative/text/fractional/overflowing/out-of-range inputs
and multiple/wrong-ID curves under B=0. VB03 compares the zero-bit translating
planar shock against the default run: 1728 curve samples and 16 emitted curve
events, exact after the checked ID remap, with regular particle records unchanged.
VB04 compares zero-bit live AMR: 9644 curve samples and 63 emitted events, with all
records and geometry exact after ID remap. VB05 restarts that active run at cycle
13; all histories and event/view identities match, and a split override rejects.
VB06 rejects inconsistent geometry/thermo, empty-journal and within-thermo-stream
layouts. VB07 verifies compatible default8 geometry v3 reading/restart. VB08
rejects a repaired-checksum B=31 checkpoint and invalid thermo-header bit metadata.

The boundary suites allocate exactly 4294967295 in both default8 and B=0 cases.
The zero-bit counter then persists 4294967296. Exhausted allocators and split
batches larger than their remaining namespace refuse insertion without wrap/reuse.
These are synthetic near-limit counter fixtures on accepted uniform flow, not runs
containing billions of nodes. Mutated checkpoints are retained in checkpoint_used/.

Existing CT/CSH physical oracles and their fixed $2\times10^{-12}$ budgets remain
unchanged. No tolerance was widened. Inside-born nodes retain unknown entry and
fresh histories; no parent encounters or exact crossing times are invented.
Human scientific sign-off remains pending.

## Reproduction commands

From `/Users/hitesh/hitesh/git/athenak_nyu/.worktrees/shock-lagrangian`, use fresh
output destinations when repeating these retained runs:

```sh
cmake --build build-lagrangian -j8
cmake --build build-lagrangian-prescribed -j8
cmake --build build-lagrangian-mpi -j8
python3 tst/tracer-output/check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs build-lagrangian/evidence/native-output-qualified-final-20261023 \
  --workdir build-lagrangian/evidence/variable-id-qualified-output-20261024
python3 tst/curve-thermo/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/variable-id-qualified-CT-20261024
python3 tst/shocks-curves/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/variable-id-qualified-CSH-20261024
python3 tst/shocks-lagrangian/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/variable-id-qualified-LT-20261024
python3 tst/tracer-output/id_check.py \
  --executable build-lagrangian-prescribed/src/athena --curve-bits 0 \
  --historical-run build-lagrangian/evidence/existing-output-final-20261024/planar \
  --workdir build-lagrangian/evidence/variable-id-final-boundary0-20261024
python3 tst/tracer-output/id_check.py \
  --executable build-lagrangian-prescribed/src/athena --curve-bits 8 \
  --historical-run build-lagrangian/evidence/existing-output-final-20261024/planar \
  --workdir build-lagrangian/evidence/variable-id-final-boundary8-20261024
python3 tst/tracer-output/variable_id_check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs build-lagrangian/evidence/variable-id-qualified-output-20261024 \
  --workdir build-lagrangian/evidence/variable-id-final-cases-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-setup/check.py \
  --executable build-lagrangian-prescribed/src/athena \
  --baseline /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-thermo-synergy-20261022/final-reviewed/runtime/build-lagrangian-prescribed/athena \
  --workdir build-lagrangian/evidence/variable-id-qualified-setup-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-output/curve_reader_check.py \
  --runs build-lagrangian/evidence/variable-id-qualified-output-20261024 \
  --three-d-run build-lagrangian/evidence/variable-id-qualified-CT-20261024/forced_both \
  --output build-lagrangian/evidence/variable-id-qualified-reader-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-output/visuals.py \
  --run build-lagrangian/evidence/variable-id-qualified-output-20261024 \
  --id-boundaries build-lagrangian/evidence/variable-id-final-boundary0-20261024 \
  --variable-layout-run build-lagrangian/evidence/variable-id-final-cases-20261024 \
  --output build-lagrangian/evidence/variable-id-final-plots-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python vis/python/plot_curves.py \
  build-lagrangian/evidence/variable-id-final-cases-20261024/zero_amr_shock \
  --node-id 33 --field pressure \
  --output build-lagrangian/evidence/variable-id-final-zero-node33-20261024.png
python3 tst/scripts/material_curves_io/test_reader.py
python3 tst/shocks-curves/io_reader.py
python3 tst/scripts/material_curves_io/test_cpp.py \
  --build build-lagrangian --kokkos-library build-lagrangian/kokkos/core/src/libkokkoscore.a
./build-lagrangian/shock_detector_test --self-check
./build-lagrangian/shock_encounter_test
```

## Retained findings and remaining gates

Initial native development suites passed. Independent review found two strict
reader bugs before final handoff: absent trace metadata could mislabel B=4 IDs,
and empty direct journals skipped layout validation. Both fixed with regressions;
full source-bound suites passed. Strict bit parsing avoided atoi aliasing before
initial builds. Line wrapping was followed by fresh builds and final native runs
so source fingerprints bind the final code. Earlier reports/figures remain retained.

Current Serial and MPI builds pass. Only Serial/double runtime is claimed. The
known local PMIx/PRRTE mismatch from the preceding compact-ID task was not retried;
MPI runtime, RUSTY, Frontier, GPU/threading/single precision, production budgets
and human scientific approval remain pending. No environmental repair, new physics
or deployment was attempted. No push performed.
