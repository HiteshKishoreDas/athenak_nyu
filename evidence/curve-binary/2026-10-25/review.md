# Independent reviews — unified CRV

Two existing reviewers continued the standing timing/persistence review gate,
read-only and using ponytail-review. No overlapping native builds or source edits.

## Timing / interaction

- Missing crv registration in mesh-variable exclusions was found and fixed in
  outputs.cpp and basetype_output.cpp; explicit file_type=crv now runs without
  supplying a fictitious mesh variable.
- Sample(fields) and Sample(observations) return before advection-only sample and
  migration counters. Geometry diagnostics are pure and shared with ASCII.
- Completed-step adaptation/shock observation precedes diagnostics. Existing
  AdvanceTracerClock preserves duplicate/final/dt restart semantics. MPI collective
  order is unchanged in the shared sampler; no runtime MPI claim is made.
- Automatic ASCII suppression is confined to unified-output mode. Explicit ASCII
  and encounter journals remain available. Unique output ids avoid path collisions.
- Final native follow-up only corrected the startup thermo-enabled summary.
- Plotting caught the test's false active-cycle13 assumption. The final test derives
  actual first inside cycle37 and asserts phase=active. Reviewer also requested
  explicit zero-step geometry comparison; final CRV08 compares both zero-step and
  continued geometry, in addition to all node/particle sample records.
- Reviewer reports no production blocker; native runtime evidence is coordinator's.

## Persistence / reader

- Reproduced cross-clock domain mixing with a repaired checksum; fixed immutable
  stream/curve metadata, counter monotonicity and node-set retention checks.
- Fixed unrelated CRV stream shadowing a specifically selected historical THP/CSV
  stream. Rank decomposition may change across complete clocks for analysis;
  consistency still required within each rank set.
- Reviewed scalar widths/order, checksum coverage, signed IDs, empty rank metadata
  and immutable replay; no remaining writer or persistence blocker.
- Independently ran 28 reader groups successfully. Final missing-rank fixture uses
  a valid rank1of2 home ownership and asserts the intended missing-rank guard.
  Separate fixture preserves signed -1 and exhausted full-width allocator exactly.
- Independent report retained from /private/tmp/crv-reader-final28-20261025/report.json.

Ponytail assessment: no useful additional abstraction or dependency needed.
