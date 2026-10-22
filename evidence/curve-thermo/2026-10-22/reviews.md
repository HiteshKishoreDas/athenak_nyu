# Independent reviews — 2026-10-22

Two read-only reviewers used the requested ponytail-review instructions. Shared
interfaces frozen: one existing thermo output class; source=particles|curves|both;
unchanged particle v2; per-rank immutable curve v1; shared field evaluator and
current leaf-owner routing. Reviewers edited no files.

## review_curve_thermo_persistence

Reviewed typed ID return path, field-scaled counts/displacements, rank ownership,
node metadata and births, migration counters, immutable commits, restart and
reader. Re-read ten native evidence streams, including IDs above2^53. No actionable
bug found. Thermo returns before geometry counters, preserving curve state. No
MPI/GPU runtime qualification inferred. Minimalism verdict: lean; no new abstraction
recommended.

## review_curve_thermo_timing

Confirmed native float-threshold comparisons, restart output ordering (thermo
before checkpoint), initialization skipping output on restart, and CT16's genuine
off-cadence failure. Found P2: a repeated Finalize at an overdue output threshold
could advance last_time twice. Suggested distinguish repeated cycle or consume
all overdue thresholds. Chose persisted last_cycle so native catch-up semantics
remain unchanged. Reviewer rechecked revised new_cycle guard and zero-step/stale
parameter cases; no remaining actionable supported-path issue.

Reviewer also found CT17's incorrect *final.rst assumption; fixed to numbered
checkpoints. Actual retained pre-fix native .00001/.00002 clocks .09/.095 prove the
numerical defect separately from the test's StopIteration. Final native CT17 checks
equal.09 and zero-step restart. Added CT18's empty populations, addressing the
review's bounded missing-check observation. Older checkpoints lacking last_cycle
remain an explicit zero-step overdue-clock compatibility limitation.

These are engineering reviews. Human scientific sign-off remains pending.
