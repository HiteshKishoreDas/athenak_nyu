# Sampled shock output versus complete encounters — 2026-10-26

## Accepted correction

User clarified that sampled CRV/THP shock properties must not imply every-step
shock finding. Evaluate sampled properties only at those outputs' saved clocks.
A separately declared complete-history output block explicitly enables every-step
encounter tracking for its selected carrier and writes complete records in .crv
or .thp. User accepted event-only binary output blocks and replacement of completed CSV
journals. Sampled properties and complete episodes are separate explicit products.

Starting branch codex/shock-lagrangian, clean at4e797769. Preserve source branch
codex/shock-multiphase and existing unrelated work. Ponytail remains active.

## Inspection

Current code couples detector configuration, per-step observation, encounter state,
full-mesh particle cache and journal/checkpoint guards through shock_enabled.
CRV BinarySnapshot demands the encounter generation/time be current. Thus output
cadence alone currently cannot suppress timestep detection. Separate configuration,
sampled-signal freshness and encounter generation are needed; do not merely skip
ObserveShocks or weaken stale-generation guards. Curves can reuse native touched
stencils; particle sampling must avoid refreshing a full-mesh cache at every step
when only sampled output is requested. Encounter updates must never occur as a
side effect of snapshot output or a repeated same-cycle Finalize call.

## Pending gates

- [x] Check branch/dirty state and trace native output/driver/configuration coupling.
- [x] Resolve binary history content and CSV replacement decisions.
- [x] Freeze minimal input/state/persistence interface, retaining typed identities.
- [x] Independent sampled/full modes, opt-in carrier scope and cadence tests.
- [x] No lost/duplicate encounters across tested buffering, migration and restart cases.
- [x] Independent timing/persistence review; actual plots and source-bound evidence.

No source behavior changed yet. No new tests claimed. Hardware/scientific gates
remain separately pending.

## Implemented interface and review gates

Default mode=snapshot preserves CRV/sample-THP schemas and evaluates requested
shock fields only at saved clocks. Explicit mode=encounters enables one carrier's
complete tracking and typed event-only .crv/.thp batches. Native negative dcycle
initial/final-only output is supported; zero disables. Global enabled is detector
configuration, not implicit history activation. Optional qualification traces stay
CSV; completed-event CSV is gone. Sampled MC properties are point samples only;
complete MC histories remain rejected. No source thresholds/physics changed.

Shared bounded HistoryStore reuses immutable commits. Sorted events keep original
emission generations and typed IDs/sequences. Unflushed records persist inside
the checksummed fluid-bound paired sidecars. Particle sidecar v3 and curve sidecar
v2 explicitly reject old incompatible tracked restart files. Full-prefix links,
hashes and orphan-tail guards prevent missing/modified/branched streams.

Independent review found terminal views were written before final flush, producing
non-identical zero-step replay. Finalize now flushes event outputs, then writes
views, then the last-position paired restart. Review also found ordinary and final
checkpoints can share a clock: test selection now chooses the terminal file for
terminal replay, and separately verifies ordinary checkpoints with their exact
clean prefix. A non-prefix tail fails before resumed integration and is not deleted.

Reader review preserved committed domain exits carrying both right-censored and
domain-exit flags; they must not be discarded as merely terminal views. Native
protocol tests preserve an already-destroyed identity, invalid-fluid endpoint
clock, signed -1 and sequence above2^53 exactly; these are unit protocol evidence,
not live outflow qualification. Known-magics dispatch supports custom directories
sharing snapshot/event extensions, while unknown/corrupt candidates stay strict.

Development gate20 groups passed after fixing checkpoint selection. Final source-
bound gate adds missing-lineage and mixed-root cases. Local prescribed thermo
oracle retains CT2e-12 budget (observed3.55e-15). Reader author independently ran
42 native/synthetic corruption/compatibility cases. Final runs/plots follow.

## Final local gate

All21 native groups pass on final source-bound executables. Final47 event-reader
checks and28 unchanged CRV snapshot-reader checks pass (96groups total); native
detector8checks and EVT01–08/11–14 checks also pass. Native protocol modes are
included in SM20; do not double-count them as additional runtime groups. Exact
run executable hashes were rechecked (39native process records). Built-in and
prescribed Serial/double builds and MPI compile pass.

Measured planar behavior: 107timesteps,23saved clocks; sampled mode has23particle
cache calls and23curve stencil batches, no histories. No requested shock fields
means0detector calls. Complete mode has108particle cache calls and131curve batches
(108history +23snapshot queries). Full histories retain15particle/16curve events
and6particle censored views exactly. Coarse CRV only saves0/107 with zero curve
inside flags, while explicit complete tracking retains all16 between-snapshot
passages. SnapshotMC mode is independently checked against native mesh samples.

Cycle80 checkpoint preserves15unflushed completed curve records. Cycle106 particle
checkpoint preserves4records before bounded overflow at107; increasing operational
capacity restores all15 without identity loss/duplication. Cycle37 active and
zero-step replay, ordinary clean-prefix restart, fork-tail rejection and repaired-
checksum pending semantics all pass. Negative IDs and the full-width curve layout
are exercised; source thresholds/initial physics/tolerances are unchanged.

Nine actual figures explain call frequencies, between-snapshot passages, original
completion versus flush clocks, pending checkpoint recovery, measured encounter
files, sampled curves, AMR, accepted short forcing and independent thermo error
(3.55e-15 within the existing2e-12 budget). Encounter file counts change from110 to12
for particles and110 to9for curves; bytes112296→31836 and135312→26218 respectively.
This combines format/cadence changes and is not a wall-time or memory benchmark.

Retained development failures include terminal-prefix ordering, ambiguous same-
cycle checkpoint selection and the query fixture's terminal-ancestor dependency.
The final47reader suite separately tests uncovered queries and corrupted terminal
ancestors. Source, runtime, raw hashes and compact reference inputs/data are being
sealed. No push. Hardware/production/scientific gates remain pending.

Final timing validation also rejects invalid-fluid endpoints from a generation after emission, while preserving unknown-entry sentinel(-1,0). All47 event-reader cases and final21 native groups pass on the exact final source.
