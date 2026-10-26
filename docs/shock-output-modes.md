# Sampled shock properties and opt-in complete encounters

User-directed behavior, 2026-10-26. Normal CRV/particle THP snapshots evaluate
requested shock properties only at saved clocks. An explicit separate
`mode=encounters` block enables every-completed-step tracking for its carrier.
That block writes event records only, in binary, at its native output cadence.
It does not duplicate snapshots. Completed-event CSV journals are replaced.

## Input

Keep existing qualified hydro, particles/curves and detector thresholds.
`shock_tracking/enabled=true` makes detection available; it no longer enables
complete histories on its own. Use distinct output ids and unused numbers. `file_type=curve` writes binary `.crv`;
`crv` remains a compatibility alias. The older `curves` text output is unchanged.

```ini
<shock_tracking>
enabled = true
# Existing width/min_mach/floors remain your qualified settings.
# Required for complete histories; omit for sampled-only output if desired.
lineage = my_unique_run
# Bounds accumulated completed records as well as per-step/terminal records.
event_capacity = 4096
trace = false

# Ordinary curve snapshots: geometry + thermo + sampled shocks.
<output3>
file_type = curve
id = curve_samples
variables = default,vmag,cs,mach,shock
dcycle = 10

# Optional complete curve histories. No geometry/thermo snapshots in this stream.
<output4>
file_type = curve
mode = encounters
variables = shock
id = curve_encounters
dcycle = 50

# Ordinary particle thermo snapshots remain native THP v2.
<output5>
file_type = prtcl_thermo_history
id = particle_samples
variables = default,mach,shock
dcycle = 10

# Optional complete particle histories. No duplicated thermo snapshots.
<output6>
file_type = prtcl_thermo_history
mode = encounters
variables = shock
id = particle_encounters
dcycle = 50
```

Delete/disable blocks4/6 for sampled-only work. Add either independently: a curve
history block does not enable particle histories or vice versa. Exactly one active
encounter block per carrier is supported; duplicates reject. Encounter THP blocks
use the existing `source=particles` default; omit it. The source selector was added during curve/particle thermo sharing; it is not
a pre-feature multiphase_turb option. The THP snapshot writer also accepts `source=curves|both` for existing workflows. `file_type=curve` already selects
curves and needs no `source` setting.

Every active `mode=encounters` block must explicitly set `variables` containing
`shock` or `mach`. The existing parser also accepts individual `shock_*` fields,
`mach_number`, case-insensitive names and comma/space separators. A missing list,
`default` alone or density/pressure alone fails before tracking starts. Shock fields
in another snapshot block do not satisfy this requirement. Disabled encounter
blocks require no fields and enable no tracking. These names are opt-in selectors:
event-only packets still contain the complete fixed encounter evidence schema.
Ordinary sampled `mach` remains local flow speed divided by sound speed; it does
not imply shock detection unless complete encounters are explicitly enabled.

Use `dt` instead of `dcycle` for native time thresholds. Presence of `dcycle` takes
precedence. Zero disables that block; negative `dcycle` retains native initial/final
only writing, while the explicit encounter block still enables complete tracking.
The default mode is `snapshot`. Output cadence controls saving/flush, never exact
interpolation of a crossing time.

| Configuration | Detection and state | Output |
| --- | --- | --- |
| No requested shock fields; no encounter block | No shock evaluation | Normal chosen fields |
| Shock fields in snapshots only | Evaluate at those saved clocks; no encounter state | Flags/strength/validity in sampled CRV/THP |
| Explicit curve encounter block | Every completed step and actual new-node birth | Complete curve events in binary CRV batches |
| Explicit particle encounter block | Every completed step, independent of snapshots | Complete particle events in binary THP batches |

Particle mesh caches are refreshed lazily and shared by same-clock consumers.
Curve snapshots evaluate native node stencils without advancing encounters.
An explicitly requested mesh shock BIN output also requires detection at its own
output clocks. Snapshot-only mode cannot claim complete histories: entire passages
can occur between saved samples. Native MC fluid tracers may sample instantaneous
properties; complete histories still require true lagrangian_tracer particles.
All existing detector/backend/physics qualification limits remain.

## Files and readers

Sample file paths and binary layout stay unchanged from the unified CRV work.
Complete histories use `shock_tracking/directory` (default `shock_events`):

- Particles: `<directory>/<basename>.<id>.cycle<G>.rank<R>of<N>.thp`.
- Curves: `<directory>/curves/<basename>.<id>.cycle<G>.rank<R>of<N>.crv`.
- Terminal censored views use `.suspended<G>` instead of `.cycle<G>`.

These are typed event packets, not THP-v2 sample records or geometry CRV records.
Distinct magics prevent accidental interpretation by old snapshot readers.
Existing `shock_history.read_history/read_suspended`, curve equivalents and
`tracer_history.read_run` discover binary history automatically; historical CSV
remains readable. Public `shock_binary.read` decodes a single batch or checkpoint
pending packet. Whole-run readers check complete rank sets, duplicate identities,
strict prefix links and original emission generations. A through_cycle query
requires a committed batch covering that clock, even if events are empty.

Known magics distinguish snapshots/events when custom directory paths share a
folder/extension. Unknown or corrupt candidates still fail. Keep output ids
unique; immutable files never silently overwrite a colliding stream. One history
stream per carrier/directory is supported. Optional `trace=true` qualification
observations remain CSV, with their metadata; these are not completed-event
journals and are not written in sampled-only mode.

## Buffering, restart and censoring

Completed events retain their original emission generation, timing brackets,
strength evidence, validity/censoring flags and exact identities while waiting
for the next flush. The bounded queue fails explicitly on overflow. Replay a
prior paired checkpoint with a larger event_capacity; increasing this operational
capacity does not alter the detector contract. Choose cadence/capacity for the
number of events that can accumulate, not only one timestep's events.

Paired checkpoints include both per-carrier encounter state and unflushed completed
records, including records whose particles no longer exist. Pending bytes are
inside the checksummed fluid-bound sidecar. New particle sidecars are version3;
new curve-shock sidecars version2. Source/ABI/configuration guards still apply;
old tracked checkpoints are not silently converted or reinitialized. Complete
history activation/id/lineage cannot be dropped or changed during continuation.
The output's original history stem persists independently of a new checkpoint job
basename. Directory relocation requires the complete referenced binary prefix.

Restart rejects committed tail/branch files outside the checkpoint's exact saved
prefix. For an older checkpoint, copy only that saved prefix (and matching views)
to a new history directory; preserve the original run. Do not delete files in
place to force continuation. Pending-packet headers declare each carrier's last
flush; different carrier cadences may have different boundaries. Terminal censored
views are separate, non-consuming records. Final flushing precedes terminal-view
creation and paired checkpoint writing, making zero-step replay byte-identical.

Binary packet layout is specified in
`evidence/shock-output-modes/2026-10-26/interface.md`. All identities/sequence and
clock counters remain integers; packet checksums detect accidental corruption,
not malicious forgery. Output buffers are bounded; production memory/I/O budgets
still need deployment qualification.

## Validation and local reproduction

See `evidence/shock-output-modes/2026-10-26/reproduce.md` for source-bound commands,
reports, retained failures and actual plots. The canonical new native suite is
`tst/tracer-output/shock_output_modes_check.py`, with strict reader checks in
`shock_binary_reader_check.py`. Old per-step CSV-specific acceptance scripts are
historical format evidence; they must not be cited as passing the new binary mode.
Physical initial conditions and previous numerical budgets are retained. Local
Serial/double checks, MPI compilation, production performance and human scientific
approval are distinct evidence levels.
