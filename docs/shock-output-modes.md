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
encounter block per carrier is supported; duplicates reject.
`file_type=prtcl_thermo_history` selects regular particles and `file_type=curve`
selects curve nodes. The former shared `source=` selector has been removed.
Remove `source=particles` from particle blocks. Replace `source=curves` by a
`file_type=curve` block; replace `source=both` by two separate blocks with the
same desired fields/cadence. Each block controls its own output clock.
All active tracer blocks with any `source=` value reject with migration guidance;
it is never silently ignored. Historical files remain readable.

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

## Meaning of each user-facing option

These descriptions were checked against the current native parsers and detector
at code commit d0359aab. Values below are implementation defaults, not scientific
approval of a detector configuration. Fields and thresholds use simulation code units.
Particle/curve initialization stays in its existing blocks.

### Output blocks

| Option | Meaning |
| --- | --- |
| `file_type=curve` | Binary curve geometry, connectivity, diagnostics and selected samples in `.crv`; `mode=encounters` selects event-only CRV instead. |
| `file_type=prtcl_thermo_history` | Ordinary particle thermo snapshots in `.thp`; `mode=encounters` selects event-only THP instead. |
| `id` | Filename identifier, defaulting to the file type. Give each stream a distinct explicit id. |
| `mode=snapshot` | Default: sample only at saved clocks. Requested shock fields cause detection at those clocks. Does not retain passages between clocks. |
| `mode=encounters` | Enable complete tracking every completed step (and curve-node birth). Write event records only at this block's cadence. Requires explicit shock or Mach fields in this same block; their selection does not trim the fixed encounter evidence schema. |
| `variables` | Comma/space-separated tracer fields, case-insensitive, aliases supported and duplicates removed. It is not named `fields`. |
| `dcycle` | Native cycle cadence. Positive: every N cycles plus native initial/final behavior. Zero: disable block. Negative: initial/final-only writing, while an encounter block still tracks every completed step. |
| `dt` | Alternative output interval in code time, at the first completed step reaching each threshold. Positive required to activate this alternative. If dcycle exists it takes precedence, so remove dcycle to use dt. |

Normal particle snapshots without variables use particles/track_variables, default
`default`. Curve snapshots without variables use `default,shock` when the detector
is enabled, otherwise `default`. Explicit field lists make both intentions clear.
`source=` is removed; active old settings reject with migration guidance.

### Selected fields

| Name/group | Meaning |
| --- | --- |
| `default` | Density, pressure, temperature proxy, entropy proxy, internal energy density and three velocity components. |
| `rho` / `density` | Fluid mass density. |
| `p` / `pressure` | Ideal-gas pressure, \(P=(\gamma-1)e_{\rm int}\). |
| `t` / `temperature` | Temperature proxy \(P/\rho\), not a direct Kelvin conversion. |
| `s` / `entropy` | Entropy proxy \(\ln(P/\rho^\gamma)\), not physical entropy units. |
| `eint` / `internal_energy` | Internal energy per unit volume, not per unit mass. |
| `v1,v2,v3` | Coordinate velocity components. |
| `vmag` | Velocity magnitude. |
| `cs` / `sound_speed` | Sound speed \(\sqrt{\gamma P/\rho}\). |
| `mach` / `mach_number` | Local flow speed divided by sound speed; different from a shock-jump Mach estimate. Alone in a snapshot, this does not invoke the shock detector. |
| `scalarN` | Existing passive scalar indexed by N; rejects indices outside configured nscalars. |
| `shock` | Expands into all six shock fields below. |
| `shock_valid` | Detector support/state valid for this sampled point. Zero means unknown/invalid, not proof that no shock exists. |
| `shock_region` | Membership in the unreduced numerical shock region; interpret with shock_valid. |
| `shock_mach_pressure` | Shock Mach estimated from the selected pressure jump. |
| `shock_mach_temperature` | Shock Mach estimated from the selected temperature jump. |
| `shock_mach_pressure_valid`, `shock_mach_temperature_valid` | Separate validity flags for the corresponding estimates. A valid detector sample can still have an unavailable strength estimate. |

These are sampled containing-cell diagnostics. Curve advection retains its existing
interpolated velocity; these selections do not modify transport. Always-present
curve geometry is not controlled by this thermo field list.

### Shared shock_tracking configuration

| Option (default) | Meaning |
| --- | --- |
| `enabled=false` | Make the detector available when true. Does not by itself enable complete histories or force evaluation every timestep. Selected shock fields/encounter outputs require it. |
| `lineage` (required for complete tracking) | User-chosen run identity retained across restart. Use a different identity for an independent run. Safe letters/digits/underscore/hyphen, maximum 128 characters. Sample-only mode can omit it. |
| `event_capacity=4096` | Bounds per-rank, per-population event buffers, including completed records awaiting flush and terminal views. Not the lifetime encounter count. Overflow fails explicitly; replay the prior paired checkpoint with a larger capacity. |
| `trace=false` | Optional qualification/debug observation CSV and metadata during complete tracking. Does not enable tracking itself. Keep false for normal binary output. |
| `directory=shock_events` | Root for encounter batches and optional trace data; curves use its curves subdirectory. Ordinary CRV/particle THP snapshot directories are unchanged. |
| `width=1` | Upstream/downstream comparison offset in cells on either side along the dominant temperature-gradient axis. Not a smoothing radius or physical shock width. Must be between 1 and available ghost width. |
| `min_mach=1.1` | Minimum pressure-OR-temperature jump Mach needed for shock-region membership after the other shock criteria pass. Does not threshold local flow Mach. |
| `mach_max=1000` | Upper admissible Mach for the temperature-jump inversion. A larger implied temperature-jump Mach invalidates that estimate. This is not a universal cap: the pressure-jump estimate is not capped by this parameter. |
| `grad_floor=1e-30` | Candidate requires temperature-proxy gradient magnitude above this threshold. Rejects flat/very weak-gradient candidates. |
| `div_limit=0` | Candidate requires velocity divergence below this threshold. Zero selects compression; a negative value demands stronger compression. Positive values are unsupported by this feature. |
| `density_floor=0` | Minimum sampled upstream density for shock membership. Zero disables this additional cutoff. Does not change hydro density. |
| `temperature_floor=0` | Minimum sampled upstream \(P/\rho\) for shock membership. Zero disables this additional cutoff. Does not change hydro temperature. |

The floors filter detection; they do not repair invalid/nonpositive fluid states
or modify the fluid equations. Gamma comes from the fluid EOS, not a separate
shock_tracking gamma input. Preserve qualified thresholds when comparing runs.

### Optional mesh shock BIN

`file_type=bin` uses AthenaK's existing mesh binary writer. `variable=shock`
(singular) requests the mesh shock field group. `data_precision=real` writes the
build's Real precision (64-bit for the currently tested build); the native default
`float32` writes 32-bit values. The same id and cadence conventions apply. This is
independent of particle/curve outputs and requests detection at its saved clocks.

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
