# Particle and material-curve startup summary

Existing inputs now print a summary automatically, once after startup seeding and
adaptation, or after checkpoint populations are restored. Only rank zero writes
the summary. Counts are global (integer reduction in MPI builds).

Example from the accepted adaptive laptop fixture:

```text
Tracer setup summary: fresh, time=0 (global counts)
Particles: type=lagrangian_tracer, requested_at_start=128, actual=128
  placement=eligible cell centers, draws with replacement; weight controls birth selection
  tracer_seed_initial: requested/event=128, region=box, weight=volume, seed=903, start=0, end=0, cadence=once
Curves: requested=1, actual=1, nodes_requested=8, nodes_actual=15
  curve1: type=line, closed=false, placement=uniform spacing
  adaptation=true, max_nodes=1000000, spacing_cells=1, max_spacing=auto(mesh), max_angle=0.2 rad, chord_error=auto(spacing), unresolved_curves=0
  geometry snapshots=initial/final, every 1 cycles
Histories: particle_shocks=false, curve_shocks=false, particle_shock_trace=false, curve_shock_trace=false
  output3: prtcl_thermo_history, source=both, every 1 cycles, id=prtcl_thermo_history
End tracer setup summary
```

`requested_at_start` counts births due at this initial endpoint, not future
scheduled injections. Each seed block separately reports configured count/event,
region, mass/volume weighting, seed, start/end/cadence and optional field filters.
A region with no eligible cells can therefore show 128 requested and 0 actual;
a future-only injection shows 0 requested now while retaining its count/event.
Particle placement remains eligible cell centers, sampled with replacement.

Curve `nodes_requested` counts original seed vertices before adaptation, and
`nodes_actual` includes nodes inserted during startup. Initial anchors already
carry this distinction, so reporting adds no persistent state. Placement is
uniform angle for circles, uniform spacing for lines, or explicit lifted vertices
for polylines. Automatic spacing/chord settings are labeled `auto(...)`; the
summary does not claim these labels are resolved numerical tolerances. Adaptation
settings remain visible when adaptation is off, explicitly labeled inactive.
`unresolved_curves` reports existing geometry-quality flags.

On restart, populations are labeled `restored` and seed requests are marked
`not reapplied`. History status comes from enabled runtime modules and constructed
output streams, so disabled output blocks are not shown as active. No tracer
summary is printed when neither fluid tracers nor curves exist.

The change only reads resolved settings/counts and prints them. It changes no
input syntax, defaults, seeding/advection, adaptation criteria, output cadence,
file format or restart state. Explicit pusher selection and `ncurves` bookkeeping
remain unchanged, as requested.

## Reproduce the narrow laptop checks

From the isolated `codex/shock-lagrangian` worktree:

```sh
cmake -S . -B build-lagrangian-prescribed -DPROBLEM=material_curve_test -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian-prescribed -j8
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
 /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-setup/check.py \
 --executable build-lagrangian-prescribed/src/athena \
 --baseline /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-thermo-synergy-20261022/final-reviewed/runtime/build-lagrangian-prescribed/athena \
 --workdir build-lagrangian-prescribed/evidence/startup-summary-NEW
```

Use a new workdir. The baseline is the preserved previously reviewed executable.
The runner saves inputs, commands, summaries, logs, source/binary fingerprints,
actual count plots and an HTML viewer. Nine checks cover fresh/adapted/empty and
future births, circles/polylines, restored populations, enabled/disabled histories,
separate populations and absent modules. Counts are compared with native records;
the four-step paired run compares particle-history bytes and initial/final curve
state exactly. A separate built-in hydro smoke check is retained in the ledger.
These are local Serial/double checks of reporting, not additional scientific,
multi-rank, GPU or performance qualification.
