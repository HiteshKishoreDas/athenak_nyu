# Combined material-curve and shock-history showcase — 2026-10-21

## Request and scope

User requested an actual showcase of both features working together. Show native
fluid evolution, material-curve deformation, unreduced shock membership, and
histories on the same nodes at synchronized times. Continuous standalone particles
also remain present. Use the already accepted physical setup rather than inventing
new initial conditions, forcing, boundaries, floors or resolution choices.

Final viewer and outputs:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/combined-curve-shock-showcase-20261021/showcase-reviewed/`.

- `index.html`: playback and per-frame slider, explanations beside figures.
- `showcase.gif` / `.mp4`: nine actual stored fluid frames, no time interpolation.
- `showcase.png`: a synchronized frame at cycle 76, with active encounters and
  previously completed passages visible together.
- `curve_evolution.png`: four states of the same curve; both inside and later
  outside samples, with persistent episode counts and measured polyline length.
- `node_passages.png`: all 28 retained node episodes, sequence numbers, sampled
  entry/exit brackets, representative strength evidence and per-node counts.
- `report.json`, `media-check.json`, `manifest.json`: exact bindings and checks.

Older `visuals/` and `showcase/` render directories are superseded drafts. The
reviewed viewer makes raw flags suppressed by the episode filter visible in purple.
The earlier sealed 2026-10-20 implementation/evidence bundle remains untouched.

## Native run and source

Native executable SHA256:
`a543768ceb0a7a16f3bf6ffff53adaf6b932af5065bdf382f2471d3f5e106421`.
Scientific implementation commit `1c325e5a59d64f743ffdaa4f8678a4805b12608f`;
compiled source fingerprint:
`1a8ee8ed78fb74a35234b3ba0417e6e8aba1d5bad37c109992404406bd6e9c56`.
This task adds analysis/documentation, with no C++/physics changes.

Replayed the accepted native forced \(32^3\), periodic, Newtonian ideal-hydro
fixture, with 4096 continuous RK2 particles and a fixed 16-node material curve,
through \(t=0.2\). Initial velocity RMS 3 is an initializer, not a claim of evolved
stationary Mach-3 turbulence. Exact physical input remains
`inputs/tests/shocks-curves/forced.athinput` and the parent acceptance-run input.

Only changes: output basename/directory/lineage and `output3/variable=hydro_w`
instead of `hydro_u`. Primitive output permits direct supplied-state detector
checks without reconstructing internal energy from total energy in Python.
Distinct lineage avoids mixing separate output stores. Native run completed with
exit 0 in 5.673 local wall seconds. Its full input/command/log and hashes reside in
`../native-run/execution.json` relative to the final viewer.

## Reused machinery and checks

Used the running-sim workflow and Ponytail: existing `bin_convert`, strict curve
and particle readers, and the compiled `shock_detector_test` adapter. No new
solver, physics generator, detector implementation or dependency was added.

1. Read nine native primitive snapshots at cycles 0,21,41,59,76,93,110,127,144.
   Join to curve geometry and histories by exact cycle; the fluid header time is
   rounded, so it is checked against the exact native geometry/trace clock.
2. Feed the exact saved density, velocity and internal-energy primitives to the
   existing C++ detector, retaining width 4 and the same complete configuration.
   Periodic adapter ghosts are appropriate to this uniform periodic fixture.
3. Compare validity and inside flags at every displayed carrier position against
   its containing detector cell. **37008 exact comparisons pass**: nine times
   4096 particles plus 16 curve nodes. Positions match native curve geometry at
   the unchanged \(2\times10^{-12}\) reconstruction budget.
4. Replaying with primitive output leaves all 145 particle trace files byte equal
   to the conserved-output parent: 593920 particle observations. All 2320 curve
   observations equal. Complete particle/node events and terminal views match
   except the deliberately distinct lineage/identity namespace.
5. Counts remain 6904 complete particle events,28 complete node events and517
   separate terminal censored particle views. No terminal curve censoring.
6. GIF and MP4 each contain exactly nine frames. Inspected poster, curve evolution
   and node passage plots visually. Source/script/media hashes are retained.

The curve's native 3D polyline length reaches \(L/L_0=7.4121\). This is a measured
coarse polyline diagnostic, not proof of geometric/trajectory convergence.

## Visible scientific limits

Raw unreduced region membership is distinct from accepted encounter episodes.
Of 775 raw inside node samples,548 are suppressed by the inherited simultaneous
upstream/downstream ordering check; one completed node event has a direction-
ambiguity flag. These facts are shown, not hidden. Purple node samples are inside
but suppressed; orange are inside an active episode; blue are outside. This
showcase does not scientifically validate that suppression rule in turbulence.

Gray extensions on the passage plot are actual stored timing brackets. Some are
broad because earlier raw flags were suppressed. Squares show single-inside-
sample episodes; no artificial nonzero duration was invented. Per-node episodes
do not count unique global shock surfaces or whole-curve shock intersections.

Fluid density and flagged-cell fractions are projected along \(z\); apparent
projected overlaps may occur at different depths. Node flags themselves are
checked in full 3D. Orange contours represent fractions 0.1,0.3,0.5 of line-of-
sight cells flagged, not 2D physical shock surfaces. Only 1/16 of standalone
particles are drawn in the spatial panel; all4096 enter counts and native checks.

Fixed uniform mesh and fixed16-node curve in this showcase: no AMR/adaptive-node
qualification claimed here. Prior engineering AMR evidence stays separate.
Short local Serial/double demonstration; no converged turbulence, cluster/GPU
qualification or human scientific approval. The inherited forcing thermal-oracle
issue remains unresolved. No physical choices or numerical tolerances changed.

## Reproduce the renderer

Use the existing NumPy/Matplotlib/Pillow environment and installed FFmpeg:

```sh
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/shocks-lagrangian/showcase.py \
  --run /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/combined-curve-shock-showcase-20261021/native-run \
  --detector build-lagrangian/shock_detector_test \
  --output build-lagrangian/evidence/showcase-UNIQUE
```

A fresh output destination is required. Exact detector input fixtures, commands,
signals and reconstructed masks are retained per frame. Native source/state is
never modified. The original run command is stored in `execution.json`; give
fresh output paths/lineage before repeating it to avoid immutable-store conflicts.

## Inspection issues retained in this record

Initial exploratory inspection treated the binary reader's meshblock list as an
array; corrected by using its existing NumPy conversion pattern. First media check
used system Python, which lacks Pillow; reran in the already installed visualization
environment and passed. Neither issue changed numerical data or tolerances.

## Curve thermodynamic capability inspection (2026-10-21)

User asked whether curves can track density/pressure like ordinary tracers.
Inspected `material_curves.hpp`, `curve_shocks.cpp`, `curve_shocks_io.cpp`,
`tracer_fields.cpp`, `prtcl_thermo_history.cpp` and the strict history readers.
Current optional `<shock_tracking>/trace=true` already writes node-local density,
pressure and temperature proxy P/rho each completed-step observation, including
outside shocks. Event-only records contain encounter snapshots, not continuous
thermodynamic histories. Generic selectable fields/cadence independent of shock
tracking are available for standalone fluid tracers, not wired to curves.
The minimal future extension reuses the existing field evaluator and leaf-owner
sampling; no advection/connectivity change is needed. Inserted nodes start at birth.
These are sampled local fluid values, not independently conserved curve variables.

Plotted existing native data, no solver change or simulation rerun:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/combined-curve-shock-showcase-20261021/thermo-demonstration/`.
Run `plot.py` there with `vis/.venv/bin/python` (set writable MPLCONFIGDIR).
Result: 145 native trace files, 16 nodes, 2320 observations; strict reader,
unique identities/generations, increasing sample times, validity, birth bounds,
and P/rho consistency (relative tolerance 2e-14) passed. Three predetermined
nodes shown (first/middle/last sorted IDs), containing-cell samples in code units.
Source SHA-256 hashes and selected IDs retained in report.json. This is evidence
of existing output contents, not a new independent physical validation. Existing
forcing/convergence/scientific-review limitations remain as above.
