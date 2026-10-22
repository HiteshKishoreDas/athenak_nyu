# Shared material-curve and particle thermodynamic histories

Material-curve nodes and ordinary fluid tracers now use one thermo output class,
one field parser/evaluator, and the same native output clock. No change to curve
advection, node insertion, particle transport, detector thresholds, or fluid equations.
This is an optional diagnostic. It works without shock tracking or standalone particles.

## Configure

Edit the existing thermo output block (or choose an unused block and unique id):

```ini
<output2>
file_type = prtcl_thermo_history
source = both
id = shared_thermo
variables = rho,p,t,s,eint,v1,v2,v3,vmag,cs,mach
dcycle = 10
```

`source=particles` is the backward-compatible default. `source=curves` requires
only material curves. `source=both` requires curves plus either continuous
`lagrangian_tracer` particles or native `lagrangian_mc` tracers. Shock histories
still require continuous Lagrangian transport; selecting thermo does not relax
that guard. Fixed-mesh MC+curves thermo and live-AMR continuous-particle+curves
thermo have separate checks. No new live-AMR MC qualification is claimed.

Use `dt` instead of `dcycle` for time-based cadence. If `dcycle` exists, native
output parsing takes precedence over `dt`; do not leave `dcycle=0` when switching
to `dt`. Initial and terminal outputs follow the native scheduler, including an
extra terminal sample when stopping between scheduled outputs. New curve nodes
have no inherited history. Their first thermo sample is the first selected output
at or after their actual birth; every-step output captures birth endpoints.

The exact shared field parser accepts the same aliases and removes duplicate
canonical names. Hydro fields: density, pressure, temperature, entropy,
internal_energy, v1/v2/v3, vmag, sound_speed, mach and available scalarN fields.
`default` selects the existing eight-field particle default. Missing `variables`
uses `particles/track_variables` for particles/both, and
`material_curves/track_variables` for curves. Missing fallback selects `default`.
Magnetic quantities remain unavailable for curves: curve MHD is unsupported.

Quantities are containing-cell fluid diagnostics, not independently carried or
conserved curve variables. In ideal-gas code units,
\(P=(\gamma-1)e_{\rm int}\), the temperature proxy is \(P/\rho\), entropy proxy is
\(\ln(P/\rho^\gamma)\), and \(c_s=\sqrt{\gamma P/\rho}\). Physical temperature or
entropy units need the simulation's unit and composition conventions. Diagnostic
velocity uses the containing cell; advection retains its existing interpolated
primitive velocity. Curves retain the current periodic Newtonian ideal-gas hydro,
RK2, double-precision support restrictions.

## Files, identities and analysis

Particle output stays `.thp` v2, with unchanged record layout. Curve output is:

`curve_thermo_history/<basename>.<id>.cycle<C>.rank<R>.csv`

Each curve chunk declares version, cycle, time, row count, rank and rank count.
Columns are `node_id,curve_id,birth,a,x1,x2,x3,gid,level`, followed by canonical
field names. Integers retain full decimal width; never load node IDs through a
float. `a` is the immutable material label, positions are continuous periodic
lifts, and gid/level describe the sampled current leaf. Particle positions are
wrapped. Join curve geometry and thermo using `(cycle,node_id)`; join shock
observations using their generation/node ID and the same run lineage. Particle
tags and curve node IDs are separate namespaces, even when numbers coincide.

Both populations have the same analysis layout (`time`, `x`, `gid`, `fields`):

```python
from tracer_thermo_history import read_particles, read_curves
particles, particle_times = read_particles('prtcl_thermo_history/run.shared_thermo.thp',
                                            allow_identical_replay=True)
curves, curve_times = read_curves('curve_thermo_history', stem='run.shared_thermo')
assert particle_times == curve_times
# curves[cycle, node_id]['fields']['pressure']
```

The shared stdlib reader lives in `vis/python/tracer_thermo_history.py`.
It preserves exact identities and checks schema, row counts, finite values,
immutable node metadata, duplicate observations, rank completeness within each
represented cycle and birth/time consistency. Detecting a completely missing
cycle additionally requires the input cadence/expected-clock audit; the reader
does not invent a missing sample. Select `stem` when multiple output streams share
a directory. Use a fresh run directory and retain input/source provenance.

Curve files use the existing durable immutable commit helper. Identical final
replay is idempotent; a conflicting chunk fails instead of overwriting data.
Particle append behavior remains compatible, including native identical terminal
rows (the reader can explicitly deduplicate those). There is no atomic joint
commit across the particle and curve output files. No additional per-node thermo
state is carried in checkpoints: persisted geometry, IDs/births, fluid state and
native output scheduling are sufficient. The thermo writer saves `last_cycle`
beside the existing `last_time` in the native output parameter block. A repeated
final write cannot advance the clock twice, and an off-cadence terminal write
cannot skip the next scheduled sample. Catch-up still follows native scheduling;
no intermediate fluid states or samples are fabricated. Older checkpoints without
`last_cycle` do not guarantee zero-step replay of an overdue output clock; that
edge is qualified only for newly written checkpoints. Resume checks compare exact field values
and identities. Keep field selection fixed during a stream; use a new output id
for a changed schema. Existing exact-source shock checkpoint checks still apply.

## Reproduce laptop checks

From `/Users/hitesh/hitesh/git/athenak_nyu/.worktrees/shock-lagrangian`:

```sh
cmake -S . -B build-lagrangian -DPROBLEM=built_in_pgens -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian -j8
cmake -S . -B build-lagrangian-prescribed -DPROBLEM=material_curve_test -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian-prescribed -j8
python3 tst/curve-thermo/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/curve-thermo-NEW
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/curve-thermo/visuals.py \
  --run build-lagrangian/evidence/curve-thermo-NEW \
  --output build-lagrangian/evidence/curve-thermo-plots-NEW
```

Destinations must be new. The runner retains exact commands, inputs, logs, failing
runs, binary hashes and compiled-source fingerprints. The renderer writes an HTML
viewer and four explained scientific figures from actual native records. Add
`--cadence-failure <retained-CT16-failure-directory>` for the fifth figure comparing
actual failed and corrected restart clocks. NumPy
and Matplotlib use the existing environment; no dependency was added.

The gates cover curve-only large uint64 IDs, aliases/shared clocks, unchanged
particle bytes and geometry, MC coexistence without shocks, live refinement and
coarsening, fresh node births, exact restart, simultaneous shock sampling, cycle
and time cadence, unavailable fields, corrupt chunks, an independent spatial cell
oracle, immutable replay, native forced flow, passive scalars, overdue-threshold replay and empty populations. The fixed
\(2\times10^{-12}\) comparison budget is inherited from LT/CSH. Neither tolerances
nor physical defaults were tuned to make these checks pass.

## Evidence and limits

`evidence/curve-thermo/2026-10-22/ledger.md` is the living evidence ledger. Exact
numbers, source binding and artifacts are recorded there after final verification.
Host primitive mirrors follow ordinary tracer output; selecting both currently
makes a separate curve mirror. Request batches obey existing max_nodes and MPI
int-count bounds. Files grow with selected cadence and node count. These local
checks are not memory/I/O/runtime-budget qualification for a production run.

RUSTY CPU/MPI, Frontier, GPU/threaded execution and multi-rank persistence remain
separate deployment gates. The short forced box is an engineering demonstration,
not converged turbulence or scientific approval. The inherited forcing thermal
oracle issue remains open. Human scientific sign-off stays pending.
