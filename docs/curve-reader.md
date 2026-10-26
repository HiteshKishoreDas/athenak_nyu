# Read and plot material curves

Current curve output is unified `.crv`; see [input, contents and reader usage](curve-binary-output.md).
Historical THP/CSV and ASCII reader support remains. The details below record the earlier split output.

`vis/python/plot_curves.py` uses the supplied existing readers to draw:

1. Initial, middle and final stored curve geometry, with 3D axes for 3D runs.
2. The selected thermo field versus immutable material label at its latest saved clock.
3. One exact node's saved thermo history and actual birth time.
4. That node's emitted shock episodes and separate censored terminal views,
   retaining validity flags and available entry/exit brackets.

The plotting script launches no simulation and changes no runtime state. Existing
reader calls retain their defaults; this script selects `include_particles=False`
to avoid loading unrelated particle streams. Current THP curve samples, retained
CSV samples and text geometry are read directly, with no conversion step.

## Laptop example

From `/Users/hitesh/hitesh/git/athenak_nyu/.worktrees/shock-lagrangian`, using the
existing scientific environment:

```sh
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  vis/python/plot_curves.py \
  build-lagrangian/evidence/variable-id-qualified-output-20261024/planar \
  --curve-id 7 --field pressure --output /private/tmp/curve7-pressure.png
```

Choose a new output path when rerunning. Existing figures are never overwritten.
For your own run, replace the run-directory argument. `--field` accepts the
canonical saved field name, such as `density`, `pressure`, `temperature`, `entropy`,
`v1` or `mach`; it does not compute fields that were not saved. Missing fields
produce an error listing available names. Proxy quantities retain the existing
code-unit conventions described in `curve-thermo-history.md`.

`--node-id 117440545` selects a particular node and infers its curve ID. Enter
IDs as decimal integers, never floats. Current files declare their configured
curve-bit split, including zero-bit single-curve runs; the CLI does not guess
that split from the number. Without a node selection the smallest
recorded node ID is used, normally the first anchor. A single curve is selected
automatically; multiple curves require `--curve-id`. Conflicting selections fail.

Use `--thermo-stem basename.output_id` when several curve thermo streams exist.
Use `--shock-directory PATH` for a custom journal location, including absolute
paths. The script prints a compact JSON summary containing the exact IDs, sample
counts, birth/first-sample times and selected-node event/view counts. Redirect
stdout to retain that summary if desired.

## Use the selected records in Python

```python
from tracer_history import read_run
from plot_curves import select_curve

run = read_run('/absolute/path/to/run', include_particles=False)
curve = select_curve(run, curve_id=7, field='pressure')
node_history = curve['history']       # saved observations, ordered by cycle
latest_profile = curve['profile']     # current saved values, ordered by material label
snapshots = curve['frames']           # independent geometry clocks and connectivity
```

Importing these functions does not import Matplotlib; only rendering needs it.

## Interpretation

Each panel uses its own saved clock. Geometry is not inferred from thermo samples;
thermo is not interpolated onto geometry output times. Geometry-only and thermo-only
inputs produce explicit unavailable panels. Line segments guide the eye between
saved samples; markers are actual records. A node can first be saved after its
birth when the thermo output cadence is coarser than insertion timing.

Closed curves use stored winding and domain lengths for their closure endpoint,
not nearest-image wrapping. Positions remain continuous lifts. The current profile
orders nodes by their immutable material label; labels are not renumbered after
insertion. No parent history is copied into a newborn's history.

The event panel's thick segment spans first/last inside evidence. Dotted segments
show available outside/inside timing brackets; they are not exact physical crossing
times. Incomplete entry is labeled unknown. Partial emitted events and suspended
terminal views are distinct. Numeric flags preserve the existing shock schema.

The script loads curve data into memory using the existing reader caps; it is not
a streaming analysis tool or a production memory qualification. Optional raw
per-step flags remain accessible through `curve_shock_history.read_trace`. Keep
the rank manifest alongside a trace for automatic layout decoding. For an isolated
32-bit trace, explicitly pass `curve_id_bits=B` and `signed=True` for current signed IDs;
missing layout is rejected. `id_layout=(32,B,True)` declares signed storage; older
unsigned layouts retain `(32,B)`. Never discard negative IDs as missing.

## Reproduce checks

```sh
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/curve_reader_check.py \
  --runs build-lagrangian/evidence/variable-id-qualified-output-20261024 \
  --three-d-run build-lagrangian/evidence/variable-id-qualified-CT-20261024/forced_both \
  --output build-lagrangian/evidence/curve-reader-NEW
```

This checks actual saved planar/AMR/3D datasets, scheduled/control equality, exact
current signed32, previous unsigned32 and historical uint64 CLI IDs, optional missing diagnostics, particle IO opt-out with unchanged
default strictness, winding closure, invalid selections, input corruption and
figure preservation. No simulation is launched by the reader checks. See
`evidence/tracer-output/2026-10-24/curve-reader-ledger.md` for the original reader results; `existing-formats-ledger.md` records the format-removal checks and new plots.

## Current binary curve thermo

`read_run` and `read_curves` discover native THP v2 chunks automatically:
`curve_thermo_history/<basename>.<id>.cycle<C>.rank<R>of<N>.bits<B>.thp`.
The same four signed integer slots hold cycle/node ID/curve ID/GID; birth, a and
refinement level are named metadata fields. The filename declares layout and the
rank set even for empty chunks. Historical CSV is still read directly.

The returned curve observation interface and plot_curves CLI are unchanged.
When shock properties were selected, they are available under fields using the
same canonical names as regular particle THP and mesh BIN output. Only the two
named shock Mach estimates may be NaN; flags and clocks are validated. Plots use
actual saved clocks/positions and do not invent pre-birth samples or crossing times.
