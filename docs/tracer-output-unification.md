# Existing particle, curve and shock outputs

The output-format experiment was withdrawn on 2026-10-24 at the user's request.
No new diagnostic format or format selector is supported. Keep existing input
files and readers. The shared analysis entry point is the integration layer.
Historical evidence below and under `evidence/tracer-output/2026-10-23/` describes
the withdrawn experiment, not the current interface.

## Existing files stay existing files

| Product | Existing format | Reader |
| --- | --- | --- |
| Particle thermo | `.thp` v2, unchanged binary records | `tracer_thermo_history.read_particles` |
| Curve thermo | `.thp` v2 cycle/rank chunks (historical CSV readable) | `tracer_thermo_history.read_curves` |
| Curve geometry/connectivity | ASCII v5 `.txt` snapshots | `material_curves.read` |
| Particle encounters and optional traces | `.csv` journals | `shock_history` |
| Curve encounters and optional traces/birth traces | `.csv` journals | `curve_shock_history` |

Curve geometry is ASCII v5. Set `material_curves/curve_id_bits` (default 8,
allowed 0–30) to choose the split of the signed 32-bit node identity. Zero means
one logical curve, with all 32 bits used for node sequence. Current files declare
the layout; the supplied reader and curve CLI use it automatically. Existing
particle THP bytes are unchanged. Negative curve node IDs retain all 32 packed bits; zero alone is missing.
Curve thermo now uses the same THP v2 layout; metadata is declared in native
cycle/rank/bit-split filenames. Complete event journals remain separate CSV. Historical geometry/thermo/event analysis retains
its original IDs, and fixed-8-bit geometry v3 remains readable for compatible
geometry restarts. Tracked restarts retain exact-source guards. Changing the
split on restart is rejected.

## Same configuration pattern

Existing particle output blocks continue to work without edits. To sample curve
thermo in the same output block, add `source=both` (or `source=curves` for a curve-only
run). Fields and cadence use the existing settings:

```ini
<output2>
file_type = prtcl_thermo_history
source = both
variables = rho,p,t,s,v1,v2,v3
dcycle = 10
```

The optional `file_type=curves` output block schedules existing ASCII geometry
through AthenaK's normal `id`, `dcycle`/`dt` and file numbering. It introduces no
new representation. Existing `material_curves/output_dcycle` and initial/final
snapshots continue to work without adding this block:

```ini
<output3>
file_type = curves
id = shape
dcycle = 10
```

Scheduled files are `curves/<basename>.<id>.<number padded to 8 digits>.txt`.
Set `material_curves/output_dcycle=0` if only scheduled geometry is wanted;
otherwise both existing cycle snapshots and scheduled snapshots are written.
Shock journals retain their every-completed-step drain independently of these
output cadences. Checkpoint layouts and encounter rules are unchanged.
Source fingerprints still guard tracked restarts; do not assume a historical
checkpoint from different code is compatible merely because layouts match.

## One supplied analysis entry point

```python
from tracer_history import read_run

data = read_run('/absolute/path/to/run')
particles = data['particles']
nodes = data['curve_nodes']
events = data['curve_events']
geometry = data['geometry']
```

This calls the existing parsers directly. No conversion or extra data file is
required. Namespaces and actual sample clocks remain separate; no resampling,
exact crossing time or inherited newborn history is invented. Missing optional
streams remain empty; corrupt input and orphan journal chunks raise errors.
Use `thermo_stem` for multiple thermo streams and `shock_directory` for a custom
journal directory. Use `include_particles=False` for curve-only analysis.

[Read and plot curves](curve-reader.md) gives the ready-to-run plotting command.
The script accepts a run directory, exact curve/node IDs and a saved field.

## Reproduce laptop checks

```sh
cmake --build build-lagrangian -j8
cmake --build build-lagrangian-prescribed -j8
python3 tst/tracer-output/check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/existing-output-NEW
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/curve_reader_check.py \
  --runs build-lagrangian/evidence/existing-output-NEW \
  --output build-lagrangian/evidence/existing-output-reader-NEW
```

The output checks retain numerical, AMR birth and restart gates; only checks
specific to the removed codec are withdrawn. Existing CT/CSH/LT suites provide
independent physical oracles. Reports distinguish local Serial/double execution
from MPI/GPU/cluster or scientific approval, which remain separate qualifications.

---

# Historical inspection record (superseded)

2026-10-23. Inspected `codex/shock-lagrangian` at
`66e436e130f39399a201539b3cc108d3cacc126c`; clean before this record.
User request: synergise output files and formats. Writer changes have not started.
The previous startup-summary change is committed locally, one commit ahead of
the remote-tracking branch; no push was performed during this inspection.

## Files and readers at the 2026-10-23 inspection

| Data | Native representation | Existing reader |
| --- | --- | --- |
| Particle thermo samples | Append-only `.thp` binary v2 | `tracer_thermo_history.read_particles` |
| Curve thermo samples | Immutable CSV v1, cycle/rank chunks | `tracer_thermo_history.read_curves` |
| Curve geometry/connectivity | ASCII v2 snapshots | `material_curves.read` |
| Particle encounters | CSV v1 events, optional traces, separate suspended views | `shock_history` |
| Curve encounters | CSV v1 extended events, optional traces/birth traces, suspended views | `curve_shock_history` |

Sources inspected: `src/outputs/prtcl_thermo_history.cpp`,
`src/material_curves/material_curves_io.cpp`,
`src/particles/shock_particles.cpp`,
`src/material_curves/curve_shocks_io.cpp`, and the readers above under
`vis/python/`. Thermo already shares field selection/evaluation and scheduling;
its two on-disk representations differ. Existing fluid-output conversion in
`vis/python/bin_convert.py` provides a precedent for native output plus an
analysis conversion step. Existing `trk` output rejects fluid tracers and is
not a substitute for these histories.

## Facts a common representation must preserve

- Separate particle and curve-node identity namespaces; exact integer IDs,
  including node IDs above floating-point exact-integer range.
- Event lineage, sequence, quality flags, five evidence observations and timing
  brackets. Suspended terminal views remain distinguishable from committed events.
- Actual sample cycle/time and independent stream cadence. No interpolation or
  invented samples to make geometry, thermo and shock clocks appear identical.
- Particle wrapped coordinates versus curve lifted coordinates, with explicit
  position conventions. Geometry/connectivity is sampled at its own cadence.
- Actual node birth and immutable material label. No copied parent histories or
  pre-birth observations for inserted nodes.
- Input/source provenance, original schema versions and explicit unavailable
  data. No silent omission of malformed or incomplete source chunks.
- Existing strict readers and replay checks should be reused before introducing
  another parser. Checkpoint formats are distinct from analysis products.

## Choices pending at initial inspection (superseded)

Two UserQuestions submitted: conversion of existing output versus common native
writers; HDF5 versus CSV versus a common reader alone. Recommendation is an
analysis conversion first, retaining native outputs and checkpoint persistence.
This was a recommendation at inspection time. The later user request to align
with pre-feature multiphase_turb supersedes these pending engineering choices;
the subsequently withdrawn experiment is retained in Git history. The current interface is documented above.

### Follow-up: what pre-feature AthenaK does

The user asked how vanilla handles outputs before selecting a format. Inspected
local Git objects without switching branches or fetching remote changes:

- Cached upstream/main and local main: `5de631377124b2696201425c2675f8d46b45d7a2`.
  This is the inspected upstream snapshot, not a claim about today's remote HEAD.
- Pre-shock development base multiphase_turb:
  `526a187f938e38634e89a88192ebe9731c4d36bb`.

Upstream uses independent `<outputN>` blocks and specialized output classes,
selecting `file_type`, `dt` or `dcycle`, basename and applicable variable/id.
There is no single container for all mesh, particle, history and restart data.
Mesh `bin` output has a versioned text header (1.1 in this snapshot), followed by
binary MeshBlock records; shared-file and per-rank options exist. `trk` appends
tracked particle positions and velocities in binary to a separate `.trk` file.
`pvtk` supplies particle visualization, `hst` formatted global histories and
`rst` a distinct checkpoint format. Python `bin_convert.py` reads mesh binary
into a dictionary and optionally writes `.athdf` plus XDMF. Thus native runtime
output plus optional analysis conversion is already an AthenaK pattern; native
HDF5 for all outputs is not the inspected baseline.

The user's multiphase_turb base already adds `prtcl_thermo_history`: binary `.thp`
v2 with named fields, typed integer/Real records and MPI gather to rank zero.
That feature is absent from the inspected upstream/main tree. It must not be
attributed to stock upstream AthenaK, nor replaced with `trk`, which assumes
drift-particle velocity slots and explicitly rejects the MC fluid tracers in
multiphase_turb.

Verification commands used `git ls-tree <ref>:src/outputs` and `git show` on
`outputs.cpp`, `binary.cpp`, `track_prtcl.cpp`, `prtcl_thermo_history.cpp` and
`vis/python/bin_convert.py`. Recommendation now stated in baseline terms: retain
native scheduled output and separate checkpoint files, harmonize the analysis
interface and schema/identity metadata, and make a combined HDF5 export optional.
The user's question alone does not select a new native binary schema.

Local availability check:

```sh
/Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  -c 'import h5py; print(h5py.__version__)'
# 3.16.0
```

This establishes an existing optional Python dependency only. It does not
qualify native HDF5, parallel HDF5, RUSTY, Frontier, throughput or memory use.
No simulation, new implementation test or scientific plot was produced in this
inspection. Next: implement the selected scope, round-trip actual native
records through it, test exact IDs/flags/clock and replay semantics, and plot
the same particle and curve signals from the common representation against
the native reference.

## Native binary shock properties

With existing shock tracking enabled, add `shock` to the normal thermo variables
for particles, curves, or both. Keep the existing particle/curve setup and lineage.

```ini
<output2>
file_type = prtcl_thermo_history
source = both
id = thermo
variables = rho,p,mach,shock
dcycle = 1

<output5>
file_type = bin
variable = shock
id = shock
data_precision = real
dcycle = 1
```

The thermo group expands to the same six named fields written by the mesh group:

| Field | Meaning |
| --- | --- |
| shock_valid | Native touched-state/stencil and physical-endpoint validity |
| shock_region | Unreduced accepted shock-region membership |
| shock_mach_pressure | Native pressure-jump Mach estimate |
| shock_mach_temperature | Native temperature-jump Mach estimate |
| shock_mach_pressure_valid | Valid native sample and finite pressure estimate |
| shock_mach_temperature_valid | Valid native sample and finite temperature estimate |

`mach` retains its original fluid speed/sound-speed definition. Raw jump estimates
can exist on candidates outside the accepted region; use shock_region when
selecting accepted shock samples. An unavailable estimate is NaN and its validity
flag is zero. Other thermo fields must remain finite. Partial field selections
still validate every available prerequisite. No episode ID or sequence is stored
in a floating-point field.

Particle THP v2 retains its file path, headers and record layout. Curve THP uses
that same layout in immutable cycle/rank chunks under curve_thermo_history; the
supplied read_curves/read_run and plot_curves interfaces load them directly.
Historical curve CSV remains readable. A changed variable list needs a new output
id/directory, as before. No format selector or conversion command is required.

Mesh diagnostics use the existing binary writer and reader, including its default
float32 storage. Set data_precision=real for native Real precision. The shock mesh
group currently requires file_type=bin, ghost_zones=false and enabled qualified
shock tracking. It uses the synchronized particle detector cache, or the same
native detector for the qualified periodic curve-only case. No ghost-zone shock
values are fabricated. Zero-step restarts rebuild this derived cache after native
primitive initialization, without advancing encounters or journals.

Thermo fields follow their configured output cadence. Complete per-passage events,
timing brackets, censoring/validity flags and restart state remain in the durable
encounter journals/checkpoint sidecars. Curve connectivity/ancestry remain in the
geometry snapshots and sidecars. These independent products are joined by exact
population/node identity and generation; sparse thermo samples do not replace
complete encounter records.

Supported runtime evidence remains local Serial/double ideal Newtonian Cartesian
hydro on qualified uniform/live-AMR fixtures. MPI compile is checked separately;
MPI runtime, GPU/cluster performance and human scientific approval remain pending.
