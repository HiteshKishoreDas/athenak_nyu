# Native particle, curve and shock output conventions

The feature now follows the pre-feature `multiphase_turb` architecture: native
scheduled output classes, typed binary diagnostics and Python analysis readers.
There is no new HDF5 dependency and no change to the fluid equations, tracer
pushers, detector, encounter transitions or checkpoint state layouts.

## Enable the common native format

```ini
<job>
basename = run
tracer_output_format = binary

<output2>
file_type = prtcl_thermo_history
source = both
id = thermo
variables = rho,p,t,s,v1,v2,v3
dcycle = 10

<output3>
file_type = curves
id = shape
dcycle = 10
```

Use unused output block numbers in an existing input. Existing particle seeding,
curve definitions and shock configuration are still required for those modules.
`source=curves` supports curve-only runs. The accepted scope is to align with
`multiphase_turb`; `tracer_output_format=binary` is an engineering opt-in,
not a scientific parameter or a request to change physics.

`legacy` remains the default, preserving existing files/readers. One job setting
selects the new representation for curve geometry, curve thermo and both shock
carriers. Particle thermo remains **byte-compatible `.thp` v2** from
`multiphase_turb`; its int32 particle tags cannot safely represent uint64 curve
node IDs, so curve IDs are not squeezed into that layout.

Geometry now also accepts an ordinary `<outputN>` block with `file_type=curves`.
It uses `id`, native file numbering and either `dcycle` or `dt`. As with existing
outputs, a present `dcycle` takes precedence over `dt`. Geometry and thermo reuse
the same existing tracer-output clock advancement, including terminal/restart
handling. Set `material_curves/output_dcycle=0` if only scheduled geometry is
wanted; leaving it positive additionally retains the existing cycle snapshots.
Initial/final geometry snapshots continue to be written.

## Files and clocks

| Product | Binary mode | Clock |
| --- | --- | --- |
| Particle thermo | `prtcl_thermo_history/<base>.<id>.thp`, unchanged v2 | Native thermo output block |
| Curve thermo | `curve_thermo_history/<base>.<id>.cycle<C>.rank<R>.tbin` | Same selected thermo output block |
| Scheduled geometry | `curves/<base>.<id>.<number padded to 8 digits>.rank<R>.tbin` | Native geometry output block |
| Existing geometry snapshots | `<base>.curves.<cycle>.rank<R>.tbin`; initial/final names also include their actual cycle | Existing geometry controls |
| Particle shock events | `<shock_tracking/directory>/cycle<C>.rank<R>.tbin` | Every completed step, including empty generations |
| Curve shock events | Same journal directory under `curves/` | Every completed step |
| Optional shock traces | `trace<C>.rank<R>.tbin`, plus curve `birthtrace<C>.rank<R>.tbin` | Existing trace/birth observation controls |
| Censored terminal views | `suspended<C>.rank<R>.tbin` | End of run; separate from committed events |

Shock journals deliberately keep their every-step durable drain. A coarse
geometry/thermo output cadence cannot discard encounters or change detection
cadence. No dense sample grid, interpolated path or exact crossing time is
invented to join different streams. Initial/final geometry includes the actual
cycle in its native filename, so a resumed run does not overwrite an earlier
terminal snapshot. All native chunks use the existing durable immutable commit
helper: byte-identical replay succeeds; conflicts fail.

## Common native table protocol

`.tbin` is a distinct diagnostic format, not a mesh `.bin` file. Like AthenaK's
mesh output it has a versioned named header followed by numeric binary records.
All new diagnostic kinds share **one** table codec and reader:

- 32-byte NUL-padded magic `Athena tracer binary output`.
- Explicitly serialized uint32 version (1), endian marker (`0x12345678`),
  Real byte width, and UTF-8 metadata byte count.
- uint64 record count and completed cycle, Real sample time, int32 rank/count.
- Named typed columns, kind, compiled-source fingerprint, code-unit declaration,
  position convention and applicable lineage/contract metadata.
- Row-major values with declared `u64`, `i64`, `u32`, `i32` or `real` widths.
  No C++ struct padding and no integer IDs converted to floating point.
- uint64 FNV-1a footer covering all preceding bytes. This detects accidental
  corruption; it is not an authenticity/security signature.

Geometry stores exact IDs/links/parents/material labels/births, lifted positions,
images, curve status/counters/winding and the same existing geometry diagnostics.
Curve metadata is repeated per node, keeping one simple named-table protocol.
The header includes domain, boundaries, dimensions, adaptation thresholds and
leaf count. Full mesh topology and AMR refinement clocks remain in the native
fluid output/checkpoints and the existing checkpoint sidecar; a binary geometry
analysis file is **not** a restart file or a replacement mesh dump.

Thermo keeps the same canonical field names and containing-cell sampling as
`multiphase_turb`. Geometry/curve thermo positions are lifted; shock observations
and particle thermo positions are wrapped, with shock periodic images retained.
Shock records retain the complete five-observation event evidence, status,
validity flags and sequence, plus the existing curve birth/refinement extension.
Raw membership flags and encounter phase are separate trace quantities.

## Read everything through one entry point

Add `vis/python` to the Python import path:

```python
from tracer_history import read_run
from shock_history import read_trace as particle_trace
from curve_shock_history import read_trace as curve_trace

data = read_run('/absolute/path/to/run', thermo_stem='run.thermo')
particles = data['particles']       # (cycle, exact particle tag) -> observation
nodes = data['curve_nodes']         # (cycle, exact node ID) -> observation
node_pressure = nodes[cycle, node_id]['fields']['pressure']
events = data['curve_events']       # existing node/lineage/sequence identities
geometry = data['geometry']        # independently clocked snapshots, exact links
```

The entry point handles both legacy and native files with the existing semantic
readers. `available` identifies present products; absent optional products remain
empty. Orphan shock chunks, missing represented ranks, incompatible schemas,
truncation, checksum errors and duplicate event/generation identities fail.
Choose `thermo_stem` when multiple thermo streams exist. For a custom journal
location pass `shock_directory`; absolute paths are supported as in existing
restart tests. Traces are read individually to avoid loading every optional
per-step observation by default. `tracer_binary.read(path)` also exposes any
native table's `schema`, `metadata` and typed `rows` directly.

Particle tags and curve-node IDs remain separate namespaces. Join each population
by exact identity plus cycle/time and run lineage as applicable; never join
particle and node IDs merely because numbers match. The reader preserves actual
clocks and births, and never copies a parent node's history.

For curve-focused figures and record selection, use `vis/python/plot_curves.py`.
[Curve reader usage](curve-reader.md) covers legacy/native inputs, exact node IDs,
3D geometry, field profiles, node histories and shock evidence.

## Restart and compatibility

Fluid `.rst`, particle restart payloads, curve `.rst.curves` and paired shock
state layouts are unchanged. The selected binary journal representation is bound
into the existing shock contract. Restart validates all committed event bytes
through its generation using the existing digest; mixing representations during
a tracked restart is rejected. Source/config compatibility remains strict.
Historical checkpoints do not automatically become compatible with this changed
source. A binary final geometry snapshot is not accepted as restart input.

Native thermo/geometry files are analysis products; they are not an atomic joint
transaction with each other or with the event journal. Particle `.thp` preserves
its existing append semantics, including explicitly deduplicable identical
terminal rows. Existing crash/restart guarantees and limits are not broadened.

## Validation and reproduction

```sh
cmake -S . -B build-lagrangian -DPROBLEM=built_in_pgens -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian -j8
cmake -S . -B build-lagrangian-prescribed -DPROBLEM=material_curve_test -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian-prescribed -j8
python3 tst/tracer-output/check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/native-output-NEW
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/visuals.py \
  --run build-lagrangian/evidence/native-output-NEW \
  --output build-lagrangian/evidence/native-output-plots-NEW
```

Use new destinations. The ledger in `evidence/tracer-output/2026-10-23/` records
actual results, source/binary hashes, all commands, review findings and retained
failures. Six explained figures use actual native records. File-byte and single
local wall-time measurements are observations, not production performance
qualification. Small mostly empty binary journals can be larger because each
chunk repeats schema/provenance. No compression or new output framework was
added.

Supported evidence remains local Serial/double Newtonian ideal hydro, the
existing fixed-mesh and live-AMR fixtures, continuous-particle/curve coexistence,
curve-only runs and fixed-mesh MC thermo with shocks off. No new physics,
MPI/GPU/cluster qualification, production budget or scientific sign-off is claimed.
The existing Python reader byte/record caps are explicit; large analysis files
require raised bounds or incremental analysis, not silently skipped rows.

---

# Original inspection record

2026-10-23. Inspected `codex/shock-lagrangian` at
`66e436e130f39399a201539b3cc108d3cacc126c`; clean before this record.
User request: synergise output files and formats. Writer changes have not started.
The previous startup-summary change is committed locally, one commit ahead of
the remote-tracking branch; no push was performed during this inspection.

## Current files and reusable readers

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
the implemented native convention and compatibility default are documented above.

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
