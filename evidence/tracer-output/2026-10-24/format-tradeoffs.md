# Output-format tradeoffs: discussion only

2026-10-24. User asked whether text and CSV are best for curve/shock data.
Inspected current source at `b99cd8ff` on `codex/shock-lagrangian`; clean before
this note. No output format, solver or configuration changes authorized by this
question, and none made.

Text/CSV are a compatibility/readability choice, not a demonstrated production
optimum. Existing particle thermo is already binary THP v2. Its MPI implementation
gathers records to rank zero, so binary representation alone does not establish
parallel scalability. Curve geometry and thermo contain exact uint64 node IDs
and connectivity; the unchanged THP v2 particle tag layout cannot hold them.

The existing planar output fixture contains 1 particle THP file (1,455,897 bytes),
108 curve thermo CSV files (351,582 bytes), 108 particle event CSV chunks
(104,012 bytes), 108 curve event CSV chunks (133,109 bytes), and 109 cycle/initial/
final ASCII geometry snapshots (419,061 bytes), excluding scheduled geometry in
`curves/`, traces and other products. Populations and field contents differ;
these are inventory measurements, not a text/binary size or speed comparison.
The fixture is `build-lagrangian/evidence/existing-output-final-20261024/planar`.

File proliferation follows per-rank/per-cycle journal scheduling and immutable
commit/replay semantics. Changing numeric encoding alone does not remove that
architecture. Any later batching must preserve bounded event handling, exact
identities, restart prefix validation and durability; no such change is made here.

For large histories, standard HDF5 is a candidate: typed datasets, chunking,
compression and MPI I/O are documented by The HDF Group. It is not a demonstrated
winner for this code or target. Compression/chunk choices can harm performance.
Adopting HDF5 still changes the output contract/build dependencies and is not
implied by the user's preference for existing workflows.

Primary references consulted:
- https://portal.hdfgroup.org/documentation/hdf5/latest/_h5_d__u_g.html
- https://portal.hdfgroup.org/documentation/hdf5/latest/improve_compressed_perf.html

Recommendation: keep existing runtime files and common reader now. Profile actual
output time, file counts and memory on the intended RUSTY/Frontier workloads;
address output cadence/batching if measurements justify it. Consider a standard
container only with explicit agreement and measured evidence. No benchmark,
new simulation, speedup, MPI/GPU qualification or scientific approval is claimed.

## Follow-up: reuse THP and AthenaK mesh binary

The user proposed curve thermo in THP, mesh shock fields in AthenaK `.bin`, and
particle shock properties in THP. Source inspection supports this direction as
an extension of the existing output families, rather than another independent
codec. This is feasibility discussion; no implementation or tests changed here.

- Mesh shock diagnostics can be added as named output variables backed by the
  synchronized detector signals and written by the existing mesh binary writer.
  The current output-variable registry does not yet expose these shock fields.
- Existing particle THP has a named scalar-field list and can gain optional
  sampled shock diagnostics. The current tracer field parser has no shock-field
  selectors; implementing them must use synchronized encounter/detector state,
  not confuse fluid Mach with pressure/temperature-jump shock Mach.
- Curve thermo can use the THP family with separate curve files/population
  metadata and a versioned extension for exact uint64 node IDs, birth/material
  labels and refinement metadata. Current THP v2 uses four int32 fields
  (cycle, tag, seed_id, gid), so it cannot directly represent existing curve IDs.
  Do not put IDs in floating-point fields or silently reinterpret the v2 tags.
  Preserve old particle v2 writing/reading when no extension is requested.
- Sampled shock properties are not the full encounter history. Episodes can
  start/end between thermo output times. Preserve every-step durable event
  records, exact sequence/lineage identities, timing brackets, flags and censored
  outcomes even if THP is later extended to carry event records as well.
- Curve geometry/connectivity is a separate product; THP thermo adoption does
  not by itself encode ordered links, winding or geometry diagnostics.

Inspected `src/outputs/prtcl_thermo_history.cpp`, `src/outputs/outputs.hpp`,
`src/outputs/basetype_output.cpp`, `src/outputs/derived_variables.cpp`,
`src/particles/tracer_fields.{hpp,cpp}`, `src/shocks/{detector,encounter}.hpp`,
`src/material_curves/material_curves_io.cpp`, `src/driver/driver.cpp` and
`vis/python/tracer_thermo_history.py` at `a560920b`. No benchmark or new scientific
approval is implied. These are proposed mappings, not already available options.

## Why current curve-node IDs are 64-bit

Source verification: `MaterialCurves::NewID` in
`src/material_curves/material_curves.cpp:144` checks curve ID and monotonically
allocated node sequence against UINT32_MAX, then packs the two components:

\[
\mathrm{nodeID}=2^{32}\,\mathrm{curveID}+\mathrm{sequence}.
\]

This is an identity design choice, not a physical requirement or a statement
that a run contains billions of nodes. Even curve 1, sequence 1 is already
larger than an unsigned 32-bit integer. Nodes retain this identity as positions,
owning mesh leaves and curve resolution change; allocation sequence persists
across restarts. A single 32-bit global ID scheme is possible but would change
allocation, identity limits and persisted link/history contracts. Storing two
32-bit components separately preserves the information but still needs explicit
curve semantics in the writer/reader; it is not the existing particle-tag schema.
The previous statement that curves "require" 64-bit IDs should be read as:
current packed IDs must be preserved exactly, not that every curve design needs
64-bit IDs. No code changes or new tests made for this explanation.

## Where connectivity lives

Each `Node` stores exact neighbour IDs in `prev` and `next`; the owning `Curve`
also stores its node vector in immutable-material-label order. `parent_left` and
`parent_right` record insertion ancestry, not necessarily current neighbours.
`Links()` rebuilds neighbour links after insertions; `Audit()` checks reciprocal
links against curve order. Open ends use zero; closed curves wrap last to first.

Current ASCII geometry snapshots serialize these fields in each `node` record.
The paired fluid checkpoint sidecar (`<fluid checkpoint>.curves`, normally
`*.rst.curves`) serializes them too. Current curve thermo CSV has node identity,
birth, material label, positions, leaf metadata and fields, but not neighbour
links. Adopting THP for thermo alone would not replace geometry/connectivity
output. Verified in `material_curves.hpp:24`, `material_curves.cpp:150`, and
`material_curves_io.cpp:129,217,649`. No source changes or tests needed.

## Latest identity choice: configurable split

The fixed 8/24 implementation is now the default of the user-requested
`material_curves/curve_id_bits` option (0–30). Zero means one logical curve ID1
and the entire unsigned32 node ID holds the sequence. Layout metadata is explicit
and immutable on restart; current readers do not guess it from numeric IDs.
Evidence, compatibility rules and 107 passing local groups are consolidated in
`evidence/curve-ids/2026-10-24/variable-ledger.md`. Existing global node/MPI batch
limits still apply. This does not itself implement THP curve/event output.
