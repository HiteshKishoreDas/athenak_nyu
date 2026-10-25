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
