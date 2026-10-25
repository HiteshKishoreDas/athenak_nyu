# Curve shock IO contract

Implementation task: optional histories belong to immutable curve-node IDs; newborn
nodes receive fresh state, never copied parent histories. No particle schema changes.

## Reused machinery and layout

- `shocks/store.hpp` copies existing MC durable exclusive temp-write/fsync/link/
  unlink/directory-fsync helpers into one inline namespace for both carriers.
  `PutEvent(out,event,newline=true)` retains MC bytes by default; curves append
  their metadata with `newline=false`. No new dependency or generic output framework.
- MC Python parsing extracts `parse_event(row,lineage,generation)` with the same
  checks, allowing the new reader to reuse the exact124-column core.
- Curve store lives at configured `shock_tracking/directory + '/curves'`.
  Manifest: `AKCURVE-SHOCK-EVENTS 1`, lineage, contract, each on its own line.
- Each rank commits `cycleGEN.rankR.csv` every generation, including empty ones.
  Header `# AKCURVE-SHOCK-EVENTS 1 LINEAGE GEN COUNT`, contract comment and
  `# END COUNT`. Rows133: MC124fields, birth, first_level, first_spacing[3],
  last_level, last_spacing[3]. Levels are native logical refinement levels.
- `suspendedGEN.rankR.csv` uses `AKCURVE-SHOCK-SUSPENDED` and contains terminal
  right-censored views without consuming identities or changing encounter state.
- Optional raw33-column `traceGEN.rankR.csv` and `birthtraceGEN.rankR.csv`:
  MC28trace fields, birth, current_level, current_spacing[3]. Birth samples have
  their own files; they do not silently overwrite earlier same-generation traces.
- Reader functions: `read_history(directory,through_cycle=None)`,
  `read_suspended(directory,generation=None)`, `read_trace(path)` keyed by node ID.
  Events retain `particle_id` as the shared core field and add `node_id`; curve ID
  uses the declared layout of the 32-bit packed ID. Current contracts start with
  `curve-node-history-v3 B`, where B is the configured curve-bit count (zero means
  one logical curve ID 1). Historical v2 means fixed 8/24; v1 means 32/32.
  Traces use their manifests or explicit `curve_id_bits` argument; 32-bit trace
  IDs are never guessed from their numeric value. Python curve identities include the carrier:
  `('curve_node',lineage,node_id,sequence)`, distinct from the unchanged MC triple
  even when both carriers use identical numeric IDs. Censored views stay separate
  from events.

## Restart and failure boundary

Parent hooks must call `PreflightShockRestart` before opening/truncating the native
fluid checkpoint, `SaveShockRestart` after native close and geometry `Save`, and
`LoadShockRestart` after native read and geometry `Load`. These functions no-op when
histories are disabled. The parent must also reject disabling histories recorded in
older inputs; missing sidecars cannot start a fresh lineage.

Sidecar filename: `FLUID.rst.curves.shock.rankR`. Native-ABI version1 header magic
`AKCSHCK1`, version uint32 at byte8, endian marker, scalar/state/int widths, rank and
rank count, generation/time, local node count, bounded contract size, fluid and
geometry lengths/digests, committed history prefix digest. Then lineage+contract,
sorted complete `ShockNode` records, trailing FNV-1a checksum. Same ABI and rank
partition only. Different ranks, missing/corrupt/torn state or prefix, wrong source/
configuration, node IDs or birth times reject before history state is installed.
Checksums detect accidental corruption; they are not hostile-collision authentication.

Sidecar loading validates into a temporary map, checking IDs/births against loaded
geometry, observation chronology/validity, phase/flags/rearm counters and AMR
metadata. An MPI barrier precedes state installation. Fluid/geometry/sidecars are
not an all-rank atomic transaction; incomplete pairs are explicitly unrecoverable.
The output cap bounds completed and suspended records; node state and optional
trace memory are bounded by the existing configured node cap.

## Checks

`python3 tst/shocks-curves/io_reader.py` exercises exact uint64 IDs, the shared MC
core,133/33-column layouts, AMR level/spacing, separate censored and newborn traces,
duplicate/missing generations, duplicate rank aliases, malformed ranges, nonfinite
metadata and invalid terminal censoring. Native checkpoint/replay and corruption
checks belong to the coordinated full acceptance runner; no device/HPC claim.
