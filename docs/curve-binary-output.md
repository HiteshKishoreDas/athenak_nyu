# Unified curve output: `.crv`

Current shock cadence/history controls: [sampled versus explicit complete encounters](shock-output-modes.md).

User-directed consolidation, 2026-10-25. One binary snapshot contains curve
geometry, explicit connectivity, provenance, geometry diagnostics, mesh location
and selected fluid/shock samples at one completed-step clock. No curve THP or
separately scheduled geometry file is needed. Ordinary particle THP v2 remains
unchanged. Complete encounters are opt-in event-only binary CRV/THP streams; they are
deliberately not embedded in sampled snapshots. Checkpoint sidecars remain authoritative
for restart and retain their existing representation.

## Input

Keep the existing curve definitions, hydro and optional shock-tracking settings.
Use an unused output number and a unique output id:

```ini
<output3>
file_type = curve
id = curves
variables = default,vmag,cs,mach,shock
dcycle = 10
```

`file_type=curve` is the binary output name; `crv` remains a compatibility alias.
The extension stays `.crv`. `variables` uses the existing tracer field parser and aliases. `default` expands
to density, pressure, temperature, entropy, internal energy and three velocity
components. Derived speed, sound speed, Mach and passive scalars can be selected
with the same names as particle thermo. Without `variables`, CRV defaults to
`default,shock` when curve tracking is enabled, otherwise `default`. Asking for
shock fields without enabled qualified tracking rejects explicitly. Geometry,
connectivity, diagnostics and mesh location are always present.

Use `dt` instead of `dcycle` for time thresholds. Native outputs occur at the first
completed timestep reaching each threshold; no interpolated clock is invented.
Fresh initial and normal final snapshots follow the native output scheduler.
Same-cycle final output replays the identical immutable file.

`file_type=prtcl_thermo_history` now selects particles only. The former `source=`
selector and curve routing have been removed; active old blocks reject with
migration guidance. Use a separate `file_type=curve` block for curve samples.
For both populations, use two blocks with the desired fields/cadence. Historical
THP files remain readable. With any unified curve output active, the
old automatic initial/final and `material_curves/output_dcycle` ASCII snapshots
are suppressed. Explicit `file_type=curves` still requests historical ASCII
snapshots. Runs without unified curve output retain their previous ASCII behavior.

Files are `curves/<basename>.<id>.cycle<C>.rank<R>of<N>.crv`. Whole curves retain
their existing home rank. Each rank emits its own bounded file, including explicit
empty files, with no new all-rank geometry gather. Files commit durably through
the existing immutable writer; conflicting content at an existing identity fails
instead of overwriting it. A changed field selection requires a new output id.

## Contents

| Scope | Saved data |
| --- | --- |
| Snapshot | Cycle, time, precision, byte order, rank/rank count, curve-ID split, domain, BCs, mesh dimensions, adaptation settings, basename, field names, counts |
| Curve | Signed ID, open/closed flag, winding, initial material period, unresolved status/reason, allocator, insertion/migration/sampling counters, node count, current length |
| Node | Signed ID, prev/next links, insertion parents, anchor/placement flags, material label, birth, continuous position, periodic images |
| Diagnostics | Edge length, stretch, curvature and validity mask, owning leaf GID and AMR level |
| Samples | Requested thermo/scalar fields and optional six-field native shock group |

The shock group contains `shock_valid`, unreduced `shock_region`, pressure- and
 temperature-jump Mach estimates, and the corresponding two estimate-valid flags.
Unavailable estimates are NaN. Use validity and region together when selecting
shock samples. `mach` remains fluid speed divided by sound speed.

IDs, links, flags and levels never pass through floating point. Full signed32
identities, including negative values and -1, remain exact. Zero alone means a
missing link. Winding/images are signed64; allocator/counters are uint64, including
the exhausted full-width allocator sentinel. New nodes retain their own birth
and parent provenance; no parental history is copied.

## Reader and plots

```python
from material_curves import read
from tracer_history import read_run
from tracer_thermo_history import read_curves

frame = read('/absolute/run/curves/name.curves.cycle10.rank0of1.crv')
node = frame['curves'][0]['nodes'][0]
links = node['prev'], node['next']
pressure = node['fields']['pressure']

data = read_run('/absolute/run')
rows, clocks = read_curves('/absolute/run/curves')
```

`read_run` supplies both geometry and curve-node thermo tables from CRV without
resampling. It retains old particle THP, historical curve THP/CSV and ASCII
geometry reading. Explicit ASCII/CRV snapshots at the same clock must agree on
all shared geometry before the richer CRV snapshot is preferred. Multiple CRV
streams require `thermo_stem='basename.id'`. The existing `plot_curves.py` CLI
works on a CRV-only curve run, including absent event journals.

Readers validate checksums, declared widths/byte order, bounded counts, complete
rank sets, clocks, immutable node metadata, topology/provenance and diagnostic
masks. Fixed geometry settings/domain must agree across a stream; node sets and
counters cannot silently regress. Rank counts may change between complete saved
clocks for analysis; runtime MPI/repartition qualification is a separate gate.

## Packed binary v1 contract

All scalars are appended individually in the writer's declared byte order, with
no C++ struct padding. Reals have declared IEEE width4 or8; currently qualified
native runtime is Serial/double. Strings have uint32 byte lengths and UTF-8
bytes. The final uint64 FNV-1a checksum covers every preceding byte, including
names and field samples. This detects accidental corruption, not forgery.

In order:

1. Magic16 bytes (`ATHKCRV` plus9 zero bytes); uint32 version1, endian marker
   `0x01020304`, Real byte width, sampled-field count.
2. int32 cycle/rank/rank count, uint32 curve bits; Real time/domain[6]; int32
   dimensions[3]/BCs[6]; uint64 max nodes/curve count/node count; Real max spacing,
   spacing cells, max angle, chord error; int32 adaptive/output_dcycle; basename
   and newline-separated sampled-field names as two length-prefixed strings.
3. Per curve: int32 id/closed/reason; int64 winding[3]; Real period/unresolved_since;
   uint64 next_sequence/inserted/migrations/samples/count; Real current length.
4. Per node in that curve's material order: int32 id/prev/next/parent_left/
   parent_right/anchor/placement; Real a/birth/x[3]; int64 image[3]; Real edge
   length/stretch/curvature; int32 curvature_mask/gid/level; requested Real fields
   in the declared name order.
5. uint64 checksum; no trailing data.

Whole-curve restart sidecars, AMR refinement scheduling state and complete
encounter state stay in the paired native checkpoint machinery. CRV is analysis
output, not a replacement restart file.

## Qualification

Exact commands, source hashes, retained failed harness cases, reports and plots:
`evidence/curve-binary/2026-10-25/`. Native parity checks compare shared values
exactly against previous output, plus independent prescribed thermo and geometry
oracles with their pre-existing budgets. Parameter text embedded in native mesh
BIN changes when output configuration changes; mesh/fluid array comparisons stay
exact, without changing numerical tolerance. Local MPI compilation does not
qualify runtime MPI, GPU, RUSTY/Frontier performance or scientific approval.
