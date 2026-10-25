# Passive material loops and open strings

This is the bounded passive, Newtonian Cartesian material-curve implementation
ported to `multiphase_turb` at `526a187f`. The original `codex/material-curves`
implementation used a different source/Kokkos baseline; its results are historical
controls rather than qualification of this port. Target-specific builds and tests
are recorded separately. CUDA/HIP sampling and push kernels remain a candidate
path requiring explicit opt-in and external qualification. This is not a
scientifically qualified turbulence release. See
[the evidence record](material-curves-evidence/README.md) for exact revisions,
commands, measurements, failures and exclusions.

## Supported envelope

| Configuration | Current status |
| --- | --- |
| Newtonian Cartesian single-fluid HD, ideal/isothermal, RK2, double precision | Implemented; tested CPU configurations |
| 2D/3D, every active boundary periodic | Implemented; 2D controls and bounded 3D controls tested |
| Open strings, zero/nonzero-winding closed loops | Implemented; line/circle/polyline initialization |
| Spacing and bending/chord-error insertion, fixed persistent anchors | Implemented; analytic strained-circle placement and cap/history tests |
| CPU Serial, one/two/four MPI ranks, empty curve-home ranks | Tested locally |
| Dynamic AMR, leaf ownership/cache refresh, complete-step checkpoints | Tested bounded prescribed-flow controls; broader physical/interface studies remain open |
| CUDA/HIP | Candidate device sampler/pusher; default rejection until explicit unqualified opt-in; no local compile/device evidence |
| Host topology, lookup, adaptation, audit, MPI buffers and ASCII I/O on GPU builds | Intentional implementation; transfers and scaling need measurement |
| Other integrators, 1D, single precision, SYCL, physical exits, reflecting/shearing boundaries, MHD/GR/radiation/ion-neutral/ordinary particles/forcing/cooling | Explicitly rejected |
| User source hooks | Prescribed-field fixture qualified; arbitrary user source/refinement hooks require their own stage/restart audit |

The target branch's Monte Carlo particles and all feature-off fluid paths retain
their existing code. Monte Carlo cell-centre jumps do not replace continuous
velocity-following curve samples; coexistence is currently rejected. Active
standalone cooling, frame tracking, turbulence forcing, and shearing/orbital
advection are rejected using their actual module state, independently of whether
cooling-history output is enabled. Disabled cooling/frame blocks are permitted.
Curves do not write fluid arrays. Curve-enabled AMR moves regridding before curve
completion/output so a shared checkpoint contains the completed mesh/curve state.
This changes AMR/output scheduling for curve-enabled runs; feature-off scheduling
is preserved. There is no level-subcycled time interpolation policy.

## Build and small examples

Use separate build directories. Initialize the repository's pinned Kokkos
submodule if it is absent; do not substitute another Kokkos version.

```sh
cmake -S . -B /absolute/path/build-hd -DCMAKE_BUILD_TYPE=Release
cmake --build /absolute/path/build-hd -j4
/absolute/path/build-hd/src/athena -i inputs/material_curves/line.athinput \
  -d /absolute/path/line-results
/absolute/path/build-hd/src/athena -i inputs/material_curves/loop.athinput \
  -d /absolute/path/loop-results
```

These evolve actual uniform HD using the built-in shock-tube generator with
identical left/right states. The loop inserts samples on the initial analytic
circle; its original anchors remain fixed in material identity.

For a laptop-scale turbulent showcase, use the same default HD build:

```sh
/absolute/path/build-hd/src/athena \
  -i inputs/material_curves/turbulence-2d.athinput \
  -d /absolute/path/turbulence-results
```

This 128 × 128 periodic box uses the branch's native one-time Fourier velocity
perturbations, then evolves unforced compressible HD to time 1.5. It seeds three
circles and three open strings, with adaptive insertion and fluid/curve snapshots
at matching cycle intervals. The native initializer supplies three velocity
components on a 2D spatial grid; curves follow the in-plane velocity. This is a
**decaying** flow, with numerical dissipation and no physical viscosity or ongoing
forcing. It demonstrates stretching and sampling; it does not establish a
converged turbulent spectrum. No solver or source machinery changes are needed.

The AMR example uses a **prescribed** translation fixture rather than a general
Euler solution. It is useful for inspecting scheduling/ownership/restart. On
`multiphase_turb`, `<amr_criterion0>/method=user` explicitly selects the fixture's
timed refinement hook; merely enrolling `user_ref_func` does not activate it.

```sh
cmake -S . -B /absolute/path/build-fixture -DPROBLEM=material_curve_test \
  -DAthena_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build /absolute/path/build-fixture -j4
mpiexec -n 2 /absolute/path/build-fixture/src/athena \
  -i inputs/material_curves/amr-restart.athinput -d /absolute/path/amr-results
mpiexec -n 4 /absolute/path/build-fixture/src/athena \
  -r /absolute/path/amr-results/rst/material.00001.rst \
  -d /absolute/path/amr-resumed time/nlim=-1
```

Inspect actual checkpoint names and recorded times first. Restart only a
completed checkpoint; retain its `.rst.curves` sidecar. The regression harness
selects and verifies legal checkpoints automatically. On the local Mac,
`/opt/homebrew/bin/prterun` was required instead of the incompatible installed
`mpiexec` wrapper; this is an environment limitation, not a dependency to install.

## Input and defaults

Absence of `<material_curves>` disables the feature. An existing block defaults
`enabled=true`; use `enabled=false` to suppress it. Units follow the simulation.

| Parameter in `<material_curves>` | Default | Meaning |
| --- | --- | --- |
| `enabled` | true when block exists | Allocate the separate passive-curve module |
| `ncurves` | 1 | Number of `<curve1>`... blocks; zero is legal |
| `curve_id_bits` | 8 | Integer 0–30; high curve bits in a 32-bit node ID; zero selects one logical curve |
| `max_nodes` | 1000000 | Global node cap; must be in `[1,INT_MAX/3]` |
| `adaptive` | true | Insertion-only geometric refinement |
| `spacing_cells` | 1 | Maximum spacing relative to the globally finest current active-cell width |
| `max_spacing` | 0 | Positive absolute cap, combined with cell-based cap; zero uses only cell width |
| `max_angle` | 0.2 | Maximum neighboring-edge bend in radians, in `(0,pi)` |
| `chord_error` | 0 | Positive length tolerance for the midpoint cubic correction; zero uses `0.02*effective_spacing` |
| `output_dcycle` | 0 | Positive cycle cadence; initial/final snapshots also written |
| `allow_unqualified_accelerator` | false | Explicit candidate CUDA/HIP execution permission; does not qualify a device/backend |

The cell-based policy uses the smallest active width anywhere on the mesh. It is
conservative: a fine patch inside a long edge influences refinement even when
both endpoint samples remain coarse. It can oversample curves far from that patch.
Fixed sampling ignores spacing/bending thresholds; quality reason zero then means
structural/degeneracy checks succeeded, rather than certifying geometric accuracy.

Each `<curveN>` accepts a unique integer `id` within its configured curve-bit
namespace, `type`, `closed`, and
integer `winding1/2/3` (default zero). Open strings cannot have winding. In 2D, x3
is constant and transverse winding is zero. Seed coordinates are continuous lifts.

* `type=line`: `count` (default16), origin `x1/2/3`, endpoint `end1/2/3`.
  An open line includes both endpoints. A closed line excludes the duplicate
  endpoint and requires `end=origin+winding*domain_length`.
* `type=circle`: `count` (default32), center `x1/2/3`, positive `radius`, and
  nonzero `normal1/2/3` (default `(0,0,1)`). Must be closed with zero winding.
  Any 3D plane is supported. Its label period is the exact smooth circumference.
* `type=polyline`: `coordinates=lifted`, and `points=x y z; x y z; ...`.
  Labels use the initial polygon arc length. Explicit winding defines closure;
  positions alone cannot disambiguate a curve drawn across a periodic seam.

Minimum count is2 for open strings and3 for closed loops. Seeds must be finite,
nondegenerate and within the reliable lift/image range of `2^40` domain images.
True polyline kinks retain unresolved quality; linear subdivision of a kink does
not demonstrate recovery of a previously unspecified smooth curve.

## Identity, time integration and refinement

A node ID and each neighbour/parent link is `int32_t`. Set `material_curves/curve_id_bits` to $B$ (default 8,
allowed 0–30): the high $B$ bits hold the curve ID and the remaining $32-B$ bits
hold the creation sequence. For $B>0$, curve IDs range from 1 through $2^B-1$.
For $B=0$, only logical curve ID 1 is allowed and the full bit pattern holds the sequence.
Sequences range from 1 through $2^{32-B}-1$; zero remains a missing-link sentinel.
Packing/decoding uses unsigned bit arithmetic; stored IDs use the corresponding
signed representation. For packed bits $u$, the displayed ID is $u$ below $2^{31}$,
and $u-2^{32}$ otherwise. Negative IDs, including $-1$, are valid; only zero is
missing. Parent creation order is checked using sequence bits, not signed numerical
order. The full namespace is retained; no sign bit is discarded.
The separate uint64 allocator counter safely represents the exhausted sentinel,
including $2^{32}$ for $B=0$. No IDs wrap or get reused. At least two nodes are
required, so allocating 31 or 32 bits to the curve field is invalid.
The existing global `max_nodes`, MPI batch bounds and 100000-curve IO cap still
apply; ID capacity is not a resident-memory or throughput qualification.

```ini
<material_curves>
ncurves = 1
curve_id_bits = 0
```

Keep `curve_id_bits` fixed during a run and restart. Current files declare the
split; the normal readers and curve CLI decode it automatically.
Neither IDs nor winding pass through `Real`. Immutable `a`
is an initial arc-length material label. Reciprocal `prev`/`next` IDs and sorted
material labels define connectivity, independently of storage/rank ownership.
No nearest-neighbor spatial reconnection, deletion, smoothing or
adaptation-induced anchor relocation is performed. Coordinates retain a continuous lift through every RK stage; wrapping
is only for fluid lookup. Closure adds the explicit winding displacement.

A curve's home rank is `curve_id % nranks`. Nodes and connectivity reside on that
home rank, while every stage routes batched interpolation requests to the current
active fluid leaf owner. This avoids global curve replication. A single long
curve's geometry is still concentrated on one rank: it is not distributed node
storage. `migrations` counts changes in serving fluid rank, not node-memory moves.

The RK2 predictor samples initialized `u^n` at `X^n` before fluid stage1. The
corrector samples `u^(1)` after stage1 boundary fill/prolongation/primitive
conversion at the lifted predictor position. The push combines lifted coordinates
as Heun RK2. Every stage resolves current sampling ownership. Cell-centered
primitive velocities use bilinear/trilinear interpolation with valid ghost
stencils and finite positive-density checks. Shock-width/interpolation error is
part of the physical error budget, not removed by the RK2 label.

Refinement runs before advection and at completed-step mesh/curve synchronization.
Spacing, neighboring-edge bend and a limited nonuniform Hermite midpoint correction
trigger insertion. Initial analytic circles are sampled from the known seed;
otherwise a guarded cubic placement falls back to a linear midpoint when its
correction is unsafe. New nodes retain both parent IDs, actual birth time,
placement method and a new label. They have no fabricated pre-birth trajectory,
thermodynamic or shock history. Original anchors remain persistent.

The global cap rejects an entire pending split batch, preserving coherent topology.
The first unresolved time/reason persists even if later refinement resolves the
current spacing. Severe per-step loss records the beginning of the affected step;
inserted points cannot repair unobserved earlier folds. Discrete midpoint/bending
indicators are not a rigorous continuous-curve error guarantee.

| `reason` | Meaning |
| --- | --- |
| 0 | No recorded structural/geometric failure under the active policy |
| 1 | Global node cap prevents requested insertion |
| 2 | Label/ID/subdivision exhaustion |
| 3 | Degenerate/nonfinite edge |
| 4 | Initial nonsmooth kink or severe geometry loss |

## Output, restart and diagnostics

Snapshots are `BASENAME.curves.initial.txt`, `BASENAME.curves.CYCLE8.txt`, and
`BASENAME.curves.final.txt`. A restart writes `FLUID.rst.curves` alongside the
**shared** fluid checkpoint. Rank-sharded fluid files are rejected. Curve snapshots and sidecars are
atomically replaced through a temporary file; native fluid files retain their
legacy write behavior, and incomplete/mixed pairs fail validation. Checksummed
ASCII has explicit
I32 IDs/links, U64 counters, I64 winding/images and I32 flags, with declared Real width. Reals use `max_digits10`.
The sidecar fingerprints the complete fluid file plus time/cycle/domain/config;
copying an unrelated fluid restart with a matching timestamp does not make a pair.
Checksums detect accidental corruption, not malicious modification.

The current schema restores identities, winding, material labels, allocator,
ancestry/birth/placement, persistent quality and counters. Completed-step RK scratch
is reset to the restored lift; stale ownership caches are discarded. AMR scheduler
clocks are part of the synchronized pair, mapped by fluid logical leaf identity
and restored independently of rank count. Thresholds/adaptive mode and AMR
scheduler/prolongation settings must match; the node cap can increase and output
preferences can change. Arbitrary stateful user refinement hooks are unqualified.
C++ preflight streams one curve at a time before changing live state. Unknown
versions, corruption/truncation, mixed pairs, malformed integer widths, duplicate
IDs, invalid links, bad ancestry/birth/labels/geometry/config or AMR clocks fail
clearly. Checkpoint schema versions from development are not a stable compatibility
promise; current output is schema 4 with explicit curve-bit metadata. Geometry restarts
also accept fixed-8-bit schema 3 when the configured split is 8. Python analysis
retains historical schemas 2 and 3 without changing their IDs. Tracked restarts
retain their separate exact-source/contract guards; schema compatibility alone
does not permit continuation from different code. Split changes reject explicitly.

Length sums each lifted edge once, including winding closure. Edge stretch is
current chord length divided by its immutable initial material interval `da`.
For a circle this compares to initial smooth arc length, so its initial polygon
stretch is slightly below1. Curvature is the three-point unsigned circumcircle
estimator `2*|e_prev cross e_next|/(|e_prev|*|e_next|*|e_prev+e_next|)`.
Open endpoints have mask1; collapsed/reversing/nonfinite stencils have mask2 and
NaN curvature. A valid local curvature stencil does not override persistent
curve-level unresolved history. Physical shock curvature and adaptive curvature
accuracy remain unqualified.

```python
import importlib.util
spec = importlib.util.spec_from_file_location('curves', 'vis/python/material_curves.py')
reader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reader)
state = reader.read('/absolute/path/material_loop.curves.final.txt')
for curve in state['curves']:
    geometry = curve['reconstructed']
    print(curve['id'], geometry['length'], curve['reason'], curve['unresolved_since'])
    # Exclude unresolved curves for a qualified observable; keep their fraction
    # as a separate result. Never backfill newly inserted samples before birth.
    if curve['reason']:
        continue
    edge_lengths = geometry['edge_length']
    # Normalize edge integrals by current arc length when that is the declared
    # scientific weight. Uniform node averages are biased by adaptive sampling.
```

The convenience reader holds all requested nodes, with a default bounded cap.
For many saved times, stream one snapshot at a time. Endpoint masks, birth coverage,
current-arc versus initial-material weights, per-loop weights and unresolved
fractions must be selected explicitly for scientific statistics.

## Reproducible checks

```sh
python3 -B tst/scripts/material_curves/analytic.py --output /absolute/path/oracles.json
python3 -B tst/scripts/material_curves/harness.py > /absolute/path/harness.json
python3 -B tst/scripts/material_curves/production.py \
  --executable /absolute/path/build-fixture/src/athena \
  --output-dir /absolute/path/new-results --suite all --launcher mpiexec
python3 -B tst/scripts/material_curves/physical.py \
  --executable /absolute/path/build-hd/src/athena --artifacts /absolute/path/physical
python3 -B tst/scripts/material_curves/evolved_shear.py \
  --executable /absolute/path/build-fixture/src/athena --output-dir /absolute/path/shear
python3 -B tst/scripts/material_curves_io/test_reader.py
python3 -B tst/scripts/material_curves_io/test_cpp.py \
  --build /absolute/path/build-fixture \
  --kokkos-library /absolute/path/build-fixture/kokkos/core/src/libkokkoscore.a \
  --launcher mpiexec
```

Result directories must be fresh. Scripts retain exact inputs, commands,
executable hashes, measured errors, frozen thresholds, snapshots and failure logs.
Oracle checks alone do not establish production accuracy. Prescribed-flow and
actual unforced-fluid controls are labeled separately. Native fluid passivity
comparisons apply at that output format's float32 precision.

Small local timing pilots use `pilot.py`; they include MPI startup, setup and
initial/final I/O. They are not strong/weak scaling, kernel profiles, cluster/GPU
resource qualification or production workload recommendations.

## Deferred GPU qualification

There is no CUDA/HIP machine available in this session. The operator will perform
these checks later. Configure the supported AthenaK CUDA or HIP backend with its
pinned Kokkos and actual device architecture, double precision, MPI as appropriate,
and the `material_curve_test` generator. Run a bounds-checked small build first.

```sh
python3 -B tst/scripts/material_curves/production.py \
  --executable /absolute/path/gpu-build/src/athena --output-dir /absolute/path/gpu-small \
  --suite core --allow-unqualified-accelerator
```

Then run all suites with the same opt-in, one/two/four ranks as available. Run
`physical.py --allow-unqualified-accelerator` on a built-in HD GPU build and
`evolved_shear.py --allow-unqualified-accelerator` on the GPU fixture build. For manual input files,
set `material_curves/allow_unqualified_accelerator=true` explicitly. Never weaken
the CPU numerical budgets merely to pass a device test. Capture compiler/Kokkos
versions, GPU model/architecture, execution-space log, commands, inputs, rank
mapping and executable hash. Use the backend profiler/Kokkos profiling hooks to
prove both `material_curve_velocity` and `material_curve_push` kernels execute on
the device, including empty-rank batches, post-birth insertion, MPI crossings,
regridding and changed-rank restart. A successful process or backend name alone
is insufficient. Check CPU/device common-label parity and convergence, run memory
checking, and measure staging/adaptation/I/O cost before accepting this backend.

Host topology/adaptation/lookup and host-buffer MPI are explicit ceilings. GPU-aware
MPI and within-curve device/distributed topology are not implemented. The default
GPU rejection remains deliberate until external qualification establishes a
reviewed supported envelope.


## Shock-feature integration (2026-10-17)

`feature_dump` includes native MC shock encounters alongside material curves.
The existing no-particles/no-forcing curve guard remains in force. Curve nodes
have no shock flags or encounter histories. Requesting shock tracking without
particles now fails explicitly rather than silently producing no history.
Fresh curve core/regression/unsupported/persistence and native restart checks
passed locally; see [the integration record](../docs/feature_dump.md) for exact
source, commands, evidence and qualification limits. This merge does not qualify
shock histories on material curves or new GPU/AMR/forcing combinations.
