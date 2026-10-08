# Typed MC shock-state migration adaptation

Base: `origin/multiphase_turb` at
`526a187f938e38634e89a88192ebe9731c4d36bb`; development branch
`codex/shock-multiphase`. Delegated production ownership was limited to
`src/bvals/bvals.hpp` and `src/bvals/bvals_part.cpp`. No commits or full Athena
binary builds were performed by this agent.

## Changes

The existing MC arrays and disabled boundary path are retained. Enabled shocks
use separate uniform-grid GID and receive/compaction methods. The GID method
checks extent, parent block and finite active-axis positions. Supported physical
faces are periodic/outflow; unsupported faces and nonuniform/missing neighbors
fail explicitly. It checks adjacent-block bounds before indexing the neighbor
array. MC moves at most one cell along one axis in its current pusher; the
boundary check additionally prevents any beyond-adjacent-block access.

Outflow particles enter a destruction list before any missing-neighbor lookup.
Sorted original indices are passed to `EmitShockExits` before any arrays move or
shrink. The coordinator owns actual partial-event/censor semantics and callback
idempotence. Periodic coordinates use half-open upper seams (`>= max`), update
active-axis signed int64 image counters with saturation rejection, and leave
historical observation images unchanged.

The complete trivially-copyable EncounterState is packed alongside each MC
particle's real and integer fields, on separate tag2/MPI_BYTE messages using the
existing particle communicator. IDs, sequences, counters and images never pass
through Real. Native byte communication requires identical executable/layout/
ABI on all ranks. Aggregate receive/send counts, Real/int buffer products and
state byte counts are checked before multiplication/allocation or MPI calls.

Enabled receives test all three request kinds before modifying particle arrays.
The union of sent and destroyed original indices is validated for duplicates and
bounds, then survivors and arrivals are copied into complete new real/int/state
arrays and published after a fence. Send/receive clear methods wait and clear
state requests alongside the two legacy kinds. No state buffer or destruction
storage is allocated by the feature-off path.

## Collective and lifecycle audit

`CountSendsAndRecvs` retains the branch's existing particle-communicator
Allgather followed by Allgatherv. Enabled migration introduces point-to-point
messages, not another collective. Both completed receive paths call
`Mesh::UpdateParticleCounts` exactly once; the enabled method returns directly
rather than falling through to the legacy call. Pending receive retries call it
zero times. The branch's count method performs one world-communicator Allgather
and updates rank/global totals. Thus ranks with no local transfers still
participate in the same completed-path count collective as transferring ranks.

The branch's task graph is push -> newgid -> count -> irecv -> send -> receive ->
clearrecv -> clearsend -> AMR adjustment -> seeding. Send buffers are independent
of the arrays replaced during compaction. Shared AMR qualification/rejection,
state seeding/AppendParticles, native restart and observer timing are owned by
the coordinator and are not validated by this packet.

The legacy metadata derived MPI datatype is committed each exchange and is not
freed by the baseline implementation; this adaptation did not expand or repair
that inherited lifecycle. Metadata collective error handling and extremely
large rank metadata counts also retain baseline behavior. State payload
allocation/byte bounds are checked explicitly. These are distinct from bounded
event-output and restart tests, which belong to coordinator integration.

## Verification and reproduction

From this checkout, after the coordinator's shared headers are present and the
baseline CPU Kokkos build exists:

```sh
python3 evidence/shocks-mc/migration/run_boundary_fixture.py
python3 tst/scripts/style/cpplint.py --filter=-build/include_subdir \
  src/bvals/bvals_part.cpp src/bvals/bvals.hpp
git diff --check
```

The fixture compiles the actual production boundary source and globals against
actual landed shared headers. It provides minimal mesh/block/particle
constructors and count/exit callbacks, not a hydro solver. Serial and singleton
MPI builds use the baseline CPU Kokkos libraries; MPI's config.hpp overlay only
changes the MPI compile switch and adds mpi.h. Exact compiler commands, source
SHA256 values, outputs, expected failures and verdicts are in
`boundary_results.json`. The initial 13-case result is retained separately.

The initial suite passed periodic 2D/3D seams, active-axis images, byte-preserved
complete state with identities above 2^53, outflow with no neighbors, and
expected extent/displacement/image-limit rejection in Serial and singleton MPI.
The singleton MPI mixed case sends a real particle to itself while two other
particles are destroyed, exercises actual tag0/tag1/tag2 payloads, and checks
complete survivor/arrival state. Count stand-ins verify exactly one completed
count call. This is not cross-rank migration qualification or proof of collective
behavior across ranks.

First compile failed because the fixture lacked a GeodesicGrid destructor
stand-in; preserved in `compile-first-failure.log`. First periodic assertions
failed because CPU `create_mirror_view` aliased the live state expectation;
preserved in `fixture-first-failure.log`. Replacing it with independent
`create_mirror` storage fixed the fixture; production was unchanged. The project
cpplint filter, 90-column check and diff whitespace check passed for both owned
production files. Style output is retained in `style.log`.

No inter-rank run, GPU/backend run, physical MC qualification, real partial-event
callback, seeded-state persistence, AMR gate or native restart round trip was
executed by this fixture. Rusty/Frontier and scientific approval remain pending.

## Frozen production hashes

- `src/bvals/bvals.hpp`: `4022e5b850d3486e45632690b39aafca8c1aa5af798e10bfc47efcd14537e623`
- `src/bvals/bvals_part.cpp`: `65f00233b419953f15936f62c01469c090dd7e2a807aa2d197611155f4fe0a60`
- shared `src/shocks/encounter.hpp`: `84cc51a2af2b397b9f3360d34a39bcfbf0fe61aee4c059cb745b1c076cba10f4`

Next gate: coordinator exact-source full build and real MC/observer/outflow/
restart integration, followed by independently reviewed multi-rank qualification.

## Final bounded fixture result

Serial and MPI fixture compilation both passed. All 18 extended cases passed:
periodic2d/periodic3d/outflow, image maximum, adjacent-displacement and extent
rejection, feature-off zero allocation and duplicate removal rejection in both
builds; singleton MPI additionally passed mixed send/destruction loopback and
expected MPI byte-limit rejection before buffer allocation. The byte-limit
fixture deliberately supplies an oversized metadata count without allocating
that payload. Full exact commands and rejection outputs are retained in the
final boundary_results.json. No production change followed this run.
