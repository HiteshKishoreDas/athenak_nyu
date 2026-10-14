# Native shock histories for existing Monte Carlo tracers

Branch `codex/shock-multiphase`, based on `multiphase_turb` at
`526a187f938e38634e89a88192ebe9731c4d36bb`. Engineering implementation;
**human scientific sign-off and R1 release acceptance remain pending**.

Ponytail audit: core fluid solver, fluxes, mesh and task graph remain unchanged.
The general particle constructor adds one call to private ConfigureShockTracking;
feature setup lives in existing shock_particles.cpp. The guarded boundary path
and paired restart hooks are retained because histories must survive migration,
outflow deletion and replay. No new framework/dependency/tracer is introduced.
Structural cleanup commit `f0053f2e`; source fingerprint
`80da5390b07a322b83f3ce829f97975c4326b3216a764d77bf530a8f355265c4`.
Independent exact-extraction review and rerun results are in
`evidence/shocks-mc/ponytail/`. Earlier 78e5 reports/plots remain historical;
all 21 physical measurements match exactly after cleanup, including failures.

Three-review cleanup evidence: source `3de18d97`, fingerprint
`3feb3bcdf77d21296516544e67b5ad75f6ac6ca7bc17c4b48be886ddc0fb7cbc`.
Forty net production lines removed: duplicate state/header/schema code and
redundant validation/sort. Encounter state is816bytes on tested local ABI,
down from1000; emitted124-column event schema is unchanged. Old sidecars reject
explicitly by size/source contract. Three independent reviews and fresh checks
are in `evidence/shocks-mc/review3/`; earlier packets keep their source identity.

## What changed and why

| Change | Purpose |
| --- | --- |
| Native C++/Kokkos detector, frozen ShockIt contract | Evaluate the live fluid state without dumping full grids; preserve unreduced shock-region membership |
| Completed-step particle sampling | Associate actual MC parent cells with the same final fluid time, after migration and scheduled births |
| Typed encounter state and event records | Preserve entry/exit brackets, strength evidence, validity, sequence and censored outcomes across timesteps |
| Complete typed migration/compaction | Move every encounter field with its particle; preserve integer IDs and periodic images; emit exits before deletion |
| Bounded immutable event chunks and strict reader | Reject overflow, damaged output and conflicting replay explicitly |
| Paired encounter sidecar beside native MC restart | Restore particle rows, seeding progress, RNG configuration and encounter identity together; verify the checkpoint/history prefix |
| Input/restart integrity guards | Restrict unqualified physics, reject ambiguous seed IDs/unsafe counts and check writes/close before committing a sidecar |
| Independent statistical, RH, parity and persistence fixtures | Distinguish numerical port fidelity, MC transport statistics and experimental physical encounter recovery |

The existing mass-flux tracer, stateless RNG, thermo output and seeding/restart
machinery are reused. No new velocity-following tracer is introduced. Positions
are cell centers; each completed fluid step selects at most one outgoing face.
Uniform-flow qualification uses ensemble statistics, not individual analytic paths.
The legacy RNG's additive key has cross-(tag,cycle) aliases; independence at research
scales is an outstanding scientific question. Native fixed-layout replay is exact.

## Support and experimental choices

CPU Serial, double precision; RK2; Newtonian ideal-gas Cartesian hydro; uniform
2D/3D mesh, one pack per rank; periodic/outflow boundaries; full synchronized
primitive ghosts. Local 2D and thin3D cases are measured. True multi-rank runtime
qualification remains pending despite an MPI build and singleton checks.

AMR, MHD, relativistic/radiation physics, cooling/other sources, diffusion,
frame tracking, shearing/orbital operators, FOFC, RK3/RK4, threaded/GPU/single
precision and per-node restart are excluded/rejected. Rank/partition changes are
not supported for restart. Material loops/strings are excluded.

Accepted engineering requirements: exact source/config identity, unchanged fluid
and MC behavior, complete state transfer, typed identities, valid timing evidence,
bounded output and explicit persistence failures. Experimental choices: temperature
normal, nearest-axis endpoints, first valid inside strength representative,
upstream log-distance entry guard and two valid outside observations to rearm.
Width1 and minMach1.1 are constructor defaults; the supplied pilot uses width4.
Strength/recovery/timing thresholds await the owner. A completed numerical episode
is not automatically a scientifically certified physical passage.

## Build and run

From this worktree (CMake uses the pinned Kokkos submodule):

```sh
git submodule update --init kokkos
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DAthena_SHOCK_TESTS=ON
cmake --build build-debug -j8
build-debug/shock_encounter_test
python3 tst/shocks/parity.py build-debug/shock_detector_test \
  --reference evidence/shocks/reference --report /tmp/mc-parity.json
mkdir -p /tmp/athenak-mc-pilot
cd /tmp/athenak-mc-pilot
/ABS/WORKTREE/build-debug/src/athena \
  -i /ABS/WORKTREE/inputs/tests/shock_mc.athinput
```

Replace `/ABS/WORKTREE` with this checkout's absolute path. Parity and legacy
hydro regression require NumPy; the strict history reader and MC helper checks
use Python's standard library. The recorded environment uses ShockIt's existing
`.venv/bin/python`; package/version identities are in the evidence manifests.
Run directories must be new to avoid conflicting immutable output.

```sh
python3 tst/shocks-mc/integration.py --executable build-debug/src/athena \
  --workdir /tmp/mc-integration-new --report /tmp/mc-integration.json
python3 tst/shocks-mc/qualification.py --executable build-debug/src/athena \
  --workdir /tmp/mc-heldout-new --report /tmp/mc-heldout.json
```

The first suite checks engineering invariants; the second records physical
measurements without inventing individual MC crossing times or scientific tolerances.

## Settings and diagnostics

Use existing `particle_type=lagrangian_mc`, `pusher=lagrangian_mc` and
`<tracer_seed...>` blocks. Enable `<shock_tracking>/enabled=true` and supply a
safe unique `lineage` (letters/digits/underscore/hyphen, maximum128characters).
`directory` defaults to `shock_events`. Resolved configuration is printed and
stored in each rank manifest; source fingerprint and native EOS/solver/mesh/RNG/
static seeding definitions are bound to replay. Bound default EOS floors and seed
settings are serialized with max_digits10 so authentic restart inputs round-trip.

`width` must lie in1..nghost; `min_mach>=1`; `mach_max>1` bounds temperature
inversion only. `grad_floor`, `density_floor`, `temperature_floor` must be finite
and nonnegative; zero disables upstream floors. `div_limit<=0` is a signed strict
comparison `div(v)<div_limit`. Constructor default0 means compression required.
`event_capacity` defaults4096 and must be positive. It bounds completed/exit and
terminal suspended records; overflow aborts explicitly. Replay a previous paired
checkpoint with increased capacity and a new checkpoint basename. Capacity is
operational, so increasing it does not change the physical contract.

`trace=true` enables qualification observation/transfer CSVs. It is expensive and
optional; ordinary operation retains event evidence without full per-step traces.
Transfer diagnostics record source parent cell, initial density, six probabilities,
RNG draw, selected direction and actual timestep bin. They are not continuous paths.
Counters report samples, invalid observations, complete/partial/suspended episodes,
suppressed reentries, buffer high-water mark and diagnostic compute/I/O time.

## Restart and output

See schema.md for exact fields, flags, censoring and failure model. Keep the native
`.rst`, its `.shock.rankR` sidecar and the entire committed event directory together.
The native MC section saves all particle rows and seed-schedule progress/next tag;
sidecar version2 saves only typed encounter state and integrity metadata.

```sh
# Stop at a bounded generation, writing a native paired checkpoint.
/ABS/WORKTREE/build-debug/src/athena -i /ABS/WORKTREE/inputs/tests/shock_mc.athinput \
  time/nlim=25 output1/dt=1
# Use the emitted final .rst filename, retain prior event directory and lineage.
/ABS/WORKTREE/build-debug/src/athena -r /ABS/RUN/rst/shock_mc.00001.rst \
  time/nlim=-1 job/basename=shock_mc_resumed output1/dt=-1 \
  shock_tracking/directory=/ABS/RUN/shock_events
```

For actual MPI, use `single_file_per_rank=true` and the same rank count/partition.
Restart from the emitted `rst/rank_00000000/...rst`; each rank resolves its own
native payload and sidecar. Per-node layouts are rejected for every rank count.
Restart requires matching source/config/native ABI, valid native/sidecar checksums
and complete matching history through the checkpoint. Missing history is a failure,
not a fresh lineage. Old checkpoints can replay existing future chunks only when
bytes agree. Never disable histories when `history_present` is set in a restart.

## Deployment qualification

RUSTY (Flatiron) and Frontier are intended targets. Neither was accessed/tested.
Local multi-rank launch fails before AthenaK because PMIx runtime6.1.0 and PRRTE's
build5.0.9 are incompatible; logs are retained. No system installation was changed.

A CPU/MPI build template is `cmake -S . -B build-mpi -DAthena_ENABLE_MPI=ON
-DCMAKE_BUILD_TYPE=Release -DAthena_SHOCK_TESTS=ON`. Qualify native tests and the
same-partition suite on the target with its supported compiler/MPI environment.
The current feature intentionally rejects GPU/threaded modes until their gate is
qualified. Frontier HIP requires device/runtime tests, probability/state migration,
race/overflow checks and scientific review before removing that restriction.
Node/rank/cell/particle scale and acceptable compute/memory/I/O budget remain open.
Small local medians do not establish a production overhead budget.

See checklist.md, ledger.md and the independent review/evidence packet under
`evidence/shocks-mc/`. Physical reentries, weak-shock misses/bias and censored
outcomes must remain visible during scientific review. Feature branch was pushed at4aebdafb. The following local pilot correction is not yet pushed; no merge/release occurred.

## Laptop plots

See [laptop-tests.md](laptop-tests.md) for a one-command local runner, measured
costs and figure interpretation. The fresh28plot gallery and complete PDF are in
`evidence/shocks-mc/laptop-visuals/`; engineering passes and physical limits
are shown together. Core production source remains3de18d97/3feb3b.

## Decaying turbulence pilot - 2026-10-13

Reentry correction source099109b6/fingerprint1ac1558c clears stale exit evidence
after accepted reentry. Fresh native state/parity/integration/replay and held-out
checks passed; prior checkpoints require their matching old-source executable.
See evidence/shocks-mc/turbulence-pilot-20261013/README.md for the reproducible
32^3 periodic Mach3 decay with4096tracers, five plots and all event statistics.
Completed7123/open540episodes at CFL.30; completed7623/open530at CFL.15.
This bounded pilot is measured, with7.02percent joint numerical sensitivity.
It does not qualify driven turbulence, cluster performance or physical shock counts.
Old reports/PDFs retain their exact prior source identities; human review stays pending.

## Earlier opt-in forcing checks - 2026-10-13 (historical)

Historical corea03d621b/fingerprint93a861ab4638613fb94f475caeb5527e971bfdb1fc43b15baeee479a89997bc7.
Native forcing ordering/non-interference/sampling/replay checks pass for the local
experimental subset. This earlier implementation blocked forcing by default: explicit
`shock_tracking/experimental_forcing=true` permits only continuous zero-start
periodic3D untiled/unlocalized `normalization=accel_rms` driving. Static forcing
settings/seed and coupled stage/RNG sources are bound to the restart contract.
No duplicate dynamic state or event-schema change. Equivalent numeric spellings
or unused parameter edits conservatively reject restart. Old-source sidecars need
their matching old executable.

**Independent gas-energy oracle fails:** native RK2 forcing adds excess
3rho|a|²dt²/8in the shear fixture. Its energy update was left unchanged.
Forcing thermodynamic/scientific qualification remains open; do not interpret the
engineering checks as physical strength accuracy. See forcing-checks.md and
evidence/shocks-mc/forcing-checks-20261013/ for exact evidence and rerun commands.


## Forcing enabled - 2026-10-13

Current source `ad4225ed`, fingerprint
`f6195adafb7e6a554c23a7537053bb00e6ae8c618a10ff1951a487d6981cb70b`.
Native forcing now works with ordinary shock tracking; no experimental opt-in.
Both accel_rms and edot pass exact local on/off and native restart checks. Seven
short forcing variants pass compatibility tests; native fluid/driver/MC machinery
is reused unchanged. See [forcing-checks.md](forcing-checks.md) and the exact-code
packet `evidence/shocks-mc/forcing-enabled-20261013/`. The known native RK2
thermal bias is documented separately and does not block forcing. Scientific
strength calibration and RUSTY/Frontier/MPI/GPU qualification remain pending.
The new source/evidence commits are local, not pushed.


## Forced supersonic paper showcase - 2026-10-13

See `evidence/shocks-mc/supersonic-showcase-20261013/`:12-page PDF with adjacent
explanations,10figure families(PDF/SVG/PNG), episode/evolution CSV and selected
flow/mask arrays. Baseline64³/8192particles/initialMach5/nativeedot50/t.3 yields
57013complete+1212open, finalrmsMach1.316. Five simulation comparisons, including
a matched physical detector-span32³case, pass data/identity/non-interference checks.
The portable input and `tst/shocks-mc/supersonic_showcase.py` reproduce the packet.
No production machinery changed. These are evolving numerical encounters, with
source thermodynamic bias and physical/statistical/hardware gates still pending.
