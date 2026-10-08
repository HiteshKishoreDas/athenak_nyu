# Exact-source MC engineering evidence

Current source is `3de18d97`, fingerprint3feb; fresh suites are in `../review3/`.
`../ponytail/` retains the previous80da constructor extraction. This packet documents
pre-cleanup78e5 source and remains intact.
All21physical measurement records match exactly after cleanup.

Final production commit `38c55570`; compiled source fingerprint
`78e5cb1d425daaa8ab2e06deaa40db6c8249de79cd1480c7d0b069fd1587f7c3`.
Release executable SHA-256
`ccb1c87b2e909869fefd1b7e4d07188b7642c857071b839ced52fbe8795e9fe9`.
The authoritative reports are `report_stylefinal.json`, `rejects_stylefinal.json`
and `replay_stylefinal.json`. They preserve exact commands/configuration/log/
checkpoint/sidecar/manifest hashes and raw run paths. These serial CPU checks
qualify engineering behavior; human scientific approval remains pending.

20 integration processes passed: fixed-seed uniform and planar on/off exact
thermo observations and byte-identical17-digit conserved slices; zero particles;
zero-event moving contact;128 completed M2 histories;128 IDs retained over four
local blocks with113 GID changes; initial and delayed birth each128 left-censored
partial events; synchronized exact sampled time/x/density. TRC11 independently
checks119808 transfer records: source parent and timestep bins exact, physical
uniform probability error max6.10623e-16 below predeclared1e-12 double roundoff,
logged draw/direction and realized wrapped cell hop consistent. It does not
invert the RNG or infer continuous physical particle trajectories.

Actual OUTSIDE/ACTIVE/PENDING checkpoints at generations1/61/68 resume with
exact canonical event content and exact same-partition particle observations.
Corrupt committed event prefix, missing/corrupt sidecar, explicit overflow and
unsupportedRK3 rejection pass. The separate12-process rejection packet passes
FOFC, duplicate schedule IDs, per-node restart, shearing-box, native-tag overflow,
oversized contract with recomputed FNV footer, future/truncated sidecar, corrupted
fluid binding and authentic-checkpoint changedEOSfloor/RNG contract rejection.
The four-process replay packet proves an old checkpoint collision preserves
existing paired checkpoint bytes and canonical event history; a new basename
allows deterministic replay. Existing legacy thermo output may append identical
terminal/replayed rows; its entire file is not claimed immutable.

Earlier records are retained as evidence:

- `report_v1.json`: startupSIGSEGV caused by nullable Hydro::psrc dereference,
  fixed by parent. Debugger launch hung and was stopped; no stack obtained.
- `report_v2.json`: analyzer used rho instead of canonical .thp field density.
  The analyzer lookup was corrected; no physical bound changed.
- `report_v3.json`: repaired Debug suite passed before final binding changes.
- `report_final.json`/`rejects_final.json`: strict authentic restart failed after
  EOSfloor binding exposed default serialization rounding1.1754943508222875e-38
  to1.17549e-38. Production now uses SetRealPrecise; equality was not relaxed.
- `report_reviewed.json`/`rejects_reviewed.json`: complete pass at earlier reviewed
  fingerprint26d4; final style-only source revision was rerun separately.

Reproduce from this checkout using fresh work directories:

```sh
python3 tst/shocks-mc/integration.py --executable build-release/src/athena \
  --workdir /private/tmp/shock-mc-integration-reproduce \
  --report evidence/shocks-mc/integration/reproduce.json
python3 tst/shocks-mc/rejects.py --executable build-release/src/athena \
  --workdir /private/tmp/shock-mc-rejects-reproduce \
  --report evidence/shocks-mc/integration/rejects_reproduce.json
python3 tst/shocks-mc/replay.py --executable build-release/src/athena \
  --workdir /private/tmp/shock-mc-replay-reproduce \
  --report evidence/shocks-mc/integration/replay_reproduce.json
```

Each process is limited to30s, integration suite300s. No multirank/GPU,
repartitioned restart, AMR, cooling, additional physics, FOFC orRK3 qualification
is claimed. Rusty node, Frontier backend/scale and performance budgets remain
owner decisions. Parent records actual active domain-exit/system/performance/
style/clean-checkout checks separately.
