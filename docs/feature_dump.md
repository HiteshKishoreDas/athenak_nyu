# feature_dump integration record

## 2026-10-17 — material curves

- Request: create `feature_dump` from `multiphase_turb`; merge
  `codex/multiphase-material-curves`; preserve the source branch.
- Worked directly in `/Users/hitesh/hitesh/git/athenak_nyu`. Initial checkout was clean;
  previous branch `codex/staircase-ordered-commits` remains intact.
- Refreshed `origin` with `git fetch origin`. Local and remote base/source tips agree.
- Base: `526a187f938e38634e89a88192ebe9731c4d36bb`.
- Source: `a0744a0d6dbc93364201bd8f5248e6be35b938b9`.
- Created with `git switch -c feature_dump multiphase_turb`.
- Merged with `git merge --no-ff codex/multiphase-material-curves`.
  Explicit merge records the feature integration and its original base.
- Merge commit: `b13ca7d65fc821653d04cf006b6bfcf6e4b61e9e`. No conflicts or code edits.
- Verified merge parents and exact tree equality with the source feature tip
  before adding this integration record. Source branch remains at its original tip.
- Check: `PYTHONDONTWRITEBYTECODE=1 python3
  tst/scripts/material_curves_io/test_reader.py` — all 10 tests passed.
- Existing numerical evidence: [material curves](../doc/material_curves.md) and
  [qualification record](../doc/material-curves-evidence/README.md).
  No new build, simulation, GPU qualification, push or PR performed in this merge.

## 2026-10-17 — publication

- User authorized publishing `feature_dump` to `origin`.
- Integration/documentation tip before publication: `430c21e9`.
- Publish with `git push --set-upstream origin feature_dump`;
  verify the remote branch SHA against local `HEAD` afterward.
- Normal branch push; preserve `codex/multiphase-material-curves` unchanged.
- Existing 10 reader tests passed during integration; publication adds no code changes.

## 2026-10-17 — shock-encounter merge into feature_dump (validated)

- User authorized merging all of `codex/shock-multiphase` into `feature_dump`,
  explicitly asking that material-curve interactions be assessed.
- Initial main checkout clean, target `cea1dafcf6830bfecba3e2aea5622b5c943a761c`;
  source `2450e1c21b0a97c0f1635da75a20701c7e0852a4` (12 local commits ahead of tracking ref).
  Shared base `526a187f938e38634e89a88192ebe9731c4d36bb`. No fetch/push this step.
- Isolated temporary integration branch `codex/merge-shocks-feature-dump` at
  `/private/tmp/athenak-feature-dump-shocks-20261017`; main target stays untouched
  until testing. Source branch/worktree preserved, including unrelated untracked AGENTS.md.
- `git merge --no-ff --no-commit codex/shock-multiphase` found two conflicts:
  `src/CMakeLists.txt` and `src/main.cpp`. Preserved both compile source lists
  and both independently conditional restart sidecar loaders. Driver/restart
  output edits auto-merged; reviewed ordering explicitly.
- First resolution helper's marker assertion matched AthenaK's decorative
  equals-sign comment banner; stopped before writing main.cpp. Corrected check
  to line-start conflict markers. No behavioral workaround or tolerance change.
- Existing curve constructor rejects particles/forcing; curve Node and MC tracer
  state are separate. This merge preserves that supported envelope. Curve shock
  sampling/history is a future qualification, not inferred from coexisting source.
  Newly inserted geometry nodes must not receive invented past shock histories.
- Build uses local exact Kokkos `08ceff92bcf3a828844480bc1e6137eb74028517` through
  a temporary worktree symlink; no dependency install/network. `git submodule status`
  also reports the preexisting unmapped sundials gitlink; no sundials change.
- Local Release built-in-pgen + existing supplied-state shock tests are being
  configured out of source. Existing curve prescribed-pgen regression gets its
  own build. Preserve original scientific/hardware limits and frozen tolerances.


### Final interaction fix and validation

The pre-fix candidate silently accepted `shock_tracking/enabled=true` with curves
but no `<particles>` block and wrote no shock history. Retained this failing input,
log and binary hash. Added one seven-line guard in `MeshBlockPack::AddPhysics`
after particle assembly. Enabled tracking now requires particles explicitly;
existing particle configuration enforces `lagrangian_mc`. This prevents silently
misreading curve output as shock history. No curve/MC integrator or physics changed.

Resolved build-list conflict by compiling both features. Resolved restart-loader
conflict by preserving both independently conditional loaders. Auto-merged driver
and checkpoint-write changes retain both parents' timing and disabled paths.
Three independent timing/integration/persistence reviews found no production blocker;
one reviewer requested the new guard. Review test suggestions were applied: require
512 actual shock observations and label unused restart fixture input accurately.

Checks on final compiled fingerprint
`48c53dc5c167604269d121412d4316215f09be8224e053d98af2eae9908a2b7b`:

- Both local Serial/double Release builds pass: built-in problems and
  `PROBLEM=material_curve_test`; encounter state checks and 53 exact-mask parity
  fixtures pass, unchanged floating tolerance 1e-12.
- MC integration: 20 native processes; immutable paired-checkpoint replay: four
  processes. Supported timing, migration, identities, bounded output and restart
  checks pass. New 11-process interaction test passes (eight successes, three
  expected rejections), including exact float64 fluid passivity, disabled-feature
  equivalence, native paired curve restart, missing-particle and curve/MC rejection.
- Curve reader: ten tests; production core: 31 records; regression: seven records;
  unsupported configurations: 18 processes; native C++ serial persistence checks
  pass. No MPI/GPU/AMR runtime requalification performed in this merge.
- All 21 held-out physical measurements exactly match the prior source's
  measurement dictionaries. This is preserved measured behavior, not scientific
  acceptance. Thirty forcing processes pass compatibility/on-off/replay checks;
  inherited native RK2 gas-energy physical oracle remains FAILED, unchanged.
- Old-source shock checkpoint is explicitly rejected by the changed exact-code
  contract. Its original files remain byte-unchanged. Continue old runs with their
  matching old executable; do not remove the source branch or relabel old evidence.

A first interaction helper wrongly treated the native reader's list of meshblock
arrays as a single ndarray. Fixed representation handling with `np.asarray`, keeping
exact float64 array equality and all gates. Original failed report/log retained.
The earlier 53-case parity and existing feature tests were unchanged. Nonfatal
font-cache warnings during the forcing plot helper did not affect its passed checks.

### Material-curve shocks: what the merge does and does not provide

Both features are available in one executable. Existing curve constructor still
rejects particles and continuous forcing. Curve nodes have ordered persistent IDs,
RK2 mesh-velocity sampling and adaptive insertion; MC tracers have mass-flux motion
and their own encounter state. No node-level shock history is added by this merge.

A future curve-history extension can reuse the detector kernel and encounter engine,
but needs a curve-owned sampling/state adapter and persistence contract. Sample at
completed synchronized fluid/curve times; preserve timing brackets and quality flags.
Key histories by stable node identity. Inserted geometry nodes start at their actual
birth time; never copy or interpolate fabricated earlier encounters. Qualify uniform
periodic hydro first; AMR, forcing and other curve physics require separate gates.
Scientific calibration and human approval remain pending for the shock feature.

Exact inputs, commands, source/binary hashes, build logs, failure fixtures, independent
reviews and restart evidence are in
[evidence/feature-dump-shocks-20261017](../evidence/feature-dump-shocks-20261017/README.md).
Full 38,724-file native archive and both exact executables are retained at
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/feature-dump-merge-20261017`.
All changes through source tip `2450e1c2` are included, not just its last pushed tip.
`codex/shock-multiphase` must remain at its original tip. No push authorized in this step.


### Merge completed

Merge commit `0b063c85e9f3f7f0857647ce944cc0b686b2f362` has both original tips
as parents: `cea1dafcf6830bfecba3e2aea5622b5c943a761c` and
`2450e1c21b0a97c0f1635da75a20701c7e0852a4`. Updated the main checkout's
`feature_dump` with `git merge --ff-only codex/merge-shocks-feature-dump`.
Verified both histories are ancestors and every tested source-file hash matches
this checkout. Main checkout was clean after fast-forward. The source branch
`codex/shock-multiphase` remains exactly `2450e1c2`; its worktree is preserved.
No push performed. This final completion entry is documentation/evidence metadata only.
