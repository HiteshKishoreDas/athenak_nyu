# Shared particle/material-curve thermo histories — 2026-10-22

Request: make material-curve histories and regular tracer thermo histories work
together. Branch codex/shock-lagrangian at 54800cff; existing documentation-only
dirty change in docs/curve-shock-showcase-20261021.md preserved. Primary checkout
feature_dump and codex/shock-multiphase remain untouched.

## Implementation decisions

- Extend existing prtcl_thermo_history output with source=particles|curves|both;
  default particles preserves existing behavior/schema.
- Reuse ParseTracerFieldList/EvaluateTracerFieldHost; same variables and native
  output scheduler for both populations. No shock-tracking dependency for thermo.
- Reuse curve leaf-owner request routing. Containing-cell sampling, matching
  particle diagnostics; continuous velocity interpolation/advection unchanged.
- Preserve exact uint64 node identity, curve identity, material label, birth,
  lifted positions and sampled leaf metadata. Never copy parent histories.
- Companion per-rank immutable CSV chunks; byte-identical final/replay accepted,
  conflicting output rejected using existing durable commit helper. No new
  particle binary or curve checkpoint schema. Restart resumes from persisted
  geometry/IDs and native output scheduler; this is sampled output, not state.
- Bound returned field batches by existing max_nodes and MPI int-count limits.
  Host evaluation matches existing tracer implementation; GPU/performance gates
  remain separate. Curves retain Newtonian periodic RK2/double support envelope.

## Gates

- [x] Shared output source parsing, curve-only and coexistence with shocks off.
- [x] Independent thermo/clock/position oracles; shared canonical field aliases.
- [x] Adaptive insertion and live AMR, exact identity/birth and restart checks.
- [x] Shock+thermo simultaneous clocks/values; passivity and legacy regression.
- [x] Local plots with explanations; reproduce commands, exact source/build hashes.
- [x] Commit coherent implementation/evidence; scientific sign-off pending.

No numerical tolerances or physical defaults will be adjusted to obtain passes.
Existing accepted fixtures supply physical parameters; no new science setup needed.

## Final supported engineering gates

All implementation decisions above are implemented. Eighteen CT checks pass in
`build-lagrangian/evidence/curve-thermo-reviewed-code/report.json`; fourteen LT
checks pass in `build-lagrangian/evidence/thermo-reviewed-LT/report.json`; twenty
CSH checks pass in `build-lagrangian/evidence/thermo-reviewed-CSH/report.json`.
The exact detector still passes all 53 parity cases at unchanged rtol=atol1e-12.
Existing standalone encounter and adversarial particle-reader self-checks pass.
These are local Serial/double results, not scientific approval or cluster tests.

- [x] Curve-only, particles-only/default, both, and MC+curve thermo with shocks off.
- [x] Eleven shared canonical fields; alias/order/dedup; passive scalar0 availability.
- [x] Independent prescribed acceleration and spatial containing-cell oracles.
- [x] Live refine/coarsen: leaf count \(4\to16\to4\), fourteen positive-time inserted nodes.
- [x] Exact AMR restart records, no duplicate new curve identities.
- [x] Simultaneous shock/thermo values:1728 curve observations, exact rho/P/T parity.
- [x] `dcycle` and `dt`; unchanged legacy particle bytes and final curve state.
- [x] Six explicit unavailable-source/field rejections; strict reader corruptions.
- [x] Immutable replay/conflict guards, exact IDs above2^53, empty populations.
- [x] Forced \(32^3\) coexistence:4096 particles, 16 nodes, 16 common sample clocks;
  65536 particle and 256 node observations. Particle bytes, curve geometry and
  curve shock events equal the legacy-output control.
- [x] CT16 off/on-cadence stop/resume; CT17 overdue thresholds and zero-step resume.
- [x] Independent timing and persistence/minimalism reviews (see reviews.md).
- [x] Explained native-data plots and reproducible run/analysis documentation.
- [ ] RUSTY/Frontier deployment, genuine multi-rank, GPU/threaded, production budgets.
- [ ] Human scientific sign-off; inherited forced thermal issue remains open.

Maximum independent thermo error 3.552713678800501e-15 against fixed 2e-12 budget.
CT11 checks 5039 same-cell comparisons. Curve/particle clocks match exactly.
CT01 preserves node ID 9223372032559808513 and all other exact uint64 identities.
CT06 compares 1179 curve and 5760 particle AMR records exactly across continuation.

## Failures retained, why changes were necessary

1. Initial attempt passed9 of10 checks; CT07 incorrectly compared lifted curve
   coordinates with wrapped shock coordinates bitwise. Corrected the test to
   compare periodic displacement using the already fixed 2e-12 budget. Correct
   final discrepancy6.94e-18. Initial exploratory attempt also used the prescribed
   executable for the planar slot; accepted campaigns use the built-in shock code
   and require actual shock membership/events.
2. CT16 found a genuine inherited thermo clock bug. At t=.225, terminal output
   advanced a .05 cadence from .2 to .25; restart then skipped the t=.25 sample.
   Raw failure: `build-lagrangian/evidence/curve-thermo-dt-repro/`.
3. Independent timing review found a second edge: same-cycle final replay could
   advance another overdue threshold when output cadence is shorter than a fluid
   timestep. Kept native catch-up semantics; persisted `last_cycle` in the existing
   output parameter block prevents advancing twice. No driver/solver change.
4. First CT17 expected a *final.rst filename, but native restarts are numbered.
   Corrected its fixture to compare the last two numbered files. The retained
   `curve-thermo-overdue-repro/overdue-native-diagnosis.json` confirms genuine
   pre-fix last_time=.09 versus.095. Final CT17 gives.09 in both, and identical
   clock after zero-step restart. Initial fixture StopIteration remains recorded.
5. A parity help command used system Python, which lacks NumPy. Actual parity
   uses the already installed ShockIt environment; no dependency was installed.

No tolerances, initial fluid states, solver settings or detector/encounter rules
were changed to obtain passes. The shared thermo clock fix intentionally corrects
continuation cadence for ordinary particles as well as curves; .thp schema stays v2.
Old checkpoints without last_cycle do not guarantee zero-step overdue-clock
idempotence. New checkpoints are qualified by the native CT17 check.

## Plots and next action

Final viewer (five figure families, explanations adjacent):
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-thermo-synergy-20261022/final-reviewed/index.html`.
Earlier viewer folders are retained drafts; final-reviewed binds final code.
Source/configuration/binary bindings and complete commands are sealed alongside
that viewer. Next external gate: human scientific review and RUSTY CPU/MPI
qualification, with node/particle scale and performance budget still unspecified.

## Sealed exact-code bundle

Implementation commit `b0d38d28ba6d290bc0c2d1b29b1e7fb600d0df79` binds all 40 compiled fingerprint sources.
Archive contains 2211 source files including Kokkos and complete local fixtures,
plus 25236 verified raw files across final campaigns and both timing failure reproducers.
`bundle-verification.json` independently verifies every archive source member and
all artifact file hashes. Manifest SHA256:
`b76d28f162cc6be324a56257f7e59b5c225879d1f01fa1af82e47c026173cd4c`.
Build fingerprints and exact binary hashes are in that JSON and the sealed packet.
Earlier unsealed draft directories remain for traceability. No source changes or
new performance/scientific claims follow from packaging.
