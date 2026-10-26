# Living engineering checklist

Task: implementation plan in das-wiki, 2026-10-08. This branch replaces the
previous velocity-following experiment as the user's selected implementation.
Checked boxes refer to the named engineering evidence, not scientific approval.

- [x] M0.1–6: repository instructions, live branch inspection, isolated worktree,
  exact selected base/submodule, baseline debug build, frozen ShockIt snapshot.
- [x] M0.7: existing MC uniform-flow ensemble, local migration, zero-particle and
  same-partition restart baseline; inherited identical terminal thermo replay noted.
- [x] M0.8–9: CPU resources and bounded evidence retention recorded.
- [x] M1: experimental RK2 MC contract and independent fixtures frozen.
  User accepted branch/tracer; strength/timing/recovery release thresholds pending.
- [x] M2: native detector reused and rechecked against frozen ShockIt: 53 cases,
  exact masks, rtol=atol=1e-12; independent RH/smooth oracles separate from parity.
- [x] M3: final live timing/feature-off/boundary integration run recorded.
- [x] M4 baseline: TRC09 uniform-flow statistical checks, fixed-layout native restart.
- [x] M4 integrated: TRC10–11 parent/time-bin evidence, birth/exit/sampler checks.
- [x] M5 engineering pilot: state invariants and stochastic translating-shock histories recorded;
  one physical episode per passage remains a scientific acceptance gate.
- [x] M6 supported local scope: migration, bounded drains and all-phase replay/
  corruption/overflow checks. Cross-rank/storage-fault qualification remains open.
- [x] M7 supported local scope: debug/release/MPI-singleton, non-interference,
  inherited-only style findings and ASan/UBSan encounter-unit evidence.
- [ ] M7: true multi-rank CPU qualification; GPU/threaded/single precision qualification.
- [ ] M8: held-out physical measurements/performance; owner tolerance/budget acceptance.
- [x] M8.5/PHY16 bounded pilot:32^3 periodic Mach3 decay,4096MC particles,
  CFL.30/.15 statistics, exact non-interference and conserved fields measured.
- [ ] M8: turbulence scientific acceptance and production-scale cluster qualification.
- [x] M9 engineering evidence: reproducible docs, independent integrated review
  and exact-code bundle; clean archived full Release build plus state, 53-case parity
  and 20-process serial integration reproduction passes for42bcf12c (clean-reproduction/).
  Subsequent constructor extraction: independent exact review, three rebuilds and
  repeated core/persistence/physical checks recorded in ponytail/.
  Three subsequent complexity reviews implemented; fresh core/restart, all18boundary
  and13production exits, native ABI/payload and23748CSV equivalence checks in review3/.
- [ ] M9.6: human review of exact code/evidence. Scientific sign-off remains pending.
- [ ] M9.7–8: PR/release action if requested; R1 is not declared release-ready.
- [ ] M10: separate AMR, cooling/source, MHD and relativistic gates. Excluded here.

See ledger.md for commands/results/failures and README.md for supported behavior.

- [x] Requested laptop visual checks: fresh20engineering/21physical runs, native
  state check, actual128×4 snapshots and18sequential small-job benchmarks;28plots
  plus portablePDF/HTML/data/input/log packet. M8/M9scientific gates remain open.

- [x] Requested PDF explanations:36-page self-contained report with test methods,
  results, limits and28 original vector plots; all pages rendered/visually reviewed.

- [x] Requested plot-local explanations: setup, reading guide, measured result and
  limits on each of28plot pages;36pages rendered and reviewed; methods in appendix.

- [x] EXT06 ordinary native forcing allowed without opt-in:30processes, accel_rms
  and edot on/off and exact native replay, zero shear events,896 synchronized
  samples, seven forcing variants, exact activation/cutoff and physical exit oracles.
  Existing runtime/physics constraints remain; known native RK2 source bias documented.
- [ ] EXT06 physical forcing qualification: native RK2 gas-energy oracle FAILS;
  resolve/qualify source work, then scientific strength/recovery/calibration gate.

- [x] Requested paper-style supersonic showcase:5actual runs,64³/8192baseline,
  10vector figure families and12-page plot-local report. Strict history/counter/
  population identities; exact saved fluid/particle non-interference; native mask
  snapshot reconstruction; brackets/censoring and matched-span sensitivity retained.
- [ ] Paper physical acceptance remains separate: native source-work accuracy,
  independent physical shock calibration, ensembles/resolution/HPC and human approval.

- [x] Larger-run preflight: width8 NUM08 exact frozen ShockIt parity and NUM10
  independent Mach3 planar jump; same 1e-12 strength tolerance. Portable128^3
  input output-cadence smoke passes initial/cycle1/cycle2.
- [x] Larger-run checkpoint replay: five evolved steps,40,960 exact particle
  trace rows,1,449 exact canonical chunks; censored identities retained separately.
- [x] Requested128^3/t1.0 run: native exit0/generation6033; all49,430,528
  observations/all8192IDs per6034frames, invalid0;367,074complete/0partial/999open.
  Exact80-step native replay crosses forcing update;655,360trace rows/all1524
  chunks match. Independent final CSV/conserved-field audit and8-page visual QA pass.
- [ ] Larger-run physical acceptance: final rmsMach0.808576 from adiabatic heating;
  evolving flow, no plateau or stationary Mach5 claim. Source-work/calibration/
  ensembles/hardware/human gates remain open.

## 2026-10-24 compact curve identity follow-up

- [x] User-selected unsigned 32-bit curve IDs: high 8 bits curve, low 24 sequence;
  zero sentinel and deterministic exhaustion, checked config/checkpoint narrowing.
- [x] Current links and insertion ancestry use 32-bit IDs; existing-node identities
  survive insertion, AMR and same-code restarts. Shared encounter rules unchanged.
- [x] Historical geometry/thermo/event analysis retains original 64-bit identities;
  old checkpoints reject explicitly rather than silently reassigning identities.
- [x] Independent persistence/simplicity review; applied contract-version fix and
  added repaired-checksum shock-state alias regression.
- [x] 93 local check groups, additional reader/C++ IO/state checks and fresh plots.
- [ ] MPI runtime: build passes; PMIx/PRRTE ABI mismatch prevents local launch.
  RUSTY/Frontier/GPU, production budgets and scientific review remain pending.

Current evidence: `evidence/curve-ids/2026-10-24/ledger.md`.

## 2026-10-24 configurable curve-bit follow-up

- [x] `material_curves/curve_id_bits` defaults to 8, accepts 0–30 with checked
  integer parsing; zero keeps logical curve 1 and uses all 32 ID bits for sequence.
- [x] ID/link storage remains uint32; the allocator counter safely holds the
  exhausted full-width sentinel. Existing node/batch caps remain in force.
- [x] Geometry, thermo and journal metadata declare the split; all supplied
  readers reject missing/inconsistent layouts rather than guessing new IDs.
- [x] Split changes reject on restart; fixed-eight-bit geometry v3 compatibility
  and old analysis decoding are retained under their declared formats.
- [x] Independent persistence/simplicity review completed; 107 local groups pass,
  including zero-bit real-shock, AMR birth, exact restart and limit fixtures.
- [ ] MPI/GPU/cluster/runtime-budget and human scientific gates remain pending.

Current evidence: `evidence/curve-ids/2026-10-24/variable-ledger.md`.

### 2026-10-25 — signed32 curve IDs (user-directed)

- [x] Curve/node/link/parent storage signed32; full32-bit namespace retained.
- [x] Negative IDs valid; zero alone missing; unsigned arithmetic for packing/ancestry.
- [x] MPI_INT32_T identity transport; shared encounter algorithm/ABI unchanged.
- [x] Versioned interpretation in existing geometry/thermo/event file families.
- [x] Historical analysis preserves unsigned identities; mixed layouts reject.
- [x] 110 local Serial/double groups plus isolated readers/IO and state checks pass.
- [x] Positive-parent/negative-child insertion and -1 exhaustion fixtures retained.
- [x] Negative-ID live AMR events and active restart match exactly after ID remap.
- [x] Actual signed-ID/pressure/restart plots produced; evidence ledger in
  `evidence/curve-ids/2026-10-25/ledger.md`.
- [x] Independent signed persistence/simplicity review; checkpoint-key conversion fixed;
  optional legacy v4/B0 full-width restart regression added and passed.
- [x] Sealed signed-ID source/runtime/plots artifact: 14 figures, 35396 raw files
  verified, 175 native executable digests; artifact.json records the manifest.
- [ ] MPI runtime/GPU/cluster budgets and human scientific approval remain pending.
- [x] Curve THP writer and mesh shock BIN variables completed in the follow-up below.

### 2026-10-25 — native THP/BIN output integration

- [x] U1: curve thermo uses literal THP v2, exact signed integer slots, metadata and immutable chunks.
- [x] U2: optional six-field shock group shared by particle and curve thermo; full journals retained.
- [x] U3: variables=shock uses native BIN, active zones, standard precision options.
- [x] Zero-step restart rebuilds derived signals independently of encounter/journal advancement.
- [x] Header clocks, partial validity flags, empty metadata, replay and historical CSV checks.
- [x] 119 local Serial/double groups plus21 THP reader cases and additional standalone checks.
- [x] Independent timing/persistence/ponytail review findings applied and source-bound suites rerun.
- [x] Accepted short32³ forced box: exact sampled property parity and passivity.
- [x] Sealed 18-figure exact-code/runtime/raw-hash bundle; artifact.json records the manifest.
- [ ] MPI runtime, GPU/cluster production budgets and scientific sign-off remain pending.

### 2026-10-25 — unified curve CRV (supersedes curve THP writing)

- [x] One `.crv` holds typed geometry/connectivity/provenance, diagnostics, current
  mesh location and selected thermo/sampled shock fields at the same saved clock.
- [x] Particle THP v2 byte compatibility; historical readers and explicit ASCII retained.
- [x] No complete encounter embedding, cursor or event-buffer change; journals retained.
- [x] 50 local native groups and28 reader groups, independent timing/persistence review.
- [x] Actual cycle37 active witness, exact zero-step topology and field continuation.
- [x] Nine explained plots; exact source/config/commands and retained failures in
  `evidence/curve-binary/2026-10-25/`.
- [x] Sealed source/runtime/raw-hash gallery and artifact manifest: nine figures,
  28,317 raw-file hashes, 89 native-run digests; artifact.json records the manifest.
- [ ] Runtime MPI/GPU/cluster budgets and human scientific approval remain pending.
