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
- [ ] M8: supported-physics turbulence pilot, production-scale cluster qualification.
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
