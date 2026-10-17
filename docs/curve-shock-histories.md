# Material-curve shock histories — requirements

## 2026-10-17: initial questions (historical)

User requested implementation, with all questions asked before coding. No feature
implementation has started. Current feature_dump base is f27e40a0. Preserve
codex/shock-multiphase unchanged.

Accepted birth rule: newly inserted nodes start their own observation/encounter
history at insertion time. Never copy parents' encounters or fabricate pre-birth
observations. Stable node identity, ancestry and birth time remain distinct.
The treatment of a node already inside a shock at birth is pending below.

Recommended approach: reuse the existing ShockIt-contract detector and encounter
engine, preserving bracketed times and validity/censoring evidence. Existing curve
advection and topology stay intact. Store complete history state through output
and restart. Proposed settings and physics support below are not accepted yet.

Questions for the user:

1. Physics: must the first version include forced supersonic turbulence? Recommend
   uniform ideal-gas Newtonian hydro with forcing included; qualify AMR, cooling
   and MHD separately unless required now.
2. Coexistence: must MC particle histories and curve histories work in the same
   run? Recommend yes, with distinct identities and passive behavior checks.
3. Product: independent per-node histories only, or whole-curve diagnostics too
   (for example instantaneous arc-length fraction in shocks)? Recommend per-node
   history first, curve summaries derived during analysis.
4. Born inside: retain a clearly incomplete/left-censored encounter starting at
   birth, or wait for exit before recording later encounters? Recommend retain
   the incomplete encounter without inventing its entry time.
5. Output: save each node's every-step shock flag as well as event histories by
   default, or make detailed traces optional? Recommend always save events and
   current restart state; every-step traces opt-in. Detection itself remains
   synchronized every complete step and at node birth.
6. Detection: reuse existing MC detector thresholds, strength evidence and
   encounter/rearming rules initially, or request changes? Recommend reuse as
   explicitly experimental defaults, with independent curve-path tests.
7. Existing restarts: must users be able to enable tracking when restarting an
   older curve-only checkpoint? Recommend fresh tracking runs first; support for
   starting new history at restart requires an explicit provenance policy.
8. Deployment: laptop CPU first, RUSTY CPU/MPI, or Frontier GPU required for the
   first supported version? Approximate curves/total nodes and acceptable runtime,
   memory and I/O costs are needed for deployment-scale qualification. Recommend
   laptop engineering gates followed by RUSTY MPI; Frontier remains untested
   until actual hardware evidence exists.

Retain the existing native RK2 forcing thermal-oracle failure as an explicit
physical limitation. Engineering implementation and tests do not supply human
scientific approval. No new physics values, floors or acceptance tolerances
were chosen in this requirements discussion.


## Accepted answers — 2026-10-17

- Include AMR initially, alongside the forced uniform-hydro scope discussed.
- MC particle and curve-node histories must coexist in one simulation.
- Per-node histories first; whole-curve summaries are derived in analysis.
- Initial/inserted nodes already inside shocks retain incomplete left-censored
  encounters beginning at initial observation/insertion. Never copy parent events.
- Always retain events and complete restart state; per-step traces are optional.
  Observe at every completed timestep and at actual insertion.
- Reuse existing detector thresholds, strength evidence and encounter/rearming
  rules initially; independently test curve crossings and AMR effects.
- Fresh tracking runs first; do not silently activate tracking from old curve-only
  checkpoints or invent earlier histories.
- Laptop checks first, followed by RUSTY CPU/MPI qualification. Node counts and
  resource budgets are not yet supplied. Frontier GPU remains outside first gate.

One AMR numerical choice is pending: width in local leaf cells (recommended,
matching existing MC definition) versus fixed physical distance. Continue storage,
birth and uniform-path work independently; do not silently choose AMR scale.

Implementation gates: freeze consumers/state interfaces; retain failing fixtures;
implement node observation/birth tracking; preserve complete curve/MC state through
regrid/migration/checkpoints; qualify uniform and actual refine/coarsen cases,
coexistence, passivity and replay with independent physical references. Existing
MC history rejects AMR and must be extended rather than bypassed. Human scientific
approval and actual RUSTY/hardware performance remain separate.
