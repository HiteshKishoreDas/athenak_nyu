# Independent reviews and applied fixes

Three read-only reviewers used ponytail-review; coordinator edited sequentially.
Reviewers did not execute tests or grant scientific approval.

- Timing reviewer: confirmed stage sampling/order and periodic RK2 increment algebra.
  Found Euler-predictor outflow deletion could fabricate completed exit evidence;
  added stage-aware fail-closed guard and native reversal regression. Found missing
  fine/coarse diagonal neighbor entries; added rare containing-leaf fallback and
  a native mixed-level test that explicitly logs all128 fallback resolutions.
  Follow-up source review found no races, duplicate sends or lost state.
- Persistence reviewer: traced every native real/int slot plus typed EncounterState
  through boundary compaction, MPI packing, AMR Alltoallv and native/paired restart.
  Confirmed transport discriminator rejects MC-to-continuous changes before arrays.
  Added remap caller to source fingerprint. Deleted unreachable MC shock-only
  instrumentation. Kept native MC RNG/draw/hops intact and tested tracking-off.
- Minimality/test reviewer: shared native infrastructure is appropriately reused;
  no speculative abstraction recommended. Strengthened RH test to require one event
  for every independently swept particle, rather than checking only emitted IDs.
  Added true3D z motion/wrapping. Restart reports hash the actual checkpoint and
  mark generated input unused instead of binding an unread input file.

Additional coordinator validation found the repeated-section input editor bug and
constructor-time versus resolved-AMR contract mismatch. Both fixed at shared entry
points. Failed runs and final passing evidence remain separate.
