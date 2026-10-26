# Independent review: sampled versus complete binary encounters

Native timing/persistence reviewer performed read-only inspection, including
source-frozen integration and final21-check report. Python reader implementation
was delegated with explicit ownership; the native reviewer independently checked
its chain/query/magic integration. Both used the requested ponytail review scope.

Concrete findings applied:

- Snapshot cache generations are distinct from encounter generations; sampled
  outputs never call ObserveShocks, allocate histories or update episodes.
- History-only pusher/seed/restart restrictions apply to complete mode; MC may
  sample point properties without continuous-history claims.
- Queue records retain original emission clocks and survive particle destruction;
  checkpoint packets are inside fluid-bound sidecar checksums. Cumulative capacity
  is bounded and failure/recovery is explicit.
- Finalize must flush all event outputs before terminal views and paired restart.
  Earlier ordering changed the terminal prefix on zero-step restart; fixed and tested.
- Scheduled and terminal checkpoints can share a clock. Test terminal replay chooses
  highest file number. Ordinary checkpoints require their exact saved committed
  prefix; native guards reject all other own-stream completed tails/branches.
- Deterministic encoding sorts buffered records. Prefix verification walks every
  committed predecessor, not only the newest file. Same-cycle replay checks the
  saved file hash and empty queue before doing nothing.
- Committed domain-exit events legitimately carry right_censored + domain_exit;
  binary readers preserve them, while rejecting terminal-only views as committed.
- Custom history directories can share extensions/folders with sampled data.
  Known magics dispatch to the correct readers; unknown/corrupt candidates fail.
- Sample-only local lineage fallback must not silently provide the other carrier's
  required history lineage. Startup summary distinguishes snapshot thermo from
  event-only output and labels complete encounters explicitly.

Native review final: no remaining blocker, all21 native checks pass including
mixed roots. Reader author independently ran final47 checks on actual native
packets and retained corrupted fixtures. The final query-coverage fixture removes
terminal views when deleting the latest completed batch, so it reaches the intended
coverage guard; separate terminal-ancestor corruptions remain tested. No production
reader change was needed for that fixture correction.

SM20's destroyed-ID/domain-exit test is protocol evidence only. No live outflow,
MPI runtime, GPU or scientific approval was supplied by review. Lean implementation
reuses the existing output scheduler, detector/encounter engine, field sampler and
immutable writer; no external dependency was added.

Final timing validation also rejects invalid-fluid endpoints from a generation after emission, while preserving unknown-entry sentinel(-1,0). All47 event-reader cases and final21 native groups pass on the exact final source.
