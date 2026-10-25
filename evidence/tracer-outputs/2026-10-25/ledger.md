# Existing-output integration

2026-10-25. User resumed output work after signed curve IDs. Starting branch
codex/shock-lagrangian, clean at a051de51. Existing worktree retained; original
branches preserved. Ponytail full and running-sim workflows apply. Existing
accepted physical fixtures/parameters are reused; no new scientific choices.

Accepted output targets: curve thermo and sampled particle/curve shock properties
through THP; mesh shock variables through the native mesh BIN writer. No separate
codec, container or new format selector. Full encounter journals and geometry
remain separate lossless products.

## Frozen first interface

U1: curve thermo uses literal THP v2, including its four signed integer slots
(cycle, node_id, curve_id, gid), standard time/x columns and named real fields.
Birth, immutable material label and sampled refinement level are additional
named fields. IDs never pass through Real. Existing particle THP bytes stay
unchanged. Retain current cycle/rank chunking and CommitFile immutable/replay
semantics; this work does not claim fewer files or append-performance gains.
A curve filename declares cycle, rank/rank-count and curve-bit split, including
empty chunks, using basename.output_id.cycleC.rankRofN.bitsB.thp. The existing
reader normalizes these records to its current curve analysis interface and
retains historical CSV reading. No user conversion step is required.

U2: optional named shock fields use the same native detector evidence and output
clock. Preserve unreduced region membership; distinguish pressure/temperature
jump estimates from fluid Mach. Missing estimates remain explicitly unavailable.
No episode identity/sequence is converted to floating point. Full every-step
events and restart state remain authoritative.

U3: variables=shock uses the existing mesh binary writer and derived-variable
registry. Reuse the synchronized particle detector array when available; the
curve-only periodic case uses the same DetectCell contract at output time.
No new time-integration/regrid ordering, detector criteria or physics support.

## Gates

- [x] U1 native THP curve round trips, exact controls and immutable replay (42 groups).
- [ ] U2 sampled fields equal raw traces, synchronized time/AMR/birth tests.
- [ ] U3 native BIN mesh fields agree with independent detector/particle samples.
- [ ] Historical-reader, empty populations, signed IDs, cadence and restart gates.
- [ ] Independent persistence/timing/simplicity review; apply concrete findings.
- [ ] Actual explained plots, exact source/configuration and commands retained.
- [ ] Coherent commits and sealed evidence; no push requested.
- [ ] MPI runtime/GPU/cluster performance and human scientific review pending.

## U1 completed

Serial/double builtin and prescribed builds pass. TO15 + CT18 + VB9 =42 groups
pass at build-lagrangian/evidence/thp-U1-{TO,CT,VB}-20261025. Complete physical
observations and events match controls; existing particle THP bytes are exact.
Signed negative IDs, zero-bit live AMR, active restart, empty populations, unchanged
cadence, corrupt-header/truncation rejection and immutable conflicting replay are
covered. Tests that targeted current CSV bytes now target native THP bytes; old
CSV decoding is retained and used by the independent stored controls. No tolerance
or physical parameter changed. Reuse of standard headers adds no full-data copy
to regular particle writes. Binary chunks retain the existing per-cycle/per-rank
file count; no throughput or file-count improvement is claimed.
