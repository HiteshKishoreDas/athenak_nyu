# Three Ponytail reviews and integrated result

User explicitly requested three differently focused subagents and implementation of
suggestions. All read applicable instructions and ponytail-review; no agent edited code.
Root owned files and integrated sequentially against8b453ed7.

- Numerical/state reviewer: remove duplicate representative state (-2lines), reuse
  ConfigValid (-3); retain all other detector/state diagnostics with real users.
- Native integration reviewer: one native restart header (-11), writer tautology (-1),
  redundant exit-index sort (-1); keep guarded migration, native transport and task graph.
- Persistence reviewer: same-file ShockError, shared PutEvent, guarded contiguous state
  append; retain immutable fsync/link, contract/checksum/replay and strict Python reader.

All findings implemented. Overlapping ConfigValid findings consolidated. Actual EOS
assignment order stays unchanged; no early gamma move. Forty net production lines
removed; fixture changes simply fill the four remaining state observations. No numerical
threshold, event schema or fluid/tracer algorithm was changed. Ten error messages gain
shared fatal prefix and MPI abort. Record format stays124columns; representative is first.

Each reviewer independently checked integrated changes. State size drops1000to816bytes
on tested ABI; event stays944, native MC header48 with all offsets unchanged. Existing
size/hash checks reject older sidecars. These are engineering reviews, not human science
sign-off. Last MC change only removed a blank line; exact reconstruction hash is recorded.

Revalidation: three builds/state tests;53parity;20integration;12reject;4replay;
18boundary fixtures;13actual domain exits; actual prior-layout rejection. All21physical
measurements match prior results, including known extra episodes. All23748heldout CSVs
are byte-identical after substituting only the implementation fingerprint. Native particle
restart tails are identical for OUTSIDE/ACTIVE/PENDING. No thresholds were loosened.
Source style check has10inherited MC warnings and0introduced findings. Initial reference
style comparison lacked repository CPPLINT.cfg; rerun with exact original config passed.

Reports retain exact commands/config/input/log/executable identities. Complete raw run
outputs remain in their named temporary directories; compact reports/inputs/logs are
committed. No new true multi-rank/GPU/cluster/performance/scientific approval is claimed.
