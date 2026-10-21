# Continuous Lagrangian shock-carrier qualification

35 native checks pass: 14 standalone-particle,20 curve-history,1 native forcing.
53 exact ShockIt parity cases pass. Encounter core and strict readers pass.
MPI compilation and single-rank runtime pass; multi-rank launcher fails before
AthenaK (PMIx runtime0x60100 vs PRRTE build0x50009). No cluster/GPU science claim.

Implementation, accepted scope, limits, failures and reproduction commands:
`docs/shock-lagrangian.md` (repository root).

Full source/binary/configuration/raw evidence and seven explained plots:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/lagrangian-shock-carriers-20261020-final/`.
`manifest.json` there binds the implementation commit and archive/binary hashes.

All particle shock histories now require `particle_type=lagrangian_tracer` and
`pusher=lagrangian_tracer`. Legacy MC remains available with tracking disabled.
RK predictor outflow fails explicitly; full outflow trajectory support unqualified.
Scientific approval remains pending until the human supplies it.
