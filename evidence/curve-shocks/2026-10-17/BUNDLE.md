# Exact-code engineering evidence bundle

Code commit: `6f1f1486ac671496fcb7b2f7812a87f9a8517f14`.

Bundle: `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-shock-histories-20261017`.

- `source.tar.gz`: committed source/tests/evidence at that code commit, including the final carrier-aware reader.
- `runtime/`: exact three macOS executables, build caches, flags, configuration and generated source fingerprints; binaries depend on their original platform/runtime.
- `raw-final.tar.gz`: complete latest native runs, checkpoints, histories, inputs, commands and logs. Original absolute paths are recorded; use fresh directories when rerunning.
- `historical-failures.tar.gz`: retained native development failures, separate from final evidence.
- `manifest.json`: sizes and SHA256 hashes; every listed file reverified after archiving. `verification.json` records that verification separately. Selected source archive files additionally compared byte-for-byte to the committed worktree.
- `implementation.md`, `event-restart-format.md`: run instructions, supported scope, event/restart documentation and open gates.

Status: **19/20 curve gates passed; required MC+curve AMR incomplete**. MC occupancy fails independently; its startup exclusion remains active. Physical AMR-width/allocation UserQuestions, compatible multi-rank checkpoint layout, actual RUSTY execution/resource budgets and human scientific approval remain pending. No push or merge into feature_dump; codex/shock-multiphase remains 2450e1c2.

Manifest SHA256: `ecab97a7c23190c354c2aeb5a87646971a37393c79a47be492206104e1218a07`.
