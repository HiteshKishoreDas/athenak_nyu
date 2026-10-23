# Tracer output alignment — living evidence ledger

2026-10-23. Branch `codex/shock-lagrangian`, isolated worktree
`/Users/hitesh/hitesh/git/athenak_nyu/.worktrees/shock-lagrangian`.
Starting HEAD `66e436e130f39399a201539b3cc108d3cacc126c`; main checkout and
`codex/shock-multiphase` preserved. No push authorized for this increment.

## Scope and gate checklist

- [x] Inspect upstream/pre-feature output conventions and existing shared thermo.
  Cached upstream `5de631377124b2696201425c2675f8d46b45d7a2`; multiphase_turb
  `526a187f938e38634e89a88192ebe9731c4d36bb`. No live upstream freshness claim.
- [x] Freeze one optional native diagnostic table codec, exact widths/names,
  identities and independent sampling/journal clocks. Legacy remains default.
- [x] Reuse native output blocks for geometry and existing thermo output class.
  Keep event drain every completed step; no numerical/pusher/state-machine edits.
- [x] Preserve particle THP bytes and checkpoint state layouts. Bind journal
  representation into the existing restart contract/digest checks.
- [x] Independent persistence and reader/clock reviews; apply both correctness
  fixes and removal of discarded legacy formatting. See `reviews.md`.
- [x] Final local Serial/double builds and 82 native check groups pass:
  21 output, 18 CT, 20 CSH, 14 LT and 9 startup summary groups.
- [x] Direct encounter state executable and offline curve IO schema/adversarial
  reader checks pass. `git diff --check` and changed Python compilation pass.
- [x] Render six explained figures from actual stored records; inspect all six.
- [x] Seal source/runtime/raw-output evidence bundle after the implementation commit.
- [ ] Human scientific approval; RUSTY/Frontier/MPI/GPU and production budgets.
  These remain separate existing gates, not a reason to claim this slice untested.

## Commands and actual final outputs

All commands run from the worktree above. Both builds use the existing local
Serial CPU/double caches and initialized Kokkos. Exact caches, generated source
headers, executable hashes and recomputed source lists enter the evidence bundle.

```sh
cmake -S . -B build-lagrangian
cmake --build build-lagrangian -j8
cmake -S . -B build-lagrangian-prescribed
cmake --build build-lagrangian-prescribed -j8
python3 tst/tracer-output/check.py --executable build-lagrangian/src/athena --prescribed-executable build-lagrangian-prescribed/src/athena --workdir build-lagrangian/evidence/native-output-qualified-final-20261023
python3 tst/curve-thermo/acceptance.py --executable build-lagrangian/src/athena --prescribed-executable build-lagrangian-prescribed/src/athena --workdir build-lagrangian/evidence/native-output-frozen-CT
python3 tst/shocks-curves/acceptance.py --executable build-lagrangian/src/athena --prescribed-executable build-lagrangian-prescribed/src/athena --workdir build-lagrangian/evidence/native-output-frozen-CSH
python3 tst/shocks-lagrangian/acceptance.py --executable build-lagrangian/src/athena --prescribed-executable build-lagrangian-prescribed/src/athena --workdir build-lagrangian/evidence/native-output-frozen-LT
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-setup/check.py --executable build-lagrangian-prescribed/src/athena --baseline /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-thermo-synergy-20261022/final-reviewed/runtime/build-lagrangian-prescribed/athena --workdir build-lagrangian-prescribed/evidence/native-output-frozen-setup
build-lagrangian/shock_encounter_test
python3 tst/shocks-curves/io_reader.py
```

The raw runners record per-run exact commands, source/binary/input/checkpoint/log
hashes, durations and expected rejection exits. Final reports are copied here.

Final plot command:

```sh
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-output/visuals.py --run build-lagrangian/evidence/native-output-qualified-final-20261023 --output /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/tracer-native-outputs-20261023/qualified
```

## Measurements, fixed oracles and limits

The planar pair has 13,824 unique particle thermo samples and 1,728 node samples
over 108 clocks, plus 15 emitted particle and 16 emitted curve events. All
corresponding thermo/geometry/event/trace values, flags and identities are exact
between modes; NaNs compare by their explicit invalid/masked meaning. THP bytes
and all planar conserved-fluid text tables are byte equal.

The AMR pair has 5,760 particle and 1,179 node samples. Leaves change 4→16→4;
14 nodes start at positive insertion time and have no earlier saved history.
The independent manufactured thermo oracle has maximum error zero, within the
unchanged 2e-12 CT/LT budget. The large-ID fixture spans
18446744069414584321–18446744069414584352, retained exactly. Hand-packed opposite
endian/Real-width fixtures qualify the file reader only, not execution on those
hardware/precision targets.

Native small journals are larger in this fixture: per-chunk schema/provenance
cost dominates sparse/empty events. Figure 4 reports actual bytes including
unchanged THP, and one local run wall time per mode. No production performance,
uncertainty/scaling or memory-budget pass is inferred from those measurements.

## Retained failures and corrections

- First Python import smoke hit an indentation error introduced while connecting
  the curve reader. Corrected; changed Python modules compile and all native/
  legacy reader checks now pass. No native simulation had launched at that point.
- `build-lagrangian/evidence/native-output-smoke-20261023/` retains the initial
  native input/log rejected by BaseTypeOutput variable validation. Root cause:
  new specialized `curves` output needed the same early-return path as existing
  specialized history outputs. Fixed in the shared base; TO01/TO10/TO16/TO20
  cover normal, cadence, missing-module and shared-directory behavior. Early
  development executable/source bytes were not sealed; this is retained failure
  evidence, not an exact-code failure reproducer.
- `ctest -R shock_` discovered no registered CTest tests; it is not counted as a
  test pass. Ran the actual encounter executable and standalone IO check instead.
- `native-output-final-setup/` ran all nine native groups successfully but its
  plotting tail failed under system Python (Matplotlib unavailable). Retained the
  report and log; reran with the existing scientific environment at
  `native-output-frozen-setup/`, which exits zero and produces its plot.
- Earlier native-output attempt/reviewed/final and pre-freeze CT/CSH/LT directories
  remain in ignored build evidence. Final qualification uses only frozen reports.
- Reviewers' two valid-but-incomplete reader reproductions remain in TO12, with
  corrupt/alias fixtures under the frozen native test directory. No tolerances
  were loosened and no detector thresholds/physics values were retuned.

### Final coupled AMR addition

TO21 reuses the accepted CSH live-shock refinement/coarsening input with both
continuous particles and curves. It finds 17 positive-time inside births and 22
outside births. A checkpoint at an actual inside birth (cycle 13) resumes thermo,
geometry, emitted events and censored terminal views exactly. First observations
have sequence zero, actual birth time and the correct incomplete-entry/clean-start
flags. The sixth figure shows those insertion observations and restarted pressure.

The first TO21 attempt in `native-output-qualified-20261023/` failed a harness
assumption that births must appear in separate `birthtrace` files. Source tracing
confirmed the unchanged native behavior: completed-step adaptation is observed
in `trace`; stage-one adaptation uses `birthtrace`. The fixture creates 39 births
in completed-step traces and zero separate birthtrace files. Corrected the
harness to audit the earliest real observation across both streams, as the
existing CSH oracle does; added exact sequence/censoring checks. No simulation
parameter or tolerance changed. The failed report/input/raw records remain.
Native positive-time *separate* birthtrace emission is not claimed tested here.

### Evidence assembly control binding

The first assembler stopped before sealing because its current-build assertion
also encountered the intentionally historical startup passivity baseline. All
other native runs matched the frozen current executable hashes. The assembler
now binds that one named baseline separately, checks its exact bytes against
the prior sealed source binding, and retains its runtime/provenance. Current
runs still must match a freshly recomputed current-source fingerprint. This is
an evidence-accounting correction, not a relaxed numerical or runtime gate.

## Sealed exact-code bundle

Implementation commit: `a8e19214`. Bundle source commit:
`a91b01c23e9649d9e5f08856a4b53f967ff2d2c3`. No numerical/runtime C++ changes
followed the frozen builds and final checks.

Artifact directory:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/tracer-native-outputs-20261023/qualified`

Manifest SHA256:
`7ba7763a498e8cf97a4b89794a1f49cf824045213741f618c5229bf19d9a681e`.
The assembler decompressed and individually verified all 38,387 raw files and
2,217 source files, recomputed both compiled source fingerprints, checked every
current run executable hash, and bound the single historical startup control
separately. A subsequent independent pass verified every manifest file hash.
Six figures, full native raw outputs (including retained failures), current and
historical-control runtimes, build configuration, exact source, reviews and
usage documentation are present. `artifact.json` records counts and binding.

Human scientific approval remains pending. RUSTY/Frontier, multi-rank/GPU,
production scaling/I/O budgets and native separate positive-time birthtrace
emission remain unqualified. No push was performed.
