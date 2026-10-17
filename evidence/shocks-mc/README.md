# Exact-code MC engineering evidence bundle

Current source `3de18d97`; base `526a187f938e38634e89a88192ebe9731c4d36bb`;
Kokkos `08ceff92bcf3a828844480bc1e6137eb74028517` (4.4.0).
Current core fingerprint
`3feb3bcdf77d21296516544e67b5ad75f6ac6ca7bc17c4b48be886ddc0fb7cbc`.
`review3/` contains three independent reviews and integrated rechecks:53parity,
20integration,12reject,4replay,18boundary cases,13production exits, actual old-layout
rejection, native ABI/payload checks and21physical cases exactly matching prior results.
All23748heldout CSVs match after source fingerprint substitution.40net source lines
removed; state1000to816bytes; public124-column events unchanged. `ponytail/` is the
previous constructor-extraction packet. Other packets below retain earlier identities.
Scientific sign-off pending. Local commits are not remote pushes/releases.

- `baseline`: frozen criteria, existing MC statistical/restart results and failures.
- `environment`: exact tool/build/parity logs, four CMake caches/generated configs
  and selected baseline identity.
- `clean-reproduction`: clean archived full Release build, native state suite,
  53-case parity and 20-process serial integration; commands/inputs/logs retained.
- `review`: independent timing/persistence audit and26source hashes; style-only
  source equivalence recorded. Source review is separate from native runtime tests.
- `integration/report_stylefinal.json`:20successful finalRelease processes.
- `integration/rejects_stylefinal.json`:12explicit rejection checks.
- `integration/replay_stylefinal.json`:4old-checkpoint/new-basename replay checks.
- `migration/boundary_results.json`:18actual-method cases on final source,
  Serial and MPI singleton; stand-in constructors/count/exit callback are explicit.
- `system/report_style_final.json`: exact source/executable/config provenance,
  baseline/final6-case existing hydro regression, bounded local timings, singleton
  and failedactual2rank launch. `exit_report.json` exercises13actual censored exits.
- `system/sanitizer.json`: ASan/UBSan encounter unit only; full solver not sanitized.
- `system/style.json`:0introduced findings; inherited20warnings remain.
- `qualification/report_stylefinal.json`:21measured heldouts; no scientific pass.
- `qualification/packet_stylefinal`: selected sampled/transfer CSVs, event/suspended
  records, nonempty chunks, allchunk digests and4PNG/PDFfrontplots.
- `qualification/parent_front_bins.json`: independent observable cell-parent bins;
  stochastic transfer/return ambiguity retained, no interpolated exact time.

Historical reports/packets retain their own source fingerprints and earlier
failures. Their outcomes are not silently promoted to current-code acceptance.
Large complete run outputs remain in the absolute directories named in manifests;
compact selected evidence is committed. No simulation binaries/full rawgrids added.
`final-manifest.json` hashes the committed packet and source/configuration inputs.

From checkout, with NumPy installed for parity/existing wave helpers:

```sh
cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug -DAthena_SHOCK_TESTS=ON
cmake --build build-debug -j8
build-debug/shock_encounter_test
python3 tst/shocks/parity.py build-debug/shock_detector_test \
  --reference evidence/shocks/reference --report /tmp/mc-parity-new.json
python3 tst/shocks-mc/integration.py --executable build-debug/src/athena \
  --workdir /tmp/mc-integration-new --report /tmp/mc-integration-new.json
python3 tst/shocks-mc/qualification.py --executable build-debug/src/athena \
  --workdir /tmp/mc-qualification-new --report /tmp/mc-qualification-new.json
python3 tst/shocks/front_bins_mc.py --qualification /tmp/mc-qualification-new.json \
  --report /tmp/mc-parent-front-new.json
```

Use fresh directories. Detailed rejection/replay helper commands and latest packet
provenance are in integration/README.md and qualification/README.md. Documentation,
supported subset, restart recipe and unresolved gates are in docs/shocks-mc/.
