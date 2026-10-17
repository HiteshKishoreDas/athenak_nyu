# feature_dump shock merge evidence — 2026-10-17

Target before merge: `cea1dafcf6830bfecba3e2aea5622b5c943a761c`.
Source, preserved unchanged: `2450e1c21b0a97c0f1635da75a20701c7e0852a4`.
Source branch: `codex/shock-multiphase`; all its committed changes included.
Merge prepared on `codex/merge-shocks-feature-dump`, then target updated only
by a fast-forward after validation. No source deletion, force update or push.

The only added production behavior beyond preserving both parents is explicit
rejection of enabled shock tracking without particles. Original silent failure,
first helper failure and corrected 11-process interaction checks are retained.
Existing material-curve exclusion of particles/forcing remains; no curve shock
history is claimed. New geometry nodes must never receive invented past history.

`verification.json` binds final source files, binary hashes, unchanged 21 held-out
measurement dictionaries, test scope and durable raw archive. `independent-reviews.json`
records timing/integration/persistence review. Historical imported evidence retains
its original exact-source identity. New merged fingerprint:
`48c53dc5c167604269d121412d4316215f09be8224e053d98af2eae9908a2b7b`.

## Reproduce local checks

From `/Users/hitesh/hitesh/git/athenak_nyu` on the merged branch, with initialized
Kokkos. Commands use fresh build/results directories; never overwrite prior evidence.
The temporary integration build reused the existing local Kokkos commit through a
symlink, removed after compilation. Future builds should use the normal submodule.

```sh
ROOT=/Users/hitesh/hitesh/git/athenak_nyu
PY=/Users/hitesh/hitesh/git/shockit/.venv/bin/python
cmake -S "$ROOT" -B /tmp/merge-shock-build-new -DCMAKE_BUILD_TYPE=Release -DAthena_SHOCK_TESTS=ON
cmake --build /tmp/merge-shock-build-new -j8
/tmp/merge-shock-build-new/shock_encounter_test
"$PY" "$ROOT/tst/shocks/parity.py" /tmp/merge-shock-build-new/shock_detector_test --reference "$ROOT/evidence/shocks/reference" --report /tmp/merge-parity-new/report.json
"$PY" "$ROOT/tst/shocks-mc/feature_dump.py" --executable /tmp/merge-shock-build-new/src/athena --workdir /tmp/merge-interactions-new
"$PY" "$ROOT/tst/shocks-mc/integration.py" --executable /tmp/merge-shock-build-new/src/athena --workdir /tmp/merge-integration-new --report /tmp/merge-integration-new-report/report.json
"$PY" "$ROOT/tst/shocks-mc/replay.py" --executable /tmp/merge-shock-build-new/src/athena --workdir /tmp/merge-replay-new --report /tmp/merge-replay-new-report/report.json
"$PY" "$ROOT/tst/shocks-mc/qualification.py" --executable /tmp/merge-shock-build-new/src/athena --workdir /tmp/merge-physical-new --report /tmp/merge-physical-new-report/report.json
MPLCONFIGDIR=/tmp/merge-mpl-new "$PY" "$ROOT/tst/shocks-mc/forcing_checks.py" --executable /tmp/merge-shock-build-new/src/athena --workdir /tmp/merge-forcing-new --output /tmp/merge-forcing-report-new
cmake -S "$ROOT" -B /tmp/merge-curve-build-new -DCMAKE_BUILD_TYPE=Release -DPROBLEM=material_curve_test -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build /tmp/merge-curve-build-new -j8
"$PY" -B "$ROOT/tst/scripts/material_curves_io/test_reader.py"
"$PY" -B "$ROOT/tst/scripts/material_curves_io/test_cpp.py" --build /tmp/merge-curve-build-new --kokkos-library /tmp/merge-curve-build-new/kokkos/core/src/libkokkoscore.a
"$PY" -B "$ROOT/tst/scripts/material_curves/production.py" --executable /tmp/merge-curve-build-new/src/athena --output-dir /tmp/merge-curve-core-new --suite core
"$PY" -B "$ROOT/tst/scripts/material_curves/production.py" --executable /tmp/merge-curve-build-new/src/athena --output-dir /tmp/merge-curve-regression-new --suite regression
"$PY" -B "$ROOT/tst/scripts/material_curves/production.py" --executable /tmp/merge-curve-build-new/src/athena --output-dir /tmp/merge-curve-unsupported-new --suite unsupported
```

NumPy is required only by existing scientific helpers/new float64 output checks.
No new dependencies. Per-case reports preserve actual executed command arrays,
exact inputs/configuration, elapsed times and failures; build logs preserve configuration.
The core/held-out/forcing tolerances were not changed. Native forcing's independent
thermal oracle remains failed; engineering compatibility passes are separate.
GPU/MPI/AMR runtime requalification, curve-shock support and human scientific sign-off
remain pending. Old shock checkpoints deliberately require their matching old binary.

Raw directory, checkpoint pairs and both runtime binaries/source archive:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/feature-dump-merge-20261017`.
`raw-sha256.json` binds every original run file; `sha256.json` binds compact evidence
and excludes itself. Do not treat old imported shock evidence as current merged-code evidence.
