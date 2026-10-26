# Reproduce curve input and encounter-field checks

Run from the shock-lagrangian worktree. Use fresh output directories on reruns.
Current report.json records exact source/configuration and executable hashes.
Same existing Serial/double CMake configuration as the preceding output-mode gate:

```sh
cmake --build build-lagrangian -j4
cmake --build build-lagrangian-prescribed -j4
PY=/Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python
REF=build-lagrangian/evidence/crv-sealed-20261025
$PY tst/tracer-output/shock_output_modes_check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs "$REF" --workdir build-lagrangian/evidence/curve-input-MY-NEW-RUN
python3 tst/tracer-output/shock_binary_reader_check.py \
  --runs build-lagrangian/evidence/curve-input-MY-NEW-RUN \
  --reference "$REF/planar" --output build-lagrangian/evidence/curve-input-MY-NEW-READER
mkdir build-lagrangian/evidence/curve-input-MY-NEW-SNAPSHOT-FIXTURE
ln -s "$PWD/build-lagrangian/evidence/curve-input-MY-NEW-RUN/curve_snapshot" \
  build-lagrangian/evidence/curve-input-MY-NEW-SNAPSHOT-FIXTURE/planar
python3 tst/tracer-output/curve_binary_reader_check.py \
  --runs build-lagrangian/evidence/curve-input-MY-NEW-SNAPSHOT-FIXTURE \
  --historical build-lagrangian/evidence/thp-qualified-CT-20261025/planar/curve_thermo_history \
  --output build-lagrangian/evidence/curve-input-MY-NEW-SNAPSHOT-READER
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
XDG_CACHE_HOME=/private/tmp/athenak-font-cache \
$PY tst/tracer-output/shock_output_modes_plots.py \
  --runs build-lagrangian/evidence/curve-input-MY-NEW-RUN --reference-runs "$REF" \
  --output build-lagrangian/evidence/curve-input-MY-NEW-PLOTS
```

SM22 makes nine real native trials; all rejected inputs must fail with their specific
message. Histories accepted via Mach aliases must match exact event identities and
evidence from the shock-field selection. Curve snapshots must have zero particle
cache builds, no encounter files and one curve stencil batch per saved clock.
Other groups exercise restarts, bounded queues, native scheduling, live AMR and
independent thermo predictions. Output-mode reference fixtures are also retained
in the preceding sealed gallery. The historical THP directory is intentionally
separate: a CRV directory cannot substitute for the old-format reader test.

Expected results:22 native +47 event-reader +28 snapshot-reader groups, all pass.
Plots:10. No MPI/GPU runtime or production performance claim.
