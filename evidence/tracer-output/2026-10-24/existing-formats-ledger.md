# Keep existing output formats; remove the unrequested format experiment

2026-10-24. User correction: "i don't want a new format!! I want to do add minimal
friction to teh existing user experienve". Worktree: `.worktrees/shock-lagrangian`,
branch `codex/shock-lagrangian`, starting clean at `c82e3d80`. Other branches preserved.

## Accepted scope and implementation

- Removed the optional `.tbin` C++ writer, Python codec, format switch, format
  members and representation-specific journal contracts. Deleted the actual
  implementation, rather than retaining a hidden option or another default.
- Restored pre-experiment (`a8e19214^`, also `66e436e1`) particle shock,
  curve shock, geometry, thermo and paired-checkpoint writer source exactly.
  Particle THP v2 is unchanged. No pusher, detector or encounter-state changes.
- Retained ordinary AthenaK geometry scheduling (`file_type=curves`), using the
  same existing ASCII v2 snapshot writer. Existing cycle/initial/final geometry
  controls continue to work without adding output blocks.
- Retained `read_run`, exact identity handling and `plot_curves.py`, simplified
  to the supplied text/CSV/THP readers. Retained the CSV particle trace helper;
  no converter, new dependency or format-selection step is needed.
- Corrected active documentation and examples. Geometry is ASCII v2 `.txt`,
  not VTK as an earlier chat response mistakenly stated. The prior output-format
  evidence and plots remain historical; they do not describe current code.

## Commands actually run

From `/Users/hitesh/hitesh/git/athenak_nyu/.worktrees/shock-lagrangian`:

```sh
cmake --build build-lagrangian -j8
cmake --build build-lagrangian-prescribed -j8
python3 tst/tracer-output/check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs build-lagrangian/evidence/native-output-qualified-final-20261023 \
  --workdir build-lagrangian/evidence/existing-output-final-20261024
python3 tst/curve-thermo/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/existing-output-CT-20261024
python3 tst/shocks-curves/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/existing-output-CSH-20261024
python3 tst/shocks-lagrangian/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/existing-output-LT-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-setup/check.py \
  --executable build-lagrangian-prescribed/src/athena \
  --baseline /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-thermo-synergy-20261022/final-reviewed/runtime/build-lagrangian-prescribed/athena \
  --workdir build-lagrangian/evidence/existing-output-setup-final-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/curve_reader_check.py \
  --runs build-lagrangian/evidence/existing-output-final-20261024 \
  --three-d-run build-lagrangian/evidence/existing-output-CT-20261024/forced_both \
  --output build-lagrangian/evidence/existing-output-reader-final-20261024
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/visuals.py \
  --run build-lagrangian/evidence/existing-output-final-20261024 \
  --output build-lagrangian/evidence/existing-output-plots-final-20261024
./build-lagrangian/shock_detector_test --self-check
./build-lagrangian/shock_encounter_test
python3 tst/shocks-curves/io_reader.py
python3 -m py_compile vis/python/tracer_history.py vis/python/shock_history.py \
  vis/python/curve_shock_history.py vis/python/tracer_thermo_history.py \
  vis/python/plot_curves.py tst/tracer-output/check.py \
  tst/tracer-output/curve_reader_check.py tst/tracer-output/visuals.py
git diff --check
```

Both builds and all final commands above passed. The builds reused the established
local Serial/double configurations and Kokkos submodule. Compiled fingerprints,
executable/build-cache hashes and exact native commands are in each raw report.
Every LT native run's executable digest was checked against the completed builds;
none used the earlier executable during rebuild.

## Results and gates

| Suite | Passing groups |
| --- | ---: |
| Output scheduling, histories, AMR and restart | 15 |
| Shared thermo (CT01–CT18) | 18 |
| Curve shocks (core CSH) | 20 |
| Continuous particles (LT) | 14 |
| Startup summaries (SUM01–SUM09) | 9 |
| Curve reader/plotter (CR01–CR11) | 11 |
| Total | 87 |

Additional detector self-check, encounter state checks and standalone curve IO
checks passed. `ctest --test-dir build-lagrangian --output-on-failure` reported
no registered tests; it is not counted as a test pass.

TO01/TO02 compare every particle/curve thermo, event and terminal-view record,
clock and geometry against retained text/CSV controls from the pre-removal build;
particle THP bytes and planar conserved-fluid tables match exactly. Scheduled
geometry on/off is also passive. TO03's independent manufactured thermo error is
zero, below the unchanged budget of $2\times10^{-12}$; actual mesh leaf count spans 4–16.

TO21 refines/coarsens with 17 inside-born and 22 outside-born nodes. Histories begin
at actual birth, entry is incomplete when born inside, and no parent encounter is
copied. Checkpoint cycle 13 resumes to exactly matching particle/curve samples,
event identities, suspended views and final geometry. Zero-step replay also leaves
journal bytes unchanged. This run has zero separate birthtrace files; the existing
per-step traces contain the birth observations. Separate positive-time birthtrace
emission is not newly qualified.

Removed TO06/08/09/11/13/14 tests belong to the deleted codec (wire/endian/checksum,
format selection, format-specific rank/schema or binary size comparisons). Their
historical results remain retained. No physical tolerance or remaining numerical
acceptance criterion was weakened to pass the format removal.

Nine fresh plots use actual newly generated text/CSV/THP records: seven curve
reader views, shared particle/curve thermo/strength and live-AMR restart equality.
Descriptions sit beside plots in the handoff gallery. No exact physical crossing
time or dense reconstructed trajectory is inferred.

## Retained failures and corrections

- First output run (`existing-output-20261024`): 14/15 groups passed; TO15 failed
  because restoring the older reader also removed its newer CSV `read_trace`
  helper. Retained that helper without binary handling, then reran the full
  output suite in a new directory: 15/15 pass.
- First startup invocation used system Python: all nine native checks passed,
  then plotting failed because Matplotlib is absent there. Reran using the
  existing scientific environment in a fresh directory: complete exit 0.
- First diagnostic-plot invocation wrote its first plot then failed due to a
  reversed `(label,key)` unpacking in the new plot code. Corrected it and
  compared encounters by exact identity. Reran in a fresh destination: exit 0.
- A direct detector invocation without its required argument returned usage
  (exit 2). Its intended `--self-check` invocation then passed.

Raw runs, failed fixtures and logs are retained. Final summary and report hashes:
`existing-formats-report.json`; command logs: `existing-formats-logs/`.

## Remaining limits and next actions

Only local Serial/double is qualified here. MPI/GPU, RUSTY/Frontier, production
runtime/memory/output budgets and human scientific approval remain separate,
pending gates. Tracked restarts retain strict source-fingerprint guards; matching
layouts alone do not permit cross-code checkpoint continuation. The reader loads
histories into RAM. The removed binary experiment remains recoverable in Git and
sealed historical artifacts only. No push performed.
