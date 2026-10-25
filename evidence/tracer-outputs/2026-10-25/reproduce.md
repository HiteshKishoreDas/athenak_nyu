# Reproduce native output integration on the laptop

Working checkout:
`/Users/hitesh/hitesh/git/athenak_nyu/.worktrees/shock-lagrangian`.
Use fresh output destinations: the scripts deliberately preserve existing results.
Exact executed commands, input files and binary/cache fingerprints are also in
all report.json run records. The sealed runtime includes the exact three build
caches, configuration and generated source-hash headers.

## Build

The following fresh-build recipe reproduces the recorded configuration choices;
the already configured builds were rebuilt with cmake --build during this task.
Recorded compiler: /usr/bin/c++; empty CMAKE_BUILD_TYPE; Serial/double. Kokkos
commit: 08ceff92bcf3a828844480bc1e6137eb74028517.

```sh
git submodule update --init kokkos
cmake -S . -B build-lagrangian -DPROBLEM=built_in_pgens \
  -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF \
  -DAthena_ENABLE_OPENMP=OFF -DAthena_ENABLE_MPI=OFF \
  -DAthena_SINGLE_PRECISION=OFF -DAthena_SHOCK_TESTS=ON
cmake -S . -B build-lagrangian-prescribed -DPROBLEM=material_curve_test \
  -DKokkos_ENABLE_SERIAL=ON -DKokkos_ENABLE_OPENMP=OFF \
  -DAthena_ENABLE_OPENMP=OFF -DAthena_ENABLE_MPI=OFF \
  -DAthena_SINGLE_PRECISION=OFF -DAthena_SHOCK_TESTS=ON
cmake --build build-lagrangian -j8
cmake --build build-lagrangian-prescribed -j8
```

MPI compile used the same prescribed problem/Serial/double configuration with
Athena_ENABLE_MPI=ON. It compiled; no successful MPI runtime is claimed. Do not
infer GPU or cluster support from the local Serial checks.

## Quick reproducible native runs and plots

This produces the accepted translating-shock/live-AMR reference inputs and the
new THP/BIN tests. The latter includes the existing short 32³ forced-box fixture.
All initial physical values are taken from the accepted repository fixtures;
new physical initial conditions are not introduced.

```sh
thp_check_root=build-lagrangian/evidence/my-output-check
python3 tst/tracer-output/check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir "$thp_check_root/reference"
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/shock_output_check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs "$thp_check_root/reference" \
  --workdir "$thp_check_root/shocks"
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
  /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/tracer-output/shock_output_plots.py \
  --run "$thp_check_root/shocks" --thermo-runs "$thp_check_root/reference" \
  --output "$thp_check_root/plots"
python3 tst/tracer-output/thp_reader_check.py
```

Optional matched-size plot: add --csv-reference pointing to retained pre-THP,
signed-ID check.py runs. Current evidence used
build-lagrangian/evidence/signed-id-reviewed-output-20261025. That plot first
checks exact equality of all thermo records. File size measurements are not a
throughput or production-memory benchmark. The final source-bound suite directories
are thp-final-{TO,CT,CSH,LT,VB,setup,reader,SHOUT,boundary0,boundary8}-20261025.

## Full final gates

The same executable pair is accepted by these scripts; use fresh --workdir paths:

- tst/curve-thermo/acceptance.py: 18 groups.
- tst/shocks-curves/acceptance.py: 20 groups, unchanged physical timing/oracle budgets.
- tst/shocks-lagrangian/acceptance.py: 14 groups.
- tst/tracer-output/variable_id_check.py: 9 groups; --reference-runs is the fresh check.py run.
- tst/tracer-output/id_check.py: prescribed executable, --curve-bits 0 or8: 8/6 groups.
- tst/tracer-setup/check.py: 9 groups; retained old executable baseline is recorded in its report.
- tst/tracer-output/curve_reader_check.py: 11 groups; --runs is check.py output,
  --three-d-run is the curve-thermo forced_both run, --output is a fresh folder.
- tst/tracer-output/shock_output_check.py: 9 groups, including zero-step restart,
  float32 exact cast, empty populations and accepted 3D forced-box coexistence.

Additional checks:

```sh
python3 tst/scripts/material_curves_io/test_reader.py
python3 tst/shocks-curves/io_reader.py
python3 tst/scripts/material_curves_io/test_cpp.py \
  --build build-lagrangian --kokkos-library build-lagrangian/kokkos/core/src/libkokkoscore.a
./build-lagrangian/shock_detector_test --self-check
./build-lagrangian/shock_encounter_test
git diff --check
```

The isolated C++ IO harness stubs physics/audits; the native suites exercise
production integration. Complete command arrays and input/restart hashes remain
in the final run reports. No tolerance was weakened. Scientific approval and
MPI/GPU/cluster production qualification remain separate pending gates.
