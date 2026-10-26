# Reproduce unified CRV checks locally

From the codex/shock-lagrangian worktree. Reuse the qualified Kokkos checkout and
existing Serial/double build configurations; MPI build is compile-only.

```sh
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

Use the existing scientific environment at the primary checkout:

```sh
PY=/Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python
$PY tst/tracer-output/curve_binary_check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs build-lagrangian/evidence/thp-final-SHOUT-20261025 \
  --workdir build-lagrangian/evidence/crv-MY-NEW-RUN
$PY tst/tracer-output/curve_binary_reader_check.py \
  --runs build-lagrangian/evidence/crv-MY-NEW-RUN \
  --historical build-lagrangian/evidence/thp-final-SHOUT-20261025/planar/curve_thermo_history \
  --output build-lagrangian/evidence/crv-MY-NEW-READER
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
$PY tst/tracer-output/curve_binary_plots.py \
  --runs build-lagrangian/evidence/crv-MY-NEW-RUN \
  --reference-runs build-lagrangian/evidence/thp-final-SHOUT-20261025 \
  --output build-lagrangian/evidence/crv-MY-NEW-PLOTS
```

Directories must be new. Exact qualified inputs, commands, source fingerprints,
executable digests and outcomes live in the source-bound reports referenced by
report.json. The reference-runs argument uses retained pre-change native outputs
for exact parity. Without that argument, the native runner uses repository fixtures
and still checks structure/oracles/restart, but does not claim old-output parity.
The supplied old SHOUT inputs/outputs are required for the exact qualified parity
configuration and storage plots; their digests are in the evidence bundle.

Existing physics/history gates:

```sh
python3 tst/shocks-lagrangian/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/crv-MY-NEW-LT
python3 tst/shocks-curves/acceptance.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --workdir build-lagrangian/evidence/crv-MY-NEW-CSH
build-lagrangian/shock_detector_test --self-check
build-lagrangian/shock_encounter_test
```

For a saved run, the existing plot CLI works directly with CRV:

```sh
$PY vis/python/plot_curves.py /absolute/path/to/run \
  --field pressure --output /absolute/path/to/new-plot.png
```

No native hardware beyond local Serial/double is qualified. MPI compiles; runtime
MPI remains blocked by the previously recorded PMIx/PRRTE ABI mismatch. Production
budgets, GPU/RUSTY/Frontier runs and human scientific sign-off remain pending.
