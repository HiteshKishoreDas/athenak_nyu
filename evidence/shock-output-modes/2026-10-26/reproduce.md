# Reproduce the 2026-10-26 output modes

Use a fresh destination for every command. Existing data is retained deliberately.
The exact executed native inputs/commands, executable/cache/source fingerprints
and expected rejection outcomes are in the native report linked by report.json.
Current checkout is codex/shock-lagrangian in .worktrees/shock-lagrangian.

Build Serial/double with the existing Kokkos commit08ceff92bcf3a828844480bc1e6137eb74028517:

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

Use the supplied scientific environment. REF below is the previous accepted CRV
run directory in this workspace. The gallery's compact `reference-runs/` contains
identical inputs, snapshot streams and encounter records needed for these commands
and can be substituted for REF when reproducing from a fresh source checkout.
No prior executable or scientific initial-condition choices need to be invented.

```sh
PY=/Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python
REF=build-lagrangian/evidence/crv-sealed-20261025
$PY tst/tracer-output/shock_output_modes_check.py \
  --executable build-lagrangian/src/athena \
  --prescribed-executable build-lagrangian-prescribed/src/athena \
  --reference-runs "$REF" --workdir build-lagrangian/evidence/modes-MY-NEW-RUN
python3 tst/tracer-output/shock_binary_reader_check.py \
  --runs build-lagrangian/evidence/modes-MY-NEW-RUN \
  --reference "$REF/planar" --output build-lagrangian/evidence/modes-MY-NEW-READER
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl \
XDG_CACHE_HOME=/private/tmp/athenak-font-cache \
$PY tst/tracer-output/shock_output_modes_plots.py \
  --runs build-lagrangian/evidence/modes-MY-NEW-RUN --reference-runs "$REF" \
  --output build-lagrangian/evidence/modes-MY-NEW-PLOTS
build-lagrangian/shock_detector_test --self-check
build-lagrangian/shock_encounter_test
```

The main native gate runs the new HistoryStore protocol test and its expected
failure variants, including queued destroyed IDs, domain exits, signed-1 and a
sequence above2^53. These do not qualify a live outflow configuration. The snapshot
reader's28 independent corruption groups also pass against the final sampled files.

MPI compilation uses the prescribed configuration above with Athena_ENABLE_MPI=ON.
It compiled, but no runtime MPI retry was performed: the previously recorded
PMIx/PRRTE ABI mismatch remains. GPU/RUSTY/Frontier runs, production budgets and
human scientific sign-off are pending. The short32³ forced-box fixture is an
engineering coexistence check, not converged turbulence statistics.
