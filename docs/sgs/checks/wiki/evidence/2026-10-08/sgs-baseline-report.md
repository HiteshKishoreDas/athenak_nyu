# P0 / M0 untouched AthenaK baseline reproduction

Executed 2026-10-08, serial CPU on this Mac. **M0 baseline evidence is reproduced.**
No solver, input, harness, repository documentation, or existing local work was edited
by this packet. Subsequent SGS edits by other packets are outside this baseline.

## Source and tools

- AthenaK `/Users/hitesh/hitesh/git/athenak_nyu` at
  `07c529202befbeb6ef5e694d95286929618ff7c1`.
- Clean Kokkos submodule `6739bc623081648af9e752b616d9671527922cbf`,
  CMake-reported version 4.7.2.
- Initial local status: modified `.gitignore` adding `vis/.venv`; untracked
  `AGENTS.md`. Both preserved. Applicable `AGENTS.md`, wiki `agent-runbook.md`,
  and `development-plan.md` were read.
- Apple Clang 17.0.0 (`clang-1700.6.4.2`), `/usr/bin/c++`, C++17,
  arm64 Apple Darwin 25.3.0; CMake 4.1.2. Serial Kokkos, double precision,
  MPI/OpenMP disabled. Build concurrency `-j4`.
- Raw manifest, binary/input SHA-256, exact commands, return codes, timings,
  logs and outputs: `/private/tmp/sgs-baseline-runs`.
- Hardware model/memory/core enumeration via `sysctl` was denied by sandbox.
  No hardware-performance portability claim follows from this packet.

## Build results

All three builds passed. Build time figures below are configure-log completion to
build-log completion, rather than instrumented wall-clock totals.

| Build | Directory | Configuration | Approximate seconds |
|---|---|---|---:|
| Debug | `/private/tmp/athenak-sgs-baseline-debug` | Debug, built-in pgens, Kokkos debug/bounds/DualView checks | 42.8 |
| Release | `/private/tmp/athenak-sgs-baseline-release` | Release `-O3 -DNDEBUG`, built-in pgens | 83.1 |
| Turbulence | `/private/tmp/athenak-sgs-baseline-turb` | Release, `PROBLEM=turb` | 79.4 |

Default builds used the working source while its tracked solver files were untouched.
For turbulence, pinned AthenaK and Kokkos Git archives were extracted into
`/private/tmp/athenak-sgs-baseline-source`; this protects baseline reproduction from
subsequent solver edits. Build logs contain no observed error.

Reproduction configuration (use fresh output directories):

```sh
cmake -S /private/tmp/athenak-sgs-baseline-source \
  -B /private/tmp/athenak-sgs-baseline-debug \
  -DCMAKE_BUILD_TYPE=Debug -DAthena_ENABLE_MPI=OFF \
  -DAthena_ENABLE_OPENMP=OFF -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
cmake --build /private/tmp/athenak-sgs-baseline-debug -j4
```

Release substitutes `CMAKE_BUILD_TYPE=Release` and its separate directory; turbulence
also sets `PROBLEM=turb`. The manifest records actual source/build directories.

## Selected existing tests

These are selected instances, not a claim that the full regression suites passed.
Largest hydro-wave mesh was 32×16×16; turbulence was 16³; shock/scalar were 32×1×1.
All successful simulation invocations returned zero.

| Case | Configuration / oracle | Result |
|---|---|---|
| Left/right sound wave | Existing `linear_wave_hydro.athinput`, RK2/PLM/HLLC, 16×8×8 and 32×16×16, multiblock, amplitude 1e-6, t=1 | PASS: RMS-L1 6.727418e-7 → 1.875944e-7; ratio 0.278851 < 0.29; fine error < 2.5e-7; left/right printed errors identical |
| Entropy wave | Same, wave_flag=3, vflow=1 | PASS: 4.477191e-7 → 1.197284e-7; ratio 0.267419 < 0.30; fine error < 2.5e-7 |
| Release wave | Left sound, fine case | PASS: same printed RMS-L1 as debug (1.875944e-7) |
| Sod | Existing input, 32 cells in 4 blocks, HLLC, t=0.2, debug and release | PASS smoke: finite positive cell fields; density minimum 0.125 and internal-energy density minimum 0.25. No exact-Sod error oracle was executed |
| Two passive scalars | Existing advection pgen/input copied to temp, nscalars=2, advect_dens=false; sine advection, RK2/PPMX, 32 cells in 4 blocks, t=1 | PASS: scalar identities agree, final mean 1.0000000000000002; one-period L1 error 8.3155854556e-5. This is kinematic isothermal legacy scalar coverage |
| Intermediate-checkpoint restart | Restart from continuous scalar run's t≈0.5 dump with untruncated dt | PASS: final tabular scalar state agrees bit-for-bit with continuous run; max absolute difference 0 |
| Stopped-run restart | Stop scalar run at t=0.5 then restart to t=1 | Strict equivalence FAIL: max scalar difference 1.0763146867e-6. See pre-existing issue below |
| Short forced turbulence | Existing `turb.athinput`, separate `PROBLEM=turb`, fixed rseed=20261008, 16³ in eight blocks, t=0.1 | PASS smoke: finite history; printed mass 1 → 1; total E 1.78571 → 1.79572; final kinetic E 0.00811535. History is default precision; no forcing-work conservation or statistical validation claim |

Existing wave thresholds were taken directly from `tst/scripts/hydro/hydro_linwave.py`.
Scalar state was written at 17 decimal digits. Sod/turbulence smoke results have the
original inputs' default output precision. `analyze.py` reruns the recorded numerical
assertions and writes `metrics.json`.

## Failures and limitations retained in evidence

1. **Pre-existing stopped-restart dt behaviour.** The stopped t=0.5 scalar run
   takes cycle 81 with a near-zero final clipped dt. Its restart stores that dt,
   then grows it: at cycle 100, time is still ≈0.5 and dt is 1.629815e-9.
   It finishes at cycle 201 versus cycle 161 for the continuous run. The changed
   timestep schedule produces a 1.0763e-6 scalar difference. A checkpoint written
   before final-time clipping reproduces exactly. This is baseline behaviour,
   neither an SGS regression nor evidence that all restart strategies reproduce
   bitwise. No source fix was attempted in P0.
2. **Initial verification setup failed.** CLI overrides cannot add an absent
   parameter: attempted `hydro/nscalars=2` on the unmodified advection input exits
   with "Parameter 'nscalars' ... not found". That attempt and the wrapper's missing
   restart-file exception remain in logs. Corrected verification supplied a temp
   input with `nscalars` and `output4` defined; no solver change was needed. Treat
   the initial failure as verification setup failure, not a physics regression.
3. **Shared harness deliberately bypassed after audit.** `tst/run_tests.py`
   unconditionally removes `cwd/build` before/after; `utils/athena.py` hardcodes
   `cmake3` and `make -j8`. `cmake3` is absent here. Direct CMake/input runs avoided
   shared cleanup and complied with the -j4 cap. No full harness result is claimed.
4. MPI, GPU, AMR, cooling/dust-specific tests, complete solver/reconstruction
   matrices and stationary turbulence were not run. Baseline serial smoke and
   selected numerical convergence do not qualify those paths.

## Evidence and resource limits

- `manifest.json`: exact revisions/configuration, input and executable hashes.
- `runs.json`, `scalar-runs.json`, `release-runs.json`, `checkpoint-runs.json`,
  `turb-runs.json`: exact commands, working directories, elapsed times, return codes.
- `debug-*.log`, `release-*.log`, `turb-*.log`: configure/build logs.
- Named run directories: raw logs, errors, fields, histories and scalar restarts.
- `metrics.json`, `analyze.py`: metrics and reproducible assertions.
- `run_selected.py`, `run_scalar.py`: direct execution recipes; `run_selected.py`
  intentionally retains the failed first scalar attempt. Use `run_scalar.py` for
  corrected scalar verification; use fresh run directories to avoid output append.
- Raw simulation/evidence output is approximately 0.16 MB; build trees and archived
  source are separate from the requested <200 MB simulation-output cap. No long
  turbulence campaign was launched. P0 finished within its 15-minute cap.

**Handoff:** serial debug/release baseline is healthy for the selected cases. Preserve
stopped-restart timestep-schedule behaviour as a known pre-existing limitation in
future SGS comparisons. Baseline scripts/output stay outside the repository.
