# Independent clean committed-source reproduction

Executed 2026-10-08 on this Mac. All bounded checks below passed. No source,
wiki, user input, or tracked checkout files were edited by this packet.

## Exact source and build

AthenaK commit `14620bb33d580246e32d80f4bd584c6c021bbd1e` was exported with
`git archive` to `/private/tmp/athenak-sgs-clean-source`. The archive contains
only this committed tree, independent of the shared checkout's later edits.
Kokkos was copied from the clean pinned baseline archive; the commit's gitlink
is `6739bc623081648af9e752b616d9671527922cbf` (CMake version 4.7.2).

Fresh build `/private/tmp/athenak-sgs-clean-build`: CMake 4.1.2, Apple Clang
17.0.0, arm64 Darwin 25.3.0, Release C++17, CPU Kokkos Serial, double precision,
MPI/OpenMP disabled, built-in problem generators, BUILD_TESTING enabled.
Configure and build returned zero; no observed compiler warning/error.
Build concurrency was three; configure-log completion to build-log completion
was about 72 seconds, within the five-minute cap. All simulation execution
elapsed times sum to 1.20 seconds, within the 60-second test cap.

Exact configure/build/CTest commands, revisions and executable/archive hashes:
`/private/tmp/sgs-clean-reproduction-manifest.json`. Configure and build logs:
`/private/tmp/sgs-clean-configure.log`, `/private/tmp/sgs-clean-build.log`.

## Reproduced results and fixed acceptance criteria

| Check | Criterion | Observed |
|---|---|---|
| Independent closure device oracle | Frozen scaled tolerance 5e-13 | PASS: 1,025 cases, 22,575 assertions, max error 1.4521283901912292e-15 |
| Existing left sound wave, SGS absent | 16/32 fine RMS-L1 <2.5e-7 and ratio <0.29 | PASS: 6.727418e-7 to 1.875944e-7, ratio 0.27885051887663287 |
| Existing right sound wave, SGS absent | Same criteria | PASS: identical printed errors to left wave and untouched baseline |
| Existing Sod, SGS absent | All output fields finite; rho and internal energy positive | PASS: minimum rho 0.125, minimum internal-energy density 0.25 |
| Two legacy passive scalars, SGS absent | Identical scalar species, mean one within 1e-13 | PASS: both species identical, final mean 1.0000000000000002 |
| Scalar intermediate-checkpoint restart | Final tabular state agrees bit for bit with continuous run | PASS: maximum absolute difference 0 |
| Generic uniform SGS seed through built-in pgen | Native restart layout six variables, all finite/positive; initial/final field difference <1e-12 and specific K seed error <1e-12 | PASS: 1,024 active cells, six variables, scaled Linf difference 0, seed error 0, final time .01, K .03 |

The SGS smoke used the built-in linear-wave pgen with amplitude zero and generic
uniform initialization, k=.03, periodic cubic active cells, RK2/PLM/LLF and
CFL=.1. Dissipation was disabled (ceps=0) to make exact invariant preservation
the oracle. It exercises the committed generic initialization and storage without
requiring the custom verification pgen. Final density was 1 and thermal-energy
density 0.8999999999936998.

Wave inputs and Sod came directly from the exported committed source. The
corrected passive-scalar baseline input was copied unchanged into the temporary
run directory; this preserves the baseline's isothermal kinematic RK2/PPMX
advection coverage. No thresholds were changed after results.

## Evidence and limitations

- `/private/tmp/sgs-clean-ctest.log` and build `Testing/Temporary/LastTest.log`:
  CTest return status and full closure oracle output.
- `/private/tmp/sgs-clean-runs/runs.json`: exact simulation commands, directories,
  elapsed times and return codes (all eight returned zero).
- `/private/tmp/sgs-clean-runs/metrics.json`: numerical metrics above.
- `/private/tmp/sgs-clean-runs/input-sha256.json`: input fingerprints.
- `/private/tmp/sgs-clean-tests.log`: adapter assertions and metrics.
- `/private/tmp/sgs-clean-run.py`: self-contained temporary reproduction adapter;
  native restart reader imported from the exact exported committed tree.
- Named run directories preserve fields, restart files and logs.

This is a clean Release build, closure oracle, selected SGS-disabled regression
and generic seeded-SGS smoke. It is not the full regression harness, custom SGS
pgen campaign, turbulent reference calibration, GPU/MPI/AMR qualification, or a
claim that every restart strategy is bitwise reproducible. The baseline's stopped
run/final-time-clipping restart limitation was not changed or reclassified.
