# Native compression-corrected SGS heating budget

New diagnostic and verification work on `codex/multiphase-schmidt-federrath-hd`
in `/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu`.
Local source commits: `1044ddd6` (native ledger/tests), `e97db55f` (exact restart
origins) and `7bf794ad` (omit undefined fraction bars). No push or merge.

The [paper's Section 4.2](https://arxiv.org/html/1010.4492#S4.SS2) reports about
75% SGS and 25% numerical heating after subtraction of resolved compression.
This is a reference, not a pass threshold. The paper averages 2T–10T; these
one-seed solenoidal pilots cover 2T–4T. Native forcing normalization/equivalence
remains unresolved, and 32³ is below the paper's minimum 64³ grid.

## What the ledger measures

For periodic, ideal-gas HD without cooling, all quantities are volume integrals:

```
Delta U = internal thermal energy increment
C = signed applied resolved thermal pressure work
D = applied SGS dissipation counter increment
N = Delta U - C - D
SGS fraction = D / (Delta U - C)
numerical fraction = N / (Delta U - C)
```

C uses the actual reconstructed Euler LLF pressure in each native momentum
flux, Pface=(Pleft+Pright)/2, including fallback. Its global rate is the negative
pressure-force work on resolved kinetic energy:
`-sum_faces(area*Pface*(vright-vleft))`. Final-solution RK2 weights are applied
through the existing stage history mechanism. Signed residuals are retained;
fractions are undefined for nonpositive corrected heating. Ratios use integrated
increments, not arithmetic means of interval fractions. Endpoints are clipped by
linear interpolation of native histories.

N includes spatial and RK effects and shock-associated heating. It cannot be
identified with LLF alone or automatically interpreted as unphysical heating.
Reconstruction changes the coupled flow and SGS production as well as numerical
heating. This diagnostic does not isolate the entire cause of low SGS K.

## Minimal source changes

Existing SGS face storage grows by one thermal-pressure component; stage budget
storage grows by one rate. A fifth integrated ledger joins existing counters,
MPI reduction, native history and precise parameter-header restart persistence.
History/reduction capacity grows 32 to 33 to preserve the already-supported
forced twelve-scalar case. No SGS coefficient, stress, energy exchange, fluid
state update, RK weight, diffusivity or support guard changes.

The test harness records exact restart inventory/counter origins in
`history-origin.json`. The report reads the actual native checkpoint header,
not the unused input generated for a restart. Older checkpoints initialize the
previously absent pressure counter at zero at their exact restart time; only
subsequent differences are supported. The 64³ pressure ledger begins around
1.9T, before the requested 2T window. Earlier pressure-free histories are not
concatenated into a fabricated cumulative ledger.

## Executed verification

- 23/23 independent heating checks across 12 native result records pass.
  The initial compression/expansion rate has maximum relative error
  1.61099e-6; the independent DC shear numerical residual has error 1.64564e-7,
  below the rate threshold 2e-5 frozen before execution. Exactly one RK2 step
  dt=1e-7, 16³/32³ compression and 16²/32² shear, no limiting/fallback.
- Fixed-space coupled smooth pressure-work refinement at CFL .1/.05/.025
  gives order 1.87323, within the prespecified 1.5–2.5 range. This is a numerical
  reference, not an exact PDE solution or a bound on turbulent RK heating.
- Native checkpoint restart including the new counter passes. Exact native
  fields match the earlier binary in the bounded 3D control (maximum difference 0).
- Configuration 34/34, PPM qualification 16/16, forced twelve-scalar history
  capacity and native production amplitude 24/24 checks across eight cases pass.
- Independent C++ closure oracle: 1025 cases, 22575 assertions; maximum scaled
  error 1.45213e-15 below 5e-13. Reporting arithmetic and unequal-density cell
  diagnostics self-tests pass; added C++ style/whitespace and new Python E/F
  lint pass. The existing forcing runner passes F lint; its inherited E501
  diagnostics prevent a whole-file E501 claim.
- Dedicated native forcing-origin/restart checks pass 3/3 records; saved origins
  agree exactly with their original checkpoints.

No overall M0–M9 scientific-release approval is implied by these diagnostic checks.

## Completed matched turbulent measurements

Both 64³ bounded continuations reach 4T and pass native admissibility,
forcing, scalar and inventory checks. Final native fields, forcing metadata
and RNG match the preceding scheme runs EXACTLY (maximum field error 0).
The new diagnostic has not changed these turbulent trajectories.

All fractions below use 2T–4T, one solenoidal seed 42, RK2/LLF, CFL .1,
gamma 1.001, acceleration RMS 3 and unchanged closure coefficients.

| Grid/reconstruction | SGS heating | Numerical residual |
|---|---:|---:|
| 32³ PLM CFL 0.1 | 47.483% | 52.517% |
| 32³ PPM4 CFL 0.1 | 63.345% | 36.655% |
| 64³ PLM CFL 0.1 | 56.639% | 43.361% |
| 64³ PPM4 CFL 0.1 | 70.054% | 29.946% |

Integrated native budgets (code volume-energy units):

| Case | Delta U | C | D | N | P | Delta K |
|---|---:|---:|---:|---:|---:|---:|
| 32³ PLM CFL 0.1 | 1.39914229 | 0.0721693167 | 0.630091613 | 0.696881365 | 0.635385282 | 0.00529366995 |
| 32³ PPM4 CFL 0.1 | 1.4087048 | 0.0694591606 | 0.848351397 | 0.490894247 | 0.848542501 | 0.000191103884 |
| 64³ PLM CFL 0.1 | 1.39662361 | 0.0937762416 | 0.73792479 | 0.564922577 | 0.735943497 | -0.00198129218 |
| 64³ PPM4 CFL 0.1 | 1.39150644 | 0.0800361638 | 0.918733969 | 0.392736311 | 0.916361243 | -0.00237272582 |

Worst window budget residuals:

- E+K versus applied forcing: 2.2915e-13.
- Delta K versus P-D: 6.32827e-15.

A substantial competing numerical heating channel is measured. PPM4
sends a larger fraction through SGS than PLM in these coupled runs.
Low instantaneous K does not imply negligible SGS transfer: production
and dissipation can nearly balance. These results do not establish the
paper's full stationary statistics or uniquely identify a closure defect.

[64³ six-panel comparison](report64/heating-comparison.pdf),
[32³ six-panel comparison](report32/heating-comparison.pdf),
[short CFL comparison](report-cfl-checked/heating-comparison.pdf),
[native analytic/refinement controls](report-controls/native-controls.pdf).
The paper-style [PLM report](turbulence32-plm/report.pdf) and
[PPM4 report](turbulence32-ppm4/report.pdf) add 18 panels each.
[Independent review](review.md), [verification counts](verification-summary.json)
and [final provenance](manifest-final.json) distinguish tests, measurements
and limitations. Original failed probes remain visible.

## Short timestep probe and retained failures

At 32³ over 2T–2.1T, halving CFL .1 to .05 changes the SGS fraction by
+.0005419 percentage points for PLM and -.0021366 for PPM4. This is small observed
short-window sensitivity with coupled evolution and interpolated endpoints. It
is not a full 2T–4T RK contribution bound or temporal convergence at 64³.

Both original longer CFL probes hit their 120-second cap and remain FAIL in their
original records. A separately pinned shorter protocol passed; no original
threshold/status was relaxed. The first shorter-probe plot failed because the
requested start preceded its first stored history row. Adding the exact restart
origin corrected coverage. The initial Mac OpenMP configure failed discovery of
existing libomp; explicit native CMake flags corrected it without installing a
new dependency. Original run logs and failure records are retained.

## Reproduction and provenance

Raw data, checkpoints, protocols, binaries and all logs remain at
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-heating-budget`.
The pinned 64³ origin checkpoints are in the preceding scheme-recovered folders;
protocols record their SHA256. Case-specific 64³/100000-cycle settings are
explicit in the continuation protocols; generic suite-manifest control defaults
(32 cells/1000 cycles) are not these turbulent case limits. Compact inputs, histories, logs, results, protocols,
plots and source snapshots are mirrored here; raw checkpoints are outside wiki Git.

`provenance/final-source-provenance.json` distinguishes final local commits from
run-time dirty-source manifests. The pre/post-format native binaries are
byte-identical. Running jobs loaded an earlier reporting harness; their origin
records were reconstructed from the exact pinned checkpoint and independently
reviewed. Final reporting snapshots are saved separately from initial snapshots.

Mac Terminal, from this evidence directory, CPU double Release/OpenMP,
existing Homebrew libomp:

```sh
SGS_SOURCE=/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu
cmake -S "$SGS_SOURCE" -B /private/tmp/athenak-sgs-heating \
  -DPROBLEM=sgs_verify -DCMAKE_BUILD_TYPE=Release \
  -DAthena_ENABLE_MPI=OFF -DKokkos_ENABLE_OPENMP=ON \
  '-DOpenMP_CXX_FLAGS=-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include' \
  -DOpenMP_CXX_LIB_NAMES=omp \
  -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib
cmake --build /private/tmp/athenak-sgs-heating -j4
sh reproduce-heating.sh
```

The reproduction script chooses fresh directories and runs the bounded controls,
fresh matched 32³ cases and plotting. It uses the existing plotting environment
`/private/tmp/sgs-plot-venv/bin/python` (NumPy/Matplotlib) and retained baseline
binary. Source README supplies the portable commands. The recorded 64³ continuations
use 4 threads per lane, 1800-second native caps and unchanged checkpoint physics;
protocols and `continue-heating64.py` preserve exact executed commands. Replay
in fresh output directories. Multi-rank MPI/GPU and wider physics remain separately
unqualified.
