# SGS hydro qualification evidence

This snapshot preserves executed evidence from the 2026-10-08 implementation.
The [implementation report](../../implementation-progress.md) and
[coverage ledger](../../coverage-ledger.md) define claim boundaries. M0–M9 is
not complete. No raw failure is converted to a pass by this bundle.

`bundle-manifest.json` lists exact SHA-256/byte counts. Archives contain original
commands, effective inputs, hashes, metrics, histories, source/build manifests and
logs. Original temporary paths in those records describe where a run actually
happened; extract archives into a fresh directory and use the reproduction
instructions rather than writing over those directories.

| Archive | Contents |
|---|---|
| baseline | Untouched builds and selected existing hydro/scalar/turbulence checks |
| oracle | Independent 1025-case algebra metrics and three detected mutations |
| debug-qualification | Initial failures, revised bounds checks and later fallback/limiter cases; native-double checkpoints retained |
| release-qualification | Earlier31-case Release matrix and native checkpoints |
| final-release-unforced | Final46-record unforced acceptance matrix |
| configuration | 22 configuration/disabled/restart checks |
| forcing-history | Earlier forced restart/reader/diagnostic failures, fixed reruns and no-advance payload audit |
| forcing-final | Fresh18-record t=1 forcing packet, inputs/history/metrics and selected native checkpoints |
| forcing-extended | Fresh18-record t=2 extension, unchanged production source and selected native checkpoints |
| forcing-debug | Actual StageForcing oracle on Debug with bounds checking |
| wave-final and wave-debug | Independent sampled longitudinal C2 pressure branch and four reservoir convergence metrics |
| filter | Ten synthetic analysis tests, conditioned/signed metric sample and hashes |
| clean-regressions | Fresh archive build's original regressions, exact scalar restart and uniform seed test |
| mpi-direct | MPI-enabled direct single-rank actual forcing oracle; no multi-rank claim |
| performance and performance-history | Final three-repeat microbenchmark and earlier contended/RSS-denied attempts |

The original unforced archives retain checkpoints. For the larger forcing packets,
initial/final and the source restart checkpoint are retained per run; intermediate
generated arrays omitted from archives are indexed by size/hash in
`raw-file-index.json`. Ordinary float32 `.bin` arrays are excluded. Numerical budget
claims use independent native-double `.rst` inventories; raw JSON/history metrics
are retained for every result. Original unarchived generated files remain at their
recorded temporary paths and are not promised permanent storage.

The retained independent reports distinguish mathematical derivation, source review,
build/run evidence and scientific interpretation. `code-patch-series.mbox` contains
the focused local patches relative to base07c5292; `code-commits.txt` records the
full local SHA list. Patches require the original repository and its pinned Kokkos
submodule, not hidden workspace files. No compiled executable or buildtree is
bundled. Cache/configuration/build logs record the actual environment.

The clean build reproduced commit14620bb3. The final acceptance/style patch adds
wave passive/fallback checks and explicit standard includes/line wrapping; its
production formulas and execution order are unchanged. Earlier tested source hashes
remain their original values. A final rebuild and smoke result is recorded separately
instead of assigning earlier measurements to a future commit.

MPI multi-rank execution failed because this machine's launcher and PMIx library
are incompatible. CUDA/HIP are unavailable. Peak RSS measurement was denied by the
execution sandbox. No environment reinstall or system modification was attempted.
Flake8 was unavailable; focused C++ changed-line and Python syntax results do not
imply a complete lint pass. No stationary turbulence or closure-benefit claim follows
from these capped, single-seed pilots.
