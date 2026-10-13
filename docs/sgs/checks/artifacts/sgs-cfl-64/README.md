# Matched SGS CFL sensitivity, 2026-10-12

**Three 64³ PPM4 branches finish 4.5T and pass 45/45 native checks.** All start
from the identical native4T checkpoint with unchanged fields, counter origins,
OU state and RNG. Only CFL changes0.1→0.05→0.025; each ends at the exact same
physical time 2.25 (T=.5). This is a 0.5T timestep-sensitivity experiment, not a
full 2T–10T statistical convergence campaign. No source equations or coefficients
were changed, and no rebuild was needed.

| CFL | Mean K,4T–4.5T | Exact4.5T K | Native steps | Mean dt | Native/verification wall seconds |
|---|---:|---:|---:|---:|---:|
| .1 | .0378136422888 | .0367334627167 | 1597 | .000156543519 | 606.82 |
| .05 | .0378132455541 | .0367312226147 | 3193 | .0000782962731 | 1201.19 |
| .025 | .0378130786941 | .0367315831875 | 6386 | .0000391481365 | 2011.49 |

K is energy density averaged over the unit cube, normalized by rho0 V²=1.
Mean K integrates linearly between stored native histories and includes the exact
restart origin. History/checkpoint cadences are 0.1T; slight timing jitter can
contribute to extremely small mean differences. Endpoints instead use independent
native-double inventories at exactly4.5T. The restart driver recomputes the
initial timestep for the changed CFL before evolving, as verified in source.

Maximum mean difference versus CFL .025 is **.00149048%**. Maximum endpoint
absolute difference versus CFL .025 is **.00511693%**. The largest pairwise endpoint
spread (CFL .1 versus.05) is about **.00609826%**. These changes are far smaller
than the approximately 49% full-window paper SGS-energy discrepancy.
The frozen operational 1% mean/endpoint screen passes. It is not a scientific
acceptance threshold or a claim of full-window statistical convergence.

The endpoint sequence is slightly non-monotonic: .05 is below .025.
The analysis therefore refuses a Richardson-order estimate (null), rather than
claiming second-order convergence from three small differences. Shocks, variable
steps, backscatter limiting and coupled driving are present. Endpoint K-field
relative L2 differences versus.025 are .08739% and .04650%; density differences
are .06270% and .03051%. Stable means do not imply exactly identical fields.

All three final **underlying Fourier OU arrays and RNG states agree exactly**.
The realized physical acceleration includes density-dependent mean subtraction
and normalization, so it can respond to the small field differences; that is part
of the coupled timestep comparison. All executable/native source fingerprints
agree. Closure coefficients, gamma 1.001, initial c0, seed 42, fixed forcing inputs,
mesh, reconstruction and integration method are unchanged.

## Energy/limiter checks

Each native branch passes15 admissibility, mass/momentum/scalar, applied-forcing,
SGS-transfer and independent history/checkpoint checks. Worst absolute window
forcing/represented-energy residual is 3.79e−13; SGS-transfer residual is 1.95e−15.
No fallback or diffusion limiting occurs. Integrated stress-limit increments
are 1.147e−6,7.504e−7,6.206e−7, respectively; their sign is suppressed backscatter,
not removal of positive production. SGS shares of corrected heating are 69.761705%,
69.761415%,69.761396%. The mean-energy reservoir turns over approximately 6.26
SGS dissipation/depletion times during this interval (D/meanK); these are not
independent statistical samples or a claim about every cell's local time scale.

The test gives no evidence for a large CFL bias over the checked interval.
It strengthens the direct RK2-weight/production checks, but does not establish
full 10T timestep independence or solve the SGS-energy accuracy question.
Spatial flux damping persists as timestep decreases. The separate resolved
spectra lack a broad clean Kolmogorov range; matched high-resolution reference
and forcing-convention comparisons remain open.

[Six-panel CFL report](report/cfl-comparison.pdf), [numeric analysis](report/results.json),
[restart wrapper](restart-cfl.py), [analysis script](analyze-cfl.py).
Each cfl directory contains the actual native command, unused generated input,
exact restart origin, native history/log, results, manifest and a protocol written
before execution. **The actual command/header, not the unused input file, defines
the restart physics.** Native resource caps were900/1800/3600 seconds, two OpenMP
threads per branch. Every native attempt completed within its original cap.

All runs executed on the Mac with the existing serial-double, CPU/OpenMP Release
sgs_verify binary, SHA256 b5f0fc19de082e07bc45c989cc24d5cc0c46dfd4d62e7994d8dd1efaac901197.
Source worktree is /Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu,
commit f9b3d596; native C++ matches 7bf794ad. No source checkout mutation or push.
The wrapper uses existing ForcingVerification, not a new runner or dependency.
For repetition, choose fresh output directories or edit the wrapper's artifact-root
path; preserved directories intentionally fail rather than overwrite native runs.
The direct native invocation pattern is:

```sh
OMP_NUM_THREADS=2 /private/tmp/athenak-sgs-heating/src/athena -r /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-heating-budget/turbulence64-ppm4/continued/rst/continued.00040.rst job/basename=continued time/tlim=2.25 time/cfl_number=0.05
```

Run inside a fresh output directory. The existing Python environment is
/private/tmp/sgs-plot-venv/bin/python, with MPLCONFIGDIR=/private/tmp/sgs-mpl-cache.
Raw checkpoints remain external; the compact wiki packet records their paths,
byte counts and hashes. Derived mean/endpoint/forcing comparisons were checked
independently before archiving. Analysis scripts pass Ruff; native source remains
unchanged. This experiment was not used to retune any coefficients or thresholds.
