# SGS timestep sensitivity and cascade comparison

Executed locally on 2026-10-12. **The matched64³ CFL check finds very small SGS
energy changes; Kolmogorov-like two-grid scaling does not establish absolute SGS
energy accuracy.** The approximately 49% paper-energy discrepancy remains open.
No source equations, coefficients, forcing semantics or reconstruction changed.
Source worktree remains clean at f9b3d596; native C++ matches7bf794ad. No push.

## CFL results: identical starting state and stochastic realization

Three 64³ PPM4/RK2/LLF branches start at exactly4T from one pinned native checkpoint
and end exactly4.5T. Actual mean timestep approximately halves with CFL; native
step counts are 1597,3193,6386. All underlying Fourier OU arrays and RNG states
are exactly identical at the endpoint. Applied physical acceleration can respond
to the tiny density changes through its existing normalization/mean subtraction.

| CFL | Time-weighted mean K,4T–4.5T | Common-time endpoint K |
|---|---:|---:|
| .1 | .0378136422888 | .0367334627167 |
| .05 | .0378132455541 | .0367312226147 |
| .025 | .0378130786941 | .0367315831875 |

Fourfold CFL reduction changes the mean by **.00149048%** relative to the finest
run. The largest pairwise endpoint spread is **.00609826%**. All 45/45 native
budget/admissibility/history checks pass. SGS corrected-heating shares remain
69.7614–69.7617%; no fallback or diffusion limiting occurs. K-field/density
relative L2 differences versus the finest run are below .1%. The endpoint sequence
is slightly non-monotonic, so no formal Richardson convergence order is claimed.
The frozen operational 1% screen passes; it is not a scientific-release criterion.

This is a 0.5T response experiment, not a full 2T–10T statistical CFL campaign.
It covers approximately 6.26 global SGS-reservoir dissipation/depletion times;
that is not an independent-sample count. It gives no support for a large CFL
explanation of the observed factor-of-two energy gap over the checked interval.
Spatial numerical damping is not removed simply by taking smaller timesteps.
[Full protocol, costs and reproduction](cfl/README.md),
[six-panel CFL report](cfl/report/cfl-comparison.pdf),
[numeric results](cfl/report/results.json).

## Direct source-weight check

The RK2 SGS sources act inside the two stages. First-stage beta=1 is subsequently
halved by the second-stage state blend; stage 2 adds beta=.5. Final history weights
are .5/.5. A direct native uniform-decay audit matches the exact discrete Heun
recurrence over one and 11 steps to maximum absolute specific-k error 1.90e−15.
A deliberately half-rate source produces clearly different endpoints. The earlier
native mixed-production checks establish production amplitudes separately.
This rejects a uniform missing source multiplier .5 in the checked supported path,
not every possible turbulent timestep/transport error.
[Exact audit results and retained harness failures](rk2-audit/README.md).

## Kolmogorov expectation and measured scale dependence

For a homogeneous incompressible inertial range, the radial energy spectrum
E(kappa)=C_K epsilon^(2/3) kappa^(−5/3) integrates above a sharp cutoff to

```
k_sgs = (3/2) C_K epsilon^(2/3) kappa_c^(−2/3)
kappa_c = pi/Delta  =>  k_sgs proportional to (epsilon Delta)^(2/3)
```

The measured PPM4 mean K decreases from .05906799 at 32³ to.03676890 at 64³ over
2T–10T, with matching physical input parameters. The two-point exponent is .68389,
near 2/3. Anchoring a constant-flux K41 curve to our 32³ point predicts64³ K=.03721050,
about1.19% above the actual value. This is a relative-scaling check, not an
independent absolute-energy prediction:32³ is below the paper grid floor, only
two resolutions/one seed exist, and flux/density statistics change with resolution.
The [SGS paper §4.2](https://arxiv.org/html/1010.4492#S4.SS2) reports a solenoidal
energy exponent .799±.009 over 64³–512³ and warns that SGS energy is not directly
related to velocity spectra or structure functions.

Absolute K41 numbers depend on flux, density and filter conventions. Choosing
3D C_K=1.6 and a sharp Nyquist cutoff gives .06582 using our SGS dissipation,
or .08372 using the independently measured total compression-corrected-heating
proxy. Measured K is .03677. A chosen Gaussian filter gives yet another ideal
estimate. These are conditional incompressible estimates, not acceptance targets
for a Mach≈5.5 compressible finite-volume flow. The heat proxy is not a measured
inertial flux. Even the published energy is lower than the homogeneous sharp-filter
estimate based on its quoted dissipation. The one-dimensional longitudinal
constant near .5 is not the 3D radial constant; see the primary
[Sreenivasan (1995) paper](https://users.ictp.it/~krs/pdf/1995_003.pdf).

Our SGS dissipation ratio relative to the conditional density-weighted paper value
is .36944; its2/3 power .51486 is close to the energy ratio .50646.
However, SGS dissipation is computed from K using the same closure. This is
internal consistency, **not an independent explanation or validation of low K**.
The model-based homogeneous/Jensen bound is .044896; observed density/k
intermittency makes the mean smaller. The inferred moment factor 1.349 is comparable
to the paper's conditional1.316. Table4's epsilon header and density-weighted prose
remain ambiguous; neither convention is silently relabeled.
[Detailed derivation, assumptions and checked Gaussian integral](kolmogorov/README.md),
[four-panel comparison](kolmogorov/kolmogorov-comparison.pdf).

## Does the resolved flow exhibit a clean cascade?

Nine equally weighted snapshots from approximately 2T–10T were Fourier analyzed.
The u and rho^(1/3)u spectra preserve Parseval to 4.34e−16 relative error.
Shells4–10 were frozen before analysis; fits are descriptive and cover less than
half a decade above forcing support, where numerical damping can overlap.
At64³ the slopes are −2.746 for u and−2.180 for rho^(1/3)u; at 32³ they are steeper.
Compensated curves show no broad flat K41 range. This is consistent with limited
scale separation and spatial damping, but does not isolate their cause or prove
an SGS defect. Density weighting, shocks, SGS stresses and forcing can all affect
these finite-grid spectra. No tail extrapolation is used as SGS ground truth.
[Procedure, limitations and spectra](spectra/README.md),
[four-panel spectrum report](spectra/resolved-spectra.pdf).

## Scope and reproducibility

The simplest useful work reused the existing Verification/ForcingVerification,
restart reader, statistics and heating helpers. No native instrumentation,
new dependency or generic reporting framework was added. Native CFL runs used
fixed coefficients, gamma 1.001, seed 42, the exact same 4T checkpoint, two OpenMP
threads each and original 900/1800/3600-second caps; all completed within cap.
Fine-branch elapsed native/verification time was 2011.49 seconds. This remains a
Mac CPU/serial-double support-envelope test; larger campaigns are unexecuted.

All detailed protocols, native commands, exact restart origins, logs, histories,
results, manifests and analysis scripts are retained. A generated input file is
not the restart physics: the actual command plus native header determine that.
`archive-manifest.json` verifies compact copies by SHA256 and records external
checkpoint paths/byte counts/hashes. Original failed RK2-output/duplicate-row
analysis attempts are preserved separately; all three CFL runs passed initially.
No numerical tolerances were relaxed. Analysis Python lint, independent math,
forcing parity and copied-byte checks pass. Native history headers preserve
original trailing spaces and are excluded from prose/code whitespace checks.

Verdict: no large CFL effect in this matched response test; provisional2/3
relative scaling, no certified Kolmogorov inertial range, and no independent
validation of absolute SGS energy. Same-forcing high-resolution Favre-filtered
reference, forcing conventions, full statistical timestep/resolution/seed
qualification and the paper-energy discrepancy remain unresolved.
