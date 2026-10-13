# Native C2 and mixed SGS production investigation

Executed locally on 2026-10-11. **116/116 specialized checks and 720/720 native
invariant checks pass across 40 one-step 3D runs. No production solver defect
was demonstrated; the approximately 49% low turbulent SGS energy remains unresolved.**
Source commit `f9b3d596` adds one test script and local run documentation only.
The native C++ source hashes match the completed 64³ 10T campaign exactly.
No coefficients, forcing, evolution, restart format or solver machinery changed.
No push or scientific release approval.

## What was checked

The [test script](verify_mixed_production.py) reuses the existing problem generator's
3D initial-state loader and `Verification` runner. No rebuild is required for an
existing `sgs_verify` binary. It uses 16³ and 32³, PLM and PPM4, two OpenMP threads,
one complete RK2 step dt=1e-7, and the pre-existing rate tolerance 2e-5. Summed
per-case native execution plus verification took 5.21 seconds; initial-state
file generation and plotting are additional. Protocols are written before runs.

Four full velocity-gradient tensors exercise compression, expansion, trace-free
strain and both signs of C2 transfer. Uniform density/K in the central affine
core permit independently calculated local rates. The surrounding periodic
seams are intentionally discontinuous and excluded from the local oracle.
Only abs(x),abs(y),abs(z)<.125 is compared analytically. These are manufactured
amplitude tests, not physically representative turbulent solutions.

The analytic decomposition of the paper stress contracted with strain is:

```
Q = A:A, d = trace(A), S = sym(A)
P1 = 2 C1 Delta sqrt(2 rho K) [S:S - d²/3]
P2 = -2 C2 K [(A A^T):S/Q - d/3]
Pc = -(2/3) K d
Sigma = P1 + P2 + Pc
local dK/dt = Sigma - K d - Cepsilon K^(3/2)/(sqrt(rho) Delta)
```

P2 is the **deviatoric** C2 contribution, including its isotropic compensation.
Pc is the full trace compression. This avoids mixing part of C2 into a term
incorrectly labeled compression. Matching C2=0/.7 native runs isolate P2 on
identical initial states. The full local K rate also tests C1, conservative
advection and dissipation together at the published coefficients. Core rho/K
pairs are 3/.24 for compression/expansion and 1/.01 for trace-free fixtures.

A smooth periodic fixture uses u=m sin(2pi x), m=(.3,.2,.1), and
K=.02[1+.4cos(2pi x)]. Density is one. With C1=Cepsilon=Ckappa=0, its exact
shared-face global production is

```
KG = .02 * .4 * .3 * cos(pi/N) * sin(pi/N) * N
P(C2) = -2[(1-C2)/3 + C2] KG
P(.7)-P(0) = -(4/3) .7 KG
```

All three velocity components contribute traction work. Both signs occur
locally, while the nonzero signed volume transfer checks native production
counters, independently read K inventories, opposite resolved E change and
independent resolved kinetic stress work. C2=0 is a test control, not a changed
production default. Every run must complete one RK2 step without any limiter
or fallback activation; all do. Global conservation alone would not pass the
independent production-amplitude comparisons.

Maximum specialized relative rate error is **1.86452e-6**, below the frozen 2e-5.
Independent rank-one traction, coordinate rotation and density scaling reference
self-checks pass. [Final native results](final/results.json),
[specialized results](final/mixed-results.json), [summary](final/mixed-summary.json)
and [frozen protocol](final/mixed-protocol.json) preserve every value/tolerance.

![Four native verification panels](final/native-mixed-production.png)

[Printable native verification](final/native-mixed-production.pdf).

## Actual 64³ campaign suppression

The existing matched PLM/PPM4 histories give exact native applied-counter
increments over 2T–10T, with endpoints clipped linearly. This uses the same
qualified pressure-ledger histories, not a new turbulence run.

| Case | Native production P | Stress limiter counter increment | Counter / net P | Fallback stages | Diffusion limiting |
|---|---:|---:|---:|---:|---:|
| PLM | 2.89564511 | 8.18600948e-7 | 0.0000282701% | 0 | 0 |
| PPM4 | 3.65474406 | 8.45290363e-6 | 0.000231286% | 0 | 0 |

The stress limiter acts only on negative production/backscatter. Its direct
suppression **increases** net SGS production; it does not remove positive
production. These tiny integrated effects cannot account for the approximately
49% energy shortfall in this campaign. No fallback or diffusion limiter is
concealing the missing reservoir. These statements refer to measured counters
and the audited limiter sign, not to hypothetical alternative trajectories.

The [snapshot diagnostic](budget/production-budget.py) reconstructs the actual
shared-face gradient formula on nine native snapshots per scheme, at 2T through
10T. C1, C2 deviatoric and full-trace compression are all present; C2 is the
largest positive component at each sampled time. Instantaneous production and
dissipation are comparable. This rules out an absent C2 contribution in these
snapshots, but does not establish the accuracy of production relative to a
resolved turbulent reference.

Snapshot components are **not** an exact decomposition of accumulated RK-stage
counters: snapshots omit intermediate stage fields. The limiter curve uses the
saved timestep as a proxy, separately from the exact native integrated counter.
No nine-snapshot quadrature is presented as the full-window native budget.

Three independent snapshot decomposition checks pass at 1e-12: the periodic
rank-one analytic sum, simultaneous rotation of coordinates/velocity components,
and density/energy scaling. Maximum absolute error is 1.73e-17.
[Oracle checks](budget/oracle-checks.json), [snapshot/counter results](budget/results.json).

![Four production-budget panels](budget/production-budget.png)

[Printable budget diagnosis](budget/production-budget.pdf).

## Retained failed oracle and corrected execution

The first attempt had 96/116 specialized checks passing; all 40 native invariant
records passed. Its periodic expected value incorrectly retained only the normal
velocity's nonlinear work, omitting the transverse components. The error was
in our independent oracle, not in the native solver. The full 3D affine-core
checks already passed in that attempt. [Original expected values and failures](first-attempt/mixed-results.json)
and [executed original oracle](first-attempt/executed-oracle.py) remain unchanged.

The analytic periodic expression was corrected to include all three components,
without changing a tolerance or native code. A fresh corrected run passed
116/116. Independent reference self-checks were then added to prevent this error,
and a fresh final complete 40-case run also passed 116/116. All three attempts
have separate protocols, native inputs/histories/logs and results. Original
failure is not relabeled a pass. No native C++ rebuild was performed.

## Reproduce locally

Source worktree: `/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu`.
The tested existing binary is `/private/tmp/athenak-sgs-heating/src/athena`,
SHA256 `b5f0fc19de082e07bc45c989cc24d5cc0c46dfd4d62e7994d8dd1efaac901197`.
CPU double Release, OpenMP, MPI off, problem `sgs_verify`.

```sh
SGS_SOURCE=/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu
/private/tmp/sgs-plot-venv/bin/python -B "$SGS_SOURCE/tst/sgs/verify_mixed_production.py" --self-test
OMP_NUM_THREADS=2 /private/tmp/sgs-plot-venv/bin/python -B "$SGS_SOURCE/tst/sgs/verify_mixed_production.py" --exe /private/tmp/athenak-sgs-heating/src/athena --output /private/tmp/sgs-mixed-local-check
```

Use a fresh output directory. If the temporary binary/environment are gone,
build `PROBLEM=sgs_verify` and use the NumPy environment as documented in
`tst/sgs/README.md`. No new dependencies were added. `final/plot-results.py`
regenerates the native figure from its JSON with installed Matplotlib.
`budget/production-budget.py` requires the external saved 64³ checkpoints and
existing full-window reports; it is a campaign-specific artifact, not a generic
new reporting framework. Its output directory must be fresh.

`archive-manifest.json` records the local source commit, executable/native source
identity, byte-verified compact copies and external hashes for every native
checkpoint and generated initial-state array. Large raw files remain under the
artifact root `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb`.
New source Python lint and whitespace checks pass. Native history headers retain
their original trailing spaces for byte identity; they are excluded from prose/code
whitespace checks. Unrelated source/wiki edits were preserved.

## Verdict and remaining diagnosis

A missing C2 amplitude/sign, paired-exchange error, or substantial limiter/fallback
suppression is **not demonstrated by these tests**. The narrow native mixed-production
gap is closed within these manufactured fixtures. No solver correction is justified
by the findings, and no coefficient or forcing change was made to obtain agreement.
The 64³ PPM4 SGS energy is still .0367689 versus the paper's .0726. This experiment
does not validate quantitative turbulent energy, density-dependent SGS transport,
reference closure correlation, resolution scaling or scientific release gates.

The next discriminating experiments remain controlled forcing conventions and a
sufficiently resolved Gaussian/Favre turbulent reference. Native conservative
Euler dissipation can compete with SGS production; matching analytic amplitudes
cannot establish that the resolved gradients capture the proper cascade near the
cutoff. Any change to forcing or gradients must be a separate, measured control.

Primary equations: [Schmidt–Federrath (2011), Eqs. 7–10](https://arxiv.org/html/1010.4492).
[Completed 2T–10T comparison](../2026-10-11-full-window/README.md).
