---
title: AthenaK SGS hydro implementation and qualification
class: living
created: 2026-10-08
last_touched: 2026-10-08
---

# AthenaK SGS hydro implementation and qualification

The Schmidt–Federrath HD SGS model is implemented on branch
`codex/schmidt-federrath-hd` in `/Users/hitesh/hitesh/git/athenak_nyu`.
Restricted serial CPU cases pass independent conservation, source, restart and
convergence checks. **The M0–M9 baseline scientific release is incomplete.**
No AMR, MHD, cooling or dynamic-coefficient extension has been added.

The [development plan](development-plan.md), [catalogue](test-matrix.md) and
[runbook](agent-runbook.md) remain the governing planned coverage. Their original
checklists are preserved. The [family ledger](coverage-ledger.md) records every
catalogue ID; an executed instance does not close all instances in a family.
The runnable code documentation is
`/Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/README.md`.

A separate [port to multiphase_turb](multiphase-port.md) is implemented and
qualified on selected serial CPU cases. Its branch, source revision and new evidence
are recorded separately; the historical results and ledger below remain unchanged.

## Source and ownership

The confirmed base was `gotham_development` at
`07c529202befbeb6ef5e694d95286929618ff7c1`. Existing `.gitignore` modifications
and untracked `AGENTS.md` were preserved and excluded from code commits. The
coordinator owned shared hydro/EOS/task/output files; independent workers reviewed
science and lifecycle, reproduced the baseline, and owned verification packets.
Each build/output directory was separate. Unrelated wiki edits were preserved.

| Local code commit | Change |
|---|---|
| `4f5dbb33d541c603ba8c5552c5cc71800466bde4` | Raw closure helpers and independent Kokkos oracle |
| `36be17da5c3335795c287986cbd437ce786a6390` | SGS state, conservative transport, coupled RK2 flux/work/sources and deterministic runner |
| `d6c5512edbb70fb38b650c78e28044a73a564325` | Stage forcing, complete OU restart state, source and sampled pressure-wave oracles |
| `14620bb33d580246e32d80f4bd584c6c021bbd1e` | Filtering analysis, examples and support documentation |
| `6ddfedd46cca0e670e19e4d82f3725824be6bc15` | Final wave passive-fraction/fallback acceptance checks and source style; production formulas unchanged |

These are local commits. No push, PR, merge or deployment was performed. Per-run
manifests bind the tested executable to source hashes and the recorded source
revision/diff; several runs predate their final commit and are not falsely labelled
as tests of a future SHA. The original reports and raw checkpoints are retained
in [the evidence bundle](evidence/2026-10-08/README.md).

## Milestone status

| Milestone | Actual status and gate limits |
|---|---|
| M0 | Untouched Debug/Release baseline reproduced; selected hydro/scalar/restart/turbulence tests and pre-existing limitations recorded. Full regression harness was not run. |
| M1 | Independent equation/stage/energy audit; restricted numerical contract frozen. The nonlinear zero-gradient discontinuity prevents a general wave/stability proof. |
| M2 | Serial double raw closure oracle and mutation detection pass. GPU, single precision and some transformation/invalid-input families remain unqualified. |
| M3 | State/ghost/reconstruction/scalar/output/restart paths implemented; selected lifecycle tests pass. Uneven packs, repeated destruction and complete reader/format coverage remain open. |
| M4 | Coupled transport/stress/work/production/decay implemented; budgets, decay order, component time refinement, actual forcing oracle and sampled mixed pressure-wave branch pass. Full compression/expansion and manufactured PDE gates remain open. |
| M5 | Strict admissibility, shared-face caps and global donor-cell fallback implemented; tested collisions, boosts, contrasts and CFL sweeps pass. Active stress limiting/fallback observed; active diffusion limiting, rarefaction/near-vacuum and general robustness remain open. |
| M6 | Serial block/restart/history and direct one-rank MPI checks pass. MPI launcher mismatch blocks multi-rank execution; CUDA/HIP unavailable. Performance is characterized only on this Mac. |
| M7 | Independently reviewed Favre/top-hat analysis and ten synthetic tests pass. No turbulent reference with demonstrated converged scales was supplied or produced within the capped runs. Closure accuracy and uncertainty gates remain incomplete. |
| M8 | Measured-work SGS-on/off solenoidal/compressive pilots pass numerical budgets. Initial pilots are unsaturated; no justified stationary multi-resolution, multi-seed reference campaign is claimed. |
| M9 | Focused local patches, inputs, runners, docs and evidence delivered. Scientific release gate remains blocked by the preceding incomplete gates. |

## Numerical contract

The scientific baseline is [arXiv 1010.4492v1](https://arxiv.org/abs/1010.4492v1),
the only submission listed by arXiv when checked on 2026-10-08, associated with
the Schmidt & Federrath (2011) A&A article. The independent review checked the
[full text](https://arxiv.org/html/1010.4492), Eqs. (7)–(14), coefficients and
conventions. Its derivations are retained in the evidence bundle.

K is SGS kinetic energy density; primitive k=K/rho. Resolved E retains ordinary
thermal-plus-resolved-kinetic semantics and excludes K. K is appended after true
passive scalars; their indices, counts and output meanings remain intact.
The raw stress is

```
A_ij = derivative_j u_i; Q = sum A_ij^2; S* = sym(A) - trace(A) I/3
tau = 2 C1 Delta sqrt(2 rho K) S*
      - 2 C2 K A A^T/Q - (2/3)(1-C2) K I
Sigma = tau:S; R = Ceps K sqrt(K/rho)/Delta
J = -Ckappa Delta sqrt(rho K) grad(K/rho)
```

Defaults are C1=.02, C2=.7, Ceps=1.5 and Ckappa=.65. Raw trace is -2K.
At exactly Q=0 use isotropic -2K I/3; nonzero-Q limits depend on direction.
At K=0 the closure vanishes. Scaled gradient normalization avoids avoidable
underflow/overflow without a denominator epsilon. No repeated seed or SGS floor
is used. Fresh uniform initialization adds K=rho*seed once, without subtracting
thermal energy; exactly zero uniform K remains zero.

The tested subset is Newtonian ideal-gas double HD, periodic uniform cubic cells,
RK2, LLF, DC/PLM, at least two ghost zones and CFL≤.2. Dimensional reductions
retain three velocity components and the 3D trace. Unsupported physics, AMR/SMR,
GPU/single precision, other solvers/integrators and incompatible hooks fail clearly.

One shared face traction enters momentum as -tau.n and resolved work as
-u_face.tau.n. Arithmetic face velocity and the final applied traction give
`Sigma_h=div(u_face.t)-u_cell.div(t)`. E and K receive opposite stage source
exchanges. SGS diffusion is a represented-energy flux, not an artificial E sink.
The semidiscrete thermal SGS source is dissipation; finite-step kinetic partition
has convergent RK error and LLF numerical heating remains distinct.

Full-dt donor caps scale stress/work consistently and cap outgoing diffusion.
Scaling stress changes its applied trace; only the raw algebra retains -2K.
Any failed full-dt SSP predictor globally rebuilds hyperbolic fluxes and closure
data with donor-cell reconstruction. A second failure stops; EOS floors do not
silently repair an SGS run. History measures actual production, dissipation,
suppressed stress/diffusion transfer and fallback stage count.

Timestep estimates include thermal-plus-SGS regularization 6k, receiving-density
eddy/diffusion row sums and decay. This is a finite-amplitude frozen-shape estimate
with actual stage checks, not an arbitrary-flow stability theorem. An attempted
k/|A| restriction stalled contacts/diffusion and was removed after review. Closure
stencil differences below eight machine epsilons times velocity magnitudes are
treated as zero; actual differences remain in compatible production. Boost tests
characterize finite conditioning, not exact frame invariance.

Stage-consistent forcing uses old-stage rho*a and rho*u.a without repeated quadratic
energy kicks. Heun supplies the quadratic constant-acceleration work. Cadence and
thermal kick caps can change applied work relative to nominal dedt. History records
actual work and impulse separately. Marked restarts preserve raw OU memory, RNG,
exact OU counter and precise accumulated budgets. Every forcing parameter string
is frozen across restart; equivalent reformatted overrides are also rejected.
Legacy forcing behavior and payloads remain unchanged when the new mode is off.

## Executed results

| Packet | Result and interpretation |
|---|---|
| Untouched baseline | Debug/Release builds; selected linear waves, scalar conservation, normal restart and Sod/turbulence smoke. Wave 32/16 error ratio .278851<.29. |
| Kernel oracle | 1025 cases, 22575 assertions; maximum scaled error 1.45213e-15 at tolerance 5e-13. Three temporary sign/factor mutations were detected. |
| Revised bounds Debug matrix | 43 result records pass; worst E+K drift 2.94099e-16, paired-transfer residual 6.93889e-18. Later fallback restart and active stress-limit cases separately pass. |
| Final Release unforced matrix | 46 records pass, zero fail. Includes actual fallback/stress-limit and normal checkpoint restarts. |
| Configuration | 22 checks pass, including bitwise identical omitted/disabled states and early rejection tests. |
| Decay | Independently derived K(t) and deposited heat; time orders 1.95053 and 2.00308. |
| Conditional stress refinement | Full-solver field and kinetic/thermal/SGS component time orders about 2.02 and 2.01, with C2/decay/diffusion isolated off. Numerical reference, not exact full PDE. |
| Mixed pressure branch | C2=.7, C1=Ceps=Ckappa=0; b=1.6, SGS acoustic factor4.16. Phase-shifted 16/32-cell density/u/P/K errors meet normalized fine L1<.1 and ratio<.4 (observed .255–.315); scalar fractions retained and no fallback. Release and bounds Debug pass. |
| Actual forcing kernel | Constant acceleration (.7,-.2,.3), dt=.01; every field and impulse/work satisfy independent Heun expectations at1e-12; max error5.09141e-16 in Release and bounds Debug. |
| Fresh forcing packet | 18 records pass: four pilots, three continued lanes, three bitwise restart comparisons and eight contract rejections. Passive masses agree exactly in all physical lanes. |
| Clean independent reproduction | Fresh archive build at14620bb3, closure, unchanged SGS-off waves/Sod/scalars and exact uniform seed pass. |
| Filtering | Ten independent synthetic tests pass at1e-12. Signed/conditional error and sign metrics retain unfavorable results; no coefficients retuned. |
| MPI | MPI build and direct one-rank execution pass; selected active state equals serial bitwise. Multi-rank launcher fails due installed PMIx6.1 versus PRRTE's PMIx5.0.9 ABI. |

Strong collisions v=5 and v=20 activated 9 and19 fallback stages. A four-cell
high-K collision activated suppressed stress transfer .2126880755333241 with
closed budgets. Two bounded high-diffusivity attempts stayed admissible and closed
their budgets, but had exactly zero suppressed diffusion. Active diffusion limiting
is therefore still unexercised.

The four initial 16³ pilots used ordinary gamma5/3, seed42, nominal dedt=.1,
tcorr=.1 and update cadence .01 through t=1. Applied work ranged .09563–.09990.
Late kinetic fractional drift was22–26% over about t=.5–1. Seed/nominal-power
matching does not imply equal realized forcing work or replayed forcing. Those
runs cannot support stationary spectra, uncertainty or model-benefit claims.

The synthetic n8,width3 filter comparison has signed production correlation .6949,
relative RMS error1.6947 and false-sign fraction .41016. It is a constructed
trigonometric field, not validated turbulent data. Top-hat Delta=width*dx is a
declared filter convention, not calibrated equivalence to the paper's Gaussian.

## Extended pilots and cost

A separately capped extension through t=2 passes18/18 records on the identical
production source and executable. Worst applied-work E+K residual is6.384e-16;
three restart comparisons remain bitwise. Late fitted kinetic drift over t≈1–2
is12.9–17.6%. The compressive late kinetic range is91–92% of its window mean.
These strong fluctuations and remaining drift do not establish stationarity.
Only seed42 and16³ were run; no three-seed or32³ scientific campaign was launched.
The extended report retains actual work, heat, SGS and time-window definitions.

An idle-build serial Release microbenchmark uses three repeats of the same16³
50-cycle unforced jobs with outputs disabled. Median Athena loop time is .06852s
SGS off and .15666s on (about2.29 times per cell-cycle); reported throughput is
2.99 million and1.31 million zone-cycles/s. Physical timesteps differ, so this
is computational overhead characterization, not equal-physical-time accuracy
or performance portability. Earlier contended timing attempts remain recorded.
Peak RSS could not be measured because the sandbox denied the timing utility's
sysctl call; no measured memory-overhead claim follows.

## Failed attempts and remaining gates

The untouched scalar restart after a deliberately final-time-clipped run changed
the timestep schedule (201 versus161 cycles; max scalar difference1.0763e-6).
A normal intermediate checkpoint was bitwise reproducible. This pre-existing
distinction is preserved; no generic legacy timestep reset was introduced.

The first coupled matrix had20 passes and6 completion failures from the provisional
gradient-based timestep. Revised runs retained the same acceptance tolerances.
An earlier Debug box exceeded its20-second cap; a separately declared short box
passed. Include-order and reader-label setup failures remain in the archived logs.

Initial forced restarts reproduced RNG bytes but failed fields/raw OU. Native
no-advance readback proved the payload intact. Default kpeak=4pi had been saved as
12.5664, changing OU amplitudes; precise serialization fixed the regenerated
checkpoint comparison. Incorrect expected-error messages in two harness attempts
were corrected under new run names. No physical failure was relabelled a pass.

Open local numerical work includes exact full mixed manufactured solutions,
homogeneous compression/expansion PDE oracles, active diffusion limiting, additional
rarefaction/low-density and reconstruction-trigger cases, complete output-format,
pack and performance coverage. General normalized-gradient stability remains a
model/numerics question. MPI/GPU qualification and a converged turbulent reference
require additional execution resources. A statistical campaign requires demonstrated
stationarity, multiple resolutions/mixtures/seeds, common-scale reference comparison
and correlated-sample uncertainty. These gates remain explicit; M9 is not complete.

## Reproduction

From the recorded checkout/commit, use fresh build/output directories:

```sh
cmake -S /Users/hitesh/hitesh/git/athenak_nyu -B /tmp/sgs-repro \
  -DPROBLEM=sgs_verify -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=OFF
cmake --build /tmp/sgs-repro -j3
ctest --test-dir /tmp/sgs-repro --output-on-failure
python3 -B /Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/verify.py \
  --exe /tmp/sgs-repro/src/athena --output /tmp/sgs-repro-unforced
python3 -B /Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/configuration.py \
  --exe /tmp/sgs-repro/src/athena --output /tmp/sgs-repro-config
python3 -B /Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/verify_wave.py \
  --exe /tmp/sgs-repro/src/athena --output /tmp/sgs-repro-wave
python3 -B /Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/verify_forcing.py \
  --exe /tmp/sgs-repro/src/athena --output /tmp/sgs-repro-forcing
```

`inputs/hydro/sgs_verify.athinput` supplies the deterministic cases;
`inputs/hydro/sgs_box.athinput` supplies the short forced example. The independent
field reader uses native-double restart arrays: ordinary binary field output is
float32 and cannot substantiate1e-10 inventory tolerances. NumPy is needed only
for offline filtering tests; no solver dependency was added. Flake8 was unavailable;
syntax checks and focused C++ checks are recorded without claiming that lint ran.
