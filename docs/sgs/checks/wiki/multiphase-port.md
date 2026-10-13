---
title: Schmidt-Federrath HD SGS port to multiphase_turb
class: living
created: 2026-10-08
last_touched: 2026-10-12
---

# Schmidt-Federrath HD SGS port to multiphase_turb

The restricted uniform-grid HD SGS implementation is now on
`codex/multiphase-schmidt-federrath-hd`, based directly on `origin/multiphase_turb`
at `526a187f938e38634e89a88192ebe9731c4d36bb`. Local code commit:
`32c03c1a581974e5edb76e7da2f2eed29ff0f63d`.
The attached worktree is
`/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu`.
The original `codex/schmidt-federrath-hd` checkout and its user modifications were
preserved. No push, PR or merge was performed.

This is a port of the restricted model documented in
[implementation-progress](implementation-progress.md), not completion of the M0–M9
scientific gates. The previous evidence and family ledger remain frozen. Results
below were executed on this port; previous MPI/performance measurements do not
qualify this branch.

## What changed and why

| Change | Reason |
|---|---|
| Append SGS energy density K after true passive scalars; transport primitive k=K/rho | SGS energy needs its own reservoir without changing passive-scalar indices or resolved energy semantics. |
| Extend hydro arrays, boundary communication, reconstruction, stage copies and native checkpoint counts | Every state path must carry K. Mesh migration counts are consistent, but AMR/SMR remain rejected. |
| Published mixed stress, diffusion and dissipation; shared traction/work and paired E/K sources | Preserve the Schmidt–Federrath closure and conservative exchange between resolved and SGS energy. The raw closure header is identical to the preceding implementation. |
| Full-stage positivity checks, stress/diffusion donor caps, donor-cell retry and SGS timestep estimates | Prevent negative SGS or thermal reservoirs without silently repairing them with EOS floors. |
| Stage-consistent forcing adapted to the branch-native Fourier-mode OU renderer | RK2 forcing must use the same stage RHS as energy exchange, while reusing the newer forcing implementation. |
| Reuse native modal restart metadata, RNG and OU amplitudes | The older cell-force restart payload was deliberately not ported. The newer format already preserves the necessary OU state. |
| Use Driver::SourceTermHistoryWeight and SetRealPrecise | Follow the branch's source quadrature and serialization rather than duplicate them. |
| Append SGS/forcing history after existing scalar/cooling columns; capacity 32 with early bounds checks | Preserve existing labels and avoid scalar-column collisions or out-of-bounds writes. Twelve true scalars plus SGS and forcing fill all 32 columns in a tested case. |
| Append hydro_u_sgs and hydro_w_sgs selectors at indices 219/220 | Existing output validation uses hard-coded index ranges; inserting selectors earlier would break existing choices. |
| Reject newer cooling, hydro_srcterms, frame tracking and actual user source/cooling hooks | Those energy exchanges need separate qualification. Merely inspecting absent legacy cooling flags must not insert them and trigger the branch's migration checks. |
| Initialize uniform K after generic initial perturbations | Ensure K/rho equals the requested seed for the final initial density. |
| Update independent checkpoint reader, examples and verification runners | Test native modal state, source budgets and the branch's modern forcing parameters. |

The existing `hydro_sgs` diagnostic selector remains the branch's filtering-moment
output. It is distinct from evolved SGS energy. Its pre-existing energy-moment
interpretation was not corrected as part of this port; these SGS tests use
independently decoded state and explicit SGS selectors.

## Simplicity review and cleanup

Following the request to minimize changes to AthenaK's core machinery, the
[ponytail skill](/Users/hitesh/.agents/skills/ponytail/SKILL.md) was applied. The
follow-up local code commit is `6e5e2d4951126b21e03f4fe837319c7c1ee85770`.
The initial commit and qualification above remain historical records.

The integration uses existing hydro arrays and scalar transport, the unchanged
hydro task dependency graph and RK2 integrator, ordinary boundary communication,
existing native modal forcing restart state and the existing source-history weight
helper. A source-only hook cannot implement SGS momentum traction and energy work.
Treating K as an ordinary user scalar would alter scalar counts/identities and apply
scalar flooring rules. A separate SGS solver/communication/integration framework
would add machinery without satisfying a missing requirement.

The review found avoidable changes and removed them:

- Restore the entire AMR load-balancing source file to the multiphase base.
- Restore coarse-grid/shearing allocation, AMR flux-receive and RK4 register paths;
  SGS rejects those configurations before these paths execute.
- Keep the stage-forcing task-registration guard in its shared method, removing
  the duplicate check at its caller.
- Remove unused predictor arguments and an unused header dependency.
- Remove the never-read face-work storage component; work still enters the energy
  flux. Five rather than six components reduce that scratch allocation by one sixth.
- Append SGS fields after existing hydro/scalar aggregate-output columns. One check
  in the existing configuration runner verifies conserved and primitive column order
  with twelve scalars. Output filenames and explicit SGS selectors are unchanged.

Conservative stress/work/production coupling, timestep restrictions, admissibility,
restart validation and recorded numerical budgets were retained. These enforce the
requested invariants; removing them would produce a smaller but incomplete model.
There is no claim that this is a provable absolute minimum implementation.

Fresh Release checks pass: 46 unforced, 32 configuration, two acoustic grids/all eight
acceptance checks, 19 forcing/restart and one actual forcing oracle. Bounds Debug
passes the full 46 unforced matrix, both acoustic grids/all eight checks, actual
forcing oracle and closure. Six SGS-off built-in cases remain bitwise identical to
the untouched base; both generic seed checks pass. **All 395 saved native checkpoints
from the unforced, acoustic and forcing suites are byte-identical to the pre-cleanup
port**, including checkpoint headers, field data, counters and forcing state. The
aggregate binary output order intentionally changes to preserve existing columns.

[Cleanup evidence](evidence/2026-10-08-multiphase-ponytail/README.md) stores new
inputs/results/build logs/manifests and byte-level before/after checkpoint hashes.
The unchanged checkpoint bytes remain in the original port archive instead of
being duplicated. Added-line C++ checks and changed-Python parsing pass. The
qualification limits in this record are unchanged.

## Reviewer suggestions implemented

Follow-up code commit: `ac384d2d4f52c3dd3b6cb9201905f02bf69adad5` on the same multiphase SGS branch.
The three independent reviews found no justified change to the numerical core.
Their small suggestions are now applied: existing fatal-error helper reuse, one
combined guard, removal of the already-bounded final history capacity check, and
one shared input serializer. All pre-write capacity checks remain.

The additional runner defect is fixed: configuration attempts are saved before
launch, timeout/read failures are recorded, partial logs survive, and later checks
continue. A focused regression passes with the fix and fails against the previous
runner. It covers both subprocess sites, malformed checkpoints and missing states.

Fresh Release checks pass: 46 unforced, 32 configuration, both acoustic grids/all
eight acceptance checks, 19 forcing/restart, actual forcing oracle and closure.
Bounds Debug passes 46 unforced, 32 configuration, both acoustic grids/all eight
checks, actual forcing oracle and closure. Six SGS-off built-in cases remain bitwise
identical to the untouched base and both seed checks pass. All 779 directly paired
native checkpoints are byte-identical to the previous cleanup. Every one of the
795 newly saved packet checkpoints matches bytes in the original archived states.

[Reviewer-fix evidence](evidence/2026-10-08-multiphase-review-fixes/README.md) retains
the patch, results, logs, regression reproduction and fingerprint/archive checks.
No numerical equations, support rules or acceptance thresholds changed. The
qualification limits below remain unchanged.

## Laptop-sized 2D visual checks

Code commit `1c35d3eecaebae9339c7ea671b6fa05454d8925e` adds `tst/sgs/laptop_2d.py`, using the
existing problem generator and verification helpers. No C++ solver changes.
A final local serial Release run completes all 22 simulations / 25 result records
on 16² and 32² grids, producing 15 PNG figure pages with 92 panels, a PDF, and an
HTML gallery in 10.61 seconds including plotting.
The [PDF](evidence/2026-10-08-laptop-2d/report.pdf) includes analytic decay/advection,
shear/diffusion profiles, box evolution, SGS-on/off fields, applied energy budgets,
restart/decomposition checks, forcing, spectra and elapsed runtime/solver health.

Decay orders are 2.00925 and 2.00137. The contact density L1 error falls from
.00555285 to .00168528. Independent unforced E+K drift is at most 1.14492e-15 in
these plotted checks; the final forced represented-energy/work residual is at most
2.35922e-16. Restart and decomposition comparisons are exact. Reported elapsed
solver costs are about 4.0–4.8 times SGS-off for the displayed matched 32² cases;
these single short timings are not benchmark medians. Plot-only regeneration and
the calculation self-checks pass. Existing default 3D forcing still passes 19
records with all 31 saved checkpoint bytes unchanged.

[Evidence and local reproduction](evidence/2026-10-08-laptop-2d/README.md) preserve
all final raw outputs and the gallery. Build/run commands are in the repository's
`tst/sgs/README.md`. Changing resolution changes the SGS filter width; forced t=1
pilots do not demonstrate stationary turbulence or reference closure accuracy.
The model's support rules and outstanding scientific qualification are unchanged.

## Filtered high-resolution turbulence comparison

Code commit `b932230caa7ae749ae43309def77ee1232621dac` adds a matched resolution
comparison: drive a 128² SGS-off flow to t=2, then decay 128² and 64² references
and 32² SGS-off/on lanes to t=.5. Conservative block averages match coarse mass,
momentum, scalar masses and represented E+K. SGS-on stores removed kinetic
variance as K and preserves filtered thermal energy; SGS-off stores that variance
as heat, so their initial thermal pressures differ. Reference fields are filtered
to 32² before comparing diagnostics.

This case shows no overall SGS improvement. Final relative RMS density errors are
8.3407% off / 8.3133% on, velocity 40.2507% / 40.2996%, pressure 13.7325% /
13.5851%, and vorticity 41.3469% / 41.6584%. Modeled K has 83.8375% relative RMS
error against filtered 128² K. The filtered 64² and 128² velocities differ by
21.14%, so a converged reference is not established. This remains a single 2D
seed with three velocity components, followed by unforced decay.

All six comparison/initialization records pass numerical guards; that pass does
not require improved SGS accuracy. The Release and bounds Debug existing suites
pass unchanged, with 520 existing native checkpoints byte-identical. Native
initializer round-trip and eight malformed-input rejections pass in both builds;
a bounds Debug 64²/16² short comparison passes. The new source hook is confined
to the test problem generator; solver, task and restart machinery are unchanged.

The final workflow takes 35.35 seconds including nine plot pages / 41 panels.
[Report PDF](evidence/2026-10-08-resolution-comparison/report.pdf) and
[evidence/reproduction](evidence/2026-10-08-resolution-comparison/README.md) retain
all inputs, raw states, histories, comparisons, source fingerprints and checks.
The outstanding physical qualification below remains open.

## Three-dimensional resolution comparison

Code commit `eac5743afb44c53dd47dc463f4e7f014a0113ca8` extends the test-only loader
and existing comparison to 3D. Drive 64³ to t=2; initialize and decay 64³/32³
references and 16³ SGS-off/on lanes to t=.5 with fixed published coefficients.
The serial Release workflow takes 482.27 seconds including ten plot pages /47
panels. A 64²/16² control with matched parameters also runs (9.09 seconds).

SGS-on improves the 3D errors only slightly: density 9.3000% off /9.2540% on,
velocity 77.9380% /77.7077%, pressure 15.5361% /15.3569%, vorticity 91.3994% /
91.3690%. Modeled K is 41.77% of the filtered 64³ reference's mean. Reference
32³/64³ differences reach 52.48% in velocity and 80.72% in K. This is not a
converged DNS or stationary/multi-seed campaign, and it does not establish that
2D caused the prior weak benefit.

Supplying exact filtered K and velocity gives production correlation .92355 but
38.12% relative RMS error, 28.18% mean overprediction, 13.33% wrong signs and
50.04% deviatoric-stress error. This is a cell-centered mathematical diagnostic,
not native face production; coupled evolution/numerical/source causes remain open.

Both full comparisons pass six numerical records each; no limiter/fallback is
used in the 3D run. Release/bounds Debug builds and 48 initializer checks pass.
Existing 46-record unforced and 32-record configuration suites pass in both builds,
with 520 old checkpoints byte-identical. Solver physics is unchanged.
[3D report](evidence/2026-10-08-3d-resolution/report-3d.pdf) and
[evidence/reproduction](evidence/2026-10-08-3d-resolution/README.md) preserve results,
review, build/check records and hashes/locations for all raw final state.

## Larger 2D resolution comparisons

Commit `675c7460` raises only the optional comparison grid/runtime caps.
Full 256² and 512² preparations (t=2), half-resolution reference checks and
64² SGS-off/on decay lanes (t=.5) run locally with fixed coefficients and solver.
Both pass six records and exact native initialization; no limiter/fallback is used.
A separately verified eight-thread CPU build retains the solver physics and reproduces
all 244 matched regression checkpoints exactly. Default test limits and numerical
thresholds remain intact. Each report contains 9 pages /41 panels.

| Reference / coarse | Velocity off / on | Vorticity off / on | Mean K recovered | Half/full reference velocity difference |
|---|---:|---:|---:|---:|
| 128² → 32² | 40.251% / 40.300% | 41.347% / 41.658% | 36.10% | 21.14% |
| 256² → 64² | 39.619% / 39.689% | 41.020% / 41.113% | 20.45% | 22.48% |
| 512² → 64² | 50.118% / 50.100% | 50.414% / 50.455% | 7.96% | 20.69% |

Larger grids still show nearly identical SGS-on/off resolved errors and a deficit
of modeled K. Reference convergence is not established. Separate fine-grid preparations
and the old 32² versus new 64² filter width prevent a pure convergence inference;
one seed and decaying 2D flow do not complete the scientific release gates.

[Full evidence and reproduction](evidence/2026-10-08-large-2d/README.md),
[256² report](evidence/2026-10-08-large-2d/256-report.pdf),
[512² report](evidence/2026-10-08-large-2d/512-report.pdf).
Large raw output is archived externally with hashes in the evidence packet.

## First driven 3D paper-test attempt

Local commits `afb8b983` and `ab670494` add signed Fourier bounds and a bounded,
locally runnable Section 4/Table 4 attempt. The published coefficients are fixed.
A 32³ pilot informed the forcing amplitude; the 64³ case ran continuously to 10T
with a native restart at 4T and complete statistics over 2T–10T.

| Quantity | Executed 2T–10T mean ± temporal scatter | Paper 64³ |
|---|---|---|
| Resolved volume RMS Mach | 5.42481 ± 0.196 | 5.38 |
| Volume RMS SGS Mach | 1.00262 ± 0.0365 | 1.107 ± 0.053 |
| Mean K/(rho0 V²) | 0.0317117 ± 0.00211 | 0.0726 ± 0.0055 |
| Volume mean specific epsilon L/V³ | 0.188696 ± 0.0146 | 1.236 ± 0.141* |
| Volume mean rho epsilon L/(rho0 V³) | 0.362124 ± 0.0392 | 1.236 ± 0.141* |

Both 64³ segments pass numerical budgets/admissibility. Applied normalized
production/dissipation are 0.361956/0.362292; work-normalized energy residual
is 7.38e-13. The SGS reservoir is sustained,
while the normalized energy/dissipation do not quantitatively match the paper.
Native forcing normalization/envelope, Enzo-versus-AthenaK differences, and
assumed V/c0=5.5 make this an approximate reproduction. Temporal scatter is not
SEM; eight block means/trends are reported rather than declaring convergence.
Mean thermal energy and sound-speed squared rise about 16.3% over sampling,
contributing to the Mach decline. Fitted K/SGS-Mach trends are −5.64%/−10.05%;
these are regression estimates, not endpoint changes.
The table's epsilon convention is ambiguous; both explicit definitions are retained.

Only the stage-consistent signed-bound guard changes production code. Checkpoint
streaming, independent cell diagnostics and 18-panel plots reuse existing machinery.
Signed-mode 21/21 and existing-forcing 19/19 records pass; a zero-duration 4T restart
has exactly equal fields/OU/RNG/counters. A reporting assumption about the first
continued output was corrected using that probe, preserving the original attempt.

[Paper-test evidence](evidence/2026-10-09-paper-3d/README.md) records full metrics,
plots, costs, raw archive verification, failed/recovered reporting attempts,
exact fingerprints and local reproduction commands. No high-resolution filtered
reference, resolution law, pressure-work heating fraction, MPI/GPU qualification
or M0–M9 completion is claimed. Prior decay evidence remains frozen.

## Fixed-forcing reconstruction diagnostic

Source commit `7cdff361` enables existing monotone PPM4 alongside DC/PLM with
LLF retained. Independent bounded qualification passes 16/16 records, config
checks 34/34 and the independent closure oracle. Fresh matched 32³/64³ turbulence
lanes use the same forcing, seed, coefficients and three ghost zones, and run to
4T with averages over 2T–4T. Completed campaign budgets/admissibility pass.
Original timeouts are retained as failed; verified continuations are recorded
separately.

| 64³ quantity, average 2T–4T | PLM | PPM4 | Relative change |
|---|---|---|---|
| Resolved RMS Mach | 5.49394 | 5.51484 | +0.38% |
| SGS RMS Mach | 1.04634 | 1.15413 | +10.30% |
| Mean instantaneous SGS/resolved Mach ratio | 0.190448 | 0.209299 | +9.90% |
| Mean SGS energy / (rho0 V²) | 0.0321678 | 0.0371305 | +15.43% |
| Density-weighted dissipation L/(rho0 V³) | 0.36822 | 0.45931 | +24.74% |

The stationary smooth-shear oracle isolates numerical kinetic loss: at 64²,
PLM loses 0.05608% versus PPM4 0.001679% over t=1. This supports numerical-method
sensitivity without identifying turbulence's compression-corrected heating fraction.
LLF remains unchanged, so its separate dissipation is unresolved. One seed/short
window does not establish the paper discrepancy's full cause or scientific release.

The forcing audit corrects the previous conversational normalization suggestion:
actual correlation T=L/V fixes V=1. Rescaling V to .661 would change the required
correlation time, so it cannot correct the same physical run. Exact pure-solenoidal
Table 4 forcing amplitude and V/c0 remain unestablished; current forcing was held
fixed.
[Reconstruction evidence](evidence/2026-10-09-scheme-comparison/README.md) includes
plots, independent review, primary-source audit, resource caps, exact commands and
verified raw archive. Earlier evidence is preserved.

## Executed qualification

Serial double CPU builds use the exact branch Kokkos gitlink
`08ceff92bcf3a828844480bc1e6137eb74028517` (Kokkos 4.4.0).
Release and Debug with bounds checking compile successfully.

| Packet | Executed result |
|---|---|
| Independent built-in SGS-off comparison | Six cases match untouched multiphase_turb fields bitwise: wave16/32, two-scalar wave, normal scalar checkpoint restart, Sod128 and modern forced16³. |
| Built-in seed initialization | Two additional checks pass: uniform seed .03 is exact and conserved state stays constant; generic density perturbations give maximum K/rho seed error 3.46945e-18. |
| Raw closure oracle | 1025 cases, 22575 assertions; maximum scaled error 1.45213e-15 against 5e-13. Release and bounds Debug pass. Prior mutation tests are not counted as new port runs. |
| Release unforced matrix | 46/46 result records pass; worst E+K drift 2.94099e-16. Actual stress limiting and fallback are exercised. |
| Release configuration | 32/32 pass, including omitted/disabled equivalence, modern guards, output selectors and history bounds. |
| Mixed longitudinal acoustic branch | Two grids (16/32), four refinement ratios .255–.315; all eight acceptance checks pass. Passive fraction errors remain at roundoff; no fallback. |
| Actual forcing kernel | One independent Heun oracle passes; maximum error 5.09141e-16 at tolerance 1e-12. |
| Release forced/restart matrix | 19/19 records pass: solenoidal/compressive SGS-on/off pilots, three continued lanes and exact restart comparisons, eight contract rejections and a twelve-scalar history case. Fields, modal amplitudes, metadata/update count, RNG and integrated budgets match exactly on restart. |
| Bounds Debug | Full 46-record unforced matrix, both acoustic grids/eight checks, actual forcing oracle and closure pass. Full forced pilot packet was run in Release only. |
| Filtering analysis | 10/10 synthetic tests pass; this verifies analysis machinery, not turbulent closure accuracy. |
| Source checks | Added C++ lines pass style/custom checks; all six Python files parse. Whole changed-file cpplint retains pre-existing diagnostics, so no whole-repository style pass is claimed. |

Decay temporal orders are 1.95053 and 2.00308. The isolated stress full-solver
Richardson orders are about 2.02 and 2.01; that numerical reference is not an exact
mixed manufactured PDE solution. Tested hard collisions activate 9 and 19 fallback
stages, and the high-K stress case suppresses .2126880755333241 transfer. Suppressed
diffusion remains exactly zero in tested cases, so active diffusion limiting is open.

The short 16³ forced pilots reach t=1. Applied work is .09682–.09990 rather than
exactly the nominal .1. The worst measured work/represented-energy residual is
2.77556e-16. Late kinetic fractional drift is about 28% in solenoidal cases and
51–52% in compressive cases: no stationarity, converged spectrum, reference accuracy
or scientific benefit is claimed.

The supported contract remains Newtonian ideal-gas HD, uniform cubic periodic
cells, double precision, RK2, LLF, DC/PLM with at least two ghost zones,
or native monotone PPM4 with three ghost zones, and CFL≤.2.
Cooling, user source/cooling hooks, frame tracking, particles, shearing,
viscosity/conduction, AMR/SMR, MHD and relativistic/multiphysics combinations are
rejected with SGS enabled. Multi-rank MPI and GPU require separate qualification.
The existing multiphase turbulence example selects unsupported SGS combinations;
it was preserved. Use the dedicated SGS inputs instead.

## Reproduction and provenance

Runnable instructions and the model contract are in
`/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu/tst/sgs/README.md`.
Build with `-DPROBLEM=sgs_verify -DCMAKE_BUILD_TYPE=Release
-DAthena_ENABLE_MPI=OFF`, then run CTest and the four verification runners. Use new
output directories. Generic pgen comparisons use a separate built-in Release build.

[The evidence bundle](evidence/2026-10-08-multiphase/README.md) retains source,
inputs, results, build settings, run logs, checkpoint state and hashes. Run manifests
record the base SHA plus the then-uncommitted diff and per-file fingerprints;
[commit provenance](evidence/2026-10-08-multiphase/commit-provenance.json) confirms
final Release production files match the local code commit. Runs are not relabelled
as having started at a future commit. Preformat Release summaries are kept separately.
Debug's full matrix preceded two formatting-only edits; both builds were refreshed
and their closure tests rerun. The final Release packet was rerun after formatting.

Three harness mistakes are preserved in the baseline bundle: the initial restart
comparison used the wrong time, the first Sod setup placed the discontinuity outside
the domain, and the first port seed check mistook a scalar output column for K.
The corrected checks use matched time, a physical Sod state and the named SGS column.
No acceptance threshold was relaxed to turn these into passes.


## Implementation re-audit on 2026 10 09

Three independent read-only equation/integration/diagnostic reviews found no
confirmed defect explaining low SGS energy. The gap between kernel checks and
native production amplitude was closed with a test-only addition, local source
commit `3ca288aa`. Eight one-step 16³/32³ PLM/PPM4 cases pass 24/24 independent
production/dissipation and structural checks. Worst relative rate errors are
7.394e-6 and 1.017e-6, below the frozen 2e-5 tolerance. The kernel oracle again
passes 1025 cases/22575 assertions. Production source and coefficients are unchanged.

LLF numerical resolved kinetic loss can heat the gas without feeding SGS K;
its driven-box contribution remains unmeasured. Near-cutoff gradient attenuation
and incomplete forcing normalization remain numerical/benchmark limitations.
This audit does not convert the short paper comparison into scientific validation.
General mixed production, diffusion and full statistical agreement remain open.
[Audit, exact protocol, results and local reproduction](evidence/2026-10-09-implementation-reaudit/README.md).


## Compression-corrected heating measured on 2026 10 10

Local source commits `1044ddd6`, `e97db55f` and `7bf794ad` add an applied
resolved thermal-pressure-work ledger and minimal reporting/verification.
Closure coefficients, fluid equations and source exchanges are unchanged.
History capacity 33 preserves the forced 12-scalar case. All 23 independent
heating checks pass, along with 34 configuration, 16 PPM qualification and 24
production-amplitude checks. The smooth pressure ledger refinement order is 1.873.

Matched 32³ fresh runs and 64³ native continuations qualify through 4T. The 64³
final fields/forcing/RNG match the preceding scheme runs exactly. Over 2T–4T,

| Case | SGS heating | Numerical residual |
|---|---:|---:|
| 32³ PLM | 47.48% | 52.52% |
| 32³ PPM4 | 63.35% | 36.65% |
| 64³ PLM | 56.64% | 43.36% |
| 64³ PPM4 | 70.05% | 29.95% |

These quantify competing numerical heating after resolved compression is
subtracted. The numerical residual combines spatial/RK/shock-associated effects,
rather than isolating LLF. The short 32³ CFL probe changes fractions by less
than .003 percentage points, but does not bound full-window RK heating.
PPM4 approaches the paper's ~75% SGS reference; full 2T–10T statistics
are not yet established;
one-seed stationarity and forcing equivalence remain unresolved.
[Budgets, plots, independent review, failures and reproduction](evidence/2026-10-10-heating-budget/README.md).


## Full-window matched statistics on 2026 10 11

Both 64³ PLM and PPM4 native continuations now reach 10T, with the applied
pressure ledger and time-weighted 2T–10T statistics. Native equations, closure
coefficients and forcing remain unchanged from `7bf794ad`; local source commit
`808066a4` adds guarded restart-window reporting, larger-grid/seed/forcing CLI
parameters and a basename regression fix. No push or scientific approval.

| 2T–10T quantity | 64³ PLM | 64³ PPM4 | Paper 64³ solenoidal |
|---|---:|---:|---:|
| Resolved volume RMS Mach | 5.42481 | 5.47433 | 5.38 |
| Volume RMS SGS Mach | 1.00262 | 1.10987 | 1.107 |
| Mean K/(rho0 V²) | .0317117 | .0367689 | .0726 |
| Mean rho epsilon L/(rho0 V³) | .362124 | .456623 | 1.236* |
| SGS fraction of compression-corrected heating | 56.3535% | 69.7232% | ~75% |

*The paper's epsilon header differs from the nearby density-weighted prose;
both explicit dissipation definitions and temporal scatter appear in the packet.
PPM4 matches mean SGS Mach closely but has only 50.65% of the published K.
The energy/dissipation gap persists over the full paper averaging interval.
Temporal scatter and eight adjacent 1T blocks are not independent-seed uncertainty.

The operational statistical screen was frozen before outcomes: fitted absolute
drift and early/late difference ≤10%, applied |P−D|/D ≤5%. Both 64³ cases fail:
PLM SGS Mach drift −10.0524%; PPM4 SGS Mach −10.5857% and density-weighted
dissipation −12.3564% (early/late difference 11.1349%). These limits are our
screen, not published acceptance requirements. Full averages are descriptive
measurements; no stationary/ensemble-convergence certificate is claimed.

Both completed recovery segments pass 15/15 native checks; every intermediate
checkpoint was checked. Original interrupted attempts retain FAIL, with separate
verified partial segments and successful native recoveries. Exact inventories
and cumulative counters were checked at restart joins without rebasing.
Full-window forcing/represented-energy residual ≤6.04e−14 and SGS-transfer
residual ≤4.89e−14; low K is not accompanied by an accounting defect in these runs.

Two separate 32³ gamma controls finish 10T and pass 15/15 native checks each.
At fixed initial c0 and forcing, gamma 1.001→1.0001 reduces the mass-weighted cs²
rise over 2T–10T from 16.4468% to 1.6761%, and SGS Mach fitted drift from −8.6424%
to −2.9029%; K trends remain similar. Heating is an established contributor to
Mach drift in this control, not a complete attribution of the 64³ K discrepancy.

Independent forcing algebra/covariance checks pass 14/14; native instantaneous
normalization differs from fixed modal Gaussian OU forcing and its final radial
amplitude peaks at mode 1.5. The original paper Enzo patch/input remains missing.
Gaussian/Favre analytic checks pass 18/18, but no resolved turbulent reference or
paper-identical Gaussian width convention is established. The 128³ startup pilot
passes independent native-file validation; a two-thread 10T startup extrapolation
is at least 7.33h. Full 128/256 scaling, extra seeds and compressive campaigns are
prepared but not run. Cluster/original data and an explicit compute allocation
are needed for adequate reference and scaling work; no remote access was inferred.

[Full results, 10T plots, failures, protocols and reproduction](evidence/2026-10-11-full-window/README.md).


## Native mixed-production gap closed on 2026 10 11

Local source commit `f9b3d596` adds a small native C2/full-mixed production test
and run instructions, using the existing 3D initial-state loader. No native
solver, closure coefficient, forcing, restart format or equation changes.
Forty one-step 16³/32³ PLM/PPM4 runs pass 116/116 specialized amplitude/exchange
checks and 720/720 native invariant checks. Maximum relative rate error is
1.86452e−6 against the pre-existing frozen 2e−5 threshold. Full 3D compression,
expansion, trace-free strain/backscatter and smooth periodic nonzero total
production are covered. Independent rank-one, rotation and density-scaling
reference self-checks pass. An initial periodic-oracle mistake omitted transverse
velocity work; that failed attempt is retained, followed by two complete corrected
runs. The solver did not require a correction.

The actual 64³ 2T–10T stress-limiter increments are 8.1860e−7 (PLM) and
8.4529e−6 (PPM4), only .00002827% and .00023129% of net production.
Both have zero fallback and diffusion limiting. Stress limiting suppresses
negative production and directly increases net production. It cannot explain
the roughly 49% low PPM4 SGS energy. Reconstructed native faces at nine snapshots
show all components present, with C2 the largest positive contribution.
Snapshot decomposition is not a partition of accumulated RK-stage counters.
Its three independent reconstruction checks pass with error ≤1.73e−17.

The turbulent SGS-energy discrepancy remains unresolved. These fixtures close
a native amplitude gap; they do not establish filtered turbulent closure accuracy
or justify retuning coefficients. Controlled forcing/reference tests remain next.
[Protocols, retained failures, local reproduction and eight diagnostic panels](evidence/2026-10-11-mixed-production/README.md).


## CFL sensitivity and cascade diagnostics on 2026 10 12

Matched 64³ PPM4 native restarts at CFL .1/.05/.025 finish exactly 4.5T from
one 4T state. All 45 native checks pass; underlying Fourier OU arrays and RNG
states match exactly. Mean K over 4T–4.5T is .0378136423/.0378132456/.0378130787.
Maximum mean difference versus the finest run is .0014905%; the largest pairwise
endpoint spread is .0060983%. The interval covers about 6.26 global SGS depletion
times, but not a full 2T–10T statistical campaign. Slightly non-monotonic endpoints
prevent a formal Richardson-order claim. The 1% operational screen passes.
The source audit and direct uniform-decay test also reject a missing RK2 factor .5:
full Heun recurrence error ≤1.90e−15; SGS source/history weights are consistent.
Native source/equations/coefficients/forcing remain unchanged; no rebuild or push.

The existing 32³→64³ PPM4 mean K decreases .059068→.036769, giving two-point
exponent .68389, near Kolmogorov 2/3. This is only relative scaling, with 32³ below
paper floor and one seed. It does not validate absolute K. Conditional classical
sharp-filter estimates using our flux/heat proxy are larger, but assume homogeneous
incompressible turbulence and a filter convention not certified for this solver.
The SGS-dissipation ratio^(2/3) tracks the energy ratio to the paper; since
SGS dissipation is derived from K, this is not independent validation or cause.

Nine-snapshot resolved spectra lack a broad flat compensated K41 range. Frozen
shells 4–10 give 64³ slopes−2.746 (u),−2.180 (rho^(1/3)u), with Parseval error
≤4.34e−16. These short-range fits show limited scale separation and possible
spatial damping, not a uniquely identified mechanism. The roughly 49% SGS energy
shortfall remains open. Same-forcing high-resolution filtered reference and
forcing/resolution/seed qualification are still required.
[Protocols, bounded CFL test, derivations, 14 diagnostic panels and retained audit attempts](evidence/2026-10-12-cfl-kolmogorov/README.md).

## Direct spectrum-tail estimate on 2026 10 12

The operator clarified that the desired Kolmogorov estimate must fit the
measured resolved spectrum and integrate its missing tail, rather than infer
amplitude from dissipation. For E(n)=A n^(-5/3), fitting log amplitude over
shells 4–10 and integrating above spherical Nyquist nc=N/2 gives
(3/2)A nc^(-2/3). No dissipation or SGS energy sets the fitted amplitude.
Existing velocity shell sums already contain the factor 1/2; no second half
or extra 2pi is applied. The rho^(1/3)u spectrum is not treated as energy.

The 64³ primary estimate is .0353572 versus .0374237 SGS averaged over the
same nine 2T–10T snapshots (.0367689 for the full time-weighted average).
Fit-window sensitivity is .0294707–.0491582 at nc=32. These estimates encompass
the evolved SGS energy and do not independently show a factor-of-two deficit.
At 32³, the primary estimate is .0299808 versus matched SGS .0592459;
its steeper, more damped spectrum makes this extrapolation especially fragile.

There is no clean inertial range. Cutoff N/4 raises each estimate by 1.5874,
making the 64³ primary value .0561261. The spherical sharp filter is assumed,
not established for this finite-volume method. Multiplying the velocity-specific
tail by rho0 approximates density energy; it does not recover exact compressible
Favre SGS energy. The paper energy gap remains open. No native source change,
new simulation or coefficient adjustment was made. Eighteen checkpoint hashes,
times, mean masses and unit domains were checked; synthetic power-law and
independent integral self-checks pass.
[Derivation, matching averages, fit/cutoff sensitivity and four-panel plot](evidence/2026-10-12-spectrum-tail/README.md).

## Native filtered-snapshot production test on 2026 10 12

Source commit `b513bc9d` adds a Gaussian extension to the existing Favre reference
helper, two independent filter tests, a host utility calling the actual SGS
stress kernel, a snapshot analysis runner and local usage instructions. Native
solver equations, coefficients, forcing, timestep/source splitting and restart
format are unchanged; the Athena executable hash is unchanged. No new dependency
or push. The analysis follows the paper's Section 3 comparison on identical
supplied flows, with measured covariance K supplied to the closure.

Fresh SGS-off 32³/64³ RK2/LLF/PPM4 runs with gamma 1.001, CFL .1, seed 42 and fixed
acceleration RMS 3 finish 10T. Nine nominal integer-T snapshots over 2T–10T are
filtered at box/16, box/8 and box/4, using Gaussian and box filters. Gaussian
width convention Delta_G=sqrt(12)*sigma is explicit, not certified paper-identical.
The primary 64³ Gaussian box/16 comparison has mean transfer ratio .929245
(7.0755% low), mean signed correlation .994893 and mean relative local RMS error
.130715. Mean trace-free production correlation .989717 confirms that the shared
isotropic pressure term does not explain all the agreement. Spectral derivatives
give mean ratio .939719; centered differencing does not produce a factor-two bias.

All 1134 comparison arithmetic checks, 590 independent inventory checks over 82
native checkpoints, 22 continuation checks, 12 filter unit tests and the existing
native closure oracle pass. Native/Python stress disagreement ≤5.82e−16; direct
positive centered covariance, mass/momentum and kinetic-energy decomposition,
covariance positivity, trace, symmetry and Galilean invariance are checked.
Initial generic report attempts remain FAIL due to a basename/history mismatch;
native simulations completed successfully and independent corrected validation
reads the actual saved history names. Continuations use the forcing-aware harness.

This supports production closure on these flows and does not establish evolved
SGS-energy accuracy. Primary K_reference=.0824366 is measured at filter box/16,
not grid box/64, so it cannot be compared directly to the 64³ evolved K=.0367689.
Mean reference K differs 27.67% between 32³ and 64³ at that physical width, exceeding
the earlier proposed 10% sensitivity screen. The32³ Gaussian is undersampled there.
No converged reference or paper-sized scale separation is claimed. A same-scale,
adequately resolved comparison with evolved SGS energy remains necessary.
[Fixed protocols, retained failures, all metrics, reproduction and 14 diagnostic panels](evidence/2026-10-12-filtered-snapshot/README.md).
