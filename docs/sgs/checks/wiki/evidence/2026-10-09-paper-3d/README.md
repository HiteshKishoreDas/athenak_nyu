# First driven 3D paper-test attempt

Code branch `codex/multiphase-schmidt-federrath-hd`, based on
`origin/multiphase_turb` at `526a187f938e38634e89a88192ebe9731c4d36bb`.
Local commits: `afb8b983` permits signed Fourier bounds;
`ab670494c9d097e59ac1210cbc812339a7954494` adds the bounded runner, streaming
checks, statistics, plots and run instructions. No push or merge.

This attempts the solenoidal 64³ case in
[Schmidt & Federrath, Section 4/Table 4](https://arxiv.org/html/1010.4492).
The box side is 1, L=.5, V=1, T=.5, gamma=1.001, rho0=1, initial u=0,
c0=1/5.5, pressure0=c0²/gamma, k0=1e-8. Coefficients are fixed at
C1=.02, C2=.7, Ceps=1.5, Ckappa=.65. Continuous OU forcing has
correlation T, updates every .01T, pure solenoidal projection, signed modes
−3…3 in each component, radial band 1…3 and peak 2. RK2/LLF/PLM,
CFL=.1, one MPI rank, eight OpenMP threads, block 16³ and double precision.
Snapshots/histories are spaced .1T. Source and executable hashes are preserved.

A 32³, 4T pilot with acceleration RMS 2 is below the paper resolution floor.
A 64³, 4T pilot with RMS 3 produced resolved Mach 5.494 over 2T–4T, so RMS 3
was held fixed for the continuation to 10T. The forcing amplitude alone was
calibrated; closure coefficients were never tuned. Both segments carry absolute
forcing and SGS counters. The native restart at 4T truncates a timestep and is
not asserted bitwise identical to an unrun fresh 10T realization. The full
2T–10T interval is sampled. The original source runner supports a fresh 10T run.

| Quantity | Executed 2T–10T mean ± temporal scatter | Paper 64³ |
|---|---|---|
| Resolved volume RMS Mach | 5.42481 ± 0.196 | 5.38 |
| Volume RMS SGS Mach | 1.00262 ± 0.0365 | 1.107 ± 0.053 |
| Mean K/(rho0 V²) | 0.0317117 ± 0.00211 | 0.0726 ± 0.0055 |
| Volume mean specific epsilon L/V³ | 0.188696 ± 0.0146 | 1.236 ± 0.141* |
| Volume mean rho epsilon L/(rho0 V³) | 0.362124 ± 0.0392 | 1.236 ± 0.141* |

The displayed deviations are temporal scatter of the stored linear time
interpolant, not uncertainty of the mean. Eight 1T block means, early/late means,
slopes and explicit averaging definitions are in `paper64-statistics.json`.
*The table heading labels specific epsilon while adjacent prose discusses
rho epsilon; both definitions are retained without choosing a favorable one.

Applied normalized production/dissipation means are 0.36195564/0.36229247;
relative imbalance is 0.00092972.
Maximum represented-energy residual divided by accumulated actual forcing work
is 7.381e-13. Numerical source, inventory,
scalar, momentum, history and admissibility checks pass for both segments.
Fallback stages and diffusion suppression are zero. Accumulated stress suppression
is tiny but nonzero (1.48e-6), too small to explain the energy deficit. Exact
pressure-work is not recorded,
so the paper's compression-corrected heating fraction is unqualified.

The SGS reservoir is sustained under driving. Quantitative energy/dissipation
agreement remains incomplete under the declared native forcing units. This case
alone cannot distinguish the effects of forcing normalization, solver dissipation,
SGS discretization and reference convergence. Full-window sampling is not a
stationarity certificate; the measured trends and blocks remain descriptive.
Over 2T–10T, mean thermal energy and volume-mean sound-speed squared rise
16.33% and 16.27%, respectively. Heating therefore contributes to decreasing Mach.
Linear fitted fractional trends are −5.64% for K and −10.05% for SGS Mach;
these are fitted trends, not endpoint changes, and do not isolate heating from
reservoir evolution.
No cutoff-scaling fit, filtered 1024³ comparison, seed ensemble or full scientific
release qualification is claimed.

## Departures and independent checks

- Native instantaneous volume RMS acceleration normalization replaces published
  stochastic modal amplitude. V/c0=5.5 is an explicit assumption because the
  solenoidal Table 4 row does not state it. Thus comparisons in V units are conditional.
- The shell-compensated native parabolic envelope is not proven identical to the
  [forcing reference](https://arxiv.org/pdf/0905.1060). AthenaK differs from Enzo 1.5.
- The small positive SGS seed avoids the current discrete zero-K absorbing state.
- Existing positive-octant defaults are preserved. New signed bounds remove the
  missing mixed-sign directions for this input. Analytic normalized covariance
  error is 5.55e-17; native solenoidal projection residual is 1.86e-15. The 21-record
  signed/restart/bound-rejection suite and 19-record existing forcing suite pass.
- Streamed/retained checkpoint lanes have exact matching fields and OU state.
  Cellwise Mach/dissipation and irregular-time averaging selfchecks pass, as does
  the independent C++ closure oracle. Frozen numerical tolerances were retained.
- The short same-checkpoint timing pilot measured 8 threads: 9.70845 s versus
  4 threads: 12.0788 s; final fields differed by 3.36e-15. These are single short timings.
- The first continued checkpoint is at 4.1T. A zero-duration restart of 4T has
  exactly equal active fields, OU/RNG and counters, bound to the checkpoint hash.

## Failed and recovered reporting attempts

The first signed rejection harness omitted max_k keys; its failed diagnostics
remain in `signed-first-results.json`. The corrected suite passes 21/21. The 32³
simulation passed, then its first plotting invocation failed because the report
module was still being written; plot-only recovery succeeded.

A reporting assumption compared the first advanced 4.1T checkpoint with 4T as a
load check. This was discovered while the simulation continued. Process enumeration
was unavailable in the tool environment, so no signal was sent to the obsolete Python reporting wrapper. Its invalid assertion failed
after the simulator and native verification passed. Corrected postprocessing
rereads all raw states, checks cadence/ordering and counter continuity, and uses
the exact zero-duration load probe. The original script, recovery reason and corrected
analysis are preserved. No simulator rerun or numeric-threshold relaxation was used.

## Plots, cost and reproduction

The main five-page report has 18 panels: Mach/K/dissipation versus published
reference bands, applied transfers, thermal drift, energy residual, guards,
eight 1T blocks, component covariance, density/Mach tails, SGS inventory and
final density/SGS-Mach distributions. Caption qualification distinguishes
numerical PASS, incomplete pilots and the full sampling window.

![Reference statistics](paper64-page1.png)

The orange curves/bands are the published reference; blue curves are this run;
the purple dissipation curve includes density weighting. The lower panels remain
below the reference in the assumed V units. Mach values provide a direct
comparison under the stated thermal sound-speed convention.

![Applied energy transfers](paper64-page2.png)

Production/dissipation use applied source counters. A small E+K residual verifies
accounting, not closure accuracy. Heating changes shown here are uncorrected for
resolved pressure work.

[Full 18-panel PDF](paper64-report.pdf), [32³ pilot](pilot32-report.pdf),
[64³ pilot](pilot64-report.pdf), [statistics](paper64-statistics.json),
[reproduce script](reproduce.sh), [saved input](paper-input.athinput).

The 64³ pilot took 1031.10 s for
simulation/native verification; continuation took
1939.68 s, and corrected
postprocessing took 91.15 s. These costs exclude
some earlier plotting/review intervals and later archive verification; they are
not a controlled performance benchmark. Simulation caps were 1200 s / 3600 s and
100000 total cycles. Fresh-run example cap is 7200 s. All inputs, histories,
checkpoints, failed attempts, exact commands and hashes are retained.

Raw archive (outside wiki Git): `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-paper-3d-raw-evidence.tar.gz`.
Size 3438044616 bytes; SHA-256 `23298f7efa10437960b4a1bc96a23e5eabd586b523181171cc456ebbcd4dcf51`.
All 381 regular members independently match source sizes
and hashes; original bytes 3670417698. The manifest and readback result
are included here. Native raw data and the executable remain in the persistent
visualization folder as well. No large binary simulation data were added to Git.
