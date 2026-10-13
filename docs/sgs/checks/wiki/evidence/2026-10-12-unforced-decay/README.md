# Initially supersonic, unforced 3D decay — 2026-10-12

At late times, the SGS energy falls below both available filtered numerical
references without ongoing forcing. All six native runs and 488 numerical checks pass. Both reference
convergence screens fail, so the absolute model accuracy remains unqualified.
This is a new local benchmark following the comparison method of
[Semenov2025](https://arxiv.org/pdf/2410.23339v3), not an exact paper reproduction.

## Setup and provenance

- Feature worktree: `/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu`,
  branch `codex/multiphase-schmidt-federrath-hd`.
- Local source commit: `7250a4e4ce3ec30bd09db460b3c777f6dc8cd2a1`; no push.
  Adds only the decay runner and reproduction documentation. No solver,
  forcing, coefficient, timestep algorithm or reservoir equation changes.
- Unit periodic cube, mean density1, initial sound speed1 and RMS Mach9;
  fixed seed42 random Fourier modes1–4 with 2/3 solenoidal energy.
  Initial density is uniform; this flow is undeveloped, unlike the paper's
  imported Kritsuk snapshot. The first1T is treated as startup.
- Ideal EOS gamma1.001 approximates isothermality; RK2/LLF/PPM4, CFL.2,
  CPU double, one rank/two OpenMP threads. C1=.02, C2=.7, Ce=1.5, Ckappa=.65.
- Independently evolve128³/64³ SGS-off references and8³/16³ SGS-on/off lanes.
  Every coarse state is a conservative/Favre restriction of the same fine IC.
  Removed kinetic variance becomes SGS K on, or heat off. Initial K agreement
  is imposed. All lanes begin with the same represented E+K.
- No `turb_driving` block. T=1/18=box/(2*initial velocity RMS); all runs reach10T.
  Fields every nominal.5T; histories every.1T; actual timestamps are retained.

The run began before the test commit. Its original harness hash remains in
`manifest.json`; final report code hash is in `summary.json`. Subsequent report,
guard and documentation improvements did not alter native inputs or simulations.
All324 native source-file hashes match the final checkout. See
[commit verification](commits-and-verification.json).

## Measured outcome

Energy ratios below use the filtered128³ reference; the convergence screen
compares64³ with128³ at the same physical filter.

| Quantity |8³ SGS|16³ SGS|
|---|---:|---:|
| SGS/reference mean K ratio over1–10T |.74912|.78753|
| SGS/reference K ratio at10T |.45006|.43702|
| Maximum64³/128³ reference sensitivity over1–10T |17.52%|34.61%|
| Predeclared10% reference screen |Fail|Fail|
| Mean reference K over1–10T |.73789|.43240|
| Mean model K over1–10T |.55277|.34053|
| Resolved+SGS kinetic energy / reference, integrated1–10T |.55147|.76396|
| SGS-off kinetic energy / reference, integrated1–10T |.43802|.67801|
| Numerical/RK heating fraction after pressure-work subtraction |12.84%|29.46%|

Means use interpolation to the declared common nominal.5T reporting grid and
trapezoid integration. Independent integration on native timestamp knots yields
K ratios .74720/.78652, differing by less than.3%; headline values are .75/.79.
The128³ reference ends at resolved Mach.89886 and RMS sound speed1.01990.
Late underprediction therefore occurs in a flow that has become subsonic.
The numerical heating percentages are fractions of compression-corrected
thermal gain, not fractions of missing SGS energy or proof of its cause.

The low-energy tendency is present without ongoing forcing. However, failed
reference convergence, a single seed, small coarse grids, initial transients
and late Mach evolution prevent universal model-accuracy claims or a direct
comparison to the original forced-box paper energy. This does not establish a
missing arithmetic factor or justify retuning. The8³ scale has16 reference cells
per coarse cell but poorly resolves energy-containing modes;16³ has only8.

## Checks, plots and failures

- [Native results](results.json): six passes,488/488 checks,126 double snapshots.
  Checks cover exact initial readback, admissibility, mass/momentum/E+K,
  kinetic decomposition, history inventories, completion and SGS source ledger.
- [Summary](summary.json), [additional metrics](additional-metrics.json),
  [protocol](protocol.json), [manifest](manifest.json), [raw hashes](raw-manifest.json).
- [Independent reader checks](reader-checks/results.json): six nonempty initial
  and final comparisons match the independent restart reader bitwise.
- [Accepted pilot](pilot-32-verified/results.json): six native passes/136 checks.
  [Short-window regression](short-window-regression/summary.json): six passes and
  finite integrals for1.1T; the self-check also covers1.25T.
- [Five-page report](report.pdf),24 panels:
  [energy](energy-comparison.png), [decay/conservation](decay-and-conservation.png),
  [source/heating budgets](source-budgets.png), [statistics](local-statistics.png),
  [slices](final-slices.png). Shared color ranges apply across slice panels.
- Three earlier reporting/reader attempts remain failed in the original output
  tree; see [failure ledger](attempt-failures.json). They are not relabelled pass.

Independent review confirmed initialization, filter/reader conventions, budgets,
ratios and interpretation. Python lint/self-checks pass; the existing native
build is up to date. All raw states remain outside the wiki Git tree at:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-unforced-decay`.
The full reference takes about21minutes locally; all six solver/verification
lanes together take about23minutes. These are single-run costs, not a benchmark.

## Reproduce locally

From the feature checkout, use the existing NumPy/Matplotlib environment and
double CPU `-DPROBLEM=sgs_verify` binary. Output directory must be new:

```sh
python3 -B tst/sgs/decay_turbulence.py --self-test
OMP_NUM_THREADS=2 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build/src/athena --output /tmp/sgs-decay-new \
  --high 128 --low 8 16 --end-T 10 --timeout 3600
python3 -B tst/sgs/decay_turbulence.py --output /tmp/sgs-decay-new --plot-only
```

Archived source copies are for inspection; execute from the feature checkout.
The next qualification requires a finer reference within a justified resource
allocation, followed by production/dissipation diagnostics across Mach regimes.
F14/T21 scientific qualification and M9 remain open.
