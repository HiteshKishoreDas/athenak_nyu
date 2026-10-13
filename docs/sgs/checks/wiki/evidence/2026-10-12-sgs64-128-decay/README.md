# 256³ reference against64³/128³ SGS — 2026-10-12

The requested64³/128³ SGS-on/off comparison completes10T. All six case analyses
and488/488 numerical checks pass. Mean SGS K relative to the filtered256³ reference
is1.02266 at64³ and1.93662 at128³; final ratios are.45964/.83733. Total kinetic
decay is similar with SGS on and off. Reference convergence remains inadequate
for an absolute SGS-accuracy verdict. No solver equations or coefficients changed.

## Setup and reuse

The experiment reuses the completed matched256³ and128³ SGS-off references from
[the higher-reference campaign](../2026-10-12-higher-reference/README.md), then runs
new64³ SGS-off,64³ SGS-on and128³ SGS-on trajectories. The existing128³ reference
also serves as the128³ SGS-off control: six reported lanes represent five unique
native simulations, with three new solver executions. Reused outputs are hashed,
revalidated and re-filtered; they are not relabelled new simulations.

- Feature worktree `/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu`,
  branch `codex/multiphase-schmidt-federrath-hd`, local test-only commit
  `ab54c37d`. Adds guarded reference reuse and unavailable-screen reporting,
  plus runnable reuse checks and reproduction notes. Native source unchanged.
- Periodic unit cube, rho1, sound speed1, initial RMS Mach9; seed42,
  Gaussian modes1–4 with2:1 solenoidal/longitudinal energy. All coarse states
  derive from the same256³ initial field. Removed kinetic variance initializes
  SGS K on and heat off, preserving equal initial E+K. Initial K equality is imposed.
- Gamma1.001 near-isothermal ideal EOS, RK2/LLF/PPM4, CFL.2, double CPU,
  eight OpenMP threads, maximum64-cell block axes, zero ordinary passive scalars.
  The separate SGS reservoir remains active. C1=.02, C2=.7, Ce1.5, Ckappa.65.
- No forcing. T=1/18 from initial velocity RMS; all trajectories reach10T.
  Score1–10T. Fields every nominal.5T, histories every.1T, checkpoints every1T;
  actual timestamps are retained and interpolated to the common reporting grid.
- This remains a local undeveloped initial field, not the paper's imported
  Kritsuk snapshot or an exact paper reproduction. One seed and near-isothermal
  EOS remain limitations. The256³ reference ends subsonic, at Mach.89576.

[Protocol](protocol.json), [manifest](manifest.json), [case results](results.json),
[summary](summary.json), [source verification](commits-and-verification.json).
The shared manifest includes manufactured-suite defaults; actual decay limits
and checks are in the protocol/results, including50000 cycles and7200s per case.

## Energy comparison

K is volume-mean SGS energy density. Ratios use the filtered256³ numerical
reference at each SGS grid's physical filter scale. Means integrate1–10T.

| Quantity |64³ SGS|128³ SGS|
|---|---:|---:|
| Mean reference K |.126653|.038808|
| Mean evolved model K |.129523|.075156|
| Model/reference mean K ratio |1.02266|1.93662|
| Model/reference K ratio at10T |.45964|.83733|
| Final reference K |.033779|.009019|
| Final model K |.015526|.007552|
| Maximum128³/256³ reference sensitivity,1–10T |59.06%|Unavailable|
| Frozen10% convergence screen |Fail|Unavailable|

Native-time-knot integration independently gives mean K ratios1.02224/1.93689,
within.05% of the declared reporting-grid results. At64³, K is initially above
the256³ reference and crosses below it near3T; at128³ it remains above until
roughly8T. A mean near unity does not imply agreement throughout the decay.
At10T the64³ model is still1.11712 times the filtered128³ reference, although
only.45964 times the filtered256³ reference. It therefore is not below both
available references at that filter.

Reference qualification is the central limit. The256³ reference supplies only
four/two fine cells across each64³/128³ coarse cell. At the64³ filter, the two
reference K values differ by up to59.06%, failing the predeclared10% screen.
At128³, filtering the128³ control to its own mesh gives zero removed variance;
that is not an independent measurement of unresolved energy below128³. Its
reference screen is explicitly unavailable, and the zero curve is omitted.
A finer reference pair is needed for absolute SGS-energy accuracy qualification.

The new results do not support a universal factor-of-two energy deficit:
ratios vary with both time and grid size. They also do not prove the128³ model
truly overproduces unresolved energy, given the unqualified256³ reference.
F14/T21 and M9 scientific qualification remain open; no retuning is justified
by this packet alone.

## Resolved and total kinetic performance

Total means resolved+SGS kinetic energy, compared with the full256³ kinetic
energy. Resolved-only ratios compare with the256³ state filtered to the same
coarse grid. See [four performance panels](performance-ratios.png).

| Kinetic-energy comparison |64³|128³|
|---|---:|---:|
| Total kinetic, SGS on /reference |.95410|.98233|
| Total kinetic, SGS off /reference |.95532|.98176|
| Resolved kinetic, SGS on /filtered reference |.94735|.95537|
| Resolved kinetic, SGS off /filtered reference |1.04942|1.00950|
| Final resolved kinetic, SGS on /filtered reference |1.00627|1.01748|
| Final resolved kinetic, SGS off /filtered reference |1.03533|1.02903|

The combined kinetic decay barely differs between on/off. The resolved/unresolved
partition differs more strongly. SGS-on is closer in final resolved kinetic
energy, while SGS-off is closer in the128³ time-integrated resolved energy.
This is mixed performance, not uniform improvement at every time or scale.
Final spectra and density-conditioned SGS velocity are available in
[local statistics](local-statistics.png). Chaotic slices are illustrative,
not pointwise trajectory-accuracy tests.

## Numerical verification and costs

All324 native source hashes, the executable and run harness match provenance.
There are126 case snapshot inventories but105 unique native field snapshots,
because the128³ reference/control is shared. [125 unique raw hashes](raw-manifest.json).
Checks include exact initial readback, admissibility, conservation, Favre kinetic
decomposition, completion, independent binary/history inventories and paired
SGS transfers. Maximum binary represented-energy drift is4.37e-16 scaled;
maximum online history drift across all lanes is3.92e-14 relative.
Both SGS runs have zero fallback stages and zero limited diffusion integral.
Final stress-limited transfer integrals are.00112046/.00045891.

The source ledger closes to rounding accuracy. Numerical/RK heating after thermal
pressure-work subtraction is34.01%/34.52% of the corrected heat gain over1–10T.
Those are heating fractions, not missing-SGS-energy fractions or proof of cause.
[Source and heating budgets](source-budgets.png).

New native loops take42s for64³ off,448s (7.5minutes) for64³ on and5400s
(90minutes) for128³ on. Reference reanalysis takes102s/13s, and shared128³-control
reanalysis4s. [Additional verification and timings](additional-verification.json).
These are single-run laptop costs, not repeated performance benchmarks.

The small native reuse pilot passes six case analyses, exercises three reused
lanes and the unavailable same-grid reference screen, and deliberately rejects
a mismatched seed. [Pilot summary](reuse-pilot/summary.json),
[guard validation](reuse-pilot/reuse-validation.json). The initial linter call
could not write its cache in the restricted worktree; the no-cache run passes.
Self-checks and diff checks pass. Full campaign has no failed native case or
numerical check.

## Plots and reproduction

[Five-page report](report.pdf),24 panels:
[energy](energy-comparison.png), [decay/conservation](decay-and-conservation.png),
[budgets](source-budgets.png), [statistics/spectra](local-statistics.png),
[final slices](final-slices.png). Four additional normalized-energy panels:
[PNG](performance-ratios.png), [PDF](performance-ratios.pdf). All plots visually inspected.
[Additional metrics](additional-metrics.json), [independent audit](independent-audit.py),
[completed state](state.json), archived [runner](source/decay_turbulence.py) and
[source reproduction notes](source/README.md).

Raw new trajectories remain outside wiki Git at:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-unforced-decay/decay-sgs64-128`.
Reused raw references remain in the adjacent`decay-256` directory.
Execute from the feature checkout with the existing NumPy/Matplotlib environment
and double CPU `-DPROBLEM=sgs_verify` binary; output must be new, and the matching
reference campaign must already exist:

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build/src/athena --output /tmp/sgs-decay64-128-new \
  --reference-output /tmp/sgs-decay256-new --high 256 --low 64 128 \
  --block 64 --nscalars 0 --end-T 10 --timeout 7200
python3 -B tst/sgs/decay_turbulence.py --output /tmp/sgs-decay64-128-new --plot-only
```

The same64³/128³ experiment can run without reference reuse by omitting
`--reference-output`, at the cost of rerunning both references. Reuse requires
matching physical and execution parameters and intact source/reference files.
