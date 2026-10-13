# Higher reference for unforced 3D decay — 2026-10-12

The matched 256³/128³ campaign completes10T: six native runs,488/488 numerical
checks and126 double field snapshots pass. The SGS energy deficit persists
against the finer reference. The8³ filter meets the predeclared10% reference
sensitivity screen; the16³ filter still fails it. This is scoped evidence for
this local decay benchmark, not universal model qualification or an exact paper
reproduction. No model equations or coefficients were changed.

## Matched protocol and execution

- Feature branch `codex/multiphase-schmidt-federrath-hd`, worktree
  `/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu`.
  Local test-only commit `921cf83279ff61b08f9f9f75c487b20d90c565a8`; no push.
- Periodic unit cube, rho1, initial sound speed1 and RMS Mach9; seed42, Gaussian
  modes1–4 with2:1 solenoidal/longitudinal energy. Initial density is uniform.
  This is an undeveloped local field, not the paper's imported Kritsuk data.
- Gamma1.001 near-isothermal ideal EOS; CPU double, RK2/LLF/PPM4, CFL.2.
  C1=.02, C2=.7, Ce1.5 and Ckappa.65 unchanged. No forcing throughout.
- Independently evolve256³ and128³ SGS-off references plus8³/16³ SGS-on/off.
  All states are conservative/Favre restrictions of the same256³ initial field.
  Removed variance is initialized into SGS K on and heat off. Initial K equality
  is imposed; score only1–10T. T=1/18, defined from initial velocity RMS.
- Native fields every nominal.5T, histories every.1T, checkpoints every1T;
  actual field times are retained and interpolated to a common reporting grid.
  Reference sensitivity is max |K256-K128|/max(K256,K128) over1–10T.
- Execution uses8 OpenMP threads, maximum64-cell block axes and zero ordinary
  passive scalars. The separate SGS reservoir remains evolved. These resource
  choices were tested before the full run; they do not change the SGS method.

[Protocol](protocol.json) and actual per-case inputs/checks in [results](results.json)
define this decay experiment. Some generic manufactured-suite defaults remain
in the shared [manifest](manifest.json); they are not the decay grid/cycle/runtime
limits. The decay protocol specifies256³,50000 cycles and7200s per native case.

The256³ native loop completes4728 cycles in5975s (99.6minutes); including field
verification/filtering it takes6062s (101minutes). The fresh128³ solver plus
verification takes451s (7.5minutes). Total case time is108.6minutes, excluding
initial preparation and final plotting. These are single-run laptop timings.

## Measured result

All ratios below refer to the fresh256³ reference, filtered at the coarse scale.

| Quantity |8³ SGS|16³ SGS|
|---|---:|---:|
| Mean model/reference K ratio,1–10T |.69840|.66505|
| Model/reference K ratio at10T |.41885|.34980|
| Maximum128³/256³ reference sensitivity,1–10T |9.67%|21.52%|
| Frozen10% reference screen |Pass|Fail|
| Mean reference K,1–10T |.79160|.51224|
| Mean evolved model K,1–10T |.55285|.34066|
| Model/reference K at10T against128³ |.45009|.43633|
| Resolved+SGS kinetic/reference, integrated1–10T |.54122|.75001|
| SGS-off kinetic/reference, integrated1–10T |.42997|.66555|
| Numerical/RK fraction of compression-corrected heating |12.84%|29.44%|

Independent integration on native time knots gives K ratios .69643/.66378,
within.3% of the declared reporting-grid values. See [additional metrics](additional-metrics.json)
and the [independent audit](independent-audit.py). Integrated kinetic improvement
does not imply that every time or spectral scale improves with SGS.

The finer reference has more measured unresolved energy, while evolved coarse
K changes very little. The earlier128³/64³ campaign reported mean ratios .75/.79;
the new ratios .70/.67 worsen agreement. The old coarse trajectories were not
silently reused: the new256³ field changes initial coarse K by.21%/1.09%, so
all six cases were rerun. [Initial comparison](higher-reference-initial-comparison.json).

At the8³ filter, the model's mean deficit is30.2% and its final deficit58.1%,
while this pair's reference sensitivity remains below10%. This strengthens the
case that the measured deficit is not just the previous reference-resolution
uncertainty. Passing this screen is not a proof of DNS convergence. At16³,
reference uncertainty still prevents precise absolute accuracy qualification,
although the model ends below both references.

The fine flow ends at resolved Mach.89576 and RMS sound speed1.01989. The late
underprediction is in a subsonic flow. Small coarse grids, one seed, startup
transients, approximate isothermality and the undeveloped initial field remain
limitations. The8³ grid poorly resolves the original energy-containing modes,
with mode4 at its Nyquist scale. This does not establish a universal missing
one-half factor, its cause, or justify coefficient retuning. F14/T21 and M9
scientific qualification remain open.

## Numerical checks and resource choices

[Commit and source verification](commits-and-verification.json): all324 native
source hashes and the executable match the run manifest; the test harness hash
also matches. All488 checks pass: exact initial readback, admissible reservoirs,
mass/momentum/E+K conservation, Favre kinetic decomposition, completion, final
history inventories and paired SGS transfer ledger. Maximum binary represented
energy drift is4.37e-16 scaled; online history drift is3.92e-14 relative.
[Additional verification](additional-verification.json).

Both SGS lanes have zero fallback stages and zero limited diffusion integral.
Final stress-limited transfer integrals are.00433/.00171. The K change matches
applied production minus dissipation to rounding accuracy. These checks verify
numerical bookkeeping; they do not validate the closure's turbulent K amplitude.
The numerical/RK heating percentage is not a missing-K percentage or a causal
attribution.

Before the expensive run, changing blocks16→32 and threads2→8 gave four bitwise
identical physical snapshots in a64³ pilot. Changing blocks32→64 and ordinary
passive scalars2→0 gave24 bitwise identical physical/SGS snapshots across six
pilot lanes. [Block equivalence](block16-equivalence/results.json),
[physical equivalence](block64-no-passives-pilot/physical-equivalence.json).
The selected32-cycle256³ resource pilot used6.42GiB maximum RSS, fitting24GiB
RAM; the earlier block32/two-passive pilot used8.89GiB. These cycle-capped pilots
are resource measurements, not completed scientific decay tests.
[Resource protocol](higher-reference-resource-protocol.json).

One resource-pilot invocation failed before initialization because a command-line
override tried to add an absent output precision parameter. The corrected saved
input succeeds; [failed attempt](resource-pilot-256-lean/failed-attempt.json)
and both logs are retained. An optional final summary query initially used wrong
check field names and was corrected; native results were unaffected. The full
six-case campaign has no failed numerical checks or failed native runs.

## Plots and reproduction

[Five-page PDF](report.pdf),24 panels, visually inspected:
[energy and ratios](energy-comparison.png), [decay/conservation](decay-and-conservation.png),
[source/heating budgets](source-budgets.png), [local statistics and spectra](local-statistics.png),
[final slices](final-slices.png). Slice colors use common ranges; chaotic fields
are compared statistically, not through pointwise trajectory agreement.

[Summary](summary.json), [150 raw hashes](raw-manifest.json),
[completed state](higher-reference-state.json). Raw field/checkpoint/initial files
remain outside wiki Git at:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-unforced-decay/decay-256`.
Archived [runner](higher-reference-source/decay_turbulence.py) and
[reproduction notes](higher-reference-source/README.md) are inspection copies;
execute from the feature checkout with the existing NumPy/Matplotlib environment
and CPU-double `-DPROBLEM=sgs_verify` binary, using a new output directory:

```sh
OMP_NUM_THREADS=8 OMP_PROC_BIND=false python3 -B tst/sgs/decay_turbulence.py \
  --exe /absolute/build/src/athena --output /tmp/sgs-decay-256-new \
  --high 256 --low 8 16 --block 64 --nscalars 0 --end-T 10 --timeout 7200
python3 -B tst/sgs/decay_turbulence.py --output /tmp/sgs-decay-256-new --plot-only
```

The next diagnostic is the production/dissipation balance of the filtered
reference against the evolved SGS across Mach regimes. A finer reference would
still be needed to qualify the16³ filter. Do not turn these scoped results into
full scientific approval or retune before identifying the discrepancy.
