# Full-window SGS validation, 2026-10-11

**Both 64³ PLM and PPM4 complete 10T, with time-weighted statistics and applied heating budgets over 2T–10T. Numerical checks pass; quantitative paper reproduction remains unestablished.** This extends the matched pressure-ledger runs, not the first-ever PLM 10T run. The recovered PLM full-window means reproduce the earlier paper-attempt means to displayed precision.

Source worktree: `/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu`, branch `codex/multiphase-schmidt-federrath-hd`, local source commit `808066a4ab5baee49c012e80acc2a48ff6134b84`. This commit adds reporting/restart-window guards and CLI campaign parameters. The native equations, coefficients and forcing remain unchanged from `7bf794adeb45cda1e3640c52d5ae52e5f00da12f`. No push or scientific release approval.

Executable `/private/tmp/athenak-sgs-heating/src/athena`, SHA256 `b5f0fc19de082e07bc45c989cc24d5cc0c46dfd4d62e7994d8dd1efaac901197`: double precision CPU Release, OpenMP, no MPI, problem `sgs_verify`. Executed on the Mac (10 logical CPUs, 24 GiB RAM). Both 64³ lanes used four threads; thermal controls one each. Native source SHA dictionaries agree across verified and recovered segments despite the reporting commit change.

## Full-window result

Reference: [Schmidt–Federrath §4.2/Table 4](https://arxiv.org/html/1010.4492#S4.SS2), pure solenoidal 64³ row. Paper statistics average 2T–10T. All runs retain C1=.02, C2=.7, Cepsilon=1.5, Ckappa=.65. Our box has L=.5, V=1, T=.5, gamma=1.001, initial c0=1/5.5, instantaneous acceleration RMS=3, seed 42, RK2/LLF, CFL .1, block16, ng3. Forcing conventions differ; this is a controlled AthenaK comparison, not certified identical paper input.

| Quantity | 64³ PLM | 64³ PPM4 | Paper 64³ solenoidal |
|---|---:|---:|---:|
| Resolved volume RMS Mach | 5.42481 ± .19617 | 5.47433 ± .15333 | 5.38 |
| Volume RMS SGS Mach | 1.00262 ± .03654 | 1.10987 ± .03670 | 1.107 ± .053 |
| Mean K/(rho0 V²) | .0317117 ± .0021063 | .0367689 ± .0025773 | .0726 ± .0055 |
| Mean specific epsilon L/V³ | .188696 ± .014643 | .257875 ± .013766 | 1.236 ± .141* |
| Mean rho epsilon L/(rho0 V³) | .362124 ± .039233 | .456623 ± .048277 | 1.236 ± .141* |
| SGS fraction of corrected heating | 56.3535% | 69.7232% | approximately 75% |

Scatter is time-weighted temporal standard deviation, not SEM. Eight adjacent 1T blocks are descriptive and may be correlated. *Table 4 labels epsilon, while nearby prose identifies rho epsilon. Both explicit definitions are reported; neither matches quantitatively. No seed-ensemble uncertainty estimate exists.

PPM4 mean SGS Mach agrees closely, but mean K is 50.65% of the published value. SGS Mach measures volume mean 2K/(rho cs²), whereas K is energy density; agreement of Mach alone cannot establish energy-density agreement. The low mean K and dissipation remain. Neither averaging longer nor changing reconstruction removes these gaps; no closure coefficient or forcing amplitude was retuned.

[Eight convergence/block panels](averages-10T.pdf), [six heating panels](heating-10T/heating-comparison.pdf), [PLM five-page report](report64-plm/report.pdf), [PPM4 five-page report](report64-ppm4/report.pdf). PNG pages are available beside each PDF. Reports include cell distributions, density conditionals, CDFs, budgets and temporal trends.

## Heating and accounting

Signed native integrated differences use D/(Delta U−C), where C is applied resolved thermal-pressure work and D is applied SGS dissipation. N=Delta U−C−D is the numerical residual, including spatial/RK/shock-associated effects. Fractions use whole-window increments, not means of interval ratios. Endpoints are linearly interpolated in stored histories. No pre-ledger interval is included: the exact native pressure-ledger origin is around 1.9T.

| Case | Delta U | C | D | N | P | Delta K |
|---|---:|---:|---:|---:|---:|---:|
| PLM | 5.53075809 | .387616023 | 2.89833977 | 2.24480230 | 2.89564511 | −.00269465643 |
| PPM4 | 5.57510261 | .335161097 | 3.65345342 | 1.58648810 | 3.65474406 | .00129064205 |

Maximum absolute full-window represented-energy/forcing residual is 6.04e−14, and Delta K−P+D residual 4.89e−14. The full sampled-history residual normalized by initial represented energy is 1.57e−14 or smaller. These are numerical inventory checks, not proof of turbulent closure accuracy. Low K can coexist with substantial transfer because production and dissipation nearly balance.

## Frozen statistical screen

`stationarity-protocol.json` predates these full-window outcomes. Operational thresholds: absolute fitted fractional drift ≤10%, early 2T–6T versus late 6T–10T difference/abs(fullmean) ≤10%, and applied |P−D|/D ≤5%. These are our screening limits, not paper acceptance requirements. Zero/undefined normalizers cannot pass. Sound-speed growth is reported separately. The proposal in `reference-protocol` uses a different max(|P|,|D|) imbalance normalization; it was **not adopted**. The executable screen uses D throughout.

- PLM fails: SGS Mach fitted drift −10.0524%. K drift −5.6427%, density-weighted dissipation −8.1436%, and imbalance .09297% pass their screens.
- PPM4 fails: SGS Mach drift −10.5857%, density-weighted dissipation drift −12.3564%, and its early/late difference 11.1349%. K drift −8.6691% passes. Applied P/D imbalance .03533% passes.
- Both reports contain 82 samples spanning the exact ledger origin through10T. Maximum snapshot gaps are .10038T and .10032T, below the proposed .11T sampling limit. Eight full 1T block means are reported.

Do not relax the 10% limit to relabel either run stationary. Full averages remain valid descriptive measurements despite failed screens. They do not establish ensemble convergence.

## Thermal control: an established contributor

Two separate 32³ PPM4 runs finish 10T and pass 15/15 native checks each. Fixed initial c0, forcing, seed and SGS coefficients; only gamma changes 1.001→1.0001, adjusting initial pressure to preserve c0. Initial thermal reservoir U0=rho c0²/[gamma(gamma−1)] grows from 33.0248 to 330.5455.

| gamma | Mass-weighted cs² increase,2T–10T | SGS Mach fitted drift | K fitted drift |
|---|---:|---:|---:|
| 1.001 | 16.4468% | −8.6424% | +1.1612% |
| 1.0001 | 1.6761% | −2.9029% | +.9148% |

Both controls pass the operational screen. The thermal reservoir changes fractional sound-speed growth strongly, while K trends are similar. This supports heating as a contributor to Mach drift. Coupled trajectories can diverge; this one-seed 32³ control is not a quantitative attribution of the 64³ K deficit or proof that gamma alone causes it.

[Four-panel thermal comparison](thermal-sensitivity.pdf), [gamma1.001 report](thermal32-gamma1001/report.pdf), [gamma1.0001 report](thermal32-gamma10001/report.pdf).

## Forcing and reference arithmetic

[Independent forcing oracle](forcing-audit/README.md): 14/14 checks pass. 65,536 independent Gaussian fields with the same native spectrum have ensemble RMS 2.999923 and instantaneous RMS scatter .180294; the native driver fixes instantaneous RMS to 3 and removes a density-weighted mean. These are different stochastic laws. Native final radial amplitude peaks at mode 1.5, despite npeak 2. This audit establishes departures, not their causal contribution to low K. Pinned modern Enzo forcing was introduced after the paper; the original Enzo 1.5 patch/input/seed was not recovered in the bounded primary-source search. No claim that acceleration6 is the confirmed Table 4 target, and no3→6 tuning.

[Gaussian/Favre analytic pilot](reference-protocol/README.txt): 18/18 checks pass, maximum absolute error 1.30e−14. Checks include analytic Fourier attenuation, unequal-density Favre covariance, centered real-space covariance, stress symmetry/positive covariance and Galilean invariance. Chosen Gaussian width is DeltaG=sqrt(12)sigma; the paper-identical width convention remains unknown. This validates arithmetic on analytic fields only. Existing top-hat comparisons are not a substitute for the paper Gaussian reference. No adequately resolved turbulent reference was executed here.

## Execution, failures and costs

Original64³ continuations were externally interrupted at 8.40028T (PLM) and 7.70009T (PPM4). Their `results.json` retain FAIL for incomplete requested10T. Every saved checkpoint was independently validated in separate `*-verified-segment` packets (56 and49 checks). Recovery from the pinned checkpoints finishes10T, with15/15 native checks each. Assembly checks exact inventories/cumulative counter boundaries and strictly increasing times; it never rebases counters. Verified-segment plus recovery records are separate from the failed original attempt. Sources, executable and checkpoint hashes are retained.

Thermal controls contribute 30 passing native checks; the two completed64³ recovery segments 30. These 60 checks do not include older prefixes, partial-checkpoint checks, or analytic-oracle checks. Restart join/statistical fixtures have3 passing unit tests; diagnostic, time-weighted-statistics and signed-heating self-tests pass. Native 16³ basename regression and seed 43/compressive CLI smoke pass. New source Python Ruff and whitespace checks pass. Original artifact scripts and failed probes remain available; the partial validator was formatted after execution with an executed text copy retained.

The 128³ startup pilot completes .2T in 527.907 seconds (two threads, 143 cycles); each checkpoint 348,937,732 bytes. Its original wrapper FAIL (stale basename lookup) is retained; independent native-file validation passes15/15. The runner fix forces the requested basename after overrides and passes a native regression. A preceding attempt failed the unused dedt/acceleration parameter guard. No failed attempts were erased. Linear startup extrapolation gives at least 7.33h for 10T at two threads; this is not a measured developed-turbulence runtime or reason to call the pilot a scaling result.

## Reproduce and remaining work

Existing Python environment: `/private/tmp/sgs-plot-venv/bin/python` (NumPy/Matplotlib), with `MPLCONFIGDIR=/private/tmp/sgs-mpl-cache`. Use the committed `tst/sgs/README.md` for native build/tests and the restart-report CLI. Artifact `assemble_recovery.py` joins already-verified native segments; `compare_10T.py` and `plot_thermal.py` regenerate figures from stored diagnostics. To regenerate heating:

```sh
SGS_SOURCE=/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu
SGS_DATA=/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-validation-10T
MPLCONFIGDIR=/private/tmp/sgs-mpl-cache /private/tmp/sgs-plot-venv/bin/python -B "$SGS_SOURCE/tst/sgs/heating_plots.py" --runs "$SGS_DATA/report64-plm" "$SGS_DATA/report64-ppm4" --start-T 2 --end-T 10 --output /private/tmp/sgs-heating-reproduction
```

Use `--start-T 2 --end-T 10` (separate flag/value tokens). Output directories must be fresh. `campaign-commands.sh` is syntax-checked runnable preparation for seeds43/44,64³ compressive,128³/256³ resolution runs; those campaigns were **not executed**. It requires explicit executable/Python/output locations and suitable resources. At 256³,100 ghost-inclusive checkpoints would consume roughly 280GB; 1024³ active reference fields require 56GiB before FFT/tensor temporaries. Proper Gaussian turbulent reference, resolution laws, seed uncertainty and original paper forcing provenance remain open; they need cluster/original data and a specified compute budget. The earlier compute-location question has no answer; no remote execution was inferred. AMR/MHD/cooling/MPI/GPU qualification remains outside this HD campaign.

Compact wiki copies retain raw paths/hashes via `archive-manifest.json`. Raw checkpoints remain in the external artifact tree and prior prefix archives. Do not claim all paper tests or scientific release gates are established.
