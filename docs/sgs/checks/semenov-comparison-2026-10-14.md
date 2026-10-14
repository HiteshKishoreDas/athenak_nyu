# Semenov comparison: evolved initial field, 10T decay and energy integration

The new unforced tests use our own evolved turbulence snapshot, as authorized on 2026-10-14. They keep the original Schmidt–Federrath coefficients and native solver unchanged. SGS captures a useful kinetic-energy partition in this pilot, but its late unresolved-energy prediction remains low. All reference-convergence screens fail, so this is not a scientifically qualified reproduction of the paper.

## Mean unresolved energy

The reported mean ratio is the integral of volume-mean SGS energy divided by the corresponding integral of filtered reference energy, over 1–10 initial crossing times. It is not the average of instantaneous ratios. Initial K equality is imposed by initialization and is excluded from validation.

| SGS grid | SGS / 128³ reference, average 1–10T | Ratio at 10T | Maximum 64³ versus 128³ reference sensitivity, 1–10T | Fixed 10% reference screen |
|---|---:|---:|---:|---|
| 8³ | 0.8669 | 0.4814 | 22.27% | Fail |
| 16³ | 0.8305 | 0.5272 | 37.62% | Fail |
| 32³ | 1.1016 | 0.6042 | 57.90% | Fail |

The 32³ ratio is about 1.28 at 1T, 1.07 at 4T and 0.60 at 10T. This temporal dependence is inconsistent with a universal omitted factor of 0.5. The high-resolution estimate itself changes materially with resolution, including on the larger 8³ filtering grid; reducing coarse-grid size does not establish convergence in this campaign.

![Mean SGS energy and reference sensitivity](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/decay-mass-rms/energy-comparison.png)

Read the lower panels first: unity means agreement; the other curve shows how much the unresolved-energy estimate changes when the reference is reduced to 64³. The grey region excludes startup. Agreement of a time average can hide early excess and late deficit.

## What SGS improves in this pilot

Over 1–10T, the resolved kinetic-energy ratio to the correspondingly filtered 128³ reference is 1.006 with SGS versus 1.240 without SGS at 16³, and 0.993 versus 1.141 at 32³. Thus SGS improves the partition between resolved and unresolved motions in these runs. These are integrated-energy comparisons with an unconverged numerical reference, not a claim of pointwise accuracy or universal benefit.

For total resolved-plus-SGS kinetic energy, the corresponding ratios are 0.958 versus 0.903 at 16³, and 1.007 versus 0.998 at 32³. The latter illustrates why total kinetic energy alone does not establish correct energy partition.

![Kinetic-energy comparison and local run cost](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/generated-energy-and-cost.png)

The cost bars are single local measurements, including output, with eight OpenMP threads. They are not a formal performance benchmark.

## Spatial comparison and distributions

The exact filtered-out velocity RMS is sqrt(2 K_filtered/rho_filtered); the model equivalent is sqrt(2 K_SGS/rho_SGS). The reference quantity uses the density-weighted velocity fluctuations removed by the aligned block filter. Matching native snapshot hashes were checked before plotting.

At 32³, same-cell K correlations at 1T, 4T and 10T are 0.502, 0.224 and 0.265. The corresponding highest-energy-decile overlaps are 42.1%, 22.0% and 25.4% (random overlap is 10%). Model maps are smoother. Late same-cell errors include different flow trajectories and shock positions; they are not an isolated constitutive-closure test. Snapshot offsets are at most 0.0031T for these pairs.

![Removed velocity versus corresponding SGS-on energy](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/generated-removed-velocity.png)

Each column uses a common full-range velocity color scale for reference and model. The bottom row compares every positive-energy cell; the red line is equality. All-cell scalar correlations are recorded separately from the logarithmic plot.

![Density and unresolved-velocity relation](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/generated-density-sigma.png)

Curves show medians in density bins, with 16th–84th percentile bands. These compare distributions, avoiding a strict assumption of identical late trajectories.

![Unweighted velocity spectra](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/generated-velocity-spectra.png)

These additional spectra use unweighted velocity, compensated by k^(5/3). The original runner's `local-statistics.png` instead shows the shell spectrum of sqrt(rho) times velocity; the two conventions are not interchangeable.

## ART energy-update comparison

This is a comparison against published methods; a working public source checkout for Semenov's ART version was not found and no ART executable was run.

[Semenov 2025, sections II.2–II.4](https://arxiv.org/html/2410.23339v3) describes the following distinctions:

| Operation | Published ART method | Current AthenaK port |
|---|---|---|
| Adiabatic SGS transport | Advect a nonthermal entropy proxy | Advect K and apply compatible pressure work |
| Explicit production | Residual of total minus kinetic energy from added SGS fluxes; explicit source fallback for small K | Stage-consistent face-work production and paired E/K exchange |
| Decay | Analytic split update after each step | Two-stage RK2 source evaluation |
| Turbulent diffusion | Diffuse energy density; Ckappa=0.4 | Diffuse specific energy K/rho; Ckappa=0.65 |
| Thermodynamics | Isothermal decay test | Ideal gas with gamma=1.001 |

Both use C1=0.02 and C2=0.7 for the nonlinear closure. The prominent mean-energy calibration curves are for SLES, whose source is numerical dissipation, rather than the explicit LES port tested here. Cepsilon=1.6 is the published SLES calibration; it was not substituted into our Cepsilon=1.5 baseline. The paper's 1024³ reference and exact initial fixture remain unmatched.

The corresponding [entropy method paper, Semenov et al. 2022](https://arxiv.org/html/2107.14240), explains why directly evolving nonthermal energy can create nonadiabatic changes at shocks. We tested that concern independently in the existing native solver.

## Native compression test: a numerical qualification gap

These controls retain isotropic SGS pressure but set C1=C2=Cepsilon=Ckappa=0 to isolate adiabatic compression. Initial K/rho^(5/3) is spatially constant. Therefore the independent adiabatic expectation at each later cell is K_ad=S0 rho^(5/3). An ordinary passive scalar with uniform specific value supplies a same-code entropy-advection proxy; its maximum relative error is below 7e-15.

Four periodic shock-pair configurations and a smooth compression control were run on 64, 128 and 256 cells. These are related tests, not the paper's cosmic-ray shock tube: SGS gamma is fixed at 5/3, boundaries are periodic and initial SGS energy is dynamically tiny. No stress limiter, diffusion limiter or fallback activates in these controls.

At 256 cells and CFL 0.1:

| Thermal gamma | Mass-weighted SGS entropy RMS error | Whole-domain K excess over adiabatic expectation | Maximum local entropy excess |
|---|---:|---:|---:|
| 5/3, strong shock | 6.55% | 0.81% | 40.75% |
| 1.1, compressible shock | 25.85% | 3.93% | 207.46% |
| 1.001, near-isothermal shock | 33.70% | 5.52% | 251.97% |

Smooth-flow entropy RMS errors fall from 7.35e-5 to 2.61e-6 over the tested refinement. Shock errors persist when CFL decreases from 0.2 to 0.05 and do not decrease over the tested compressible-shock grids. Thus the energy formulation has a shock-compression qualification gap beyond the earlier smooth-flow and initial-rate equation checks. Its measured sign in these controls raises K; it does not explain a 50% mean-K deficit in the turbulent box. This does not establish the magnitude of the effect with the full nonlinear production and dissipation model in 3D turbulence.

![Native compression profiles](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/compression-profiles.png)

Unity is constant entropy. The ordinary passive entropy proxy follows that expectation to rounding accuracy. No ART solution is plotted.

## Analytic decay versus RK2

The uniform-decay controls retain the original coefficients and compare native K with the published analytic-decay formula evaluated using our Cepsilon=1.5. Density is constant; total evolution time is 0.25 and dx=1/16. Relative errors at CFL 0.2, 0.1, 0.05 and 0.025 are 0.005890%, 0.001460%, 0.000364% and 0.000091%. The errors shrink by approximately four on halving CFL. This isolates decay integration and rejects a universal factor-of-two source-weight error in this control; it does not compare the full coupled ART and AthenaK timestep operators.

![Compression and decay diagnostics](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/compression-and-decay.png)

The shock and decay questions have different behavior: shock entropy error persists with smaller timesteps, while uniform decay converges toward the analytic result.

## Initial conditions and interpretation limits

The parent is our 256³ SGS-off decay snapshot `reference.hydro_u.00004.bin`, from approximately 2.001 original crossing times. It contains evolved density structure and a cascade generated from the earlier low-mode random field. It is restricted to 128³, normalized to mean density 1, given zero mass-weighted bulk velocity, rescaled to global Mach 8.6 and assigned sound speed 1. The recorded volume RMS velocity is 9.3204; these Mach definitions differ because density is nonuniform. The initial kinetic energy is 36.98 and 40.3% of unweighted velocity spectral energy is above shell k=4. Density standard deviation is 1.143 and log-density standard deviation is 1.010.

This is an evolved, structured initial field. It is not certified statistically stationary at Mach 8.6: resetting temperature and rescaling a lower-Mach decay snapshot alters its original density–velocity balance. It should be described as a local Semenov-style surrogate, not an exact developed forced-state reproduction.

All decay lanes use RK2, LLF, monotone PPM4, gamma=1.001, CFL=0.2, no forcing and coefficients C1=0.02, C2=0.7, Cepsilon=1.5, Ckappa=0.65. T=1/(2*8.6). Coarse rho, momentum and represented energy are conservative averages. SGS-on stores exact removed kinetic variance in K and subtracts it from E; SGS-off leaves that variance as heat. Paired initial rho/momentum are identical and represented-energy errors are below 1e-12. Initial equality of K with the reference is imposed.

The final reference Mach is 0.987 and sound speed is 1.018. The late window therefore includes transonic/subsonic decay. This differs from a stationary supersonic calibration. There is one realization, no ensemble or confidence interval, a near-isothermal rather than exact-isothermal EOS, and only 128³/64³ numerical references. All fixed 10% reference screens fail. No coefficient was retuned to hide these limitations.

The original [Kritsuk 2011 paper](https://arxiv.org/pdf/1103.5525) identifies the historical KITP fixture and describes its STAGGER parent/regridding; Semenov describes a PPM-prepared initial state. The exact variant was not recovered. The legacy Wikispaces fixture URL was inaccessible, and the inspected [author website](https://vadimsemenov.com/index.html) and papers supplied no working fixture/source download. After this search the user explicitly authorized locally generated data.

## Validation and retained attempts

- Final compression suite: 43 native runs, 774 native checks plus 118 specialized checks, all engineering checks passed. Scientific entropy errors are reported above and are not turned into passes.
- Main six-lane mass-normalized 3D decay: 488 checks, zero failures, all lanes completed 10T.
- Supplementary bandwidth campaign: two reference runs were hash-checked/revalidated; four coarse runs were executed, including a repeat of 16³. Its 488 checks passed. All scientific reference screens failed.
- Sixteen additional paired-IC, fixed-coefficient and unforced-input checks passed.
- Snapshot normalization self-checks cover nonuniform density, global Mach, zero momentum, Galilean boosts and invalid states.
- All 324 native source hashes and the executable hash match the prior tested solver. No native source or default coefficient changed.

Executable SHA256: `b5f0fc19de082e07bc45c989cc24d5cc0c46dfd4d62e7994d8dd1efaac901197`. Starting branch commit: `21bfe7a6`.

The first six-lane pilot in `decay/` used volume RMS Mach 8.6 and consequently global Mach 7.935. It is retained with its executed harness in `decay_turbulence_volume_executed.py`; its results are not used for the corrected headline comparison. The original 34-case compression pilot in `native/` is also retained; `native-expanded/` is the final 43-case suite.

Raw evidence, source fingerprints, protocols, histories, initial states and plots are retained in [sgs-semenov-comparison-2026-10-14](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14). `decay-mass-rms/` is the corrected main campaign; `decay-bandwidth/` holds the larger filtering-scale control. `additional-verification.json`, `spatial-results.json` and `benefit-results.json` contain the extra measurements.

## Reproduce locally

Use the existing CPU-double `PROBLEM=sgs_verify` build and NumPy/Matplotlib environment. Output directories must be new. The [test README](/Users/hitesh/.codex/worktrees/multiphase-sgs/athenak_nyu/tst/sgs/README.md) documents the build and runner commands. `decay_turbulence.py --initial-snapshot` accepts a native SGS-off snapshot and records its origin, hash and normalization. `verify_adiabatic_sgs.py` runs the compression/decay controls; `plot_adiabatic_sgs.py` renders their saved JSON without launching simulations.

[All plot gallery](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/index.html) · [Main five-page decay PDF](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/decay-mass-rms/report.pdf) · [Larger filtering-scale PDF](/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-semenov-comparison-2026-10-14/decay-bandwidth/report.pdf)
