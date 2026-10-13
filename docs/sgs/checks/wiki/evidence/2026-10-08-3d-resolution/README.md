# Three-dimensional SGS turbulence comparison

Code commit `eac5743afb44c53dd47dc463f4e7f014a0113ca8`, branch
`codex/multiphase-schmidt-federrath-hd`, based on `multiphase_turb`.

## Executed protocol

A serial double CPU Release run drives an SGS-off 64³ flow to t=2 using native
modal OU forcing, seed 42, dedt=1, sol_fraction=.5, physical modes 1–3 in all three
directions. Conservative aligned averages initialize 64³ and 32³ SGS-off
references and 16³ SGS-off/on lanes. All then decay without forcing to t=.5.
RK2/LLF/PLM, CFL=.1, gamma=5/3 and published SGS coefficients are fixed.
SGS-on stores removed velocity variance as K and retains filtered thermal energy;
SGS-off stores the variance as heat. Their initial mass, momenta, scalar masses
and represented E+K match. All supplied initial states read back exactly.

The full 3D workflow takes 482.27 seconds, including 10 plotted pages /47 panels.
Solver elapsed timings: preparation 376.50s, reference 75.53s, 32³ reference 3.37s,
16³ off .192s, 16³ on 1.602s. These are single local elapsed timings with output,
not benchmark medians. The short 2D control was run during 3D preparation; none
of its solver runs overlapped 3D decay timing. An independent 64²/16² control
uses the same duration, coefficients, forcing strength and seed parameter,
but dimensionality changes the forcing realization and dynamics. It finishes
in 9.09 seconds. This is not an isolated causal dimensionality experiment.

All quantitative errors, energy means and spectra use complete 3D volumes.
Maps show central z-cell slices: native/coarse cell centers differ, while filtered
reference and coarse lanes share a slice. Vorticity errors compare all components.

## Measured outcome

| Final relative RMS error in 3D | SGS off | SGS on |
|---|---:|---:|
| Density | 9.3000% | 9.2540% |
| Velocity | 77.9380% | 77.7077% |
| Pressure | 15.5361% | 15.3569% |
| Vorticity vector | 91.3994% | 91.3690% |
| K | 100% | 75.0835% |

Going to 3D does not produce a substantial SGS benefit in this case. Mean modeled
K is .0141399 versus .0338485 in the filtered 64³ reference (41.77%). Its initial
K=.0642307, applied production=.0249539 and dissipation=.0750448 close the global
K ledger. No stress/diffusion limiter or fallback was used. In the matched 2D
control velocity errors are 47.6759% off /48.2028% on; mean modeled K is 47.57%
of its filtered 64² reference.

The 32³/64³ reference sensitivity is 52.48% for velocity and 80.72% for K. Neither
reference is established converged DNS. Preparation's late kinetic drift is
7.07%; it is developed flow, not statistically stationary qualification.
One seed and short decay do not complete M7–M9 or demonstrate general SGS benefit.

With exact filtered K AND Favre velocity supplied, the final cell-centered
stress-production diagnostic has correlation .92355 and relative RMS error 38.12%.
Mean production is overestimated by 28.18%, production sign is wrong in 13.33%
of cells and deviatoric stress error is 50.04%. Thus correlation is good, but
the closure is not established accurate. This diagnostic uses block-filtered
covariance and centered coarse gradients, not the native face production stencil.
It suggests examining coupled evolution, numerical smoothing, native discretization
and source balance; it does not identify the cause of the K deficit.

## Verification and provenance

Both full comparisons pass six numerical/initialization records each. Native
initializer checks pass 48 cases across 2D/3D, SGS on/off and Release/bounds Debug,
including exact varied multiblock readback and ten rejection cases. Independent
3D curl, FFT Parseval, averaging and stress trace checks pass. Existing Release
and bounds Debug suites pass 46 unforced +32 configuration records each; closure
checks pass. All 520 existing native checkpoints are byte-identical to the prior
results. All-cell numerical admissibility and conserved inventories pass; maximum
independent history represented-energy drift is 5.56e-14. No threshold was relaxed.
Added C++ style checks, Python syntax and new-file Ruff E/F checks pass; existing
whole-file cpplint baseline diagnostics remain.

Production physics is unchanged. Changes are restricted to the test initial-state
loader, dimension-aware comparison/reporting, a-priori diagnostic and explicit
optional runtime caps retaining existing verification-runner defaults. Source and
harness fingerprints of both full runs match the committed source. Both PDFs were
checked for expected page counts, and all 3D figure pages were visually inspected.
Independent review of the implementation and final interpretation found no blocker.

The first pilot mistakenly supplied a negative min_kz, which the native driver
rejects. It remains a failed attempt. Correcting to the driver's nonnegative mode
convention gave a successful 35.26-second pilot (all six records pass). The initial
new timeout signature caused an import error before execution; corrected before
pilot/full runs. An initializer rejection-test draft used unsupported policy
'constant'; it was corrected to supported 'uniform', then all48 checks reran.
No numerical equations or thresholds changed to obtain a pass.

## Files and reproduction

- `report-3d.pdf`: ten pages /47 panels from the full 3D run.
- `report-2d-control.pdf`: nine pages /41 panels from the matched-cell 2D control.
- `verified-*.json`, `results-*.json`, `manifest-*.json`: measured data, acceptance,
  source and executable fingerprints, applied transfers and qualification.
- `dimensional-comparison.*`: direct comparison of resolved errors and energies.
- `implementation.patch`: exact focused code commit.
- `verification.tar.gz`: regression inputs/results/logs, native initializer data,
  independent mathematical review, build caches/logs and the rejected pilot.
- `raw-data-manifest.json`: per-file hashes and absolute locations for all final
  raw inputs, checkpoints, histories, logs, plots and their complete archive.
- `SHA256SUMS.json`: hashes of the files saved here.

The full raw archive is stored beside the persistent visualization outputs to
keep hundreds of megabytes of numerical state out of the wiki Git history:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-3d-final-evidence.tar.gz`.
The hash manifest binds that archive to this record. Extract and open
`sgs-3d-resolution-comparison/index.html` or `sgs-2d-matched-control/index.html`.
The corrected pilot's raw data remain at the sibling `sgs-3d-pilot-corrected` path.
Run/build/self-test instructions are in `tst/sgs/README.md`, section '3D comparison'.
Use a new output folder; the 3D laptop runner is capped at 64³ and the full run used
an explicit 900-second cap per simulation. No GPU/MPI qualification is added.
