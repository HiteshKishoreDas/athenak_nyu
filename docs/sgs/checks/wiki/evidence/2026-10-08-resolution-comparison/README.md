# Filtered high-resolution turbulence comparison

Code commit: `b932230caa7ae749ae43309def77ee1232621dac`, branch
`codex/multiphase-schmidt-federrath-hd`, based on `multiphase_turb`.

A 128² SGS-off flow is driven with the existing modal forcing to t=2 (seed 42,
mixed solenoidal/compressive forcing, modes 1–3, kz=0). Four unforced decay lanes
then evolve to t=.5: 128² reference, 64² reference sensitivity, 32² SGS-off and
32² SGS-on. All retain three velocity components. The final workflow, including
nine figure pages / 41 panels, takes 35.35 seconds on the local serial CPU build.

Aligned conservative block averaging sets the coarse initial density, momentum,
total energy and scalar masses. Removed kinetic variance becomes SGS energy K in
the SGS-on lane; its resolved total energy is reduced by K to retain filtered
thermal energy. In the SGS-off lane this variance becomes heat. Thus the coarse
lanes share mass, momentum and represented E+K, but have different initial thermal
pressures. All saved native initial states match the supplied values exactly.
Reference diagnostics are filtered to the same 32² grid before comparison.
Intermediate curves interpolate conserved fields to common times; the final
comparison uses the exact common stopping time without interpolation.

| Final relative RMS error | 32² SGS-off | 32² SGS-on |
|---|---:|---:|
| Density | 8.3407% | 8.3133% |
| Velocity | 40.2507% | 40.2996% |
| Pressure | 13.7325% | 13.5851% |
| Vorticity | 41.3469% | 41.6584% |
| Subgrid kinetic energy K | 100% | 83.8375% |

Errors are normalized by reference RMS, not fluctuation standard deviation.
This case does not demonstrate an overall SGS benefit. Density/pressure improve
slightly; velocity/vorticity worsen slightly. Modeled K underestimates the filtered
reference. The filtered 64² and 128² velocities differ by 21.14% relative RMS,
so the reference is not demonstrated converged. This is one 2D developed-flow
seed followed by decay, not stationary turbulence, converged DNS or 3D physical
qualification. No acceptance condition requires SGS to improve the flow.

All six execution/initialization result records pass. Existing Release and bounds
Debug unforced suites (46 records each), configuration suites (32 each), closure
checks, new initializer round-trip and eight malformed-input rejection cases pass.
A bounds Debug 64²/16² short comparison also passes six records. All 520 existing
native regression checkpoints are byte-identical to the preceding results.
Added-line C++ checks and the new runner's Ruff E/F checks pass; no whole-repository
style pass is claimed. No numerical threshold was relaxed.

The source change is restricted to a test-only initial-state loader in sgs_verify,
a comparison runner, an explicit optional larger test grid cap and documentation.
Solver, task, restart and numerical physics machinery remain unchanged.

- `report.pdf`: all nine plotted pages.
- `verified-summary.json`, `comparison.json`: numerical results and qualification.
- `results.json`, `manifest.json`: acceptance records and run/source fingerprints.
- `implementation.patch`: complete committed change.
- `qualification.tar.gz`: full final inputs, initial-state payloads, native outputs,
  histories, logs, plots and HTML gallery; all fresh regression packets, Debug
  smoke outputs, build caches and development/check logs. Extract and open
  `sgs-resolution-comparison/index.html`. It also includes pilot summaries.
- `SHA256SUMS.json`: hashes of the files in this evidence folder.

Build, environment, run, plot-only and self-test commands are in the code
repository's `tst/sgs/README.md` under the high-resolution comparison section.
The full current gallery is saved at
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-resolution-comparison/index.html`.
Final main-run source/harness fingerprints match the committed source. Earlier
regression manifests retain their original then-uncommitted provenance rather
than being relabelled with a future commit. First Ruff formatting checks found
long lines, which were fixed before final verification. Numerical runs passed.
