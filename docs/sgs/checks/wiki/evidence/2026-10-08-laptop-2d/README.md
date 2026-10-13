# Local laptop 2D SGS tests and plots

Code commit: `1c35d3eecaebae9339c7ea671b6fa05454d8925e` on `codex/multiphase-schmidt-federrath-hd`.
No C++ solver changes. One visual test runner was added; the existing forcing
runner accepts dimensions/overrides while retaining its default 3D behavior.

Ran 22 serial Release simulations on 16² and 32² grids: uniform decay at three
CFLs, SGS-on/off shear and boxes at two resolutions, diffusion, isolated contact
advection, zero K, block-layout/restart lanes, and solenoidal/compressive driven
pilots. All 25 result records pass. The workflow including 15 figures/92 panels
completed in 10.61 seconds in the final local run.
The report PDF was checked for all 15 pages and interpretation captions; plot-only
regeneration was executed successfully. Grid orientation, pressure, periodic
vorticity and FFT Parseval self-checks pass. New-file Ruff E/F checks at 90 columns
pass; existing-file baseline style was not reformatted. Default 3D forcing still
passes all 19 records and its 31 checkpoints are byte-identical to pre-change runs.

Analytic decay temporal orders: 2.00925 and 2.00137. Contact density L1 errors:
.00555285 (16²) and .00168528 (32²). Independent unforced E+K drift reaches only
1.14492e-15 in the plotted checks. Final forced represented-energy/work residual
is at most 2.35922e-16. Restart and decomposition field errors are exactly zero.
Final 32² solver elapsed-time ratios (SGS-on/off) are approximately 4.0 for shear,
4.8 for the decaying box, 4.2 for solenoidal driving and 4.3 for compressive driving.
The underlying Kokkos timer measures elapsed time even though AthenaK's log labels
it CPU time. These are single short measurements with output, not controlled
benchmark medians.

- `report.pdf`: 15 plotted pages, directly viewable.
- `verified-summary.json`, `results.json`, `manifest.json`: measured diagnostics,
  acceptance results and source/executable fingerprints.
- `implementation.patch`: complete runnable workflow and README instructions.
- `laptop-2d-evidence.tar.gz`: all 22 final-case inputs, raw checkpoints, histories,
  logs, results, plots and HTML gallery; build cache and calculation/lint/regression
  logs. Extract and open `sgs-2d-local/index.html`.

The complete local gallery is at
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-2d-local/index.html`.
The repository README provides build, environment, run and plot-only commands.
Every case's saved input is linked in the gallery for independent reruns.
Source and harness fingerprints match the final code commit. The executable was
built from unchanged C++ source before adding the Python report workflow; its
production hashes match the new commit.

Initial development runs were successful and remain in the local visualization
folder. Captions and correct elapsed-time labels were added before the final run;
one formatting-only lint pass corrected long lines and an unnecessary f-string.
The first environment setup could not access the existing package cache from the
sandbox; cached packages were then installed into a temporary isolated environment
with the authorized sandbox escalation. No physical threshold was relaxed.

Full-model resolution changes also change the SGS filter width. The contact
oracle disables dissipation/diffusion. Forced 2D pilots retain three velocity
components, use kz=0 and stop at t=1. None of these short runs establish stationary
turbulence, physical closure accuracy, MPI/GPU support or completed M0–M9 gates.
