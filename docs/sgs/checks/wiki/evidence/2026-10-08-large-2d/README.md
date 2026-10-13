# Larger 2D SGS turbulence comparisons

Code commit `675c74606efe7dff75f771e6b3979a38fbed9820`, branch
`codex/multiphase-schmidt-federrath-hd`, based on `multiphase_turb`.

Two full local comparisons drive separate 256² and 512² SGS-off preparations to
t=2, then decay to t=.5. Each uses a half-resolution reference check and 64²
SGS-off/on lanes initialized through aligned conservative/Favre filtering. The
SGS-on lane stores removed kinetic variance as K and preserves filtered thermal
energy; SGS-off stores it as heat. Mass, momentum, scalar masses and represented
E+K match initially. RK2/LLF/PLM, CFL=.1, gamma=5/3, seed 42, forcing parameters,
published coefficients and intervals match previous tests. Reference-check
initializations come from their own full-resolution preparation snapshot.

| Reference / coarse | Velocity off / on | Vorticity off / on | Mean K recovered | Half/full reference velocity difference |
|---|---:|---:|---:|---:|
| 128² → 32² | 40.251% / 40.300% | 41.347% / 41.658% | 36.10% | 21.14% |
| 256² → 64² | 39.619% / 39.689% | 41.020% / 41.113% | 20.45% | 22.48% |
| 512² → 64² | 50.118% / 50.100% | 50.414% / 50.455% | 7.96% | 20.69% |

Errors use reference RMS, not fluctuation standard deviation. Separate
on-off-difference files measure SGS-on against SGS-off directly, normalized by
SGS-off RMS; their maps use individual color ranges. Mean K recovered
is modeled mean K divided by filtered reference mean K, distinct from the cellwise
K error. Final diagnostics compare at the exact common time and common 64² filter.
Larger grids do not demonstrate a substantial SGS advantage. Direct final SGS-on/off
velocity differences are 0.779%
and 1.049% of SGS-off RMS
for the 256² and 512² preparations respectively. These measure the actual field
difference, independently of errors against the reference. Modeled K drains because
its applied dissipation exceeds production; stresses depend on K and sqrt(K).
The SGS production path does not automatically collect every kinetic-energy loss
from the hydro solver. Numerical smoothing is a plausible contributor to the weak
coupling, but pressure work and implicit numerical loss have not been separately
measured, so this is not a demonstrated root-cause diagnosis. References remain
resolution sensitive and neither is established converged DNS. Across the old
128²/32² run and new tests, both coarse filter width and preparation resolution
change, so the overview is not a controlled resolution-convergence proof. These
remain one-seed developed-flow decay tests, not stationary statistical qualification.

Each new workflow passes all six acceptance/initialization records. Native
initial states read back exactly, all cells remain admissible, represented energy
and periodic conserved integrals pass the existing thresholds, and the K ledger
closes. All stress/diffusion limiter and fallback counters remain zero. Each full
run has nine plotted pages /41 panels; two reports total 18 pages /82 panels, plus
six overview panels and eight direct on/off map panels. All plotted pages were visually inspected and PDF page
counts checked. Production solver physics is unchanged: source edits only raise
optional test-grid/runtime caps and document them, preserving defaults and the
3D cap. No numerical acceptance threshold changed.

The Apple M4 laptop uses a separate double-precision Kokkos OpenMP build with
eight threads and the already-installed libomp. The independent closure oracle,
46 evolution records, 32 configuration records and 2D/3D initializer self-checks
pass. All 244 matched native regression checkpoints have exactly equal active
fields to the serial build. The 512² serial/parallel pilots pass six records each;
reported diagnostics agree within 1.8e-14. Parallel preparation takes 40.32 seconds
versus 98.63 serial in the short pilot. The full 256² workflow takes
256.02 seconds and 512² 1418.65 seconds,
including analysis/plots. These are single elapsed runs, not benchmark medians;
the end of the parallel pilot overlaps the start of 256² preparation. Full 256²
and 512² workflows do not overlap. The larger source snapshot and harness
fingerprints match the clean committed code. Build/runtime commands and binary
hash are saved in commands.txt/environment.json. No dependency was installed.

Initial Ruff invocation used its default 88-column limit and flagged five existing
89/90-column lines. A subsequent 90-column run could not write its cache in the
restricted worktree; the final run with --no-cache and the repository 90-column
limit passes (ruff-final.log). No source changes were needed. Plot page-count inspection first used the isolated plotting
Python, which lacks pypdf, then successfully used bundled Python. These tooling
attempts do not alter or invalidate numerical runs. Earlier pilot manifests retain
their actual pre-commit provenance.

The raw archive remains outside wiki Git, with its path, byte size and SHA256 in
sgs-large-2d-raw-manifest.json. It preserves both complete workflows, native
checkpoints, ASCII initial states, histories, inputs, logs, all figures and HTML
galleries; build/verification evidence; the compiled Apple-arm64 binary; and
compact pilot records. Original pilot raw checkpoints remain in their output
folders. The saved binary dynamically depends on the installed macOS libomp.

- 256-report.pdf /512-report.pdf: complete reports.
- sgs-large-2d-overview.png/pdf/json: comparison with prior 128²/32² result.
- *-verified-summary.json: independently checked energy ledgers and metrics.
- *-results.json /manifest.json: execution checks and source fingerprints.
- implementation.patch: minimal test-runner extension.
- commands.txt /environment.json: reproduction and execution environment.
- sgs-large-2d-raw-manifest.json: complete raw-file hashes and archive location.
- SHA256SUMS.json: hashes of this compact evidence packet.

Gallery paths:
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-2d-256-reference/index.html`
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbe-7042-7051-b251-90b2d08343eb/sgs-2d-512-reference/index.html`
