# Laptop tests and visual checks

Run from this development checkout, branch `codex/shock-multiphase`.
The runner reuses the existing integration/qualification tests and output readers.
It adds no production hooks, dependencies or numerical acceptance thresholds.

## One command with the existing local Release build

```sh
cd /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu
MPLCONFIGDIR=/tmp/athenak-mpl \
  /Users/hitesh/hitesh/git/shockit/.venv/bin/python tst/shocks-mc/laptop_visuals.py \
  --executable build-release/src/athena \
  --workdir /tmp/athenak-laptop-check-new \
  --output /tmp/athenak-laptop-plots-new
open /tmp/athenak-laptop-plots-new/index.html
```

Both directories must be new; change their suffix each run. `--reuse` is only for
reusing completed integration/qualification reports after an interrupted plotting
attempt; additional snapshot/benchmark directories still must be absent.
Python needs NumPy, Matplotlib and h5py (AthenaK's binary reader imports h5py).
The command above uses the already installed ShockIt environment. Tested versions:
Python3.14.0, NumPy2.5.0, Matplotlib3.11.0, h5py3.16.0; no package installation was performed.

For another laptop/check-out, replace the Python path with an environment containing
those packages, and build the local Serial CPU executable first:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DAthena_SHOCK_TESTS=ON
cmake --build build-release -j8
build-release/shock_encounter_test
```

This is a small thin-2D planar test family, not general 2D shock geometry.
Twenty physical fixtures use four y cells; one held-out fixture is thin3D.
There are20 engineering subprocesses,21 physical fixture runs,1 additional full-grid
snapshot run and18 sequential benchmark runs (3 repetitions,2 meshes,3 modes).
Fresh standalone engineering/physical suites took13.83/28.05seconds on this laptop;
raw runs occupied1.2GB. The compact exported gallery is about9MB. Build time is extra.

## What to inspect

The HTML gallery and explained36-page PDF (28 figure pages) include native128×4 density, pressure andP/rho maps;
fluid profiles against independent Rankine–Hugoniot states; Eulerian front position;
space/time membership; ten sets of sampled parent histories; three strength/state
histories; event counts and Mach errors for all21 fixtures; resolution/CFL timing
brackets; uniform-flow transport statistics; negative controls; all-phase restart
comparison; actual wall time and diagnostic output sizes/counts.

CSV data, original native2D snapshots, exact inputs/logs, JSON test measurements,
executable/core-source identity and SHA256 file digests are retained with the plots.
Full raw trace/transfer directories remain in `--workdir`; they can be inspected
with the existing strict `vis/python/shock_history.py` reader.

| Test evidence | Local result | Boundary of the claim |
| --- | --- | --- |
| CFG01/06, PHY08/09, SYS02 | On/off output exact; uniform/contact/empty emit0events | Same-partition Serial CPU |
| TRC01, TRC10/11, INT06/08 | Uniform moments and source/timestep evidence; sampling aligned | Statistical cell-hop tracer, no continuous material paths |
| EVT/RST01/07 | All three checkpoint phases replay exact observations and128events | Matching source, native ABI, rank count and partition |
| RST03/04, IO03/04/08 | Expected damaged-history/checkpoint and overflow failures verified | No storage fault injection or distributed runtime claim |
| Held-out physical cases | All21 run; results recorded without acceptance thresholds | Physical/scientific qualification remains pending |

Known outcomes remain visible: all eight Mach2 seeds produce128 complete episodes
for128 unique particles; Mach5 produces134 for128 (six additional episodes);
boosted Mach2 produces132 for127 plus one suspended; Mach1.5 produces127 plus one
suspended. First-inside endpoint strength is not plateau-certified. Recorded time
brackets describe adjacent mask observations, not exact material crossing times.

## Measured laptop costs

Release, Serial CPU,1024 particles, uniform flow,tlim0.1, three rotated-order trials.
Whole-process medians include startup, initial/final outputs and durable empty
per-generation event chunks. Other fluid/particle settings are identical across modes.

| Mesh | Tracking off | Tracking on, trace off | Tracking on, trace on |
| --- | ---: | ---: | ---: |
| 64×4 | 0.01163s | 0.04171s (3.59×) | 0.28880s (24.84×) |
| 256×4 | 0.04174s | 0.12862s (3.08×) | 1.10772s (26.54×) |

Full traces are expensive qualification diagnostics. Output-cost plots report actual
bytes/file counts; resident memory was not measured. These short local timings do
not qualify sustained production performance, MPI, GPUs, RUSTY or Frontier.

The first plotting attempt failed in Matplotlib when laying out a tiny all-zero
restart difference axis. The final display expresses differences in1e−14code units
and places labels in axis coordinates. Exact equality assertions were unchanged.
Partial first-attempt artifacts are preserved locally under
`/private/tmp/shock-mc-first-plot-attempt` and `first-attempt-snapshots` in the raw run root.
Scientific sign-off and broader hardware/performance gates remain pending.

The documented one-command workflow was independently run without `--reuse`
in a fresh directory. It passed all engineering checks, completed 21 physical
fixtures and 28 plots, and repeated 18 benchmarks. All 21 physical measurement
dictionaries and 27 PNGs match the first gallery exactly, including output cost.
Runtime variation was retained; reproduction reports are under `clean-reproduction/`.

## Explained PDF report

The gallery PDF now includes eight explanatory pages and reading notes beneath all
28 original vector plots. The original plot-only PDF is preserved as
`evidence/shocks-mc/laptop-visuals/plots-only.pdf`. The explained report is also saved
in `output/pdf/athenak-shock-test-report.pdf`. No simulator source or results changed.

To author explanations for a new plot packet, use the bundled document Python:

```sh
/Users/hitesh/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 \
  tst/shocks-mc/explained_pdf.py \
  --packet /tmp/athenak-laptop-plots-new \
  --plots /tmp/athenak-laptop-plots-new/all-plots.pdf \
  --output /tmp/athenak-laptop-explained.pdf
```

This separate documentation step uses already bundled ReportLab/pypdf; it consumes
recorded results without rerunning AthenaK. The local test environment needs no
new package. Rendering/layout verification is in evidence/shocks-mc/explained-report/.
