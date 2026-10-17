# Larger longer forced-box evidence - 2026-10-14

Actual local Serial CPU/double run:128^3,8192MC particles, initial Mach5,
continuous native edot50, gamma1.4/RK2/PLM/HLLC/CFL.30, t=1.0.
Compared with64^3/t.3:8times cells,3.333times code-time span.
Width8/nghost8 keeps physical detector half-span1/16box (earlier width4/64).
No production source changed. Every-step raw shock flags and phases saved.

Run exit0; generation6033; wall126.457minutes.
Complete367,074, partial0, open999;
complete/particle44.808838. Native sampled
49,430,528, invalid0; all8192IDs every frame.
Final volume rmsMach0.808576, supersonic volume
20.399percent. Adiabatic heating changes
sound speed; this is an evolving forced realization, not stationary Mach5.

Exact width8 NUM08 parity/NUM10 independent planarMach3 oracle retain1e-12
strength tolerance and exact masks. Strict bounded history parsing verifies
canonical identities, sequences, counters, brackets and censored records.
Five-step native same-partition quarter restart matches40,960trace rows and
all1449canonical chunks byte for byte. An independent80-step replay crosses
the scheduled forcing update: all1524canonical chunks and655,360continued
trace rows are byte-identical, with no restored-step reobservation. Terminal
censored records are checked separately, not compared with nonexistent views.
The final independent audit recomputes episode counts/timing/flags/strength
summaries and derives final Mach/energy from conserved fields, agreeing with
the primitive route at unchanged1e-12 tolerance (fractions exact).
Checkpoint recovery evolves zero steps. Full8192thermo samples agree with
original traces; native plotted z64 slices match actual live particles' masks,
strengths, normals and endpoint evidence at1e-12. Crop support applies only
to the central plane, not whole-volume mask reconstruction.

Original dcycle=-1 key presence suppressed intermediate volume output; this
frozen failing cadence configuration remains in input.athinput. Same-binary
paired checkpoints supply actual quarter-time maps. Corrected portable input
omits those keys; smoke produces0/1/2, with an exact original two-step sample
prefix check. Initial CLI nlim attempt failed before evolution; successful
variant explicitly declares nlim. SystemPython lackedNumPy on two helper
attempt; the installed scientific environment executes the actual checks.

Scientific approval is pending. The inherited RK2 source thermal oracle
remains a failure; plateau/physical strength truth, stationary ensembles,
physical convergence tolerance, MPI/GPU/RUSTY/Frontier budgets remain open.
Numerical elapsed spans retain entry/exit bounds; no exact crossing time or
continuous MC trajectory is invented. All flags and open identities remain.

## Inspect and reproduce

The PDF has explanations adjacent to each of7vector figure families. PNG/SVG/
PDF plots, episodes.csv, evolution.csv, selected-traces.csv and frames.jsonl
are included. Complete per-step all-particle histories and paired checkpoints
are preserved at `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/highres128-long-20261014/raw` (outsideGit). Exact executable, detector fixture binary and build cache are retained in `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/highres128-long-20261014/runtime`. The original `/private/tmp`
path is a symlink to this durable copy, keeping executed metadata valid.
`raw-sha256.json` binds every primary raw file; checkpoint-recoveryJSON binds
the input pair for every recovered map. The packet intentionally omits large
raw cubes/traces and native sliceCSV fixtures; these remain in the raw/plot
artifact directories. Source-verificationJSON gives the current absolute paths.

Fresh run, in an empty directory:

```sh
ATHENA=/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu/build-release/src/athena
INPUT=/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu/inputs/tests/shock_mc_highres_long.athinput
"$ATHENA" -i "$INPUT" shock_tracking/directory=/absolute/new/events
```

Use a new directory and a distinct lineage when combining separate realizations.
Allow tens ofGiB for fulltrace=true. All native restart files must remain paired
with their .shock.rankN sidecars and canonical event prefix. Samebinary and
qualified samepartition only; the feature restart contract rejects incompatible
state. Both quarter-replay/replay_check.py and quarter-replay80/replay_check.py
record concrete restart commands and their exact-source verification reports.

Analysis (existingNumPy/Matplotlib scientific environment):

```sh
/Users/hitesh/hitesh/git/shockit/.venv/bin/python /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu/tst/shocks-mc/highres_long_analysis.py \
  --root /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu --run /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/highres128-long-20261014/raw \
  --trace /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/highres128-long-20261014/trace-validation --baseline /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu/evidence/shocks-mc/supersonic-showcase-20261013 \
  --output /absolute/new/plots --pdf
```

The frozen executed-trace-validator.py provides bounded fulltrace collection
for a fresh equivalent run; it takes --run and --output. Original input/config
commands, all logs/scripts, independent reviews/checks and frozen source hashes
are retained. BundledReportLab/pypdf Python runs highres_long_pdf.py with the
plot packet and --root/Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu; PDF pages were parsed/rendered/visually checked.
New source/evidence commits are local, not pushed.
