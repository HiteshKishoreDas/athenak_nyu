# Curve shock-history acceptance evidence

## Frozen before candidate runs — 2026-10-17

Owner: bounded independent fixtures only. Production numerical/IO source belongs
to coordinator and persistence reviewer. No scientific sign-off from these tests.
Commands, exact inputs, binary/source hashes, stdout and failures remain under the
selected `--workdir`; reports are Markdown/JSON, never replace prior failing runs.
The old merged binary is a red baseline, not a candidate feature qualification.

Accepted: global tracking controls both present consumers; event output always,
trace optional; every completed timestep and insertion sampled; IDs distinguish
MC and curves; nodes never inherit encounters; inside births left-censored;
fresh starts first; actual AMR required; laptop before RUSTY MPI.
Experimental unchanged detector: width 4, minimum Mach 1.1, unreduced membership,
existing encounter/rearm rules. AMR width units await the user's answer.

## Independent checks / strict criteria

- CSH01: smooth zero, translation and shear controls emit zero events; zero flow
  and affine prescribed translation/shear positions match independent formulas
  within the existing curve roundoff budget \(2\times10^{-12}\). No physical
  convergence claim from this roundoff check. Existing prescribed fixture values
  supply velocity, shear amplitude and runtime, not invented production values.
- CSH02: periodic planar forward-RH \(M=2,\gamma=1.4\), existing accepted states,
  front \(x_f=0.3+0.5t\). Central-front analysis stops at \(t=0.08\): maximal
  characteristic speed from either initial state is below \(3.05\), boundary
  influence reaches at most \(0.244\), while target front is \(0.34\). Do not
  analyze the unrelated periodic discontinuity. Native density-profile midpoint
  crossings bracket material passage independently of mask/event bookkeeping.
  Each initial target node must have one emitted completed passage, no duplicate
  identity; report strength errors without assigning scientific tolerances.
- CSH03: sampled entry/exit are observed membership brackets. Require observation
  fields/time/generation match actual curve snapshots, valid inside/outside
  chronology, no fabricated exact crossing times. Report density-front crossing
  time brackets separately; do not demand a midpoint equals a mask boundary.
- CSH04: initial and positive-time inserted nodes inside membership must have
  left-censored first episode and no valid pre-birth observation. Outside births
  start OUTSIDE without inherited sequence/events. Require positive-time births,
  both inside/outside classes actually exercised; missing coverage is not pass.
- CSH05: fixedseed MC+curves coexist; tracking off/on and traces off/on preserve
  node geometry, material IDs, native MC trajectories and fluid arrays exactly.
  The optional trace setting does not suppress events or restart state.
- CSH06: compare uninterrupted/restarted same-partition complete event identity,
  evidence, current node state and RNG trajectory; active/pending/outside phases
  must actually be observed before claiming phase coverage. Missing, corrupted,
  future-version sidecars and incompatible prefix/configuration reject explicitly.
- CSH07: existing timed refine/coarsen hook actually changes leaf count/levels;
  preserve IDs, births, events; smooth AMR remains zero-event; current node's
  leaf-local level/spacing agrees with native saved leaf geometry. AMR width
  qualification depends on accepted units; acceptance report must name units.
- CSH08: overflow rejects explicitly without silent drop. Unsupported old curve
  checkpoints cannot silently activate history. No MPI/GPU/RUSTY claim without
  those exact processes/hardware being run.

## Scope / source of physical values

Planar states reused from `inputs/tests/shock_mc.athinput`; periodic material
geometry already required by the accepted curve implementation. Smooth/AMR values
reuse `inputs/tests/material_curves/prescribed.athinput` and
`inputs/material_curves/amr-restart.athinput`. New geometry is a test stimulus,
not a recommended scientific setup. No floor or detector threshold tuned.

## Reproduction

```sh
python3 tst/shocks-curves/acceptance.py --baseline /path/to/old/athena \
  --workdir /private/tmp/curve-shock-red-unique
python3 tst/shocks-curves/acceptance.py --executable /path/to/built-in/athena \
  --prescribed-executable /path/to/material_curve_test/athena \
  --workdir /private/tmp/curve-shock-candidate-unique
```

Each process capped at 30 seconds; suite capped at 300 seconds. Never remove a
failed fixture or weaken these bounds to claim a pass. Human scientific review,
forced thermal-oracle issue and RUSTY resource-scale qualification stay pending.

## Red baseline and oracle observations — 2026-10-17

Command: reproduction baseline form above, executable merged built-in AthenaK
`feature-dump-merge-20261017/runtime/built-in/athena`, SHA256
`bd2038a751e42335e1ef3e860d827fd060c4e272cf42a2f78374665f43de2769`.
Evidence `/private/tmp/curve-shock-red-v3-20261017/report.json`, all native logs and
inputs retained. Both red cases exit 1: MC+curves rejected by passive HD guard;
curve-only tracking explicitly reports material-curve shock histories unsupported.
The third control disables tracking and runs the identical accepted periodic RH
fluid/curve setup: exit 0, 107 actual steps, 16 anchors retained. All 16 native
node paths cross the independently reconstructed central density midpoint once;
first bracket $[0.035409375502100174,0.03615526222668155]$, last
$[0.05187278068517021,0.05262134828197653]$. These are native density-front brackets,
not exact physical crossings and not shock-mask entry times. Profile timestamps
are printed at limited precision; join to native snapshots by cycle instead.
Historical first/second red attempts retained separately (`red`, `red-v2`).

Runtime dependencies: standard-library Python, existing strict MC output helper,
existing strict material snapshot reader and the new curve-history reader. No new
Python package or test framework. Parent owns production code/builds/review.

Runner refinements before candidate execution: retain NaNs using canonical JSON
when comparing native curvature-bearing snapshots (Python NaN inequality must not
cause fake passivity failures); join formatted density profiles by native cycle;
reconstruct sampled AMR leaf level and cell width from saved logical leaves;
future sidecar fixture uses finalized native `AKCSHCK1` version field at byte 8.
Checksums test accidental corruption; they are not authentication. Source hashes
read from the current worktree are explicitly unbound to the old binary; only
compiled fingerprints present in executable bytes count as binary/source binding.

## Independent AMR MC limitation reported during implementation

Migration reviewer exercised the native MC regrid transaction: 4096 distinct
particles at coarse cell center $(0.125,0.125)$ all mapped to fine child center
$(0.1875,0.1875)$. Child counts $[0,0,0,4096]$ differ from equal-child uniform mass
fractions $1/4$. Existing cell-center remapping preserves IDs/state but does not
establish mass-following AMR. This is reported reviewer evidence, not a run by
this acceptance runner. The live RH+AMR fixture can qualify identity, phase,
observation and persistence engineering only; scientific MC AMR qualification
remains open/failing. No RNG/remap algorithm or physical threshold changed here.

## Added forced-box gate CSH09 — frozen before forced candidate runs

Use the exact physical settings of `inputs/tests/shock_mc_forced.athinput`, including
its $32^3$ mesh, 4096 MC particles, native velocity/forcing seeds, forcing amplitude,
RK2 step settings and $t=0.2$ duration. Add only the existing contained curve seed
from the planar fixture. No new forcing physics, floor, strength or runtime value.
Three fresh runs: tracing on, tracing off, global tracking off. Require exact
native fluid-array, MC-trajectory and material-geometry equality, unchanged event
identities/evidence when traces disabled, synchronized actual node sampling and
valid bounded outputs. Report counts, censoring and runtime; no physical strength
truth claims. Reuse existing native binary reader; excludes text input headers
from fluid-array comparison. Inherited RK2 forcing thermal-oracle failure remains.

```sh
/Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python \
  tst/shocks-curves/acceptance.py --suite forced --executable /path/to/athena \
  --workdir /private/tmp/curve-shock-forced-unique
```

Initial candidate report `/private/tmp/curve-shocks-acceptance-initial/report.json`
retained. Two harness defects found from output: periodic folding changes the
last floating bit of positions; compare reconstructed lifts using the already
frozen $2\times10^{-12}$ coordinate roundoff budget, keeping exact time/birth.
Restart fixture had `dt=-1`, which suppresses restart files: seed/checkpoint-off
runs now explicitly use the existing `dt=1` stream, without changing physics.
Neither is a production timing finding or a loosened physical tolerance.

MC+AMR tests are blocked pending the user's remapping decision; native startup
rejection remains enabled. Do not reinterpret that rejection as an engineering
pass for the accepted physical AMR requirement.

Read-only corrected reanalysis of that initial candidate (`bb2a70215ed7a95995b91ce4c3e93ac40fd104961577fbe5f5382c5273d5d3eb`):
MC+curve and curve-only planar runs each retain 16 completed events, 16 nodes,
zero terminal views, zero extra episodes, 1728 synchronized node observations;
all 16 native-density crossings bracket uniquely. Initial-inside run has 16
partial left-censored events and 880 synchronized observations. These are
preliminary old-code measurements, not final-code qualification. No new process
was launched for the reanalysis, and original failed reports are unchanged.

CSH03 timing check explicitly requires each independent native-density crossing
bracket to overlap its event's observed pre-to-post evidence envelope. This
checks coverage, assigns no exact crossing time, uses no adjustable tolerance,
and does not equate a density midpoint with a mask entry/exit boundary.

Final-run hardening requested by independent review, before final rebuild:
full trace generations and exact node-ID sets must equal every native snapshot;
earliest observation across normal and birth traces must equal each actual
node's birth with sequence zero and correct inside censoring. Replayed MC
thermo must equal the entire expected post-checkpoint generation/ID/value set,
not an arbitrary subset. Compare full replayed traces and terminal censored
views as well as emitted events. Add `-d` continuation with a relative configured
history prefix, and authentic checkpoint rejection after changing the saved
`problem/refine_time` callback parameter. No tolerance changed. Read-only checks
of preliminary planar output cover all 16 independent density-front brackets;
maximum lifted-coordinate fold error is $6.94\times10^{-18}$, below the frozen
roundoff budget. Original initial report remains failed and unmodified.

Negative-path isolation: CSH08 overflow uses curve-only vertical 16-node line at
$x=0.4$ and capacity 1, so a simultaneous node passage exercises the curve event
buffer rather than accidentally passing via MC overflow. Old untracked restart
activation uses a curve-only checkpoint, so a missing MC sidecar cannot mask the
curve history guard. Geometry stays inside the frozen planar stimulus range.
