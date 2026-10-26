# Curve output name and explicit encounter fields — 2026-10-26

## Accepted request

Use `file_type=curve`, require shock or Mach fields in each explicit complete-history
block, and explain existing THP `source`. Preserve `.crv` binary schemas, detector
thresholds, transport, event timing and restart algorithms. Retain `crv` as an input
compatibility alias and existing text `curves` output. Source defaults to particles
in the THP writer; no source is needed for curve output. `variables` is the existing
field-list key. Ordinary flow Mach and shock-strength estimates keep distinct meanings.

## Changes and justification

- Extend the existing native dispatch conditions through output construction,
  snapshots, encounters and startup summaries. No new writer or format.
- Enforce field opt-in in shared `shocks::HistoryOutput`, before either population
  enables history. Reuse `ParseTracerFieldList` so aliases/case/validation agree
  with snapshots. Each active block needs explicit `variables`; another block's
  fields or tracking defaults do not qualify. Disabled blocks still do no work.
- Preserve event-only fixed complete evidence; selected fields are the opt-in gate,
  not a lossy event-schema filter. Keep MC complete-history restriction.
- Extend native SM22 with expected rejections, Mach alias/case acceptance, old-name
  history parity and new-name snapshot clocks. Existing mode/restart/AMR checks
  remain. Update runnable plots to show these actual outcomes beside explanation.

## Verification pending

Base: 932c88e5ef3ff9c25ec060eff95bd0981200f98d on codex/shock-lagrangian.
Worktree: .worktrees/shock-lagrangian. Clean before edits. Primary feature_dump and
codex/shock-multiphase remain untouched. Previous sealed evidence is not modified.
Serial/double local qualification only; no new MPI runtime/GPU/performance claim.


Source provenance checked live: pre-feature `multiphase_turb` THP has no source
selector. Commit b0d38d28 introduced it for shared curve/fluid-tracer thermo output.
Its particles default keeps existing particle-only inputs working. Canonical curve
examples no longer need it; existing source=curves|both compatibility remains.

A local documentation/formatting helper first failed Python syntax parsing before
any writes; corrected immediately. This did not execute or invalidate native tests.


## Results

22 native groups, 47 binary encounter reader checks and 28 unchanged curve snapshot
reader checks pass: 97 local groups. SM22 contains nine native acceptance/rejection
trials. Every native executable digest is checked against the final actual binary;
commands/inputs/restart hashes and compiled source fingerprints are in native-report.json.

- New `curve` snapshots produce 23 saved clocks and 23 curve stencil batches,
  zero particle cache builds and no complete-history directory.
- `curve`/`crv`, shock and mixed-case Mach aliases preserve exactly 15 particle and
  16 curve completed events in the planar reference. THP source defaults correctly.
- Both populations reject missing variables, density-only and unknown `shockless`
  tokens; shock fields in the other snapshot block do not satisfy encounter opt-in.
- Unchanged 21 groups pass: sampled/no-work frequency, exact particle THP bytes,
  native mesh sampling, active/zero-step/pending/overflow/prefix restarts, live AMR
  births/negative IDs, short forced-box coexistence and independent thermo oracle.
- Independent accelerating thermo maximum error remains 3.552713678800501e-15;
  original 2e-12 budget unchanged. No detector/physics/timing tolerance changed.
- Ten plots regenerated from these actual runs, with explanations immediately nearby.

## Retained failure / correction

The first snapshot-reader command pointed --historical at CRV planar/curves,
which has no old THP and raised StopIteration after its preceding mutations.
The partially populated curve-input-snapshot-reader-20261026 directory is retained.
Corrected --historical to thp-qualified-CT-20261025/planar/curve_thermo_history and
used a fresh output directory; all28 checks then passed without source changes.
No test tolerance or expectation was relaxed.

## Limits / next actions

Local Serial/double only. Historical MPI compilation remains prior evidence;
this input cleanup has no new MPI/GPU/RUSTY/Frontier runtime qualification.
Scientific sign-off pending. No push. Original branches/evidence preserved.


## Committed source and sealed gallery

Code commit: `2caf18a09fe1bd0f7fdcea601b41db308b6ddb8f`. Gallery: `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-output-input-20261026/index.html`.
Verified 3618 packaged file hashes and 48 native-run executable digests.
Manifest SHA256: `c5794f8b630669fc07ad98ee8d7f170ea1e68b0eca0423a230cf9854ee90a0bb`.
Includes exact source archive, two actual local executable/configurations, all raw
native runs, reader mutations (including failed fixture setup), compact reference
data and historical THP fixtures. Previous galleries remain unchanged. No push.
