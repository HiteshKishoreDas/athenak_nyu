# Remove output source selector — 2026-10-26

## Accepted scope

User explicitly requested removing source= after confirming it was redundant.
Particle THP must select regular fluid particles; curve CRV selects curve nodes.
Old source=curves|both inputs must fail clearly instead of dropping curve data.
This cleanup also rejects explicit source=particles; omit the removed option.
Physical source terms and unrelated variables named source are outside scope.

## Steps / rationale

Base commit998c52ef on codex/shock-lagrangian, clean worktree before edits.
Traced the THP constructor, data sampling/write paths, member state, early shock
activation, output dispatch, startup summary and current docs before editing.
Deleted carrier-selector parsing, booleans, curve buffer, conditional gathering,
curve write branch and early source defaults. Kept particle THP writer/byte layout,
shared field parser/evaluators, curve writer, detector, transport and restart engines.
One active-output boundary guard rejects the removed setting for both populations
in snapshot/encounter modes. Startup output no longer references source.

Reference fixtures used source=both; the runnable test now explicitly splits that
historical input into particle output2 and curve output10, preserving id/fields/
cadence and original physical inputs. It strips only old output-block selectors
before applying test options; negative tests add those selectors back deliberately.
No accepted reference data is rewritten. Added SM23:12 actual native rejections
across three old values, two populations and two modes, plus startup-summary check.
Existing SM01-22 and reader tests remain. Added an explained plot of SM23 outcomes.
Updated current docs and marked older source-based configuration examples historical.

## Pending checks

Local Serial/double main/prescribed builds and mode/reader regressions pending.
No new MPI/GPU/performance/scientific qualification. Original branches preserved.


## Results

98 local groups pass:23 native +47 encounter-reader +28 curve-snapshot reader.
All12 SM23 native cases reject old source=particles/curves/both on active THP/curve
snapshot and encounter outputs with the specific migration message. Startup source
line is absent. No unexpected failing fixture or tolerance change in this gate.
The suite retains every expected negative fixture and raw run output.

Particle THP v2 bytes match accepted output exactly after splitting the previous
source=both fixture. Curve sample/evidence fields and native schedules are exact.
Planar complete histories still15particle/16curve events; sample-only23clock calls.
All bounded pending, overflow, active/zero-step and exact-prefix restart groups pass.
Live AMR new-node/negative-ID parity and short32³ forced-box6904particle/28curve
encounters remain exact. Independent thermo maximum error3.552713678800501e-15,
original2e-12 limit unchanged. Scientific approval remains pending.

Two actual Serial/double binaries rebuilt. report.json records all changed source
hashes, actual binary/cache hashes and report hashes. Native report records exact
inputs, commands, expected errors and compiled source fingerprints. Every reported
native executable digest verified against the final actual binary. Eleven plots
regenerated from these runs, with adjacent explanations. No push or branch deletion.


## Sealed evidence

Exact code commit `d0359aab90b12d69556457c4313708175b5c6141`. Gallery: `/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/output-source-removed-20261026/index.html`.
Eleven explained plots, exact source archive, actual local executable/configuration
pairs, raw native runs and reader/reference fixtures. Verified 3616 packaged
file hashes and 60 native-run binary digests. Manifest SHA256:
`063c83b2ae4fc22d6a8512414e8db80230fbf5b5c94ae2d74c4567dc39863664`. Primary feature_dump remains clean;
codex/shock-multiphase remains at2450e1c21b0a97c0f1635da75a20701c7e0852a4.
No push; hardware/performance/scientific gates remain pending.
