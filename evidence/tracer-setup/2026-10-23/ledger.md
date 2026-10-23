# Startup-summary implementation — 2026-10-23

User authorized only printing requested/actual counts, placement rules,
adaptation settings and enabled histories. Pusher defaults and curve block
numbering/discovery explicitly excluded.

## Implementation and scope

Base: clean published codex/shock-lagrangian at 89338ee96b70e0d2e0da938d8cfc109a3ee6ec2b.
Read Driver::Initialize, particles seed parsing/selection, curve initialization
and adaptation, runtime output setup, applicable repository/skill guidance.

Only runtime source edit: 120 added lines in src/driver/driver.cpp, comprising a
read-only startup formatter, its includes and one call after initial populations
and histories are ready. Existing global mesh/rank ownership stays unchanged.
Global particle/curve/node/anchor/quality counts use one integer MPI reduction;
only rank zero formats output. No default-setting calls or particle/device data
copies are introduced. Initial requested particle count uses the native due-time
criterion and resolved initial seed settings. Curve anchors count original seed
vertices independently of inserted nodes. Restart is labeled restored, not reseeded.

Output reports type, global requested/actual counts, per-seed selection/time/filter
settings, placement rules, curve shapes/closure, adaptive limits, geometry output
cadence, shock event/trace status, and constructed thermo/tracking streams.
No new user configuration is needed. Neither excluded UX change was implemented.

## Verification

- [x] Configure/build existing local Serial/double prescribed and built-in targets.
- [x] SUM01: fresh counts; particle history and initial/final curve files byte-equal
  to preserved pre-change executable across four completed timesteps.
- [x] SUM02: initial adaptation, requested 8 versus actual 15 nodes independently
  checked against native geometry; automatic settings and quality displayed.
- [x] SUM03: no eligible cells, requested 128 versus actual 0; filter shown.
- [x] SUM04: future injection, requested 0 now / actual 0; 128 per future event shown.
- [x] SUM05–06: analytic closed circle and explicit lifted polyline placement.
- [x] SUM07: restart counts labeled restored; no reseeding; geometry byte-equal.
- [x] SUM08: enabled shock histories/traces and actual thermo cadence; disabled
  output blocks omitted from active histories.
- [x] SUM09: empty, curve-only, particle-only and absent populations.
- [x] Built-in hydro smoke: exactly one summary, 128 particles /16 nodes,
  particle+curve shock histories enabled.
- [x] Source diff whitespace check and native-data plots inspected.

Final runner: tst/tracer-setup/check.py, nine passing groups. Full raw input,
command, binary/source fingerprints and logs are in report.json and native.tar.gz
in the visualization packet. Existing numerical tolerances and physical fluid
setups were not changed. Counts in plots come from actual native startup summaries;
positions come from native particle/curve output.

## Retained failing fixture

First SUM08 enabled particle shock histories on the prescribed acceleration/source
fixture; native validation rejected it because arbitrary user source stages are
unqualified for particle shocks. Corrected that check to use the established
native RH shock-tube fixture through material_curve_test's existing shock_tube
mode. No guard was weakened. The first run/report/logs remain archived.

No full science suite was repeated for this logging-only edit. CPU/MPI aggregation
is implemented, but this run qualifies local Serial/double only; actual multi-rank
execution, GPU/cluster deployment, performance budgets and scientific approval
retain their separate outstanding gates.

## Reviewable artifacts

Viewer: /Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/tracer-startup-summary-20261023/index.html

The packet includes the actual summary for every case, figure, complete raw runs
(including the rejected fixture), runner and source snapshot/diff. Source/config
and archive hashes are recorded in manifest.json. Old sealed thermo/science
packets remain untouched.
