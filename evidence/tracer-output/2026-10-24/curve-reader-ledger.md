# Curve reader and plotting script — evidence ledger

2026-10-24. Worktree `codex/shock-lagrangian`, starting HEAD
`22c36dfc7857b2e7c60e84ad81516e4b3e8f0772`. The prior compatibility audit was
already present as an untracked Markdown record; preserved and documented with
this increment. No native simulation was run and no existing evidence bundle
was modified.

## Changes and rationale

- Add `vis/python/plot_curves.py`, using existing parsers rather than another
  schema implementation. It selects exact curve/node IDs and draws stored 2D/3D
  geometry, the latest thermo profile, a node's thermo history and event brackets.
- Add optional `include_particles=False` to the existing `read_run`; default
  behavior stays unchanged. A curve-focused plot should not load or depend on
  unrelated particle files. CR06 verifies both opt-out and default strictness.
- Preserve independent clocks and actual birth times; do not interpolate a
  missing geometry stream or fill a newborn's prehistory. Handle geometry-only
  and thermo-only records with explicit unavailable panels.
- Reconstruct closure using stored integer winding and domain lengths. Keep
  all IDs as Python integers; the CLI checks IDs above \(2^{63}\) without rounding.
- Refuse nonexistent fields/IDs, ambiguous curve selection and figure overwrite;
  retain partial events, numeric flags and separate terminal views.
- Existing in-memory reader reuse is intentional. Production streaming/memory
  qualification is not claimed. No C++, format, detection, pusher or restart
  source changes; no new dependency was installed.

## Checks and exact commands

Run from `/Users/hitesh/hitesh/git/athenak_nyu/.worktrees/shock-lagrangian`:

```sh
MPLCONFIGDIR=/private/tmp/athenak-lagrangian-mpl /Users/hitesh/hitesh/git/athenak_nyu/vis/.venv/bin/python tst/tracer-output/curve_reader_check.py --runs build-lagrangian/evidence/native-output-qualified-final-20261023 --three-d-run build-lagrangian/evidence/native-output-frozen-CT/forced_both --output build-lagrangian/evidence/curve-reader-20261024-final
python3 -m py_compile vis/python/tracer_history.py vis/python/plot_curves.py tst/tracer-output/curve_reader_check.py
python3 tst/shocks-curves/io_reader.py
git diff --check
```

All 11 CR groups pass, plus the standalone existing curve IO/schema/adversarial
reader check. `curve-reader-report.json` records source hashes, selected exact
IDs, source run paths and plotted results. Native/legacy selection has 108 saved
node observations and 16 latest-profile nodes, identical across formats. The CLI
retains node ID 18446744069414584321 and curve ID 4294967295 exactly.

The AMR example selects actual inside-born node 30064771105. Its first saved
thermo sample equals its positive actual birth; all later samples follow birth,
and its emitted episode retains incomplete-entry evidence. The 3D figure reads
the previously saved forced-box run; no new turbulence simulation or physical
qualification is implied. A hand-built winding example tests closure in all
three axes separately from the actual saved-flow plots.

Checks also cover geometry-only/thermo-only inputs, corrupt geometry rejection,
invalid selections, unrelated bad particle files and preservation of an existing
figure's exact bytes. No tolerances were changed. The first development check
directory also passed all 11 groups and remains retained; the final directory
reruns after the plotting-label cleanup with source hashes matching the final code.

## Artifacts and remaining limits

Seven actual-data figures are retained in the final check directory and copied
to the reviewable gallery at
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/curve-reader-20261024`.
Descriptions sit beside each plot. Source/report hashes and gallery verification
are recorded in `curve-reader-artifact.json`; the gallery is sealed.

These are reader/plotting checks on existing local datasets. They do not rerun
the previous 82 native groups, qualify MPI/GPU/cluster execution or supply human
scientific approval. Full traces remain available through the supplied trace
reader; this script plots saved event evidence rather than reconstructing raw
membership between event observations.

## Final handoff

Implementation commit: `52cee8d1d46ffe445a95dfee364ec2fc0d3a20bb`. Gallery
manifest SHA256: `5b3a4b7048f217f1f2832069b8fc85b9ffef39bb50e24e3bb74829e23937f970`.
Verified 3,257 original native input files against the previously sealed native
qualification bundle, all gallery source/figure/report hashes, and the retained
check-fixture archive contents. Seven actual-data figures are present. No native
simulation was launched and no sealed prior artifact was overwritten.

The new script, reader option, checks and usage examples are committed locally
on `codex/shock-lagrangian`; no push was performed.
