# SGS reviewer suggestions: implementation evidence

Code commit: `ac384d2d4f52c3dd3b6cb9201905f02bf69adad5` on `codex/multiphase-schmidt-federrath-hd`.
Previous commit: `6e5e2d4951126b21e03f4fe837319c7c1ee85770`.
The [port record](../../multiphase-port.md) describes the implementation and limits.

Applied the reviewers' four small simplifications: reuse the existing turbulence
fatal-error helper, collapse one nested guard, remove a redundant final history
capacity check, and share the identical input serializer across three runners.
All pre-write history capacity checks remain. Fatal messages now use the existing
helper's standard prefix and output stream. No equations, coefficients, timestep
rules, numerical task order, supported configurations or acceptance thresholds changed.

The configuration runner also records each attempt before launch and saves it on
completion or failure. Subprocess timeouts, missing states and unreadable output
are recorded as failures while later checks continue. Partial process logs survive.
A standalone standard-library regression checks both subprocess call sites,
malformed checkpoint handling, preserved earlier results and continued later results.
It passes with this fix and fails with the previous runner because no results file
exists before the first launch. Its intentionally injected failures are not failed
physical qualification cases.

Release passes 46 unforced records, 32 configuration checks, two acoustic grids
with all eight acceptance checks, 19 forcing/restart records, the actual forcing
oracle and closure CTest. Bounds Debug passes the same unforced/configuration/wave
checks, actual forcing oracle and closure. The full forcing packet ran in Release.
The built-in comparison passes six bitwise SGS-off checks against the unchanged
multiphase base and two generic SGS seed checks.

- `reviewer-fixes.patch`: exact committed implementation and regression.
- `summary.json`: counts, executable hashes and source/harness-to-run verification.
- `qualification.tar.gz`: inputs, logs, history, results, manifests, build caches,
  regression logs, old-runner reproducer, style checks and baseline comparison.
- `before-after-checkpoints.json`: all 779 directly paired checkpoints match the
  previous cleanup byte for byte.
- `checkpoint-archive-index.json`: all 795 new packet checkpoint files have bytes
  already retained in the original port's `final-release-evidence.tar.gz`, verified
  by SHA256. No duplicate raw checkpoints are stored in this smaller archive.

Run manifests retain the prior committed revision plus the execution-time diff.
Production and harness fingerprints match the final commit. They are not rewritten
as if execution began after the commit. Added-line C++ checks and Python parsing pass;
whole-file lint still has baseline diagnostics.

Two development checks needed correction before qualification: the first regression
mock also intercepted Git provenance reads, producing a missing-results-file error;
it now forwards Git reads to the real subprocess runner. One helper call exceeded
90 columns and was wrapped. The old-runner failure retained in the archive is a
separate deliberate regression reproduction. No acceptance threshold was relaxed.
There is no added MPI/GPU, AMR, cooling or scientific turbulence-accuracy qualification.
