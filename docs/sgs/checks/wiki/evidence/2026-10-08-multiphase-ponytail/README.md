# Minimal-intrusion SGS cleanup evidence

Code commit: `6e5e2d4951126b21e03f4fe837319c7c1ee85770`.
Previous port: `32c03c1a581974e5edb76e7da2f2eed29ff0f63d`.
The [port record](../../multiphase-port.md) explains the architectural review and
changes. Only the existing SGS branch was edited; unrelated checkout/wiki work
was preserved.

- `summary.json`: packet counts, executable hashes, source-to-run fingerprint
  verification, unchanged hydro task graph and unchanged AMR source check.
- `cleanup.patch`: exact committed simplification and output compatibility test.
- `qualification.tar.gz`: inputs, logs, history, results, manifests, build caches,
  closure/Ctest/style logs, baseline comparison harness and twelve-scalar binary
  outputs demonstrating that SGS columns are appended.
- `before-after-checkpoints.json`: SHA256 for all 395 native checkpoints, all
  byte-identical to their pre-cleanup counterparts.

The native checkpoints are retained in
[the original port evidence](../2026-10-08-multiphase/README.md), including its
`final-release-evidence.tar.gz`; duplicating identical large states is unnecessary.
The temporary runs contain the same checkpoints. Run manifests retain the revision
and working diff from execution time, with per-file hashes checked against the
follow-up commit; they are not rewritten to imply execution began after the commit.

Intermediate cleanup checks preceded the output-order correction; their results
are retained in the archive. The final-* packets include the output correction.
All final packet checks pass; no thresholds were changed. There is no new AMR,
cooling, MHD, MPI/GPU or scientific turbulence-accuracy qualification.
