# multiphase_turb SGS port evidence — 2026-10-08

See [the port report](../../multiphase-port.md) for changes, measured results and
qualification limits. Code commit: `32c03c1a581974e5edb76e7da2f2eed29ff0f63d`.
Base: `526a187f938e38634e89a88192ebe9731c4d36bb`.

- `final-release-evidence.tar.gz`: all final Release inputs, logs, history,
  independent native checkpoints, thresholds, results and source/executable manifests.
  Also includes first/preformat run root summaries and manifests.
- `final-release-summary.json`: counts, worst errors, temporal orders, actual forcing
  work, pilot trends and restart equality metrics.
- `baseline-and-debug.tar.gz`: untouched baseline attempts and corrected results,
  independent built-in port comparison and bounds Debug coupled packets. Includes
  native fields/checkpoints, harnesses and manifests where produced.
- `source.tar.gz`: committed repository source (gitlinks do not embed dependencies).
  Obtain Kokkos at the exact gitlink `08ceff92bcf3a828844480bc1e6137eb74028517`.
- `port.patch`: complete local port commit including metadata.
- `commit-provenance.json`: final source fingerprint to local commit binding.
- `*-CMakeCache.txt` and build/oracle logs: compilation and configuration evidence.
- `multiphase-sgs-added-style.json`, style log and checker: added-line lint and
  Python syntax checks. Full changed-file lint retains baseline warnings.

The final Release archive has 663 independently verified members and is 25,476,478
bytes. SHA256:
`261ae56b586469a68dafba4445feff4c15f789f6e91b61420bbaf245a0ddc852`.
Final verification executable SHA256:
`0b4ac3cfaba15ef6cb6f6a8000bc07ffa104f1a10fd07a6557267b179ab8aed2`.
The executable is not archived; build configuration and committed source are.

Run manifests retain original absolute execution paths. Temporary paths identify
where execution occurred, not required permanent dependencies. Extract archives
into a fresh directory and use the recorded inputs with your rebuilt executable.
Native restarts/readers are serial, unrefined, native-double HD only. Binary layout
is platform-dependent; cross-platform restart portability is not claimed.

The first baseline restart/Sod setup mistakes and first port output-column reader
mistake remain in their original attempt directories. These are failed harness
attempts, not relabelled physical results. Original 2026-10-08 SGS evidence on the
other base is separate and unchanged. No new MPI/GPU or performance claim is made.
