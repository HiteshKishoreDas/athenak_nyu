# Particle workflow compatibility audit — 2026-10-24

Read-only code/report audit of output integration through `22c36dfc`; no simulation rerun.

- Legacy format remains default; existing particle THP v2 is byte equal in TO01.
- TO17 fixed-mesh MC thermo (shocks off) passes. TO05/TO21 same-source restart continuation is exact, including live-AMR insertion in TO21.
- The output increment does not modify particle pushers/tasks, native MC transport, trk/pvtk writers or restart.cpp (`git diff 66e436e1..a8e19214` on those paths is empty).
- Opting into binary changes shock/curve filenames/content to .tbin; scripts hard-coded to CSV/text need the supplied readers or legacy mode.
- Checkpoint state layouts are unchanged, but tracked restart contracts bind source fingerprint and format. An older build's shock checkpoint is not a guaranteed cross-build restart; retain its executable. Changing legacy/binary during a tracked restart is rejected.
- Full drift/cosmic-ray workflows, MPI/GPU/cluster execution and production scaling were not comprehensively requalified by these local tracer checks. Unchanged source is not runtime qualification.

Conclusion: no observed regression for the tested existing fluid-tracer workflows; no blanket claim of backward compatibility for every particle workflow or historical checkpoint.

Evidence: `../2026-10-23/{output-report.json,ledger.md}` and `../../../docs/tracer-output-unification.md`.
