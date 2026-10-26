# Curve THP topology and binary geometry — 2026-10-25

## Initial scope (superseded by unified CRV below)

User requests all missing curve geometry in curve THP, excluding current length,
stretch, curvature and mesh-block/AMR level. Separate binary .crv geometry includes
those diagnostics. Interpret "strentth" as stretch, matching the prior question.
Preserve ordinary particle THP v2, old geometry/CSV readers, shock encounters and
restart state. No dynamics, detector thresholds or physical tolerances change.

## Initial interface (superseded before completion)

- Curve THP v3: same named thermo fields and THP v2 sample block, fourth integer
  reserved zero, metadata birth/a; no leaf GID or level. Prefix embeds bounded,
  checksummed typed geometry metadata without excluded diagnostics.
- Shared typed geometry payload v1: signed32 IDs/links/flags, signed64 winding/images,
  uint64 counters/allocator, declared native endian/Real width, clocks/rank/split,
  domain/dimensions/BC/config, complete per-curve/per-node topology and provenance.
- file_type=crv: same payload with geometry diagnostics and current leaf GID/level,
  rank-local immutable chunks through the existing output scheduler. Existing
  file_type=curves ASCII and restart/checkpoint representation retained.
- Existing readers/CLI discover new files and retain historical support. Legacy
  particle reader explicitly rejects curve v3 rather than misreading IDs.

## Gates

- [x] Inspect clean branch/source, callers, format readers and scheduler.
- [x] Unified native implementation and reader strict validation.
- [x] Native parity, insertion, periodic geometry, AMR, signed limits, empty ranks,
  restart/replay and ordinary particle compatibility; retain failed cases.
- [x] Actual plots near explanations and independent persistence/timing review.
- [x] Exact source/config/commands/outcomes and coherent commits.

Deployment MPI runtime/GPU/cluster budgets and human scientific approval stay pending.

## Superseding accepted scope — unified CRV

User stopped the split THP-v3/geometry-only implementation before completion, then
requested one .crv carrying geometry/topology, all diagnostics, mesh location and
selected thermo/shock fields. Subsequently clarified: do NOT embed complete
encounter history. Complete event journals and restart persistence stay separate.
No cursor or event-bundling changes were implemented. The unfinished compile
failure (missing store.hpp include) is retained in build-dev.log.

Current implementation removes curve THP writing; source=curves/both writes .crv
using the same cadence and selected fields. file_type=crv is the explicit curve
output option; variables default to native default fluid fields plus shock when
tracking is enabled. Particle THP v2 bytes remain unchanged. Unified binary output
suppresses automatic initial/final/periodic ASCII geometry, while explicit
file_type=curves ASCII and checkpoint .curves files remain available.

Independent review found/fixed missing crv dispatch exclusions and reader
cross-clock fixed-metadata consistency and stem-selection issues. A native
three-step planar smoke run passed: four frames, 64 sampled nodes, full default
fluid and six shock fields. Reader/testing qualification is ongoing.

## Final engineering gate

50 native groups passed on the final source-bound Serial/double executables:
CRV01–16, LT14 and CSH20. Another28 reader groups pass locally and independently.
Detector device self-checks (8) and EVT01–08/11–14 encounter checks pass; Python
compilation and git diff --check pass. All source-bound run executable hashes
were checked against the final binaries. Built-in/prescribed and MPI compile
builds passed. No runtime MPI retry or GPU/cluster qualification is claimed.

Exact-value parity against prior accepted output: planar1,728 samples, live-AMR
negative-ID9,644 samples and accepted short32³ forced box256 samples. All old
shared thermo/shock values, geometry/diagnostics, particle THP bytes, encounters
and native fluid/mesh arrays match. BIN embedded input text changes with output
settings; the harness correctly compares arrays/clocks rather than parameter text.

Independent physical oracles retain original budgets: accelerating thermo max
error5.33e-15 and geometry reconstruction max9.95e-14 within CT2e-12; accepted
circle production budgets are2e-11 length and2e-10 curvature. Periodic translation
uses production.py translation_wrap values. The earlier ad-hoc circle/wrap smoke
combination is retained as development evidence, not the final physical fixture.

Retained harness failures: NaN direct equality, checkpoint filename digit-count
assumption, native BIN parameter-dump byte comparison, and an import-path error.
All corrected without changing physics or weakening numerical tolerances. Plotting
also exposed that initial checkpoint13 preceded entry. Final CRV08 derives the
actual first inside cycle37, asserts trace phase=active, and verifies exact active
continuation plus zero-step geometry/thermo/particle records. Earlier equality-only
cycle13 checks must not be cited as active-encounter evidence.

Nine source-bound plots use actual CRV records and the supplied plotting API.
Matched record storage: planar740099→454788 bytes,216→108 files; AMR
3956634→2318246 bytes,324→162 files; forced108983→67472 bytes,32→16 files.
These are measured file sizes/counts, not runtime/memory/cluster benchmarks.

Reader regressions include repaired-checksum semantic corruptions, exact signed
-1/full-width allocator preservation and the intended missing-rank-set rejection.
Six representative binary fixtures and expected outcomes are retained here.

Complete encounter journals/checkpoint state remain independent by explicit user
choice. No event cursor, event buffering or encounter/restart-schema change added.
No push requested/performed. RUSTY/Frontier/GPU production budgets and scientific
review remain pending. Evidence seal and coherent commit follow.

## Evidence seal

Source commit: `9816a84bf782b3c143e9e71d3e751bbddcc53ae0`. Nine explained figures, exact source and
Kokkos archives, three executables/configurations and reports are sealed in
`/Users/hitesh/.codex/visualizations/2026/10/08/01a11cbd-ecd8-7103-86ce-b245d30c9407/unified-crv-20261025`.
Verified28,317 raw-file hashes and89 native-run executable digests. Manifest SHA256:
`37e9a48beac28d7db438a9c3576856b98579b3825a377bcadaae8f9436a009b6`. The packaged manifest hash was rechecked after sealing.
Deployment and human scientific gates above remain pending.
