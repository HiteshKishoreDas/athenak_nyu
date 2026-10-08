# Independent baseline MC numerical/time/persistence audit

Status: read-only baseline review; no runtime qualification or scientific approval.
Base: `526a187f938e38634e89a88192ebe9731c4d36bb` (`multiphase_turb`), new branch
`codex/shock-multiphase`; pinned Kokkos `08ceff92bcf3a828844480bc1e6137eb74028517`.
The findings below refer to baseline Git objects; later integrated edits require
independent rereview. Only review evidence was written by this worker.

## What the existing method actually does

`particles_lagrangian_mc.cpp:500–583` samples a single outgoing face or no move
per completed fluid step. It divides saved signed mass-density transfers by the
**initial-step cell density** in `u1` and takes the positive outward contribution:

`p_xminus=max(-Q_xleft/rho_n,0)`, `p_xplus=max(Q_xright/rho_n,0)`, similarly y/z.
`Q_face=(dt/nstages) sum_stage(F_mass_face/dx_axis)` (`hydro_fluxes.cpp:370–392`).
One draw in [0,1) chooses cumulative categories x−,x+,y−,y+,z−,z+, with the
remaining probability representing no transfer. One chosen face changes one
coordinate by exactly one cell width. Particle positions are discrete cell
centers; the stored velocity slots1/3/5 are zero/unused, not a fluid path.
Physical cell mass is rho times volume; using Q and rho cancels volume on the
uniform Cartesian mesh. Fluxes are saved **before** flux-reflux communication
(`hydro_tasks.cpp:57–60`), so AMR is separately unqualified.

`PLASTMOVE` stores parent-cell parity bits32/16/8 and direction0–6;
`PLASTLEVEL` stores the parent level (`particles_lagrangian_mc.cpp:560–582`).
These fields are used by the subsequent coarse/fine adjustment. The RNG is
stateless SplitMix64, top53bits converted to Real, with key
`uint64(tag)*7919 + uint64(ncycle)*104729 + uint64(random_seed)` (lines39–50,
555–558). Rank, local array order and GID are absent from the push key. The
parent-cell probability and unchanged persistent tag determine each transfer.
AMR adjustment uses key+1 (693–696), not a new stored RNG state.

## Frozen field contract versus the previous velocity branch

No ideal Newtonian hydro pressure/axis drift was found:

- `athena.hpp:58–59`: rho0, vx/vy/vz1/2/3, primitive IEN4.
- `tracer_fields.cpp:175–181`: pressure=(gamma−1)*internal-energy density,
  temperature=pressure/rho. Detector must still use `(m,k,j,i)` coordinates and
  the physical dx/dy/dz with active dimensions.
- `athena.hpp:65–70`: MC position slots are **0/2/4**, creation time6;
  integer slots GID0/tag1/lastmove2/lastlevel3/seedid4. Seven real/five integer
  fields. Neither contiguous position slots0/1/2 nor drift velocity fields are
  valid substitutes. Migration must preserve all existing fields plus new state.
- Existing thermo-history `mach` is lab speed/sound speed, not shock Mach.
  Existing entropy evaluates log(P/rho^gamma), while the frozen detector uses
  logP−gamma*logrho; these are equivalent on ordinary states but have different
  extreme-value overflow behavior. Do not use thermo-history as an exact
  detector parity oracle or silently change the frozen detector.

## Correct observation phase

The actual completed-step order is:

1. Fluid stages: initial `u0` copied into `u1` at stage1; accumulated saved
   fluxes; RK update/source/user-source/boundaries/conserved-to-primitive;
   complete all after-stage communication (`hydro_tasks.cpp:55–78,131–146`).
2. `after_timeintegrator`: MC push → new GID → particle migration/receive/send
   completion → refinement adjustment → due seeding (`particles_tasks.cpp:31–46`).
3. Driver advances `mesh.time` and `ncycle` (driver.cpp:467–473), then outputs.

An observation after due seeding sees final `w0`, received destination ownership,
and final parent cell. Before driver counter advancement it must explicitly use
`t_n+dt`, generation `ncycle+1`. Alternatively sample immediately after counter
advancement and before output/restart. Labeling it mesh.time while still in the
post-integrator task list creates an actual one-step timestamp error. Sampling
before MC movement would pair a final fluid state with an old parent cell.
A restart's first refresh must not run the MC pusher or increment its RNG cycle.

At initialization, boundaries/primitives are refreshed and initial seeding runs
before ordinary output (`driver.cpp:345–349`). Initialize new encounter state
**after** seeding and observe once at actual mesh.time, with no invented prior
parent or upstream evidence. Existing `AppendParticles` (265–292) reallocates
only real/int arrays and copies only those; side-view shock state must be
explicitly preserved and initialized for new IDs. Due schedules may catch up
multiple events at one completed step (`376–490`): stored creation time is the
requested `sched.next_time`, while the eligibility fluid and first available
observation are at `t_n+dt`. Preserve both; do not relabel scheduled creation as
an actually sampled fluid time. Such new particles cannot be assigned a prior
transfer in that just-completed bin.

For MC histories, required transfer evidence is prior and destination parent
cell (global/logical cell plus GID provenance), move direction, time bin
[t_n,t_n+dt], creation/sample times, current seed/stream metadata and flags for
unobserved/ambiguous membership transitions. A line interpolated between random
cell centers is not a deterministic material path or physical crossing oracle.
The previous entry-upstream logdistance guard is only an experimental detector
heuristic, not a validated MC one-passage rule; stochastic returns and stays need
independent continuity tests and possibly explicit ambiguity.

## Ranked findings / required restrictions

**P1: saved-flux averaging is not the actual RK3 mass update.**
`hydro_fluxes.cpp:370–391` gives each stage dt/3; the driver's RK3 coefficients
(driver.cpp:127–143) give final flux weights1/6,1/6,2/3. Unequal stage fluxes
therefore yield biased parent transfer probabilities. RK4 additionally modifies
`u1` rather than retaining rho_n (hydro_tasks.cpp:135–146). Restrict the initial
shock-MC slice to RK2 (or RK1 after its separate qualification); do not advertise
RK3/RK4 mass-flow equivalence. This is inherited baseline behavior, not a defect
introduced by the frozen detector. Independent exact-rational counterexample is
in `arithmetic-checks.json`.

**P1: missing finite/probability-support validation permits biased sampling.**
`particles_lagrangian_mc.cpp:546–553` checks only mass<=0. NaN mass or saved flux
can silently turn a transition into no move; outward probabilities summing>1
are truncated according to face ordering (563–582) instead of representing the
requested flux transfer. Need an enabled-scope guard or demonstrated supported
probability bound with a loud error on invalid cells. Uniform CFL evidence alone
cannot prove the bound for evolving flows. Do not renormalize overoutflow to1:
that changes the mass-flow contract. Both directional and no-move probabilities
must be validated, with a fixture forcing bad support.

**P1: tag narrowing can silently destroy identity and RNG uniqueness.**
Seed generation casts int64 next_tag+q to int without a range check at460;
PTAG is int. Requested 64-bit event IDs cannot repair already wrapped native
parent tags. Until a deliberate typed native tag migration is implemented,
reject negative tags, overflow/exhaustion and duplicate seeding schedules/tags
before append; cap advertised tag support at INT_MAX. Enforce restart tag bounds
and next_tag metadata too. This is a storage contract restriction, not a reason
to round IDs through Real.

**P2: keyed RNG has known cross-tag/cycle aliases.**
The linear prehash key repeats exactly for (tag104729,cycle0) and
(tag0,cycle7919), for any seed. SplitMix cannot remove equal input keys.
This does not invalidate ordinary same-layout replay, but it forbids claiming
unique independent streams from tag/cycle alone. The same stream is deterministic
without restart state, provided tag,ncycle,seed and arithmetic match. Qualify
population statistics and disclose finite-key stream correlation; a future
key-combination change would intentionally change realizations and require new
seed/version metadata. Single-precision conversion can round the largest
53-bit variate to1, so do not claim float [0,1) semantics without qualification.

**P2: restart configuration changes are not part of a guaranteed replay.**
The MC section (855–974) stores dims, local count, next_tag, schedule progress
by schedule ID and all real/int arrays. `random_seed` is preserved indirectly
through the full ParameterInput snapshot (`restart.cpp:232–234`), not in that
particle header; changed overrides/seeds or changed schedules can alter future
population/streams. Unknown schedule IDs are silently ignored (945–952), header
counts lack negative/size sanity checks, duplicate schedule IDs can copy progress
to multiple schedules. Shock restart contract must bind resolved seeding/RNG
configuration, counts and ID mapping or reject changes. Existing MC restart
writes/reads require per-rank files for>1MPIrank and the same rank count
(`restart.cpp:680–689`, `pgen.cpp:1728–1754`); changed-rank restart is explicitly
unsupported, not statistically tested by inspection.

**P2: thermo history is not the required event persistence.**
`prtcl_thermo_history.cpp:239–255` records four native int fields
(cycle,tag,seedid,GID), time and xyz plus configured sampled fields. It has no
source-parent transfer, encounter sequence, active state, event identity,
validity/censoring brackets or event watermark. Its append-only file cannot
provide event restart reconciliation. Retain a separate strict event contract.

## Required MC qualification gates and predeclared oracle proposals

TRC09, TRC10, TRC11 and RST10 become mandatory for selected MC support. Deterministic
TRC01/02 analytic paths are substituted by mass-flow statistics; retain TRC03
geometry, TRC04 timestamps, TRC06 underresolved passage, TRC08 compaction,
EVT/IO/RST structure and disabled-dynamics comparisons.

- **TRC09 controlled uniform translation:** actual fluid solver + native pusher;
  periodic uniform rho/P, single positive/negative axis velocity; parent cell
  transition probability p=|v|dt/dx independently known. For N particles, one
  step has expected directional count Np and no-move count N(1−p). Predeclare
  population N=8192, seeds[17,101,12345,8675309], a finite simultaneous interval
  method (e.g. exact binomial tail with predeclared family error0.01 across all
  cases), and geometric exactness of one-cell move/no move. Extend to two
  nonzero axes with multinomial expected counts and mass-weighted nonuniform
  populations/known flux. Keep individual seeds and denominators, not only a
  pooled passing mean. Native RNG generation may be used for replay, not to
  generate the expected population probabilities. The interval proposal is not
  accepted or executed evidence; coordinator/owner must freeze it before results.
- **TRC10/RST10 replay:** fixed executable/backend/precision, timestep sequence,
  mesh decomposition, seed schedule and tag configuration; compare uninterrupted
  versus restart for full particle state, parent-transfer timelines, generation,
  event identities/state/watermarks and scheduled append across restart. A
  stateless RNG needs no stored draw counter, but restart must preserve ncycle
  and the exact seed/config. Rank-layout changes require independent initial
  seed selection and flux equality; floating mass-prefix accumulation can change
  seeding outcomes under rank grouping. Same rank-count restart restriction
  remains until actual new support is implemented and exercised.
- **TRC11 shock-cell transfers:** known mask cells with explicit source/destination
  parents; stay in inside cell, enter/leave, stochastic return, transfer across
  a thin moving mask that neither endpoint resolves, and newly seeded inside
  particles. Verify time bins, ambiguity/cadence flags and no invented exact
  crossing. Native finite-volume shock pilot reports ensemble encounter
  distributions and left/right censoring; a deterministic one-event per ideal
  velocity path count is not an MC acceptance oracle.
- **TRC04/INT time:** a manufactured time-varying diagnostic and source-controlled
  final primitive identifies stale-step sampling; check one observation per
  completed step, after migration and append; no extra RK-stage observation.
- **IO/RST adversarial:** active and pending episodes migrate/compact; due append
  into existing population; zero-particle rank; event saturation; malformed
  headers and unsupported changed rank/config; exact persistent sequence/
  lineage recovery. Larger event IDs than2^53 remain typed in state/readers,
  while native MC PTAG bounds must be separately enforced.

No runtime result is claimed here. No hardware beyond the available local CPU
was accessed. AMR/cooling/MHD/relativistic support and human sign-off stay outside
this audit. Next gate: integrated source review followed by bounded actual-MC
uniform-transfer, initial/due-seeding and restart tests with frozen source hashes.
