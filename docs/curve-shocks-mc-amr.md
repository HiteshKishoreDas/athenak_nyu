# MC AMR encounter-state support — living evidence

## Scope and decisions

2026-10-17. Working branch `codex/curve-shock-histories`, initial HEAD
`dd15995315cf85041740f8cefbe64e0c1bff3d87`, separate from preserved
`codex/shock-multiphase`. Delegated ownership: native MC source, particle observer
source, particle boundary source, focused tests and this document.

Accepted: initial AMR support and MC/curve histories coexisting. Detector support
width in leaf cells versus physical units is still pending. No width semantics,
threshold, floor, forcing, flux probability, RNG or refinement-adjustment rule has
been changed. Scientific approval and actual RUSTY/GPU evidence remain separate.

## Inspection and implementation

- Native MC boundary selection chooses fine subfaces or falls back to coarse
  neighbors; native `AdjustMeshRefinement` then uses destination fine-face fluxes
  and the original RNG to select child cells, or applies coarse-cell parity
  corrections. Those transport methods remain untouched.
- Enabled boundary selection now follows that native level/subface rule, with
  bounded fallback within the actual face/edge slots. Outflow still precedes
  missing-neighbor lookup and complete state still follows real/int migration.
- Native regrid remapping locates the new leaf, snaps to its cell center and
  resets native move/level metadata. Complete typed EncounterState now follows
  the identical ordering through serial copies and a checked MPI byte exchange.
  Historical observations, IDs, sequences, counters, flags and periodic images
  are preserved. State extent, destination rank, payload totals and byte
  displacements are checked before copying or MPI multiplication/allocation.
- Enabled remap rejects unexpected nonfinite/out-of-domain positions rather than
  silently wrapping/clamping without image or exit evidence. Normal boundary
  handling must have consumed physical exits before the regrid transaction.
- Detector scratch shape is checked before every observation and reallocated
  when the local block count changes. Its active cells are regenerated then.
- Startup MC AMR exclusions are currently restored: complete state preservation
  passes, but native regrid mass allocation fails the independent occupancy oracle.
  Coordinator owns completed-step observation timing after AMR and curve changes.

## Scientific blocker: inherited MC regrid mass allocation

The native-remap occupancy probe places \(4096\) distinct MC identities at a
legitimate coarse-cell center \((0.125,0.125)\), with coarse cell width
\(\Delta x=\Delta y=1/4\). Replacing its parent block with four refined leaves
gives fine-cell width \(1/8\). The actual native remap sends every tracer to
\((0.1875,0.1875)\), with child counts \([0,0,0,4096]\). For uniform fluid density,
the four equal-volume children have mass fractions \(1/4\) each. The native snap
is deterministic and does not use identity/RNG to allocate that parent-cell
population. This is an inherited transport deficiency for mass-following MC
regrid semantics, independently of complete encounter-state preservation.

`RemapAfterMeshRefinement` snaps position, sets the native move code to zero and
sets last level to the new leaf level. The later native face-crossing adjustment
therefore does not perform a coarse-to-fine split for that regrid. This packet
retains the reproducer and does not invent a new mass-allocation/RNG rule.
MC mass-following AMR qualification is **not passed**. The focused harness result
means the native behavior was reproduced; it does not label that scientific
oracle as passing. Real fluid refine/coarsen identity tests remain useful and
separate. Human owner decisions and physical qualification are required before
calling the MC tracer mass-following on adaptive regrids.

## Validation in progress

`mc_amr_state_test.cpp` executes the actual native remap and boundary source,
with clearly identified constructor/count stand-ins. It replaces leaf metadata
with refined and then coarsened layouts and asserts exact equality of every semantic state field,
large integer identities, source-associated creation time, native level metadata,
local ownership, and zero-particle ranks. This checks the native remap routine;
it does not itself execute the fluid AMR transaction or supply a physical oracle.

A real native refine/coarsen transaction, observer scratch resize, coexisting
histories, same-state restart and cross-rank migration are still required.
No numerical tolerance may be relaxed to close those gates.

## Focused checks, first closure

Actual native remap compiled and passed serial and singleton MPI refined-leaf
then coarsened-parent layouts. Every current EncounterState field is compared
exactly, including integer IDs beyond binary64's exact integer range, all four
observation slots, periodic images, counters, phase and validity flags. Invalid
state extent and unexpected outside-domain remap positions reject explicitly.
After those preservation checks passed, the startup AMR exclusion was briefly
removed during development. The independent occupancy failure then established
that preservation alone was insufficient; startup adaptive/multilevel rejection
has been restored pending the user's physical allocation decision. Typed remap
support and its tests remain present for independent engineering review. The
detector algorithm and configured width are unchanged; the pending scientific
AMR width decision has not been supplied or approved by this work.

Initial harness failure: constructor stand-in omitted the native MC type, so the
native remap correctly returned without work. A later test edit referenced an
obsolete representative state slot (current events derive it from first).
Both were fixture mistakes, corrected without changing production semantics.
The first failure report is retained as `mc_amr_initial_failure.json`.

The local MPI launcher itself fails before launch: installed PMIx runtime and
PRRTE build versions are binary-incompatible. Exact version-probe output is in
`mc_amr_mpiexec_unavailable.json`. Singleton MPI native collective code executed;
it is not evidence of inter-rank transfer or RUSTY qualification.

## Reproduction commands and evidence

From this candidate checkout:

```sh
python3 tst/shocks-curves/mc_amr_state_test.py \
  --kokkos-build /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu/build-shock \
  --output /private/tmp/mc-amr-state-serial-occupancy
python3 tst/shocks-curves/mc_amr_state_test.py \
  --kokkos-build /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu/build-mpi \
  --output /private/tmp/mc-amr-state-mpi-final --mpi
```

These compile only the bounded actual-method harness against the named existing
CPU Kokkos libraries, not the coordinator's full Athena executable. Build flags,
source SHA256 values, exact compile/run commands and outputs are in the retained
`mc_amr_*` JSON reports. `CheckNeighbors` separately executes native fine-subface
and bounded coarse-fallback selection with exact destination/state checks.
The fixture supplies minimal constructors/count callbacks; it does not execute
fluid prolongation/restriction, saved-face allocation or the full mesh transaction.
All current state fields are compared, not uninitialized C++ padding bytes.

## Pending physical repair: design constraints only

The parent has asked the user to choose mass-weighted refined-fluid redistribution
or explicitly unqualified native remapping. No physical repair is implemented
while that answer is pending. Relevant native transaction details:

- The regrid routine completes native conserved-variable restriction/prolongation
  and MPI fluid transfers before replacing leaf metadata and calling particle
  remap. New active-cell conserved density is available at that point; freshly
  recomputed primitives/ghosts are established by the caller afterward. A repair
  must use the actual resulting conserved fluid mass, not assume uniform weights
  for a nonuniform fluid or reconstruct density from old shock observations.
- Reliable old leaf/cell geometry must be captured **before** old block and
  logical-location arrays are deleted. Old PGID becomes stale during regrid.
  PLASTLEVEL alone is insufficient: native Push records the pre-move level, and
  AdjustMeshRefinement can subsequently cross levels without updating that field.
  Inferring the old cell from that stale level can choose the wrong child set.
- For refinement, candidate new active cells must partition the old tracer's
  containing cell, including any ownership changes. Candidate probabilities
  would be proportional to actual post-prolongation child mass
  \\(p_c=m_c/\\sum_j m_j\\), with \\(m_c=\\rho_c\\,\\Delta V_c\\). Reject nonfinite,
  negative or zero total mass rather than invent floors or replacement weights.
  Coarsening merges into the containing coarse cell; unchanged-grid redistribution
  must not draw a new random child or move its cell center.
- Preserve integer particle identity, complete encounter state and all prior
  observations. Regridding changes the current native MC host cell; it must not
  fabricate temporal entry/exit evidence or reset an ongoing episode.
- Use a separately identified deterministic RNG stream keyed by stable native
  integer particle ID, persistent seed, completed-step/transaction identity and
  stable old-cell logical coordinates/level. Rank, current GID, memory index and
  receive order are not stable keys. Candidate ordering must be canonical across
  partitions. Reusing an existing face-transfer draw risks unintended correlation;
  the RNG domain and any transaction counter require explicit provenance and
  restart persistence. A once-per-step transaction key is sufficient only if
  that invariant is actually guaranteed, including initialization/restarts.
- Any native MC mass-allocation repair must apply consistently with histories
  enabled **and disabled**, or the detector would change tracer transport merely
  by being enabled. Exact fixed-seed on/off passivity remains a required check.
- Required new oracles: uniform child-mass occupancy, prescribed nonuniform child
  mass fractions, unchanged/coarsened layout behavior, repeated refine/coarsen,
  restart and partition-invariant draws, zero-particle ranks, active encounter
  continuity, and actual native fluid transactions. Statistical criteria must be
  specified independently before evaluating nonuniform sampling runs; the current
  deterministic collapse fixture is retained rather than tuned away.

No source change for those proposed physical rules has been made. The restored
startup rejection prevents the current branch from advertising MC shock AMR as
mass-following while this decision and qualification remain open.

## Final packet handoff

Both focused Serial and singleton MPI builds completed. In both, native refined
and coarsened layout state preservation and native fine/coarse neighbor checks
passed; invalid extent/outside-domain remap inputs rejected. Both reproduced the
same scientific occupancy failure. Final MPI source was rebuilt after a
format-only line wrap; Serial preprocessing is unchanged by that wrap inside
the MPI-only block. No statistical tolerance was added or weakened.

`git diff --check` passed. Project cpplint on owned boundary/native MC sources
reports ten inherited diagnostics on unchanged original source lines; no new
changed-line diagnostics. Exact inherited diagnostics are retained in
`mc_amr_style.log`. This is not a claim that whole-file lint is clean.

Final native MC source SHA256:
`2f924c81622de40432e0df6080860546aff226e771f4ee5242017d05e7c787e5`.
Final boundary source SHA256:
`12b49e1eeecc868cc20a9ddd2ba839b2cadf2e5e46ee797cdc166501454bb7c1`.
The particle observer file was handed back to the coordinator immediately after
restoring startup AMR rejection; subsequent shared-store extraction is coordinator
owned. No commits, full binary builds or source-branch edits were performed by
this delegated packet.

Next gates: user's physical mass-allocation and detector-scale decisions,
independent typed-state/MPI review, actual native fluid refine/coarsen with
coexisting histories and passivity/restart checks, then working multi-rank runtime
and RUSTY qualification. Scientific approval remains pending.
