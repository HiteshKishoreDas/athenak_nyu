---
title: AthenaK SGS turbulence — agent work packets and review runbook
class: living
created: 2026-10-08
last_touched: 2026-10-08
---

# Agentic development runbook

Use this with the [development plan](development-plan.md) and
[test catalogue](test-matrix.md). This document describes future development.
No agents, jobs, schedules or AthenaK code changes were launched in this planning session.

## Implementation status on 2026 10 08

Restricted HD code and independent verification are now recorded in [implementation progress](implementation-progress.md) and the [execution ledger](coverage-ledger.md). M0–M9 scientific release is incomplete. Original planning statements and checklists below are preserved as their historical record.

## 1. Operating model

The coordinator owns integration and the evidence ledger. Specialist agents work
on bounded packets. Review is independent enough to catch equation, index and
budget errors: a reviewer must not approve a test solely because its expected
value comes from the implementation being tested.

| Role | Owns | Required output |
|---|---|---|
| Coordinator / integrator | Scope, dependencies, shared-file ownership, branch integration, milestone closure | Packet ledger, coherent diff and gate report |
| Science/numerics reviewer | Equations, signs, units, stencils, temporal treatment, admissibility and filter assumptions | Independently derived oracles and explicit findings |
| Implementer | One coherent code slice on a pinned base | Focused diff, runnable verification and observed results |
| Verification/reproduction agent | Independent execution, boundary cases, regression and defect injection | Raw evidence, pass/fail decisions and coverage gaps |
| Reference-analysis role | Filtering, a priori tests, statistics and turbulent-box protocol | Validated analysis tools, uncertainty and comparison report |
| Performance/portability role | Available CPU/GPU/MPI qualification and resource estimates | Hardware/build manifest, profiles and measured overhead |

These are roles, not six mandatory concurrent agents. A practical arrangement
uses coordinator plus up to three workers. Rotate responsibilities across phases.
An author can run their own checks, but final approval comes from another role.
Do not invent a scheduler, orchestration service or elaborate agent framework;
use available task tools, Git and the existing test harness.

## 2. What can run in parallel

Before M1, parallel work is read-only: equation derivation, lifecycle audit and
baseline/tool inventory. After M1, M2 kernels/reference calculations and M3 state
lifecycle may proceed independently if interfaces are frozen and file ownership
is disjoint. Reference-analysis preparation can run beside deterministic solver
verification using synthetic fields.

The coordinator serializes integration into shared hydro, EOS, task, configuration
and output files. Different worktrees prevent filesystem interference but not
semantic conflicts; review the combined state rather than trusting separate green
branches. Changes to equations, state layout or stage timing invalidate affected
packets and their dependent evidence.

Long statistical runs start only after deterministic gates and a resource-capped
pilot. Once a frozen campaign is running, independent analysis may proceed on
completed outputs. A failing conservation test stops dependent interpretation,
not unrelated read-only preparation.

## 3. Workspace and ownership rules

- [ ] Confirm the operator's target checkout/base once, then pin the actual commit.
- [ ] Read applicable instructions before assignments; distinguish wiki and code worktrees.
- [ ] Inventory changes without resetting anything. Preserve staged/untracked user work.
- [ ] Use isolated development branches/worktrees when needed; proposed branch prefix
  is `codex/`, with names recorded in the packet ledger.
- [ ] Give each mutable file or coherent subsystem one active owner.
- [ ] Give each agent its working directory, base revision, allowed files and forbidden scope.
- [ ] Allocate separate build/output directories and unique run basenames.
- [ ] Do not run a test harness that removes another agent's build directory.
- [ ] Integrate one reviewed packet at a time; preserve source commits and handoff evidence.
- [ ] Never include unrelated local edits in an implementation or wiki commit.
- [ ] Record local commit, push, PR, review and merge independently; do not equate them.
- [ ] Do not send messages to outside collaborators or other user-owned chats merely
  because an agent suggests it. Use the authorized task's own coordination tools.

## 4. Standard packet contract

Every assignment contains these fields in plain language:

| Field | Content |
|---|---|
| Objective | One observable deliverable, not “work on SGS” |
| Base | Exact checkout/worktree and commit |
| Dependencies | Reviewed packets and interfaces required first |
| Allowed edits | Explicit files/subsystem, with shared-file owner |
| Invariants | Equations, units, conservation, positivity and legacy meanings to preserve |
| Not in scope | Extensions or refactors to defer for this packet |
| Exit tests | Catalogue IDs and concrete instances, expected metrics and thresholds |
| Resource cap | Time, ranks, memory, output size and largest allowed grid |
| Evidence | Input/build hashes, logs, metrics, raw outputs and failure records |
| Handoff | Changed paths/commits, executed checks, open risks and next dependency |

No agent silently expands a packet from HD to MHD/AMR or changes a coefficient,
filter convention, tolerance or reservoir definition. Propose the change with a
reason; the coordinator reviews whether it fits the operator's scope and records
the decision. Routine implementation choices can be resolved within the locked
contract without asking the operator at every step.

## 5. Concrete work packets

### P0 — reproduce the untouched baseline (M0)

**Owner:** verifier; coordinator handles workspace selection.

- [ ] Pin source/build/toolchain/hardware and local changes.
- [ ] Build serial CPU debug and release in separate directories.
- [ ] Run selected existing hydro/scalar tests and a short turbulence input.
- [ ] Return baseline failures, minimal reproduction, timings and resource estimates.

**Edits:** manifests/results only; no solver code. **Exit:** reproducible baseline,
with pre-existing failures separated from environmental blockers. No full campaign.

### P1 — lock the model and numerical interfaces (M1)

**Owner:** science/numerics reviewer, with read-only integration audit from implementer.

- [ ] Independently check equations, signs, coefficients, units and omitted terms.
- [ ] Review proposed storage, SGS zero-state treatment and positivity policy.
- [ ] Derive wave/timestep, stress-work and source-stage requirements.
- [ ] Trace scalar, fallback, BC, restart and forcing paths in the confirmed revision.
- [ ] Freeze supported modes and concrete test instances/thresholds.

**Edits:** model contract and term map. **Exit:** no unresolved energy/sign/stage
issue; coordinator and independent reviewer accept interfaces. A design disagreement
is recorded and resolved before dependent code changes.

### P2 — closure kernel plus independent oracle (M2)

**Owner:** implementer writes production helper; reviewer writes/checks oracle separately.

- [ ] Implement mixed stress, production, diffusivity and dissipation primitives.
- [ ] Parameterize hand cases and deterministic random matrices.
- [ ] Cover K01–K25 as applicable, especially trace, symmetry, sign and zero norms.
- [ ] Show chosen sign/factor mutations are caught by the verification suite.

**Edits:** focused helper/test files and build hookup only. **Exit:** kernel
results agree with independent mathematics on the declared domain and backend.

### P3 — state lifecycle and compatibility (M3)

**Owner:** integration implementer, single owner of hydro/EOS/output count changes.

- [ ] Add optional SGS state and preserve ordinary scalar identities.
- [ ] Carry state through allocation, primitive conversion, reconstruction, advection,
  stage copies, ghost communication, initialization, output and restart.
- [ ] Exercise C/L cases and relevant I restart/legacy tests.

**Exit:** transported SGS reservoir works; disabled path remains compatible.
Do not call this a completed turbulence model.

### P4 — conservative coupled solver (M4)

**Owner:** implementer; numerics reviewer checks discrete work and RK coupling.

- [ ] Add stress momentum/energy fluxes, SGS pressure, SGS diffusion and paired exchanges.
- [ ] Integrate actual timestep constraints and valid stage dependency/data flow.
- [ ] Add applied-flux/source diagnostics and independent inventory checks.
- [ ] Run P01–P24 as applicable, including smooth convergence and full manufactured case.

**Exit:** component reservoirs and represented total energy behave correctly;
temporal/spatial order established. No long turbulent run used as a substitute.

### P5 — positivity and fallback (M5)

**Owner:** implementer; verifier constructs adversarial admissible/near-limit cases.

- [ ] Implement bounded transfer and consistent fallback on shared faces.
- [ ] Ensure no source survives from a stress removed by flux correction.
- [ ] Record every allowed floor/limiter energy change.
- [ ] Exercise R families and first-step/restart/global dt cases.

**Exit:** supported hard-flow cases pass; correction frequency/energy impact
quantified. Restrict unsupported regimes rather than hiding failures with looser thresholds.

### P6 — persistence, parallelism and diagnostics (M6)

**Owner:** independent verifier and available-backend performance role.

- [ ] Run I families for claimed block/rank/backend/precision paths.
- [ ] Verify forcing random state and integrated actual work across restart.
- [ ] Check history/field output against independent active-cell sums.
- [ ] Measure H01–H07; report hardware/resource gaps without simulating success.

**Exit:** a second agent reproduces qualified execution paths and energy ledger.

### P7 — a priori analysis pipeline (M7)

**Owner:** reference-analysis role; reviewer checks weighting/sign/filters independently.

- [ ] Validate filters/derivatives on synthetic data before reference turbulence.
- [ ] Run F families on a justified useful reference range with fixed coefficients.
- [ ] Report signed errors, conditional failures and uncertainty, not correlation alone.

**Exit:** reproducible analysis pipeline and closure assessment. This packet can
prepare synthetic tools earlier, but scientific results require validated source data.

### P8 — evolved-box pilot and campaign (M8)

**Owner:** coordinator approves resource-capped protocol; verifier executes; analyst evaluates.

- [ ] Pilot T01–T03 to set saturation/time/storage and acceptance protocol.
- [ ] Freeze selected forcing/resolution/seed matrix and comparison metrics.
- [ ] Run core T families with exact source/build/input provenance.
- [ ] Include null/negative outcomes and correction/numerical-method sensitivity.
- [ ] Report H08 error-versus-cost and practical limits.

**Exit:** validation report reproducible from raw metrics and inputs; no retuning
on held-out results. A crashed job stays failed in the ledger even after later fixes.

### P9 — independent release review (M9)

**Owner:** verifier not responsible for the principal implementation; coordinator integrates fixes.

- [ ] Check every equation against final code and every claimed path against evidence.
- [ ] Reproduce minimal example, conservation case and one held-out box comparison.
- [ ] Inspect compatibility, limitations, documentation and clean-checkout dependencies.
- [ ] Prepare final diff/PR and attach test/performance evidence when that workflow is used.

**Exit:** baseline HD release with honest scope and durable reproduction package.

### PA0–PA2 — separate refinement programme

Repeat P1-style contract review for inter-level partition/positivity, then implement
remap arithmetic, fixed-hierarchy transport/reflux and dynamic regridding in that
order. Use A01–A32. Add no MHD/cooling closure under the cover of AMR implementation.

## 6. Ready-to-use assignment prompts

Fill the base/path/resource placeholders from the packet ledger before dispatch.
These prompts are templates, not messages sent during this session.

**Science review:**

> Independently derive the Schmidt–Federrath HD SGS stress, production, dissipation,
> transport and represented-energy budget using the pinned model contract. Audit
> signs, energy density versus specific energy, exact-zero gradient/SGS behaviour,
> compression and positivity. Inspect the proposed numerical treatment in the pinned
> AthenaK revision. Return concrete findings and independent hand/ODE/PDE oracles.
> Do not edit solver files, tune coefficients or introduce AMR/MHD/cooling extensions.

**Implementation:**

> Implement only packet [ID] in [WORKTREE] at [BASE]. Follow its allowed files,
> locked interfaces, support matrix and invariants. Preserve existing user changes
> and passive scalar meanings. Leave the smallest focused runnable verification
> using the existing harness. Run [TEST INSTANCES] within [RESOURCE CAP]. Return
> exact changed paths/commits, actual results, failures and remaining concerns.
> Do not reinterpret a planned test as passed or change thresholds silently.

**Independent verification:**

> Review packet [ID] without relying on the production helper as its oracle. Run
> the declared tests from the recorded build and inputs, including the listed
> boundary/admissibility/decomposition cases. Confirm component energy partition
> as well as total conservation. Record raw evidence and unsupported/unavailable
> paths. Return pass/fail findings and the smallest reproduction for each defect.
> Do not rewrite unrelated code or launch the full turbulence campaign.

**Campaign analysis:**

> Use the frozen campaign/filter protocol and recorded raw runs. Verify density
> weighting, effective scales, actual forcing work, correction budgets, stationarity
> and uncertainty before comparing model-on/off cases. Hold coefficients fixed and
> retain null/negative outcomes. Return metrics and reproducible figures with each
> claim tied to its run manifest. Do not infer universal accuracy from one spectrum.

## 7. Evidence and handoff formats

Packet status: `planned → ready → in progress → review → verified → integrated`.
Failures and resource blockers are explicit flags, not a misleading success state.
Do not create a second evidence platform; plain manifests, CSV/JSON metrics and
Markdown summaries beside existing test outputs are sufficient.

Each run manifest records:

- Source commit plus dirty-patch hash, input hash and actual command/overrides.
- Compiler/Kokkos versions, build type/options, precision and execution backend.
- Hardware/job allocation, ranks/threads, block/pack/grid layout and boundary conditions.
- Model coefficients, filter convention, SGS seed/initial energy convention.
- Solver/reconstruction/integrator/CFL, physical units and forcing protocol/seed/state.
- Run duration/steps, output cadence and analysis version/filter/sampling window.
- Exit status, conservation/admissibility metrics, correction budgets and resource cost.
- Raw log/output locations and test thresholds with rationale.

Handoff example structure:

1. Packet, base and submitted commits/changed paths.
2. What changed and which invariant it implements.
3. Tests actually executed, exact metrics and evidence locations.
4. Failures, not-run/unsupported paths and limitations.
5. Shared-interface changes and dependent evidence needing rerun.
6. Proposed next packet; no automatic expansion of scope.

The coordinator checks the handoff against actual diffs and evidence. An agent's
confident summary does not replace run logs or independently recomputed budgets.

## 8. Baseline build/run recipes

**Execution environment: confirmed AthenaK worktree on a machine with the required
compiler and Kokkos dependencies.** These are proposed commands derived from the
inspected checkout; revalidate them at M0. The wiki directory is not a solver build.
Set `ATHENAK_SGS_WORKTREE` to the confirmed isolated source directory.

```sh
cmake -S "$ATHENAK_SGS_WORKTREE" -B "$ATHENAK_SGS_WORKTREE/build-sgs-debug" \
  -DCMAKE_BUILD_TYPE=Debug -DKokkos_ENABLE_SERIAL=ON
cmake --build "$ATHENAK_SGS_WORKTREE/build-sgs-debug" --parallel 8

cmake -S "$ATHENAK_SGS_WORKTREE" -B "$ATHENAK_SGS_WORKTREE/build-sgs-release" \
  -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_SERIAL=ON
cmake --build "$ATHENAK_SGS_WORKTREE/build-sgs-release" --parallel 8
```

Expected result: separate debug/release builds containing the target's AthenaK
executable. Use the actual discovered executable path, not an assumed output location.
Run inputs in unique output directories and record all overrides.

**Execution environment: target worktree's `tst/`, isolated from other agents.**
The inspected harness exposes:

```sh
python3 run_tests.py hydro
```

It builds/cleans its own test `build/` directory. Confirm dependencies and cleanup
before invocation; choose specific tests for iteration. GPU/MPI commands are written
only after actual compiler/hardware/module availability is known. No invented
cluster login, allocation or CUDA architecture belongs in a reproducible recipe.

## 9. Failure, review and completion rules

- Never change expected values to match code without re-deriving the oracle.
- Reduce a failure to the smallest case before adding more physics or larger grids.
- Preserve failing evidence; add a fix/retest entry with commit and rationale.
- Fix sign/count/stage bugs at their shared cause, not by special-casing one test.
- If a gate fails, stop dependent claims/runs; continue independent useful preparation.
- If a resource is missing, mark its coverage blocked and limit support claims accordingly.
- Rerun affected tests after fixes; repeat an expensive campaign only when new changes
  invalidate it, unresolved sensitivity requires it or the reference protocol changes.
- Do not conceal energy defects by subtracting corrections without also reporting their cost.
- Distinguish a correct implementation with limited scientific benefit from a broken implementation.
- Close M9 only when a clean independent reproduction succeeds and all mandatory
  baseline families and claimed support paths have evidence.

## 10. Initial dispatch checklist

- [ ] Operator's target checkout/base known.
- [ ] Existing changes preserved and isolated workspace attached/recorded as appropriate.
- [ ] Coordinator assigns P0 and read-only P1 preparation with separate output locations.
- [ ] No agents currently share mutable files or a test build directory.
- [ ] Baseline resource cap and available backends recorded.
- [ ] Model/source references and this plan linked in each assignment.
- [ ] Packet/evidence ledger initialized with all tests not run.
- [ ] Review gates and integration ownership understood before solver edits.
