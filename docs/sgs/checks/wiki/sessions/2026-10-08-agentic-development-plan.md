---
title: AthenaK SGS turbulence — comprehensive agentic development planning
class: frozen
created: 2026-10-08
---

# 2026-10-08 — agentic development plan

## Request

> Make a full thorough, exhaustive plan with todos, milestones, all possible tests fro implementing this feature using an agentic developemnet

Interpreted as a detailed planning/documentation task in the existing wiki context.
The preceding operator decision remains authoritative: implement Schmidt–Federrath
in AthenaK before extensions. This session produces a plan, not the implementation
or a dispatched agent team.

## Outputs

- [Development plan](../development-plan.md): equations and numerical decisions,
  integration map, 10 baseline milestones M0–M9, later A0–A2 AMR milestones,
  extension boundary X, dependency graph, checkbox TODOs, campaign tiers,
  threshold protocol, risks and release inventory.
- [Test catalogue](../test-matrix.md): 230 unique test-family IDs with setup,
  independent oracle/outcome and baseline/support/diagnostic/AMR/future labels.
- [Agent runbook](../agent-runbook.md): responsibilities, file ownership,
  parallelism/dependency rules, concrete work packets, assignment templates,
  handoffs, run manifests, build recipes and gate/failure policy.
- Updated living subject/index and appended reference/science logs.

The plan and runbook contain 212 unchecked TODOs. Catalogue coverage is 151 baseline,
28 conditional support-path, 11 diagnostic, 32 AMR and 8 future-physics families.
An interaction table specifies critical FOFC/MPI/GPU/restart/scalar/AMR combinations.

The tests span configuration, algebra, state lifecycle, smooth PDE/ODE verification,
energy partition, positivity/fallback, shocks, timesteps, boundary/parallel paths,
restart/output, filtering, evolved boxes, cost and later AMR/physical extensions.
No finite catalogue literally covers every possible future flow. New support paths
must expand the coverage contract. No family was executed or marked passed.

## Evidence revisited

Reopened the primary 2011 text, especially §§2–4, for mixed-stress normalization,
coefficients, SGS transport, paired energy equations and pseudo-isothermal evolved
tests. Reopened the 2014 AMR text and Appendix A.2 for later compensation/remap
work; did not adopt its algorithm without a separate AMR contract.

Read-only local investigation of `/Users/hitesh/hitesh/git/athenak_nyu` at
`07c529202befbeb6ef5e694d95286929618ff7c1` reconfirmed the candidate revision.
Read repository guidance, hydro task/fallback paths, restart state counts, output
configuration, global timestep selection, forcing configuration and an existing
hydro regression script.

Particular integration risks recorded: the fallback predictor covers hydro rather
than all SGS state; corrected fluxes can invalidate previously calculated stress
work; ordinary scalar floors can alter energy without paired compensation; first
passive scalar has a metallicity meaning; forcing restart state requires explicit
verification. These are planning findings from inspected paths, not proven failures
of an SGS implementation that does not yet exist.

## Decisions and open choices

The plan proposes initial Newtonian ideal-gas HD on uniform Cartesian grids,
then parallel/restart/backend qualification and scientific validation. AMR follows
the baseline release. A pseudo-isothermal reproduction lane is distinguished from
ordinary-gamma heating verification. Published coefficients stay fixed in the
baseline; calibration is a separate study.

The target checkout/base is still unconfirmed. State layout, pressure/solver
treatment, positivity/startup policy, stencil, support matrix, exact tolerances,
reference data and resource budget remain contract-review items. Proposed thresholds
are not measured achievements. Hardware/backend support requires actual access and
tests. No schedule or allocation was created and no agents were spawned.

Stage remains framed rather than claiming an approved numerical specification,
scientific verdict or implemented feature. Earlier frozen sessions are unchanged.
No AthenaK edits, solver build, simulation, test campaign or performance measurement
were performed. File/link/coverage bookkeeping checks verify the plan, not the physics.
