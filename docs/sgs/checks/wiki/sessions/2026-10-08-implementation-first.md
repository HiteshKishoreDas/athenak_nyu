---
title: Subgrid turbulence in AthenaK — implement the published method first
class: frozen
created: 2026-10-08
---

# 2026-10-08 — implementation-first decision

## Operator decision and documentation request

After discussing the Schmidt–Federrath SGS model and its validation, the operator said:

> Before any extension, I would first liek to implement this method n athenaK

The next instruction was:

> document this in the wiki

“This method” refers to Schmidt & Federrath (2011),
[A Fluid-Dynamical Subgrid Scale Model for Highly Compressible Astrophysical Turbulence](https://arxiv.org/abs/1010.4492),
A&A 528, A106, DOI 10.1051/0004-6361/201015630. The agreed priority is a faithful
AthenaK implementation before extensions. AMR remains part of the original ambition.
The assistant proposed establishing the hydrodynamic model on uniform grids first.
This sequencing is a proposed implementation brief, not an adopted detailed run plan.

## Proposed first deliverable

A working implementation should evolve SGS kinetic energy, apply the published
mixed stress and turbulent pressure, include transport and production, and return
SGS dissipation to heat. Resolved-plus-SGS energy exchanges must be consistent.
Forcing, conservative shock heating and numerical corrections need separate accounting.

Initial checks proposed in the discussion were uniform-flow preservation, SGS decay
and heating, energy conservation, then driven turbulent boxes across resolutions.
The living subject expands these into a proposed test sequence. No test has passed
or been run. AMR energy compensation, MHD and cooling/multiphase extensions are later
work, rather than prerequisites for implementing the original HD model.

## Read-only checkout investigation

Candidate: `/Users/hitesh/hitesh/git/athenak_nyu`.

- Branch: `gotham_development`.
- HEAD: `07c529202befbeb6ef5e694d95286929618ff7c1`.
- An older checkout, `/Users/hitesh/hitesh/git/athenak_nyu_old`, was also found.
- The candidate had an existing tracked `.gitignore` change and an untracked
  `AGENTS.md`; these were not edited or adopted into an implementation commit.
- Read its repository instructions, hydro allocation and task flow, scalar fluxes,
  conservative-to-primitive conversion, timestep logic, diffusion/source examples,
  turbulence problem/input, output hooks and build/test infrastructure.

Relevant files include `src/hydro/hydro.cpp`, `hydro_tasks.cpp`,
`hydro_fluxes.cpp`, `hydro_newdt.cpp`, `src/eos/ideal_hyd.cpp`,
`src/diffusion/viscosity.cpp`, `src/srcterms/srcterms.cpp`,
`src/pgen/turb.cpp` and `inputs/hydro/turb.athinput`.

Existing scalar transport and boundary machinery appear useful for carrying the
SGS reservoir. This is an integration opportunity, not proof that an ordinary
passive scalar alone implements the closure. The existing scalar floor can change
energy without compensating the resolved reservoir; positivity needs deliberate
treatment. Existing passive-scalar indices must retain their meaning, including
the cooling code's use of the first scalar as metallicity.

The assistant asked whether to use this checkout/branch or another target. No answer
was received before the documentation request. The checkout is therefore a candidate,
not a confirmed implementation destination. No SGS code was written, no checkout
was modified, and no build, simulation or performance test was performed.

## Record changes and limits

Updated the existing living subject and science index, appended the science log and
reference record, and saved this session. The frozen 2026-10-06 digest remains intact.
Stage stays `framed`, with no formal verdict or floor. The new priority supersedes
the earlier statement that all model choices were open; detailed numerical and run
choices still are open.

This filing adds no new literature verification or validation measurements.
The earlier discussion of how thoroughly the model was tested remains a literature
assessment, not evidence of an AthenaK reproduction.
