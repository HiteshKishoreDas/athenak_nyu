---
title: Subgrid turbulence in AthenaK
aliases: [athenak-subgrid-turbulence, AthenaK SGS turbulence, subgrid turbulence with AMR]
class: living
stage: framed
created: 2026-10-06
last_touched: 2026-10-12
floor: ""
ceiling: "Proposed: conservative turbulence modelling in AthenaK whose resolved-plus-subgrid predictions remain consistent across changing AMR levels, with validated heating and mixing in radiatively cooling multiphase gas."
domains: [compressible turbulence, large eddy simulation, adaptive mesh refinement, numerical dissipation, radiative cooling]
expertise_confidence: "High for compressible-flow energy bookkeeping; medium for SGS closure selection and statistical validation; limited for AthenaK integration internals. Literature and selected upstream pages inspected; restricted uniform HD implementation and deterministic/pilot verification recorded on 2026-10-08; full scientific validation remains incomplete."
prior_art_checked: ""
related: [amr-criteria-multiphase-turbulence, cascade-under-cooling, numerical-dissipation-density-bias]
source: Direct user idea, captured and digested on request 2026-10-06; implementation-first decision documented on request 2026-10-08
---

# Subgrid turbulence in AthenaK

Add subgrid turbulence to AthenaK, then repeat the turbulent-box tests with
proper validation, including tests with adaptive mesh refinement (AMR).

## Implementation status on 2026 10 08

Restricted HD code and independent verification are now recorded in [implementation progress](implementation-progress.md) and the [execution ledger](coverage-ledger.md). M0–M9 scientific release is incomplete. Original planning statements and checklists below are preserved as their historical record.

## Full-window evidence on 2026 10 11

Matched 64³ PLM/PPM4 runs finish 10T with 2T–10T statistics and applied heating
budgets. PPM4 mean SGS Mach agrees closely with the paper, but mean SGS energy
remains about half its value; both cases fail the frozen operational drift screen.
Gamma controls establish a heating contribution to Mach drift. Forcing algebra
and Gaussian analytic checks pass within their stated scopes. Exact paper forcing,
turbulent filtered reference, resolution scaling and independent-seed uncertainty
remain open. See [current port evidence](multiphase-port.md) and the
[full-window packet](evidence/2026-10-11-full-window/README.md).

## Timestep and cascade evidence on 2026 10 12

Matched 64³ CFL .1/.05/.025 restarts pass native checks with mean SGS-energy
changes below .002% over 4T–4.5T. A direct uniform-decay source-weight audit
matches full RK2. Preliminary 32³–64³ energy scaling is near Kolmogorov 2/3,
but resolved spectra lack a clean inertial range and the absolute energy gap
remains unvalidated. See the [CFL/cascade packet](evidence/2026-10-12-cfl-kolmogorov/README.md).

The subsequent direct measured-spectrum extrapolation gives 64³ unresolved
energy .03536 versus .03742 SGS using matching nine-snapshot averages.
Fit-window estimates .02947–.04916 include the evolved SGS energy under an
assumed spherical Nyquist cutoff. No clean inertial range or exact compressible
Favre/filter equivalence is established; this does not settle the paper gap.
See the [direct spectrum-tail estimate](evidence/2026-10-12-spectrum-tail/README.md).

The subsequent native filtered-snapshot test uses fresh SGS-off 32³/64³ flows
through 10T with fixed forcing. Over nine 2T–10T snapshots, the primary 64³ Gaussian
comparison gives mean production 7.1% low, correlation .995 and local RMS error
13.1%; trace-free correlation remains .990. Arithmetic and inventory checks pass.
Reference energy still changes materially with resolution, and this test supplies
measured K rather than evolving it: the absolute SGS-energy gap remains open.
See the [paper-style production test](evidence/2026-10-12-filtered-snapshot/README.md).

## Unforced decay evidence on 2026 10 12

A matched initially Mach9, unforced3D experiment reaches10T with128³/64³
references and8³/16³ SGS-on/off lanes. All six runs and488 numerical checks pass.
Relative to the filtered128³ reference, mean SGS K over1–10T is .75/.79
and falls to .45/.44 at10T. The
low-energy tendency therefore also appears without ongoing forcing. However,
the reference energy changes17.5%/34.6% between grids, failing the frozen10%
screen, and the fine flow becomes subsonic (Mach.90) by10T. This neither
qualifies absolute model accuracy nor reproduces the original paper benchmark.
No solver equations or coefficients changed. See the
[unforced decay packet](evidence/2026-10-12-unforced-decay/README.md).

## Higher unforced reference on 2026 10 12

A fresh matched256³/128³ campaign reaches10T with all six runs and488 checks
passing. The evolved SGS/reference mean K ratios are .70/.67 at8³/16³ and
.42/.35 at10T. Reference sensitivity falls to9.67% at the8³ filter, meeting
the frozen10% screen, but remains21.52% at16³. The energy deficit therefore
persists against the finer reference; the8³ result is better constrained.
The fine flow ends subsonic at Mach.896, and coarse-scale, initial-state and
single-seed limitations remain. No solver equations or coefficients changed.
This is scoped local evidence, not universal accuracy or paper reproduction.
See the [higher-reference packet](evidence/2026-10-12-higher-reference/README.md).

## Matched64³/128³ SGS on 2026 10 12

The256³ reference is now compared with the requested64³/128³ SGS-on/off grids
through10T. All six analyses and488 checks pass; three native runs are new,
with the existing128³ reference also serving as the128³ SGS-off control.
Mean model/reference K ratios are1.02/1.94, and final ratios .46/.84. Total
kinetic decay is nearly the same on/off. Thus a universal low-energy factor is
not supported by these runs. Reference qualification remains inadequate:
64³-filter sensitivity is59.06%;128³ has no independent finer reference pair.
Coefficients and native source remain unchanged. See the
[64³/128³ packet and28 plot panels](evidence/2026-10-12-sgs64-128-decay/README.md).

## Original capture

> Adding subgrid turbulence to athenaK
> And test again turbulent boxes properly, and with AMR

## Status

Framed on 2026-10-06. On 2026-10-08 the operator chose the immediate priority:
implement the Schmidt–Federrath method in AthenaK before pursuing extensions.
The original AMR ambition remains a later step. The baseline implementation and
test sequence below are a proposed brief; the target checkout is still unconfirmed.
The [detailed agentic development plan](development-plan.md),
[230-family test catalogue](test-matrix.md) and [agent runbook](agent-runbook.md)
were added on request on 2026-10-08. Their TODOs are unchecked and tests are not run;
numerical choices must be frozen at the first review gates before implementation.
No SGS code changes, builds or simulation tests have been performed. Stage remains
`framed`: this is not a completed implementation specification or scientific verdict.
`floor` and `prior_art_checked` remain blank under the tree's judgment convention.

## Implementation first — decision recorded 2026-10-08

Operator wording:

> Before any extension, I would first liek to implement this method n athenaK

In the preceding discussion, “this method” refers to Schmidt & Federrath's
[2011 compressible hydrodynamic SGS model](https://arxiv.org/abs/1010.4492).
The first deliverable is a faithful, working AthenaK implementation of that
published method. New closure physics, MHD, cooling/multiphase modelling and
AMR transfer algorithms come after the baseline implementation and validation.
This sets the immediate priority; it does not abandon the broader research idea.

Proposed baseline components:

- Evolve SGS kinetic energy, distinguishing energy density from specific energy.
- Implement the published mixed stress closure and its turbulent-pressure term.
- Include SGS advection, modelled transport, production and dissipation into heat.
- Couple the SGS and resolved energy equations consistently, so internal exchanges
  conserve their sum. Account separately for forcing, shock heating and numerical
  corrections; do not count them again as SGS injection.
- Reuse AthenaK's existing transport and boundary machinery where appropriate,
  while preserving the meaning and indices of existing passive scalars.

Proposed validation order:

1. Uniform flow: preserve a spatially uniform state and avoid spurious SGS production.
2. SGS decay and heating: check the decay law and the paired gain in thermal energy.
3. Stress and energy accounting: check the stress trace, compression work and
   resolved-plus-SGS conservation, including any positivity corrections.
4. Driven uniform-grid turbulent boxes: compare model-on/model-off runs across
   resolutions and against suitable filtered references; track energy budgets
   as well as spectra and statistical uncertainty.
5. Once the baseline passes, specify and test static refinement and dynamic AMR
   energy transfers separately, then consider additional physical extensions.

These are proposed checks, not completed tests or a fixed run campaign. The
thermodynamic setup, initial SGS energy, numerical treatment, tolerances, forcing
parameters and reference data still need a concrete implementation specification.
An energy-evolving run is needed to check deposited heat.

Read-only investigation identified `/Users/hitesh/hitesh/git/athenak_nyu` as a
candidate checkout, on `gotham_development` at commit
`07c529202befbeb6ef5e694d95286929618ff7c1`. Its scalar transport, EOS conversion,
flux/source hooks, timestep logic, turbulence generator and test infrastructure
were inspected. **The operator has not confirmed this checkout or branch as the
implementation target.** No SGS implementation was made there.

See the [dated decision and investigation record](sessions/2026-10-08-implementation-first.md).

For execution planning, use the [milestones and TODO ledger](development-plan.md),
the [test catalogue](test-matrix.md) and the [agent work-packet runbook](agent-runbook.md).
The [planning session record](sessions/2026-10-08-agentic-development-plan.md)
preserves the request, evidence and scope limits.

## Research question

**Can an explicit subgrid turbulence model make AthenaK's predictions consistent
across resolution and changing AMR levels, while correctly accounting for energy
transfer, heating and mixing?**

Two connected outcomes are possible: a reusable numerical capability and a physics
result about which turbulent observables that capability recovers. The box tests
should establish where the model helps, where numerical dissipation dominates, and
where its assumptions fail. Agreement with a spectral slope alone cannot answer
those questions.

## What adding subgrid turbulence means

The model represents the effect of fluctuations below a stated filter scale on the
resolved flow. It does not reconstruct the missing eddies. Numerical dissipation
already supplies an implicit treatment of unresolved motions in an ordinary
shock-capturing calculation; an explicit model must be evaluated alongside that
treatment. [Schmidt-Brückner's 2025 review](https://doi.org/10.1007/s41115-025-00023-1),
retrieved 2026-10-06.

| Candidate direction | What it supplies | What needs validating |
|---|---|---|
| Algebraic eddy viscosity | A simple stress closure for unresolved transport | Whether it adds useful transfer or excessive damping on top of the solver |
| Evolved SGS kinetic energy | An unresolved energy reservoir, stress/pressure and dissipation closures | Conservative exchanges, compression response, transport and positivity |
| Structural or mixed closure | An approximation to filtered stresses and transfers from resolved gradients | Filter dependence, transfer direction and stable evolution |
| MHD extension | Unresolved magnetic stresses and electromotive force, alongside velocity effects | Induction, kinetic–magnetic transfer and magnetic-field constraints |

These were alternatives in the original 2026-10-06 framing. The 2026-10-08
implementation-first decision selects the Schmidt–Federrath HD model as the
immediate baseline; the other directions remain possible later comparisons or extensions.
Schmidt & Federrath (2011) provide a compressible kinetic-energy model and a mixed
closure tested on filtered isothermal and adiabatic boxes. Vlaykov et al. and Grete
et al. (2016) provide compressible MHD structural closures and filtered-data tests.
Their work supplies foundations for the HD/MHD branches; neither establishes
performance in the operator's AthenaK setup.
[HD model](https://arxiv.org/abs/1010.4492),
[MHD derivation](https://arxiv.org/abs/1606.02064),
[MHD tests](https://arxiv.org/abs/1606.01573); retrieved 2026-10-06.

For the energy-reservoir branch, a normalized, nonnegative spatial filter and
density weighting define a specific unresolved kinetic energy

$$
k_{\rm sgs}=\frac{1}{2}\left(\widetilde{|\boldsymbol u|^2}
-|\widetilde{\boldsymbol u}|^2\right),\qquad
\widetilde f=\frac{\overline{\rho f}}{\bar\rho}.
$$

In a Newtonian HD formulation the represented total energy density is

$$
E_{\rm represented}=E_{\rm internal}+E_{\rm kinetic,resolved}
+\bar\rho k_{\rm sgs}.
$$

The implementation must define which stored variable contains each reservoir.
SGS production transfers energy from resolved motions; SGS dissipation transfers
it to internal energy. Neither is extra external injection. The isotropic part of
the SGS stress and compression work must also be consistent. Existing conservative
shock heating must be retained without adding it a second time. Magnetic and
relativistic branches require their own energy accounting.

## Why AMR is central

Refinement changes the filter scale and therefore changes what counts as resolved
or subgrid energy. Conservative averaging of density and momentum does not preserve
resolved kinetic energy by itself, because kinetic energy is quadratic in momentum.
Schmidt et al. (2014) address resolved/SGS energy compensation and test adaptively
refined LES, including forced turbulence. [Paper, Section 2](https://arxiv.org/html/1309.3996),
retrieved 2026-10-06.

For restriction of equal-volume fine cells to one coarse cell, the kinetic-energy
difference is

$$
\delta K=\frac{1}{N}\sum_n\frac12\rho_n|\boldsymbol u_n|^2
-\frac{|\boldsymbol m_c|^2}{2\rho_c}
=\frac{1}{2N}\sum_n\rho_n|\boldsymbol u_n-\boldsymbol u_c|^2\geq0,
$$

where $\rho_c$ and $\boldsymbol m_c$ are the averaged density and momentum, and
$\boldsymbol u_c=\boldsymbol m_c/\rho_c$. A compatible HD energy-reservoir model
can assign this difference to SGS energy during restriction. This preserves the
represented kinetic-energy sum while retaining the averaged internal energy.
The [session record](sessions/2026-10-06-digest.md) reconstructs and checks the identity.

That identity does **not** prove that every lost velocity variation is turbulence.
Coherent shear, shocks and bulk gradients need separate attention. Refinement is
also underdetermined: coarse SGS energy does not specify the fine velocity field.
The proposed method must state its reconstruction assumptions, prevent negative
energy, and record any nonconservative corrections. Repeated coarsen/refine cycles
should not manufacture heat, kinetic energy or turbulence. A universal inertial-range
rescaling should be tested rather than presumed valid in shocks, cooling interfaces
or magnetically anisotropic flow.

## What testing the boxes properly could mean

This is a proposed validation design, not a run specification or a floor.

| Test | Comparison | Question answered |
|---|---|---|
| A priori closure test | Filter a better-resolved reference at several scales; compare exact filtered stresses/transfer with the model | Does the closure predict the missing term, including its magnitude and sign? |
| A posteriori uniform boxes | Evolve matched boxes with and without SGS; vary resolution and check against filtered references | Does the model improve the evolved statistics? |
| Static refinement | Matched uniform and fixed-hierarchy boxes | Do coarse–fine interfaces introduce bias before regridding is added? |
| Dynamic AMR | Matched physics with changing hierarchies; deliberately exercise refinement and restriction | Are regridding errors controlled and energy partitions consistent? |
| Cooling and multiphase branch | Matched cooling boxes with closure, resolution and mesh changes | Are heating, mixing and phase statistics improved? |

Proposed controls and diagnostics:

- **Driving:** match forcing band, correlation time and compressive/solenoidal
  content; measure actual injected work. Matching an input amplitude or seed alone
  does not guarantee identical forcing histories or energy injection on different meshes.
- **Thermodynamics:** separate isothermal from energy-evolving runs. An isothermal
  box cannot validate deposited heat. In a driven adiabatic box, dissipative heating
  can change sound speed and Mach number; use a declared time window or an explicitly
  stated heat-removal prescription rather than assume stationary thermodynamics.
- **Observables:** measure spectra over justified common scales, structure functions,
  density/velocity statistics, resolved/SGS energy, transfer and dissipation budgets.
  For cooling runs add temperature–density distributions, cold-gas mass and exchange
  rates, and spatial heating/mixing diagnostics with explicit definitions.
- **Statistics:** use independent realizations and time blocks long enough to estimate
  correlated sampling uncertainty. Do not demand cell-by-cell agreement after
  chaotic trajectories diverge, or use one realization to establish improvement.
  Establish convergence of the reference on the reported scales; a better-resolved
  run is not automatically an exact solution or a direct numerical simulation.
- **AMR analysis:** weight active leaf cells by physical volume or mass as appropriate;
  exclude covered coarse cells and ghosts. Use a declared common filter for comparisons.
  Interpolation onto a fine plotting grid does not supply unresolved information.
- **Cost:** compare error at matched compute cost as well as matched nominal finest
  resolution. Include refined-volume fraction, cell updates, timestep restrictions,
  closure overhead and output cost; an AMR hierarchy may cover most of a box.
- **Numerical effects:** fix the numerical method within matched comparisons, then
  check transferability to another solver/reconstruction choice. Record floors,
  limiters and energy corrections instead of letting them hide closure failures.

A relevant caution is already experimental precedent: Grete et al. (2017) found
that SGS effects in their decaying MHD tests depended strongly on the explicit
filter; models evaluated directly at the grid scale had no statistical influence
in that setup. Their result motivates checking the effective filter and numerical
damping rather than expecting adding SGS terms to guarantee improvement.
[A posteriori comparison](https://arxiv.org/abs/1703.00858), retrieved 2026-10-06.

## Directions and possible ceilings

All are exploratory ambitions, not rankings or agreed deliverables.

1. **A conservative AthenaK capability and benchmark.** Resolve the interaction
   among closure, filter, solver and AMR transfer. Produce a reproducible picture of
   the regimes where resolved-plus-SGS predictions are stable across mesh choices,
   including regimes where an explicit model does not help.
2. **A turbulence–cooling closure.** Ask whether unresolved turbulent energy,
   dissipation location and scalar/enthalpy transport alter condensation or cold-gas
   survival. Momentum closure alone does not specify thermal mixing. Nonlinear
   radiative losses evaluated at a coarse mean need not equal filtered losses;
   additional thermodynamic closure may be required. This connects to
   [cascade under cooling](../cascade-under-cooling/subject.md).
3. **Joint refinement and modelling.** Choose where to refine according to estimated
   error in an observable and the unresolved energy fraction, then measure the
   resulting error-versus-cost frontier. Compare with established refinement
   criteria in [the AMR subject](../amr-criteria-multiphase-turbulence/subject.md).
   Whether an SGS indicator itself selects the right cells is a question to test.
4. **Use AMR as a filtering laboratory.** Explore whether overlapping parent/fine
   representations can constrain scale-dependent transfer or closure coefficients
   during evolution. This is an ambitious diagnostic idea: mesh selection biases the
   available samples, and learning a closure only from refined regions could fail
   elsewhere. Dynamic filtering methods are established; novelty is unassessed.
5. **Observable closure.** Link resolved and unresolved turbulent velocities to
   synthetic line profiles, and test whether those predictions remain consistent
   across AMR hierarchies. Energy alone fixes neither anisotropy nor velocity tails;
   density, temperature and tracer weighting would need to be modelled and checked.

## Scope and parked branches

The immediate implementation baseline is the published Schmidt–Federrath HD model.
The proposed first test sequence uses uniform grids before refinement. MHD and
cooling remain possible later branches; driven/decaying setups and detailed run
parameters have not been finalized.

- Park a universal closure spanning HD, MHD and relativistic fluids until the
  filtered equations and transfer budgets for each branch are specified.
- Park unconstrained learned closures until reference coverage and stable,
  conservative evolution tests exist; keep them as a possible later comparison.
- Park claims of a universal spectral law or automatic compute savings; both need
  a stated physical regime and measured comparisons.

## AthenaK evidence and limits

The upstream README identifies Kokkos and block-based AMR. The inspected upstream
pages include a turbulent-box problem, hydro input and stochastic forcing driver
with hydrodynamic and MHD task integration. These are starting infrastructure to
examine, not proof of SGS support or AMR-forcing correctness.
[AthenaK](https://github.com/IAS-Astrophysics/athenak),
[forcing driver](https://github.com/IAS-Astrophysics/athenak/blob/main/src/srcterms/turb_driver.cpp);
retrieved 2026-10-06.

Exact-name searches did not establish an existing AthenaK SGS implementation.
Only selected upstream pages were inspected in the 2026-10-06 digest; no exhaustive
repository, branch, fork or unpublished-work audit was performed. The later local
inspection identified a candidate checkout, but the operator's implementation
target remains unconfirmed. No claim of an unoccupied niche follows from this search.

The [reference log](references.md) records retrieved sources, reading depth and
coverage. No simulation, code build or performance measurement was performed.

## Related work in the wiki

- [AMR criteria for multiphase turbulence](../amr-criteria-multiphase-turbulence/subject.md)
  — related refinement and turbulent-box benchmarking subject; its record also
  mentions AMR with subgrid turbulence models.
- [Cascade under cooling](../cascade-under-cooling/subject.md) — includes a parked
  cooling-regulated subgrid turbulence closure direction.
- [Numerical dissipation density bias](../numerical-dissipation-density-bias/subject.md)
  — related diagnostics for numerical dissipation in turbulent simulations.
