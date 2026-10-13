---
title: AthenaK SGS turbulence — agentic development plan
class: living
created: 2026-10-08
last_touched: 2026-10-08
---

# AthenaK SGS turbulence: agentic development plan

Implement the Schmidt–Federrath hydrodynamic SGS model, establish its numerical
correctness and turbulent-box behaviour, then address refinement as a separate
extension. This plan defines work packages, dependencies, test gates and evidence
requirements. It authorizes no claim that code exists or tests have passed.

Companion documents:

- [Test catalogue](test-matrix.md): individually identifiable checks and their oracles.
- [Agent runbook](agent-runbook.md): assignments, handoffs, review and execution rules.
- [Living research subject](subject.md): scientific motivation and broader directions.
- [Implementation-first decision](sessions/2026-10-08-implementation-first.md).

## Implementation status on 2026 10 08

Restricted HD code and independent verification are now recorded in [implementation progress](implementation-progress.md) and the [execution ledger](coverage-ledger.md). M0–M9 scientific release is incomplete. Original planning statements and checklists below are preserved as their historical record.

## 1. Decisions, assumptions and completion boundaries

**Operator decision:** implement the published Schmidt–Federrath method in AthenaK
before extensions. The original idea includes proper turbulent-box testing and AMR.

**Proposed implementation sequence:** Newtonian HD, ideal-gas energy equation,
Cartesian uniform grids, initially periodic boundaries; establish serial CPU
correctness, then decomposition/MPI, restart and an available GPU backend. Start
with RK2, PLM and HLLC if the numerical-contract audit supports them; this is a
starting configuration, not a promise of support for every solver combination.
Pseudo-isothermal boxes are energy-evolving runs, not AthenaK's strict isothermal EOS.

**Not yet settled:** target checkout/branch, complete support matrix, stencil,
pressure treatment in the Riemann solver, positivity mechanism, exact-zero SGS
initialization policy, timestep coefficients, reference dataset and compute budget.
M0–M1 settle these in a recorded contract. They do not block writing this plan.

Candidate inspected on 2026-10-08:
`/Users/hitesh/hitesh/git/athenak_nyu`, `gotham_development`,
`07c529202befbeb6ef5e694d95286929618ff7c1`. The operator has not confirmed it as
the development target. Recheck its state before implementation. Existing local
changes belong to their author; do not reset, stage or overwrite them.

Three distinct completion claims:

| Claim | Required evidence | What it does not establish |
|---|---|---|
| Numerically implemented | Equations mapped to code; transport, energy exchange, positivity, restart and supported execution paths pass deterministic gates | Accuracy of the closure in a new turbulent regime |
| Baseline HD validated | Filtered-reference checks and evolved uniform-grid boxes with budgets, convergence, uncertainty and limitations | Universality, automatic savings or AMR correctness |
| AMR validated | Separate inter-level contract, static/dynamic hierarchy tests and statistical comparisons | New MHD, cooling or thermodynamic closures |

All checkboxes start unchecked. A proposal is not an accepted numerical contract;
a test existing in code is not a passing result. The catalogue is comprehensive
for these declared contracts. New supported physics or execution paths must add
their own cases; a finite catalogue cannot literally enumerate every possible flow.

The initial plan and runbook contain 212 unchecked TODOs. The catalogue identifies
230 test families: 151 baseline, 28 conditional support-path, 11 characterization,
32 later AMR and 8 future-physics families. Required interaction combinations are
listed separately in the catalogue; one family can generate multiple concrete runs.

## 2. Scientific contract to freeze before solver integration

Use $\rho$ for filtered density, $\boldsymbol u$ for Favre velocity,

$$
K\equiv K_{\rm sgs},\qquad k=K/\rho,\qquad
E=e_{\rm int}+\frac12\rho|\boldsymbol u|^2,\qquad E_{\rm represented}=E+K.
$$

Here both $K$ and $e_{\rm int}$ are energy densities. Primitive SGS storage
may contain $k$, but the conservative transport equation is for $K$.

For $A_{ij}=\partial_j u_i$, define

$$
S_{ij}=\frac12(A_{ij}+A_{ji}),\quad d=A_{ii},\quad
S^*_{ij}=S_{ij}-\frac13d\delta_{ij},\quad Q=\sum_{i,j}A_{ij}^2.
$$

For $Q>0$, the algebraically normalized mixed closure is

$$
\tau_{ij}=2C_1\Delta\sqrt{2\rho K}\,S^*_{ij}
-2C_2K\frac{\sum_l A_{il}A_{jl}}{Q}
-\frac23(1-C_2)K\delta_{ij}.
$$

$$
\operatorname{tr}\tau=-2K,\qquad P_{\rm sgs}=\frac23K,\qquad
\tau^*_{ij}=\tau_{ij}+P_{\rm sgs}\delta_{ij},\qquad
\Sigma=\tau_{ij}S_{ij}=\tau^*_{ij}S_{ij}-P_{\rm sgs}d.
$$

$$
\mathcal{E}_{\rm diss}=\rho\epsilon
=C_\epsilon\frac{K^{3/2}}{\sqrt\rho\Delta},\qquad
\kappa= C_\kappa\Delta\sqrt{\rho K},\qquad
\boldsymbol J=-\kappa\nabla k.
$$

Baseline coefficients: $C_1=0.02$, $C_2=0.7$, $C_\epsilon=1.5$,
$C_\kappa=0.65$. Neglect SGS gravitational production and SGS thermal
pressure-dilatation as in this baseline. These equations and coefficient choices
come from [Schmidt & Federrath (2011), §§2–3](https://arxiv.org/html/1010.4492).

The closed equations, excluding external sources for clarity, are

$$
\partial_t\rho+\nabla\cdot(\rho\boldsymbol u)=0,
$$

$$
\partial_t(\rho u_i)+\partial_j(\rho u_i u_j+P\delta_{ij}-\tau_{ij})=0,
$$

$$
\partial_t E+\nabla\cdot\left[(E+P)\boldsymbol u-\boldsymbol u\cdot\tau\right]
=-\Sigma+\mathcal{E}_{\rm diss},
$$

$$
\partial_t K+\nabla\cdot(K\boldsymbol u+\boldsymbol J)
=\Sigma-\mathcal{E}_{\rm diss}.
$$

Their sum supplies the independent budget oracle:

$$
\partial_t(E+K)+\nabla\cdot\left[(E+K+P)\boldsymbol u
-\boldsymbol u\cdot\tau+\boldsymbol J\right]=0.
$$

Forcing adds its actual resolved work; radiative losses, if later supported, add
their measured sink. SGS diffusion transports represented total energy. It is not
an internal conversion that should automatically be subtracted from $E$.
Production and dissipation, in contrast, must appear as paired exchanges.

### Required numerical decisions

- [ ] Specify $\Delta$: begin with cubic cells and $\Delta=\Delta x$.
  Reject or separately qualify stretched cells until their filter convention is tested.
- [ ] Define the $Q=0$ branch. A proposed exact-zero convention is isotropic
  stress $\tau=-2K I/3$. The normalized nonlinear tensor has a direction-dependent
  limit as the gradient vanishes; this is a numerical convention, not a continuous
  limit derived from the formula. Avoid an untracked denominator epsilon that
  destroys the trace identity.
- [ ] Specify small positive $Q$, underflow, overflow and scaling-safe norm evaluation.
- [ ] Define zero and tiny $K$: no NaNs, negative square roots or hidden energy injection.
  Exact $K=0$ gives zero production in the $K$ formulation. Decide explicitly
  between a documented initial seed and a derived alternative evolution variable;
  do not introduce a repeated floor and call it physical production.
- [ ] Choose face/cell gradient stencils and interpolation consistent with the ghost budget.
- [ ] Derive the discrete relation between stress work, resolved kinetic work,
  SGS production and internal heating. Global $E+K$ conservation alone does not
  prove the kinetic/thermal partition is correct.
- [ ] Choose isotropic-pressure handling: integrate it into a consistent hyperbolic
  solver, or use a documented split stress update with demonstrated stability.
  Do not merely change a pressure argument while leaving incompatible energy fluxes.
- [ ] Derive hyperbolic wave bounds and explicit stress/diffusion/source restrictions.
  SGS pressure has an isotropic compression response with exponent $5/3$,
  independently of thermal $\gamma$. A thermal sound-speed check alone is insufficient.
- [ ] Select a positivity strategy for $K$, density and internal energy, including
  bounded backscatter. A paired energy transfer can conserve the sum and still make
  either reservoir negative; check available energy and stage weights.
- [ ] Decide whether dissipation is stage-consistent explicit integration, bounded
  integration or a split exact local solve. Verify the resulting temporal order;
  repeatedly applying an exact solve with RK weights is not automatically correct.
- [ ] State supported dimensional reductions: retain the three-component stress trace
  in 1D/2D tests; do not silently replace $2/3$ by $2/D$.
- [ ] Freeze the numerical contract and its deviations from the continuum model.

Useful independent analytic oracles:

$$
K(t)=\frac{K_0}{[1+C_\epsilon\sqrt{K_0}\,t/(2\sqrt\rho\Delta)]^2}
\quad\text{(constant density, dissipation only)},
$$

$$
e_{\rm int}(t)=e_{\rm int,0}+K_0-K(t),
$$

$$
K\propto\rho^{5/3},\qquad k\propto\rho^{2/3}
\quad\text{(homogeneous isotropic compression, no diffusion/dissipation)}.
$$

These are derived test oracles, not substitute closure calibrations.

## 3. AthenaK integration map and design constraints

The following paths describe the inspected candidate revision. Revalidate after
target selection; upstream layouts may differ.

| Area | Candidate files | Required work |
|---|---|---|
| State and configuration | `src/hydro/hydro.hpp`, `hydro.cpp` | Optional SGS state, coefficients, indices, counts, allocation and explicit support checks |
| Advection | `hydro_fluxes.cpp`, reconstruction helpers | Consistent advected $K/\rho$, upwind mass flux, nonnegative reconstructed states |
| Stress and transport | New focused `src/hydro/hydro_sgs.*` if justified; existing diffusion files as patterns | Device-callable closure, gradients, conservative face fluxes and stage data |
| Integration | `hydro_tasks.cpp`, `hydro_update.cpp` | Paired stage sources, valid primitive state, no stale gradients, correct dependencies |
| Robustness | `hydro_fofc.cpp`, `src/eos/ideal_hyd.cpp` | SGS-aware fallback, admissibility checks, explicit correction accounting |
| Timestep | `hydro_newdt.cpp`, `src/mesh/mesh.cpp` | Local constraints and global/rank minimum; correct first step and restart step |
| Boundary/communication | `src/bvals/`, hydro task receive/send calls | All needed state/ghosts, one shared face flux, correct extents and BC semantics |
| Output/restart | `src/outputs/`, restart input path | SGS names, history budgets, state and parameter compatibility |
| Initial conditions | `src/pgen/turb.cpp`, dedicated small verification pgen if needed | Explicit seed and energy convention; no reinitialization on restart |
| Build/tests | `src/CMakeLists.txt`, `inputs/`, `tst/scripts/hydro/`, existing readers | Focused tests using existing harness and Kokkos backend |
| Later refinement | mesh refinement, restriction/prolongation, flux correction, load balance | Separate level-transfer contract and tests after HD gates |

### State layout

Prefer conservative SGS storage alongside communicated hydro state if the audit
confirms every lifecycle path can carry it correctly. Append the slot after true
passive scalars, and distinguish the user-specified passive count from total
advected/storage counts. Record a single SGS index; do not scatter magic offsets.
An alternative separate array must justify its extra communication/restart lifecycle.

- [ ] Audit every `nhydro+nscalars` assumption, scalar loop, allocation and output label.
- [ ] Preserve existing passive indices and meanings. The inspected cooling code
  assumes metallicity is the first passive scalar at `nhydro`.
- [ ] Keep ordinary thermal $E$ semantics or explicitly document a fully audited
  alternative. Never let EOS thermal pressure accidentally include $K$.
- [ ] Preserve the disabled feature's layout, numerical path and old input behaviour.
- [ ] Avoid an extensible closure factory, generic turbulence framework or new scheduler.
  One optional model and focused kernels are sufficient for this deliverable.

### Stage and fallback coupling

The inspected task graph computes fluxes, corrects fine/coarse fluxes, updates
conserved variables, applies sources, exchanges state, applies BCs/prolongation,
then recovers primitives and computes the next timestep. Sources are documented
to use the old stage primitives because conserved state has already been updated.

- [ ] Write the exact execution order on paper, including forcing/user hooks and
  Kokkos kernel ordering. Audit every supported RK stage.
- [ ] Compute closure/source data from the intended stage and avoid reading partially
  updated conserved values as if they were the old state.
- [ ] Decide what is stored versus cheaply recomputed; measure before adding persistent
  3×3 tensor arrays over the whole domain.
- [ ] Trace first-order flux correction (FOFC). Its trial state currently covers
  `nhydro`, and its fallback replaces hydro fluxes; SGS admissibility/stress work
  cannot be assumed to survive this path automatically.
- [ ] Use the final applied face stresses when constructing compatible work/production.
  Never keep a source computed from a stress that fallback removed.
- [ ] Ensure every shared face uses the same limiter and flux for both neighbours.
- [ ] Include SGS effects in the first timestep, subsequent stages and restart timestep.
- [ ] Audit task dependencies rather than depending on incidental serial execution.

## 4. Milestones and dependency graph

| Milestone | Deliverable | Dependencies | Exit gate |
|---|---|---|---|
| M0 | Confirmed workspace, reproducible untouched baseline | Target choice | Baseline build/test report and preserved local changes |
| M1 | Reviewed equation and numerical contract | M0 | Independent science/numerics review; support matrix frozen |
| M2 | Independent closure oracle and device kernels | M1 | Algebra, units, edge-state and backend kernel tests |
| M3 | SGS state lifecycle and disabled-path compatibility | M1 | Advection, scalar coexistence, allocation, restart scaffolding |
| M4 | Fully coupled uniform-grid HD solver | M2, M3 | Stress/work/exchange/transport tests and stage-order checks |
| M5 | Positivity, fallback, timestep and hard-flow robustness | M4 | No hidden corrections; supported robustness cases pass |
| M6 | Parallel, restart and diagnostics qualification | M5 | Decomposition/restart/budget and available GPU checks |
| M7 | Filtered-reference closure assessment | M2, diagnostics from M6 | Independent filter/oracle review; signed accuracy report |
| M8 | Evolved turbulent-box validation | M6, M7 | Statistical campaign, converged budgets, documented limits |
| M9 | Baseline HD release package | M8 | Reproducible examples, review, compatibility and results |
| A0 | AMR numerical contract and isolated remap tests | M9 | Chosen transfer rule, invariants and positivity proved/tested |
| A1 | Static-refinement implementation and validation | A0 | Coarse/fine transport and conservation gates |
| A2 | Dynamic AMR, restart and statistical qualification | A1 | Regridding/load-balance/cost tests and comparison report |
| X | Additional physics or closure research | M9; AMR where relevant | A separate approved scientific/numerical contract |

Critical path: `M0 → M1 → (M2 and M3) → M4 → M5 → M6 → M8 → M9`.
M7 can start when M2 and sufficient validated reference tools exist, but M8 cannot
claim closure accuracy without it. AMR starts after M9; X is a separate programme.
Agent roles can overlap in preparation, not by bypassing these gates.

## 5. Milestone TODO ledger

### M0 — workspace and baseline

- [ ] Confirm the target checkout and base branch once; do not repeatedly ask.
- [ ] Read target/parent instructions and inventory tracked, staged and untracked changes.
- [ ] Create an isolated branch/worktree where appropriate, using `codex/` names.
- [ ] Record exact commit, compiler, Kokkos, architecture, precision and build options.
- [ ] Build serial CPU debug and release configurations in separate build directories.
- [ ] Run selected existing hydro/scalar/restart tests before SGS edits.
- [ ] Run a short existing turbulent-box input; save actual parameters and baseline output.
- [ ] Check the test harness's cleanup behaviour before running it in a shared workspace.
- [ ] Classify pre-existing failures and environmental limitations; preserve raw logs.
- [ ] Record baseline memory use, updates/second and representative per-step timing.
- [ ] Gate: another agent can reproduce the baseline and identify any pre-existing failures.

### M1 — equations, support and numerical design

- [ ] Pin the paper version and independently check equations, signs, units and coefficients.
- [ ] Produce the continuum-to-discrete term map, including stress work and compression.
- [ ] Select state layout and migration/restart policy; enumerate every affected loop.
- [ ] Choose the default supported solver/integrator and an unsupported-combination policy.
- [ ] Specify boundaries, dimensions, filter width and stretched-cell handling.
- [ ] Derive source, advection, diffusion/stress and wave-speed restrictions.
- [ ] Specify initialization, positivity, bounded backscatter and correction diagnostics.
- [ ] Define deterministic/stochastic acceptance thresholds before generating release results.
- [ ] Review the design with a science reviewer and an integration reviewer.
- [ ] Gate: no unresolved energy sign, variable unit, wave-speed or positivity ambiguity.

### M2 — mathematical oracle and kernels

- [ ] Write an independent small reference calculation from the locked equations.
- [ ] Derive hand cases for zero gradient, strain, compression, shear and rotation.
- [ ] Implement the smallest device-callable closure helpers; do not duplicate production
  formulas in several places.
- [ ] Use deterministic random tensors plus hand cases to test symmetry and trace.
- [ ] Test coefficient endpoints, scaling and extreme positive inputs.
- [ ] Check $Q=0$, near-zero norms, zero $K$, tiny density and overflow boundaries.
- [ ] Verify CPU/Kokkos outputs against the independent oracle; GPU when available.
- [ ] Deliberately change a sign, factor or contraction in a temporary test branch and
  ensure the relevant test fails; revert the deliberate defect immediately.
- [ ] Gate: closure tests detect wrong signs, factors, dimensions and contractions.

### M3 — state lifecycle

- [ ] Add one optional SGS state with explicit activation/configuration checks.
- [ ] Initialize conservative/primitive state with documented energy bookkeeping.
- [ ] Carry it through reconstruction, advection, stage copies and ghost exchange.
- [ ] Preserve true scalar count, names, indices and meanings.
- [ ] Add output fields and restart serialization using existing machinery.
- [ ] Reject incompatible old/new restart combinations clearly; test old restarts with SGS off.
- [ ] Verify no SGS allocation/kernel cost when disabled beyond a small configuration branch.
- [ ] Exercise one and many MeshBlocks, 1D/2D test reductions and periodic crossings.
- [ ] Gate: passive SGS transport works; no claim of full turbulence closure yet.

### M4 — coupled fluxes and sources

- [ ] Add face stress to momentum and its work to resolved energy consistently.
- [ ] Include turbulent pressure once; avoid duplicate isotropic compression work.
- [ ] Add conservative SGS diffusion of specific energy, with density-weighted coefficient.
- [ ] Pair production and dissipation between $E$ and $K$ at each applicable stage.
- [ ] Verify compatible resolved kinetic/thermal accounting, not merely total cancellation.
- [ ] Implement the selected temporal treatment and prove/test order on smooth cases.
- [ ] Integrate timestep estimates into the global minimum, including the initial step.
- [ ] Instrument actual integrated applied exchanges and fluxes.
- [ ] Run isolated decay, pressure, shear, diffusion and combined manufactured cases.
- [ ] Gate: conservative coupled solver and expected smooth-case convergence.

### M5 — admissibility and robustness

- [ ] Add SGS admissibility to relevant predictors and source-update checks.
- [ ] Make fallback fluxes and closure work consistent on corrected shared faces.
- [ ] Test negative-transfer states with insufficient SGS energy and intense forward transfer.
- [ ] Track any correction, its energy/momentum change and its location/count.
- [ ] Demonstrate stability over a declared CFL range and timestep convergence.
- [ ] Exercise shocks, contacts, strong rarefactions, low density and high SGS pressure.
- [ ] Test reconstruction-generated negative SGS states and ghost/corner stencil access.
- [ ] Ensure unsupported modes fail at initialization with an actionable message.
- [ ] Gate: robustness cases remain admissible; corrections are small/quantified or
  the support restriction is tightened before moving forward.

### M6 — execution, persistence and budgets

- [ ] Test block layouts, MeshBlockPack sizes, MPI rank counts and available threads/backends.
- [ ] Compare short deterministic runs before chaotic divergence with precision-aware norms.
- [ ] Check uninterrupted versus restarted runs, including forcing random state and counters.
- [ ] Verify new output/readers and history integration against independently summed fields.
- [ ] Confirm periodic conservation and open-boundary flux accounting where supported.
- [ ] Qualify the available GPU backend; record unavailable CUDA/HIP coverage honestly.
- [ ] Check style/build diagnostics and applicable existing regressions.
- [ ] Benchmark disabled/on performance and memory using matched jobs and repeated timing.
- [ ] Gate: runnable qualification bundle with no unaccounted source or flux residual.

### M7 — a priori closure tests

- [ ] Select or generate a better-resolved reference and establish its useful converged scales.
- [ ] Freeze filter type/width, Favre weighting, derivative scheme and sampling windows.
- [ ] Independently compute exact filtered $K$, stress and production.
- [ ] Apply the model to filtered fields at the same stated scale.
- [ ] Test several filter widths, forcing types and Mach/thermodynamic regimes.
- [ ] Report signed correlations, biases, magnitude errors, conditional statistics and
  forward/backscatter fractions with sampling uncertainty.
- [ ] Hold published coefficients fixed for the baseline; reserve any refit for a separate study.
- [ ] Review filtering artefacts, density weighting and effective versus nominal grid scale.
- [ ] Gate: trustworthy closure-accuracy report, even if it identifies limitations.

### M8 — a posteriori turbulent boxes

- [ ] Run pilot boxes to estimate saturation, timestep, correction frequency and total cost.
- [ ] Freeze the campaign protocol before looking at final comparison outcomes.
- [ ] Run matched model-off/on cases, multiple resolutions and forcing mixtures.
- [ ] Include pseudo-isothermal conservation tests and ordinary energy-evolving heating tests.
- [ ] Include decay after driving stops and initial-seed sensitivity.
- [ ] Measure actual forcing work; do not infer it from the input injection setting.
- [ ] Compare common filtered resolved-plus-SGS quantities and stated resolved observables.
- [ ] Estimate uncertainty with independent realizations and correlated-time-aware sampling.
- [ ] Repeat selected cases at smaller timestep and a second qualified numerical method.
- [ ] Compare error versus cost as well as resolution; quantify limiter and numerical effects.
- [ ] Gate: results distinguish numerical correctness, reproduction and scientific benefit.
  A null improvement is a valid result if correctly measured and reported.

### M9 — baseline HD release

- [ ] Provide one minimal example, one conservation verification and the box-campaign recipe.
- [ ] Document equations, conventions, coefficients, support restrictions and restart policy.
- [ ] Bundle input hashes, code/build manifests, raw metrics and scripts for every release claim.
- [ ] Run independent checkout reproduction and final science/code review.
- [ ] Fix material findings; rerun impacted gates rather than all expensive jobs by reflex.
- [ ] Keep long scientific campaigns outside the fast CI suite.
- [ ] Prepare a focused PR or local patch series when repository workflow permits.
- [ ] Record exactly what was committed, pushed, reviewed and merged as separate statuses.
- [ ] Gate: reproducible supported HD capability, with limitations and measured results.

### A0–A2 — refinement extension after baseline release

- [ ] Choose and document an inter-level energy compensation rule; do not silently
  identify every lost resolved gradient with turbulence.
- [ ] Separate conserved represented energy from desired thermal/kinetic/SGS partition.
- [ ] Prove/test restriction identities and positivity-constrained prolongation.
- [ ] Audit coarse ghosts separately from permanent cell remapping.
- [ ] Include SGS stress and transport in refluxing and the timestep logic.
- [ ] Test fixed nested grids before enabling regridding.
- [ ] Force repeated refine/coarsen cycles in uniform, shear, shock and random fields.
- [ ] Check regridding followed by load balance, MPI exchange and restart.
- [ ] Compare changing hierarchies with uniform runs using a common filter and leaf weights.
- [ ] Quantify active-volume fraction, finest-level occupancy, cell updates and cost.
- [ ] Gate: separate AMR review/report; do not relabel ordinary passive restriction as LES AMR.

## 6. Box campaign specification

Proposed sizes are scheduling tiers, not guaranteed allocations. Verify memory and
wall-time from a pilot before launching a full tier. A better-resolved ILES reference
is not automatically DNS.

| Tier | Proposed scope | Purpose |
|---|---|---|
| Smoke | Small deterministic grids; 16³–32³ short boxes | Kernel/lifecycle crashes, budgets, rapid iteration |
| Pilot | 32³ and 64³; one seed per selected forcing mode | Cost, saturation, correction rate, forcing check |
| Validation | 64³, 128³, 256³; model off/on; solenoidal/compressive/mixed; initially 3 independent seeds | Resolution trends and uncertainty; extend seed/time sample if confidence is inadequate |
| Reference/reproduction | 512³ or higher only with demonstrated resources and need | Common-filter reference and closer paper comparison |

Do not take the Cartesian product of every option. First run the mandatory
core, then pairwise numerical coverage, then targeted stress tests and the campaign.

- [ ] Choose the box size, density, thermal sound speed, $\Delta$, forcing band,
  injection method, autocorrelation time and mode-mixture convention.
- [ ] Define $T=L/V$ and the Mach number weighting; report measured values.
- [ ] Pilot until stationary kinetic/SGS statistics are demonstrated, then freeze burn-in
  and sampling windows. Proposed start: assess saturation over at least several $T$,
  then sample multiple correlation times; do not equate a fixed duration with equilibrium.
- [ ] Use published coefficients without retuning in the reproduction lane.
- [ ] Distinguish replayed forcing from merely repeated seed, especially across dt/rank changes.
- [ ] Use a separate ordinary-$\gamma$ run to expose heating that a large internal reservoir masks.
- [ ] Define the comparison filter and compute both volume- and mass-weighted quantities
  where physically appropriate; do not mix their Mach/energy definitions.
- [ ] Report $\langle K\rangle$, resolved kinetic energy, internal energy, SGS Mach,
  production, dissipation, diffusion flux, injected work and every correction budget.
- [ ] Include spectra, structure functions, density PDFs, compressive/solenoidal fractions
  and conditional dissipation statistics as declared diagnostics, with scale limits.
- [ ] Estimate uncertainty in fitted resolution scalings; do not hard-code literature
  exponents as exact truth for a different solver/filter combination.

The original evolved-model tests used stochastic forcing and an energy equation
with $\gamma=1.001$ to approximate isothermality while checking conservation.
That is a reproduction lane, not a heating-only verification strategy.
[Schmidt & Federrath (2011), §4](https://arxiv.org/html/1010.4492).

## 7. Acceptance thresholds and evidence policy

Use two normalized budget measures. With boundary flux $F_{\rm out}$ positive
outward and prescribed external work/source $W_{\rm ext}$,

$$
R_E=\mathcal{E}(t)-\mathcal{E}(0)-W_{\rm ext}(t)
+\int_0^t F_{\rm out}\,dt-\Delta\mathcal{E}_{\rm recorded\ corrections},
\quad \mathcal E=\sum_{\rm active}V(E+K).
$$

Normalize once by a declared represented-energy scale and once by the cumulative
absolute transferred/injected energy scale with a fixed small denominator guard.
The latter prevents pseudo-isothermal internal energy from hiding large errors in
the turbulent energy budget. Count correction energy as a reported defect cost,
not something whose subtraction makes a run physically conservative.

Proposed threshold protocol:

- Algebra: scaled double-precision errors near a few hundred machine epsilons for
  moderate inputs; adjust only for explained conditioning, not to rescue a bad factor.
- Smooth PDE/ODE: observed order from at least three refinements in time or space,
  keeping the other error source small; fit the asymptotic range before setting regression bounds.
- Conservative deterministic boxes: start with a double-precision target of relative
  energy residual $10^{-10}$ over a modest short verification run, plus a transfer-budget
  target. This is a design target pending conditioning and baseline characterization;
  not a universal promised precision. Do not silently borrow the paper's number.
- Positivity: no NaN/Inf; $K\geq0$ and thermal/density admissibility after each checked
  stage. Log any allowed tiny roundoff repair and its paired effects.
- Decomposition/backends: short-run norm agreement within characterized floating-point
  differences; turbulent statistics within declared uncertainty, not bitwise chaos matching.
- Performance: publish measured overhead and slowest path; set a resource budget after
  the pilot. There is no promised percent-level overhead or automatic speedup.
- Statistical benefit: predeclare observable error metric, reference uncertainty and a
  meaningful improvement threshold. Failure to improve is not permission to tune
  coefficients against the held-out validation sample.

Every threshold file must state precision, input range, duration, expected order,
normalization, baseline measurement and rationale. Freeze it before release results.
Failed runs remain in the ledger. Review any threshold change as a code/design change.

## 8. Agentic execution strategy

Use the [agent runbook](agent-runbook.md). Roles are responsibilities, not a
requirement to create a separate persistent agent or new chat for each one.

1. Coordinator locks a small work packet with explicit scope, invariant and exit test.
2. Science/numerics reviewer checks the equation mapping before the implementer codes it.
3. Implementer changes one coherent slice and runs its focused verification.
4. Independent verifier runs the oracle and tries the documented failure modes.
5. Integrator reviews the diff, incorporates it, runs affected integration tests and updates evidence.
6. Coordinator closes the packet only after the recorded exit gate passes.

Independent preparation can run in parallel after M1: reference calculations,
read-only lifecycle audit and filter-analysis tools. Assign one owner to shared
hydro/EOS/task files. Isolate overlapping code changes in separate branches/worktrees;
integrate sequentially. Long jobs are dispatched only after pilot validation and a
resource cap; do not repeatedly launch expensive parameter sweeps on failure.

This request produces the plan only. It does not launch agents, modify AthenaK,
schedule future runs or allocate cluster time.

## 9. Risks, stop conditions and fallback decisions

| Risk | Early signal | Required response |
|---|---|---|
| Wrong variable/sign convention | Trace, decay or energy-source test fails | Stop coupled work; fix contract/oracle before more coding |
| Scalar floor hides SGS loss | Total budget jumps at C2P | Implement explicit admissibility/correction policy; do not loosen conservation threshold |
| FOFC removes stress but not production | Errors localize on fallback faces | Make flux/source coupling consistent; rerun targeted shocks and energy tests |
| Missing pressure/timestep response | Large-$K$ wave instability or CFL dependence | Derive/update wave bounds and timestep; restrict unsupported solvers if necessary |
| Zero SGS reservoir never activates | Seed-dependent or identically-zero evolution | Resolve startup formulation; label seed energy explicitly |
| Nonlinear closure backscatter exceeds reservoir | Negative $K$ or thermal energy | Bound applied exchange consistently; quantify changed transfer |
| Numerical damping overwhelms model | No measurable change at grid scale | Measure effective filter/solver effects; report null result before proposing altered filter |
| Apparent improvement from parameter tuning | Same data used for fit and evaluation | Separate calibration lane and held-out assessment |
| GPU/MPI bug masked by serial tests | Pack/rank/backend disagreement | Minimize deterministic reproducer; stop long campaigns |
| False stationary state | Secular thermal/Mach drift or correlated samples | Revise protocol before interpreting statistics |
| AMR manufacture of turbulence/heat | Regrid-only energy partition changes | Fix/test remap contract; separate coherent gradients from turbulence |
| Excessive memory/compute | Pilot exceeds budget | Profile/recompute/reduce diagnostic cadence; shrink campaign transparently |

Stop the affected milestone when its gate fails; unrelated preparation may continue.
Do not substitute a more elaborate closure to conceal a numerical implementation bug.

## 10. Delivery inventory and definition of done

- [ ] Model contract with equation-to-code map and implementation deviations.
- [ ] Focused patch series and supported-configuration matrix.
- [ ] Runnable deterministic verification suite and independent reference calculations.
- [ ] SGS state/output/restart documentation and minimal input examples.
- [ ] Test ledger with every required ID marked pass/fail/blocked/not run and evidence.
- [ ] Filtered-reference and evolved-box reports with uncertainty and full budgets.
- [ ] Performance/memory report on actually tested hardware.
- [ ] Known limitations, unsupported combinations and deferred work.
- [ ] Reproduction instructions from a clean checkout and final independent review.
- [ ] Separate AMR deliverable when A0–A2 are undertaken.

A baseline release can finish without AMR, MHD or cooling support, provided those
limits are enforced and documented. GPU support is claimed only for tested backends.
The larger scientific question remains open even after correct implementation.

## 11. First executable work packet after target confirmation

**Packet P0:** pin the confirmed checkout, preserve its existing changes, build the
untouched serial CPU debug/release configurations, run selected existing hydro and
short turbulence cases, and return a baseline manifest with actual failures and cost.

**Packet P1:** freeze the SGS numerical contract and independently derive the closure,
decay, compression and energy-budget oracles. Return a reviewed term map and support
matrix. Only then start solver code.

No fixed calendar or person-day estimate is supplied before P0: baseline build health,
backend access and turbulence reference cost have not yet been measured.
