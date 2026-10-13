# P4 independent integration review

Reviewed current dirty `codex/schmidt-federrath-hd` implementation on
2026-10-08, based on `07c529202befbeb6ef5e694d95286929618ff7c1`.
Allowed source inspection only; no solver edits, new builds or simulation runs
in this packet. Source changed during review; the coordinator owns integration.
This report is not M4/M5 approval.

## Material findings

1. **Resolved during review: legacy single-precision build regression.**
   The new `sgs_closure_oracle` target originally built whenever BUILD_TESTING,
   while its compile-time assertion requires double precision. CTest defaults
   BUILD_TESTING on, so a legacy single-precision build could fail even with SGS
   disabled. Reported to coordinator. Final readback shows the target gated by
   `BUILD_TESTING AND NOT Athena_SINGLE_PRECISION`. Rebuild evidence remains
   the coordinator/verifier's responsibility.

2. **Resolved in source during review: invisible diffusion limiting.**
   The initial `sgs-limited-integral` recorded suppressed stress production only;
   donor clipping of J was unrecorded. A density-contrast diffusion run could
   heavily clip transport while this diagnostic stayed zero. Reported to coordinator.
   Final source readback shows `sgs-stress-limited-integral` and separate
   `sgs-diffusion-limited-integral`, with face deficit `abs(raw_J-J)/dx` and half
   allocation to each adjacent cell, integrated with Heun quadrature. The half
   allocation counts each periodic shared face once after volume summation.
   These are suppression diagnostics, not lost represented-energy corrections.

3. **Open practical blocker: nonlinear timestep collapses at ordinary smooth
   shear extrema.** The unconditional local estimate `8 C2 k/norm(A)` is
   singular as the gradient goes to zero, and an exact-zero branch skips it.
   Therefore floating-point remnants of a theoretically zero centered derivative
   can make the timestep roughly fifteen orders of magnitude smaller than at
   a bitwise-zero derivative. This is stronger than a generic high-gradient
   stiffness caveat and undermines practical box execution.

   Minimal setup for independent reproduction: 2D shear verification pgen,
   one 10x10 block/domain, x/y in [-0.5,0.5], vx=0.3 sin(2 pi y), k=0.01,
   C2=0.7 and CFL=0.1. Direct evaluation using the same cell-center arithmetic
   gives a centered derivative of approximately 4.163336342344337e-16 at a
   shear extremum. The nonlinear restriction alone then gives dt approximately
   1.858632295689437e-18, compared with ordinary hydro timesteps near 1e-3.
   This is an independently evaluated stencil estimate, not a simulation result.

   In coherent shear, the normalized A A^T/Q shape is constant on either side
   of a sign change, so a worst-case arbitrary-direction derivative bound can
   be grossly overconservative for that perturbation. A reviewed regularization,
   directional stencil bound or controlled finite-step strategy needs separate
   evidence; do not silently insert a norm epsilon or remove the restriction
   while claiming a proven full-closure wave bound. Add an explicit no-progress
   diagnostic and resource-capped regression for this case. Smaller CFL does
   not solve practical usability. This remains an open implementation limitation.

## Coupling inspected and found consistent at source level

- Full SGS traction is added once to momentum flux, and its arithmetic-face
  velocity work once to E flux. SGS diffusion goes only into K flux.
- `SGSRates` uses the stored final limited traction, with the derived compatible
  production from both adjacent velocity jumps. It does not keep raw production
  from a removed or scaled stress.
- `SGSSources` applies one paired exchange to K and E using beta=(1,1/2),
  while cumulative production/dissipation use correct Heun weights=(1/2,1/2).
- Stage rates are assembled from unchanged stage w0 before RKUpdate, so the
  source does not accidentally use partially updated conservative state.
- `SGSTrial` checks the hidden full-dt FE substep at both SSP stages. It includes
  actual momentum and E/K sources, then the final stage state is checked again
  before EOS conversion. Ordinary scalar floors cannot touch the separate SGS slot.
- On a failing trial, the globally reduced failure count makes all ranks rebuild
  donor-cell hydro fluxes, then recompute SGS traction/work/production/diffusion.
  Failure after the second trial exits explicitly instead of silently flooring.
  AthenaK presently has one MeshBlockPack for the mesh/rank, so this is globally
  selected for the current architecture, not independently chosen per block.
- Backscatter caps use full dt, the minimum adjacent donor K and the same face
  traction for both neighbours. Diffusion donor caps also use full dt.
- E retains thermal-plus-resolved-kinetic semantics. EOS primitive SGS is k=K/rho;
  ordinary passive scalar indices remain nhydro ... nhydro+nscalars-1.
- Communication, RK registers and restart payload counts use `nvars`, including K.
  Coefficient/mode/scalar-layout changes across restart are explicitly rejected.
  The pgen initialization is skipped on restart; counters are saved to restart
  parameters rather than being reconstructed from current output.
- First and resumed timestep calls occur after ghost/primitives initialization and
  feed the ordinary global minimum. CFL is applied once by Mesh::NewTimeStep.
- Output fields distinguish `sgs_energy` density and `sgs_specific` specific energy.
  History independently sums SGS, represented and thermal inventories; restart
  preserves the cumulative source/suppression diagnostics.

## Remaining qualification gaps and numerical limits

- Compatible thermal partition is semidiscrete. The FE stress substep changes
  resolved kinetic energy quadratically; explicit work can produce a thermal
  residual of order dt^2 per substep. Trial admissibility does not validate this
  partition or prove RK2 convergence. P17/P19 need separate component inventories
  and at least three temporal refinements with spatial error suppressed.
- LLF speed squared `gamma P/rho+6k` is a frozen-shape stabilization estimate,
  not a universal eigenvalue bound for the normalized-gradient model. The exact
  zero-gradient convention remains directionally discontinuous.
- The local diffusion estimate uses cell k, whereas face kappa divided by a low
  receiving density can be substantially larger at strong density contrasts.
  Donor caps/trial checks protect nonnegativity but do not prove accurate unmodified
  diffusion there. New diffusion-suppression diagnostics must be exercised.
- Scaling full stress changes its applied trace from -2K to -2 theta K. Raw
  closure identities remain valid; the limited applied operator is a declared
  numerical modification. Quantify its use before interpreting closure accuracy.
- Thermal/density floor prerequisites are checked before conversion. At exact
  entropy-floor equality the predictor accepts >= while underlying EOS floors <=;
  this may cause an explicit abort after a passed predictor. Use matching strict
  equality semantics if floor-boundary states are intended to be supported.
- The wave initial condition uses isotropic K proportional to rho^(5/3) and
  c squared gamma P/rho+10k/9, whereas full mixed longitudinal shape gives
  b=2/3+4C2/3 and c squared gamma P/rho+b(1+b)k (4.16k at C2=0.7),
  away from gradient-zero points and excluding eddy damping. The current wave
  remains a smoke configuration, not an exact full-model propagation oracle.
- Ordinary mass-flux-upwind SGS advection is retained. Its positivity/convergence
  and interaction with stress must be tested under nonuniform rho/k and restart;
  source algebra or a uniform decay result cannot establish this lifecycle gate.

No executed PDE, fallback, decomposition, restart or output result was independently
reproduced in this read-only review packet. Numerical correctness should remain
qualified by the actual deterministic evidence ledger; smoke tests do not close M4.

## Reference preparation source inspection

Read `tst/sgs/filter_reference.py` after the lifecycle agent's handoff. Its positive
separable periodic filter, Favre density weighting, negative covariance stress,
K=-trace(tau)/2 and production=tau:A conventions agree with the scientific
baseline. The raw Python matrix formula is independent of the C++ helper.
Even-width half-cell displacement is explicitly documented and applied equally
to all fields; there is no downsampling. Width*dx is the declared top-hat
cutoff convention, not a proven mapping to the paper's Gaussian filtering scale.
No new tests executed in this packet. The direct covariance subtraction can
lose variance under very large boosts; rejecting negative K is honest, but a
stable covariance calculation and Galilean conditioning tests are needed if
such data are to be compared. The synthetic-only/no-validated-reference label
correctly prevents pipeline preparation from being presented as M7 validation.

## Focused timestep decision after prototype failures

Coordinator reported actual prototype contact advancement near 1e-15 after
1,000 steps and diffusion advancement 0.00485 after 1,000 steps with the
unconditional `8 C2 k/norm(A)` restriction. These run results are coordinator
evidence, not simulations executed independently in this review packet.

**Decision: removing this unqualified tangent estimate is defensible for continuing
deterministic verification**, while retaining the raw published closure unchanged.
Replace it with the declared finite-amplitude LLF regularization speed squared
`gamma P/rho+6k`, explicit eddy/diffusion/decay restrictions, shared-face donor
limiting of negative production and diffusion, and full-dt SSP trial admissibility.
Record the failed estimate and revised contract before rerunning; do not change
test thresholds to hide results. This does not close M1 or prove general stability.

The normalized nonlinear tensor has bounded finite amplitude proportional to K,
so momentum kicks scale as dt*K/(rho*dx)=dt*k/dx. A timestep proportional to
dx/sqrt(k) bounds these finite kicks even when an arbitrary-direction tangent
is singular. The frozen-shape pressure contribution is also bounded by the
chosen 6k stabilization term. Neither fact proves accuracy in a stiff direction-
changing perturbation: timestep refinement, active limiter diagnostics and spatial
refinement must characterize that case. In particular, full-stage admissibility
cannot substitute for temporal order or exclude bounded-but-wrong evolution.

A stencil roundoff convention `abs(du)<=8 eps (abs(u_L)+abs(u_R)) => du=0`
is defensible when constructing the closure's velocity-gradient tensor. It makes
explicit that velocity recovery/arithmetic can leave uncertainty of order machine
precision times bulk velocity, which would otherwise select an O(K) anisotropic
stress from a roundoff gradient. Apply the same convention to the paired
differences used in transverse gradients; retain the raw helper's tiny-gradient
behaviour and its existing oracle tests.

**Do not independently threshold the velocity jumps in compatible production.**
Sigma_h must still use actual unthresholded neighbour velocities and the final
applied traction, so Sigma_h=W-u_c dot F remains the exact semidiscrete identity.
Thresholding production jumps without changing work/force would break that
identity. At a noise-level gradient, isotropic closure traction can legitimately
have a correspondingly tiny compatible production residual.

The magnitude-based roundoff convention depends on bulk-velocity conditioning.
Test constant-velocity contacts plus controlled Galilean boosts and finite gradients
just above/below the declared threshold, preserving explicit limitations for
physical gradients smaller than the velocity differences reliably resolved by
floating-point arithmetic. No repeated K floor, raw norm epsilon or stress-formula
change is authorized by this decision.
