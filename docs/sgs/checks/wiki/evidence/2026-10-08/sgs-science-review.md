# P1 independent science/numerics review

Reviewed 2026-10-08 against AthenaK base
`07c529202befbeb6ef5e694d95286929618ff7c1` in
`/Users/hitesh/hitesh/git/athenak_nyu`, branch `gotham_development`.
Read-only packet: no solver changes, builds or simulation runs. This report is
an independent derivation, not evidence that milestone gates passed.

Read the requested development-plan.md, test-matrix.md, agent-runbook.md and
subject.md under `/Users/hitesh/Documents/das-wiki/science/athenak-subgrid-turbulence`.
Inspected hydro tasks, update, FOFC, flux/timestep and driver RK weights, with the
existing viscosity face stencil as the integration precedent.

Primary scientific source retrieved directly:
[Schmidt & Federrath (2011)](https://arxiv.org/html/1010.4492), especially
Eqs. 5–14, Secs. 3.2, 3.4 and 4. The arXiv HTML is the retrieved article;
an explicit archival version/hash is still needed in the final release manifest.

## 1. Independently checked continuum convention

K is SGS kinetic energy density, k=K/rho is specific SGS energy, and
E=e_int+|m|^2/(2rho) excludes K. Pressure uses thermal e_int only. The
represented energy is E+K.

Set A_ij=partial_j u_i, S=(A+A^T)/2, d=tr(A), S*=S-d I/3,
Q=sum_ij A_ij^2. For Q>0 the paper's Eq. 8 becomes

tau_ij = 2 C1 Delta sqrt(2 rho K) S*_ij
         -2 C2 K (A A^T)_ij/Q -(2/3)(1-C2)K delta_ij.

The factor two in the paper's gradient-norm definition cancels the numerator
factor two. There must not be an extra factor two in the normalized tensor.
The trace is -2K, turbulent pressure is 2K/3 independent of gas gamma, and
Sigma=tau:S=tau*:S-(2K/3)d. Defaults C1=0.02, C2=0.7, Ceps=1.5 and
Ckappa=0.65 agree with Sec. 3.2, Sec. 3.4 and Eq. 11 discussion.

R_diss = Ceps K^(3/2)/(sqrt(rho) Delta);
kappa = Ckappa Delta sqrt(rho K);
J = -kappa grad(k).

Momentum flux contribution is -tau. Resolved-energy stress-work flux is
-u dot tau. Resolved energy source is -Sigma+R_diss; SGS source is its
negative. SGS diffusion belongs to K's flux and therefore represented energy's
flux. Subtracting SGS diffusion from E would be an additional, incorrect model.
Neglect gravity production and thermal SGS pressure dilatation for this baseline,
as the source does; this is a stated omission, not a new pressure closure.

At Q=0 choose tau=-2K I/3, Sigma=0. The Q>0 normalized tensor has a
direction-dependent limit as A approaches zero; the exact-zero branch is a
numerical convention. An additive denominator epsilon violates the exact trace.
Evaluate normalized A using its largest component to avoid squared-norm overflow
and underflow; finite input checks remain necessary.

At K=0 all raw stress, dissipation and diffusion coefficients vanish. A conservative
K equation with explicit RHS stays zero for a spatially uniform zero SGS field.
Paper Sec. 3.2 footnote discusses cancellation to an SGS-speed equation that can
activate from zero. That is a different startup formulation; do not silently emulate
it with repeated floors. A supplied positive seed must have explicit initial-energy
bookkeeping. Conservative diffusion can legitimately populate a zero cell.

## 2. Independent hand and ODE oracles

For A=a I: tau=-2K I/3, Sigma=-2Ka. Compression a<0 produces K.
Along homogeneous isotropic compression, K proportional to rho^(5/3), and
k proportional to rho^(2/3), with diffusion/dissipation disabled.

For simple shear A12=a only:
tau12=C1 Delta sqrt(2rhoK) a;
tau11=-[2C2+(2/3)(1-C2)]K;
tau22=tau33=-(2/3)(1-C2)K;
Sigma=C1 Delta sqrt(2rhoK) a^2.
The nonlinear term produces no simple-shear transfer despite anisotropic stress.

For solid rotation A12=-omega, A21=omega, S=0 and Sigma=0.
The nonlinear stress has diagonal contributions (-C2 K,-C2 K,0) plus
the isotropic term; stress is not required to be isotropic.

For A=diag(a,-a,0):
Sigma=4 C1 Delta sqrt(2rhoK) a^2; the nonlinear contraction cancels.

For A=diag(-2a,a,a):
Sigma=12 C1 Delta sqrt(2rhoK) a^2 +2 C2 K a.
At rho=K=Delta=1, a=-0.1 and published coefficients,
Sigma=0.00339411255-0.14=-0.13660588745, a concrete backscatter oracle.

For stationary constant-density dissipation only:
K(t)=K0/[1+Ceps sqrt(K0)t/(2sqrt(rho)Delta)]^2,
e_int(t)=e_int0+K0-K(t).
This tests deposited heat separately from represented-energy conservation.

## 3. Compatible finite-volume stress/work stencil

On each Cartesian face, orient normal in the positive coordinate direction and
use one shared traction t_f=tau_f dot n and u_f=(u_L+u_R)/2.
Add -t_f to momentum flux and -u_f dot t_f to resolved-energy flux.
Use normal velocity differences for normal gradients; average centered transverse
gradients from adjacent cells. Face rho and K may use arithmetic means.
Ghost requirements include transverse corner data; absent 1D/2D derivatives
remain zero while the tensor is still three-dimensional.

Write F_c=D(t), W_c=D(u_f dot t_f), and choose applied production
Sigma_h,c=W_c-u_c dot F_c. Equivalently per direction,

Sigma_h,c = [t_R dot (u_R-u_c)+t_L dot (u_c-u_L)]/(2 dx).

Then resolved SGS RHS is W-Sigma_h+R_diss=u_c dot F_c+R_diss,
and SGS RHS is -D(J)+Sigma_h-R_diss. Consequently their sum contains
only stress work and diffusion divergence, while the semidiscrete internal-energy
SGS derivative is exactly R_diss. This is a second-order smooth-space stencil.
It does not require nine persistent cell tensor arrays; only final applied face
tractions/work and stage source data are needed.

All fallback/limiting must act on the shared face traction, work and compatible
production together. A cell-centered raw tau:S source computed independently
of final face traction would conserve E+K but spoil the thermal/kinetic partition.

## 4. RK2 and finite-step partition audit

Driver RK2 is SSP Heun:
U1=Un+dt L(Un),
Unext=0.5 Un+0.5[U1+dt L(U1)].
The source is inserted with beta dt (dt, then dt/2) after flux update. Stage
w0 still holds the intended RHS primitive state; u0 has already been changed.
Save K/density/velocity for the RHS or reconstruct them from unchanged w0;
never read the updated u0 as the source state.

An explicit decay RHS with these coefficients is second order and paired.
Applying an exact decay map repeatedly with beta weights is not automatically
Heun or second order. Keep explicit decay for the initial implementation with
its declared timestep restriction, or independently derive a Strang split.

The stencil above is semidiscrete, not an exact finite-step kinetic identity.
For a stress-only FE substep at fixed density with momentum increment dm,
kinetic change=u_old dot dm+|dm|^2/(2rho), while its energy work gives the
linear term. Thus thermal residual is -|dm|^2/(2rho) before heat. RK2 should
converge away this temporal error but near-cold states need actual trial checks.
The SSP convex average contributes its own kinetic variance. Tests must assess
thermal, resolved kinetic and SGS changes individually and temporal convergence.

Existing FOFC predicts only nhydro and calls EOS floor detection; it does not
cover K or the coupled SGS sources. It runs after extra flux additions and can
replace hydro fluxes without their sources. For SGS either extend the whole
trial/fallback contract or explicitly reject the legacy FOFC mode. A full SGS
trial predictor must be checked before any floor-based C2P conversion.

## 5. Timestep and positivity findings

Isotropic SGS compression gives c_iso^2=gamma P/rho+(10/9)k. This is not
a general wave bound for the full mixed closure. For a one-dimensional
longitudinal A=diag(a,0,0), excluding eddy stress, tau11=-bK where
b=2/3+4C2/3. Its K compression equation has effective exponent 1+b and
frozen-shape acoustic speed c^2=gamma P/rho+b(1+b)k. At C2=0.7 this
SGS coefficient is 4.16, not 10/9. For C2 in [0,1], b<=2 and coefficient
6 is a conservative frozen-shape stabilization estimate. It remains an estimate,
because the normalized-gradient closure is directionally discontinuous at A=0.

For frozen coefficients, the exact nonnegative diffusion FE rate bound is
dt/rho_c *sum_faces[kappa_f/dx_f^2] <=1. The density is the receiving
cell density; replacing this with a bound using only average density can miss
strong contrasts. Dissipation only requires dt Ceps sqrt(k)/Delta<=1.
For negative local production the FE donor condition is dt |Sigma^-|<=K.
Advection, diffusion, dissipation and backscatter depletions must be combined;
independent individually-safe bounds do not prove their sum admissible.

Eddy kinematic viscosity is nu=C1 Delta sqrt(2k). Include a conservative
dx^2 explicit stress limit, with compressive/shear tensor factors and the active
dimension sum, not only hydro CFL. Normalized nonlinear stress has derivative
scale C2 k/||A|| (units of diffusivity), so an 8 C2 k/||A|| estimate for
Q>0 is a sensible provisional conservative local restriction. It cannot bound
the exact-zero discontinuity or remove spatial antidiffusive modes. Trial checks
establish admissibility, not stability of the underlying spatial operator.

RK2 limiter bounds should use full dt even on stage 2: SSP admissibility follows
the hidden full-dt FE substep, not just its final beta dt=dt/2 increment.

For c_f=t_f dot (u_R-u_L)/(2dx), both adjacent cells receive the same
signed production contribution. If c_f<0, scaling full traction by
min[1,min(K_L,K_R)/(8 D dt |c_f|)] limits each donor depletion to K/(8D)
per face. It handles exact-zero local K, where dt reduction alone cannot prevent
a negative compatible production increment. A donor SGS diffusion cap
|J_f|<=K_donor dx/(8 D dt) similarly limits depletion. Over 2D faces each
mechanism consumes at most K/4; advection and decay still need the combined
trial check. Recompute work and production from scaled traction.

Scaling full raw stress changes applied trace to -2 theta K. Keep the raw
algebraic closure identity separate from the limited applied flux. Document the
limiter as a numerical model deviation and measure frequency, suppressed transfer
and affected energy. Do not claim the limited stress still satisfies -2K.

## 6. Provisional contract assessment and required exit evidence

The coordinator proposed periodic uniform Cartesian cubic grids, Newtonian
ideal-gas double-precision serial CPU, RK2 and LLF with DC/PLM; full stress
face addition, compatible production, explicit paired dissipation, donor diffusion
and backscatter caps, c_regularization^2=gamma P/rho+6k and explicit
eddy/diffusion/nonlinear/decay limits. Full stage trial includes E and K sources;
globally rebuild first-order faces/stress/source on failure, then fail clearly if
still inadmissible. No floors on SGS runs; forcing awaits its own audit.

This is acceptable as a restricted provisional implementation contract. It is
not a proof of arbitrary-flow stability and does not close full M1/M5 by itself.
First-order reconstruction alone cannot ensure K>=0 under negative compatible
production; the face donor limit is necessary. Density and e_int checks must
include final paired sources and actual trial momentum, with no scalar/EOS floor.
Same-level faces must be numerically identical across blocks; fallback selection
must not depend on a per-block flag if neighbouring blocks could disagree.

Required independent evidence before claiming numerical qualification:

- K01–K25 raw stress/production algebra including zero norm and extreme inputs.
- P03/P04 independent decay plus heat and temporal order.
- P06/P07 compression/expansion, keeping 3D trace convention in reductions.
- P08/P17 separate stress kinetic/thermal/SGS work and convergence.
- P11/P12 diffusion of k rather than K under density contrast.
- P15/P16 longitudinal pressure waves; use mixed b response where applicable.
- R01/R04/R07–R10 zero/tiny donor and fallback on different stages/shared faces.
- R18–R24 first-step/restart/global combined timestep constraints.
- Negative production must remain signed when the donor reservoir permits it.
- No long-box result may replace these deterministic gates.

Reference/turbulence campaign, MPI/GPU, restart path and release coverage require
separate actual execution. This review did not run them or mark their IDs passed.
