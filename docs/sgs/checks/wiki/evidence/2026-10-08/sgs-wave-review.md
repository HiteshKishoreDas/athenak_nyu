# Mixed longitudinal SGS wave oracle (bounded subset, prepared before testing)

For C1=Cepsilon=Ckappa=0, C2=0.7, and only Axx nonzero,
tau_xx=-b K with b=2/3+4 C2/3=1.6; tau_yy=tau_zz=-0.2 K.
The centered-gradient norm cancels from the nonlinear term for every nonzero
Axx, including either sign. All off-diagonal stresses vanish.

The one-dimensional continuum equations imply
D rho/Dt=-rho du/dx,
D P/Dt=-gamma P du/dx,
D K/Dt=-(1+b) K du/dx,
and rho Du/Dt=-d(P+b K)/dx.
Thus K is proportional to rho^(1+b)=rho^2.6, while k=K/rho is
proportional to rho^b=rho^1.6. The distinction is essential.

Linearize about rho0,P0,K0,u0=0. The right-going eigenvector is
(delta rho, delta u, delta P, delta K)
=epsilon sin(theta) (rho0,c,gamma P0,(1+b)K0),
with c^2=gamma P0/rho0+b(1+b)K0/rho0
=gamma P0/rho0+4.16 k0.
Resolved total-energy perturbation is delta P/(gamma-1) at u0=0;
SGS energy K is stored separately. Including K in resolved IEN would double count.
The production is tau:A=-b K du/dx, and the opposite resolved-energy source
cancels SGS transfer, leaving the usual gas adiabatic pressure equation.

Suggested initialization: theta=2*pi*(x-xmin)/L+0.123 radians;
rho=rho0*(1+epsilon sin(theta)); P=P0*(rho/rho0)^gamma;
K=K0*(rho/rho0)^2.6; ux=epsilon*c*sin(theta); uy=uz=0.
The exact comparison is the *linearized* right-going wave with
 theta(t)=theta(0)-2*pi*c*t/L. Nonlinear truncation scales as epsilon after
normalization, so use epsilon=1e-5, rho0=P0=1, gamma=5/3, k0=0.3.
This k0 makes the SGS contribution substantial enough to reject the isotropic
speed coefficient. Profile powers preserve the same linear eigenvector.

Three-dimensional replication is valid: all fields depend only on x, transverse
velocities vanish, and transverse diagonal stresses have no transverse divergence.
Passive scalar fractions should be constant distinct values for this oracle,
so conserved scalar densities follow rho and fractions remain unchanged. Existing
spatially varying tracer fractions have a separate advection response and must not
be compared to a stationary initial profile. Their global masses remain usable.

## Zero-gradient convention and domain of the oracle

At exactly A=0 the implementation deliberately returns isotropic tau=-2 K/3 I.
This differs from the longitudinal nonzero-gradient limit (-1.6 K in xx), so the
mixed model has no unique direction-independent linearization at the uniform
state. The sinusoid has isolated continuum extrema with zero Axx; choosing phase
0.123 avoids aligning them exactly with initial cell/face gradient sample points
on the proposed grids. It does not remove continuum extrema or prove that no
later numerical sample hits zero. The test is a longitudinal sampled branch
oracle, not a general small-perturbation oracle. Do not add a gradient epsilon or
replace the mandated exact Q=0 branch to make this test pass. If a sample hits the
zero branch and creates a large discontinuity, report that limitation explicitly.

## Frozen proposed acceptance criteria

Run periodic replicated waves with nx1=16 and 32 (nx2=nx3 may match nx1 or be
one for a strictly 1D test), fixed-grid ideal HD RK2/PLM/LLF, CFL=0.05,
C1=Cepsilon=Ckappa=0, C2=0.7, no forcing, no floors triggered, and end time
L/(4*c). Use exact actual final time in the analytic phase, no fitted phase speed.
Compute active-cell volume-weighted L1 errors for rho, ux, P, K, each divided by
its eigenvector amplitude: epsilon*rho0, epsilon*c, epsilon*gamma*P0,
epsilon*(1+b)*K0, respectively. Freeze BEFORE results:
- all fields/state/errors finite and rho,P,K positive;
- each normalized L1 error at nx1=32 < 0.1;
- each E32/E16 < 0.4 (when E16 is above 1e-10; otherwise require E32<=1e-10);
- total E_resolved+K conserved to the existing periodic-budget threshold;
- constant passive fractions agree within 1e-12;
- no FOFC/floor corrections in this smooth low-amplitude case.

This closes only pressure-wave speed/eigenvector and convergence of the
longitudinal mixed branch. It is not a full manufactured-solution campaign,
validation of C2 on turbulent reference data, or proof at Q=0. Criteria above
are proposed and unexecuted; no scientific pass is claimed.
