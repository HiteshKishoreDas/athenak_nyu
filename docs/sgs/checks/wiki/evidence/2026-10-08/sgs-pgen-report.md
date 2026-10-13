# SGS verification pgen packet

Added only `src/pgen/sgs_verify.cpp` and `inputs/hydro/sgs_verify.athinput`.
No CMake edits/builds/runs by this owner. Header packet separately owns
`src/hydro/sgs_closure.hpp`; lifecycle report is `/private/tmp/sgs-lifecycle-audit.md`.

## Interface

Pgen uses Hydro `nhydro`, true `nscalars`, appended `isgs`, and
`sgs_pgen_initialization` bool configured by parent from `<sgs>/initialization`.
When true, pgen writes K; parent must skip generic seed assignment but validate
the initialized state. When false, pgen initializes only thermal/hydro/passive state.
Restart returns without reseeding or modifying loaded values; no hooks enrolled.

Parameters: `problem/mode`, `density0`, `pressure0`, `k0`, `amplitude`, `velocity`.
Require finite rho0,p0>0, k0>=0, 0<=amplitude<1 and finite velocity. Unknown modes
fail clearly. `shear` requires at least 2D. Diffusion pressure balance must keep
thermal pressure positive. Domain phases use global physical bounds and local
cell centers, making the same physical profile independent of block partition.

## Modes

| mode | initialization | intended verification |
|---|---|---|
| uniform | rho/P/K/velocity constant | exact homogeneous decay, heat pairing, moving uniform invariance, exact zero K when k0=0 |
| shear | vx=v0 sin(2pi y/L), rho/P/K constant | nonzero shear production, stress symmetry/energy accounting |
| diffusion | K=rho0*k0(1+a cos(2pi x/L)), v=0, P adjusted so P+2K/3 constant | nonlinear specific-energy diffusion and periodic conservation in coupled equations |
| wave | rho/rho0=1+a sin(2pi x/L), P/P0=(rho/rho0)^gamma, K/K0=(rho/rho0)^(5/3), acoustic vx | SGS pressure/compression response; small a for linearized comparison |
| contact | rho varies as cosine, P/K density/velocity constant | density contact with constant effective pressure and nonuniform specific K |
| collision | vx=−v0 sin(2pi x/L), constant rho/P/K | compression/rarefaction, shock robustness at increased v0 |
| box | periodic cyclic velocity sines, deterministic density/K modulation | multidimensional bounded pilot, block seams and energy budgets |

All modes initialize true passive n at `nhydro+n` with positive specific profile
`0.2+0.05*n+0.1*sin((n+1)*phase_x)`, distinct from K. IEN is thermal plus resolved
kinetic only; initial represented energy adds initialized K. Passive concentration
is not constrained <=1 by this pgen; intended small count is two or three tracers.

Wave speed initialization uses `c_eff²=gamma*p0/rho0+(10/9)*k0`, combining thermal
and SGS compression pressure. The finite-amplitude field is an initial condition,
not an exact traveling-wave solution. Box sine flow is deterministic; no stochastic
forcing, fitted coefficients, new closure physics or promised statistical benefit.

Input is 16^3, eight 8^3 MeshBlocks, two ghosts, ideal gamma5/3,
RK2/PLM/LLF cfl=.1, tlim=.05, two true passives, full published coefficients,
history every .01, conservative/primitive binary output every .025, restart .025.
Resource cap for this verification packet is <=32^3 when parent overrides grid.

## Checks and reproduction

Static style checks: no >90-column C++ lines, tabs or trailing whitespace. Build
and execution remain parent's integration packet. Build with `-DPROBLEM=sgs_verify`;
run supplied input and override `problem/mode` to each of seven modes. Useful
variants: `problem/k0=0`, `problem/amplitude=1e-4`, `problem/velocity=3`, one block
versus eight blocks, and a restart from the intermediate .025 checkpoint.

No claim of tests passed until executable evidence is returned. Parent should check
the initial K profiles in binary output independently of initialization code and
verify all true passive indices remain independent across output/restart.
