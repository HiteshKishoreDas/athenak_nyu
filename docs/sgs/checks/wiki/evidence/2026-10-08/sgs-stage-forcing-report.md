# Stage-consistent HD forcing packet

Owned edits only `src/srcterms/turb_driver.hpp` and `turb_driver.cpp`.
Legacy `AddForcing` implementation unchanged. Parent owns task routing,
restart payload selection/compatibility and shared SGS trial/timestep integration.

## Public API and saved parameters

- `bool stage_consistent`: default true with enabled Hydro SGS, false otherwise.
- `Real integrated_work`, `Real integrated_momentum[3]`: global cumulative physical
  forcing exchanges, restored from input/restart parameters.
- `TaskStatus StageForcing(Driver*,int)`: old-primitive RHS momentum rho*a,
  resolved energy rho*u dot a, beta*dt updates. No per-stage quadratic work.
- `static Real StageNormalization(Real m0,Real m1,Real target)` supplies the stable
  nonnegative scale for m0*s²+m1*s=target. target=0 deliberately gives s=0.
- Pin keys: `turb_driving/stage_consistent`, `integrated_work`,
  `integrated_momentum1`, `integrated_momentum2`, `integrated_momentum3`,
  `ou_state_version=1`, `ou_updates`.

Explicit stage_consistent=false with SGS enabled fails. New mode requires
Newtonian ideal-gas HD, RK2, no MHD/radiation/ion-neutral and uniform refinement
level. Validates finite forcing settings, positive correlation/update times,
nonnegative dedt, admissible solenoidal fraction/bands/spectral options and
finite restored counters. Legacy mode retains its original validation/behavior.

## RK and forcing normalization

StageForcing computes RHS from old-stage primitive density/velocity and fixed
acceleration field. It updates u0 with beta(stage)*dt; counters use half the full
dt for **each** RHS, the exact Heun quadrature. Global MPI sums use MPI_ATHENA_REAL.
The active time/start/duration predicate matches legacy AddForcing.

UpdateForcing applies new-mode cadence cap dt<=dt_turb_update before calculating
normalization moments. Positive root evaluation uses scaled hypot and a
cancellation-free positive-m1 branch; fixes original negative-m1 sign bug only
in new mode. Removes original artificial m0 epsilon only in new mode. Positive
target with zero norm fails explicitly; no hidden injected forcing.

New mode enables UpdateForcing at ncycle0; legacy's ncycle>=1 condition is retained
when stage_consistent=false. Therefore active source behavior differs intentionally
only under the opt-in/default-SGS new mode.

After normalization, SGS-enabled new mode caps dt by the global minimum
`0.5*sqrt(2*e_internal/rho)/|acceleration|` to protect the initial Euler predictor
from an oversized newly generated kick. The parent handles the trial-state forcing RHS. There is no previous-acceleration
SGSNewDt cap: this bound uses the newly generated acceleration every cycle. A thermal cap can reduce dt after
normalization, so dedt is a **nominal target**, not a promise of exact applied rate.
Cumulative actual work is authoritative; do not equate target*elapsed time to it.

## OU restart review

Ctor calls Initialize(), which zeros OU arrays and seeds RNG before pgen reader
restores raw forcing state/RNG. InitializeModes does not reset force_tmp1; it
applies a legitimate recurrence only when exact saved ou_updates is behind the
current required OU tick. Ctor restores ou_updates after old time-derived count,
so restart count is not overwritten. EnsureBasisSize can rebuild raw arrays if
mesh topology changes; new mode rejects multilevel, so no AMR/regridding claim.

Parent must serialize raw force_tmp1 in the existing 3-component forcing payload
for new mode, preserve RNG and exact update count, save counters, reject mode
mismatch/forcing parameter overrides, and skip legacy IncludeAddForcingTask when
using StageForcing through HydroSrcTerms.

## Executed checks and limits

1. Current turb_driver.cpp passed C++17 `-fsyntax-only` with baseline generated
   configuration/Kokkos include paths. No entire solver link/run by this owner.
2. Standalone executable compiled the actual StageNormalization function body and
   compared seven cases against independently evaluated long-double quadratic
   roots at frozen relative tolerance1e-12. Included positive/negative/zero m1,
   zero target with negative m1, tiny target1e-30, and m0/target1e±300. All passed.
3. Independent scalar constant-acceleration Heun oracle checked final momentum,
   retained thermal energy and exact work quadrature at tolerance1e-12. Passed.
   This is a mathematical RK oracle, not an invocation of StageForcing's Kokkos
   kernel; parent must execute actual forced solver energy/restart cases.
4. Diff whitespace check passed. Added source/header lines meet90 columns.

Temporary normalization executable `/private/tmp/sgs-normalization-oracle`.
MPI reduction, actual kernel forcing and OU restart continuity are not claimed
passed here. Parent's independent integration verifier supplies those checks.

Restart precision correction: stage-consistent construction now serializes every
finite runtime Real forcing configuration and ledger value with max_digits10.
GetOrAddReal formerly saved default kpeak=4*pi as 12.5664 instead of the runtime
12.566370614359172; at the parabolic endpoint this changes raw OU amplitudes
although the RNG state and exact OU update counts agree. Legacy serialization
remains unchanged. A new continuation comparison is required to confirm this fix.
