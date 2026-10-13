# Independent serial-double SGS verification: P4–P6 bounded packet

**Revised frozen numerical source: 43 executed result records passed, zero failed,
zero blocked.** An additional normal checkpoint restart after actual fallback
passed both active-cell state and integrated diagnostics exactly. This verifies the
listed deterministic instances; it does not close the full wiki catalogue or establish
statistical turbulence validation, MPI/GPU support, forcing, or every hard-flow regime.

## Files and reproduction

- Production source belongs to the coordinator; this packet edited only
  `tst/sgs/verify.py` and `tst/sgs/thresholds.json`.
- Build `/private/tmp/athenak-sgs-verify-debug`, `PROBLEM=sgs_verify`, Debug,
  Kokkos Serial and bounds checking, MPI/OpenMP disabled, double precision,
  Apple Clang 17 / CMake 4.1.2. Build concurrency `-j3`.
- Main passing evidence: `/private/tmp/sgs-verification/qualification-02`.
  `manifest.json` includes executable SHA-256, base commit, tracked diff digest
  and every source-file hash. `source.patch`, `new-source/` and `CMakeCache.txt`
  preserve the source delta/new files/build configuration used for reproduction.
- Every input has its hash; every invocation has command, return code, runtime,
  actual final time/cycle, minima and measured-versus-threshold checks in `results.json`.
- Additional fallback restart: `/private/tmp/sgs-verification/fallback-restart-01`.
- Additional active stress limiter: `/private/tmp/sgs-verification/limiter-02`.
- Initial failures remain under `build.log`, `rebuild.log`, `pilot-02` and
  `qualification-01`; passing reruns never replace those files.

```sh
cmake -S /Users/hitesh/hitesh/git/athenak_nyu \
  -B /private/tmp/athenak-sgs-verify-debug -DPROBLEM=sgs_verify \
  -DCMAKE_BUILD_TYPE=Debug -DAthena_ENABLE_MPI=OFF \
  -DAthena_ENABLE_OPENMP=OFF -DKokkos_ENABLE_SERIAL=ON \
  -DKokkos_ENABLE_DEBUG_BOUNDS_CHECK=ON
cmake --build /private/tmp/athenak-sgs-verify-debug -j3
python3 -B /Users/hitesh/hitesh/git/athenak_nyu/tst/sgs/verify.py \
  --exe /private/tmp/athenak-sgs-verify-debug/src/athena \
  --output /private/tmp/sgs-verification/new-unique-run
```

Use the recorded source delta and unique output directories. The current runner also
includes the later active-limiter and fallback-restart cases; these were executed in
separate evidence directories after the 43-record matrix. `--pilot` runs seven modes;
`--stress-only` executes the conditional component time-refinement packet.

## Independent measurement and frozen limits

AthenaK binary field output stores float32 even in a double build; formatted table
output permits only 1D slices. Neither supports independent full-3D 1e-10 inventories.
The runner instead reads native-double serial restart files with Python's standard
library. It validates metadata, byte count, variable count, refinement level, active-cell
count and unique global cell coordinates; excludes ghosts; and sums with `math.fsum`.
Initial uniform-state inventory independently recovers density 1, momentum 0.3,
resolved energy 1.545, distinct passive scalar masses 0.2/0.25 and SGS energy 0.01.
High-precision history is inserted into verification-only input copies. Existing history
labels truncate at 10 characters; the reader explicitly resolves those known labels.

Frozen scaled absolute tolerances use `max(1,abs(a),abs(b))`:

| Measurement | Tolerance |
|---|---:|
| Periodic represented E+K | 1e-10 |
| Paired applied production/dissipation versus K change | 1e-8 |
| Mass, momentum and ordinary scalar mass | 1e-11 |
| Block decomposition cellwise scaled Linf | 1e-12 |
| Normal intermediate restart fields/diagnostics | 1e-12 |
| Online history versus independent native-double inventories | 1e-10 |
| Observed temporal order interval | 1.7–2.3 |

Density and thermal energy must be strictly positive; SGS energy nonnegative; all
stored state finite at every available checkpoint. Reaching nlim before requested
tlim fails. Per-run timeout is 20 seconds, maximum nlim 1000, maximum grid axis 32.
The revised 43-record matrix used 16³ or smaller and approximately 13.8 seconds
aggregate measured simulation runtime. It did not start a long statistical campaign.

## Passing revised matrix

- Seven modes: uniform, shear, nonlinear diffusion, pressure wave, material contact,
  collision and deterministic unforced box. Applicable instances cover 1D/2D/3D;
  shear begins in 2D. Three-dimensional short cases end at t=0.004, other default
  cases at t=0.01. Exact command/input files distinguish the concrete instances.
- All completed states pass periodic represented-energy, mass/momentum/scalar,
  admissibility, applied-source ledger and independent history-inventory checks.
- Worst scaled E+K drift: **2.94099e-16**. Worst paired-transfer residual:
  **6.93889e-18**. Worst mass drift: **4.44089e-16**.
- Worst thermal-history inventory residual: **4.39614e-14**; kinetic-history
  residual **2.78944e-15**.
- 16³ box, blocks 8³ versus single 16³ block: state scaled Linf **0**.
- Normal intermediate box checkpoint restart: state and production/dissipation,
  stress-limiter and fallback counters **exactly reproduced**.
- Diffusion/contact at CFL 0.1/0.05/0.025: completion, conservation and positivity pass.
- Tiny K=1e-12 contacts with boosts 0, 0.3 and 100 pass. At boost 100, minimum thermal
  energy was 1.499999999985448: frame conditioning is visible, not treated as an exact
  Galilean-invariance proof. It used 163 cycles.
- SGS k0=0.1, thermal pressure 0.001 wave passes the declared stability/budget smoke.
  Its nonlinear-wave initialization is not an exact wave-speed/partition oracle.
- Hard collision v=5, pressure=0.01 activates **9 fallback stages** and passes all
  frozen checks. Hard collision v=20, pressure=1e-4, k0=1e-10, CFL=0.2 activates
  **19 fallback stages**, with all reservoirs admissible and budget checks passing.
- Three-CFL uniform dissipation: K errors 9.88706e-8, 2.55799e-8, 6.38135e-9;
  temporal orders **1.95053, 2.00308** against the independent law
  `K(t)=rho*k0/[1+ceps*sqrt(k0)*t/(2*dx)]^2`. Thermal oracle is
  `pressure0/(gamma-1)+rho*k0-K(t)`; both reservoir errors are recorded.

## Conditional kinetic/thermal component convergence

Smooth 16² transverse shear, t=0.04, C2=ceps=ckappa=0, CFL 0.1/0.05/0.025/0.0125
passes the full-solver Richardson check. Revised-source field orders are
**2.01919, 2.01173**; kinetic, thermal and SGS reservoir orders are approximately
2.02 and 2.01. The finest CFL is a numerical reference, not an exact PDE solution.
LLF numerical heating remains present. This establishes time convergence of the
resolved kinetic/thermal/SGS components under that conditional coefficient isolation;
it does not independently validate exact stress-work partition or the full mixed
closure's manufactured PDE. Those remaining catalogue gates stay open.

## Additional active-event checks

**I22 normal restart after actual fallback:** the v=5 collision checkpoint was taken
with 3 prior fallback stages, before final-time clipping. Evolution to t=0.01 reproduced
all active fields and the production/dissipation/stress-limiter/fallback counters
exactly (scaled Linf 0). The final total was 9 fallback stages. This is separate from
restarting a run that terminates with a clipped final dt, which P0 identified as a
pre-existing timestep-schedule difference.

**R01 active stress limiter:** a deliberately coarse four-cell collision with v=20,
k0=10, pressure0=0.1, CFL=0.2 activated integrated rejected stress transfer
**0.2126880755333241** and passed admissibility, conservation, source-ledger and
history checks. This is an adversarial arithmetic/limiter case, not a turbulence
accuracy result. Minimum SGS energy 4.539993220673556; no fallback was needed.
Earlier coarse low-K and smooth high-contrast diffusion/contact attempts are retained
as passing characterization cases but did not activate their requested limiter.
**Active diffusion-limiter behaviour remains unexercised.**

## Initial failures and their disposition

1. Two initial build attempts failed at `mesh.hpp:209`, member access into incomplete
   `MeshBlock`, through new pgen include order. The source owner fixed include order;
   rebuild-02 and rebuild-03 passed. This verifier edited no production file.
2. Initial runner history parsing expected untruncated labels; successful physical runs
   were logged as parser failures. Reader correction resolved the documented 10-character
   labels. Those setup failures remain in pilot-02.
3. Before the approved timestep revision, contact 1D stalled at t≈1.29e-15 after
   1000 cycles; diffusion 1D stopped at t≈0.004853742 of requested 0.01. Other diffusion
   and contact dimensions hit the cycle/timeout cap. Qualification-01 records
   **20 pass / 6 fail**. No completion failure was reclassified as a pass.
4. An initial 16³ box with t=0.01 exceeded the 20-second pilot cap. A separately declared
   shorter t=0.004 box passed. No stationary turbulence interpretation followed.
5. Coordinator/reviewer approved removal of the failed unconditional k/|A| timestep
   estimate, roundoff-only closure-gradient convention, receiving-density diffusion
   row bound and zero-K timestep corner fix. Revised qualification-02 reran the affected
   matrix with identical tolerances. This resolves those observed cases, not a proof of
   unconditional stability of the nonlinear normalized stress.

## Scope still open

No exact full nonlinear manufactured PDE, exact compression/expansion oracle, long-time
statistics, forcing-work validation, GPU, MPI, AMR, single precision or unsupported
solver/mode qualification is claimed here. Pressure-wave cases are smoke/budget cases;
conditional stress temporal convergence is not an exact partition oracle. Active
fallback and stress limiter are observed; active diffusion limiter is not. The full
catalogue remains a ledger of pass/partial/not-run cases, not 43 fully closed families.

## Final bounded diffusion-limiter characterization

Two additional approved term-isolation attempts are retained in
`/private/tmp/sgs-verification/diffusion-limiter-03`: c1=c2=ceps=0, ckappa=100,
amplitude=0.999999, k0=0.01, pressure0=1, velocity=0, CFL=0.2, tlim=0.02.
The requested four-cell case completed in 9 cycles. A six-cell case samples the
sinusoidal minimum at initial K=1e-8 and completed in 15 cycles. Both passed strict
admissibility, periodic represented-energy budget, momentum/scalar inventories and
SGS source ledger checks; the four-cell represented-energy residual was
2.94099e-16. Both reported zero integrated rejected diffusion. These are passing
characterization cases, and **active diffusion-limiter behaviour stays unexercised**.
No threshold was loosened and no further activation tuning was performed.

Stage forcing was subsequently verified in a separate bounded packet. Its exact
physical results, retained restart defect and correction, actual-work budgets and
unsaturated pilot limitations are documented in
`/private/tmp/sgs-forcing-verification-report.md`.
