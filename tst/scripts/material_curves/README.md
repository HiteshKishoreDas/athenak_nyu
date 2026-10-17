# Material-curve validation

`analytic.py` provides independent exact position/tangent and geometry oracles.
It does not import or call any production curve sampler, integrator or geometry
routine. Its RK implementations only verify reference convergence and the error
gate; their success is not evidence for the production integrator.

Standalone checks (no build required):

```
python3 -B tst/scripts/material_curves/analytic.py --output /absolute/path/oracles.json
python3 -B tst/scripts/material_curves/harness.py > /absolute/path/harness.json
```

Both modules also implement the existing regression runner's `run()`/`analyze()`
convention. The runner builds AthenaK even for oracle-only tests and removes its
build directory. Store result artifacts outside that directory and run it only
from an isolated checkout. These modules require only Python's standard library.

## Frozen gates

These thresholds were chosen before examining production curve output:

* Roundoff checks: absolute `2e-12` on unit-scale fixture quantities.
* Independent Simpson quadrature: absolute `1e-12`, verified against a further
  factor-16 tolerance reduction for a strained circle.
* Smooth RK2 position gate: four timestep levels with factor-two refinement;
  both final successive slopes must be at least `1.8`; final maximum material
  position error must be at most `1e-4` per unit characteristic length.
* Dense RK4 reversing-shear reference: both final successive slopes at least
  `3.8` and final maximum intermediate position error at most `1e-8`.
* Independent polygon sample refinement: circle/helix relative length error
  and helix absolute curvature error at most `5e-4` at 64 intervals, with
  both final successive slopes at least `1.8`.

The RK2 oracle selfcheck uses uniform `u_x=sin(t)` to `t=0.8` with 8/16/32/64
steps. The reversing reference uses `u_x=A sin(k*y) cos(omega*t)`, with
`A=0.8`, `k=2.3`, `omega=1.7`. Its exact map is
`x=x0+A sin(k*y0) sin(omega*t)/omega`, `y=y0`, `z=z0`; intermediate errors
are evaluated at one quarter period with 4/8/16/32 steps. The full-period return
is checked separately and cannot replace the intermediate error/convergence test.

Positions use continuous lifts and common immutable initial material labels.
Normalization uses a declared characteristic length rather than distance from
the origin. Circle polygon accounting is checked against
`2*N*R*sin(pi/N)`; comparison to the smooth circumference is a separate question.
Curvature is unsigned inverse circumradius, with `None` for coincident triples.
Helix smooth length/curvature are analytic; discretized helix geometry requires
an independent sampling sweep when used to assess production output.

`geometric_error` measures symmetric sampled distance to line segments. Refine
its subdivision count independently before claiming a geometric error bound.
For nonzero winding, pass the continuous displaced closure endpoint explicitly
and treat that lifted polyline as open for distance measurements. Length uses
an explicit closure displacement. Trapezoidal arc weights are current arc-length
weights; local stretch retains its initial material-tangent denominator.

Oracle reports explicitly declare `implementation_execution=false`. The launcher
report uses mocked commands and does not claim a real MPI process executed.
Malformed/NaN/empty output, stale stage velocities, first-order convergence,
one favourable slope, and excessively large final error must fail the analyzer.
Actual insertion, migration, AMR, restart, GPU and fluid-coupling counters belong
to executable-path tests and cannot be inferred from these reference selfchecks.

## Production-path fixtures

Build the dedicated prescribed-field problem in a separate build directory:

```
cmake -S . -B /absolute/path/build -DPROBLEM=material_curve_test \
  -DAthena_ENABLE_MPI=ON -DCMAKE_BUILD_TYPE=Debug
cmake --build /absolute/path/build -j8
python3 -B tst/scripts/material_curves/production.py \
  --executable /absolute/path/build/src/athena \
  --output-dir /absolute/path/new-results --suite all
```

Use `--launcher 'launcher --options'` or `CURVES_MPI_LAUNCHER` when the local MPI
installation requires another launcher. Individual suites are `core`, `extended`,
`mpi`, `amr`, `adaptation`, `regression`, and `amr_restart`. The existing runner can run
`material_curves/production` with `--cmake=-DPROBLEM=material_curve_test`; it runs
the core suite. Results and per-case inputs, hashes, executable hash, commands,
stdout/stderr, snapshots and checkpoints remain outside the deleted build tree.
Existing case directories are rejected; `--replay` reanalyzes matching recorded
inputs and commands without executing AthenaK, writing `replay-SUITE.json`.

The pgen imposes analytic density, momentum and energy through the ordinary user
source hook after each RK2 conserved update. Boundary fills and primitive
conversion subsequently populate the actual mesh velocity stencil. No analytic
advection bypass exists in the production curve sampler. Uniform, sinusoidal,
affine and reversing nonlinear fields are explicitly prescribed kinematic flows.
The physical acceleration control disables that hook and evolves AthenaK's
ordinary constant-acceleration source. The separate `evolved_shear.py` test disables
the hook, sets `constant_pressure=1`, and evolves an exact stationary ideal-gas
sinusoidal shear on four mesh levels using HLLC. It checks the evolved fluid
profile and independently measures curve spatial and geometric errors. The passivity control compares complete
native binary fluid payloads with curves enabled/disabled, at that format's float32
output precision.

On an isotropic ideal-gas mesh, pressure uses
`cs=signal_speed-max(abs(v_i))`, `P=rho*cs^2/gamma`, so the maximum cell signal
speed is constant. Fixture choices keep `cs` positive. Four CFL-number halvings
therefore isolate temporal convergence. Affine fields stay away from their
periodic seams. The nonlinear spatial sweep uses a periodic steady sinusoidal
shear and a dense vertical set of fixed anchors; its time integration is exact,
so mesh interpolation error is isolated.

The additional strained-circle placement control freezes the original `1e-4`
material-position budget and `1e-3` geometric-distance/relative-length budgets
before candidate execution. It compares actual post-birth cubic nodes to the exact
affine circle at their immutable labels and to linear interpolation of the same
material parents. Dense exact-ellipse geometry and independent derivative
quadrature supply references. Curvature discrepancy against exact same-label
circumcircles is reported separately, without an adaptive-curvature accuracy claim.
The deliberately kinked control must preserve its initial unresolved quality after
insertion. Cap exhaustion must retain coherent IDs/anchors and persistent status.

MPI counters measure changes in the fluid rank serving a velocity sample. Whole
curves keep their deterministic owner; these tests do not claim node-storage
migration. The rank-crossing seeds start less than one predictor displacement from
a block/rank interface. The AMR control requires positive actual creation/deletion
counts and an additional case requires new curve births after mesh refinement.
Restart comparison uses real paired fluid/curve checkpoints from two ranks and
continues on four, auditing identities, ancestry, placement, births and geometry.
The AMR restart interaction uses refinement interval 5 and check cadence 3,
checkpoints cycles 30/33 immediately before/after coarsening eligibility, and
resumes on 2→4 and 1→1 ranks. Every subsequent cycle must match physical time,
logical leaves, per-leaf refinement clocks, identities, ancestry, placements,
births and coordinates. A delayed coarsening step must fail even if the final
mesh agrees. `--allow-unqualified-accelerator` appends the explicit accelerator
opt-in to every case input while retaining all numerical gates. Accelerator
execution remains experimental until independently qualified.
GPU qualification, open boundaries, other integrators, long turbulence and statistical
qualification remain outside these bounded fixtures.
