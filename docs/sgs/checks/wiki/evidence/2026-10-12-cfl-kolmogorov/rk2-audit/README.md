# Direct RK2 factor audit, 2026-10-12

Read-only native solver audit plus two small uniform-decay runs. No source edits.
The current SGS state update uses beta*dt (1,.5), and the second-stage state blend
halves the first-stage increment. Final source/history weights are (.5,.5).
Sources/rates are computed from the same pre-update primitives as stage fluxes,
then added inside each RK stage. This SGS path is not a separate Strang half-step.

Native 16³, rho=3, initial specific k=.04, Cepsilon=1.5, zero velocity, all
published SGS defaults, two OpenMP threads, existing sgs_verify Release binary.
A one-step dt=.0001 run and a .025-time run (11 steps) were checked against the
exact discrete Heun recurrence using the saved native history times:

```
f(k)=-(Cepsilon/dx)*k^(3/2)
k*=k+dt*f(k)
k_next=.5*k+.5*k*+.5*dt*f(k*)
```

Maximum absolute specific-k recurrence errors:8.535e-16 and1.894e-15,
below the prespecified5e-13. Deliberately multiplying f by .5 gives endpoint
errors9.595e-6 and.002104, readily distinguishable. The independently read final
restart inventories agree with the full recurrence. Both native verification
records pass. Previous native manufactured production tests separately establish
C1/C2 production amplitudes at tiny timestep; conservation alone is not the basis
of the conclusion. No global factor .5 is detected in the supported RK2 source path.
This does not prove small finite-timestep bias in turbulent runs or paper agreement.

The first attempt disabled history by setting dt=0 without dcycle. Its native
one-step run completed, but the runner retained FAIL because the expected history
was absent. The next attempt used dcycle=1 and passed native verification, but
analysis encountered an identical repeated final history row. Both attempts remain
in sibling folders sgs-rk2-factor-audit and sgs-rk2-factor-audit-cycle-output.
The checked attempt deduplicates only same-time identical SGS-energy rows; strictly
positive time increments drive the recurrence. Tolerances and native code unchanged.

See audit-results.json, results.json, protocol.json, manifest.json and audit.py.
The script uses the existing Verification runner and retained native binary at
/private/tmp/athenak-sgs-heating/src/athena. Its output directory must be fresh.
The test executable/source hashes are recorded in manifest.json.
