# Native forcing compatibility checks

Engineering integration checks pass for the explicitly experimental subset.
**The independent forcing thermodynamics oracle fails. Scientific qualification
and use of strengths as calibrated physical measurements remain pending.**

The new shock_tracking/experimental_forcing flag defaults false. Explicit true
permits native continuous(turb_flag2), zero-start, untiled/unlocalized, periodic3D
accel_rms driving with Serial CPU/double/RK2 ideal hydro and existing MC tracers.
Edot, tiling, localized forcing, outflow driving and other physics are rejected.
No solver, forcing energy update, density flux, task graph, event schema or particle
state layout was changed. Only qualification guards/static restart binding and
fingerprint coverage changed in production code.

## What passed

15processes:8normal successes and7expected rejects. Smooth single-mode transverse
forcing emitted zero complete/partial/suspended shock episodes.896forced samples
(128tracers×7generations) match full-double native thermo time/position/density
exactly, and completed native primitive density/pressure; pressure error is0at the
predeclared1e-12 check. Binary time metadata is explicitly only6significant digits,
so its formatted value is checked separately instead of inventing exact precision.
Tracking on/off saved double fluid fields, rendered force and particle records
match exactly. All4096particles survive. Native restart from generation30 preserves
all49152expected resumed rows, canonical committed events, all556terminal censored
records, and final fluid/OU force arrays exactly. Altered seed/acceleration reject;
failed paths preserve committed history. Native OU/modal RNG state is reused.
Static parameter value text is deliberately strict: even1versus1.0 or an unused
parameter change rejects restart. New source fingerprints reject old-source sidecars.

Fresh Release build/state/offline checks,53ShockIt parity fixtures (exact masks,
rtol=atol=1e-12),20unforced integration processes and21held-out physical fixtures
passed. All21measurement dictionaries exactly match the preceding packet. Three
independent read-only reviews cover timing, persistence and evidence interpretation.
No new Debug/MPI/GPU/cluster result or human scientific approval is claimed.

## What failed physically

A homogeneous gas with divergence-free transverse acceleration has no compressive
work and should retain its internal energy. In the first active native RK2 step,
at rms acceleration1 and dt=.01875, mean gas energy instead grows by
0.00013183593749990785. Momentum is exact. Pointwise excess matches
3rho|a|²dt²/8 with maximum error2.36e-16. This comes from native quadratic kick
energy inside the two RK stages; the source implementation is unchanged here.
The gas-energy physical oracle explicitly remains FAILED. The first full timestep
has zero OU force, so the oracle uses the second full timestep/first active kick.
This is a concrete source bias, not a shock-detector or clock mismatch. Zero
false detections does not establish thermodynamic neutrality. The next physical
gate is to resolve/qualify native RK2 forcing work before trusting shock strengths.

## Driven demonstration

32^3periodic cells,4096MC particles, gamma1.4, initial rms Mach3, forcing rms10,
OU tcorr.1/update.01/modes1..3/seed2718, PLM/HLLC, code tlim.2; no cooling.
CFL.30:7063completed/0partial/556suspended episodes; CFL.15:7744/0/520.
Completed count changes+9.64percent. CFL changes fluid integration, random draws
and detector sampling jointly; this is not isolated cadence or convergence evidence.
Final rms local Mach1.706/1.707; total energies7.01337/7.01053. Driven total energy
contains physical forcing work and its native numerical bias; conservative shock
heating and force bias must not be equated. This is a short realization, not steady
supersonic turbulence or a truth shock catalogue. Fixed-endpoint strengths remain
uncertified and open episodes remain right-censored.

## Reproduce

From /Users/hitesh/.codex/worktrees/shock-multiphase/athenak_nyu:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DAthena_SHOCK_TESTS=ON
cmake --build build-release -j8
MPLCONFIGDIR=/tmp/shock-turb-mpl \
/Users/hitesh/hitesh/git/shockit/.venv/bin/python tst/shocks-mc/forcing_checks.py \
  --executable build-release/src/athena \
  --workdir /tmp/shock-forcing-new \
  --output /tmp/shock-forcing-plots-new
```

Use new work/output directories. Each process is capped at180s. report.json records
exact commands/inputs/logs/source/executable hashes, successes, expected rejects
and the explicitly failing physical oracle. Plot driven_checks.png shows fluid
Mach/total-energy evolution and uncalibrated first-inside strength distributions.
Native timing/shear fixtures and checkpoint/event histories are in the two ZIPs;
preflight/ preserves failed fixture/parser assumptions and their corrections.
Full raw final runs remain/private/tmp/shock-forcing-checks-final-20261013.
