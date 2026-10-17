# Pinned MC baseline evidence

Base `526a187f938e38634e89a88192ebe9731c4d36bb`, branch
`codex/shock-multiphase`, Kokkos
`08ceff92bcf3a828844480bc1e6137eb74028517`. Baseline executable SHA-256
`8e23c3a37b0a784191f1b99e1e02209a1208207eac0612be2336ef4d4759fd18`.
Source-file hashes are computed from the pinned Git commit because production
integration proceeded after the completed baseline build. Cache/input/runner/
criteria hashes, commands and exact outcomes are in `report.json`; compact
inputs/logs/run manifests are under `runs` and raw output/checkpoints remain in
the report's absolute run directories.

Frozen pre-runtime criteria: `predeclared.json`. Eight fixed movement seeds,
4096 tracers each, uniform positive translation across a periodic seam. The
independent physical probability is `p_n=u*dt_n/dx`; the oracle uses recorded
step times and forward Bernoulli moments, without copying the hash generator.
All per-step cell hops are0 or1. Mean z scores [-1.064,1.068], variance z
scores [-1.229,1.166] pass the predeclared six-standard-error diagnostic. This
qualifies the selected ensemble test; it does not provide a deterministic
physical particle path or certify all stochastic regimes.

Four-block run:1024/1024 exact tags retained and1024 changedGID; mean/variance
z scores-.891/-2.403. Same-partition restart at generation17:102400 continued
particle observations exactly match uninterrupted output. Zero particles run
successfully and emit no thermo observations. There is no multirank/GPU claim.

The first process succeeded, but strict thermo parsing found identical terminal
rows duplicated by the pre-existing `dcycle=1` output's final append. The
original failure is preserved in `report_initial_strict_thermo_duplicate.json`.
Legacy thermo analysis permits only exactly identical replay rows and records
4096 replay rows per seed; conflicting duplicates still reject. The new shock
event reader's canonical identities remain strict. No physical bounds changed.

Reproduce from this checkout, using a fresh directory:

```sh
python3 tst/shocks-mc/baseline.py --executable build-baseline/src/athena \
  --workdir /private/tmp/shock-mc-baseline-reproduce \
  --report evidence/shocks-mc/baseline/reproduce.json
```

Each process has a30s ceiling, the suite300s. Integration against the modified
binary and human scientific approval are separate gates. Rusty node, Frontier
backend, intended scale and performance budgets remain unspecified.
