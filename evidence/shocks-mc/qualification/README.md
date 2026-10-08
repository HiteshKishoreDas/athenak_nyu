# Exact-source MC scientific-review measurements

`report_stylefinal.json` and `packet_stylefinal` bind production commit38c55570,
compiled fingerprint
`78e5cb1d425daaa8ab2e06deaa40db6c8249de79cd1480c7d0b069fd1587f7c3`, and Release
executable SHA-256
`ccb1c87b2e909869fefd1b7e4d07188b7642c857071b839ced52fbe8795e9fe9`.
21 fixtures were prepared before modified outcomes (`prepared.json`). Forward
RH states are independent of detector inversions; ideal mass, momentum and
energy conservation is checked for both directions. Native defaults width4 and
minMach1.1 remain fixed. No deterministic particle crossing-time estimator or
physical recovery/strength/timing tolerance was invented. Human review remains
pending; observed extra episodes remain unresolved scientific evidence.

| Case | Completed events | Unique emitted / raw-region particles | Terminal suspended | Pressure-Mach error range |
|---|---:|---:|---:|---|
| Eight independent M2 movement seeds |128 each|128/128 each|0|overall[-.0164371,.0119524]|
| M1.2 |128|128/128|0|[-.103117,-.0156752]|
| M1.5 |127|127/128|1|[-.00927597,.00992491]|
| M3 |128|128/128|0|[-.0166745,.0102860]|
| M5 |134|128/128|0|[-3.96565,.0593853]|
| Boost2 |132|127/128|1|[-.109141,-.0771184]|
| gamma5/3, stationary, reflected, nx64/256, CFL.15/.45, live3D |128 each|128/128 each|0|per-case ranges inJSON|

M5 emits two episodes for IDs12,43,47,57,74,115. Boost2 emits two for
IDs3,52,67,80,102 and one other particle remains suspended. Counts describe
actual stochastic cell-center paths and full-region observations; they are not
forced against deterministic ideal material-particle crossing times. Weak,
reentry and censored histories require scientific interpretation. The selected
packet includes full emitted/censored observation rows, sequence-aware sampled
CSV, conditional source-parent/probability/draw/direction transferCSV, nonempty
immutable event chunks, manifest and SHA-256 for every committed chunk. Empty
chunks and full128-particle traces remain in the report's raw run directories.

Four PNG/PDF figures overlay known Eulerian planar fronts with selected MC
cell-center samples; dots indicate unreduced region membership. They contain
no interpolation-based crossing estimator. Figures were visually inspected.
`packet_prefinal_floor_serialization` retains the earlier source version where
restart-floor serialization failed; it is historical evidence. Neither that
packet nor the reviewed26d4 reports are advertised as latest source.

Reproduce measurements and selected scientific figures:

```sh
python3 tst/shocks-mc/qualification.py --executable build-release/src/athena \
  --workdir /private/tmp/shock-mc-qualification-reproduce \
  --report evidence/shocks-mc/qualification/reproduce.json
python3 tst/shocks-mc/packet.py --runs /private/tmp/shock-mc-qualification-reproduce \
  --report evidence/shocks-mc/qualification/reproduce.json \
  --output evidence/shocks-mc/qualification/packet_reproduce --plot
```

The plot command requires Matplotlib; measurements need only the standard
library and strict history reader. Each process has a30s ceiling and the suite
300s. This is local Serial CPU/double, RK2, Newtonian ideal-gas Cartesian hydro,
uniform mesh, cooling/sources/transport/FOFC disabled. Rusty and Frontier,
threaded/GPU/multirank execution, AMR and additional physics remain unqualified.
Physical accuracy/stochastic recovery criteria and human scientific signoff
remain open decisions.
