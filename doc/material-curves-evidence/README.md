# Material curves on multiphase_turb: port evidence

Base: `526a187f938e38634e89a88192ebe9731c4d36bb` (`multiphase_turb`).
Branch: `codex/multiphase-material-curves`.
Original implementation source: `d854fe64c0aafb05db5938c80619b20105b10609`
on `codex/material-curves`, based on `07c52920`. Old build/test artifacts have
not been relabelled or copied as qualification of this port.

## M0 and integration contracts

The local branch and cached origin ref both identify the base above. Initial
unprivileged remote verification failed DNS; no remote branch update was assumed.
Pinned target Kokkos is `08ceff92bcf3a828844480bc1e6137eb74028517` (4.4.00).
Separate baseline Release/MPI and port Debug/MPI builds use that exact pin.
Hardware is the local Apple M4 CPU; Kokkos Serial, AppleClang 17, Open MPI 5.0.9.
Compiler parallelism is bounded at four jobs per build; local tests use 1/2/4 ranks.

The support scope remains Newtonian Cartesian HD, ideal/isothermal, double,
2D/3D periodic boundaries, RK2, insertion-only, and shared native checkpoints.
The RK2 coefficients are unchanged. The predictor samples valid initialized
primitives before stage1; the corrector samples refreshed stage1 primitives and
ghosts before stage2. Curves write no fluid arrays or timestep constraints.
Curve-enabled AMR is synchronized before adaptation/checkpoint output; native
feature-off output-before-AMR ordering is preserved.

Target changes require null-safe source validation, actual standalone cooling
status independent of history output, and separate shearing/orbital/frame guards.
Forcing, MC/ordinary particles, cooling, frame tracking, RK3/other integrators,
rank/node-sharded checkpoints, MHD/relativity and physical exits remain rejected.
The input interface for actual constant acceleration is now `hydro_srcterms`.

Schema2 and its integer/identity/birth/ancestry/quality/AMR-clock semantics are
unchanged. The paired native file checksum binds dt/time/cycle as before.
Target `dtold`/`dt_last_completed` are outside the curve runtime and are not
new curve state; MC/cooling-history configurations using them are unsupported.
Existing target restart payloads for MC/forcing/frame remain unchanged feature-off.

## Fresh target qualification

Runtime/examples: `313977256c95e0aee17433afed27379f0b5579bc`.
Tests/reader/native AMR selector: `15b40a7ceb35521651244326fa7cab92b87552aa`.
The 323 runtime source files are byte-identical between those revisions and the
archived build source. Documentation/evidence commits through `dda93e67` do not change runtime.
The reviewer cleanup below modifies two runtime files and has fresh validation.
The generated Kokkos git-state file identifies Kokkos, not the application revision.

| Gate | Executed result |
| --- | --- |
| Builds R07 | Untouched Release/MPI baseline; port Debug/MPI fixture with bounds/DualView checks; port Release/MPI builtin; port Release without MPI fixture: pass |
| Independent reference/harness M1/R04/R05 | Analytic/order/geometry and wrong-result/failure detection checks: pass |
| Final integrated suite | 93 passing records from 78 fresh executable cases; `implementation_execution=true`, `analysis_replay=false` |
| Serial execution | 31 core records on the non-MPI Release fixture: pass |
| E02 temporal order | Sinusoidal, rotation and strain final slopes approximately 2; frozen minimum 1.8 and final position budget 1e-4 retained |
| E03/evolved Euler shear | Spatial slopes 2.00044/2.01492/1.98214; finest position error 2.40381e-5; fluid shear profile error 0 |
| D04 insertion/geometry | Cubic error 3.76993e-5 versus linear 1.90680e-3; geometric error 5.02422e-4 within 1e-3; 96 insertions with 64 post-birth cubic samples |
| MPI/periodic/AMR | One/two/four ranks, empty curve ranks, actual owner migration; refinement/coarsening creates/deletes 12 blocks; AMR inserts 21 samples |
| Q10 restart | Checkpoint cycles 30/33 with check interval3/refinement interval5; 2 to4 and1 to1 ranks; every-cycle leaves/clocks/births exact and coordinate error0 |
| O02 ordinary rank-change restart | 2 to4 ranks; position error 3.33067e-16 within 2e-11 |
| Physical fluid controls | Six actual Euler Sod/contact/uniform cases: pass; finest Sod trajectory error 2.02668e-4 within 2e-3; contact 6.20634e-17 and uniform 6.66712e-16 |
| R01 passivity | Two curves-off/on conserved payloads byte-identical at native float32 output precision |
| Feature-off native branch parity | Six probes: plain HD, active standalone cooling, actual frame boosts, forcing, MC serial and two MPI ranks; all 45 native outputs byte-identical to the untouched target |
| Native hydro regression | Six RK2/PLM/HLLC sound/entropy launches; frozen errors/convergence/symmetry gates pass and every error file is byte-identical to target baseline |
| R03 exclusions | Nine native-valid curves-off initialization controls paired with nine exact curve rejections; all pass within final integrated run |
| I/O | 10 independent reader tests; native C++ exact metadata/clock restore and corruption/mixed-pair rejection; MPI2 to4/1/2 including empty curve ranks: pass |
| Documented examples | Line/loop/AMR examples reach requested final times with quality0; loop inserts112 samples; two-rank AMR inserts21 samples |

[Final qualification summary](qualification-summary.json),
[production records](production-final/results.json),
[feature-off comparisons](branch-feature-off/results.json),
[manifest](manifest.json), and the adjacent cohort inputs/logs/configuration record
exact commands, executable/source/input hashes, ranks, measurements and tolerances.
The final Debug executable SHA256 is
`1386ca430052fa412628f5f823412973606678bb8232bc09c5dbfbf184c1ded6`;
builtin Release is
`14461bc1c008dcfa3738efe5c171ff4e703fb1ff47929a52e4ae681cbaec7f39`.

## Reviewer cleanup qualification

The six accepted simplicity findings remove 21 net runtime lines: reuse the
validated seed count and home-rank decision, use one ID map for duplicate and
parent checks, remove the repeated end-step audit and unused include, and share
curve/node decoding between file reads and MPI restore. The defensive native
restart write guard remains. Solver hooks, numerical gates and schemas are unchanged.

Fresh Debug/MPI and Release/non-MPI builds pass. The full integrated suite passes
93 records from 78 actual cases; the serial core passes 31 records. All 390 curve
output/checkpoint files match the pre-cleanup cohort byte for byte. Ten independent
reader tests and both serial/MPI C++ persistence checks pass, including exact
metadata/AMR-clock restore, empty ranks, rank changes and corrupt/mixed-pair
rejection. Both changed C++ files pass cpplint and length/tab/whitespace checks.

[Cleanup summary](reviewer-cleanup/summary.json),
[production results](reviewer-cleanup/production-results.json),
[serial results](reviewer-cleanup/serial-core-results.json),
[comparison](reviewer-cleanup/before-after-comparison.json),
[commands](reviewer-cleanup/commands.json), and
[manifest](reviewer-cleanup/manifest.json) record the new checks separately from
initial port evidence. Source and rebuilt executable hashes accompany the records.
Raw per-case artifacts remain in
`/private/tmp/athenak-multiphase-curves-cleanup-20261008`.
No CUDA/HIP tests ran; the pending device qualification remains unchanged.

## Preserved failures and review

Earlier cohorts remain explicitly separate from the successful final cohort:

* `unsupported-pilot`: native kinematic initialization used LLF, whereas this
  branch requires `advect`; corrected input selects that solver before testing
  the curve rejection. No numerical tolerance changed.
* `production-all`: native AMR setup failed because no criterion block was
  supplied. Target requires `amr_criterion0/method=user` to invoke the enrolled
  fixture hook. AMR pilot and a complete fresh integrated cohort then pass.
* `branch-feature-off-sandbox-failure`: local MPI socket creation was denied by
  the sandbox. The bounded local run succeeded through authorized execution.
* `branch-feature-off-frame-pilot`: the short run reached only filter priming.
  An extended eight-cycle budget exercises actual boosts; unchanged-output
  checks are from that later run, not the incomplete coverage pilot.

Independent source review traced stage/ghost timing, AMR completion, null-source
handling, actual cooling status, frame/shearing exclusions and sidecar ordering.
Independent evidence review rehashed inputs/executables, all45 feature-off outputs
and six native-wave comparisons, and checked actual regrids/restarts. No material
source/evidence defect remains in the bounded scope. New C++ files pass the repo's
cpplint policy and custom length/tab/whitespace checks; Python syntax passes.
Flake8 is unavailable locally and was not reported passed. Native-tree pre-existing
lint issues were not changed.

## Retention and limits

Captured generated headers/build logs retain their original whitespace for byte
provenance. A directory-local Git whitespace attribute applies only to this
evidence bundle; runtime/test source checks remain unchanged.

Exact text inputs, commands, logs, numerical result records and configuration/
provenance hashes are committed here. Heavy native payloads/curve snapshots and
duplicate source tarballs remain under
`/private/tmp/athenak-multiphase-curves-validation`; the committed tests/examples
regenerate them. The complete small I/O source/config bundle is retained. Original
absolute artifact paths in raw records are preserved as execution provenance.

These are bounded software controls and actual unforced-fluid checks, not a
scientifically qualified turbulence release. Home-rank curve storage, host
geometry/adaptation/diagnostics, host-buffer MPI and ASCII I/O remain explicit
implementation ceilings; no new cluster scaling claim is made.
CUDA/HIP hardware is unavailable; the operator will perform device qualification
later using the opt-in procedure in [the feature guide](../material_curves.md).
Forcing/cooling/frame/MC coexistence and RK3 remain rejected extensions.
No merge, push, publication or substantial remote compute was performed.
